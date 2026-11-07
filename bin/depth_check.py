#!/usr/bin/env python3
"""Read-depth consistency check for DEL / DUP calls in a <sample>.filtered.tsv.

Annotation only: no row is dropped.  Appends seven columns after the existing ones:
  DEPTH_IN            mean depth inside the SV, breakpoints excluded (EDGE_EXCLUDE bp each side)
  DEPTH_FLANK         mean depth of the two flanks, each min(SV length, FLANK_MAX) bp
  DEPTH_GLOBAL        the sample's autosomal median depth (same value on every row)
  DEPTH_RATIO         DEPTH_IN / DEPTH_FLANK, 2 decimals; "." when the check does not apply
  DEPTH_RATIO_GLOBAL  DEPTH_IN / DEPTH_GLOBAL, 2 decimals
  DEPTH_SUPPORT       supported | refuted | inconclusive | na
  DEPTH_NOTE          "flank inflated" / "flank depressed" (segdup signature) or "."

Two baselines, because a segmental duplication inflates the local flank (a normal-depth
region then looks deleted).  "supported" needs BOTH ratios in the supported band; "refuted"
needs EITHER in the refuted band; otherwise "inconclusive".  If the flank is too shallow
for a ratio, the global ratio alone decides.

A caller-supported DEL that reads at ~1.0x of its flanks is not a deletion (het DEL ~0.5x,
hom ~0); a DUP at ~1.0x is not a duplication (het DUP ~1.5x).  "refuted" means the read
depth contradicts the call; it is evidence against the call, not proof it is false (e.g.
segmental-duplication regions give ratios between 0.6 and 0.9).

Depth comes from `samtools coverage -r` (no new dependency; samtools is in the default
utils container).  Default samtools filters apply (unmapped, secondary, QC-fail, duplicates
are excluded).

Usage: depth_check.py <filtered.tsv> <bam> <out.tsv> [--threads N]
"""
import argparse
import csv
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

# ---- thresholds: the single config block ------------------------------------------------
# Calls shorter than MIN_SV_LEN are within depth noise at ~30x, so the check does not apply.
# Flank depth below MIN_FLANK_DEPTH is too sparse for a ratio.  Between the "supported" and
# "refuted" cut-offs the ratio is "inconclusive" (het events and segdups land here).
MIN_SV_LEN = 1000
EDGE_EXCLUDE = 500
FLANK_MAX = 100_000
MIN_FLANK_DEPTH = 5.0
DEL_SUPPORTED_MAX, DEL_REFUTED_MIN = 0.65, 0.85
DUP_SUPPORTED_MIN, DUP_REFUTED_MAX = 1.35, 1.15
# Global baseline: median of samtools-coverage mean depth over GLOBAL_WINDOWS windows of
# GLOBAL_WINDOW_LEN bp at fixed, evenly spaced positions along the concatenated chr1-22
# (deterministic, not random); zero-coverage windows are skipped.  Cached in a sidecar
# "<out_tsv>.global_depth" so re-annotation does not recompute.
GLOBAL_WINDOWS = 300
GLOBAL_WINDOW_LEN = 20_000
# Local flank vs global depth: outside this band the flank is a segdup-like outlier.
FLANK_INFLATED_MIN, FLANK_DEPRESSED_MAX = 1.3, 0.7
# -----------------------------------------------------------------------------------------

NEW_COLS = ["DEPTH_IN", "DEPTH_FLANK", "DEPTH_GLOBAL", "DEPTH_RATIO",
            "DEPTH_RATIO_GLOBAL", "DEPTH_SUPPORT", "DEPTH_NOTE"]


def band(sv_type, ratio):
    """One ratio -> supported/refuted/inconclusive (None for non DEL/DUP or missing)."""
    if ratio is None:
        return None
    if sv_type == "DEL":
        if ratio <= DEL_SUPPORTED_MAX:
            return "supported"
        return "refuted" if ratio >= DEL_REFUTED_MIN else "inconclusive"
    if sv_type == "DUP":
        if ratio >= DUP_SUPPORTED_MIN:
            return "supported"
        return "refuted" if ratio <= DUP_REFUTED_MAX else "inconclusive"
    return None


def depth_support(sv_type, ratio, ratio_global=None):
    """Flank ratio + global ratio -> supported/refuted/inconclusive/na.

    supported only if every available ratio is supported; refuted if any is refuted."""
    bands = [x for x in (band(sv_type, ratio), band(sv_type, ratio_global)) if x]
    if not bands:
        return "na"
    if "refuted" in bands:
        return "refuted"
    return "supported" if all(x == "supported" for x in bands) else "inconclusive"


def depth_note(d_flank, d_global):
    if d_flank is None or not d_global:
        return "."
    r = d_flank / d_global
    if r >= FLANK_INFLATED_MIN:
        return "flank inflated"
    if r <= FLANK_DEPRESSED_MAX:
        return "flank depressed"
    return "."


def global_depth(bam, lens):
    """Median autosomal depth over fixed evenly spaced windows (see config block)."""
    autos = []
    for c, n in lens.items():
        base = c[3:] if c.startswith("chr") else c
        if base.isdigit() and 1 <= int(base) <= 22:
            autos.append((int(base), c, n))
    autos.sort()
    total = sum(n for _, _, n in autos)
    vals = []
    for i in range(GLOBAL_WINDOWS):
        pos = int((i + 0.5) * total / GLOBAL_WINDOWS)
        for _, c, n in autos:
            if pos < n:
                s = max(1, min(pos, n - GLOBAL_WINDOW_LEN + 1))
                ws, wl = mean_depth(bam, c, s, s + GLOBAL_WINDOW_LEN - 1)
                if wl and ws > 0:
                    vals.append(ws / wl)
                break
            pos -= n
    if not vals:
        raise SystemExit("depth_check: no covered autosomal windows; cannot estimate global depth")
    vals.sort()
    m = len(vals) // 2
    return vals[m] if len(vals) % 2 else (vals[m - 1] + vals[m]) / 2


def depth_ratio(d_in, d_flank):
    """Ratio or None when the flank is too shallow / a depth is missing."""
    if d_in is None or d_flank is None or d_flank < MIN_FLANK_DEPTH:
        return None
    return d_in / d_flank


def regions(start, end, chrom_len):
    """(interior, [flanks]) as 1-based inclusive (s, e) tuples; interior None if empty."""
    length = end - start + 1
    flank = min(length, FLANK_MAX)
    interior = (start + EDGE_EXCLUDE, end - EDGE_EXCLUDE)
    if interior[0] > interior[1]:
        interior = None
    left = (max(1, start - flank), start - 1)
    right = (end + 1, min(chrom_len, end + flank)) if chrom_len else (end + 1, end + flank)
    return interior, [r for r in (left, right) if r[0] <= r[1]]


def mean_depth(bam, contig, s, e):
    """(length-weighted sum of depth, length) for a region via samtools coverage."""
    out = subprocess.run(["samtools", "coverage", "-r", f"{contig}:{s}-{e}", bam],
                         capture_output=True, text=True, check=True).stdout
    for line in out.splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        span = int(f[2]) - int(f[1]) + 1
        return float(f[6]) * span, span
    return 0.0, 0


def contig_lengths(bam):
    out = subprocess.run(["samtools", "view", "-H", bam], capture_output=True,
                         text=True, check=True).stdout
    lens = {}
    for line in out.splitlines():
        if line.startswith("@SQ"):
            kv = dict(x.split(":", 1) for x in line.split("\t")[1:])
            lens[kv["SN"]] = int(kv["LN"])
    return lens


def resolve_contig(chrom, lens):
    for c in (chrom, "chr" + chrom, chrom[3:] if chrom.startswith("chr") else None):
        if c and c in lens:
            return c
    return None


def check_one(bam, contig, start, end, chrom_len):
    """-> (depth_in, depth_flank) floats or None."""
    interior, flanks = regions(start, end, chrom_len)
    if interior is None or not flanks:
        return None, None
    ws, wl = mean_depth(bam, contig, *interior)
    d_in = ws / wl if wl else None
    fs = fl = 0
    for r in flanks:
        s, l = mean_depth(bam, contig, *r)
        fs += s
        fl += l
    return d_in, (fs / fl if fl else None)


def annotate(rows, bam, threads=4, d_global=None):
    """rows: list of dict (DictReader rows). Adds NEW_COLS in place; returns rows."""
    lens = contig_lengths(bam)
    if d_global is None:
        d_global = global_depth(bam, lens)
    keys = {}
    for r in rows:
        if r["SV_type"] not in ("DEL", "DUP"):
            continue
        try:
            s, e = int(r["SV_start"]), int(r["SV_end"])
        except ValueError:
            continue
        if e - s + 1 < MIN_SV_LEN:
            continue
        c = resolve_contig(r["SV_chrom"], lens)
        if c:
            keys[(r["SV_chrom"], s, e)] = c

    def work(k):
        chrom, s, e = k
        return k, check_one(bam, keys[k], s, e, lens[keys[k]])

    with ThreadPoolExecutor(max(1, threads)) as ex:
        res = dict(ex.map(work, list(keys)))

    for r in rows:
        r.update({c: "." for c in NEW_COLS})
        r["DEPTH_SUPPORT"] = "na"
        r["DEPTH_GLOBAL"] = f"{d_global:.2f}"
        try:
            k = (r["SV_chrom"], int(r["SV_start"]), int(r["SV_end"]))
        except ValueError:
            continue
        if k not in res:
            continue
        d_in, d_fl = res[k]
        if d_in is not None:
            r["DEPTH_IN"] = f"{d_in:.2f}"
        if d_fl is not None:
            r["DEPTH_FLANK"] = f"{d_fl:.2f}"
        ratio = depth_ratio(d_in, d_fl)
        ratio_g = d_in / d_global if d_in is not None and d_global > 0 else None
        if ratio is not None:
            r["DEPTH_RATIO"] = f"{ratio:.2f}"
        if ratio_g is not None:
            r["DEPTH_RATIO_GLOBAL"] = f"{ratio_g:.2f}"
        r["DEPTH_SUPPORT"] = depth_support(r["SV_type"], ratio, ratio_g)
        r["DEPTH_NOTE"] = depth_note(d_fl, d_global)
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("filtered_tsv")
    ap.add_argument("bam")
    ap.add_argument("out_tsv")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--global-depth", type=float, default=None,
                    help="override/skip the autosomal median estimate")
    a = ap.parse_args(argv)
    with open(a.filtered_tsv, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        fields = list(rd.fieldnames or [])
        rows = list(rd)
    if "DEPTH_SUPPORT" in fields:   # re-annotation: replace, keep column order
        fields = [c for c in fields if c not in NEW_COLS]
    sidecar = a.out_tsv + ".global_depth"
    d_global = a.global_depth
    if d_global is None:
        try:
            d_global = float(open(sidecar).read().split()[0])
        except (OSError, ValueError, IndexError):
            d_global = global_depth(a.bam, contig_lengths(a.bam))
            with open(sidecar, "w") as fh:
                fh.write(f"{d_global:.4f}\n")
    print(f"depth_check: DEPTH_GLOBAL (autosomal median) = {d_global:.2f}x", file=sys.stderr)
    annotate(rows, a.bam, a.threads, d_global)
    with open(a.out_tsv, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=fields + NEW_COLS, delimiter="\t",
                           extrasaction="ignore", lineterminator="\r\n")
        w.writeheader()
        w.writerows(rows)
    n = {}
    for r in rows:
        n[r["DEPTH_SUPPORT"]] = n.get(r["DEPTH_SUPPORT"], 0) + 1
    print("depth_check:", ", ".join(f"{k}={v}" for k, v in sorted(n.items())), file=sys.stderr)


if __name__ == "__main__":
    main()
