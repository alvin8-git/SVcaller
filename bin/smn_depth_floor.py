#!/usr/bin/env python3
"""Depth floor for the SMN1/SMN2 copy-number result.

Below `--min-depth` mean coverage the SMNCopyNumberCaller can return a wrong SMN1/SMN2
split that still carries a PASS tag (a GIAB HG002 downsampling series gave 3/1 for a
true 2/2 at 14.7x).  A result that is not confident is not reported as a call, so under
the floor the row becomes an explicit no-call in the shape of the caller's own low-depth
no-call: isSMA, isCarrier, SMN1_CN and SMN2_CN = None, Info = NoCall:BelowValidatedDepth.
The other columns (Total_CN_raw, SMN1_CN_raw, ...) are left as the caller wrote them.

  at/above the floor            output is a byte-identical copy of the input
  below the floor               floored row written; the untouched input is saved to --raw
  caller's own `NoCall:*` row   passed through unchanged (never re-labelled)
  our NoCall + raw copy exists  re-evaluated from the raw copy (fix a wrong depth / new floor
                                by re-running; above the floor the raw row is restored)
  bad depth / bad TSV           error, exit 2, nothing written (no depth, nan/inf, malformed
                                `total` row, no SMN1_CN/SMN2_CN column)

Depth: `--depth X`, or `--summary <sample>.mosdepth.summary.txt` (the file the pipeline's
MOSDEPTH step already writes; mean = column 4 of the first `total` row).

Usage: smn_depth_floor.py <in.smn.tsv> (--depth X | --summary S) [--min-depth 18]
                          [--out OUT] [--raw RAW]
  --out defaults to rewriting <in> in place (--out must not resolve to <in> through a
  symlink in a pipeline: stage the input under another name). --raw defaults to
  <S>.smn.raw.tsv next to <in> (only written when the floor applies, never overwritten).
"""
import argparse
import math
import os
import shutil
import sys

NOCALL = "NoCall:BelowValidatedDepth"
BLANK = ("isSMA", "isCarrier", "SMN1_CN", "SMN2_CN")


def mean_depth(summary):
    with open(summary) as fh:
        for line in fh:
            f = line.rstrip("\r\n").split("\t")
            if f and f[0] == "total":
                if len(f) < 4:
                    raise ValueError(f"malformed 'total' row in {summary}")
                return float(f[3])
    raise ValueError(f"no 'total' row in {summary}")


def _rows(lines):
    idx = [i for i, l in enumerate(lines) if l.strip() and not l.startswith("#")]
    if len(idx) < 2:
        raise ValueError("SMN TSV has no header + sample row")
    header = lines[idx[0]].split("\t")
    for c in ("SMN1_CN", "SMN2_CN"):
        if c not in header:
            raise ValueError(f"SMN TSV has no {c} column")
    return idx, header


def info_of(text):
    lines = text.replace("\r\n", "\n").split("\n")
    idx, header = _rows(lines)
    if "Info" not in header:
        return ""
    row = lines[idx[1]].split("\t")
    k = header.index("Info")
    return row[k] if k < len(row) else ""


def floor_text(text, depth, min_depth):
    """Return (new_text, floored). Pure function on the TSV text."""
    if not math.isfinite(depth):
        raise ValueError(f"non-finite depth {depth}")
    crlf = "\r\n" in text
    lines = text.replace("\r\n", "\n").split("\n")
    idx, header = _rows(lines)
    if depth >= min_depth:
        return text, False
    changed = False
    for i in idx[1:]:
        row = lines[i].split("\t")
        row += [""] * (len(header) - len(row))
        d = dict(zip(header, row))
        if d.get("Info", "").startswith("NoCall:"):
            continue  # already a no-call: keep the caller's own reason
        for k in BLANK:
            if k in header:
                row[header.index(k)] = "None"
        if "Info" in header:
            row[header.index("Info")] = NOCALL
        else:
            row.append(NOCALL)
        lines[i] = "\t".join(row)
        changed = True
    if changed and "Info" not in header:
        lines[idx[0]] = lines[idx[0]] + "\tInfo"
    out = "\n".join(lines)
    return (out.replace("\n", "\r\n") if crlf else out), changed


def default_raw(tsv):
    return (tsv[:-len(".smn.tsv")] if tsv.endswith(".smn.tsv") else tsv) + ".smn.raw.tsv"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("tsv")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--depth", type=float)
    g.add_argument("--summary")
    ap.add_argument("--min-depth", type=float, default=18.0)
    ap.add_argument("--out")
    ap.add_argument("--raw")
    a = ap.parse_args(argv)
    out = a.out or a.tsv
    raw = a.raw or default_raw(a.tsv)
    try:
        depth = a.depth if a.depth is not None else mean_depth(a.summary)
        with open(a.tsv, newline="") as fh:
            text = fh.read()
        if info_of(text) == NOCALL and os.path.exists(raw):
            with open(raw, newline="") as fh:
                text = fh.read()  # re-evaluate from the caller's own row
        new, floored = floor_text(text, depth, a.min_depth)
    except (OSError, ValueError) as e:
        print(f"ERROR: smn_depth_floor: {e}", file=sys.stderr)
        return 2
    if floored:
        if not os.path.exists(raw):
            with open(raw, "w", newline="") as fh:
                fh.write(text)
        with open(out, "w", newline="") as fh:
            fh.write(new)
        print(f"SMN depth floor: mean depth {depth:g}x < {a.min_depth:g}x -> {NOCALL}")
    else:
        with open(out, "w", newline="") as fh:
            fh.write(new)
        print(f"SMN depth floor: mean depth {depth:g}x, row unchanged")
    return 0


if __name__ == "__main__":
    sys.exit(main())
