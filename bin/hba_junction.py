#!/usr/bin/env python3
"""Channel 3 (HBA_JUNCTION) — targeted split-read / discordant-pair junction
detection in the alpha-globin cluster.

WHAT THIS MEASURES
------------------
Given a BAM and a small region, find deletion junctions by stacking soft-clip
positions, pair a left breakend with a right breakend, and — critically —
report the ALLELE BALANCE at the junction so the call can be typed het vs hom.

    split    = reads clipped exactly at a breakend (variant-allele evidence)
    spanning = unclipped reads crossing the left breakend (reference-allele
               evidence: they can only come from the INTACT chromosome)
    vaf      = split / (split + spanning)

A caller that only counts split reads cannot tell a heterozygous 20 kb
deletion from a homozygous one, and on this locus that difference is the whole
clinical question. Spanning reads are not an optional extra; they ARE the
zygosity measurement.

WHAT THIS DOES NOT DO
---------------------
It does NOT name an allele. It reports coordinates and support. Mapping a
junction onto `--SEA` vs `--MED` requires published GRCh38 breakpoints for both,
which this repo does not have (see the header of
`assets/hba_deletion_alleles.tsv`: "THESE ARE NOT BREAKPOINTS"). Naming an
allele from a bare coordinate — or from the sample's likely ancestry — is a
population inference dressed as a measurement. Don't.

USE `confidence == "high"` ONLY. Both real samples carry recurrent soft-clip
stacks in this locus (chr16 ~182000, ~188350, ~188365) that look like junctions
on split reads alone and appear even in THAL2, which has no deletion at all.
Discordant-pair corroboration is what separates signal from those artifacts, and
the confidence tier encodes it. Rows are still emitted for transparency; do not
promote a `medium`/`low` row into an allele call.

It also does not realign the clipped sequence. The clipped bases of a real
junction read are reference sequence from the far side of the deletion, so
realigning them would confirm the partner breakend independently. Here the
partner is inferred only from the opposing clip stack.

Pure stdlib + subprocess calls to `samtools`. No pysam: the `svcaller/utils:1.2`
container does not have it.

Usage:
    hba_junction.py --bam S.bam --sample S --out S.alpha_junction.tsv \\
                    [--region chr16:1-250000] [--min-split 3] \\
                    [--min-clip 15] [--min-mapq 1] [--window 10]

An output file with only the header row is a valid NEGATIVE result ("no
junction detected"), not an error. Exit status stays 0.
"""
import argparse
import os
import re
import subprocess
import sys
from collections import Counter

# ---------------------------------------------------------------------------
# Defaults, and why they are what they are.
# ---------------------------------------------------------------------------

# Breakend clustering window (bp). True junction reads clip at the SAME base
# (the fixture's six left splits all end at exactly one position). Real BAMs
# jitter a few bp: alpha-cluster deletions arise by NAHR inside near-identical
# homology boxes, so the aligner has several equally-good places to put the
# breakend. 10 bp absorbs that jitter while staying two orders of magnitude
# below the smallest allele of interest (-a3.7 = 3.7 kb), so two genuinely
# distinct junctions can never be merged into one.
DEFAULT_WINDOW = 10

# Minimum clipped bases for a clip to count as junction evidence. Clips of a
# few bp are produced constantly by adapter remnants and terminal mismatches.
# 15 bp is well clear of that noise floor and long enough that the clipped
# sequence is, in principle, realignable.
DEFAULT_MIN_CLIP = 15

# Minimum reads in a breakend cluster. At 30x with 150 bp reads a heterozygous
# junction is covered by ~15x on the variant allele, and a read votes for the
# breakend if its start falls in the (read_len - 2*min_clip) window around it,
# so E[split] ~ 15 * 120/150 = 12 per side. 3 is therefore far below what a
# real het deletion produces at target coverage, while still being above the
# 1-2 stray clipped reads that mismapping piles up at an arbitrary base.
DEFAULT_MIN_SPLIT = 3

# Minimum MAPQ. MAPQ 0 (multi-mapping) reads are the dominant artifact source
# in this segmentally-duplicated locus and are dropped. Anything above 0 is
# kept: demanding MAPQ>=20 here would discard genuine reads, because the alpha
# cluster is repetitive by nature.
DEFAULT_MIN_MAPQ = 1

# Discordant-pair insert cutoff. Derived per-run from the data (see
# `insert_cutoff`), never hardcoded to a library-specific number.
MIN_INSERT_CUTOFF = 1000     # absolute floor, see insert_cutoff()
INSERT_CUTOFF_MULT = 3       # multiples of the observed median proper insert

# Allele-balance thresholds. With split counted on BOTH breakends and spanning
# counted at the left breakend only, a true het at 30x gives
# E[split] ~ 24, E[spanning] ~ 12 -> vaf ~ 0.67; a true hom gives spanning ~ 0
# -> vaf ~ 1.0. 0.85 sits between them with a wide margin on both sides
# (equivalently: to call hom you may tolerate at most ~1 stray reference read
# per 6 junction reads).
HOM_VAF = 0.85
# Below this the balance is inconsistent with a germline heterozygote (e.g.
# 3 split against 27 spanning is ~10% allele fraction: mosaicism or, far more
# likely here, mismapping). Refuse to type it rather than guess.
MIN_HET_VAF = 0.15

# Total junction+reference depth below which the zygosity call is not trusted
# enough to be called "high" confidence.
MIN_DEPTH_FOR_HIGH = 10

COLUMNS = ["sample", "chrom", "left_bp", "right_bp", "size", "split_left",
           "split_right", "discordant", "spanning", "vaf", "zygosity",
           "confidence"]
HEADER = "#" + "\t".join(COLUMNS)

# SAM flag bits we refuse outright.
FLAG_UNMAPPED = 0x4
FLAG_SECONDARY = 0x100
FLAG_DUPLICATE = 0x400
FLAG_SUPPLEMENTARY = 0x800
# Supplementary alignments are the *other* half of a split read. Their primary
# already carries the soft clip and already votes, so counting them too would
# double-count one fragment as two independent observations.
FLAG_EXCLUDE = FLAG_UNMAPPED | FLAG_SECONDARY | FLAG_DUPLICATE | FLAG_SUPPLEMENTARY

FLAG_PAIRED = 0x1
FLAG_PROPER = 0x2
FLAG_MATE_UNMAPPED = 0x8

CIGAR_OP = re.compile(r"(\d+)([MIDNSHP=X])")
# Operations that advance the position on the REFERENCE.
REF_CONSUMING = frozenset("MDN=X")
REGION_RE = re.compile(r"^([^:]+)(?::([\d,]+)(?:-([\d,]+))?)?$")


# ---------------------------------------------------------------------------
# Pure functions (unit-testable without a BAM)
# ---------------------------------------------------------------------------

def parse_cigar(cigar):
    """['150M'] -> [(150, 'M')]. Raises ValueError on a malformed CIGAR."""
    if not cigar or cigar == "*":
        return []
    ops = [(int(n), op) for n, op in CIGAR_OP.findall(cigar)]
    if "".join(f"{n}{op}" for n, op in ops) != cigar:
        raise ValueError(f"malformed CIGAR: {cigar!r}")
    return ops


def ref_span(cigar):
    """Bases of REFERENCE consumed by an alignment.

    M/D/N/=/X consume reference; I/S/H/P do not. Summing every number in the
    CIGAR is a real and common bug: '10S100M20I30D' spans 130, not 160, and a
    caller using the naive sum places the breakend 30 bp off.
    """
    return sum(n for n, op in parse_cigar(cigar) if op in REF_CONSUMING)


def aln_end(pos, cigar):
    """1-based coordinate of the LAST aligned base. 0 if nothing aligns."""
    span = ref_span(cigar)
    return pos + span - 1 if span else 0


def clip_ends(pos, cigar, min_clip=1):
    """Breakend candidates implied by this read's soft clips.

    Returns ``(leading_clip_bp, trailing_clip_bp)``, either of which may be
    ``None``. Read the names as "the breakpoint implied by the clip at that END
    OF THE READ" — note the deliberate inversion:

      * a LEADING ``<n>S`` means the read's start is unanchored, so the
        alignment begins abruptly at ``POS``: that is the RIGHT breakend of a
        deletion (the distal side).
      * a TRAILING ``<n>S`` means the alignment stops abruptly at the last
        aligned base: that is the LEFT breakend (the proximal side).

    Hard clips flanking the soft clip are tolerated (bwa emits ``30H20S100M``).
    A read that is entirely clipped, or aligns nothing, yields nothing.
    """
    ops = parse_cigar(cigar)
    span = sum(n for n, op in ops if op in REF_CONSUMING)
    if span <= 0:
        return None, None

    i = 0
    while i < len(ops) and ops[i][1] == "H":
        i += 1
    j = len(ops) - 1
    while j >= 0 and ops[j][1] == "H":
        j -= 1

    leading = None
    if i <= j and ops[i][1] == "S" and ops[i][0] >= min_clip:
        leading = pos
    trailing = None
    if j > i and ops[j][1] == "S" and ops[j][0] >= min_clip:
        trailing = pos + span - 1
    return leading, trailing


def cluster(positions, window=DEFAULT_WINDOW):
    """Group nearby breakend positions.

    Single-linkage on the sorted positions: a gap larger than `window` starts a
    new cluster. Each cluster is reported as ``(representative, count)`` where
    the representative is the MODE, not the mean. That is what makes
    base-resolution recovery possible: the true breakpoint is the position most
    reads clip at, whereas a mean is dragged off-base by jitter and is not even
    guaranteed to be an integer.

    Caveat: single linkage can chain (positions 0,10,20,... all merge at
    window=10). Harmless at window=10 given real alleles are kilobases apart,
    but do not raise `window` into the hundreds without revisiting this.
    """
    positions = sorted(positions)
    if not positions:
        return []
    groups, cur = [], [positions[0]]
    for p in positions[1:]:
        if p - cur[-1] <= window:
            cur.append(p)
        else:
            groups.append(cur)
            cur = [p]
    groups.append(cur)

    out = []
    for g in groups:
        counts = Counter(g)
        # highest count wins; ties break to the leftmost position for determinism
        rep = min(counts, key=lambda p: (-counts[p], p))
        out.append((rep, len(g)))
    return out


def zygosity_from_balance(split, spanning):
    """Allele balance -> (vaf, call). call is 'het' | 'hom' | 'NA'.

    THE POINT OF THIS MODULE. `spanning` reads cross the left breakend without
    a clip, which is only possible on a chromosome that still has the deleted
    interval. spanning ~ split  => het.  spanning ~ 0  => hom.
    """
    total = split + spanning
    if total <= 0:
        return None, "NA"
    vaf = round(split / total, 3)
    if vaf >= HOM_VAF:
        return vaf, "hom"
    if vaf < MIN_HET_VAF:
        return vaf, "NA"
    return vaf, "het"


def insert_cutoff(proper_inserts):
    """|TLEN| above which a pair is treated as a discordant candidate.

    Derived from the sample's own proper pairs rather than hardcoded, because
    libraries differ. A PCR-free Illumina insert distribution has sd ~ 0.15 x
    median, so 3 x median is roughly median + 13 sd — no normal pair reaches it.
    The 1000 bp floor stops a degenerate estimate (a region where every TLEN is
    0, or a handful of tiny inserts) from producing an absurdly low cutoff.
    """
    if proper_inserts:
        vals = sorted(proper_inserts)
        n = len(vals)
        median = vals[n // 2] if n % 2 else (vals[n // 2 - 1] + vals[n // 2]) / 2
    else:
        median = 500          # nothing to measure; assume a typical WGS library
    return max(MIN_INSERT_CUTOFF, int(INSERT_CUTOFF_MULT * median))


def confidence(split_left, split_right, discordant, spanning, zyg,
               min_split=DEFAULT_MIN_SPLIT):
    """Confidence tier. ONLY `high` should be used to collapse a degenerate
    allele group (see module docstring).

    The ladder is built around discordant corroboration because that is what
    empirically separates real junctions from artifacts at this locus. In
    THAL1 the true ~19.3 kb deletion carries 10+ discordant fragments whose
    inserts match its size, while the recurrent clip stacks at chr16:182xxx /
    188350 / 188365 — which appear in THAL2, a sample with NO deletion —
    carry none. Split reads alone do not distinguish a junction from a
    reference-mismatch pileup in a segmental duplication.
    """
    if zyg == "NA":
        return "low"
    both_sides = min(split_left, split_right) >= min_split
    deep_enough = (split_left + split_right + spanning) >= MIN_DEPTH_FOR_HIGH
    if not both_sides:
        return "low"
    if discordant >= min_split and deep_enough:
        return "high"
    if discordant >= 1:
        return "medium"
    return "low"          # clip stack with no paired-end corroboration


def parse_region(region):
    """'chr16:1-250000' -> ('chr16', 1, 250000). End may be None."""
    m = REGION_RE.match(region.strip())
    if not m:
        raise ValueError(f"unparseable region: {region!r}")
    chrom, start, end = m.group(1), m.group(2), m.group(3)
    return (chrom,
            int(start.replace(",", "")) if start else 1,
            int(end.replace(",", "")) if end else None)


def format_row(sample, chrom, left_bp, right_bp, split_left, split_right,
               discordant, spanning, vaf, zyg, conf):
    size = right_bp - left_bp - 1
    return "\t".join([
        sample, chrom, str(left_bp), str(right_bp), str(size),
        str(split_left), str(split_right), str(discordant), str(spanning),
        "NA" if vaf is None else f"{vaf:.3f}", zyg, conf,
    ])


# ---------------------------------------------------------------------------
# BAM access
# ---------------------------------------------------------------------------

def read_sam_records(bam, region, min_mapq=DEFAULT_MIN_MAPQ):
    """Stream `samtools view bam region` once and keep only what we need.

    ALWAYS pass a region: this is a targeted caller, and a full BAM scan on a
    real 30x WGS file is tens of minutes of pointless IO.
    """
    cmd = ["samtools", "view", bam]
    if region:
        cmd.append(region)
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except FileNotFoundError:
        raise SystemExit("hba_junction: samtools not found on PATH")
    except subprocess.CalledProcessError as exc:
        raise SystemExit(
            f"hba_junction: samtools view failed on {bam} {region or ''}: "
            f"{exc.stderr.strip()[:400]}")

    records = []
    for line in proc.stdout.splitlines():
        if not line or line.startswith("@"):
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        flag = int(f[1])
        if flag & FLAG_EXCLUDE:
            continue
        if int(f[4]) < min_mapq:
            continue
        cigar = f[5]
        if cigar == "*":
            continue
        try:
            span = ref_span(cigar)
        except ValueError:
            continue
        if span <= 0:
            continue
        pos = int(f[3])
        records.append({
            "qname": f[0], "flag": flag, "chrom": f[2], "pos": pos,
            "cigar": cigar, "rnext": f[6], "pnext": int(f[7]),
            "tlen": int(f[8]), "end": pos + span - 1,
        })
    return records


def discordant_fragments(records, chrom, cutoff):
    """qname -> (fragment_start, fragment_end) for far-apart same-chrom pairs.

    Keyed by QNAME so a pair whose two mates both fall inside the region is ONE
    piece of evidence, not two.
    """
    frags = {}
    for r in records:
        flag = r["flag"]
        if not flag & FLAG_PAIRED:
            continue
        if flag & FLAG_PROPER:          # aligner already says this is normal
            continue
        if flag & FLAG_MATE_UNMAPPED:
            continue
        mate_chrom = r["chrom"] if r["rnext"] == "=" else r["rnext"]
        if mate_chrom != chrom:         # interchromosomal: not a deletion
            continue
        tlen = abs(r["tlen"])
        if tlen < cutoff:
            continue
        start = min(r["pos"], r["pnext"])
        frags[r["qname"]] = (start, start + tlen - 1)
    return frags


# ---------------------------------------------------------------------------
# Calling
# ---------------------------------------------------------------------------

def call_junctions(records, chrom, min_split=DEFAULT_MIN_SPLIT,
                   min_clip=DEFAULT_MIN_CLIP, window=DEFAULT_WINDOW):
    """records -> list of junction dicts. Empty list is a valid answer."""
    left_bps, right_bps, unclipped = [], [], []
    for r in records:
        lead, trail = clip_ends(r["pos"], r["cigar"], min_clip)
        if lead is not None:
            right_bps.append(lead)      # leading clip  -> RIGHT breakend
        if trail is not None:
            left_bps.append(trail)      # trailing clip -> LEFT breakend
        if lead is None and trail is None:
            unclipped.append(r)

    left_clusters = [(p, n) for p, n in cluster(left_bps, window) if n >= min_split]
    right_clusters = [(p, n) for p, n in cluster(right_bps, window) if n >= min_split]
    if not left_clusters or not right_clusters:
        return []

    cutoff = insert_cutoff([abs(r["tlen"]) for r in records
                            if r["flag"] & FLAG_PROPER and r["tlen"] != 0])
    frags = discordant_fragments(records, chrom, cutoff)

    # Enumerate every left/right pairing that describes a real deletion, score
    # it, then take them greedily so each breakend cluster is used at most once.
    candidates = []
    for li, (lbp, ln) in enumerate(left_clusters):
        for ri, (rbp, rn) in enumerate(right_clusters):
            if rbp <= lbp + 1:          # need at least one deleted base
                continue
            size = rbp - lbp - 1
            # A pair supports THIS junction only if it straddles both breakends
            # AND its insert is consistent with this deletion's size: collapsing
            # the deleted interval must leave a normal-looking fragment. Without
            # the size check a single chimeric pair with a 50 Mb TLEN (THAL1 has
            # one) "supports" every candidate junction in the region.
            disc = sum(1 for s, e in frags.values()
                       if s <= lbp and e >= rbp
                       and (e - s + 1) - size <= cutoff)
            candidates.append({
                "li": li, "ri": ri, "left_bp": lbp, "right_bp": rbp,
                "split_left": ln, "split_right": rn, "discordant": disc,
            })
    # Best-supported first: discordant corroboration, then split depth, then the
    # tighter (smaller) deletion, then leftmost — the last two only for
    # determinism when everything else ties.
    candidates.sort(key=lambda c: (-c["discordant"],
                                   -min(c["split_left"], c["split_right"]),
                                   c["right_bp"] - c["left_bp"],
                                   c["left_bp"]))

    used_l, used_r, out = set(), set(), []
    for c in candidates:
        if c["li"] in used_l or c["ri"] in used_r:
            continue
        used_l.add(c["li"])
        used_r.add(c["ri"])
        lbp = c["left_bp"]
        # Reference-allele evidence: unclipped reads whose aligned block
        # STRICTLY crosses the left breakend (start < bp < end). A read merely
        # abutting the breakend proves nothing about the intact allele.
        spanning = sum(1 for r in unclipped if r["pos"] < lbp < r["end"])
        split = c["split_left"] + c["split_right"]
        vaf, zyg = zygosity_from_balance(split, spanning)
        c["spanning"] = spanning
        c["vaf"] = vaf
        c["zygosity"] = zyg
        c["confidence"] = confidence(c["split_left"], c["split_right"],
                                     c["discordant"], spanning, zyg, min_split)
        out.append(c)

    out.sort(key=lambda c: (c["left_bp"], c["right_bp"]))
    return out


def run(bam, sample, out_path, region=None, min_split=DEFAULT_MIN_SPLIT,
        min_clip=DEFAULT_MIN_CLIP, min_mapq=DEFAULT_MIN_MAPQ,
        window=DEFAULT_WINDOW):
    """Full pipeline. Returns the list of junction dicts written."""
    chrom = parse_region(region)[0] if region else None
    records = read_sam_records(bam, region, min_mapq)
    if chrom is None:
        chroms = Counter(r["chrom"] for r in records)
        chrom = chroms.most_common(1)[0][0] if chroms else "NA"

    junctions = call_junctions(records, chrom, min_split, min_clip, window)

    out_dir = os.path.dirname(os.path.abspath(out_path))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write(HEADER + "\n")
        for j in junctions:
            fh.write(format_row(sample, chrom, j["left_bp"], j["right_bp"],
                                j["split_left"], j["split_right"],
                                j["discordant"], j["spanning"], j["vaf"],
                                j["zygosity"], j["confidence"]) + "\n")
    return junctions


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Targeted split-read/discordant-pair junction caller for "
                    "the alpha-globin cluster. Reports coordinates and allele "
                    "balance; it does NOT name a deletion allele.")
    ap.add_argument("--bam", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--region", default="chr16:1-250000",
                    help="always targeted; never scan the whole BAM "
                         "(default: %(default)s)")
    ap.add_argument("--min-split", type=int, default=DEFAULT_MIN_SPLIT,
                    dest="min_split")
    ap.add_argument("--min-clip", type=int, default=DEFAULT_MIN_CLIP,
                    dest="min_clip")
    ap.add_argument("--min-mapq", type=int, default=DEFAULT_MIN_MAPQ,
                    dest="min_mapq")
    ap.add_argument("--window", type=int, default=DEFAULT_WINDOW,
                    help="breakend clustering window in bp (default: %(default)s)")
    args = ap.parse_args(argv)

    junctions = run(args.bam, args.sample, args.out, args.region,
                    args.min_split, args.min_clip, args.min_mapq, args.window)
    sys.stderr.write(
        f"hba_junction: {args.sample} {args.region} -> {len(junctions)} "
        f"junction(s); header-only output means no junction detected\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
