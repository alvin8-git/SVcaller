#!/usr/bin/env python3
"""Generate the alpha-globin deletion-allele definitions.

WHY THIS LOOKS DIFFERENT FROM make_globin_panels.py
---------------------------------------------------
Site coordinates are derivable: an HGVS position plus a RefSeq gene model gives
exactly one answer, and the FASTA confirms it. Deletion breakpoints are not.
They come from the literature, are usually published against hg19 or as relative
distances, and the alpha cluster's NAHR breakpoints sit inside near-identical
homology boxes where "the" breakpoint is not even well defined for a given
allele. Hand-typing GRCh38 coordinates for them would be exactly the failure
mode the site-panel generator exists to prevent.

So this does NOT define alleles by breakpoint. It defines them by their
**copy-number signature over diagnostic segments** — which genes each allele
removes — because that is what a depth caller can actually measure, and because
the segment boundaries ARE derivable from the same RefSeq models.

Three provenance classes, carried in the `basis` column and never blurred:

  derived      segment coordinates, from RefSeq gene models (as reliable as the site panel)
  observed     measured in this project's own data, cited inline
  literature   allele composition (which genes are lost) — textbook, but the
               SIZES are approximate and must not be used as coordinates

    python3 bin/make_hba_deletion_alleles.py --outdir assets/
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hgvs_map import load_models, DEFAULT_BED  # noqa: E402

# Diagnostic segments. Boundaries are derived from the gene models below; the
# inter-genic segments are what distinguish a 2-gene deletion that spares HBZ
# (--SEA) from one that does not (--FIL/--THAI).
#                                              reliability, measured 2026-07-22
SEGMENTS = [
    ("HBZ",        "HBZ",  "needs_own_baseline",
     "zeta-globin; spared by --SEA/--MED, lost in --FIL/--THAI. GC-rich and "
     "subtelomeric, so depth runs LOW even when intact: measured 0.86 (THAL1) "
     "and 0.71 (THAL2) against a chr2 control, and THAL2 has NO deletion. A "
     "global 0.8 threshold FALSE-POSITIVES here. Needs a per-segment baseline "
     "from known-normal samples, never the genome-wide control ratio."),

    ("INTER_Z_A",  None,   "do_not_average",
     "between HBZ and HBA2. DO NOT use as a single averaged segment. It spans "
     "18 kb of partly-duplicated sequence where mapping piles up: 1 kb windows "
     "in THAL1 read 1.4-1.9 across chr16:155000-162000, which offsets the real "
     "deletion over 164000-172875 and averages to 0.99 - i.e. it reports "
     "'intact' while CONTAINING half of a --SEA deletion. Use mappable "
     "sub-windows or exclude it."),

    ("HBA2",       "HBA2", "good",
     "alpha-2. Measured 0.37 (THAL1, het loss) vs 0.81 (THAL2, intact)."),

    ("INTER_A2_A1", None,  "good",
     "between HBA2 and HBA1; deleted in -a3.7. Measured 0.50 (THAL1) vs 0.99 (THAL2)."),

    ("HBA1",       "HBA1", "good",
     "alpha-1. Measured 0.40 (THAL1, het loss) vs 0.90 (THAL2, intact)."),
]

# allele, class, approx size, (HBZ, HBA2, HBA1) copy change on the affected
# chromosome, net functional alpha genes lost, populations, basis, note.
#   0 = intact, -1 = lost, "h" = disrupted/hybrid (partial signal)
ALLELES = [
    ("--SEA",     "deletional", "~20 kb",
     (0, -1, -1), 2, "EAS,SEA", "observed",
     "OBSERVED in THAL1: normalized depth 0.45-0.60 across chr16:164000-186000, "
     "HBZ 1.00, flanks 1.05. Spares HBZ. Commonest 2-gene allele in SE Asia."),

    ("--MED",     "deletional", "~17 kb",
     (0, -1, -1), 2, "MED", "literature",
     "Removes both alpha genes, spares HBZ. Depth-degenerate with --SEA; "
     "separating them needs the extent or a junction read."),

    ("--FIL",     "deletional", "~31 kb",
     (-1, -1, -1), 2, "SEA(PHL)", "literature",
     "Larger; extends into HBZ. Loss of HBZ signal is the discriminator vs --SEA/--MED."),

    ("--THAI",    "deletional", "~34 kb",
     (-1, -1, -1), 2, "SEA(THA)", "literature",
     "As --FIL, also removes HBZ. Depth-degenerate with --FIL."),

    ("-a3.7",     "deletional", "~3.7 kb",
     (0, "h", "h"), 1, "global,AFR,MED,SEA", "literature",
     "Rightward NAHR between Z boxes; fuses HBA2 5' to HBA1 3' into ONE hybrid "
     "gene. Neither gene body vanishes cleanly - the deleted 3.7 kb lies between "
     "them, so INTER_A2_A1 is the diagnostic segment, not HBA1/HBA2. Most common "
     "alpha-thal deletion worldwide."),

    ("-a4.2",     "deletional", "~4.2 kb",
     (0, -1, 0), 1, "SEA,PAC", "literature",
     "Leftward NAHR between X boxes; removes HBA2, leaves HBA1 intact."),

    ("anti-3.7",  "triplication", "+3.7 kb",
     (0, "+", "+"), -1, "global", "literature",
     "Reciprocal product of the -a3.7 NAHR: alpha-gene TRIPLICATION (aaa/). Gains "
     "a gene rather than losing one. Benign alone, but modifies beta-thal severity, "
     "and a caller that only looks for losses will silently miss it."),
]

HEADER = ["allele", "class", "approx_size", "d_HBZ", "d_HBA2", "d_HBA1",
          "alpha_genes_lost", "population", "basis", "depth_distinguishable", "note"]


def build_segments(models):
    z, a2, a1 = models["HBZ"], models["HBA2"], models["HBA1"]
    g = {"HBZ": (z.tx_start + 1, z.tx_end),
         "HBA2": (a2.tx_start + 1, a2.tx_end),
         "HBA1": (a1.tx_start + 1, a1.tx_end)}
    rows = []
    for name, gene, reliability, note in SEGMENTS:
        if gene:
            s, e = g[gene]
        elif name == "INTER_Z_A":
            s, e = g["HBZ"][1] + 1, g["HBA2"][0] - 1
        elif name == "INTER_A2_A1":
            s, e = g["HBA2"][1] + 1, g["HBA1"][0] - 1
        rows.append((z.chrom, s, e, name, reliability, note))
    return sorted(rows, key=lambda r: r[1])


def distinguishability():
    """Which alleles a DEPTH-only signature can separate. Groups sharing a
    signature are degenerate: the caller must report the group, not pick one."""
    groups = {}
    for a in ALLELES:
        groups.setdefault((a[3], a[4]), []).append(a[0])
    return {name: ("yes" if len(members) == 1 else "no:" + "|".join(members))
            for members in groups.values() for name in members}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", default=DEFAULT_BED)
    ap.add_argument("--outdir", default="assets")
    a = ap.parse_args()

    models = load_models(a.bed, {"HBZ", "HBA1", "HBA2"})
    segs = build_segments(models)

    for chrom, s, e, name, _, _ in segs:
        if e <= s:
            raise SystemExit(f"FATAL: segment {name} is empty or inverted ({s}-{e})")

    bed = os.path.join(a.outdir, "hba_segments.bed")
    with open(bed, "w") as fh:
        fh.write("# Diagnostic segments for alpha-globin copy-number calling — GRCh38\n")
        fh.write("# GENERATED by bin/make_hba_deletion_alleles.py; boundaries DERIVED from\n")
        fh.write("# RefSeq gene models (AnnotSV bundle). Do not hand-edit.\n")
        fh.write("#\n")
        fh.write("# col5 = reliability. NOT ALL SEGMENTS ARE USABLE THE SAME WAY — measured\n")
        fh.write("# on THAL1 (--SEA het) and THAL2 (no deletion) 2026-07-22:\n")
        fh.write("#   good               ratio vs a genome-wide control separates het loss cleanly\n")
        fh.write("#   needs_own_baseline intact depth is systematically below 1.0; a global\n")
        fh.write("#                      threshold false-positives. Calibrate per segment.\n")
        fh.write("#   do_not_average     averaging the whole segment is misleading; mapping\n")
        fh.write("#                      inflation can cancel a real deletion out.\n")
        fh.write("#\n")
        fh.write("# A caller that treats all five identically WILL emit wrong calls. See the\n")
        fh.write("# per-segment notes.\n")
        for chrom, s, e, name, reliability, note in segs:
            fh.write(f"{chrom}\t{s-1}\t{e}\t{name}\t{reliability}\t# {note}\n")

    dist = distinguishability()
    tsv = os.path.join(a.outdir, "hba_deletion_alleles.tsv")
    with open(tsv, "w") as fh:
        fh.write("# Alpha-globin deletion/triplication alleles — copy-number signatures\n")
        fh.write("# GENERATED by bin/make_hba_deletion_alleles.py. Do not hand-edit.\n")
        fh.write("#\n")
        fh.write("# THESE ARE NOT BREAKPOINTS. Alpha-cluster NAHR breakpoints sit inside\n")
        fh.write("# near-identical homology boxes and are published against varying builds;\n")
        fh.write("# hand-typing GRCh38 coordinates for them would be guesswork. Alleles are\n")
        fh.write("# therefore defined by WHICH SEGMENTS THEY REMOVE (see hba_segments.bed),\n")
        fh.write("# which is what a depth caller measures anyway.\n")
        fh.write("#\n")
        fh.write("# approx_size is DOCUMENTARY ONLY. Never use it as a coordinate.\n")
        fh.write("#\n")
        fh.write("# basis:  observed   = measured in this project's data, cited in note\n")
        fh.write("#         literature = allele composition is textbook; SIZE is approximate\n")
        fh.write("#\n")
        fh.write("# d_* columns are the copy change on the AFFECTED chromosome:\n")
        fh.write("#   0 intact · -1 lost · h disrupted/hybrid (partial depth signal) · + gained\n")
        fh.write("#\n")
        fh.write("# depth_distinguishable=no:A|B means depth ALONE cannot separate those\n")
        fh.write("# alleles. The caller MUST report the group, not arbitrarily pick one.\n")
        fh.write("# Resolving within a group needs the deletion extent or a junction read.\n")
        fh.write("#\n")
        fh.write("\t".join(HEADER) + "\n")
        for allele, cls, size, (dz, d2, d1), lost, pop, basis, note in ALLELES:
            fh.write("\t".join([allele, cls, size, str(dz), str(d2), str(d1),
                                str(lost), pop, basis, dist[allele], note]) + "\n")

    print(f"hba_segments.bed          {len(segs)} segments (derived)")
    for chrom, s, e, name, rel, _ in segs:
        flag = "" if rel == "good" else f"   <-- {rel}"
        print(f"  {name:<12} {chrom}:{s}-{e}  ({e-s+1:,} bp){flag}")
    print(f"\nhba_deletion_alleles.tsv  {len(ALLELES)} alleles")
    degenerate = sorted({v for v in dist.values() if v.startswith("no:")})
    print(f"  depth-degenerate groups: {len(degenerate)}")
    for d in degenerate:
        print(f"    {d[3:]}")


if __name__ == "__main__":
    main()
