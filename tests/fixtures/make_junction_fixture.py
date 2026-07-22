#!/usr/bin/env python3
"""Build the synthetic alpha-globin junction BAM used by channel-3 tests.

WHY SYNTHETIC RATHER THAN A SLICE OF THAL1
------------------------------------------
A slice of a real sample carries a real --SEA junction, which is more realistic,
but its exact breakpoint is unknown — so a test could only assert "a junction was
found somewhere", which is barely a test. Here the breakpoint is a constant, so
a caller can be asserted to find THE RIGHT ONE, to the base. Keep THAL1 for the
integration check; this is the unit fixture.

WHAT IT CONTAINS, AND WHY EACH PART EARNS ITS PLACE
---------------------------------------------------
A heterozygous ~20 kb deletion at chr16:165001-185000 (approximating --SEA,
whose real extent we measured in THAL1 as roughly 164000-186000).

  split reads, left      soft-clipped AT the breakpoint; the clipped bases are
                         the real reference sequence from the far side, so a
                         caller that realigns the clip lands on 185001 exactly
  split reads, right     mirror image, clipped bases from 164951-165000
  discordant pairs       ~20 kb apparent insert, flagged not-proper
  REFERENCE-SPANNING     reads crossing 165000 with no clip at all. THIS IS THE
                         POINT: the deletion is HET, so the intact allele is
                         present too. A fixture of pure junction reads would let
                         a caller that ignores allele balance call a HOMOZYGOUS
                         deletion and still pass.
  distal normal pairs    away from the locus, so "finds a junction everywhere"
                         fails

Soft-clipped sequence is taken from the real GRCh38 reference at generation
time, so the clips genuinely match the far side. The reference is needed to
BUILD the fixture, never to use it — the BAM is self-contained.

    python3 tests/fixtures/make_junction_fixture.py --ref /path/hg38.fa
"""
import argparse
import os
import subprocess
import sys

CHROM = "chr16"
CHROM_LEN = 90338345          # GRCh38 chr16, as it appears in the real BAM headers
DEL_START = 165001            # first deleted base
DEL_END = 185000              # last deleted base
# on the deleted allele, base DEL_START-1 is joined to DEL_END+1
LEFT_BP = DEL_START - 1       # 165000
RIGHT_BP = DEL_END + 1        # 185001
READ_LEN = 150
QUAL = "I" * READ_LEN         # Q40 flat; nothing here tests quality handling


def fetch(ref, start, end):
    """1-based inclusive reference slice, uppercased."""
    out = subprocess.run(["samtools", "faidx", ref, f"{CHROM}:{start}-{end}"],
                         capture_output=True, text=True, check=True).stdout
    seq = "".join(l.strip() for l in out.splitlines()[1:]).upper()
    assert len(seq) == end - start + 1, f"short fetch at {start}-{end}"
    return seq


def rc(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", default="/data/alvin/ref/GRCh38/hg38.fa")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                  "alpha_junction.bam"))
    a = ap.parse_args()
    if not os.path.exists(a.ref):
        sys.exit(f"FATAL: reference not found: {a.ref}\n"
                 "The fixture is committed; you only need the reference to REBUILD it.")

    recs = []

    def add(name, flag, pos, cigar, seq, mate_pos=0, tlen=0, mapq=60):
        assert len(seq) == len(QUAL), f"{name}: seq length {len(seq)} != {len(QUAL)}"
        recs.append((name, flag, pos, mapq, cigar, "=" if mate_pos else "*",
                     mate_pos, tlen, seq, QUAL))

    # ---- split reads spanning the junction, left side -----------------------
    # aligned block ends exactly at LEFT_BP; the soft-clipped tail is real
    # sequence from RIGHT_BP onward, so realignment resolves to 185001.
    for i, m in enumerate((110, 100, 95, 90, 85, 80)):
        s = READ_LEN - m                       # clipped length
        aligned = fetch(a.ref, LEFT_BP - m + 1, LEFT_BP)
        clipped = fetch(a.ref, RIGHT_BP, RIGHT_BP + s - 1)
        add(f"split_left_{i}", 99, LEFT_BP - m + 1, f"{m}M{s}S", aligned + clipped,
            mate_pos=RIGHT_BP + 200, tlen=0)

    # ---- split reads spanning the junction, right side ----------------------
    for i, m in enumerate((110, 100, 95, 90, 85, 80)):
        s = READ_LEN - m
        clipped = fetch(a.ref, LEFT_BP - s + 1, LEFT_BP)
        aligned = fetch(a.ref, RIGHT_BP, RIGHT_BP + m - 1)
        add(f"split_right_{i}", 147, RIGHT_BP, f"{s}S{m}M", clipped + aligned,
            mate_pos=LEFT_BP - 200, tlen=0)

    # ---- discordant pairs: ~20 kb apparent insert, NOT flagged proper -------
    for i in range(6):
        p1 = LEFT_BP - 400 - i * 30
        p2 = RIGHT_BP + 100 + i * 30
        span = (p2 + READ_LEN - 1) - p1 + 1
        add(f"disc_{i}", 97, p1, f"{READ_LEN}M", fetch(a.ref, p1, p1 + READ_LEN - 1),
            mate_pos=p2, tlen=span)
        add(f"disc_{i}", 145, p2, f"{READ_LEN}M", fetch(a.ref, p2, p2 + READ_LEN - 1),
            mate_pos=p1, tlen=-span)

    # ---- REFERENCE-SPANNING reads: the intact allele of a HET deletion ------
    # No clip, straight across LEFT_BP. Roughly balanced with the split reads so
    # a caller must conclude het, not hom.
    for i in range(10):
        p = LEFT_BP - 120 + i * 12
        add(f"refspan_{i}", 99, p, f"{READ_LEN}M", fetch(a.ref, p, p + READ_LEN - 1),
            mate_pos=p + 300, tlen=450)

    # ---- distal normal pairs: a caller must NOT fire here -------------------
    for i, p in enumerate((155000, 158000, 190000, 193000)):
        add(f"normal_{i}", 99, p, f"{READ_LEN}M", fetch(a.ref, p, p + READ_LEN - 1),
            mate_pos=p + 300, tlen=450)
        add(f"normal_{i}", 147, p + 300, f"{READ_LEN}M",
            fetch(a.ref, p + 300, p + 300 + READ_LEN - 1), mate_pos=p, tlen=-450)

    recs.sort(key=lambda r: r[2])

    sam = a.out.replace(".bam", ".sam")
    with open(sam, "w") as fh:
        fh.write("@HD\tVN:1.6\tSO:coordinate\n")
        fh.write(f"@SQ\tSN:{CHROM}\tLN:{CHROM_LEN}\n")
        fh.write("@RG\tID:synthetic\tSM:JUNCTION_FIXTURE\tPL:ILLUMINA\tLB:synthetic\n")
        fh.write(f"@CO\tSYNTHETIC FIXTURE - not real data. Heterozygous deletion "
                 f"{CHROM}:{DEL_START}-{DEL_END}; junction joins {LEFT_BP} to {RIGHT_BP}.\n")
        fh.write("@CO\tGenerated by tests/fixtures/make_junction_fixture.py. "
                 "Regenerate rather than hand-editing.\n")
        for name, flag, pos, mapq, cigar, rnext, pnext, tlen, seq, qual in recs:
            fh.write(f"{name}\t{flag}\t{CHROM}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t"
                     f"{pnext}\t{tlen}\t{seq}\t{qual}\tRG:Z:synthetic\n")

    subprocess.run(["samtools", "view", "-b", "-o", a.out, sam], check=True)
    subprocess.run(["samtools", "index", a.out], check=True)
    os.remove(sam)

    n = subprocess.run(["samtools", "view", "-c", a.out],
                       capture_output=True, text=True, check=True).stdout.strip()
    print(f"{a.out}  {n} reads  {os.path.getsize(a.out):,} bytes")
    print(f"  deletion   {CHROM}:{DEL_START}-{DEL_END}  ({DEL_END-DEL_START+1:,} bp)")
    print(f"  junction   {LEFT_BP} -> {RIGHT_BP}")
    print(f"  split L/R  6 / 6      discordant pairs 6")
    print(f"  ref-spanning 10       distal normal pairs 4")


if __name__ == "__main__":
    main()
