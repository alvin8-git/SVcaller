#!/usr/bin/env python3
"""Map HGVS coding coordinates (c.N) to GRCh38 genomic positions.

The globin site panels must not be built by hand. Alpha- and beta-globin sit on
opposite strands, both genes are three-exon, and a panel keyed on the wrong base
calls nothing while looking like it works. This derives every coordinate from a
RefSeq gene model and refuses to run if it cannot reproduce known answers.

Gene models come from the AnnotSV annotation bundle (BED12/genePred layout):
    ${annotsv_db}/Genes/GRCh38/genes.RefSeq.sorted.bed
    chrom  txStart  txEnd  strand  gene  tx  cdsStart  cdsEnd  exonStarts  exonEnds
All starts are 0-based half-open; this module returns 1-based genomic positions.

    python3 bin/hgvs_map.py --selftest
    python3 bin/hgvs_map.py --gene HBA2 --cpos 377
"""
import argparse
import re
import sys

DEFAULT_BED = "/data/alvin/ref/annotsv/Annotations_Human/Genes/GRCh38/genes.RefSeq.sorted.bed"

# (gene, c.pos, expected 1-based genomic) — every one of these was established
# independently of this code, so agreement is real evidence the mapper is right:
#   HBB c.20   HbS/rs334      1000G chr11:5227002 T>A, AF_AFR=0.100
#   HBB c.79   HbE            1000G chr11:5226943 C>T, EAS+SAS only
#   HBB c.126  CD41-42        1000G chr11:5226762 CAAAG>C (anchor base 5226762)
#   HBB c.52   CD17           1000G chr11:5226970 T>A, EAS only
#   HBA2 c.377 Hb Quong Sze   observed in THAL2 reads, chr16:173548 T>C VAF 0.55
#   HBB c.316-197  IVS-II-654  chr11:5225923 — cross-checked two ways: as an
#                  offset from exon 3, and by counting 654 nt into IVS-2 from
#                  its 5' end (IVS-2 is 850 nt; 850-654+1 = 197). Both agree.
KNOWN = [
    ("HBB", "20", 5227002),
    ("HBB", "79", 5226943),
    ("HBB", "126", 5226766),  # first deleted base; VCF anchors one base 3' on the genome
    ("HBB", "52", 5226970),
    ("HBA2", "377", 173548),
    ("HBB", "316-197", 5225923),
]


class GeneModel:
    def __init__(self, row):
        f = row.rstrip("\n").split("\t")
        if len(f) < 10:
            raise ValueError(f"expected 10 fields, got {len(f)}: {row[:80]}")
        self.chrom = f[0] if f[0].startswith("chr") else "chr" + f[0]
        self.strand = f[3]
        self.gene = f[4]
        self.tx = f[5]
        self.tx_start = int(f[1])           # 0-based
        self.tx_end = int(f[2])             # 1-based inclusive
        self.cds_start = int(f[6])          # 0-based
        self.cds_end = int(f[7])            # 1-based inclusive end of CDS
        starts = [int(x) for x in f[8].rstrip(",").split(",")]
        ends = [int(x) for x in f[9].rstrip(",").split(",")]
        self.exons = list(zip(starts, ends))  # 0-based half-open

    def coding_positions(self):
        """1-based genomic positions of the CDS, in translation order."""
        pos = []
        for s, e in self.exons:                      # exons are ascending on the genome
            lo = max(s, self.cds_start) + 1          # -> 1-based
            hi = min(e, self.cds_end)
            if lo <= hi:
                pos.extend(range(lo, hi + 1))
        if self.strand == "-":
            pos.reverse()
        return pos

    def c_to_g(self, cpos):
        """c.N (1-based coding) -> 1-based genomic position."""
        coding = self.coding_positions()
        if not 1 <= cpos <= len(coding):
            raise ValueError(
                f"{self.gene}: c.{cpos} outside CDS (length {len(coding)})")
        return coding[cpos - 1]

    def cds_length(self):
        return len(self.coding_positions())

    def _step(self, gpos, n):
        """Move n bases in transcript direction from a genomic position."""
        return gpos + n if self.strand == "+" else gpos - n

    def _in_exon(self, gpos):
        return any(s < gpos <= e for s, e in self.exons)

    def resolve(self, hgvs):
        """Resolve 'c.' notation to a 1-based genomic position.

        Handles  c.377          coding
                 c.92+1         intronic, offset 3' of the anchor coding base
                 c.316-197      intronic, offset 5' of the anchor coding base
                 c.-28          5' UTR, upstream of the start codon

        Offsets step along the genome from the anchor base. That is exact for
        standard IVS nomenclature (the offset lies in the intron immediately
        flanking the anchor exon) and for UTR positions in the same exon, but it
        would be wrong if the offset crossed another exon. Rather than return a
        plausible wrong answer, that case raises.
        """
        s = hgvs.strip()
        if s.startswith("c."):
            s = s[2:]
        m = re.fullmatch(r'(-?\d+)(?:([+-])(\d+))?', s)
        if not m:
            raise ValueError(f"{self.gene}: cannot parse HGVS position {hgvs!r}")
        base, sign, off = m.group(1), m.group(2), m.group(3)
        base = int(base)

        if base < 0:                       # c.-N : 5' UTR or promoter
            if sign:
                raise ValueError(f"{self.gene}: offsets on UTR positions unsupported ({hgvs})")
            g = self._step(self.c_to_g(1), base)   # base is negative
            # Inside an exon (5' UTR) is fine, and so is upstream of the whole
            # transcript (promoter — nothing is spliced there). Landing inside
            # the transcript but between exons means the 5' UTR is spliced and
            # contiguous stepping would give the wrong base.
            outside_tx = g <= self.tx_start or g > self.tx_end
            if not (self._in_exon(g) or outside_tx):
                raise ValueError(
                    f"{self.gene}: {hgvs} lands in an intron — the 5' UTR is "
                    "spliced here, resolve it manually")
            return g
        if base == 0:
            raise ValueError(f"{self.gene}: c.0 is not a valid position")

        anchor = self.c_to_g(base)
        if not sign:
            return anchor
        g = self._step(anchor, int(off) if sign == "+" else -int(off))
        if self._in_exon(g):
            raise ValueError(
                f"{self.gene}: {hgvs} resolves into an exon, not an intron — "
                "the offset crosses a splice boundary and this mapping would be wrong")
        return g


def load_models(bed_path, genes):
    """Return {gene: GeneModel}. Picks the NM_ transcript; errors on ambiguity."""
    found = {}
    with open(bed_path) as fh:
        for line in fh:
            f = line.split("\t", 6)
            if len(f) < 6 or f[4] not in genes:
                continue
            if not f[5].startswith("NM_"):
                continue
            m = GeneModel(line)
            if m.gene in found and found[m.gene].tx != m.tx:
                raise ValueError(
                    f"{m.gene}: multiple RefSeq transcripts "
                    f"({found[m.gene].tx}, {m.tx}) — pin one explicitly")
            found[m.gene] = m
    missing = set(genes) - set(found)
    if missing:
        raise SystemExit(f"FATAL: gene model(s) not found in {bed_path}: {sorted(missing)}")
    return found


def selftest(bed_path):
    models = load_models(bed_path, {g for g, _, _ in KNOWN})
    ok = True
    for gene, cpos, expect in KNOWN:
        m = models[gene]
        got = m.resolve(cpos)
        flag = "ok  " if got == expect else "FAIL"
        if got != expect:
            ok = False
        print(f"  {flag} {gene} c.{cpos:<9} -> {m.chrom}:{got:<9} expected {expect}")
    for gene, m in sorted(models.items()):
        n = m.cds_length()
        print(f"  info {gene} ({m.tx}, {m.strand}) CDS={n} nt = {n/3:.0f} codons")
        if n % 3:
            print(f"  FAIL {gene} CDS length {n} is not a multiple of 3")
            ok = False
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bed", default=DEFAULT_BED, help="RefSeq gene model BED (AnnotSV bundle)")
    ap.add_argument("--gene")
    ap.add_argument("--cpos", help="HGVS coding position, e.g. 377, 92+1, 316-197, -28")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()

    if a.selftest:
        print("hgvs_map selftest — known coordinates established independently of this code:")
        sys.exit(0 if selftest(a.bed) else 1)

    if not (a.gene and a.cpos):
        ap.error("need --gene and --cpos, or --selftest")
    m = load_models(a.bed, {a.gene})[a.gene]
    print(f"{m.chrom}:{m.resolve(a.cpos)}\t{a.gene}\tc.{a.cpos}\t{m.tx}\t{m.strand}")


if __name__ == "__main__":
    main()
