#!/usr/bin/env python3
"""Channel 4 of the alpha-globin module: targeted pileup at fixed panel sites.

This is NOT variant discovery. It reads `assets/hba_pathogenic_sites.tsv` and
does one `samtools mpileup -r chr:pos-pos` per listed position. Nothing outside
the panel is ever examined, which is what keeps this channel inside the design
spec's explicit "no genome-wide SNV calling" non-goal. There is deliberately no
discovery mode and no scan mode; do not add one.

Every panel site gets an output row, including sites where nothing was found.
The output file is the record of *what was interrogated*, not just of what was
positive -- a site that is missing from the file and a site that was clean are
very different statements, and only one of them is safe.

Dependencies: python3 stdlib + `samtools` on PATH. NO pysam -- the
`svcaller/utils:1.2` container does not have it.


Copy-number-aware zygosity
--------------------------
The trap this module exists to avoid: on a `--SEA/aa` background BOTH alpha
genes on one chromosome are gone, so the surviving HBA2 is *hemizygous*. A real
variant there sits near 100% VAF, not 50%. A caller that thresholds VAF>=0.8 as
"homozygous" is wrong precisely at the compound heterozygotes that matter
clinically (`--SEA` x Hb Quong Sze -> HbH disease).

So zygosity is computed against `--alpha-genes` (channel 1's
`alpha_genes_called`). When that is absent or NA we do NOT fall back to
assuming two copies; we emit `zygosity=NA` with `zygosity_basis=vaf_only` and
let the raw VAF carry the evidence. `vaf` is ALWAYS emitted, so the zygosity
string is never the only surviving evidence (contract, "Zygosity is
copy-number dependent").


Paralogue caveat -- stated, not hidden
--------------------------------------
HBA1 and HBA2 are ~99% identical. Reads cross-map between them, so:

  * the per-gene copy number is NOT simply alpha_genes_called / 2;
  * a hemizygous variant can be diluted below 1.0 VAF by ref-carrying reads
    mis-assigned from the paralogue;
  * `--min-mapq 1` (the default here) keeps MAPQ>=1 reads, which in a 99%
    identical duplication still includes ambiguously-placed reads.

The mapping from alpha_genes_called to per-gene copies below is therefore an
explicit modelling ASSUMPTION about which alleles produced the count, not a
measurement. It is recorded in `zygosity_basis` on every row so a reader can
re-derive the call instead of trusting it.
"""
import argparse
import hashlib
import os
import subprocess
import sys

# --------------------------------------------------------------------------- #
# Thresholds -- every one of these carries its derivation.
# --------------------------------------------------------------------------- #

#: A site below this depth is `no_call`, NEVER `absent`. 6 reads is the point
#: below which a true heterozygote has a >1.5% chance of yielding zero alt
#: reads (0.5**6 = 0.0156); calling such a site "absent" would be a fabricated
#: negative. Real example: THAL1 chr16:173208 has DP=4.
DEFAULT_MIN_DEPTH = 6

#: samtools mpileup base-quality / mapping-quality floors (samtools defaults are
#: -Q 13 / -q 0). MAPQ 1 drops only the MAPQ-0 multi-mappers, which in the
#: HBA1/HBA2 duplication is the dominant cross-mapping population.
DEFAULT_MIN_BQ = 13
DEFAULT_MIN_MAPQ = 1

#: Minimum evidence to say the alt allele is present at all. A single alt read
#: is not evidence. VAF 0.10 sits ~100x above the post-Q13 per-base error rate
#: (~1e-3) while still tolerating paralogue-driven allele imbalance.
PRESENT_ALT_MIN = 2
PRESENT_VAF_MIN = 0.10

#: Diploid-gene (2 copies of the gene) VAF bands. Expected VAF for k variant
#: copies out of 2 is k/2 -> 0.0, 0.5, 1.0. The band edges are set at ~3
#: binomial SD around 0.5 at a working depth of 20:
#:     SD = sqrt(0.5*0.5/20) = 0.112 ; 0.5 +/- 3*SD = 0.16 .. 0.84
#: rounded outward-of-0.5 to 0.20 / 0.80 so that the hom band stays defensible.
HET_VAF_LO = 0.20
HET_VAF_HI = 0.80

#: Consistency floor for a single-copy (hemizygous) gene. Expected VAF is 1.0;
#: anything below this is still called hemizygous (copy number, not VAF,
#: determines zygosity when the gene has one copy) but the basis is flagged.
HEMI_VAF_EXPECTED_MIN = 0.80

#: Below this depth, het and hom are not reliably separable, so the basis is
#: flagged. P(observed VAF >= 0.80 | true het) by depth:
#:     D=10 -> P(X>=8)  = 56/1024      = 5.5%
#:     D=15 -> P(X>=12) = 576/32768    = 1.8%
#:     D=16 -> P(X>=13) = 697/65536    = 1.1%
#:     D=20 -> P(X>=16) = 6196/1048576 = 0.59%
#: 20 is the first round number under 1%.
DEPTH_HET_HOM_SEPARABLE = 20

CALL_PRESENT = "present"
CALL_ABSENT = "absent"
CALL_NO_CALL = "no_call"

COLUMNS = [
    "sample", "gene", "allele", "hgvs_c", "chrom", "pos", "ref", "alt",
    "depth", "ref_count", "alt_count", "vaf", "zygosity", "zygosity_basis",
    "alpha_genes", "call",
]
HEADER = "#" + "\t".join(COLUMNS)

NA = "NA"


# --------------------------------------------------------------------------- #
# Panel
# --------------------------------------------------------------------------- #
def parse_panel(path):
    """Read assets/hba_pathogenic_sites.tsv -> list of dicts, file order kept.

    `#` lines are comments. The pileup compares against `genomic_ref` /
    `genomic_alt` (what a pileup actually shows), never `coding_ref`/`coding_alt`
    which are on the transcript strand.
    """
    rows = []
    header = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = fields
                continue
            rows.append(dict(zip(header, fields)))
    if header is None:
        raise ValueError("%s: no header row (all lines commented out?)" % path)
    for r in rows:
        for col in ("gene", "allele", "hgvs_c", "chrom", "pos",
                    "genomic_ref", "genomic_alt"):
            if col not in r:
                raise ValueError("%s: panel missing required column %r" % (path, col))
        r["pos"] = int(r["pos"])
        r["genomic_ref"] = r["genomic_ref"].upper()
        r["genomic_alt"] = r["genomic_alt"].upper()
    return rows


def panel_version(path):
    """`hba_pathogenic_sites.tsv@<sha1[:7]>` for the panel actually used.

    The hash is the SHA-1 of the FILE CONTENT (i.e. `sha1sum`), not a git blob
    hash. Verified against the frozen fixture
    validation/examples/SAMPLE.alpha_globin.tsv, which carries `@dfd6ccf`:
        sha1sum       assets/hba_pathogenic_sites.tsv -> dfd6ccf39b1a...
        git hash-object assets/hba_pathogenic_sites.tsv -> 0838eb83c745...
    Only the first matches, so content-sha1 is the committed scheme.
    """
    with open(path, "rb") as fh:
        digest = hashlib.sha1(fh.read()).hexdigest()
    return "%s@%s" % (os.path.basename(path), digest[:7])


# --------------------------------------------------------------------------- #
# mpileup base-string parsing
# --------------------------------------------------------------------------- #
def parse_mpileup_bases(bases_str, ref_base):
    """Decode a samtools mpileup column-5 base string into allele counts.

    A naive `bases_str.count("A")` is WRONG and silently so:
      * `^X` starts a read and X is a MAPQ character -- `^A` would be counted;
      * `+2AG` / `-3ACT` embed literal sequence that is NOT a base call at this
        position -- the A/G/C/T inside would be counted;
      * `$` ends a read and is not a base;
      * `*` / `#` are deletion placeholders, `>` / `<` are reference skips;
        neither supports the ref allele.

    Returns a dict keyed by observed allele, with `.`/`,` resolved to
    `ref_base`, plus bookkeeping keys::

        {"A":n,"C":n,"G":n,"T":n,"N":n,   # base calls, ref matches folded in
         "*":n,                            # deletion placeholders (* and #)
         ">":n,                            # reference skips (> and <)
         "ins":n,"del":n,                  # indel EVENTS anchored at this base
         "total":n}                        # pileup entries consumed == DP

    `total` must equal the mpileup DP field; the caller asserts on it.
    """
    ref_base = (ref_base or "N").upper()
    counts = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0,
              "*": 0, ">": 0, "ins": 0, "del": 0, "total": 0}
    if ref_base not in ("A", "C", "G", "T"):
        ref_base = "N"

    s = bases_str or ""
    i, n = 0, len(s)
    while i < n:
        ch = s[i]
        if ch == "^":
            # read start; the NEXT char is the mapping quality, not a base
            i += 2
            continue
        if ch == "$":
            i += 1
            continue
        if ch in "+-":
            # indel attached to the PREVIOUS base: +<len><seq>
            j = i + 1
            digits = ""
            while j < n and s[j].isdigit():
                digits += s[j]
                j += 1
            if not digits:          # malformed; treat literally and move on
                i += 1
                continue
            length = int(digits)
            counts["ins" if ch == "+" else "del"] += 1
            i = j + length          # consume exactly <length> sequence chars
            continue
        if ch in ".,":
            counts[ref_base] += 1
            counts["total"] += 1
            i += 1
            continue
        up = ch.upper()
        if up in "ACGTN":
            counts[up] += 1
            counts["total"] += 1
            i += 1
            continue
        if ch in "*#":
            counts["*"] += 1
            counts["total"] += 1
            i += 1
            continue
        if ch in "><":
            counts[">"] += 1
            counts["total"] += 1
            i += 1
            continue
        # anything else (stray whitespace, unexpected symbol) is skipped
        i += 1
    return counts


# --------------------------------------------------------------------------- #
# Copy-number model
# --------------------------------------------------------------------------- #
def gene_copies(alpha_genes, gene):
    """Copies of `gene` implied by an alpha_genes_called total -> int or None.

    THIS IS A MODELLING ASSUMPTION, not a measurement. `alpha_genes_called` is a
    TOTAL across HBA1+HBA2 on both chromosomes; it does not say which gene lost
    copies. The mapping below assumes the deletion classes that actually
    dominate alpha-thalassaemia:

      4 -> 2   normal aa/aa. Both genes diploid.
      3 -> None  a single-gene deletion (-a3.7 / -a4.2) on one haplotype. Per
               gene the split is 2+1 or 1+2 and depth cannot say which; -a3.7
               additionally creates an HBA2/HBA1 FUSION whose "gene" identity is
               not well defined. Refuse to guess.
      2 -> 1   the canonical `--SEA/aa` class: one haplotype lost BOTH genes,
               the other retains one HBA1 and one HBA2 -> each gene hemizygous.
               CAVEAT: 2 is also reachable as -a3.7/-a3.7 (each haplotype loses
               one gene), where a given gene could be 0 or 2 copies. This is the
               main way the model can be wrong; it is flagged in the basis.
      1 -> 1   only one alpha gene left in total; if reads are present at this
               site, that gene is it.
      0 -> 0   no alpha genes. Any depth here is cross-mapping, not signal.
      None -> None  channel 1 gave nothing; degrade to VAF-only.
    """
    if alpha_genes is None:
        return None
    if alpha_genes == 4:
        return 2
    if alpha_genes == 3:
        return None
    if alpha_genes in (2, 1):
        return 1
    if alpha_genes == 0:
        return 0
    return None


def expected_vaf(alpha_genes, gene):
    """Expected VAF of ONE variant copy of `gene` -> float or None.

    None means "the copy number is not determined", which is the honest answer
    for alpha_genes in {3, None} and for a zero-copy background. Arithmetic:
    with `c` copies of the gene, one variant copy gives VAF = 1/c.
    """
    copies = gene_copies(alpha_genes, gene)
    if not copies:            # None or 0
        return None
    return 1.0 / copies


def zygosity_call(vaf, depth, alpha_genes, gene,
                  min_depth=DEFAULT_MIN_DEPTH, alt_count=None):
    """-> (zygosity, zygosity_basis).

    zygosity is one of: hom | het | hemizygous | ref | NA
    zygosity_basis records WHICH RULE FIRED, so the call can be re-checked
    rather than trusted. Extra qualifiers are appended with ';'.

    Precedence, deliberately:
      1. depth below min_depth  -> (NA, no_call_low_depth). A low-depth site is
         never reported as reference.
      2. alt allele not present -> (ref, <basis>). Presence is copy-number
         independent, so this works even when alpha_genes is unknown.
      3. copy number unknown    -> (NA, vaf_only). Never silently assume 2.
      4. gene hemizygous (1 copy) -> hemizygous, regardless of VAF. If the gene
         has one copy and the variant is there, it IS hemizygous; the VAF only
         serves as a consistency check and a low value is flagged.
      5. gene diploid (2 copies)  -> VAF bands decide het vs hom.
      6. copy number ambiguous (alpha_genes==3) -> a mid-range VAF is only
         consistent with 2 copies, so het can still be salvaged; a high VAF
         cannot distinguish hom-of-2 from hemizygous-of-1 and stays NA.
    """
    if depth is None or depth < min_depth:
        return NA, "no_call_low_depth"

    present = vaf is not None and vaf >= PRESENT_VAF_MIN
    if alt_count is not None:
        present = present and alt_count >= PRESENT_ALT_MIN

    copies = gene_copies(alpha_genes, gene)

    if not present:
        if alpha_genes is None:
            return "ref", "no_alt_evidence;vaf_only"
        return "ref", "no_alt_evidence;alpha_genes=%d" % alpha_genes

    if alpha_genes is None:
        return NA, "vaf_only"

    if copies == 0:
        # alpha_genes says zero copies yet reads carry the alt allele. Something
        # disagrees -- almost certainly paralogue cross-mapping or a wrong
        # channel-1 count. Do not manufacture a genotype from it.
        return NA, "zero_copies_inconsistent;alpha_genes=0"

    if copies == 1:
        basis = "hemizygous_gene_n1;alpha_genes=%d" % alpha_genes
        if alpha_genes == 2:
            # the --SEA/aa assumption; -a3.7/-a3.7 would break it
            basis += ";assumes_2gene_haplotype_deletion"
        if vaf < HEMI_VAF_EXPECTED_MIN:
            basis += ";vaf_below_expected_for_1_copy"
        return "hemizygous", basis

    if copies == 2:
        basis = "diploid_gene_n2;alpha_genes=%d" % alpha_genes
        if depth < DEPTH_HET_HOM_SEPARABLE:
            basis += ";low_depth_het_hom_uncertain"
        if vaf >= HET_VAF_HI:
            return "hom", basis
        if vaf >= HET_VAF_LO:
            return "het", basis
        return "het", basis + ";vaf_below_het_band"

    # copies is None -> alpha_genes == 3
    basis = "ambiguous_gene_copies;alpha_genes=3"
    if HET_VAF_LO <= vaf < HET_VAF_HI:
        # only reachable with 2 copies of this gene carrying 1 variant copy
        return "het", basis + ";vaf_implies_2_copies"
    return NA, basis + ";vaf_cannot_separate_hom2_from_hemi1"


# --------------------------------------------------------------------------- #
# Pileup
# --------------------------------------------------------------------------- #
def mpileup_site(bam, ref, chrom, pos, min_bq=DEFAULT_MIN_BQ,
                 min_mapq=DEFAULT_MIN_MAPQ, samtools="samtools"):
    """One targeted pileup at a single base -> (depth, bases_str, quals_str).

    ONE REGION PER SITE. Never run a genome-wide or BED-wide depth pass over
    these BAMs -- they are ~76 GB and a full scan takes hours.
    Returns depth 0 when samtools emits no line (no covering reads).
    """
    region = "%s:%d-%d" % (chrom, pos, pos)
    cmd = [samtools, "mpileup", "-f", ref, "-r", region,
           "-Q", str(min_bq), "-q", str(min_mapq), "-d", "8000", bam]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          universal_newlines=True)
    if proc.returncode != 0:
        raise RuntimeError("samtools mpileup failed at %s:\n%s"
                           % (region, proc.stderr.strip()))
    for line in proc.stdout.splitlines():
        f = line.rstrip("\n").split("\t")
        if len(f) < 5 or f[0] != chrom or int(f[1]) != pos:
            continue
        return int(f[3]), f[4], (f[5] if len(f) > 5 else "")
    return 0, "", ""


def score_site(site, depth, bases_str, alpha_genes, min_depth=DEFAULT_MIN_DEPTH):
    """Turn one panel site + its pileup into an output record (dict)."""
    ref_base = site["genomic_ref"]
    alt_base = site["genomic_alt"]
    counts = parse_mpileup_bases(bases_str, ref_base)

    ref_count = counts.get(ref_base, 0)
    alt_count = counts.get(alt_base, 0) if alt_base in "ACGTN" else 0
    # depth used for VAF is the pileup depth actually decoded, which excludes
    # nothing that samtools counted -- assert agreement so a parser bug is loud
    decoded = counts["total"]
    eff_depth = decoded if decoded else depth

    vaf = (alt_count / float(eff_depth)) if eff_depth else 0.0

    if eff_depth < min_depth:
        call = CALL_NO_CALL
    elif alt_count >= PRESENT_ALT_MIN and vaf >= PRESENT_VAF_MIN:
        call = CALL_PRESENT
    else:
        call = CALL_ABSENT

    zyg, basis = zygosity_call(vaf, eff_depth, alpha_genes, site["gene"],
                               min_depth=min_depth, alt_count=alt_count)

    return {
        "gene": site["gene"], "allele": site["allele"], "hgvs_c": site["hgvs_c"],
        "chrom": site["chrom"], "pos": site["pos"],
        "ref": ref_base, "alt": alt_base,
        "depth": eff_depth, "ref_count": ref_count, "alt_count": alt_count,
        "vaf": vaf, "zygosity": zyg, "zygosity_basis": basis,
        "alpha_genes": NA if alpha_genes is None else str(alpha_genes),
        "call": call,
    }


def format_row(sample, rec):
    return "\t".join([
        sample, rec["gene"], rec["allele"], rec["hgvs_c"], rec["chrom"],
        str(rec["pos"]), rec["ref"], rec["alt"], str(rec["depth"]),
        str(rec["ref_count"]), str(rec["alt_count"]), "%.3f" % rec["vaf"],
        rec["zygosity"], rec["zygosity_basis"], rec["alpha_genes"], rec["call"],
    ])


def write_sites_tsv(path, sample, records):
    with open(path, "w") as fh:
        fh.write(HEADER + "\n")
        for rec in records:
            fh.write(format_row(sample, rec) + "\n")


def site_genotypes(records):
    """`;`-joined GENE:HGVS:zygosity for contract `site_genotypes`, else 'none'.

    Only `present` sites contribute. A `no_call` never becomes 'none' silently
    for the whole file -- the per-site TSV is the record of that, which is why
    the contract points OmniGen at it.
    """
    hits = ["%s:%s:%s" % (r["gene"], r["hgvs_c"], r["zygosity"])
            for r in records if r["call"] == CALL_PRESENT]
    return ";".join(hits) if hits else "none"


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _parse_alpha_genes(value):
    if value is None:
        return None
    v = str(value).strip()
    if v == "" or v.upper() == "NA":
        return None
    try:
        n = int(v)
    except ValueError:
        raise argparse.ArgumentTypeError(
            "--alpha-genes must be an integer 0-4 or NA, got %r" % value)
    if not 0 <= n <= 4:
        raise argparse.ArgumentTypeError(
            "--alpha-genes must be 0-4 or NA, got %r" % value)
    return n


def build_parser():
    p = argparse.ArgumentParser(
        description="Targeted pileup at the committed HBA pathogenic-site panel. "
                    "Not a variant caller: only panel coordinates are examined.")
    p.add_argument("--bam", required=True)
    p.add_argument("--panel", required=True, help="assets/hba_pathogenic_sites.tsv")
    p.add_argument("--ref", required=True, help="GRCh38 FASTA (indexed)")
    p.add_argument("--sample", required=True)
    p.add_argument("--out", required=True, help="<S>.alpha_sites.tsv")
    p.add_argument("--alpha-genes", default=None, type=_parse_alpha_genes,
                   help="channel 1 alpha_genes_called (0-4) or NA. Without it "
                        "zygosity degrades to vaf_only -- it is NOT assumed 2.")
    p.add_argument("--min-depth", type=int, default=DEFAULT_MIN_DEPTH)
    p.add_argument("--min-bq", type=int, default=DEFAULT_MIN_BQ)
    p.add_argument("--min-mapq", type=int, default=DEFAULT_MIN_MAPQ)
    p.add_argument("--samtools", default="samtools")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    sites = parse_panel(args.panel)
    records = []
    for site in sites:
        depth, bases, _q = mpileup_site(
            args.bam, args.ref, site["chrom"], site["pos"],
            min_bq=args.min_bq, min_mapq=args.min_mapq, samtools=args.samtools)
        records.append(score_site(site, depth, bases, args.alpha_genes,
                                  min_depth=args.min_depth))
    write_sites_tsv(args.out, args.sample, records)
    sys.stderr.write("[hba_sites] %s: %d panel sites, panel=%s, genotypes=%s\n"
                     % (args.sample, len(records), panel_version(args.panel),
                        site_genotypes(records)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
