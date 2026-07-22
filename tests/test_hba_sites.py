"""Deterministic unit tests for bin/hba_sites.py (alpha-globin channel 4).

Pure and fast: no containers, no network. The only tests that touch a BAM are
the two at the bottom, which run a *targeted* single-base pileup against the
THAL validation BAMs and skip when those are absent.

The central case here has no real sample and never will in this project: a
compound heterozygote (a deletion allele plus a point mutation on the surviving
gene). Neither THAL1 nor THAL2 is one. The synthetic fixture below is therefore
the only validation that path will ever get, so it is written to fail if the
copy-number correction is removed rather than to restate it.
"""
import os
import shutil
import subprocess
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(REPO, "bin")
sys.path.insert(0, BIN)

import hba_sites as hs  # noqa: E402

PANEL = os.path.join(REPO, "assets", "hba_pathogenic_sites.tsv")
REF_FASTA = "/data/alvin/ref/GRCh38/hg38.fa"
THAL1_BAM = "/data/alvin/ref/THAL/THAL1_30X.bwa.sortdup.bqsr.bam"
THAL2_BAM = "/data/alvin/ref/THAL/THAL2_30X.bwa.sortdup.bqsr.bam"


# --------------------------------------------------------------------------- #
# mpileup base-string decoding
# --------------------------------------------------------------------------- #
def test_read_start_mapq_char_is_not_a_base():
    """`^A` is a read start whose second char is a MAPQ, not an A."""
    s = "^A.^].,"
    counts = hs.parse_mpileup_bases(s, "T")
    assert counts["A"] == 0, "the MAPQ char of ^A was miscounted as an A base"
    assert s.count("A") == 1, "fixture must actually contain a naive-trap 'A'"
    assert counts["T"] == 3          # . . ,
    assert counts["total"] == 3


def test_insertion_sequence_is_not_counted_as_bases():
    """`+2AG` carries literal sequence that is NOT a call at this position."""
    s = ".,+2AG.,"
    counts = hs.parse_mpileup_bases(s, "C")
    assert counts["A"] == 0 and counts["G"] == 0
    assert s.count("A") == 1 and s.count("G") == 1, "fixture lost its naive trap"
    assert counts["C"] == 4          # . , . ,
    assert counts["ins"] == 1
    assert counts["total"] == 4


def test_deletion_run_is_consumed_exactly():
    """`-3ACT` consumes exactly three sequence characters, no more, no less."""
    counts = hs.parse_mpileup_bases("A-3ACTG", "C")
    # first A is a real mismatch, ACT belongs to the deletion, G is a real call
    assert counts["A"] == 1
    assert counts["G"] == 1
    assert counts["C"] == 0
    assert counts["del"] == 1
    assert counts["total"] == 2


def test_multi_digit_indel_length():
    """A 12 bp insertion must consume 12 chars, not 1 (naive single-digit bug)."""
    counts = hs.parse_mpileup_bases(".+12AAAAAAAAAAAA.", "G")
    assert counts["A"] == 0
    assert counts["G"] == 2
    assert counts["total"] == 2


def test_deletion_placeholder_and_refskip_are_not_ref_support():
    counts = hs.parse_mpileup_bases(".*#.><", "G")
    assert counts["G"] == 2, "* # > < must not be folded into ref support"
    assert counts["*"] == 2
    assert counts[">"] == 2
    assert counts["total"] == 6


def test_end_of_read_marker_ignored():
    counts = hs.parse_mpileup_bases(".$,$C$", "T")
    assert counts["T"] == 2 and counts["C"] == 1 and counts["total"] == 3


def test_naive_count_disagrees_on_a_realistic_string():
    """One string exercising ^X, $, indels and * at once."""
    s = "^A.,C$.^].-2AGc,*+3CCCA"
    counts = hs.parse_mpileup_bases(s, "T")
    naive_c = s.count("C") + s.count("c")
    assert counts["C"] == 2, "expected exactly the two real C calls"
    assert naive_c == 5, "fixture no longer traps the naive counter"
    assert counts["total"] == 9      # . , C . . c , * A


def test_total_matches_real_thal2_pileup_string():
    """The literal column-5 string measured at chr16:173548 in THAL2."""
    s = ".,CC,.c.Cc.Cc,.c,.CCCC"
    counts = hs.parse_mpileup_bases(s, "T")
    assert counts["total"] == 22, "must equal the mpileup DP field"
    assert counts["C"] == 12 and counts["T"] == 10


# --------------------------------------------------------------------------- #
# Panel + version
# --------------------------------------------------------------------------- #
def test_parse_panel_reads_all_four_alpha_sites():
    sites = hs.parse_panel(PANEL)
    assert len(sites) == 4, "the committed alpha panel has 4 sites"
    assert {s["gene"] for s in sites} == {"HBA1", "HBA2"}
    qs = [s for s in sites if s["allele"] == "Hb Quong Sze"][0]
    assert (qs["chrom"], qs["pos"]) == ("chr16", 173548)
    assert (qs["genomic_ref"], qs["genomic_alt"]) == ("T", "C")
    assert isinstance(qs["pos"], int)


def test_panel_version_matches_frozen_fixture():
    """The contract requires `site_panel_version`, and the frozen example
    validation/examples/SAMPLE.alpha_globin.tsv pins it to @dfd6ccf. That is the
    sha1 of the FILE CONTENT (`sha1sum`), not the git blob hash (0838eb8...).
    If this fails, either the panel changed or the scheme was reinvented."""
    assert hs.panel_version(PANEL) == "hba_pathogenic_sites.tsv@dfd6ccf"


def test_panel_version_in_fixture_agrees():
    fixture = os.path.join(REPO, "validation/examples/SAMPLE.alpha_globin.tsv")
    with open(fixture) as fh:
        header, row = [l.split("\t") for l in fh.read().splitlines()]
    assert dict(zip(header, row))["site_panel_version"] == hs.panel_version(PANEL)


# --------------------------------------------------------------------------- #
# THE CENTRAL TEST -- copy-number-aware zygosity
# --------------------------------------------------------------------------- #
def _naive_zygosity(vaf):
    """What a copy-number-blind caller does: threshold the VAF and stop."""
    return "hom" if vaf >= 0.8 else ("het" if vaf >= 0.2 else "ref")


def test_same_vaf_means_different_zygosity_by_alpha_gene_count():
    """0.97 VAF at HBA2 is 'hom' on 4 genes but HEMIZYGOUS on a --SEA background.

    On --SEA/aa the surviving HBA2 is single-copy, so a variant on it sits near
    100% VAF while being ONE variant allele -- i.e. --SEA/aQS, a compound
    heterozygote (HbH disease), not a homozygote. A caller that reads 0.97 and
    says 'hom' has made the clinically wrong call.
    """
    vaf, depth = 0.97, 60

    z4, b4 = hs.zygosity_call(vaf, depth, 4, "HBA2")
    z2, b2 = hs.zygosity_call(vaf, depth, 2, "HBA2")

    assert z4 == "hom"
    assert z2 == "hemizygous"
    assert z4 != z2, "copy-number correction had no effect -- it was removed"

    # the naive caller collapses both to the same answer; that is the bug
    assert _naive_zygosity(vaf) == "hom"
    assert _naive_zygosity(vaf) == z4 and _naive_zygosity(vaf) != z2

    # the basis must say which rule fired, in both directions
    assert "diploid_gene_n2" in b4 and "alpha_genes=4" in b4
    assert "hemizygous_gene_n1" in b2 and "alpha_genes=2" in b2
    assert "assumes_2gene_haplotype_deletion" in b2, \
        "the -a3.7/-a3.7 caveat must be recorded, not hidden"


def test_half_vaf_on_sea_background_is_flagged_not_called_het():
    """0.50 at a single-copy gene is physically inconsistent; say so.

    A naive caller reports a clean 0/1 here. With one gene copy there is no
    such thing as 0/1, so the call stays hemizygous and the basis carries the
    disagreement for a human to resolve.
    """
    z, b = hs.zygosity_call(0.50, 40, 2, "HBA2")
    assert z == "hemizygous"
    assert "vaf_below_expected_for_1_copy" in b
    assert _naive_zygosity(0.50) == "het"


def test_expected_vaf_arithmetic():
    """One variant copy of a gene present in c copies gives VAF 1/c."""
    assert hs.expected_vaf(4, "HBA2") == pytest.approx(0.5)
    assert hs.expected_vaf(2, "HBA2") == pytest.approx(1.0)
    assert hs.expected_vaf(1, "HBA1") == pytest.approx(1.0)
    assert hs.expected_vaf(3, "HBA2") is None, "3 genes cannot be split per gene"
    assert hs.expected_vaf(0, "HBA2") is None
    assert hs.expected_vaf(None, "HBA2") is None


def test_three_genes_is_ambiguous_but_a_mid_vaf_still_resolves_het():
    """With 3 alpha genes the per-gene split is 2+1 or 1+2 and depth cannot say.
    A mid-band VAF is only reachable with 2 copies, so het survives; a high VAF
    cannot separate hom-of-2 from hemizygous-of-1 and must stay NA."""
    z_mid, b_mid = hs.zygosity_call(0.48, 40, 3, "HBA2")
    assert z_mid == "het" and "vaf_implies_2_copies" in b_mid

    z_hi, b_hi = hs.zygosity_call(0.96, 40, 3, "HBA2")
    assert z_hi == "NA", "must not pick between hom-of-2 and hemizygous-of-1"
    assert "ambiguous_gene_copies" in b_hi


def test_missing_alpha_genes_degrades_to_vaf_only_not_to_diploid():
    """Absent/NA copy number must NOT silently mean 'assume 2'."""
    for missing in (None,):
        z, b = hs.zygosity_call(0.97, 60, missing, "HBA2")
        assert b == "vaf_only"
        assert z == "NA"
        assert z != hs.zygosity_call(0.97, 60, 4, "HBA2")[0], \
            "unknown copy number produced the same answer as a known diploid"


def test_zero_alpha_genes_with_reads_is_inconsistent_not_a_genotype():
    z, b = hs.zygosity_call(0.95, 40, 0, "HBA2")
    assert z == "NA" and "zero_copies_inconsistent" in b


def test_low_depth_het_hom_uncertainty_is_flagged():
    """At DP=12 a true het yields VAF>=0.8 about 2% of the time, so a 'hom'
    call there carries a caveat. At DP=40 it does not."""
    _, b_low = hs.zygosity_call(0.90, 12, 4, "HBA2")
    _, b_ok = hs.zygosity_call(0.90, 40, 4, "HBA2")
    assert "low_depth_het_hom_uncertain" in b_low
    assert "low_depth_het_hom_uncertain" not in b_ok


# --------------------------------------------------------------------------- #
# call = present | absent | no_call
# --------------------------------------------------------------------------- #
def _site(gene="HBA2", pos=173548, ref="T", alt="C"):
    return {"gene": gene, "allele": "Hb Quong Sze", "hgvs_c": "c.377",
            "chrom": "chr16", "pos": pos,
            "genomic_ref": ref, "genomic_alt": alt}


def test_low_depth_is_no_call_never_absent():
    """The failure mode this guards: a site with 4 reads reported as clean.
    THAL1 chr16:173208 really is DP=4."""
    rec = hs.score_site(_site(pos=173208, ref="G", alt="A"), 4, "....", 2)
    assert rec["call"] == "no_call"
    assert rec["call"] != "absent"
    assert rec["zygosity"] == "NA"
    assert rec["zygosity_basis"] == "no_call_low_depth"
    assert rec["depth"] == 4


def test_depth_exactly_at_threshold_is_callable():
    rec = hs.score_site(_site(), hs.DEFAULT_MIN_DEPTH, "." * hs.DEFAULT_MIN_DEPTH, 4)
    assert rec["call"] == "absent" and rec["zygosity"] == "ref"


def test_single_alt_read_is_not_presence():
    """1 alt in 40 is sequencing noise / cross-mapping, not an allele."""
    rec = hs.score_site(_site(), 40, "." * 39 + "C", 4)
    assert rec["call"] == "absent"
    assert rec["alt_count"] == 1
    assert rec["vaf"] == pytest.approx(0.025)


def test_raw_vaf_is_always_emitted():
    """The contract forbids the zygosity string being the only evidence."""
    rec = hs.score_site(_site(), 22, ".,CC,.c.Cc.Cc,.c,.CCCC", 4)
    assert rec["vaf"] == pytest.approx(12 / 22.0, abs=1e-6)
    row = hs.format_row("S", rec).split("\t")
    assert row[hs.COLUMNS.index("vaf")] == "0.545"
    # present even when no zygosity could be assigned
    rec_na = hs.score_site(_site(), 22, ".,CC,.c.Cc.Cc,.c,.CCCC", None)
    assert rec_na["zygosity"] == "NA"
    assert hs.format_row("S", rec_na).split("\t")[hs.COLUMNS.index("vaf")] == "0.545"


# --------------------------------------------------------------------------- #
# Output file shape
# --------------------------------------------------------------------------- #
EXPECTED_HEADER = (
    "#sample\tgene\tallele\thgvs_c\tchrom\tpos\tref\talt\tdepth\tref_count\t"
    "alt_count\tvaf\tzygosity\tzygosity_basis\talpha_genes\tcall")


def test_output_header_is_exact():
    """Schema pin. A rename here is an interface change, not a refactor."""
    assert hs.HEADER == EXPECTED_HEADER


def _run_with_fake_pileup(tmp_path, pileups, alpha_genes, sample="S"):
    """Drive main() with a stubbed pileup so the whole CLI path is exercised."""
    out = tmp_path / ("%s.alpha_sites.tsv" % sample)

    def fake(bam, ref, chrom, pos, **kw):
        return pileups.get(pos, (0, "", ""))

    real, hs.mpileup_site = hs.mpileup_site, fake
    try:
        argv = ["--bam", "x.bam", "--panel", PANEL, "--ref", "r.fa",
                "--sample", sample, "--out", str(out)]
        if alpha_genes is not None:
            argv += ["--alpha-genes", alpha_genes]
        assert hs.main(argv) == 0
    finally:
        hs.mpileup_site = real
    return out.read_text().splitlines()


def test_every_panel_row_appears_even_with_no_coverage(tmp_path):
    """A site missing from the file and a site that was clean are different
    statements. Zero coverage everywhere must still produce four rows."""
    lines = _run_with_fake_pileup(tmp_path, {}, "4")
    assert lines[0] == EXPECTED_HEADER
    assert len(lines) == 1 + 4
    panel = hs.parse_panel(PANEL)
    got = [(l.split("\t")[1], int(l.split("\t")[5])) for l in lines[1:]]
    assert got == [(s["gene"], s["pos"]) for s in panel], "row order/identity drifted"
    for l in lines[1:]:
        f = l.split("\t")
        assert f[hs.COLUMNS.index("call")] == "no_call"
        assert f[hs.COLUMNS.index("depth")] == "0"


def test_row_width_matches_header(tmp_path):
    lines = _run_with_fake_pileup(tmp_path, {173548: (22, ".,CC,.c.Cc.Cc,.c,.CCCC", "")}, "4")
    width = len(lines[0].split("\t"))
    for l in lines[1:]:
        assert len(l.split("\t")) == width


def test_alpha_genes_na_is_accepted_and_recorded(tmp_path):
    lines = _run_with_fake_pileup(
        tmp_path, {173548: (22, ".,CC,.c.Cc.Cc,.c,.CCCC", "")}, "NA")
    qs = [l for l in lines if "\t173548\t" in l][0].split("\t")
    assert qs[hs.COLUMNS.index("alpha_genes")] == "NA"
    assert qs[hs.COLUMNS.index("zygosity")] == "NA"
    assert qs[hs.COLUMNS.index("zygosity_basis")] == "vaf_only"
    assert qs[hs.COLUMNS.index("call")] == "present"     # presence is CN-free
    assert qs[hs.COLUMNS.index("vaf")] == "0.545"


def test_alpha_genes_omitted_behaves_like_na(tmp_path):
    lines = _run_with_fake_pileup(
        tmp_path, {173548: (22, ".,CC,.c.Cc.Cc,.c,.CCCC", "")}, None)
    qs = [l for l in lines if "\t173548\t" in l][0].split("\t")
    assert qs[hs.COLUMNS.index("zygosity_basis")] == "vaf_only"


def test_alpha_genes_out_of_range_rejected():
    with pytest.raises(SystemExit):
        hs.build_parser().parse_args(["--bam", "b", "--panel", PANEL, "--ref", "r",
                                      "--sample", "S", "--out", "o",
                                      "--alpha-genes", "7"])


def test_site_genotypes_string_for_the_contract(tmp_path):
    recs = [hs.score_site(_site(), 22, ".,CC,.c.Cc.Cc,.c,.CCCC", 2)]
    assert hs.site_genotypes(recs) == "HBA2:c.377:hemizygous"
    assert hs.site_genotypes([hs.score_site(_site(), 40, "." * 40, 4)]) == "none"


def test_module_import_has_no_side_effects():
    """Importing must not shell out, read a BAM, or write anything."""
    r = subprocess.run([sys.executable, "-c", "import hba_sites"],
                       cwd=BIN, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert r.stdout == "" and r.stderr == ""


# --------------------------------------------------------------------------- #
# Real BAMs -- targeted single-base pileups only.
# NEVER run a BED-wide depth pass over these; they are ~76 GB.
# --------------------------------------------------------------------------- #
needs_thal = pytest.mark.skipif(
    not (shutil.which("samtools") and os.path.exists(REF_FASTA)
         and os.path.exists(THAL1_BAM) and os.path.exists(THAL2_BAM)),
    reason="THAL BAMs / GRCh38 FASTA / samtools not present")


@needs_thal
def test_thal2_is_a_quong_sze_heterozygote():
    """validation/thal_truth_table.tsv: THAL2 alpha = aQSa/aa, 4 genes, het,
    chr16:173548 T>C. Measured here: DP=22, 12 C / 10 T, VAF 0.545."""
    depth, bases, _ = hs.mpileup_site(THAL2_BAM, REF_FASTA, "chr16", 173548)
    rec = hs.score_site(_site(), depth, bases, 4)
    assert rec["call"] == "present"
    assert rec["zygosity"] == "het"
    assert 0.40 <= rec["vaf"] <= 0.70, rec
    assert rec["depth"] >= 15


@needs_thal
def test_quong_sze_is_absent_in_thal1():
    depth, bases, _ = hs.mpileup_site(THAL1_BAM, REF_FASTA, "chr16", 173548)
    rec = hs.score_site(_site(), depth, bases, 2)
    assert rec["call"] == "absent"
    assert rec["alt_count"] == 0
    assert rec["zygosity"] == "ref"


@needs_thal
def test_thal1_hba2_adana_site_is_no_call_not_absent():
    """Real evidence for the no_call path: THAL1 is a --SEA het, so HBA2 is
    single-copy and chr16:173208 falls to DP=4 -- below --min-depth. Reporting
    that as 'absent' would be a fabricated negative on a severe allele."""
    depth, bases, _ = hs.mpileup_site(THAL1_BAM, REF_FASTA, "chr16", 173208)
    assert depth < hs.DEFAULT_MIN_DEPTH, "coverage changed; re-derive this test"
    rec = hs.score_site(_site(pos=173208, ref="G", alt="A"), depth, bases, 2)
    assert rec["call"] == "no_call"
