"""Channel 2 (allele naming) and the alpha-globin contract emitter.

These tests are pure: no containers, no network, no BAMs, no reference. Every
observed-copy dict below is what bin/hba_depth.py's `call` column maps to, so a
change to the depth caller's thresholds cannot silently change the genotype a
signature produces without a test here moving too.
"""
import os
import subprocess
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))

import alpha_globin as ag  # noqa: E402

ALLELE_TSV = os.path.join(REPO, "assets", "hba_deletion_alleles.tsv")
PANEL = os.path.join(REPO, "assets", "hba_pathogenic_sites.tsv")
FIXTURE = os.path.join(REPO, "validation", "examples", "SAMPLE.alpha_globin.tsv")
CONTRACT = os.path.join(REPO, "docs", "contracts", "alpha_globin_contract.md")


@pytest.fixture(scope="module")
def alleles():
    return ag.parse_alleles(ALLELE_TSV)


# Observed copies per segment, as bin/hba_depth.py's `call` column maps them:
#   intact -> 2 · het_loss -> 1 · hom_loss -> 0 · gain -> 3
THAL1 = {"HBZ": 2, "HBA2": 1, "INTER_A2_A1": 1, "HBA1": 1}   # --SEA het, measured
NORMAL = {"HBZ": 2, "HBA2": 2, "INTER_A2_A1": 2, "HBA1": 2}  # THAL2 + 6 GIAB


# --------------------------------------------------------------------------- #
# channel 2 — allele naming
# --------------------------------------------------------------------------- #
def test_thal1_is_named_as_the_group_never_as_sea_alone(alleles):
    """THE headline rule. --SEA and --MED have identical depth signatures, so
    naming --SEA from depth invents precision. Choosing it because the sample
    looks SE Asian is a population inference dressed as a measurement."""
    called, genes, note = ag.name_alleles(THAL1, alleles)
    assert called == "--MED|--SEA/aa", called
    assert genes == "2"
    assert "degenerate" in note
    # the failure this guards: a bare allele name with no group
    assert called != "--SEA/aa"


def test_normal_sample_is_not_ambiguous(alleles):
    """A --SEA on one chromosome and an anti-3.7 triplication on the other
    restore every segment to 2 copies, so they are formally indistinguishable
    from aa/aa. Without the parsimony rule EVERY normal sample would report as
    ambiguous — technically true, operationally useless."""
    called, genes, note = ag.name_alleles(NORMAL, alleles)
    assert called == "none"
    assert genes == "4"
    assert note == ""


def test_a37_is_matched_on_inter_a2_a1_alone(alleles):
    """hba_segments.bed states INTER_A2_A1, not HBA1/HBA2, is the diagnostic
    segment for -a3.7 — the commonest alpha-thal allele worldwide."""
    obs = {"HBZ": 2, "HBA2": 2, "INTER_A2_A1": 1, "HBA1": 2}
    called, genes, _ = ag.name_alleles(obs, alleles)
    assert called == "-a3.7/aa"
    assert genes == "3"


def test_a37_does_not_swallow_the_sea_signature(alleles):
    """REGRESSION. 'h' (hybrid) in the allele table once meant "skip this
    segment", making -a3.7 a near-wildcard that matched THAL1's --SEA signature
    and named the wrong allele on the only real sample in the project."""
    hits = ag.match_genotypes(THAL1, alleles)
    assert "-a3.7" not in {n for h in hits for n in (h[0], h[1])}


def test_a42_and_a37_are_not_confused(alleles):
    """Both lose one gene; they differ in WHICH segment goes."""
    a42 = {"HBZ": 2, "HBA2": 1, "INTER_A2_A1": 2, "HBA1": 2}
    assert ag.name_alleles(a42, alleles)[0] == "-a4.2/aa"


def test_hbz_loss_switches_to_the_other_degenerate_group(alleles):
    """--SEA|--MED vs --FIL|--THAI is decided ENTIRELY by whether HBZ is lost."""
    obs = {"HBZ": 1, "HBA2": 1, "INTER_A2_A1": 1, "HBA1": 1}
    called, genes, _ = ag.name_alleles(obs, alleles)
    assert called == "--FIL|--THAI/aa"
    assert genes == "2"


def test_unmeasured_hbz_widens_the_group_rather_than_guessing(alleles):
    """If HBZ was not scored, --FIL/--THAI cannot be excluded. The caller must
    widen the reported group, not quietly assume HBZ is intact."""
    obs = {"HBA2": 1, "INTER_A2_A1": 1, "HBA1": 1}
    called, _, _ = ag.name_alleles(obs, alleles)
    assert set(called.split("/")[0].split("|")) == {"--SEA", "--MED", "--FIL", "--THAI"}


def test_compound_heterozygote_is_resolved(alleles):
    """--SEA/-a3.7 is the clinically important genotype and no sample has it."""
    obs = {"HBZ": 2, "HBA2": 1, "INTER_A2_A1": 0, "HBA1": 1}
    called, genes, _ = ag.name_alleles(obs, alleles)
    assert called == "--MED|--SEA/-a3.7"
    assert genes == "1"


def test_homozygous_deletion_gives_zero_genes(alleles):
    obs = {"HBZ": 2, "HBA2": 0, "INTER_A2_A1": 0, "HBA1": 0}
    called, genes, _ = ag.name_alleles(obs, alleles)
    assert genes == "0"
    assert called == "--MED|--SEA/--MED|--SEA"


def test_signature_matching_nothing_is_NA_not_the_nearest_allele(alleles):
    """A caller that rounds an unexplained signature onto the closest allele
    reports a genotype it did not measure."""
    obs = {"HBZ": 1, "HBA2": 2, "INTER_A2_A1": 2, "HBA1": 1}
    called, genes, note = ag.name_alleles(obs, alleles)
    assert called == "NA" and genes == "NA"
    assert "match no allele" in note


def test_triplication_is_reported_but_gene_count_is_NA(alleles):
    """ESCALATION, pinned as a test. A triplication carrier has 5 alpha genes;
    the frozen contract declares alpha_genes_called as 0-4. Widening the
    contract unilaterally is forbidden, so the count is NA and the allele is
    carried in deletion_alleles instead — no information is lost."""
    obs = {"HBZ": 2, "HBA2": 3, "INTER_A2_A1": 3, "HBA1": 3}
    called, genes, note = ag.name_alleles(obs, alleles)
    assert called == "anti-3.7/aa"
    assert genes == "NA"
    assert "0-4" in note


def test_inter_z_a_is_never_evidence():
    """do_not_average: INTER_Z_A reads 0.99 'intact' in a sample where a --SEA
    deletion covers half of it. It must not become a vote for intact."""
    rows = [{"segment": "INTER_Z_A", "call": "intact"},
            {"segment": "HBA2", "call": "het_loss"}]
    assert ag.observed_copies(rows) == {"HBA2": 1}


def test_uncalibrated_segment_is_omitted_not_assumed_intact():
    rows = [{"segment": "HBZ", "call": "uncalibrated"},
            {"segment": "HBA1", "call": "intact"}]
    assert ag.observed_copies(rows) == {"HBA1": 2}


def test_hybrid_marker_is_never_coerced_to_a_number(alleles):
    """'h' is a marker, not a measurement."""
    a37 = next(a for a in alleles if a["allele"] == "-a3.7")
    assert a37["d_HBA2"] == "h" and a37["d_HBA1"] == "h"
    assert a37["d_INTER_A2_A1"] == "-1"
    assert ag.has_hybrid(a37, "HBA2") and not ag.has_hybrid(a37, "INTER_A2_A1")


def test_allele_table_without_the_inter_segment_is_rejected(tmp_path):
    """The original table carried only (HBZ, HBA2, HBA1). That cannot express
    -a3.7 — the commonest alpha-thal allele worldwide — because its only
    diagnostic segment is INTER_A2_A1. Loading such a table must fail loudly,
    not silently lose the allele."""
    old = tmp_path / "old.tsv"
    old.write_text("allele\td_HBZ\td_HBA2\td_HBA1\talpha_genes_lost\n"
                   "--SEA\t0\t-1\t-1\t2\n")
    with pytest.raises(ag.AlphaGlobinInputError, match="d_INTER_A2_A1"):
        ag.parse_alleles(str(old))


# --------------------------------------------------------------------------- #
# never collapse a degenerate group
# --------------------------------------------------------------------------- #
def test_group_is_never_collapsed_even_with_a_junction(alleles):
    """The contract permits collapsing only on a junction read or an extent that
    EXCLUDES the alternative. No breakpoint table exists (deliberately: the
    alleles are defined by signature, and approx_size is documentary only), so
    there is nothing to compare an extent against. Collapsing on anything else
    is a population inference."""
    junction = [{"left_bp": "165000", "right_bp": "185001", "size": "20000",
                 "zygosity": "het"}]
    called, _, _ = ag.name_alleles(THAL1, alleles, junction_rows=junction)
    assert called == "--MED|--SEA/aa"
    assert ag._collapse_group(["--SEA", "--MED"], junction) is None


# --------------------------------------------------------------------------- #
# contract emission
# --------------------------------------------------------------------------- #
def _depth_rows(obs, marginal=()):
    return [{"segment": s, "call": {2: "intact", 1: "het_loss", 0: "hom_loss",
                                    3: "gain"}[c],
             "marginal": "True" if s in marginal else "False"}
            for s, c in obs.items()]


def _site_row(**kw):
    row = {"gene": "HBA2", "allele": "Hb Quong Sze", "hgvs_c": "c.377",
           "call": "absent", "zygosity": "NA"}
    row.update(kw)
    return row


def test_contract_columns_match_the_frozen_fixture():
    with open(FIXTURE) as fh:
        assert fh.readline().rstrip("\n").split("\t") == ag.CONTRACT_COLUMNS


def test_emitted_row_validates_against_the_fixture_shape(tmp_path, alleles):
    row, _ = ag.build_row("S", _depth_rows(THAL1), None,
                          [_site_row()], alleles, PANEL)
    out = tmp_path / "S.alpha_globin.tsv"
    ag.write_contract(out, row)
    lines = out.read_text().splitlines()
    assert len(lines) == 2
    assert lines[0].split("\t") == ag.CONTRACT_COLUMNS
    assert len(lines[1].split("\t")) == len(ag.CONTRACT_COLUMNS)


def test_interpretation_complete_is_structurally_false(tmp_path, alleles):
    """SVcaller measures; it must never be able to claim it interpreted."""
    row, _ = ag.build_row("S", _depth_rows(NORMAL), None, [_site_row()],
                          alleles, PANEL)
    assert row["interpretation_complete"] == "false"


def test_beta_globin_is_always_declared_not_screened(alleles):
    """not_screened is load-bearing: a consumer must never be able to gate on
    os.path.exists() and render 'Clear'."""
    row, _ = ag.build_row("S", _depth_rows(NORMAL), None, [_site_row()],
                          alleles, PANEL)
    assert "beta_globin" in row["not_screened"]
    assert "alpha_nondeletional_outside_panel" in row["not_screened"]


def test_a_channel_that_did_not_run_moves_its_tier_to_not_screened(alleles):
    row, _ = ag.build_row("S", _depth_rows(NORMAL), None, None, alleles, PANEL)
    assert ag.TIER_SITES in row["not_screened"].split(",")
    assert ag.TIER_SITES not in row["screened"].split(",")
    assert row["site_genotypes"] == "NA"      # not 'none' — nothing was looked at


def test_site_panel_version_pins_the_panel_that_ran():
    """'no pathogenic site found' means nothing without knowing what was
    interrogated."""
    v = ag.panel_version(PANEL)
    assert v.startswith("hba_pathogenic_sites.tsv@")
    assert len(v.split("@")[1]) == 7
    with open(FIXTURE) as fh:
        fh.readline()
        fixture_version = dict(zip(ag.CONTRACT_COLUMNS,
                                   fh.readline().rstrip("\n").split("\t")))
    assert v == fixture_version["site_panel_version"], (
        "the committed panel no longer hashes to the version in the frozen "
        "example fixture — the panel changed, or the hashing scheme did")


def test_present_site_reaches_site_genotypes(alleles):
    row, _ = ag.build_row("S", _depth_rows(NORMAL), None,
                          [_site_row(call="present", zygosity="het")],
                          alleles, PANEL)
    assert row["site_genotypes"] == "HBA2:c.377:het"
    assert row["genotype"].startswith("aa/aa +HBA2:c.377:het")


def test_a_low_depth_site_does_not_become_a_negative(alleles):
    """no_call must not be laundered into 'none', which reads as "we looked and
    it is absent" — the exact false-negative shape this module exists to remove.
    THAL1 really does hit this: its HBA2 is single-copy, so DP at chr16:173208
    is 4."""
    row, _ = ag.build_row("S", _depth_rows(THAL1), None,
                          [_site_row(call="no_call")], alleles, PANEL)
    assert row["site_genotypes"] == "none"
    assert row["alpha_genes_confidence"] != "high"


def test_empty_channel_file_raises_rather_than_reporting_normal(tmp_path):
    """A present-but-empty channel file is a crashed caller. Rendering it as a
    4-gene normal result is the SMN incident repeated."""
    p = tmp_path / "S.alpha_depth.tsv"
    p.write_text("")
    with pytest.raises(ag.AlphaGlobinInputError):
        ag.read_tsv(str(p), "alpha depth")


def test_absent_channel_is_a_skip_not_a_failure(tmp_path):
    assert ag.read_tsv(str(tmp_path / "NO_FILE"), "alpha depth") is None
    assert ag.read_tsv("NO_DEPTH", "alpha depth") is None


def test_confidence_is_medium_without_a_junction_channel(alleles):
    row, _ = ag.build_row("S", _depth_rows(THAL1), None, [_site_row()],
                          alleles, PANEL)
    assert row["alpha_genes_confidence"] == "medium"
    assert row["deletion_evidence"] == "depth"


def test_confidence_is_high_when_junction_agrees(alleles):
    junction = [{"left_bp": "165000", "right_bp": "185001", "zygosity": "het"}]
    row, _ = ag.build_row("S", _depth_rows(THAL1), junction, [_site_row()],
                          alleles, PANEL)
    assert row["alpha_genes_confidence"] == "high"
    assert row["deletion_evidence"] == "both"


def test_confidence_is_low_when_the_channels_disagree(alleles):
    """Depth says no deletion, the junction caller found one."""
    junction = [{"left_bp": "165000", "right_bp": "185001", "zygosity": "het"}]
    row, _ = ag.build_row("S", _depth_rows(NORMAL), junction, [_site_row()],
                          alleles, PANEL)
    assert row["alpha_genes_confidence"] == "low"
    assert row["deletion_evidence"] == "junction"


def test_a_marginal_segment_downgrades_confidence(alleles):
    junction = [{"left_bp": "165000", "right_bp": "185001", "zygosity": "het"}]
    row, _ = ag.build_row("S", _depth_rows(THAL1, marginal=("HBA2",)),
                          junction, [_site_row()], alleles, PANEL)
    assert row["alpha_genes_confidence"] == "low"


def test_genotype_does_not_assert_phase(alleles):
    """The contract ILLUSTRATES genotype as `--SEA/aQSa`, writing the site
    variant into a haplotype. Short reads do not tell us which chromosome a
    site variant is on, and on a deletion background that placement is the whole
    clinical question — so we keep the two measurements separable."""
    row, _ = ag.build_row("S", _depth_rows(THAL1), None,
                          [_site_row(call="present", zygosity="hemizygous")],
                          alleles, PANEL)
    assert row["genotype"] == "--MED|--SEA/aa +HBA2:c.377:hemizygous"
    assert "aQS" not in row["genotype"]


def test_cli_end_to_end(tmp_path, alleles):
    depth = tmp_path / "S.alpha_depth.tsv"
    depth.write_text(
        "#sample\tsegment\tcall\tmarginal\n"
        + "".join("S\t{}\t{}\tFalse\n".format(
            s, {2: "intact", 1: "het_loss"}[c]) for s, c in THAL1.items()))
    out = tmp_path / "S.alpha_globin.tsv"
    r = subprocess.run(
        [sys.executable, os.path.join(REPO, "bin", "alpha_globin.py"),
         "--sample", "S", "--depth", str(depth), "--panel", PANEL,
         "--alleles", ALLELE_TSV, "--out", str(out)],
        capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    header, row = [l.split("\t") for l in out.read_text().splitlines()]
    d = dict(zip(header, row))
    assert d["deletion_alleles"] == "--MED|--SEA/aa"
    assert d["alpha_genes_called"] == "2"
    assert d["interpretation_complete"] == "false"


def test_cli_refuses_when_no_channel_produced_anything(tmp_path):
    """Nothing to report must not become an empty-but-successful contract row."""
    r = subprocess.run(
        [sys.executable, os.path.join(REPO, "bin", "alpha_globin.py"),
         "--sample", "S", "--out", str(tmp_path / "o.tsv")],
        capture_output=True, text=True)
    assert r.returncode != 0
    assert not (tmp_path / "o.tsv").exists()


@pytest.mark.parametrize("obs", [THAL1, NORMAL,
                                 {"HBZ": 2, "HBA2": 2, "INTER_A2_A1": 1, "HBA1": 2}])
def test_emitted_values_never_interpret(obs, alleles):
    """SVcaller measures; OmniGen owns the clinical narrative. No EMITTED value
    may carry a diagnosis, a risk word, or a reassurance — checked on the row,
    not on the source, because the docstrings legitimately discuss these terms
    in order to forbid them."""
    for sites in (None, [_site_row(call="present", zygosity="het")]):
        row, _ = ag.build_row("S", _depth_rows(obs), None, sites, alleles, PANEL)
        blob = " ".join(str(v) for v in row.values()).lower()
        for token in ("thalass", "carrier", "hbh", "bart", "clear", "normal",
                      "negative", "affected", "risk", "disease"):
            assert token not in blob, "{!r} leaked into the contract row".format(token)


def _run_cli(*args):
    return subprocess.run(
        [sys.executable, os.path.join(REPO, "bin", "alpha_globin.py")] + list(args),
        capture_output=True, text=True)


def test_print_alpha_genes_mode(tmp_path):
    """Channel 4 needs the alpha-gene count to turn a VAF into a zygosity, but
    that count comes from channel 2 (here), not from channel 1's per-segment
    TSV. subworkflows/alpha_globin.nf calls this mode rather than re-deriving
    the count, so the two cannot drift apart."""
    depth = tmp_path / "d.tsv"
    depth.write_text(
        "#sample\tsegment\tcall\tmarginal\n"
        + "".join("S\t{}\t{}\tFalse\n".format(
            s, {2: "intact", 1: "het_loss"}[c]) for s, c in THAL1.items()))
    r = _run_cli("--sample", "S", "--depth", str(depth),
                 "--alleles", ALLELE_TSV, "--print-alpha-genes")
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip() == "2"


def test_print_alpha_genes_degrades_to_NA_not_to_2(tmp_path):
    """A missing depth channel must NOT silently become a 2-copy diploid
    assumption downstream — that is the wrong denominator for every VAF on a
    deletion background."""
    r = _run_cli("--sample", "S", "--depth", "NO_FILE",
                 "--alleles", ALLELE_TSV, "--print-alpha-genes")
    assert r.returncode == 0
    assert r.stdout.strip() == "NA"


def test_out_is_required_in_normal_mode():
    r = _run_cli("--sample", "S", "--alleles", ALLELE_TSV)
    assert r.returncode != 0
    assert "--out" in r.stderr
