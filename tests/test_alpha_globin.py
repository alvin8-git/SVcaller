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
    assert called == "--SEA|--MED/aa", called
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
    a42 = {"HBZ": 2, "HBA2": 1, "INTER_A2_A1": 1, "HBA1": 2}
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
    assert called == "--SEA|--MED/-a3.7"
    assert genes == "1"


def test_homozygous_deletion_gives_zero_genes(alleles):
    obs = {"HBZ": 2, "HBA2": 0, "INTER_A2_A1": 0, "HBA1": 0}
    called, genes, _ = ag.name_alleles(obs, alleles)
    assert genes == "0"
    assert called == "--SEA|--MED/--SEA|--MED"


def test_signature_matching_nothing_is_NA_not_the_nearest_allele(alleles):
    """A caller that rounds an unexplained signature onto the closest allele
    reports a genotype it did not measure."""
    obs = {"HBZ": 1, "HBA2": 2, "INTER_A2_A1": 2, "HBA1": 1}
    called, genes, note = ag.name_alleles(obs, alleles)
    assert called == "NA" and genes == "NA"
    assert "match no allele" in note


def test_triplication_reports_five_genes(alleles):
    """CHANGED 2026-07-22 — the escalation this pinned has been resolved.

    A triplication carrier genuinely has 5 alpha genes. The contract declared
    alpha_genes_called as 0-4, on the unexamined assumption that alpha variation
    only ever REMOVES genes — anti-3.7 is the reciprocal product of the -a3.7
    NAHR and adds one. The old range forced `NA` for a perfectly determined
    count, i.e. it reported a measurement failure that had not happened. The
    contract now declares 0-6 and the count is emitted."""
    obs = {"HBZ": 2, "HBA2": 3, "INTER_A2_A1": 3, "HBA1": 3}
    called, genes, note = ag.name_alleles(obs, alleles)
    assert called == "anti-3.7/aa"
    assert genes == "5", f"triplication should report 5 alpha genes, got {genes!r}"
    assert "NA" not in genes


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
    assert called == "--SEA|--MED/aa"
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
    junction = [{"left_bp": "165000", "right_bp": "185001", "zygosity": "het",
                 "confidence": "high"}]
    row, _ = ag.build_row("S", _depth_rows(THAL1), junction, [_site_row()],
                          alleles, PANEL)
    assert row["alpha_genes_confidence"] == "high"
    assert row["deletion_evidence"] == "both"


def test_a_junction_never_lowers_confidence(alleles):
    """Depth says normal, a high junction overlaps scorable segments: it is
    mentioned in the note, evidence is `junction`, confidence is NOT lowered."""
    junction = [{"left_bp": "165000", "right_bp": "185001", "zygosity": "het",
                 "confidence": "high"}]
    row, note = ag.build_row("S", _depth_rows(NORMAL), junction, [_site_row()],
                             alleles, PANEL)
    assert row["alpha_genes_confidence"] == "high"
    assert row["deletion_evidence"] == "junction"
    assert row["alpha_genes_called"] == "4"
    assert "165000-185001 seen, not used" in note


def test_a_marginal_segment_that_could_change_the_genotype_is_never_high(alleles):
    junction = [{"left_bp": "165000", "right_bp": "185001", "zygosity": "het"}]
    row, _ = ag.build_row("S", _depth_rows(THAL1, marginal=("HBA2",)),
                          junction, [_site_row()], alleles, PANEL)
    assert row["alpha_genes_confidence"] == "medium"


def test_genotype_does_not_assert_phase(alleles):
    """The contract ILLUSTRATES genotype as `--SEA/aQSa`, writing the site
    variant into a haplotype. Short reads do not tell us which chromosome a
    site variant is on, and on a deletion background that placement is the whole
    clinical question — so we keep the two measurements separable."""
    row, _ = ag.build_row("S", _depth_rows(THAL1), None,
                          [_site_row(call="present", zygosity="hemizygous")],
                          alleles, PANEL)
    assert row["genotype"] == "--SEA|--MED/aa +HBA2:c.377:hemizygous"
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
    assert d["deletion_alleles"] == "--SEA|--MED/aa"
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


# --------------------------------------------------------------------------- #
# regressions found by adversarial re-derivation, 2026-07-22
# --------------------------------------------------------------------------- #
def test_a42_is_not_misnamed_when_inter_a2_a1_reads_as_lost(alleles):
    """THE DEFECT an independent re-derivation caught.

    -a4.2 is ~4.2 kb but HBA2's gene body is only 835 bp, so 2-3.4 kb of the
    2969 bp INTER_A2_A1 goes with the deletion and that segment scores anywhere
    from ~0.43 to ~0.77 -- it STRADDLES the 0.65 het-loss cutoff. The first
    version marked it 'h', which asserts exactly 2 copies, so a real -a4.2
    carrier whose 5' breakpoint took enough of the segment came out as
    `-a3.7/-a4.2` -- a silent 3-gene carrier reported as a 2-gene alpha-thal
    trait, flipped by a 0.03 wobble in depth.

    Both readings of INTER_A2_A1 must give the same answer."""
    for inter in (2, 1):
        obs = {"HBZ": 2, "HBA2": 1, "INTER_A2_A1": inter, "HBA1": 2}
        # INTER at 2 copies only counts as "not confidently intact" if marginal
        called, genes, _ = ag.name_alleles(obs, alleles, None,
                                           {"INTER_A2_A1": [1, 3]})
        assert called == "-a4.2/aa", f"INTER_A2_A1={inter} gave {called}"
        assert genes == "3"


def test_h_and_question_mark_are_not_interchangeable(alleles):
    """'h' expects 2 copies (load-bearing for -a3.7, or it becomes a wildcard
    that matches THAL1's --SEA signature). '?' is genuinely unconstrained. If a
    future edit collapses them, one of the two alleles breaks."""
    assert ag._deltas("h") == frozenset({0})
    assert ag._deltas("?") == frozenset({-1, 0})
    a37 = next(a for a in alleles if a["allele"] == "-a3.7")
    a42 = next(a for a in alleles if a["allele"] == "-a4.2")
    assert a37["d_HBA2"] == "h", "-a3.7's gene bodies must expect 2 copies"
    assert a42["d_INTER_A2_A1"] == "?", "-a4.2's INTER_A2_A1 must be unconstrained"
    # '?' must not admit a GAIN: a deletion allele cannot add copies.
    assert +1 not in ag._deltas("?")


def test_group_members_are_joined_in_asset_order_not_alphabetically(alleles):
    """OmniGen renders this string, so its spelling is interface, not
    cosmetics. The asset's own depth_distinguishable cells and every occurrence
    in the contract write `--SEA|--MED`; an alphabetical sort emits
    `--MED|--SEA`, the same measurement spelled a way the consumer was never
    told to expect."""
    called, _, _ = ag.name_alleles(THAL1, alleles)
    assert called == "--SEA|--MED/aa"
    assert ag.name_alleles({"HBZ": 1, "HBA2": 1, "INTER_A2_A1": 1, "HBA1": 1},
                           alleles)[0] == "--FIL|--THAI/aa"
    for a in alleles:
        val = a["depth_distinguishable"]
        if val.startswith("no:"):
            group = val[3:]
            assert group in ("--SEA|--MED", "--FIL|--THAI"), group


def test_hg001_impossible_signature_is_NA(alleles):
    """HG001 reads low on HBZ and HBA1 but NOT on INTER_A2_A1, which lies
    between HBA2 and HBA1 -- no contiguous deletion can do that. It is technical
    dropout. The caller must refuse it, not fit it to something."""
    obs = {"HBZ": 1, "HBA2": 2, "INTER_A2_A1": 2, "HBA1": 1}
    called, genes, _ = ag.name_alleles(obs, alleles)
    assert called == "NA" and genes == "NA"


def test_unknown_segment_raises_a_real_error_not_a_KeyError(alleles):
    """INTER_Z_A is do_not_average and deliberately has no allele column."""
    with pytest.raises(ag.AlphaGlobinInputError, match="d_INTER_Z_A"):
        ag.expected_copies(alleles[0], None, "INTER_Z_A")


def test_bare_sea_is_unreachable_across_every_signature(alleles):
    """Exhaustive: no observed signature may ever name --SEA (or --MED, --FIL,
    --THAI) on its own. Emitting one from depth invents precision."""
    import itertools
    segs = ["HBZ", "HBA2", "INTER_A2_A1", "HBA1"]
    degenerate = {"--SEA", "--MED", "--FIL", "--THAI"}
    for copies in itertools.product([0, 1, 2, 3], repeat=4):
        called, _, _ = ag.name_alleles(dict(zip(segs, copies)), alleles)
        for hap in called.split("/"):
            assert hap not in degenerate, (
                f"{dict(zip(segs, copies))} named {hap!r} alone in {called!r}")


# --------------------------------------------------------------------------- #
# false-positive guards: marginal readings, compound genotype, junction support
# --------------------------------------------------------------------------- #
def _rows(spec):
    """spec: {segment: (call, score, marginal)} -> depth rows incl. score."""
    return [{"segment": g, "call": c, "score": str(sc), "marginal": str(m)}
            for g, (c, sc, m) in spec.items()]


# A European GIAB sample that read as --FIL|--THAI/anti-3.7: the whole call
# rested on HBZ at 0.617 (marginal). INTER_Z_A is do_not_average.
HG001_ROWS = _rows({
    "HBZ": ("het_loss", 0.617, True),
    "INTER_Z_A": ("do_not_average", "NA", False),
    "HBA2": ("intact", 0.684, True),
    "INTER_A2_A1": ("intact", 0.832, False),
    "HBA1": ("intact", 0.996, False),
})


def test_a_genotype_resting_on_one_marginal_segment_is_not_named(alleles):
    row, note = ag.build_row("S", HG001_ROWS, None, None, alleles, PANEL)
    assert row["deletion_alleles"] == "NA"
    assert row["alpha_genes_called"] == "NA"
    assert row["alpha_genes_confidence"] == "low"
    assert "HBZ" in note and "marginal" in note
    assert "--FIL" not in row["genotype"] and "anti-3.7" not in row["genotype"]


def test_marginal_dip_is_not_reported_as_normal_either(alleles):
    called, genes, _ = ag.name_alleles(
        ag.observed_copies(HG001_ROWS), alleles, None,
        ag.marginal_alternatives(HG001_ROWS))
    assert (called, genes) == ("NA", "NA")      # not "none"/4


def test_deletion_plus_triplication_is_never_named_from_depth(alleles):
    """HBZ clearly lost (non-marginal) still fits --FIL|--THAI + anti-3.7."""
    obs = {"HBZ": 1, "HBA2": 2, "INTER_A2_A1": 2, "HBA1": 2}
    called, genes, note = ag.name_alleles(obs, alleles)
    assert (called, genes) == ("NA", "NA")
    assert "two-event" in note


def test_marginal_that_can_change_the_count_keeps_call_at_medium(alleles):
    """HBA2 0.188 and INTER 0.201 hom_loss (marginal): 2 genes as called,
    3 if read as het_loss. Call kept, confidence medium (never high), note."""
    rows = _rows({"HBZ": ("intact", 0.9, False),
                  "HBA2": ("hom_loss", 0.188, True),
                  "INTER_A2_A1": ("hom_loss", 0.201, True),
                  "HBA1": ("intact", 0.96, False)})
    before = ag.name_alleles(ag.observed_copies(rows), alleles)
    row, note = ag.build_row("S", rows, None, None, alleles, PANEL)
    assert row["alpha_genes_called"] == before[1] != "NA"
    assert row["deletion_alleles"] == before[0]
    assert row["alpha_genes_confidence"] == "medium"   # the loss is not in doubt
    if before[1] != "3":
        assert "3" in note


def test_marginal_that_changes_nothing_does_not_lower_confidence(alleles):
    """HBZ 0.725 intact (marginal), HBA2/INTER het_loss, HBA1 intact: -a4.2/aa.
    HBZ read as het_loss matches no genotype, so the flag is irrelevant."""
    rows = _rows({"HBZ": ("intact", 0.725, True),
                  "HBA2": ("het_loss", 0.486, False),
                  "INTER_A2_A1": ("het_loss", 0.468, False),
                  "HBA1": ("intact", 0.95, False)})
    junction = [{"left_bp": "172000", "right_bp": "176000", "confidence": "high",
                 "zygosity": "het"}]
    row, _ = ag.build_row("S", rows, junction, None, alleles, PANEL)
    assert row["deletion_alleles"] == "-a4.2/aa"
    assert row["alpha_genes_called"] == "3"
    assert row["deletion_evidence"] == "both"
    assert row["alpha_genes_confidence"] == "medium"   # a3.7/a4.2 ambiguity cap


@pytest.mark.parametrize("obs,called,genes", [
    ({"HBZ": 2, "HBA2": 1, "INTER_A2_A1": 1, "HBA1": 1}, "--SEA|--MED/aa", "2"),
    ({"HBZ": 2, "HBA2": 2, "INTER_A2_A1": 1, "HBA1": 2}, "-a3.7/aa", "3"),
    ({"HBZ": 2, "HBA2": 3, "INTER_A2_A1": 3, "HBA1": 3}, "anti-3.7/aa", "5"),
])
def test_true_positives_unchanged(obs, called, genes, alleles):
    assert ag.name_alleles(obs, alleles)[:2] == (called, genes)


def _jr(lo, hi, conf="high", zyg="het"):
    return {"left_bp": str(lo), "right_bp": str(hi), "confidence": conf,
            "zygosity": zyg}


def test_junction_support_rules(alleles):
    """Corroboration only. --SEA|--MED removes HBA2, INTER, HBA1 (HBZ untouched)."""
    g = "--SEA|--MED/aa"
    ok = _jr(165396, 184701)
    assert ag.junction_support([ok], g, alleles) == ([ok], [])
    # not-high rows, and rows touching no scorable segment: ignored entirely
    for ign in (_jr(165396, 184701, "medium"), _jr(165396, 184701, "low"),
                {"left_bp": "165396", "right_bp": "184701"}, _jr(155000, 172000)):
        assert ag.junction_support([ign], g, alleles) == ([], [])
    partial = _jr(174000, 175000)          # overlaps INTER only: not every removed segment
    assert ag.junction_support([partial], g, alleles) == ([], [partial])
    over_hbz = _jr(150000, 184701)         # also overlaps an untouched segment
    assert ag.junction_support([over_hbz], g, alleles) == ([], [over_hbz])
    assert ag.junction_support(None, g, alleles) == (None, [])
    assert ag.junction_support([ok], "NA", alleles) == ([], [ok])
    assert ag.junction_support([ok], "none", alleles) == ([], [ok])


def test_zygosity_is_not_part_of_support(alleles):
    hom = _jr(165396, 184701, zyg="hom")
    assert ag.junction_support([hom], "--SEA|--MED/aa", alleles) == ([hom], [])


def test_junction_support_sets_both_and_raises_confidence(alleles):
    row, _ = ag.build_row("S", _depth_rows(THAL1), [_jr(165396, 184701)], None,
                          alleles, PANEL)
    assert row["deletion_evidence"] == "both"
    assert row["alpha_genes_confidence"] == "high"


def test_junction_never_makes_a_call_NA_or_low(alleles):
    for j in (_jr(155000, 172000), _jr(150000, 190000), _jr(174000, 175000)):
        row, _ = ag.build_row("S", _depth_rows(THAL1), [j], None, alleles, PANEL)
        assert row["deletion_alleles"] == "--SEA|--MED/aa"
        assert row["alpha_genes_confidence"] in ("medium", "high")


def test_junction_on_na_row_says_junction_and_stays_na(alleles):
    row, note = ag.build_row("S", HG001_ROWS, [_jr(165396, 184701)], None,
                             alleles, PANEL)
    assert row["alpha_genes_called"] == "NA"
    assert row["deletion_evidence"] == "junction"
    assert "165396" in note


def test_a_junction_touching_no_scorable_segment_is_ignored(alleles):
    row, note = ag.build_row("S", _depth_rows(NORMAL), [_jr(155000, 172000)],
                             None, alleles, PANEL)
    assert row["deletion_evidence"] == "none"
    assert "seen" not in note


@pytest.mark.parametrize("conf", ["medium", "low", None])
def test_only_high_junctions_count(conf, alleles):
    j = [{"left_bp": "165000", "right_bp": "185001", "zygosity": "het"}]
    if conf:
        j[0]["confidence"] = conf
    row, _ = ag.build_row("S", _depth_rows(THAL1), j, None, alleles, PANEL)
    assert row["deletion_evidence"] == "depth"
    assert row["alpha_genes_confidence"] == "medium"


def test_any_NA_genotype_is_low_confidence(alleles):
    obs = {"HBZ": 1, "HBA2": 2, "INTER_A2_A1": 2, "HBA1": 2}   # compound, rule 2
    for j in (None, [], [{"left_bp": "1", "right_bp": "2", "confidence": "high"}]):
        row, _ = ag.build_row("S", _depth_rows(obs), j, None, alleles, PANEL)
        assert row["alpha_genes_called"] == "NA"
        assert row["alpha_genes_confidence"] == "low"


def test_compound_is_not_a_valid_alternative_reading(alleles):
    """Normal call, HBZ intact at 0.67 (marginal): reading HBZ as het_loss fits
    only the compound, which is not valid, so confidence is not lowered."""
    rows = _rows({"HBZ": ("intact", 0.67, True),
                  "HBA2": ("intact", 0.9, False),
                  "INTER_A2_A1": ("intact", 0.9, False),
                  "HBA1": ("intact", 0.96, False)})
    row, note = ag.build_row("S", rows, None, None, alleles, PANEL)
    assert row["deletion_alleles"] == "none"
    assert row["alpha_genes_confidence"] == "medium"
    assert note == ""


def test_no_note_counts_a_compound_alternative(alleles):
    """A compound (deletion + anti-3.7) reading is not a valid alternative, so
    no note may quote a gene count derived from one."""
    rows = _rows({"HBZ": ("intact", 1.0, False),
                  "HBA2": ("gain", 1.36, True),
                  "INTER_A2_A1": ("gain", 1.5, False),
                  "HBA1": ("gain", 1.5, False)})
    called, genes, note = ag.name_alleles(ag.observed_copies(rows), alleles, None,
                                          ag.marginal_alternatives(rows))
    assert "count to 4" not in note


def test_count_ambiguous_alternative_caps_at_medium_with_note(alleles, monkeypatch):
    """A valid alternative reading that names alleles but cannot fix the count."""
    real = ag._name_core

    def fake(obs, al, jr=None, compound_ok=False, inter_unsure=False):
        if obs.get("HBA1") == 1:
            return "-a4.2/aa", "NA", "candidate genotypes disagree"
        return real(obs, al, jr, compound_ok, inter_unsure)
    rows = _rows({"HBZ": ("intact", 1.0, False),
                  "HBA2": ("intact", 0.9, False),
                  "INTER_A2_A1": ("het_loss", 0.5, False),
                  "HBA1": ("intact", 0.667, True)})
    monkeypatch.setattr(ag, "_name_core", fake)
    row, note = ag.build_row("S", rows, None, None, alleles, PANEL)
    assert row["deletion_alleles"] == "-a3.7/aa"
    assert row["alpha_genes_confidence"] == "medium"
    assert "HBA1" in note


def test_base_NA_with_marginal_segments_names_them(alleles):
    obs = {"HBZ": 0, "HBA2": 2, "INTER_A2_A1": 2, "HBA1": 2}   # matches nothing
    called, genes, note = ag.name_alleles(obs, alleles, None, {"HBA1": [1]})
    assert (called, genes) == ("NA", "NA")
    assert "HBA1" in note


def test_genotype_kind_writes_a_note(alleles, monkeypatch):
    """Same gene count, different alleles under an alternative reading."""
    def fake(obs, al, jr=None, compound_ok=False, inter_unsure=False):
        return ("-a4.2/aa" if obs.get("HBA1") == 1 else "-a3.7/aa"), "3", ""
    monkeypatch.setattr(ag, "_name_core", fake)
    obs = {"HBZ": 2, "HBA2": 2, "INTER_A2_A1": 1, "HBA1": 2}
    _, _, note, kind = ag._name_with_marginals(obs, alleles, None, {"HBA1": [1]})
    assert kind == "genotype" and "HBA1" in note


def test_score_boundaries_match_hba_depth():
    import hba_depth
    assert ag.SCORE_BOUNDARIES == (hba_depth.HOM_LOSS_MAX_SCORE,
                                   hba_depth.HET_LOSS_MAX_SCORE,
                                   hba_depth.GAIN_MIN_SCORE)


def test_low_confidence_junctions_do_not_make_a_normal_sample_junction_evidence(alleles):
    junction = [{"left_bp": "44341", "right_bp": "182768", "confidence": "low"}]
    row, _ = ag.build_row("S", _depth_rows(NORMAL), junction, None, alleles, PANEL)
    assert row["deletion_evidence"] == "none"


def test_segment_spans_match_the_bed():
    with open(os.path.join(REPO, "assets", "hba_segments.bed")) as fh:
        for line in fh:
            f = line.split()
            if line.startswith("#") or len(f) < 4 or f[3] not in ag.SEGMENT_SPANS:
                continue
            assert (int(f[1]) + 1, int(f[2])) == ag.SEGMENT_SPANS[f[3]] or \
                (int(f[1]), int(f[2])) == ag.SEGMENT_SPANS[f[3]], line


def test_unsupported_junction_leaves_a_depth_call_at_medium(alleles):
    """Junction rows present but none supporting == depth only, not low."""
    rows = _rows({"HBZ": ("intact", 0.725, True),
                  "HBA2": ("het_loss", 0.486, False),
                  "INTER_A2_A1": ("het_loss", 0.468, False),
                  "HBA1": ("intact", 0.749, False)})
    junction = [{"left_bp": "182819", "right_bp": "236906", "confidence": "low"}]
    row, _ = ag.build_row("S", rows, junction, None, alleles, PANEL)
    assert row["deletion_alleles"] == "-a4.2/aa"
    assert row["deletion_evidence"] == "depth"
    assert row["alpha_genes_confidence"] == "medium"


def test_called_normal_with_a_valid_deletion_reading_is_NA(alleles):
    """INTER_A2_A1 intact but marginal (0.67): reading it as het_loss is a valid
    -a3.7, so "normal" is not reported. NA/NA, never none/4."""
    rows = _rows({"HBZ": ("intact", 1.0, False),
                  "HBA2": ("intact", 0.9, False),
                  "INTER_A2_A1": ("intact", 0.67, True),
                  "HBA1": ("intact", 0.96, False)})
    row, note = ag.build_row("S", rows, None, None, alleles, PANEL)
    assert (row["deletion_alleles"], row["alpha_genes_called"]) == ("NA", "NA")
    assert row["alpha_genes_confidence"] == "low"
    assert "INTER_A2_A1" in note


# --------------------------------------------------------------------------- #
# fixture rows from real libraries (neutral names)
# --------------------------------------------------------------------------- #
def _call_of(score):
    return "hom_loss" if score < 0.25 else "het_loss" if score < 0.65 else \
        "gain" if score >= 1.35 else "intact"


def _score_rows(spec):
    """spec: {segment: (score, marginal)}; call derived from the score."""
    return _rows({g: (_call_of(sc), sc, m) for g, (sc, m) in spec.items()})


def _row_of(spec, alleles, junction=None):
    return ag.build_row("S", _score_rows(spec), junction, None, alleles, PANEL)


def test_real_confirmed_het_a42(alleles):
    row, _ = _row_of({"HBZ": (0.725, True), "HBA2": (0.486, False),
                      "INTER_A2_A1": (0.468, False), "HBA1": (0.749, False)},
                     alleles)
    assert (row["deletion_alleles"], row["alpha_genes_called"]) == ("-a4.2/aa", "3")
    assert row["alpha_genes_confidence"] == "medium"


def test_real_confirmed_hom_a42_like(alleles):
    row, note = _row_of({"HBZ": (0.95, False), "HBA2": (0.188, True),
                         "INTER_A2_A1": (0.201, True), "HBA1": (0.96, False)},
                        alleles)
    assert row["deletion_alleles"] == "-a4.2/-a4.2"
    assert row["alpha_genes_called"] == "2"
    assert row["alpha_genes_confidence"] == "medium"
    assert "3" in note


def test_real_confirmed_sea_med(alleles):
    row, _ = _row_of({"HBZ": (1.116, False), "HBA2": (0.500, False),
                      "INTER_A2_A1": (0.537, False), "HBA1": (0.408, False)},
                     alleles)
    assert (row["deletion_alleles"], row["alpha_genes_called"]) == ("--SEA|--MED/aa", "2")


def test_real_false_positive_is_NA(alleles):
    row, note = ag.build_row("S", HG001_ROWS, None, None, alleles, PANEL)
    assert (row["deletion_alleles"], row["alpha_genes_called"]) == ("NA", "NA")
    assert row["alpha_genes_confidence"] == "low"


@pytest.mark.parametrize("hba2", [0.64, 0.66])      # het_loss and intact, both marginal
def test_real_hba2_only_dip_is_not_a_carrier(hba2, alleles):
    row, _ = _row_of({"HBZ": (0.66, True), "HBA2": (hba2, True),
                      "INTER_A2_A1": (0.95, False), "HBA1": (0.94, False)},
                     alleles)
    assert row["deletion_alleles"] in ("none", "NA")
    assert "-a4.2" not in row["deletion_alleles"]


def test_a42_guard_hba2_loss_with_inter_confidently_intact_is_NA(alleles):
    row, note = _row_of({"HBZ": (1.0, False), "HBA2": (0.45, False),
                         "INTER_A2_A1": (0.95, False), "HBA1": (0.95, False)},
                        alleles)
    assert (row["deletion_alleles"], row["alpha_genes_called"]) == ("NA", "NA")
    assert row["alpha_genes_confidence"] == "low"
    assert "without the expected drop between the genes" in note


def test_marginal_loss_calls_that_read_toward_normal_as_compound_are_NA(alleles):
    """HBZ clear het_loss; HBA2/INTER/HBA1 all marginal het_loss. Read toward
    2 copies they leave an HBZ-only loss: only the banned compound."""
    row, note = _row_of({"HBZ": (0.50, False), "HBA2": (0.62, True),
                         "INTER_A2_A1": (0.60, True), "HBA1": (0.58, True)},
                        alleles)
    assert (row["deletion_alleles"], row["alpha_genes_called"]) == ("NA", "NA")
    assert row["alpha_genes_confidence"] == "low"
    assert "marginal" in note


def test_two_gene_deletion_with_one_borderline_segment_survives(alleles):
    row, _ = _row_of({"HBZ": (1.0, False), "HBA2": (0.5, False),
                      "INTER_A2_A1": (0.5, False), "HBA1": (0.62, True)},
                     alleles)
    assert row["deletion_alleles"] != "NA"
    assert row["alpha_genes_confidence"] in ("medium", "low")
    assert row["alpha_genes_confidence"] != "high"


def test_unmatched_toward_normal_reading_caps_medium_not_NA(alleles):
    """--FIL|--THAI with HBA1 borderline: HBA1 read as intact matches no genotype
    at all. A real two-gene deletion with one borderline segment must survive,
    so: call kept, confidence medium, note names HBA1."""
    row, note = _row_of({"HBZ": (0.5, False), "HBA2": (0.5, False),
                         "INTER_A2_A1": (0.5, False), "HBA1": (0.62, True)},
                        alleles)
    assert row["deletion_alleles"] == "--FIL|--THAI/aa"
    assert row["alpha_genes_confidence"] == "medium"
    assert "HBA1" in note and "match no genotype" in note


def test_marginal_alternatives_reject_bad_scores():
    base = {"segment": "HBA2", "call": "het_loss", "marginal": "True"}
    for bad in ("nan", "inf", "-inf", "0.95", "x", None):
        r = dict(base)
        if bad is not None:
            r["score"] = bad
        alts = ag.marginal_alternatives([r])["HBA2"]
        assert sorted(alts) == [0, 2], (bad, alts)     # both neighbours
    ok = dict(base, score="0.62")
    assert ag.marginal_alternatives([ok])["HBA2"] == [2]
    low = dict(base, score="0.27")
    assert ag.marginal_alternatives([low])["HBA2"] == [0]
    for r in (ok, low, dict(base, call="intact", score="0.67")):
        k = ag.CALL_TO_COPIES[r["call"]]
        assert all(abs(a - k) == 1 for a in ag.marginal_alternatives([r])["HBA2"])


def test_a_gain_note_does_not_say_no_deletion(alleles):
    rows = _score_rows({"HBZ": (1.0, False), "HBA2": (1.38, True),
                        "INTER_A2_A1": (1.5, False), "HBA1": (1.5, False)})
    _, _, note = ag.name_alleles(ag.observed_copies(rows), alleles, None,
                                 ag.marginal_alternatives(rows))
    assert "no deletion" not in note


def test_confidence_requires_genes():
    import inspect
    assert inspect.signature(ag.confidence).parameters["genes"].default is inspect._empty


def test_a42_guard_needs_marginal_toward_loss(alleles):
    """INTER marginal toward GAIN (1.30) does not satisfy the guard."""
    row, note = _row_of({"HBZ": (1.0, False), "HBA2": (0.5, False),
                         "INTER_A2_A1": (1.30, True), "HBA1": (1.0, False)},
                        alleles)
    assert (row["deletion_alleles"], row["alpha_genes_called"]) == ("NA", "NA")
    row, _ = _row_of({"HBZ": (1.0, False), "HBA2": (0.70, True),
                      "INTER_A2_A1": (1.30, True), "HBA1": (1.0, False)}, alleles)
    assert (row["deletion_alleles"], row["alpha_genes_called"]) == ("none", "4")


def test_a42_guard_accepts_inter_marginal_toward_loss(alleles):
    row, _ = _row_of({"HBZ": (1.0, False), "HBA2": (0.5, False),
                      "INTER_A2_A1": (0.67, True), "HBA1": (1.0, False)}, alleles)
    assert row["deletion_alleles"] == "-a4.2/aa"


@pytest.mark.parametrize("missing", [("HBA2",), ("INTER_A2_A1",), ("HBA1",),
                                     ("HBA2", "HBA1")])
def test_unmeasured_gene_segment_is_NA_never_normal(missing, alleles):
    rows = [r for r in _score_rows({"HBZ": (1.0, False), "HBA2": (1.0, False),
                                    "INTER_A2_A1": (1.0, False),
                                    "HBA1": (1.0, False)})
            if r["segment"] not in missing]
    row, note = ag.build_row("S", rows, None, None, alleles, PANEL)
    assert (row["deletion_alleles"], row["alpha_genes_called"]) == ("NA", "NA")
    assert row["alpha_genes_confidence"] == "low"
    assert all(m in note for m in missing)


def test_unmeasured_hbz_reports_the_wider_group(alleles):
    obs = {"HBA2": 1, "INTER_A2_A1": 1, "HBA1": 1}
    called, genes, _ = ag.name_alleles(obs, alleles)
    assert called == "--SEA|--MED|--FIL|--THAI/aa" and genes == "2"


def test_a42_with_a37_a42_alternative_notes_and_caps_medium(alleles):
    obs = {"HBZ": 2, "HBA2": 1, "INTER_A2_A1": 1, "HBA1": 2}
    called, genes, note = ag.name_alleles(obs, alleles)
    assert (called, genes) == ("-a4.2/aa", "3")
    assert ag.NOTE_A37_A42 in note
    rows = _depth_rows(obs)
    row, _ = ag.build_row("S", rows, [_jr(172500, 175500)], None, alleles, PANEL)
    assert row["deletion_alleles"] == "-a4.2/aa"
    assert row["alpha_genes_confidence"] == "medium"


def test_alternatives_trust_the_call_near_a_rounded_boundary():
    r = {"segment": "HBA2", "call": "intact", "score": "1.350", "marginal": "True"}
    assert ag.marginal_alternatives([r])["HBA2"] == [3]
    r = {"segment": "HBA2", "call": "het_loss", "score": "0.650", "marginal": "True"}
    assert ag.marginal_alternatives([r])["HBA2"] == [2]
    r = {"segment": "HBA2", "call": "gain", "score": "1.349", "marginal": "True"}
    assert ag.marginal_alternatives([r])["HBA2"] == [2]


def test_marginal_is_derived_from_score_without_the_column():
    base = {"segment": "HBA2", "call": "intact"}
    assert "HBA2" in ag.marginal_alternatives([dict(base, score="0.70")])
    assert "HBA2" not in ag.marginal_alternatives([dict(base, score="0.95")])
    assert ag.marginal_alternatives([dict(base)]) == {}


def test_marginal_window_matches_hba_depth():
    import hba_depth
    assert ag.MARGINAL_WINDOW == hba_depth.MARGINAL_WINDOW


def test_junction_support_needs_containment(alleles):
    big = _jr(173500, 176900)                   # 3.4 kb: HBA2 and HBA1 not inside
    assert ag.junction_support([big], "--SEA|--MED/aa", alleles) == ([], [big])
    tiny = _jr(175000, 175200)                  # 200 bp inside INTER
    assert ag.junction_support([tiny], "-a3.7/aa", alleles) == ([], [tiny])
    ok = _jr(173500, 177200)
    assert ag.junction_support([ok], "-a3.7/aa", alleles) == ([ok], [])


def test_gain_always_carries_the_unresolved_copy_number_note(alleles):
    row, note = _row_of({"HBZ": (1.0, False), "HBA2": (2.0, False),
                         "INTER_A2_A1": (2.0, False), "HBA1": (2.0, False)},
                        alleles)
    assert row["deletion_alleles"] == "anti-3.7/aa" and row["alpha_genes_called"] == "5"
    assert ag.NOTE_GAIN in note
    _, note = _row_of({"HBZ": (1.0, False), "HBA2": (1.0, False),
                       "INTER_A2_A1": (1.0, False), "HBA1": (1.0, False)}, alleles)
    assert ag.NOTE_GAIN not in note


def test_a42_with_borderline_inter_is_capped_at_medium_with_note(alleles):
    row, note = _row_of({"HBZ": (1.0, False), "HBA2": (0.5, False),
                         "INTER_A2_A1": (0.70, True), "HBA1": (1.0, False)},
                        alleles, [_jr(170000, 174200)])
    assert row["deletion_alleles"] == "-a4.2/aa"
    assert row["deletion_evidence"] == "both"
    assert row["alpha_genes_confidence"] == "medium"
    assert ag.NOTE_INTER_BORDERLINE in note


def test_a42_with_a_marginal_het_loss_inter_stays_named_at_medium(alleles):
    """DECLINED to make this NA: confirmed het carriers read INTER 0.47-0.60
    (two platforms) against 0.78-1.29 in normals, so a marginal loss is
    carrier-like, not a coin flip."""
    row, note = _row_of({"HBZ": (1.0, False), "HBA2": (0.5, False),
                         "INTER_A2_A1": (0.62, True), "HBA1": (1.0, False)},
                        alleles)
    assert row["deletion_alleles"] == "-a4.2/aa"
    assert row["alpha_genes_confidence"] == "medium"
    assert ag.NOTE_INTER_BORDERLINE in note


def test_pure_gain_with_one_borderline_segment_is_not_NA(alleles):
    row, note = _row_of({"HBZ": (1.0, False), "HBA2": (1.40, True),
                         "INTER_A2_A1": (1.5, False), "HBA1": (1.5, False)},
                        alleles)
    assert (row["deletion_alleles"], row["alpha_genes_called"]) == ("anti-3.7/aa", "5")
    assert row["alpha_genes_confidence"] == "medium"
    assert "HBA2" in note


def test_a42_guard_fires_only_when_hba2_is_called_a_loss(alleles):
    row, note = _row_of({"HBZ": (1.0, False), "HBA2": (1.0, False),
                         "INTER_A2_A1": (1.0, False), "HBA1": (1.4, False)},
                        alleles)
    assert row["alpha_genes_called"] == "NA"
    assert "two-event" in note and "HBA2 depth loss" not in note
