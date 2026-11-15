"""The HTML report must not contradict OmniGen's repeat-expansion grading (owner decisions 2026-11-15).

The curation tables are carried as copies (assets/str_*.tsv). When $OMNIGEN_PANEL_SOURCES points at OmniGen's
prototype/panel_sources directory the copies are compared with the originals; otherwise that comparison is skipped."""
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import pytest
import html_report as hr
import str_grading_scope as sgs

SRC = os.environ.get("OMNIGEN_PANEL_SOURCES", "")
need_src = pytest.mark.skipif(not SRC or not os.path.isdir(SRC), reason="OMNIGEN_PANEL_SOURCES not set")
FX = os.path.join(REPO, "tests", "fixtures")


def _body(path):
    return [l for l in open(path).read().splitlines() if l.strip() and not l.startswith("#")]


@need_src
def test_retired_loci_are_the_audit_retire_verdicts():
    rows = [l.split("\t") for l in _body(os.path.join(SRC, "repeat_audit.tsv"))]
    h = rows[0]
    want = {r[h.index("gene")] for r in rows[1:] if r[h.index("verdict")] == "RETIRE"}
    assert set(sgs.retired()) == want and want


@need_src
def test_convention_and_not_tested_copies_are_verbatim():
    assert _body(sgs.CONVENTION_FILE) == _body(os.path.join(SRC, "str_count_convention.tsv"))
    assert _body(sgs.NOT_TESTED_FILE) == _body(os.path.join(SRC, "str_not_tested.tsv"))


def test_narrow_gap_is_computed_from_the_table_and_pabpn1_has_its_convention():
    rng = hr._load_str_ranges()
    ng = sgs.narrow_gap(rng)
    assert "PABPN1" in ng and "PABPN1" in sgs.convention_loci() and "HTT" not in ng
    sc = sgs.load_scope(rng, True)
    assert sgs.not_graded_reason("PABPN1", sc) is None            # narrow gap, but its convention is cited
    assert sgs.not_graded_reason("EIF4A3", sc) == sgs.CONVENTION_REASON
    assert sgs.not_graded_reason("NIPA1", sc) == sgs.RETIRED_REASON
    assert sgs.not_graded_reason("EIF4A3", sgs.load_scope(hr._load_str_ranges(hr.STR_RANGES_LEGACY_PATH), False)).endswith("not tested")


def test_bound_wording():
    assert sgs.unresolved_bound(">=52") == 52 and sgs.unresolved_bound("17") is None and sgs.unresolved_bound(">=0") == 0
    assert sgs.bound_text(0) == "not measurable" and sgs.bound_text(52) == "at least 52 repeats"


def _summary(tmp_path, rows):
    hdr = ["locus", "genotype_display", "allele1_sizing", "allele2_sizing", "variant_id"]
    p = tmp_path / "s.tsv"
    p.write_text("\t".join(hdr) + "\n" + "".join("\t".join(r) + "\n" for r in rows))
    return hr._load_str_summary(str(p))[0]


def _grade(tmp_path, rows):
    g = hr._load_str_graded()
    return {l["locus"]: l for l in hr.parse_str_loci("", hr._load_str_ranges(), g, _summary(tmp_path, rows))}


def test_summary_counts_replace_raw_vcf_counts_pabpn1_offset(tmp_path):
    L = _grade(tmp_path, [["PABPN1", "10/10", "spanned", "spanned", "PABPN1"]])
    assert L["PABPN1"]["repcn"] == "10/10" and L["PABPN1"]["status"] == "NORMAL"


def test_vcf_without_summary_applies_the_cited_offset():
    g = hr._load_str_graded()
    vcf = os.path.join(FX, "hg002_eh_new_catalogue.str.vcf")
    off = hr.eh_sizing.load_offsets(str(hr.STR_OFFSETS_PATH))
    assert off[("PABPN1", "PABPN1")] == -1
    assert hr._offset_repcn("11/11", -1) == "10/10" and hr._offset_repcn("33/.", -1) == "32/."
    assert {l["locus"] for l in hr.parse_str_loci(vcf, hr._load_str_ranges(), g, None, off)} == {"CNBP", "FMR1", "HTT"}


def test_unresolved_and_unmeasurable_alleles_are_never_graded(tmp_path):
    L = _grade(tmp_path, [["HTT", "17/>=60", "spanned", "unresolved", "HTT"],
                          ["ATXN3", ">=0/>=0", "unresolved", "unresolved", "ATXN3"],
                          ["ATXN2", "23/>=31", "spanned", "unresolved", "ATXN2"]])
    assert L["HTT"]["status"] == "UNRESOLVED" and L["HTT"]["repcn"] == "17/>=60 (unresolved)"
    assert L["ATXN3"]["status"] == "NOT MEASURABLE" and L["ATXN3"]["repcn"] == "not measurable/not measurable"
    assert L["ATXN2"]["status"] == "UNRESOLVED"


def test_retired_and_unconfirmed_loci_are_listed_not_graded(tmp_path):
    L = _grade(tmp_path, [["NIPA1", "4/18", "spanned", "spanned", "NIPA1"], ["EIF4A3", "20/20", "spanned", "spanned", "EIF4A3"]])
    for k in ("NIPA1", "EIF4A3"):
        assert L[k]["status"] == "NOT GRADED" and L[k]["status_note"] and L[k]["path_min"] == ""
    assert hr.build_known_diseases(list(L.values()), "", "", cnv_syndromes=[]) == []        # never a Finding


def test_a_resolved_expanded_allele_is_still_expanded(tmp_path):
    L = _grade(tmp_path, [["HTT", "17/60", "spanned", "spanned", "HTT"]])
    assert L["HTT"]["status"] == "EXPANDED"


def test_sv_not_assessed_states(tmp_path):
    assert hr.load_sv_not_assessed(None) == ([], "unrecorded")
    p = tmp_path / "n.tsv"
    p.write_text("sample\tcaller\tsvtype\treason\tlimit\n")
    assert hr.load_sv_not_assessed(str(p)) == ([], "complete")
    p.write_text("sample\tcaller\tsvtype\treason\tlimit\nS\tdelly\tINV\ttime_limit_exceeded\t7200s\n")
    rows, st = hr.load_sv_not_assessed(str(p))
    assert st == "rows" and rows[0]["caller"] == "delly"
    p.write_text("")
    with pytest.raises(SystemExit):
        hr.load_sv_not_assessed(str(p))


def test_smn_cell_without_a_count_is_never_a_carrier(tmp_path):
    p = tmp_path / "smn.tsv"
    p.write_text("Sample\tSMN1_CN\tSMN2_CN\nS\t\t2\n")
    f = hr.build_known_diseases([], str(p), "", cnv_syndromes=[])
    assert [x["status"] for x in f] == ["NOT TYPED"]          # shown as not typed, never CARRIER/AFFECTED and never silent


def test_low_depth_note_between_20x_and_30x_only():
    for d, note in (("22.0", True), ("29.9", True), ("30.0", False), ("N/A", False)):
        q = {"mean_depth": d}
        hr._low_depth_note(q)
        assert bool(q.get("low_depth_note")) is note
    assert "2-4 points fewer deletions" in hr.LOW_DEPTH_NOTE


def test_min_depth_default_is_20x():
    import re
    cfg = open(os.path.join(REPO, "nextflow.config")).read()
    assert re.search(r"^\s*min_depth\s*=\s*20\b", cfg, re.M)


def test_ci_follows_the_same_offset_as_the_count(tmp_path):
    assert hr._offset_repci("11-11/12-13", -1) == "10-10/11-12" and hr._offset_repci(".", -1) == "."
    g = hr._load_str_graded()
    vcf = os.path.join(FX, "hg002_eh_new_catalogue.str.vcf")
    off = hr.eh_sizing.load_offsets(str(hr.STR_OFFSETS_PATH))
    raw = {l["locus"]: l for l in hr.parse_str_loci(vcf, hr._load_str_ranges(), g, None, {})}
    cor = {l["locus"]: l for l in hr.parse_str_loci(vcf, hr._load_str_ranges(), g, None, off)}
    assert raw["HTT"]["repci"] == cor["HTT"]["repci"]                 # no offset, nothing moves
    s = tmp_path / "s.tsv"
    s.write_text("locus\tgenotype_display\tallele1_sizing\tallele2_sizing\tvariant_id\tci\nPABPN1\t10/10\tspanned\tspanned\tPABPN1\t10-10/10-10\n")
    L = {l["locus"]: l for l in hr.parse_str_loci("", hr._load_str_ranges(), g, hr._load_str_summary(str(s))[0])}
    assert L["PABPN1"]["repci"] == "10-10/10-10"


def test_legacy_summary_without_vcf_still_lists_the_loci(tmp_path):
    s = tmp_path / "s.tsv"
    s.write_text("locus\tgenotype_display\tallele1_sizing\tallele2_sizing\nHTT\t17/18\tspanned\tspanned\n")
    rng = hr._load_str_ranges(hr.STR_RANGES_LEGACY_PATH)
    L = hr.parse_str_loci("", rng, None, hr._load_str_summary(str(s))[0])
    assert [l["locus"] for l in L] == ["HTT"] and L[0]["repcn"] == "17/18"      # never an empty section


def test_empty_smn1_is_not_typed_never_clean(tmp_path):
    t = tmp_path / "smn.tsv"
    t.write_text("Sample\tisSMA\tisCarrier\tSMN1_CN\tSMN2_CN\tInfo\nx\tFalse\tFalse\t\t\tNoCall:BelowValidatedDepth\n")
    f = hr.build_known_diseases([], str(t), "", cnv_syndromes=[])
    assert [(x["source"], x["status"]) for x in f] == [("SMN", "NOT TYPED")] and "BelowValidatedDepth" in f[0]["detail"]


def test_scope_tables_are_read_once_per_call_site(monkeypatch):
    calls = []
    real = sgs.narrow_gap
    monkeypatch.setattr(sgs, "narrow_gap", lambda r: calls.append(1) or real(r))
    sc = sgs.load_scope(hr._load_str_ranges(), True)
    for _ in range(50):
        sgs.not_graded_reason("HTT", sc)
    assert len(calls) == 1
