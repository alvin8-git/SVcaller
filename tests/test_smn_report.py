import os
import re
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "bin"))
import pytest

def test_classify_affected():
    from smn_report import classify_sma
    result = classify_sma(smn1_cn=0, smn2_cn=2)
    assert result["status"] == "Affected"
    assert result["smn1_cn"] == 0

def test_classify_carrier():
    from smn_report import classify_sma
    result = classify_sma(smn1_cn=1, smn2_cn=1)
    assert result["status"] == "Carrier"
    assert result["badge_class"] == "badge-warning"
    assert result["smn1_cn"] == 1
    assert result["smn2_cn"] == 1
    assert "carrier" in result["interpretation"].lower()

def test_classify_normal():
    from smn_report import classify_sma
    result = classify_sma(smn1_cn=2, smn2_cn=2)
    assert result["status"] == "Normal"
    assert result["badge_class"] == "badge-success"
    assert result["smn1_cn"] == 2
    assert result["smn2_cn"] == 2
    assert "normal" in result["interpretation"].lower() or "2" in result["interpretation"]

def test_two_plus_zero_flagged():
    from smn_report import detect_two_plus_zero
    assert detect_two_plus_zero(smn1_cn=2, smn1_allele1=2, smn1_allele2=0) is True

def test_two_plus_zero_not_flagged_when_balanced():
    from smn_report import detect_two_plus_zero
    assert detect_two_plus_zero(smn1_cn=2, smn1_allele1=1, smn1_allele2=1) is False


# --- Silent-failure guard (OmniGen clean-bill-of-health regression) ------------
# A crashed SMN caller must never be rendered as a default "Normal (CN=2)" result.
# absent (NO_* sentinel) != empty (0-byte crash placeholder) != populated.

def test_smn_empty_file_raises_not_normal(tmp_path):
    """A 0-byte SMN TSV is a crashed caller, not a Normal result -> must fail loud."""
    from smn_report import parse_smn_tsv, SmnInputError
    empty = tmp_path / "HG.smn.tsv"; empty.write_text("")
    with pytest.raises(SmnInputError):
        parse_smn_tsv(str(empty))


def test_smn_header_only_raises(tmp_path):
    """A header with no data row means the caller produced no measurement -> fail loud."""
    from smn_report import parse_smn_tsv, SmnInputError
    hdr = tmp_path / "HG.smn.tsv"
    hdr.write_text("Sample\tSMN1_CN\tSMN2_CN\n")
    with pytest.raises(SmnInputError):
        parse_smn_tsv(str(hdr))


def test_smn_render_empty_file_raises(tmp_path):
    """End-to-end: the section renderer must not emit 'Normal' for an empty input."""
    from smn_report import render_html_section, SmnInputError
    empty = tmp_path / "HG.smn.tsv"; empty.write_text("")
    with pytest.raises(SmnInputError):
        render_html_section("HG_TEST", str(empty))


def test_smn_healthy_file_parses(tmp_path):
    """Header + one data row (the healthy shape) parses cleanly."""
    from smn_report import parse_smn_tsv
    good = tmp_path / "HG.smn.tsv"
    good.write_text("Sample\tSMN1_CN\tSMN2_CN\tConfidence\n"
                    "HG_TEST\t1\t3\tHIGH\n")
    parsed = parse_smn_tsv(str(good))
    assert parsed["smn1_cn"] == 1 and parsed["smn2_cn"] == 3


def test_smn_absent_sentinel_is_skip(tmp_path):
    """A Nextflow NO_* sentinel is a legitimate skip, not a failure -> no raise."""
    from smn_report import parse_smn_tsv
    sentinel = tmp_path / "NO_SMN"; sentinel.write_text("")
    parsed = parse_smn_tsv(str(sentinel))
    assert parsed["smn1_cn"] is None  # neutral/unknown, not a fabricated CN=2


# --- Low-depth no-call (SMN caller patch, svcaller/smncopynum:1.2) -------------
NOCALL_TSV = ("Sample\tisSMA\tisCarrier\tSMN1_CN\tSMN2_CN\tSMN2delta7-8_CN\tTotal_CN_raw\t"
              "Full_length_CN_raw\tg.27134T>G_CN\tSMN1_CN_raw\tInfo\n"
              "HG_TEST\tNone\tNone\tNone\tNone\t0\t4.1\t3.9\t0\t1.6,2.5\tNoCall:InsufficientDepth\n")


def test_smn_nocall_parses_as_no_copy_number(tmp_path):
    from smn_report import parse_smn_tsv
    f = tmp_path / "HG.smn.tsv"; f.write_text(NOCALL_TSV)
    parsed = parse_smn_tsv(str(f))
    assert parsed["smn1_cn"] is None and parsed["smn2_cn"] is None
    assert parsed["info"] == "NoCall:InsufficientDepth"


def test_smn_nocall_render_is_not_normal(tmp_path):
    """A no-call must render as 'No call', never the default 'Normal / CN=2'."""
    from smn_report import render_html_section
    f = tmp_path / "HG.smn.tsv"; f.write_text(NOCALL_TSV)
    html = render_html_section("HG_TEST", str(f))
    assert "No call" in html and "NoCall:InsufficientDepth" in html
    assert "Normal" not in html and "badge-success" not in html



def test_smn_nocall_render_escapes_html(tmp_path):
    from smn_report import render_html_section
    f = tmp_path / "HG.smn.tsv"; f.write_text(NOCALL_TSV.replace("NoCall:InsufficientDepth", "<b>x</b>"))
    html = render_html_section("<i>S</i>", str(f))
    assert "<b>x</b>" not in html and "<i>S</i>" not in html and "&lt;b&gt;" in html


def test_xlsx_smn_sheet_nocall_interpretation(tmp_path):
    openpyxl = pytest.importorskip("openpyxl")
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "bin"))
    from html_report import export_xls
    tsv = tmp_path / "HG.smn.tsv"; tsv.write_text(NOCALL_TSV)
    out = tmp_path / "o.xlsx"
    export_xls([], "", [], str(tsv), str(out))
    rows = list(openpyxl.load_workbook(out)["SMN"].iter_rows(values_only=True))
    assert rows[1][-1] == "No call (NoCall:InsufficientDepth)"
    good = tmp_path / "G.smn.tsv"; good.write_text("Sample\tSMN1_CN\tSMN2_CN\nS\t1\t3\n")
    export_xls([], "", [], str(good), str(out))
    assert list(openpyxl.load_workbook(out)["SMN"].iter_rows(values_only=True))[1][-1] == "SMA Carrier"


# --- SMN depth floor (bin/smn_depth_floor.py, params.smn_min_depth) -------------
HDR = ("Sample\tisSMA\tisCarrier\tSMN1_CN\tSMN2_CN\tSMN2delta7-8_CN\tTotal_CN_raw\t"
       "Full_length_CN_raw\tg.27134T>G_CN\tSMN1_CN_raw\tInfo\n")
WRONG_SPLIT_3_1 = HDR + "S\tFalse\tFalse\t3\t1\t0\t4.06\t3.9\t0\t1.6,2.5\tPASS:AgreeWithSum\n"
SMA_AFF = HDR + "S\tTrue\tFalse\t0\t3\t0\t3.0\t2.9\t0\t0.1,0.2\tPASS:Majority\n"
CARRIER = HDR + "S\tFalse\tTrue\t1\t1\t0\t2.0\t1.9\t0\t0.9,1.1\tPASS:Majority\n"


def _run(tmp_path, text, depth, min_depth=18):
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "bin"))
    import smn_depth_floor as m
    f = tmp_path / "S.smn.tsv"; f.write_bytes(text.encode())
    rc = m.main([str(f), "--depth", str(depth), "--min-depth", str(min_depth)])
    return rc, f


def test_floor_below_is_nocall_and_raw_kept(tmp_path):
    from smn_report import parse_smn_tsv
    rc, f = _run(tmp_path, WRONG_SPLIT_3_1, 14.7)
    assert rc == 0
    row = f.read_text().split("\n")[1].split("\t")
    assert row[1:5] == ["None"] * 4 and row[-1] == "NoCall:BelowValidatedDepth"
    assert row[6:10] == ["4.06", "3.9", "0", "1.6,2.5"]          # raw-derived columns untouched
    assert (tmp_path / "S.smn.raw.tsv").read_text() == WRONG_SPLIT_3_1     # caller row preserved
    p = parse_smn_tsv(str(f))
    assert p["smn1_cn"] is None and p["info"] == "NoCall:BelowValidatedDepth"


def test_floor_nocall_renders_no_call_not_normal(tmp_path):
    from smn_report import render_html_section
    _, f = _run(tmp_path, WRONG_SPLIT_3_1, 10)
    html = render_html_section("S", str(f))
    assert "No call" in html and "NoCall:BelowValidatedDepth" in html and "Normal" not in html


def test_floor_at_or_above_is_byte_identical(tmp_path):
    for depth in (18, 30.9):
        rc, f = _run(tmp_path, WRONG_SPLIT_3_1, depth)
        assert rc == 0 and f.read_bytes() == WRONG_SPLIT_3_1.encode()
        assert not (tmp_path / "S.smn.raw.tsv").exists()


def test_floor_real_sma_and_carrier_above_floor_unchanged(tmp_path):
    for txt in (SMA_AFF, CARRIER):
        _, f = _run(tmp_path, txt, 25)
        assert f.read_bytes() == txt.encode()


def test_floor_existing_nocall_passes_through(tmp_path):
    _, f = _run(tmp_path, NOCALL_TSV, 9.8)
    assert f.read_bytes() == NOCALL_TSV.encode()
    assert not (tmp_path / "S.smn.raw.tsv").exists()


def test_floor_is_idempotent_and_keeps_first_raw(tmp_path):
    rc, f = _run(tmp_path, WRONG_SPLIT_3_1, 12)
    first = f.read_bytes()
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "bin"))
    import smn_depth_floor as m
    m.main([str(f), "--depth", "12"])
    assert f.read_bytes() == first and (tmp_path / "S.smn.raw.tsv").read_text() == WRONG_SPLIT_3_1


def test_floor_appends_info_column_for_old_tsv_and_reads_mosdepth(tmp_path):
    import smn_depth_floor as m
    old = "Sample\tisSMA\tisCarrier\tSMN1_CN\tSMN2_CN\nS\tFalse\tFalse\t2\t2\n"
    new, floored = m.floor_text(old, 10, 18)
    assert floored and new == "Sample\tisSMA\tisCarrier\tSMN1_CN\tSMN2_CN\tInfo\nS\tNone\tNone\tNone\tNone\tNoCall:BelowValidatedDepth\n"
    s = tmp_path / "S.mosdepth.summary.txt"
    s.write_text("chrom\tlength\tbases\tmean\tmin\tmax\ntotal\t100\t1470\t14.7\t0\t9\ntotal_region\t1\t1\t99\t0\t1\n")
    assert m.mean_depth(str(s)) == 14.7


def test_floor_missing_depth_fails_loud(tmp_path):
    import smn_depth_floor as m
    f = tmp_path / "S.smn.tsv"; f.write_text(WRONG_SPLIT_3_1)
    assert m.main([str(f), "--summary", str(tmp_path / "absent.txt")]) == 2
    assert f.read_text() == WRONG_SPLIT_3_1


def test_floor_param_and_wiring_separate_from_min_depth():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    cfg = open(os.path.join(root, "nextflow.config")).read()
    assert re.search(r"^\s*smn_min_depth\s*=\s*18\b", cfg, re.M) and re.search(r"^\s*min_depth\s*=\s*25\b", cfg, re.M)
    assert "params.smn_min_depth" in open(os.path.join(root, "subworkflows/smn_calling.nf")).read()


def test_floor_symlinked_input_same_basename_out_leaves_source_untouched(tmp_path):
    import smn_depth_floor as m
    src = tmp_path / "cache"; src.mkdir()
    (src / "S.smn.tsv").write_text(WRONG_SPLIT_3_1)
    stage = tmp_path / "work" / "in"; stage.mkdir(parents=True)
    os.symlink(src / "S.smn.tsv", stage / "S.smn.tsv")
    out = tmp_path / "work" / "S.smn.tsv"
    assert m.main([str(stage / "S.smn.tsv"), "--depth", "10", "--out", str(out),
                   "--raw", str(tmp_path / "work" / "S.smn.raw.tsv")]) == 0
    assert (src / "S.smn.tsv").read_text() == WRONG_SPLIT_3_1 and not os.path.islink(out)
    assert "NoCall:BelowValidatedDepth" in out.read_text()


def test_floor_rerun_reevaluates_from_raw(tmp_path):
    import smn_depth_floor as m
    f = tmp_path / "S.smn.tsv"; f.write_text(WRONG_SPLIT_3_1)
    m.main([str(f), "--depth", "10"])
    m.main([str(f), "--depth", "30"])           # wrong depth corrected: raw row restored
    assert f.read_text() == WRONG_SPLIT_3_1
    m.main([str(f), "--depth", "20", "--min-depth", "25"])  # changed floor
    assert "NoCall:BelowValidatedDepth" in f.read_text()


def test_floor_bad_inputs_exit2_file_untouched(tmp_path):
    import smn_depth_floor as m
    f = tmp_path / "S.smn.tsv"; f.write_text(WRONG_SPLIT_3_1)
    bad = tmp_path / "b.txt"; bad.write_text("total\t1\t2\n")
    assert m.main([str(f), "--summary", str(bad)]) == 2
    for d in ("nan", "inf"):
        assert m.main([str(f), "--depth", d]) == 2
    nocols = tmp_path / "N.smn.tsv"; nocols.write_text("Sample\tfoo\nS\t1\n")
    assert m.main([str(nocols), "--depth", "5"]) == 2
    assert f.read_text() == WRONG_SPLIT_3_1 and nocols.read_text() == "Sample\tfoo\nS\t1\n"


def test_floor_crlf_single_info_column(tmp_path):
    import smn_depth_floor as m
    text = WRONG_SPLIT_3_1.replace("\n", "\r\n")
    new, floored = m.floor_text(text, 5, 18)
    assert floored and new.count("Info") == 1 and new.count("\r\n") == 2
    assert new.split("\r\n")[1].split("\t")[-1] == "NoCall:BelowValidatedDepth"
