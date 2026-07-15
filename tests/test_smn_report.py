import sys
sys.path.insert(0, "/data/alvin/SVcaller/bin")
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
