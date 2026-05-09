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

def test_classify_normal():
    from smn_report import classify_sma
    result = classify_sma(smn1_cn=2, smn2_cn=2)
    assert result["status"] == "Normal"

def test_two_plus_zero_flagged():
    from smn_report import detect_two_plus_zero
    assert detect_two_plus_zero(smn1_cn=2, smn1_allele1=2, smn1_allele2=0) is True

def test_two_plus_zero_not_flagged_when_balanced():
    from smn_report import detect_two_plus_zero
    assert detect_two_plus_zero(smn1_cn=2, smn1_allele1=1, smn1_allele2=1) is False
