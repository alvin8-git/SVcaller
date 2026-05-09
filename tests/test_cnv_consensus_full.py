import sys
sys.path.insert(0, "/data/alvin/SVcaller/bin")
from cnv_consensus import CNVSegment, merge, load_cnvpytor, load_gatk
import tempfile, textwrap
from pathlib import Path
import pytest

def make_seg(chrom="chr1", start=1000, end=5000, cn=1, svtype="DEL", caller="CNVpytor", qual=None):
    return CNVSegment(chrom, start, end, cn, svtype, caller, qual)

def test_both_callers_gives_high_confidence():
    cnv = [make_seg()]
    gatk = [make_seg(caller="GATK_gCNV", qual=50.0)]
    result = merge(cnv, gatk)
    assert len(result) == 1
    assert result[0]["confidence"] == "HIGH"
    assert result[0]["caller_support"] == "BOTH"

def test_gatk_only_high_quality_included():
    cnv = []
    gatk = [make_seg(caller="GATK_gCNV", qual=35.0)]
    result = merge(cnv, gatk)
    assert len(result) == 1
    assert result[0]["confidence"] == "MEDIUM"

def test_gatk_only_low_quality_excluded():
    cnv = []
    gatk = [make_seg(caller="GATK_gCNV", qual=20.0)]
    result = merge(cnv, gatk)
    assert len(result) == 0

def test_cnvpytor_only_excluded():
    cnv = [make_seg()]
    gatk = []
    result = merge(cnv, gatk)
    assert len(result) == 0

def test_different_svtype_not_merged():
    cnv  = [make_seg(svtype="DEL")]
    gatk = [make_seg(svtype="DUP", caller="GATK_gCNV", qual=50.0)]
    result = merge(cnv, gatk)
    assert any(r["svtype"] == "DUP" for r in result)
    assert not any(r["svtype"] == "DEL" for r in result)
