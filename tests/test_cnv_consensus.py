"""Tests for cnv_consensus: overlap logic, GATK parsing, and the merge tiers.

Regression coverage for the 2026-06 CNV-empty bug: the converted GATK TSV uses
CALL_COPY_NUMBER/QUALITY columns; if the raw .seg (CONTIG/.../CALL) is passed
instead, every row must be dropped AND a warning emitted — never a silent 0.
"""
import sys
sys.path.insert(0, "/data/alvin/SVcaller/bin")
import cnv_consensus as c


# --- reciprocal overlap -------------------------------------------------------

def test_full_overlap():
    assert c.reciprocal_overlap(100, 200, 100, 200) == 1.0

def test_no_overlap():
    assert c.reciprocal_overlap(100, 200, 300, 400) == 0.0

def test_partial_overlap():
    assert abs(c.reciprocal_overlap(100, 200, 150, 250) - 0.5) < 1e-9

def test_contained():
    assert c.reciprocal_overlap(100, 400, 150, 250) == 1.0


# --- load_gatk: converted TSV (correct contract) ------------------------------

def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return str(p)

def test_load_gatk_parses_converted_tsv(tmp_path):
    tsv = ("CONTIG\tSTART\tEND\tCALL_COPY_NUMBER\tQUALITY\n"
           "chr1\t1000\t5000\t1\t50\n"      # deletion
           "chr1\t8000\t9000\t3\t50\n"      # duplication
           "chr2\t2000\t3000\t2\t50\n")     # neutral -> skipped
    segs = c.load_gatk(_write(tmp_path, "g.tsv", tsv))
    assert len(segs) == 2
    assert {s.svtype for s in segs} == {"DEL", "DUP"}
    assert all(s.cn != 2 for s in segs)

def test_load_gatk_raw_seg_fails_loud(tmp_path, capsys):
    # Regression: the raw CallCopyRatioSegments .seg has no CALL_COPY_NUMBER column.
    # Every row must drop AND a WARNING must be emitted (not a silent empty result).
    seg = ("CONTIG\tSTART\tEND\tNUM_POINTS_COPY_RATIO\tMEAN_LOG2_COPY_RATIO\tCALL\n"
           "chr1\t1000\t5000\t50\t-0.8\t-\n"
           "chr1\t8000\t9000\t40\t0.7\t+\n")
    segs = c.load_gatk(_write(tmp_path, "raw.seg", seg))
    assert segs == []
    assert "WARNING" in capsys.readouterr().err


# --- merge: three confidence tiers --------------------------------------------

def _seg(chrom, start, end, cn, svtype, caller, q=None):
    return c.CNVSegment(chrom, start, end, cn, svtype, caller, q)

def test_merge_both_high_on_overlap():
    cp = [_seg("chr1", 1000, 5000, 1, "DEL", "CNVpytor")]
    ga = [_seg("chr1", 1000, 5000, 1, "DEL", "GATK_gCNV", 50)]
    out = c.merge(cp, ga)
    assert any(r["caller_support"] == "BOTH" and r["confidence"] == "HIGH" for r in out)

def test_merge_gatk_only_medium():
    cp = []
    ga = [_seg("chr1", 1000, 5000, 1, "DEL", "GATK_gCNV", 50)]
    out = c.merge(cp, ga)
    assert [r["confidence"] for r in out] == ["MEDIUM"]

def test_merge_cnvpytor_only_low_size_gated():
    # small unmatched CNVpytor call is suppressed; large one surfaces as LOW
    cp = [_seg("chr1", 1000, 5000, 1, "DEL", "CNVpytor"),               # 4 kb -> dropped
          _seg("chr3", 0, 2_000_000, 3, "DUP", "CNVpytor")]            # 2 Mb -> LOW
    out = c.merge(cp, [], cnvpytor_only_min_bp=1_000_000)
    lows = [r for r in out if r["confidence"] == "LOW"]
    assert len(lows) == 1 and lows[0]["chrom"] == "chr3"
