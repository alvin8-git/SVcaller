"""Regression tests for CNV report display (2026-06): syndrome-gated Clinical Findings,
top-N-by-size CNV list, and the copy-neutral (CN=2) CNVpytor skip."""
import sys
sys.path.insert(0, "/data/alvin/SVcaller/bin")
import html_report as H
import cnv_consensus as C

# A CNV bed with: a 22q11.2-region DEL (syndrome), a generic large non-syndrome DUP,
# and a tiny CNV clipping the edge of the Williams region (should NOT match).
_BED = """\
#chrom\tstart\tend\tcn\tsvtype\tcaller_support\tconfidence\tquality\tsample
chr22\t18900000\t21500000\t1\tDEL\tBOTH\tHIGH\t50.0\tS1
chr8\t40000000\t90000000\t3\tDUP\tGATK_only\tMEDIUM\t40.0\tS1
chr7\t74100000\t74200000\t3\tDUP\tGATK_only\tMEDIUM\t30.0\tS1
"""


def _bed(tmp_path):
    p = tmp_path / "cnv.bed"
    p.write_text(_BED)
    return str(p)


def test_load_cnv_syndromes_nonempty():
    syn = H.load_cnv_syndromes()
    assert len(syn) >= 20
    assert any("22q11.2" in s["syndrome"] for s in syn)


def test_match_cnv_syndrome_hits_and_misses():
    syn = H.load_cnv_syndromes()
    # 22q11.2 DEL fully covering the region -> match
    assert H.match_cnv_syndrome("chr22", 18900000, 21500000, "DEL", syn) is not None
    # Same region but wrong type for a DEL-only... 22q11.2 is BOTH, so DUP also matches.
    assert H.match_cnv_syndrome("chr22", 18900000, 21500000, "DUP", syn) is not None
    # Tiny 100kb clip of the Williams 7q11.23 region (<50% of region) -> no match
    assert H.match_cnv_syndrome("chr7", 74100000, 74200000, "DUP", syn) is None
    # NF1 region is DEL-only; a DUP there must NOT match
    assert H.match_cnv_syndrome("chr17", 29100000, 30300000, "DUP", syn) is None


def test_known_diseases_only_syndrome_cnvs(tmp_path):
    findings = H.build_known_diseases([], "", _bed(tmp_path))
    cnv = [f for f in findings if f["source"] == "CNV"]
    # Only the 22q11.2 DEL qualifies; the generic 50Mb DUP and the tiny clip do not.
    assert len(cnv) == 1
    assert "22q11.2" in cnv[0]["disease"]
    assert "Copy Number" not in cnv[0]["disease"]   # no generic label


def test_top_cnvs_by_size_ranks_and_flags(tmp_path):
    top, total = H.top_cnvs_by_size(_bed(tmp_path), n=10)
    assert total == 3
    # Largest first: the 50Mb chr8 DUP
    assert top[0]["pos"].startswith("chr8:")
    assert top[0]["size"].endswith("Mb")
    # The 22q11.2 DEL carries a syndrome label
    syn_rows = [c for c in top if c["syndrome"]]
    assert any("22q11.2" in c["syndrome"] for c in syn_rows)


def test_cnvpytor_skips_copy_neutral(tmp_path):
    # CNVpytor TSV: a real dup (cn 3), a copy-neutral call (cn 2 -> skip), a del (cn 1)
    p = tmp_path / "pytor.tsv"
    p.write_text(
        "duplication\tchr1:1000-5000\t4000\t3\t0.01\n"
        "duplication\tchr2:1000-5000\t4000\t2\t0.5\n"
        "deletion\tchr3:1000-5000\t4000\t1\t0.01\n"
    )
    segs = C.load_cnvpytor(str(p))
    cns = sorted(s.cn for s in segs)
    assert 2 not in cns          # copy-neutral dropped
    assert cns == [1, 3]
