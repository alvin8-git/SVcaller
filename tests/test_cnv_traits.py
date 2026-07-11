"""Deterministic unit tests for the CNV / blood-group trait interpreters.

All of these run in-session WITHOUT a pipeline run: interpreter arithmetic uses
synthetic mosdepth region fixtures, the consensus cross-check uses the real
results/HG002/HG002.cnv_consensus.bed, and HTML rendering uses fixture TSVs.
"""
import gzip
import sys
from pathlib import Path

import pytest

REPO = Path("/data/alvin/SVcaller")
sys.path.insert(0, str(REPO / "bin"))

import cnv_traits_common as c  # noqa: E402
import rh_status  # noqa: E402
import amy1_cn  # noqa: E402
import gst_null  # noqa: E402
import lpa_kiv2  # noqa: E402
import html_report  # noqa: E402

HG002_BED = REPO / "results" / "HG002" / "HG002.cnv_consensus.bed"


# --------------------------------------------------------------------------- #
# Fixtures
# --------------------------------------------------------------------------- #
def write_depth(tmp_path, region_means, control_depth=30.0, n_ctrl=8):
    """Write a mosdepth-style regions.bed.gz {label: mean}. Adds CTRL_* controls."""
    depths = dict(region_means)
    for i in range(1, n_ctrl + 1):
        depths.setdefault(f"CTRL_{i}", control_depth)
    coords = {
        "RHD": ("chr1", 25272393, 25330445),
        "AMY1_CLUSTER": ("chr1", 103655000, 103760000),
        "GSTM1": ("chr1", 109687814, 109693020),
        "GSTT1": ("chr22", 24376133, 24384680),
        "LPA_KIV2": ("chr6", 160605000, 160650000),
    }
    p = tmp_path / "s.trait_depth.regions.bed.gz"
    with gzip.open(p, "wt") as fh:
        for label, mean in depths.items():
            chrom, start, end = coords.get(label, ("chr9", 1000, 2000))
            fh.write(f"{chrom}\t{start}\t{end}\t{label}\t{mean}\n")
    return p


# --------------------------------------------------------------------------- #
# common helpers
# --------------------------------------------------------------------------- #
def test_read_region_depths_and_baseline(tmp_path):
    p = write_depth(tmp_path, {"RHD": 30.0}, control_depth=30.0)
    depths = c.read_region_depths(p)
    assert depths["RHD"] == 30.0
    assert c.control_baseline(depths) == 30.0
    assert c.estimate_copies(depths["RHD"], 30.0) == pytest.approx(2.0)


def test_control_baseline_is_median(tmp_path):
    # one wild control must not skew the baseline
    p = write_depth(tmp_path, {"RHD": 30.0}, control_depth=30.0)
    depths = c.read_region_depths(p)
    depths["CTRL_1"] = 300.0
    assert c.control_baseline(depths) == 30.0


def test_control_baseline_missing_raises():
    with pytest.raises(ValueError):
        c.control_baseline({"RHD": 30.0})


# --------------------------------------------------------------------------- #
# Rh / RHD
# --------------------------------------------------------------------------- #
def test_rh_positive_normal_depth(tmp_path):
    p = write_depth(tmp_path, {"RHD": 30.0}, control_depth=30.0)
    depths = c.read_region_depths(p)
    copies, status, conf = rh_status.call_rh(depths, [])
    assert copies == 2 and status == "pos" and conf == "HIGH"


def test_rh_negative_homozygous_deletion(tmp_path):
    p = write_depth(tmp_path, {"RHD": 0.2}, control_depth=30.0)
    depths = c.read_region_depths(p)
    copies, status, conf = rh_status.call_rh(depths, [])
    assert copies == 0 and status == "neg"
    # no consensus support -> MEDIUM
    assert conf == "MEDIUM"


def test_rh_negative_with_consensus_del_is_high(tmp_path):
    p = write_depth(tmp_path, {"RHD": 0.1}, control_depth=30.0)
    depths = c.read_region_depths(p)
    cnv_rows = [{"chrom": "chr1", "start": 25270000, "end": 25335000,
                 "cn": "0", "svtype": "DEL", "caller_support": "BOTH",
                 "confidence": "HIGH"}]
    copies, status, conf = rh_status.call_rh(depths, cnv_rows)
    assert status == "neg" and conf == "HIGH"


def test_rh_unknown_without_controls():
    copies, status, conf = rh_status.call_rh({"RHD": 30.0}, [])
    assert copies == "NA" and status == "unknown" and conf == "LOW"


# --------------------------------------------------------------------------- #
# AMY1
# --------------------------------------------------------------------------- #
def test_amy1_copy_number(tmp_path):
    # AMY1 depth 6x the diploid baseline -> 12 copies
    p = write_depth(tmp_path, {"AMY1_CLUSTER": 180.0}, control_depth=30.0)
    depths = c.read_region_depths(p)
    assert amy1_cn.call_amy1(depths) == 12


def test_amy1_na_without_region(tmp_path):
    p = write_depth(tmp_path, {}, control_depth=30.0)
    depths = c.read_region_depths(p)
    assert amy1_cn.call_amy1(depths) == "NA"


# --------------------------------------------------------------------------- #
# GSTM1 / GSTT1 null
# --------------------------------------------------------------------------- #
def test_gst_both_present(tmp_path):
    p = write_depth(tmp_path, {"GSTM1": 30.0, "GSTT1": 30.0}, control_depth=30.0)
    depths = c.read_region_depths(p)
    base = c.control_baseline(depths)
    assert gst_null.call_gene("GSTM1", depths, base, []) == "present"
    assert gst_null.call_gene("GSTT1", depths, base, []) == "present"


def test_gst_null_homozygous_deletion(tmp_path):
    p = write_depth(tmp_path, {"GSTM1": 1.0, "GSTT1": 30.0}, control_depth=30.0)
    depths = c.read_region_depths(p)
    base = c.control_baseline(depths)
    assert gst_null.call_gene("GSTM1", depths, base, []) == "null"     # ratio .033 < .15
    assert gst_null.call_gene("GSTT1", depths, base, []) == "present"


def test_gst_heterozygous_is_present(tmp_path):
    p = write_depth(tmp_path, {"GSTM1": 15.0}, control_depth=30.0)
    depths = c.read_region_depths(p)
    base = c.control_baseline(depths)
    assert gst_null.call_gene("GSTM1", depths, base, []) == "present"  # ratio .5


# --------------------------------------------------------------------------- #
# LPA KIV-2
# --------------------------------------------------------------------------- #
def test_lpa_kiv2_copies(tmp_path):
    # KIV-2 window at 15x the diploid baseline -> 30 diploid copies
    p = write_depth(tmp_path, {"LPA_KIV2": 450.0}, control_depth=30.0)
    depths = c.read_region_depths(p)
    assert lpa_kiv2.call_kiv2(depths) == 30


# --------------------------------------------------------------------------- #
# Consensus cross-check against the REAL HG002 consensus BED
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(not HG002_BED.exists(), reason="HG002 consensus bed absent")
def test_hg002_consensus_no_deletion_at_deletable_loci():
    rows = c.parse_cnv_bed(HG002_BED)
    assert rows, "expected consensus rows to parse"
    for label in ("RHD", "GSTM1", "GSTT1"):
        chrom, start, end = c.LOCI[label]
        dels = c.consensus_overlaps(rows, chrom, start, end, svtype="DEL")
        # Consensus under-calls these homozygous deletions -> depth must lead.
        assert dels == [], f"unexpected consensus DEL at {label}"


@pytest.mark.skipif(not HG002_BED.exists(), reason="HG002 consensus bed absent")
def test_hg002_rh_positive_from_consensus_and_depth(tmp_path):
    # HG002/NA24385 is Rh D positive: normal RHD depth + no consensus DEL -> pos
    rows = c.parse_cnv_bed(HG002_BED)
    p = write_depth(tmp_path, {"RHD": 30.0}, control_depth=30.0)
    depths = c.read_region_depths(p)
    _, status, _ = rh_status.call_rh(depths, rows)
    assert status == "pos"


@pytest.mark.skipif(not HG002_BED.exists(), reason="HG002 consensus bed absent")
def test_hg002_lpa_locus_has_overlap():
    rows = c.parse_cnv_bed(HG002_BED)
    chrom, start, end = c.LOCI["LPA_KIV2"]
    # LPA region overlaps at least one consensus segment in HG002
    assert c.consensus_overlaps(rows, chrom, 160531483, 160664275) != []


# --------------------------------------------------------------------------- #
# Contract-schema round-trips (end-to-end via each script's main writer)
# --------------------------------------------------------------------------- #
def test_contract_schemas(tmp_path):
    depth = write_depth(tmp_path, {
        "RHD": 30.0, "AMY1_CLUSTER": 90.0, "GSTM1": 0.5,
        "GSTT1": 30.0, "LPA_KIV2": 300.0}, control_depth=30.0)
    depths = c.read_region_depths(depth)

    rh_out = tmp_path / "s.rh_status.tsv"
    c.write_tsv(rh_out, ["sample", "RHD_copies", "Rh_status", "confidence"],
                ["S", *rh_status.call_rh(depths, [])])
    hdr = rh_out.read_text().splitlines()[0]
    assert hdr == "#sample\tRHD_copies\tRh_status\tconfidence"

    amy_out = tmp_path / "s.amy1.tsv"
    c.write_tsv(amy_out, ["sample", "AMY1_copies", "method"],
                ["S", amy1_cn.call_amy1(depths), amy1_cn.METHOD])
    assert amy_out.read_text().splitlines()[0] == "#sample\tAMY1_copies\tmethod"

    base = c.control_baseline(depths)
    gst_out = tmp_path / "s.gst_null.tsv"
    c.write_tsv(gst_out, ["sample", "GSTM1", "GSTT1"],
                ["S", gst_null.call_gene("GSTM1", depths, base, []),
                 gst_null.call_gene("GSTT1", depths, base, [])])
    assert gst_out.read_text().splitlines()[0] == "#sample\tGSTM1\tGSTT1"

    lpa_out = tmp_path / "s.lpa_kiv2.tsv"
    c.write_tsv(lpa_out, ["sample", "KIV2_copies", "method"],
                ["S", lpa_kiv2.call_kiv2(depths), lpa_kiv2.METHOD])
    assert lpa_out.read_text().splitlines()[0] == "#sample\tKIV2_copies\tmethod"


# --------------------------------------------------------------------------- #
# HTML rendering
# --------------------------------------------------------------------------- #
def _fixture_tsvs(tmp_path):
    (tmp_path / "rh.tsv").write_text(
        "#sample\tRHD_copies\tRh_status\tconfidence\nHG002\t2\tpos\tHIGH\n")
    (tmp_path / "amy1.tsv").write_text(
        "#sample\tAMY1_copies\tmethod\nHG002\t8\tread-depth-normalized\n")
    (tmp_path / "gst.tsv").write_text(
        "#sample\tGSTM1\tGSTT1\nHG002\tnull\tpresent\n")
    (tmp_path / "lpa.tsv").write_text(
        "#sample\tKIV2_copies\tmethod\nHG002\t26\tread-depth-ratio\n")
    return tmp_path


def test_build_cnv_traits_section_renders(tmp_path):
    d = _fixture_tsvs(tmp_path)
    html = html_report.build_cnv_traits_section(
        str(d / "rh.tsv"), str(d / "amy1.tsv"),
        str(d / "gst.tsv"), str(d / "lpa.tsv"))
    assert "Blood Group &amp; Copy-Number Traits" in html
    assert "Rh(D) positive" in html
    assert "8 copies" in html          # AMY1
    assert "GSTM1" in html and "null" in html
    assert "26 repeat copies" in html  # LPA KIV-2


def test_build_cnv_traits_section_empty_when_absent():
    assert html_report.build_cnv_traits_section(None, None, None, None) == ""
    assert html_report.build_cnv_traits_section("NO_FILE", "NO_FILE",
                                                "NO_FILE", "NO_FILE") == ""


def test_build_cnv_traits_section_partial(tmp_path):
    d = _fixture_tsvs(tmp_path)
    html = html_report.build_cnv_traits_section(str(d / "rh.tsv"), None, None, None)
    assert "Rh(D) positive" in html
    assert "AMY1 (salivary amylase)" not in html   # AMY1 row absent (caption still lists loci)
    assert "repeat copies" not in html              # LPA row absent
    assert "8 copies" not in html                   # AMY1 numeric row absent
