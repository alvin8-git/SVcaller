import sys, json, tempfile
from pathlib import Path
sys.path.insert(0, "/data/alvin/SVcaller/bin")
import pytest

def test_report_contains_sample_id(tmp_path):
    from html_report import render_report
    smn_html = tmp_path / "smn.html"
    smn_html.write_text("<p>SMN stub</p>")
    cnv_bed = tmp_path / "cnv.bed"
    cnv_bed.write_text("#chrom\tstart\tend\tcn\tsvtype\tcaller_support\tconfidence\tquality\tsample\n")
    sv_tsv = tmp_path / "sv.tsv"
    sv_tsv.write_text("AnnotSV_ID\tSV_chrom\tSV_start\tSV_end\tSV_type\tAnnotSV_ranking_score\n")
    circos_svg = tmp_path / "circos.svg"
    circos_svg.write_text("<svg/>")
    out = tmp_path / "report.html"
    render_report(
        sample_id="HG002_TEST",
        smn_html_path=str(smn_html),
        cnv_bed_path=str(cnv_bed),
        sv_tsv_path=str(sv_tsv),
        circos_svg_path=str(circos_svg),
        out_path=str(out),
        pipeline_version="1.0.0",
    )
    content = out.read_text()
    assert "HG002_TEST" in content
    assert "<svg" in content
    assert "SVcaller" in content
