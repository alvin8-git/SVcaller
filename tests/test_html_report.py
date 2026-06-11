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


# --- SV merged-VCF fallback (2026-06 empty-SV-sheet regression) ---------------

_MERGED_VCF = """\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t1000\tsv1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=5000;SVLEN=-4000;SUPP=2;SUPP_VEC=110
chr2\t2000\tsv2\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP;END=3000;SVLEN=1000;SUPP=1;SUPP_VEC=100
chr3\t10\tsv3\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=900000000;SVLEN=900000000;SUPP=1
"""

def _merged(tmp_path):
    p = tmp_path / "merged.vcf"; p.write_text(_MERGED_VCF); return str(p)

def test_parse_sv_from_vcf_basic(tmp_path):
    from html_report import parse_sv_from_vcf
    rows = parse_sv_from_vcf(_merged(tmp_path))
    assert len(rows) == 2                       # sv3 exceeds max artifact size, filtered
    by_type = {r["svtype"]: r for r in rows}
    assert set(by_type) == {"DEL", "DUP"}
    assert by_type["DEL"]["size_bp"] == 4000
    assert by_type["DEL"]["supp_n"] == 2
    assert rows[0]["tier"] == 3 and rows[0]["gene"] == "—"

def test_parse_sv_from_vcf_missing_file():
    from html_report import parse_sv_from_vcf
    assert parse_sv_from_vcf("/no/such/file.vcf") == []

def test_sv_sheet_falls_back_to_vcf_when_tsv_empty(tmp_path):
    """End-to-end: empty AnnotSV TSV + merged VCF -> SV sheet populated, banner shown."""
    from html_report import render_report
    import openpyxl
    smn_html = tmp_path / "smn.html"; smn_html.write_text("<p>SMN</p>")
    cnv_bed = tmp_path / "cnv.bed"
    cnv_bed.write_text("#chrom\tstart\tend\tcn\tsvtype\tcaller_support\tconfidence\tquality\tsample\n")
    sv_tsv = tmp_path / "sv.tsv"      # header only -> empty -> triggers fallback
    sv_tsv.write_text("AnnotSV_ID\tSV_chrom\tSV_start\tSV_end\tSV_type\tAnnotSV_ranking_score\n")
    circos_svg = tmp_path / "circos.svg"; circos_svg.write_text("<svg/>")
    out = tmp_path / "FALLBACK_TEST.report.html"   # name must end .report.html for xlsx derivation
    render_report(
        sample_id="FALLBACK_TEST",
        smn_html_path=str(smn_html), cnv_bed_path=str(cnv_bed),
        sv_tsv_path=str(sv_tsv), circos_svg_path=str(circos_svg),
        out_path=str(out), pipeline_version="1.0.0",
        sv_vcf_path=_merged(tmp_path),
    )
    html = out.read_text()
    assert "Annotation unavailable" in html        # banner rendered
    xls = openpyxl.load_workbook(str(out).replace(".report.html", ".variants.xlsx"))
    sv_sheet = xls["SV"]
    assert sv_sheet.max_row - 1 == 2               # 2 data rows from the merged VCF
