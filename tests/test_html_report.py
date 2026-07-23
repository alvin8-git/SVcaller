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


# --- alpha-globin card wiring (M8 -> report) ----------------------------------

def _minimal_inputs(tmp_path):
    smn_html = tmp_path / "smn.html"; smn_html.write_text("<p>SMN stub</p>")
    cnv_bed = tmp_path / "cnv.bed"
    cnv_bed.write_text("#chrom\tstart\tend\tcn\tsvtype\tcaller_support\tconfidence\tquality\tsample\n")
    sv_tsv = tmp_path / "sv.tsv"
    sv_tsv.write_text("AnnotSV_ID\tSV_chrom\tSV_start\tSV_end\tSV_type\tAnnotSV_ranking_score\n")
    circos_svg = tmp_path / "circos.svg"; circos_svg.write_text("<svg/>")
    return dict(smn_html_path=str(smn_html), cnv_bed_path=str(cnv_bed),
                sv_tsv_path=str(sv_tsv), circos_svg_path=str(circos_svg))


def test_alpha_card_is_rendered_when_supplied(tmp_path):
    """The whole point of wiring: an alpha-globin fragment must reach the report.
    For months the card generator (hba_report.py) existed but was never connected,
    so a --SEA carrier's result was invisible in the HTML."""
    from html_report import render_report
    alpha = tmp_path / "alpha.html"
    alpha.write_text('<div class="card"><h5>Alpha-globin (HBA1/HBA2)</h5>'
                     '<td>--SEA|--MED/aa</td></div>')
    out = tmp_path / "report.html"
    render_report(sample_id="AG_TEST", out_path=str(out), pipeline_version="1.0.0",
                  alpha_html_path=str(alpha), **_minimal_inputs(tmp_path))
    content = out.read_text()
    assert "Alpha-globin (HBA1/HBA2)" in content, "alpha card did not reach the report"
    assert "--SEA|--MED/aa" in content


def test_alpha_section_absent_when_not_supplied(tmp_path):
    """No fragment -> the {% if alpha_html %} guard drops the section, rather than
    emitting an empty card that reads as 'nothing found at this locus'."""
    from html_report import render_report
    out = tmp_path / "report.html"
    render_report(sample_id="NO_AG", out_path=str(out), pipeline_version="1.0.0",
                  alpha_html_path=None, **_minimal_inputs(tmp_path))
    assert "Alpha-globin (HBA1/HBA2)" not in out.read_text()


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


# --- Silent-failure guard: a 0-byte required input is a crash, not "no findings" ---
# A header-only file is a legitimate negative result (see the fallback test above);
# a truly EMPTY (0-byte / whitespace) required input is a crashed upstream stage and
# must raise UpstreamEmptyError instead of rendering a clean "no findings" report.

def _guard_inputs(tmp_path):
    smn_html = tmp_path / "smn.html"; smn_html.write_text("<p>SMN</p>")
    cnv_bed = tmp_path / "cnv.bed"
    cnv_bed.write_text("#chrom\tstart\tend\tcn\tsvtype\tcaller_support\tconfidence\tquality\tsample\n")
    sv_tsv = tmp_path / "sv.tsv"
    sv_tsv.write_text("AnnotSV_ID\tSV_chrom\tSV_start\tSV_end\tSV_type\tAnnotSV_ranking_score\n")
    circos_svg = tmp_path / "circos.svg"; circos_svg.write_text("<svg/>")
    return smn_html, cnv_bed, sv_tsv, circos_svg


@pytest.mark.parametrize("victim", ["sv_tsv", "smn_html", "circos_svg", "cnv_bed"])
def test_empty_required_input_raises(tmp_path, victim):
    from html_report import render_report, UpstreamEmptyError
    smn_html, cnv_bed, sv_tsv, circos_svg = _guard_inputs(tmp_path)
    paths = {"smn_html": smn_html, "cnv_bed": cnv_bed, "sv_tsv": sv_tsv, "circos_svg": circos_svg}
    paths[victim].write_text("")   # crashed upstream stage -> 0-byte placeholder
    out = tmp_path / "GUARD_TEST.report.html"
    with pytest.raises(UpstreamEmptyError):
        render_report(
            sample_id="GUARD_TEST",
            smn_html_path=str(smn_html), cnv_bed_path=str(cnv_bed),
            sv_tsv_path=str(sv_tsv), circos_svg_path=str(circos_svg),
            out_path=str(out), pipeline_version="1.0.0",
        )
    assert not out.exists()         # no clean-looking report left behind


def test_absent_cnv_sentinel_is_not_a_failure(tmp_path):
    """A NO_FILE sentinel (sample legitimately has no CNV data) must NOT raise."""
    from html_report import render_report
    smn_html, cnv_bed, sv_tsv, circos_svg = _guard_inputs(tmp_path)
    sentinel = tmp_path / "NO_FILE"; sentinel.write_text("")   # sentinel is empty by design
    out = tmp_path / "SENTINEL_TEST.report.html"
    render_report(
        sample_id="SENTINEL_TEST",
        smn_html_path=str(smn_html), cnv_bed_path=str(sentinel),
        sv_tsv_path=str(sv_tsv), circos_svg_path=str(circos_svg),
        out_path=str(out), pipeline_version="1.0.0",
    )
    assert out.exists() and "SENTINEL_TEST" in out.read_text()


# ---- SUPP_VEC caller attribution must not name a caller it cannot identify ----
# JASMINE_MERGE writes the core three unconditionally but appends Scramble and
# MELT ONLY when each has calls. The old code indexed a fixed six-name list by bit
# position, so a sample with no Scramble calls but MELT calls (a 4-bit vector) had
# MELT reported as "Scramble" — the report named a caller that contributed nothing
# to that variant. Found 2026-07-22 while wiring SvABA in.
def test_supp_vec_never_misattributes_optional_callers():
    from html_report import _parse_supp_vec, _supp_vec_names

    # positions 1-3 are guaranteed by merge.nf's unconditional printf
    assert _supp_vec_names(3) == ["Manta", "Delly", "GRIDSS"]

    # 4 bits: exactly one optional caller merged, and SUPP_VEC cannot say which
    four = _supp_vec_names(4)
    assert four[:3] == ["Manta", "Delly", "GRIDSS"]
    assert four[3] == "Scramble or MELT", (
        f"a 4-bit vector names {four[3]!r}; it is unknowable which optional "
        "caller occupies that bit, so it must not be claimed")

    # 5 bits: both present, so the order IS unambiguous
    assert _supp_vec_names(5) == ["Manta", "Delly", "GRIDSS", "Scramble", "MELT"]

    # the specific regression: bit 4 alone must not read as "Scramble"
    got = _parse_supp_vec("SUPP=1;SUPP_VEC=0001")["callers"]
    assert got != "Scramble", "the original mislabel is back"
    assert got == "Scramble or MELT"

    # a wider vector degrades positionally rather than reusing a name
    assert _supp_vec_names(6)[5].startswith("caller")

    # degenerate input stays honest
    assert _parse_supp_vec("SUPP=0;SUPP_VEC=000")["callers"] == "?"
    assert _parse_supp_vec("")["callers"] == "?"
