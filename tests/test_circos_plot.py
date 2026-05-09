import sys, tempfile
from pathlib import Path
sys.path.insert(0, "/data/alvin/SVcaller/bin")
import pytest

def test_sv_colour_mapping():
    from circos_plot import sv_colour
    assert sv_colour("DEL") == "#1F77B4"
    assert sv_colour("DUP") == "#D62728"
    assert sv_colour("INV") == "#9467BD"
    assert sv_colour("BND") == "#FF7F0E"
    assert sv_colour("TRA") == "#FF7F0E"

def test_parse_cnv_bed_gains_and_losses(tmp_path):
    from circos_plot import parse_cnv_bed
    bed = tmp_path / "cnv.bed"
    bed.write_text(
        "#chrom\tstart\tend\tcn\tsvtype\tcaller_support\tconfidence\tquality\tsample\n"
        "chr1\t1000000\t5000000\t3\tDUP\tBOTH\tHIGH\t50\tS1\n"
        "chr2\t2000000\t4000000\t1\tDEL\tBOTH\tHIGH\t50\tS1\n"
    )
    gains, losses = parse_cnv_bed(str(bed))
    assert len(gains)  == 1
    assert len(losses) == 1
    assert gains[0]["chrom"]  == "chr1"
    assert losses[0]["chrom"] == "chr2"

def test_parse_sv_vcf_extracts_links(tmp_path):
    from circos_plot import parse_sv_vcf_links
    vcf = tmp_path / "sv.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t1000000\t.\tN\t<INV>\t.\tPASS\tSVTYPE=INV;END=2000000\n"
        "chr3\t5000000\t.\tN\tN[chr7:8000000[\t.\tPASS\tSVTYPE=BND;MATEID=.\n"
    )
    links = parse_sv_vcf_links(str(vcf))
    types = {l["svtype"] for l in links}
    assert "INV" in types
    assert "BND" in types
