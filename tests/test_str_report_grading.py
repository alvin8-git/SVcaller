"""html_report grades STR calls against the table of the catalogue that produced them (version gate), and shows only the
graded block of a compound locus.  Fixtures are real ExpansionHunter VCF records for GIAB HG002 (old and STRchive catalogue)."""
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import pytest
import html_report as hr

FX = os.path.join(REPO, "tests", "fixtures")
NEW = os.path.join(FX, "hg002_eh_new_catalogue.str.vcf")
OLD = os.path.join(FX, "hg002_eh_old_catalogue.str.vcf")
graded = hr._load_str_graded()


def grade(vcf):
    cat = hr.choose_str_catalogue(vcf, graded)
    rng = hr._load_str_ranges(hr.STR_RANGES_PATH if cat == "new" else hr.STR_RANGES_LEGACY_PATH)
    return cat, {l["locus"]: l for l in hr.parse_str_loci(vcf, rng, graded if cat == "new" else None)}


def test_gate_picks_the_table_of_the_catalogue_that_made_the_vcf():
    assert grade(NEW)[0] == "new" and grade(OLD)[0] == "legacy"


def test_new_catalogue_shows_only_the_graded_block_and_grades_it():
    _, loci = grade(NEW)
    assert sorted(loci) == ["CNBP", "FMR1", "HTT"]            # no DM2_CNBP_CA / HD_HTT_CCG rows, no UNKNOWN
    assert loci["CNBP"]["repcn"] == "15/16" and loci["CNBP"]["status"] == "NORMAL" and loci["CNBP"]["normal_max"] == "26"
    assert "DM2_CNBP_CAGA 7/9" in loci["CNBP"]["other_blocks"] and "DM2_CNBP_CA 22/26" in loci["CNBP"]["other_blocks"]
    assert loci["HTT"]["repcn"] == "17/24" and loci["HTT"]["status"] == "NORMAL" and "HD_HTT_CCG" in loci["HTT"]["other_blocks"]
    assert loci["FMR1"]["status"] == "NORMAL"
    assert all(l["status"] != "UNKNOWN" for l in loci.values())


def test_old_catalogue_output_is_graded_with_the_legacy_table():
    _, loci = grade(OLD)
    assert loci["CNBP"]["normal_max"] == "30" and loci["CNBP"]["status"] == "INTERMEDIATE"   # 13/51 whole-tract CCTG, legacy limits (the old artefact)
    assert loci["HTT"]["normal_max"] == "35" and loci["HTT"]["status"] == "NORMAL"
    assert loci["FMR1"]["normal_max"] == "44" and loci["FMR1"]["status"] == "NORMAL"


def test_premutation_label_is_restored_where_strchive_defines_one():
    rng = hr._load_str_ranges(hr.STR_RANGES_PATH)
    st = lambda loc, cn: hr._str_status(str(cn), rng[loc])
    assert st("FMR1", 120) == "PREMUTATION" and st("FMR1", 50) == "INTERMEDIATE" and st("FMR1", 201) == "EXPANDED"
    assert st("DMPK", 40) == "PREMUTATION" and st("CNBP", 45) == "PREMUTATION" and st("CNBP", 28) == "INTERMEDIATE"
    assert st("HTT", 30) == "PREMUTATION" and st("ATXN7", 34) == "INTERMEDIATE"
    assert st("TBP", 45) == "INTERMEDIATE"                    # STRchive: reduced penetrance, not a premutation


def test_mixed_catalogues_abort(tmp_path):
    mixed = tmp_path / "m.vcf"
    mixed.write_text(open(NEW).read() + "".join(l for l in open(OLD) if "VARID=HTT;" in l))
    # a record whose VARID is not a new graded variant is ignored; make one that IS graded but on the old region
    bad = [l for l in open(NEW) if "VARID=FMR1;" in l][0].replace("147912049", "147912050")
    mixed.write_text(open(NEW).read() + bad)
    with pytest.raises(SystemExit):
        hr.choose_str_catalogue(str(mixed), graded)
