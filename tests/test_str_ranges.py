"""assets/str_disease_ranges.tsv is generated from STRchive (bin/make_str_ranges.py); every row must equal the source."""
import csv
import importlib.util
import json
import os
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("mkr", os.path.join(REPO, "bin", "make_str_ranges.py"))
mkr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mkr)
rows = list(csv.DictReader(open(mkr.OUT), delimiter="\t"))
loci = {x["id"]: x for x in json.load(open(mkr.LOCI))}
stranger = {x["LocusId"]: x for x in json.load(open(os.path.join(REPO, "assets", "sources",
                                                                 "STRchive-disease-loci.hg38.stranger.json")))}


def test_committed_ranges_equal_the_generator_output():
    assert subprocess.run([sys.executable, os.path.join(REPO, "bin", "make_str_ranges.py"), "--check"]).returncode == 0


def test_every_row_equals_strchive_and_covers_the_catalog():
    cat = [e["LocusId"] for e in json.load(open(os.path.join(REPO, "assets", "eh_catalog.json")))]
    assert [r["locus_id"] for r in rows] == cat
    for r in rows:
        x = loci[r["strchive_id"]]
        assert r["path_min"] == str(x["pathogenic_min"]) == str(stranger[r["strchive_id"]]["PathologicMin"])
        assert r["normal_max"] == ("" if x["benign_max"] is None else str(x["benign_max"]))
        if x["benign_max"] is not None:
            assert x["benign_max"] == stranger[r["strchive_id"]]["NormalMax"]
        if x["intermediate_min"] is None:
            assert r["intermediate_min"] == ""
        else:       # STRchive's range, clamped (see test_intermediate_is_clamped...)
            assert int(r["intermediate_min"]) >= x["intermediate_min"] and int(r["intermediate_max"]) <= x["intermediate_max"]
        assert "STRchive v2.28.0" in r["source"] and r["strchive_id"] in r["source"]


def test_blank_means_strchive_gives_none_and_says_so():
    for r in rows:
        if not r["normal_max"]:
            assert "no benign/normal maximum" in r["notes"]
        if not r["intermediate_min"]:
            assert "no intermediate range" in r["notes"]
    assert {r["locus_id"] for r in rows if not r["normal_max"]} == {"BEAN1", "MARCHF6", "SAMD12", "STARD7"}


def test_disease_labels_come_from_strchive():
    by = {r["locus_id"]: r for r in rows}
    assert "Oculopharyngodistal myopathy type 4" in by["RILPL1"]["disease"]
    assert "myoclonic epilepsy 2" in by["STARD7"]["disease"]
    assert by["CNBP"]["graded_variant"] == "DM2_CNBP" and by["RFC1"]["graded_variant"] == "CANVAS_RFC1"


def test_intermediate_is_clamped_strictly_between_normal_max_and_path_min():
    for r in rows:
        if not r["intermediate_min"]:
            continue
        lo, hi = int(r["intermediate_min"]), int(r["intermediate_max"])
        if r["normal_max"]:
            assert lo > int(r["normal_max"]), r["locus_id"]
        assert hi < int(r["path_min"]), r["locus_id"]
    by = {r["locus_id"]: r for r in rows}
    assert (by["ATXN8OS"]["intermediate_min"], by["RFC1"]["intermediate_min"], by["ATXN10"]["intermediate_max"]) == ("51", "12", "799")


def test_premutation_only_where_strchive_text_says_so():
    details = {x["id"]: x["details"] for x in loci.values()}
    for r in rows:
        if r["premut_max"]:
            assert r["locus_id"] in mkr.PREMUTATION
            assert "premutation" in details[r["strchive_id"]].lower() or r["locus_id"] == "FMR1", r["locus_id"]
            assert "premutation range per STRchive text" in r["notes"]
            assert int(r["premut_max"]) < int(r["path_min"]) + 1
    by = {r["locus_id"]: r for r in rows}
    assert "gray zone" in details["FXS_FMR1"] and "premutation range" in details["FXS_FMR1"]
    assert (by["FMR1"]["premut_min"], by["FMR1"]["premut_max"]) == ("55", "200")
    assert (by["CNBP"]["premut_min"], by["CNBP"]["premut_max"]) == ("30", "74")
    assert not by["TBP"]["premut_max"]                  # STRchive: reduced penetrance, not a premutation


def test_repeat_unit_is_the_motif_eh_counts_and_inheritance_codes_are_the_old_ones():
    by = {r["locus_id"]: r for r in rows}
    assert by["CNBP"]["repeat_unit"] == "CAGG" and by["CNBP"]["disease_motif"] == "CCTG"
    assert by["RFC1"]["repeat_unit"] == "AAGGG" and by["FGF14"]["repeat_unit"] == "AAG"
    cat = {e["LocusId"]: e for e in json.load(open(os.path.join(REPO, "assets", "eh_catalog.json")))}
    for r in rows:
        assert r["repeat_unit"] == mkr.counted_motif(cat[r["locus_id"]], r["graded_variant"])
    assert by["FMR1"]["inheritance"] == "XLD" and by["AFF2"]["inheritance"] == "XLR" and by["AR"]["inheritance"] == "XLR"
    assert not any(r["inheritance"] in ("XD", "XR") for r in rows)


def test_legacy_table_is_the_previous_one_unchanged():
    import hashlib
    txt = open(os.path.join(REPO, "assets", "str_disease_ranges.legacy.tsv"), "rb").read()
    assert hashlib.sha256(txt).hexdigest() == "b188d8fcec1dded4a8091f02c15ae2d85ed177edb5cc09998547a7c40fcd8a14"
    assert txt.startswith(b"locus_id\tgene\tdisease\tinheritance\trepeat_unit\tnormal_max\tpremut_max\tpath_min\tnotes\n")
