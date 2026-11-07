"""ExpansionHunter per-allele sizing (bin/eh_sizing.py): spanned / flanked / unresolved.

EH prints a point estimate even for alleles longer than the read (HG002 ATXN2 22/62, CI 52-102;
Q100 truth is ~22/29).  Those must come out as '>=<CI lower>' with conf=unresolved, never 62.
"""
import importlib.util
import json
import os
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("eh_sizing", os.path.join(REPO, "bin", "eh_sizing.py"))
es = importlib.util.module_from_spec(spec)
spec.loader.exec_module(es)


def variant(gt, ci, span="()", flank="()", inrep="()", ru="CAG"):
    return {"Genotype": gt, "GenotypeConfidenceInterval": ci, "RepeatUnit": ru,
            "CountsOfSpanningReads": span, "CountsOfFlankingReads": flank, "CountsOfInrepeatReads": inrep}


def profile(**loci):
    return {"LocusResults": {k: {"ReadLength": 150, "Variants": {k: v}} for k, v in loci.items()}}


def test_spanned_needs_two_reads_at_that_size():
    r = es.classify_variant(variant("14/14", "11-14/14-42", span="(14, 4)"), 150)
    assert [a["sizing"] for a in r["alleles"]] == ["spanned", "spanned"]
    assert r["conf"] == "high" and r["genotype_display"] == "14/14"
    one = es.classify_variant(variant("10/20", "10-10/20-20", span="(10, 5), (20, 1)"), 150)
    assert [a["sizing"] for a in one["alleles"]] == ["spanned", "flanked"]   # 1 read is not enough


def test_flanked_short_allele_without_spanning_reads():
    r = es.classify_variant(variant("5/19", "5-5/19-22", flank="(5, 3), (19, 2)"), 150)
    assert [a["sizing"] for a in r["alleles"]] == ["flanked", "flanked"]     # 19*3=57 bp < 150
    assert r["conf"] in ("high", "low")


def test_long_allele_without_spanning_is_lower_bound_not_point_estimate():
    """HG002 ATXN2: 62*3 = 186 bp > 150 bp read, no spanning read."""
    r = es.classify_variant(variant("22/62", "20-22/52-102", span="(22, 6)", flank="(60, 2)", inrep="(62, 3)"), 150)
    assert [a["sizing"] for a in r["alleles"]] == ["spanned", "unresolved"]
    assert r["genotype_display"] == "22/>=52" and r["alleles"][1]["display"] == ">=52"
    assert r["conf"] == "unresolved"
    assert "62" not in r["genotype_display"]


def test_wide_ci_alone_is_unresolved_even_if_it_fits_the_read():
    r = es.classify_variant(variant("20/20", "20-20/20-60", ru="A"), 150)    # 20 bp, CI width 40 > 30
    assert r["alleles"][1]["sizing"] == "unresolved" and r["alleles"][1]["display"] == ">=20"


def test_spanning_support_beats_length_rule():
    r = es.classify_variant(variant("40/40", "40-40/40-40", span="(40, 3)", ru="GAA"), 150)  # 120 bp
    assert r["alleles"][0]["sizing"] == "spanned"


def test_read_length_comes_from_json_and_haploid_is_handled():
    p = profile(X=variant("60", "60-60", ru="CTG"))
    p["LocusResults"]["X"]["ReadLength"] = 250       # 180 bp fits a 250 bp read
    row = next(es.summarise(p))
    assert row["allele1_sizing"] == "flanked" and row["allele2_sizing"] == "." and row["read_length"] == 250


def test_cli_writes_tsv_and_leaves_json_untouched(tmp_path):
    p = profile(ATXN2=variant("22/62", "20-22/52-102", span="(22, 6)"), AR=variant("27/28", "27-28/28-28", span="(27, 3), (28, 3)"))
    src, out = tmp_path / "S.str_profile.json", tmp_path / "S.str_summary.tsv"
    src.write_text(json.dumps(p))
    before = src.read_bytes()
    subprocess.run([sys.executable, os.path.join(REPO, "bin", "eh_sizing.py"), str(src), str(out)], check=True)
    assert src.read_bytes() == before
    lines = out.read_text().splitlines()
    assert lines[0].startswith("# str_summary.tsv")
    hdr = [l for l in lines if not l.startswith("#")]
    cols = hdr[0].split("\t")
    for c in ("allele1_sizing", "allele2_sizing", "allele2_display", "genotype_display", "conf"):
        assert c in cols
    rows = {r.split("\t")[0]: dict(zip(cols, r.split("\t"))) for r in hdr[1:]}
    assert rows["ATXN2"]["conf"] == "unresolved" and rows["ATXN2"]["genotype_display"] == "22/>=52"
    assert rows["AR"]["conf"] == "high"


def test_empty_profile_fails_loudly(tmp_path):
    src = tmp_path / "e.json"
    src.write_text("{}")
    r = subprocess.run([sys.executable, os.path.join(REPO, "bin", "eh_sizing.py"), str(src), str(tmp_path / "o.tsv")])
    assert r.returncode != 0


def test_male_chrx_is_hemizygous_even_for_a_legacy_diploid_call():
    """HG002 was run without --sex (EH default female): FMR1 31/31, AR 27/28. Male chrX has one allele."""
    p = {"SampleParameters": {"Sex": "male"}, "LocusResults": {
        "FMR1": {"ReadLength": 150, "Variants": {"FMR1": dict(variant("31/31", "27-31/31-37", span="(31, 4)", ru="CGG"),
                                                              ReferenceRegion="chrX:147912050-147912110")}},
        "AR": {"ReadLength": 150, "Variants": {"AR": dict(variant("27/28", "27-28/28-28", span="(27, 3), (28, 3)"),
                                                          ReferenceRegion="chrX:67545315-67545384")}},
        "HTT": {"ReadLength": 150, "Variants": {"HTT": dict(variant("17/24", "17-17/24-24", span="(17, 3), (24, 3)"),
                                                            ReferenceRegion="chr4:3074876-3074932")}}}}
    rows = {r["locus"]: r for r in es.summarise(p)}
    assert rows["FMR1"]["allele1_est"] == 31 and rows["FMR1"]["allele2_est"] == "."
    assert rows["FMR1"]["genotype"] == "31" and rows["FMR1"]["genotype_display"] == "31"
    assert rows["AR"]["genotype"] == "28" and rows["AR"]["allele2_display"] == "."
    assert rows["HTT"]["genotype"] == "17/24"                    # autosome untouched
    p["SampleParameters"]["Sex"] = "female"
    assert next(r for r in es.summarise(p) if r["locus"] == "FMR1")["genotype"] == "31/31"


def _load_infer_sex():
    sp = importlib.util.spec_from_file_location("infer_sex", os.path.join(REPO, "bin", "infer_sex.py"))
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


def test_sex_rule_matches_omnigen_thresholds():
    inf = _load_infer_sex()
    assert inf.decide(0.012, 0.986) == "female"       # HG001
    assert inf.decide(0.556, 0.509) == "male"         # HG002
    assert inf.decide(0.10, 0.9) == "unknown"         # Y at the female ceiling
    assert inf.decide(0.30, 0.80) == "unknown"        # X too high for male
    assert inf.decide(0.05, 0.50) == "unknown"        # XO-like: neither rule


def test_cli_header_carries_sex_and_method(tmp_path):
    p = profile(AR=dict(variant("27/28", "27-28/28-28", span="(27, 3), (28, 3)"), ReferenceRegion="chrX:67545315-67545384"))
    p["SampleParameters"] = {"Sex": "male"}
    sj = tmp_path / "S.sex.json"
    sj.write_text(json.dumps({"inferred_sex": "male", "chrX_ratio": 0.5, "chrY_ratio": 0.55, "method": "depth windows", "sex_flag": ""}))
    src, out = tmp_path / "S.str_profile.json", tmp_path / "o.tsv"
    src.write_text(json.dumps(p))
    subprocess.run([sys.executable, os.path.join(REPO, "bin", "eh_sizing.py"), str(src), str(out), "--sex-json", str(sj)], check=True)
    txt = out.read_text()
    assert "# sex inferred: male" in txt and "method: depth windows" in txt and "SampleParameters.Sex): male" in txt
