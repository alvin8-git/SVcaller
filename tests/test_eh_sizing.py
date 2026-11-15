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


def _compound_profile():
    """CNBP as EH reports it with the STRchive catalog: three blocks, each counted in its own motif."""
    def v(gt, ru, region):
        return dict(variant(gt, "%s-%s/%s-%s" % (gt.split("/")[0], gt.split("/")[0], gt.split("/")[1], gt.split("/")[1]),
                            span="(%s, 5), (%s, 5)" % tuple(gt.split("/")), ru=ru), ReferenceRegion=region)
    return {"LocusResults": {"CNBP": {"ReadLength": 150, "Variants": {
        "DM2_CNBP": v("15/16", "CAGG", "chr3:129172576-129172656"),
        "DM2_CNBP_CAGA": v("7/9", "CAGA", "chr3:129172656-129172696"),
        "DM2_CNBP_CA": v("22/26", "CA", "chr3:129172696-129172734")}}}}


def test_compound_locus_reports_only_the_graded_block_in_its_own_motif(tmp_path):
    graded = es.load_graded(os.path.join(REPO, "assets", "eh_graded_variants.tsv"))
    rows = list(es.summarise(_compound_profile(), graded))
    assert len(rows) == 1                                   # one row per locus, contract with OmniGen unchanged
    r = rows[0]
    assert r["locus"] == "CNBP" and r["variant_id"] == "DM2_CNBP" and r["genotype"] == "15/16"
    assert r["repeat_unit"] == "CAGG" and r["disease_motif"] == "CAGG"
    assert "26" not in r["genotype"] and "22" not in r["genotype"]    # never the CA block or the whole tract


def test_compound_locus_without_a_graded_table_is_refused():
    try:
        list(es.summarise(_compound_profile(), {}))
    except SystemExit as ex:
        assert "not in the graded-variant table" in str(ex)
    else:
        raise AssertionError("a multi-block locus must not be guessed")


def test_stale_profile_missing_the_graded_variant_is_refused():
    p = {"LocusResults": {"CNBP": {"ReadLength": 150, "Variants": {"CNBP": variant("13/46", "20-20/46-90", ru="CCTG")}}}}
    graded = es.load_graded(os.path.join(REPO, "assets", "eh_graded_variants.tsv"))
    try:
        list(es.summarise(p, graded))
    except SystemExit as ex:
        assert "different catalog" in str(ex)
    else:
        raise AssertionError("old-catalog profile with the new graded table must fail loudly")


# ---- PABPN1 count offset (owner decision 2026-11-15) -----------------------------------------------------------
OFFSETS = os.path.join(REPO, "assets", "eh_count_offsets.tsv")
FASTA = os.environ.get("REF_ROOT", "/path/to/ref") + "/GRCh38/hg38.canonical.fa"


def test_offsets_file_lists_pabpn1_only_with_minus_one_and_evidence():
    off = es.load_offsets(OFFSETS)
    assert off == {("PABPN1", "PABPN1"): -1}


def _fetch(chrom, start0, end):
    fai = {l.split("\t")[0]: l.split("\t")[1:5] for l in open(FASTA + ".fai")}
    ln, off, bases, width = map(int, fai[chrom])
    def pos(i):
        return off + (i // bases) * width + i % bases
    with open(FASTA, "rb") as fh:
        fh.seek(pos(start0))
        raw = fh.read(pos(end - 1) - pos(start0) + 1)
    return raw.replace(b"\n", b"").decode().upper()


def test_grch38_reference_pabpn1_ten_codons():
    if not os.path.exists(FASTA):
        import pytest
        pytest.skip("GRCh38 fasta absent")
    import json as _j
    cat = _j.load(open(os.path.join(REPO, "assets", "sources", "STRchive-loci.json")))
    rec = [r for r in cat if r.get("id") == "OPMD_PABPN1"][0]
    assert rec["chrom"] == "chr14" and (rec["start_hg38"], rec["stop_hg38"]) == (23321472, 23321503)
    st = rec["start_hg38"]                    # 0-based; the 31 bp interval = 10 codons + the first base of the next GGG
    seq = _fetch("chr14", st - 3, st + 30 + 9)    # ATG + 10 codons + GGG GCT GCG
    assert seq[:3] == "ATG" and seq.endswith("GGGGCTGCG"), seq
    tract = seq[3:-9]
    assert len(tract) % 3 == 0 and all(tract[i:i + 3].startswith("GC") for i in range(0, len(tract), 3))
    assert len(tract) // 3 == 10


def test_pabpn1_eleven_eleven_becomes_ten_ten_normal_and_thirteen_stays_pathogenic_range():
    off = {("PABPN1", "PABPN1"): -1}
    graded = {"PABPN1": {"variant_id": "PABPN1", "disease_motif": "GCN"}}
    def run(gt, ci, span):
        p = profile(PABPN1=variant(gt, ci, span=span, ru="GCN"))
        return next(es.summarise(p, graded, off))
    r = run("11/11", "11-11/11-11", "(10, 1), (11, 22)")
    assert r["genotype"] == "10/10" and r["genotype_display"] == "10/10" and r["ci"] == "10-10/10-10"
    assert r["eh_raw_genotype"] == "11/11" and r["count_offset"] == -1 and r["allele1_est"] == 10
    r = run("11/13", "11-11/13-13", "(11, 5), (13, 4)")
    assert r["genotype"] == "10/12" and r["allele2_sizing"] == "spanned"          # 12 = STRchive pathogenic minimum
    ranges = {l.split("\t")[0]: l.split("\t") for l in open(os.path.join(REPO, "assets", "str_disease_ranges.tsv"))
              if l.strip() and not l.startswith("locus_id")}
    assert ranges["PABPN1"][5] == "10" and ranges["PABPN1"][7] == "12"
    assert int(r["allele2_est"]) >= int(ranges["PABPN1"][7])


def test_no_offset_changes_other_loci_and_raw_is_kept():
    graded = {"NIPA1": {"variant_id": "NIPA1", "disease_motif": "GCG"}}
    r = next(es.summarise(profile(NIPA1=variant("9/9", "9-9/9-9", span="(9, 7)", ru="GCG")), graded,
                          {("PABPN1", "PABPN1"): -1}))
    assert r["genotype"] == "9/9" and r["eh_raw_genotype"] == "9/9" and r["count_offset"] == 0


def test_offset_without_an_eh_ci_keeps_the_empty_ci_text_not_none_none():
    graded = {"PABPN1": {"variant_id": "PABPN1", "disease_motif": "GCN"}}
    p = profile(PABPN1=variant("11/11", "", span="(10, 1), (11, 22)", ru="GCN"))
    r = next(es.summarise(p, graded, {("PABPN1", "PABPN1"): -1}))
    assert r["genotype"] == "10/10" and "None" not in r["ci"]
