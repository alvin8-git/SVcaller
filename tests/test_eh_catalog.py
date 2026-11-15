"""assets/eh_catalog.json is generated from the pinned STRchive hg38 catalogs (assets/sources/, release v2.28.0,
commit e964fa88df3b2189ee478cbd6eaad15bb28e3a0a); every entry must equal its source, and EH_SIZING must grade the
disease-relevant block of a compound locus (assets/eh_graded_variants.tsv)."""
import hashlib
import importlib.util
import json
import os
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
A = os.path.join(REPO, "assets")
SRC = os.path.join(A, "sources")
spec = importlib.util.spec_from_file_location("make_eh_catalog", os.path.join(REPO, "bin", "make_eh_catalog.py"))
mk = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mk)

SHA256 = {   # upstream files at STRchive e964fa88df3b2189ee478cbd6eaad15bb28e3a0a (v2.28.0)
    "STRchive-disease-loci.hg38.expansionhunter.json": "e87177dcae3e09acaa3fe766818e16d334b1f206f29f714bcf6b42fd2f5a5c9a",
    "STRchive-disease-loci.hg38.stranger.json": "f0475cd8803f4802d592edd86e42716728d3410f6952fb3dedcd3d8906478877",
    "STRchive-loci.json": "16301ecb00239f0e3c684913aefc2f48c6f5a7c2a481958aab4db55d282bca86",
    # Illumina/ExpansionHunter variant_catalog/grch38/variant_catalog.json @ 7ff9aa5335fdb5b3124cc0f365a64bb7a185803a
    "Illumina-ExpansionHunter.grch38.variant_catalog.json": "41043d6cac8a535922c3e29960991e09aaad3b24b14595a7d864232c455781a9",
}
OVERRIDES = {"AR", "AFF2"}      # owner decision 2026-11-15: regions from the Illumina EH catalog
cat = json.load(open(os.path.join(A, "eh_catalog.json")))
eh_src = {e["LocusId"]: e for e in json.load(open(mk.EH_SRC))}
LOCI_32 = ["FMR1", "AFF2", "AR", "ATN1", "ATXN1", "ATXN2", "ATXN3", "ATXN7", "ATXN8OS", "ATXN10", "CACNA1A", "CNBP", "C9ORF72",
           "DMPK", "FXN", "HTT", "NIPA1", "NOTCH2NLC", "NOP56", "PABPN1", "PPP2R2B", "RFC1", "TBP", "BEAN1", "EIF4A3",
           "MARCHF6", "SAMD12", "RILPL1", "STARD7", "CSTB", "FGF14", "DAB1"]


def test_source_files_are_the_pinned_upstream_bytes():
    for name, digest in SHA256.items():
        assert hashlib.sha256(open(os.path.join(SRC, name), "rb").read()).hexdigest() == digest, name


def test_every_repo_entry_equals_its_strchive_source():
    m = {a: b for a, b, _ in mk.read_map()}
    assert [e["LocusId"] for e in cat] == LOCI_32            # LocusIds are the names OmniGen reads
    for e in cat:
        if e["LocusId"] in OVERRIDES:
            continue
        s = eh_src[m[e["LocusId"]]]
        assert e["LocusStructure"] == s["LocusStructure"], e["LocusId"]
        assert e["ReferenceRegion"] == s["ReferenceRegion"], e["LocusId"]
        assert e["VariantType"] == s["VariantType"], e["LocusId"]


def test_committed_outputs_equal_what_the_generator_makes():
    r = subprocess.run([sys.executable, os.path.join(REPO, "bin", "make_eh_catalog.py"), "--check"])
    assert r.returncode == 0


def test_compound_loci_carry_variant_ids_and_a_graded_block():
    graded = {r["locus"]: r for r in (dict(zip(mk.GRADED_COLS, l.rstrip("\n").split("\t")))
                                      for l in open(mk.OUT_GRADED) if not l.startswith(("#", "locus")))}
    assert set(graded) == set(LOCI_32)
    for e in cat:
        g = graded[e["LocusId"]]
        if isinstance(e["ReferenceRegion"], list):
            assert len(e["VariantId"]) == len(e["ReferenceRegion"]) == len(e["VariantType"]) == int(g["n_variants"])
            assert g["variant_id"] in e["VariantId"]
        else:
            assert g["variant_id"] == e["LocusId"]       # EH names a single variant after its locus
    assert graded["CNBP"]["variant_id"] == "DM2_CNBP" and graded["RFC1"]["variant_id"] == "CANVAS_RFC1"
    assert graded["RFC1"]["disease_motif"] == "AAGGG" and graded["HTT"]["variant_id"] == "HD_HTT"


def test_no_catalog_region_is_the_old_wrong_build_or_position():
    old_wrong = {"ATXN1": "chr6:16327864-16327916", "ATXN8OS": "chr13:70713515-70713570", "CACNA1A": "chr19:13206436-13206455",
                 "PPP2R2B": "chr5:146258289-146258327", "NIPA1": "chr15:22787428-22787497", "SAMD12": "chr8:119379048-119379081"}
    by = {e["LocusId"]: e for e in cat}
    for loc, region in old_wrong.items():
        assert by[loc]["ReferenceRegion"] != region


def test_ar_and_aff2_regions_are_the_illumina_ones_and_nothing_else_is_overridden():
    ill = {e["LocusId"]: e for e in json.load(open(mk.ILL_SRC))}
    over = {a for a, _, src in mk.read_map() if src}
    assert over == OVERRIDES and all(src == "illumina" for _, _, src in mk.read_map() if src)
    by = {e["LocusId"]: e for e in cat}
    for loc in OVERRIDES:
        assert by[loc]["ReferenceRegion"] == "chr" + ill[loc]["ReferenceRegion"]
        assert by[loc]["LocusStructure"] == ill[loc]["LocusStructure"]
    assert by["AR"]["ReferenceRegion"] == "chrX:67545316-67545385" and by["AFF2"]["ReferenceRegion"] == "chrX:148500631-148500691"
