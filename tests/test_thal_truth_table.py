"""Consume validation/thal_truth_table.tsv.

`smn_truth_table.tsv` is read by nothing — not bin/, not tests/, not validation/
— so SMN accuracy rests on a human eyeballing TSVs and a stale coordinate there
would never be noticed. The plan calls that out explicitly and says not to repeat
it. This is the not-repeating.

Two jobs:
  1. schema and enum validation, so the table stays machine-readable;
  2. COUPLING — every coordinate the table cites as evidence is checked against
     the generated panels in assets/. Regenerate a panel and change a position,
     and this fails rather than leaving the truth table quietly wrong.
"""
import os
import re

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TRUTH = os.path.join(REPO, "validation", "thal_truth_table.tsv")

LOCI = {"alpha", "beta"}
ZYG = {"het", "hom", "compound_het", "none"}
CONF = {"high", "medium", "low"}
METHODS = {"depth", "vcf", "pileup+vcf", "reads", "panel_genotype"}


def _rows(path):
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if not l.startswith("#") and l.strip()]
    header = lines[0].split("\t")
    return header, [dict(zip(header, l.split("\t"))) for l in lines[1:]]


def _panel(name):
    with open(os.path.join(REPO, "assets", name)) as fh:
        lines = [l.rstrip("\n") for l in fh if not l.startswith("#") and l.strip()]
    header = lines[0].split("\t")
    return [dict(zip(header, l.split("\t"))) for l in lines[1:]]


@pytest.fixture(scope="module")
def truth():
    assert os.path.exists(TRUTH), "validation/thal_truth_table.tsv is missing"
    return _rows(TRUTH)


def test_schema(truth):
    header, rows = truth
    expected = ["sample", "locus", "genotype", "alpha_genes", "zygosity", "panel_ref",
                "method", "evidence", "orthogonal", "confidence", "note"]
    assert header == expected, f"header drift: {header}"
    assert rows, "truth table has no rows"
    for r in rows:
        assert r["locus"] in LOCI, f"{r['sample']}: bad locus {r['locus']}"
        assert r["zygosity"] in ZYG, f"{r['sample']}: bad zygosity {r['zygosity']}"
        assert r["confidence"] in CONF, f"{r['sample']}: bad confidence {r['confidence']}"
        assert r["method"] in METHODS, f"{r['sample']}: bad method {r['method']}"
        assert r["evidence"].strip(), f"{r['sample']}: evidence must be stated"
        if r["locus"] == "alpha":
            assert r["alpha_genes"].isdigit(), f"{r['sample']}: alpha row needs a gene count"
            assert 0 <= int(r["alpha_genes"]) <= 4
        else:
            assert r["alpha_genes"] == "NA", f"{r['sample']}: beta row must not claim a gene count"


def test_nothing_claims_orthogonal_confirmation(truth):
    """Every genotype here came from the same data a caller would be scored on.
    If this ever flips to `yes`, the note must say which assay - otherwise the
    table starts asserting an accuracy it cannot support."""
    _, rows = truth
    for r in rows:
        assert r["orthogonal"] in ("yes", "no")
        if r["orthogonal"] == "yes":
            assert re.search(r"gap-PCR|MLPA|Sanger|ddPCR", r["note"], re.I), \
                f"{r['sample']}: claims orthogonal confirmation without naming the assay"


def test_panel_refs_resolve(truth):
    """Every panel_ref must exist in a committed panel — this is the coupling."""
    _, rows = truth
    hba = {r["allele"] for r in _panel("hba_pathogenic_sites.tsv")}
    hbb = {r["allele"] for r in _panel("hbb_pathogenic_sites.tsv")}
    known = hba | hbb
    for r in rows:
        ref = r["panel_ref"]
        if ref == "-":
            continue
        assert ref in known, f"{r['sample']}: panel_ref {ref!r} is in no assets/ panel"


@pytest.mark.parametrize("sample,locus,allele,panel,coord", [
    ("THAL2", "alpha", "Hb Quong Sze", "hba_pathogenic_sites.tsv", "173548"),
    ("THAL1", "beta", "IVS-II-654", "hbb_pathogenic_sites.tsv", "5225923"),
    ("HG02379", "beta", "HbE", "hbb_pathogenic_sites.tsv", "5226943"),
    ("HG00583", "beta", "CD41-42", "hbb_pathogenic_sites.tsv", "5226762"),
])
def test_cited_evidence_matches_panel_coordinate(truth, sample, locus, allele, panel, coord):
    """The coordinate quoted in `evidence` must agree with the panel's derived
    position. A regenerated panel that moved a site would otherwise leave this
    table asserting a genotype at a position the caller no longer looks at.

    SNVs must match exactly. INDELS MUST NOT be compared naively: the panel
    records the first affected CODING base, while a VCF anchors one base away on
    the genome per its own left-alignment convention. Both are right, and
    assert-equal here would be asserting that two different conventions are the
    same convention. For indels this checks proximity instead, which still
    catches a panel regenerating to a different site.
    """
    _, rows = truth
    row = next((r for r in rows if r["sample"] == sample and r["locus"] == locus
                and r["panel_ref"] == allele), None)
    assert row, f"no {locus} row for {sample} with panel_ref {allele}"
    assert coord in row["evidence"], f"{sample}: evidence does not cite {coord}"
    site = next((p for p in _panel(panel) if p["allele"] == allele), None)
    assert site, f"{allele} missing from {panel}"

    if site["event"] == "snv":
        assert site["pos"] == coord, \
            f"{allele}: panel says {site['pos']}, truth table cites {coord}"
    else:
        delta = abs(int(site["pos"]) - int(coord))
        assert delta <= 10, (
            f"{allele}: panel {site['pos']} vs cited {coord} differ by {delta} bp — "
            "more than an anchoring convention can explain")


def test_degenerate_groups_agree_with_allele_table(truth):
    """A genotype written `--SEA|--MED` must match a group the allele table
    actually declares degenerate. If channel 3 later resolves the group, both
    files have to change together."""
    _, rows = truth
    alleles = {r["allele"]: r for r in _panel("hba_deletion_alleles.tsv")}
    for r in rows:
        for hap in r["genotype"].split("/"):
            if "|" not in hap:
                continue
            members = hap.split("|")
            for m in members:
                assert m in alleles, f"{r['sample']}: unknown allele {m!r} in genotype"
                dist = alleles[m]["depth_distinguishable"]
                assert dist.startswith("no:"), \
                    f"{r['sample']}: {m} is written as degenerate but the allele table says {dist}"
                for other in members:
                    assert other in dist, \
                        f"{r['sample']}: group {hap} disagrees with declared group {dist}"


def test_giab_not_recorded_as_alpha_negative(truth):
    """The pipeline has never produced a chr16 call below 14.6 Mb. That is the
    caller's blindness, not a genotype. Recording GIAB as alpha-negative would
    repeat the exact false negative this project exists to fix."""
    _, rows = truth
    giab = {f"HG00{i}" for i in range(1, 8)}
    for r in rows:
        if r["sample"] in giab and r["locus"] == "alpha":
            pytest.fail(f"{r['sample']} recorded with an alpha genotype — "
                        "GIAB alpha status is unknown, not negative")


def test_unmeasured_alleles_are_declared():
    """The table must say what it does NOT cover; a truth table read as complete
    silently becomes a sensitivity claim."""
    text = open(TRUTH).read()
    assert "UNMEASURED" in text
    for allele in ("-a3.7", "--FIL", "--THAI", "anti-3.7", "Hb Constant Spring"):
        assert allele in text, f"{allele} should be declared unmeasured"
