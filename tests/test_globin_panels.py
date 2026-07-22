"""Guard the globin site panels and the SVcaller->OmniGen contract.

Two independent tracks build against these artifacts, so drift here surfaces as
a silent integration failure late. These tests are cheap and catch it early.

The gene-model tests need the AnnotSV bundle and are skipped without it; the
contract tests are pure file checks and always run.
"""
import os
import subprocess
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(REPO, "bin")
ANNOTSV_BED = "/data/alvin/ref/annotsv/Annotations_Human/Genes/GRCh38/genes.RefSeq.sorted.bed"
CONTRACT = os.path.join(REPO, "docs/contracts/alpha_globin_contract.md")
FIXTURE = os.path.join(REPO, "validation/examples/SAMPLE.alpha_globin.tsv")

needs_annotsv = pytest.mark.skipif(
    not os.path.exists(ANNOTSV_BED), reason="AnnotSV gene-model bundle not present")


@needs_annotsv
def test_hgvs_map_selftest():
    """Known coordinates, each established independently of the mapper."""
    r = subprocess.run([sys.executable, os.path.join(BIN, "hgvs_map.py"), "--selftest"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "FAIL" not in r.stdout


@needs_annotsv
@pytest.mark.parametrize("gene,cpos,expect", [
    ("HBB", "20", 5227002),        # HbS
    ("HBB", "79", 5226943),        # HbE
    ("HBA2", "377", 173548),       # Hb Quong Sze
    ("HBB", "316-197", 5225923),   # IVS-II-654, intronic
    ("HBB", "-78", 5227099),       # -28 TATA, promoter (upstream of the transcript)
])
def test_known_coordinates(gene, cpos, expect):
    r = subprocess.run([sys.executable, os.path.join(BIN, "hgvs_map.py"),
                        "--gene", gene, "--cpos", cpos],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert f":{expect}\t" in r.stdout, r.stdout


@needs_annotsv
def test_mapper_refuses_to_guess():
    """A splice-crossing offset must raise, not return a plausible wrong base."""
    r = subprocess.run([sys.executable, os.path.join(BIN, "hgvs_map.py"),
                        "--gene", "HBB", "--cpos", "20+1"],
                       capture_output=True, text=True)
    assert r.returncode != 0
    assert "crosses a splice boundary" in (r.stdout + r.stderr)


def _panel_rows(name):
    path = os.path.join(REPO, "assets", name)
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if not l.startswith("#") and l.strip()]
    header = lines[0].split("\t")
    return header, [dict(zip(header, l.split("\t"))) for l in lines[1:]]


@pytest.mark.parametrize("panel", ["hba_pathogenic_sites.tsv", "hbb_pathogenic_sites.tsv"])
def test_panel_wellformed(panel):
    header, rows = _panel_rows(panel)
    assert rows, f"{panel} has no sites"
    for col in ("gene", "allele", "hgvs_c", "chrom", "pos", "strand",
                "coding_ref", "coding_alt", "genomic_ref", "genomic_alt"):
        assert col in header, f"{panel} missing column {col}"
    seen = set()
    for r in rows:
        key = (r["gene"], r["hgvs_c"])
        assert key not in seen, f"{panel}: duplicate site {key}"
        seen.add(key)
        assert r["pos"].isdigit(), f"{panel}: non-numeric pos {r['pos']}"
        assert r["chrom"].startswith("chr")


def test_minus_strand_bases_are_complemented():
    """HBB is minus-strand; a panel that forgets this calls nothing at all."""
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    _, rows = _panel_rows("hbb_pathogenic_sites.tsv")
    checked = 0
    for r in rows:
        if r["coding_ref"] == "-":       # indel rows carry only an anchor base
            continue
        assert r["strand"] == "-", "HBB should be annotated minus-strand"
        assert r["genomic_ref"] == comp[r["coding_ref"]], f"{r['allele']}: ref not complemented"
        assert r["genomic_alt"] == comp[r["coding_alt"]], f"{r['allele']}: alt not complemented"
        checked += 1
    assert checked >= 8, "expected most HBB sites to be SNVs"


def test_contract_columns_match_fixture():
    """The contract's declared columns and the example file must not drift apart."""
    with open(CONTRACT) as fh:
        declared = [l.split("`")[1] for l in fh
                    if l.startswith("| `") and "`" in l[3:]]
    with open(FIXTURE) as fh:
        lines = fh.read().splitlines()
    fixture_cols = lines[0].split("\t")
    assert declared == fixture_cols, (
        f"contract declares {declared}\nfixture has  {fixture_cols}")
    assert len(lines) == 2, "fixture must be exactly one header + one data row"
    assert len(lines[1].split("\t")) == len(fixture_cols)


def test_fixture_states_interpretation_incomplete():
    """SVcaller measures; it must never claim to have interpreted."""
    with open(FIXTURE) as fh:
        header, row = [l.split("\t") for l in fh.read().splitlines()]
    assert dict(zip(header, row))["interpretation_complete"] == "false"


def test_fixture_declares_beta_not_screened():
    """The commonest false-reassurance path: beta-thal silently implied covered."""
    with open(FIXTURE) as fh:
        header, row = [l.split("\t") for l in fh.read().splitlines()]
    assert "beta_globin" in dict(zip(header, row))["not_screened"]
