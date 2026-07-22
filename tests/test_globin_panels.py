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


CODON = {
 'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
 'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
 'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
 'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
 'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
 'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
 'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
 'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
REF_FASTA = "/data/alvin/ref/GRCh38/hg38.fa"

needs_ref = pytest.mark.skipif(
    not os.path.exists(REF_FASTA), reason="GRCh38 FASTA not present")


def _coding_base(model, hgvs):
    import subprocess as sp
    g = model.resolve(hgvs)
    out = sp.run(["samtools", "faidx", REF_FASTA, f"{model.chrom}:{g}-{g}"],
                 capture_output=True, text=True, check=True).stdout
    b = "".join(l.strip() for l in out.splitlines()[1:]).upper()
    return b if model.strand == "+" else COMP[b]


def _models():
    sys.path.insert(0, BIN)
    from hgvs_map import load_models
    return load_models(ANNOTSV_BED, {"HBB", "HBA1", "HBA2"})


@needs_annotsv
@needs_ref
@pytest.mark.parametrize("gene,intron_end,exon_start", [
    ("HBB", "92+1", "93-1"),     # IVS-1
    ("HBB", "315+1", "316-1"),   # IVS-2
])
def test_splice_sites_obey_gt_ag(gene, intron_end, exon_start):
    """A wrong exon boundary shows up immediately as a broken GT..AG rule."""
    m = _models()[gene]
    donor = _coding_base(m, intron_end) + _coding_base(m, intron_end[:-1] + "2")
    acc = _coding_base(m, exon_start[:-1] + "2") + _coding_base(m, exon_start)
    assert donor == "GT", f"{gene} {intron_end} donor is {donor}, not GT"
    assert acc == "AG", f"{gene} {exon_start} acceptor is {acc}, not AG"


@needs_annotsv
@needs_ref
@pytest.mark.parametrize("allele,gene,cpos,cref,calt,legacy,consequence", [
    ("HbS",                "HBB",  20,  "A", "T", 6,   "E6V"),
    ("HbC",                "HBB",  19,  "G", "A", 6,   "E6K"),
    ("HbE",                "HBB",  79,  "G", "A", 26,  "E26K"),
    ("CD17",               "HBB",  52,  "A", "T", 17,  "K17*"),
    ("CD39",               "HBB",  118, "C", "T", 39,  "Q39*"),
    ("Hb Quong Sze",       "HBA2", 377, "T", "C", 125, "L125P"),
    ("Hb Adana",           "HBA2", 179, "G", "A", 59,  "G59D"),
    ("Hb Constant Spring", "HBA2", 427, "T", "C", 142, "*142Q"),
])
def test_allele_name_matches_hgvs(allele, gene, cpos, cref, calt, legacy, consequence):
    """Legacy globin names number codons from the MATURE protein (no initiator
    Met), so legacy codon N is c.(3N+1)..c.(3N+3). If the curated name and the
    curated HGVS disagree, the translated consequence exposes it."""
    m = _models()[gene]
    start = 3 * legacy + 1
    idx = cpos - start
    assert 0 <= idx <= 2, f"{allele}: c.{cpos} is not inside legacy codon {legacy}"
    codon = "".join(_coding_base(m, str(start + i)) for i in range(3))
    assert codon[idx] == cref, f"{allele}: FASTA has {codon[idx]}, curated ref is {cref}"
    mut = codon[:idx] + calt + codon[idx + 1:]
    got = f"{CODON[codon]}{legacy}{CODON[mut]}"
    assert got == consequence, f"{allele}: consequence is {got}, expected {consequence}"


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


def _tsv(name):
    path = os.path.join(REPO, "assets", name)
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if not l.startswith("#") and l.strip()]
    header = lines[0].split("\t")
    return [dict(zip(header, l.split("\t"))) for l in lines[1:]]


def test_hba_segments_are_contiguous_and_ordered():
    """Segments tile the cluster; a gap or overlap would double-count or lose depth."""
    path = os.path.join(REPO, "assets", "hba_segments.bed")
    rows = []
    with open(path) as fh:
        for l in fh:
            if l.startswith("#") or not l.strip():
                continue
            f = l.split("\t")
            rows.append((f[0], int(f[1]), int(f[2]), f[3]))
    assert len(rows) == 5, "expected HBZ, INTER_Z_A, HBA2, INTER_A2_A1, HBA1"
    for i in range(1, len(rows)):
        prev_end, start = rows[i - 1][2], rows[i][1]
        assert start == prev_end, f"segments not contiguous at {rows[i][3]}: {prev_end} -> {start}"
    for _, s, e, name in rows:
        assert e > s, f"{name} is empty or inverted"


def test_segments_declare_reliability():
    """Measured on THAL1/THAL2 2026-07-22: the five segments are NOT
    interchangeable, and a caller that averages all of them against one
    genome-wide control emits wrong calls in two distinct ways.

      HBZ        reads 0.71 in THAL2, which has NO deletion -> a global 0.8
                 threshold false-POSITIVES. Needs its own baseline.
      INTER_Z_A  mapping inflation (1.4-1.9 over chr16:155000-162000) cancels
                 out the real deletion at 164000-172875, averaging to 0.99 -
                 reporting 'intact' while containing half a --SEA deletion.
    """
    path = os.path.join(REPO, "assets", "hba_segments.bed")
    rel = {}
    with open(path) as fh:
        for l in fh:
            if l.startswith("#") or not l.strip():
                continue
            f = l.split("\t")
            assert len(f) >= 5, f"segment row missing the reliability column: {l[:60]}"
            rel[f[3]] = f[4]
    assert rel["HBZ"] == "needs_own_baseline"
    assert rel["INTER_Z_A"] == "do_not_average"
    for seg in ("HBA2", "INTER_A2_A1", "HBA1"):
        assert rel[seg] == "good", f"{seg} should be a clean diagnostic segment"


def test_hbz_dependent_discrimination_is_flagged():
    """--SEA|--MED vs --FIL|--THAI is decided ENTIRELY by whether HBZ is lost.
    HBZ is the segment that needs its own baseline, so that discrimination is
    uncalibrated until a known-normal cohort supplies one. If someone marks HBZ
    'good' without doing that, this fails."""
    path = os.path.join(REPO, "assets", "hba_segments.bed")
    hbz = [l for l in open(path) if not l.startswith("#") and "\tHBZ\t" in l]
    assert hbz, "HBZ segment missing"
    assert "needs_own_baseline" in hbz[0], (
        "HBZ is flagged usable, but it is the sole discriminator between the two "
        "degenerate 2-gene-deletion groups — it must be calibrated first")


def test_deletion_alleles_declare_degeneracy():
    """Depth alone cannot separate --SEA from --MED, nor --FIL from --THAI. A
    caller that names one of them from depth is overclaiming, so the file must
    say so and the agent must gate on it."""
    rows = _tsv("hba_deletion_alleles.tsv")
    by = {r["allele"]: r for r in rows}
    for pair in (("--SEA", "--MED"), ("--FIL", "--THAI")):
        for a in pair:
            assert a in by, f"{a} missing from the allele table"
            val = by[a]["depth_distinguishable"]
            assert val.startswith("no:"), f"{a} should be flagged depth-degenerate, got {val!r}"
            assert all(p in val for p in pair), f"{a} degeneracy group should name {pair}"
    # a signature unique to one allele must NOT be flagged degenerate
    assert by["-a4.2"]["depth_distinguishable"] == "yes"


def test_deletion_alleles_have_provenance():
    rows = _tsv("hba_deletion_alleles.tsv")
    assert rows, "no alleles defined"
    for r in rows:
        assert r["basis"] in ("observed", "literature"), \
            f"{r['allele']}: basis must state where the definition came from"
        assert r["note"].strip(), f"{r['allele']}: needs a note"
    # the only 'observed' allele is the one we actually measured
    obs = [r["allele"] for r in rows if r["basis"] == "observed"]
    assert obs == ["--SEA"], f"unexpected observed alleles: {obs}"


def test_triplication_is_present():
    """A caller that only looks for losses silently misses anti-3.7."""
    rows = _tsv("hba_deletion_alleles.tsv")
    trip = [r for r in rows if r["class"] == "triplication"]
    assert trip, "anti-3.7 triplication missing — a loss-only caller would never report it"
    assert int(trip[0]["alpha_genes_lost"]) < 0, "triplication should GAIN a gene"


def test_sizes_are_documentary_not_coordinates():
    """approx_size must never be mistakable for a usable coordinate."""
    for r in _tsv("hba_deletion_alleles.tsv"):
        s = r["approx_size"]
        assert s.startswith(("~", "+")), f"{r['allele']}: size {s!r} should be marked approximate"
        assert "kb" in s, f"{r['allele']}: size {s!r} should carry units"


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
