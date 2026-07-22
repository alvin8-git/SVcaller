"""Guard the synthetic junction fixture used by channel 3 (HBA_JUNCTION).

These tests do NOT test a junction caller — none exists yet. They assert the
FIXTURE still contains what channel 3 will be written against, so that when the
caller arrives it is developed against a known-good input rather than one that
quietly rotted.

The fixture is committed (~2 KB); the reference is only needed to rebuild it.
"""
import os
import subprocess

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BAM = os.path.join(REPO, "tests", "fixtures", "alpha_junction.bam")

CHROM = "chr16"
DEL_START, DEL_END = 165001, 185000
LEFT_BP, RIGHT_BP = DEL_START - 1, DEL_END + 1     # 165000 -> 185001

pytestmark = pytest.mark.skipif(
    subprocess.run(["which", "samtools"], capture_output=True).returncode != 0,
    reason="samtools not on PATH")


def _view(*args):
    r = subprocess.run(["samtools", "view", *args, BAM],
                       capture_output=True, text=True, check=True)
    return [l.split("\t") for l in r.stdout.splitlines() if l.strip()]


@pytest.fixture(scope="module")
def reads():
    assert os.path.exists(BAM), (
        "tests/fixtures/alpha_junction.bam is missing. If git dropped it, check "
        "the *.bam rule in .gitignore still has the tests/fixtures/ exception.")
    return _view()


def test_fixture_is_valid_bam():
    assert subprocess.run(["samtools", "quickcheck", BAM]).returncode == 0
    assert os.path.exists(BAM + ".bai"), "fixture must ship its index"
    assert os.path.getsize(BAM) < 200_000, "fixture should stay tiny enough to commit"


def test_header_declares_it_is_synthetic():
    """Nobody should ever mistake this for real sample data."""
    hdr = subprocess.run(["samtools", "view", "-H", BAM],
                         capture_output=True, text=True, check=True).stdout
    assert "SYNTHETIC" in hdr.upper()
    assert str(DEL_START) in hdr and str(DEL_END) in hdr, \
        "header must record the ground-truth breakpoint"


def test_split_reads_clip_exactly_at_the_breakpoint(reads):
    """The whole point: a caller must be able to recover 165000 and 185001 to
    the base. Aligned blocks must END at LEFT_BP / START at RIGHT_BP."""
    import re
    left = [r for r in reads if r[0].startswith("split_left_")]
    right = [r for r in reads if r[0].startswith("split_right_")]
    assert len(left) >= 5 and len(right) >= 5, "need split support on both sides"

    for r in left:
        pos, cigar = int(r[3]), r[5]
        m = re.fullmatch(r"(\d+)M(\d+)S", cigar)
        assert m, f"{r[0]}: expected <n>M<n>S, got {cigar}"
        assert pos + int(m.group(1)) - 1 == LEFT_BP, \
            f"{r[0]}: aligned block ends at {pos + int(m.group(1)) - 1}, not {LEFT_BP}"

    for r in right:
        pos, cigar = int(r[3]), r[5]
        m = re.fullmatch(r"(\d+)S(\d+)M", cigar)
        assert m, f"{r[0]}: expected <n>S<n>M, got {cigar}"
        assert pos == RIGHT_BP, f"{r[0]}: aligned block starts at {pos}, not {RIGHT_BP}"


def test_clipped_sequence_matches_the_far_side(reads):
    """Soft-clipped bases are real reference from ACROSS the deletion, so a
    caller that realigns the clip resolves the partner breakpoint. If this ever
    became random sequence the fixture would still look right but would silently
    stop testing realignment."""
    import re
    ref = "/data/alvin/ref/GRCh38/hg38.fa"
    if not os.path.exists(ref):
        pytest.skip("reference not present; clip-content check needs it")

    def fetch(s, e):
        out = subprocess.run(["samtools", "faidx", ref, f"{CHROM}:{s}-{e}"],
                             capture_output=True, text=True, check=True).stdout
        return "".join(l.strip() for l in out.splitlines()[1:]).upper()

    r = next(x for x in reads if x[0] == "split_left_0")
    m = re.fullmatch(r"(\d+)M(\d+)S", r[5])
    clip = r[9][int(m.group(1)):]
    assert clip == fetch(RIGHT_BP, RIGHT_BP + len(clip) - 1), \
        "left split's clip does not match the sequence at the right breakpoint"

    r = next(x for x in reads if x[0] == "split_right_0")
    m = re.fullmatch(r"(\d+)S(\d+)M", r[5])
    n = int(m.group(1))
    assert r[9][:n] == fetch(LEFT_BP - n + 1, LEFT_BP), \
        "right split's clip does not match the sequence at the left breakpoint"


def test_has_reference_spanning_reads(reads):
    """THE FIXTURE'S MOST IMPORTANT PROPERTY. The deletion is heterozygous, so
    the intact allele is present: reads cross LEFT_BP with no clip at all. A
    fixture of pure junction reads would let a caller that ignores allele
    balance report a HOMOZYGOUS deletion and still pass every other test here."""
    span = [r for r in reads
            if r[0].startswith("refspan_")
            and "S" not in r[5]
            and int(r[3]) < LEFT_BP < int(r[3]) + 150]
    assert len(span) >= 8, f"only {len(span)} unclipped reads cross the breakpoint"

    split = [r for r in reads if r[0].startswith("split_")]
    ratio = len(split) / (len(split) + len(span))
    assert 0.3 < ratio < 0.7, (
        f"split fraction {ratio:.2f} does not look heterozygous — the fixture "
        "should not be callable as hom-del or as reference")


def test_discordant_pairs_span_the_deletion(reads):
    disc = [r for r in reads if r[0].startswith("disc_")]
    assert len(disc) >= 10, "need discordant pair support"
    for r in disc:
        assert not (int(r[1]) & 0x2), f"{r[0]}: discordant read flagged proper-pair"
    tlens = [abs(int(r[8])) for r in disc if int(r[8]) != 0]
    assert tlens and min(tlens) > 15_000, \
        "discordant inserts should reflect the ~20 kb deletion"


def test_distal_reads_are_clean(reads):
    """A caller must not fire away from the locus."""
    distal = [r for r in reads if r[0].startswith("normal_")]
    assert len(distal) >= 6
    for r in distal:
        assert "S" not in r[5], f"{r[0]}: distal read should not be clipped"
        assert int(r[1]) & 0x2, f"{r[0]}: distal read should be a proper pair"
