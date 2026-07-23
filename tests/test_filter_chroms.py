"""Guard bin/filter_chroms.sh — the canonical-chromosome BAM filter.

WHY. FILTER_CHROMS was a single-threaded awk that pinned one core for ~70 min on
a 79 GB BAM. bin/filter_chroms.sh parallelizes it by chromosome. The parallel
version MUST produce content-identical output, because every SV caller runs on
its result: a read dropped wrong, a header contig kept wrong, and Manta either
crashes or silently assembles nothing.

This test drives the script on a tiny synthetic BAM that exercises all four
transforms and asserts the exact output. It is the contract, not a smoke test.

The four transforms (see the script header):
  1. header keeps only canonical @SQ, in FAI order;
  2. non-canonical entries stripped from SA:Z tags;
  3. reads whose mate (RNEXT) is on a non-canonical contig are dropped;
  4. RNEXT=="*" with PNEXT!=0 has PNEXT zeroed.
Reads on a non-canonical contig are dropped by virtue of the per-chromosome
region queries (unmapped reads too), exactly as the old `view -h <regions>` did.
"""
import os
import shutil
import subprocess

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO, "bin", "filter_chroms.sh")

pytestmark = pytest.mark.skipif(
    not (shutil.which("samtools") and shutil.which("bash")),
    reason="samtools and bash required")

ALT = "chr1_KI270706v1_random"

# QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL [TAGS]
# r1  chr1, mate on chr1 (=)                         -> kept
# r2  chr1, mate on the alt contig                   -> DROPPED (mate non-canonical)
# r3  chr2, SA:Z listing chr1 AND the alt contig     -> kept, SA stripped to chr1
# r4  chr1, RNEXT="*", PNEXT=500                      -> kept, PNEXT zeroed
# r5  on the alt contig                              -> DROPPED (non-canonical RNAME)
_SAM = f"""\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:10000
@SQ\tSN:chr2\tLN:10000
@SQ\tSN:{ALT}\tLN:5000
r1\t99\tchr1\t100\t60\t10M\t=\t300\t210\tACGTACGTAC\tIIIIIIIIII
r4\t89\tchr1\t400\t60\t10M\t*\t500\t0\tACGTACGTAC\tIIIIIIIIII
r2\t97\tchr1\t200\t60\t10M\t{ALT}\t50\t0\tACGTACGTAC\tIIIIIIIIII
r3\t0\tchr2\t300\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tSA:Z:chr1,500,+,10M,60,0;{ALT},10,+,10M,60,0;
r5\t0\t{ALT}\t50\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII
"""


def _sh(*args):
    r = subprocess.run(args, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r.stdout


@pytest.fixture
def filtered(tmp_path):
    sam = tmp_path / "in.sam"
    sam.write_text(_SAM)
    unsorted = tmp_path / "u.bam"
    bam = tmp_path / "in.bam"
    with open(unsorted, "wb") as fh:
        fh.write(subprocess.run(["samtools", "view", "-bS", str(sam)],
                                capture_output=True).stdout)
    _sh("samtools", "sort", "-o", str(bam), str(unsorted))
    _sh("samtools", "index", str(bam))
    # FAI: only column 1 (name) is read; canonical order here is chr1, chr2.
    fai = tmp_path / "ref.fa.fai"
    fai.write_text(f"chr1\t10000\t6\t60\t61\nchr2\t10000\t10200\t60\t61\n"
                   f"{ALT}\t5000\t20400\t60\t61\n")
    out = tmp_path / "out.bam"
    _sh("bash", SCRIPT, str(bam), str(fai), "2", str(out))
    return out


def _reads(bam):
    return {l.split("\t")[0]: l.split("\t")
            for l in _sh("samtools", "view", str(bam)).splitlines()}


def _sq_names(bam):
    return [f.split("SN:")[1] for l in _sh("samtools", "view", "-H", str(bam)).splitlines()
            if l.startswith("@SQ") for f in l.split("\t") if f.startswith("SN:")]


def test_header_keeps_only_canonical_sq_in_fai_order(filtered):
    assert _sq_names(filtered) == ["chr1", "chr2"], "alt contig leaked into header"


def test_reads_with_noncanonical_mate_or_rname_are_dropped(filtered):
    reads = _reads(filtered)
    assert set(reads) == {"r1", "r3", "r4"}, \
        f"expected r1,r3,r4; got {sorted(reads)} (r2 mate-on-alt / r5 on-alt must go)"


def test_sa_tag_stripped_of_noncanonical(filtered):
    r3 = _reads(filtered)["r3"]
    sa = [f for f in r3 if f.startswith("SA:Z:")]
    assert sa == ["SA:Z:chr1,500,+,10M,60,0;"], f"SA:Z not stripped to canonical: {sa}"


def test_pnext_zeroed_when_rnext_star(filtered):
    r4 = _reads(filtered)["r4"]
    assert r4[6] == "*" and r4[7] == "0", f"PNEXT not zeroed for RNEXT=*: {r4[6:8]}"


def test_output_is_coordinate_sorted_and_indexable(filtered, tmp_path):
    # samtools index only succeeds on a coordinate-sorted BAM; the script indexes
    # it, so a stale/missing index would already have failed. Re-check explicitly.
    assert (str(filtered) + ".bai") and os.path.exists(str(filtered) + ".bai")
    hdr = _sh("samtools", "view", "-H", str(filtered))
    assert "SO:coordinate" in hdr
