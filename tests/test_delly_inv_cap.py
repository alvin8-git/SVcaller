"""delly INV time cap (bin/delly_call_capped.sh): the sample continues, and says so.

A fake `delly` on PATH stands in for the real one; the real process was also run in its
container (see docs/CHANGES.md). What is asserted: INV is the only type that is capped; a cap
leaves a valid header-only BGZF VCF plus a one-row not-assessed record; a completed or
failed delly never produces that row.
"""
import gzip
import struct
import os
import shutil
import stat
import subprocess

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO, "bin", "delly_call_capped.sh")
HDR = "sample\tcaller\tsvtype\treason\tlimit\n"

pytestmark = pytest.mark.skipif(shutil.which("bgzip") is None, reason="bgzip not on PATH")


def _run(tmp_path, svtype, limit, fake_body, bam="none.bam"):
    bindir = tmp_path / "fakebin"
    bindir.mkdir(exist_ok=True)
    d = bindir / "delly"
    d.write_text("#!/bin/sh\n" + fake_body)
    d.chmod(d.stat().st_mode | stat.S_IEXEC)
    (tmp_path / "ref.fa.fai").write_text("chr1\t1000\t6\t50\t51\nchr2\t500\t2000\t50\t51\n")
    env = dict(os.environ, PATH=f"{bindir}:{os.environ['PATH']}")
    return subprocess.run(
        ["sh", SCRIPT, "S1", svtype, str(limit), "ref.fa", "ref.fa.fai", bam],
        cwd=tmp_path, env=env, capture_output=True, text=True, timeout=60)


OK = 'while [ $# -gt 0 ]; do [ "$1" = -o ] && o=$2; shift; done; echo "##ok" > "$o"\n'
HANG = "exec sleep 30\n"


def test_cap_fires_for_inv(tmp_path):
    r = _run(tmp_path, "INV", 1, HANG)
    assert r.returncode == 0, r.stderr
    assert "NOT ASSESSED" in r.stderr and "WARN" in r.stderr
    rows = (tmp_path / "S1.sv_not_assessed.tsv").read_text().splitlines(True)
    assert rows[0] == HDR
    assert rows[1:] == ["S1\tdelly\tINV\ttime_limit_exceeded\t1s\n"]
    lines = gzip.open(tmp_path / "S1.delly.INV.vcf", "rt").read().splitlines()
    assert lines[0] == "##fileformat=VCFv4.2"
    assert lines[-1].startswith("#CHROM") and not any(l for l in lines if not l.startswith("#"))
    assert "##contig=<ID=chr2,length=500>" in lines


def test_inv_that_finishes_has_header_only_record(tmp_path):
    r = _run(tmp_path, "INV", 30, OK)
    assert r.returncode == 0, r.stderr
    assert (tmp_path / "S1.sv_not_assessed.tsv").read_text() == HDR
    assert (tmp_path / "S1.delly.INV.vcf").read_text() == "##ok\n"


def test_inv_without_cap_is_not_capped(tmp_path):
    r = _run(tmp_path, "INV", 0, OK)
    assert r.returncode == 0 and (tmp_path / "S1.sv_not_assessed.tsv").read_text() == HDR


def test_real_failure_fails_and_is_not_recorded_as_a_cap(tmp_path):
    r = _run(tmp_path, "INV", 30, "exit 3\n")
    assert r.returncode == 3
    assert "NOT ASSESSED" not in r.stderr
    assert (tmp_path / "S1.sv_not_assessed.tsv").read_text() == HDR   # header only, no row


def test_other_types_ignore_the_limit_and_write_no_record(tmp_path):
    r = _run(tmp_path, "DEL", 1, "sleep 3\n" + OK)
    assert r.returncode == 0, r.stderr
    assert (tmp_path / "S1.delly.DEL.vcf").read_text() == "##ok\n"
    assert not (tmp_path / "S1.sv_not_assessed.tsv").exists()


def test_cap_still_fires_when_delly_ignores_term(tmp_path):
    """The watchdog escalates to KILL after its grace period; the capped path is unchanged."""
    import time
    t0 = time.time()
    r = _run(tmp_path, "INV", 1, "trap '' TERM\nwhile :; do sleep 1; done\n")
    assert r.returncode == 0, r.stderr
    assert time.time() - t0 < 30
    assert (tmp_path / "S1.sv_not_assessed.tsv").read_text() == \
        HDR + "S1\tdelly\tINV\ttime_limit_exceeded\t1s\n"


def _fake_bam(path, header_text):
    """A gzip stream holding a BAM header (enough for zcat; delly itself is faked)."""
    t = header_text.encode()
    with gzip.open(path, "wb") as f:
        f.write(b"BAM\x01" + struct.pack("<i", len(t)) + t + struct.pack("<i", 0))


def _sample_col(tmp_path):
    lines = gzip.open(tmp_path / "S1.delly.INV.vcf", "rt").read().splitlines()
    return lines[-1].split("\t")[-1]


def test_placeholder_sample_is_the_first_rg_sm_not_a_pg_decoy(tmp_path):
    _fake_bam(tmp_path / "x.bam",
              "@HD\tVN:1.6\n@PG\tID:a\tPN:t\tCL:run SM:WRONG -R @RG\\tSM:WRONG2\n"
              "@RG\tID:g1\tLB:l\tSM:REAL\n@RG\tID:g2\tSM:OTHER\n")
    assert _run(tmp_path, "INV", 1, HANG, bam="x.bam").returncode == 0
    assert _sample_col(tmp_path) == "REAL"


def test_placeholder_sample_without_rg_is_the_bam_file_stem(tmp_path):
    _fake_bam(tmp_path / "my_sample.bam", "@HD\tVN:1.6\n@PG\tID:a\tCL:x SM:WRONG\n")
    assert _run(tmp_path, "INV", 1, HANG, bam="my_sample.bam").returncode == 0
    assert _sample_col(tmp_path) == "my_sample"


def test_sentinel_with_an_unrelated_exit_status_is_a_failure(tmp_path):
    """Only 143/137 count as killed by the watchdog; anything else fails the task."""
    r = _run(tmp_path, "INV", 30, ": > inv_cap_fired\nexit 3\n")
    assert r.returncode == 3
    assert "NOT ASSESSED" not in r.stderr
    assert (tmp_path / "S1.sv_not_assessed.tsv").read_text() == HDR
