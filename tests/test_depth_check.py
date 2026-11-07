"""Read-depth consistency check (bin/depth_check.py): ratio -> status, and a synthetic BAM.

A DEL that reads at ~1.0x of its flanks is refuted; a real (het) deletion at ~0.5x is
supported.  This is the contract that stops depth-contradicted calls (e.g. the HG001
VHL/PPARG/NSF deletions) from looking like findings downstream.
"""
import csv
import importlib.util
import os
import shutil
import subprocess

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("depth_check", os.path.join(REPO, "bin", "depth_check.py"))
dc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dc)


@pytest.mark.parametrize("t,ratio,want", [
    ("DEL", 0.50, "supported"), ("DEL", 0.65, "supported"), ("DEL", 0.75, "inconclusive"),
    ("DEL", 0.85, "refuted"), ("DEL", 1.00, "refuted"),
    ("DUP", 1.50, "supported"), ("DUP", 1.35, "supported"), ("DUP", 1.25, "inconclusive"),
    ("DUP", 1.15, "refuted"), ("DUP", 1.00, "refuted"),
    ("INV", 1.00, "na"), ("DEL", None, "na"),
])
def test_depth_support(t, ratio, want):
    assert dc.depth_support(t, ratio) == want


def test_inflated_flank_is_refuted_by_global_baseline():
    """NSF case: in 29, flank 46, global 30. Flank ratio 0.63 alone says supported; global 0.97 refutes."""
    rf, rg = 29 / 46, 29 / 30
    assert dc.depth_support("DEL", rf, rg) == "refuted"
    assert dc.depth_note(46, 30) == "flank inflated"
    assert dc.depth_note(15, 30) == "flank depressed"
    assert dc.depth_note(30, 30) == "."
    assert dc.depth_support("DEL", 0.5, 0.5) == "supported"      # both agree
    assert dc.depth_support("DEL", 0.5, 0.75) == "inconclusive"  # one not in supported band
    assert dc.depth_support("DUP", 1.5, 1.5) == "supported"


def test_ratio_na_on_shallow_flank():
    assert dc.depth_ratio(10.0, 4.9) is None
    assert dc.depth_ratio(15.0, 30.0) == 0.5


@pytest.mark.skipif(not (shutil.which("samtools")), reason="samtools not on PATH")
def test_two_dels_on_synthetic_bam(tmp_path):
    """chr1: 30x everywhere except 15x over 20,000-30,999 (a het DEL). Second DEL sits on flat 30x."""
    L = 60000
    sam = tmp_path / "t.sam"
    with open(sam, "w") as f:
        f.write(f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:{L}\n")
        n = 0
        for pos in range(1, L - 99, 10):          # 10 bp step x 100 bp reads = 10x per layer
            layers = 3
            if 20000 <= pos <= 30900:             # drop to 1.5 layers = 15x
                layers = 1 if (pos // 10) % 2 else 2
            for _ in range(layers):
                n += 1
                f.write(f"r{n}\t0\tchr1\t{pos}\t60\t100M\t*\t0\t0\t{'A'*100}\t{'I'*100}\n")
    bam = tmp_path / "t.bam"
    subprocess.run(["samtools", "sort", "-o", str(bam), str(sam)], check=True)
    subprocess.run(["samtools", "index", str(bam)], check=True)

    cols = ["SV_chrom", "SV_start", "SV_end", "SV_length", "SV_type"]
    rows = [dict(zip(cols, ["1", "20000", "31000", "-11001", "DEL"])),   # real
            dict(zip(cols, ["1", "40000", "51000", "-11001", "DEL"])),   # flat -> refuted
            dict(zip(cols, ["1", "100", "600", "-501", "DEL"])),         # <1 kb -> na
            dict(zip(cols, ["1", "40000", "51000", "11001", "INV"]))]    # not DEL/DUP -> na
    out = dc.annotate(rows, str(bam), threads=2, d_global=30.0)
    assert out[0]["DEPTH_SUPPORT"] == "supported" and float(out[0]["DEPTH_RATIO"]) < 0.65
    assert out[1]["DEPTH_SUPPORT"] == "refuted" and 0.9 < float(out[1]["DEPTH_RATIO"]) < 1.1
    assert out[2]["DEPTH_SUPPORT"] == "na" and out[2]["DEPTH_RATIO"] == "."
    assert out[3]["DEPTH_SUPPORT"] == "na"
    assert out[0]["DEPTH_GLOBAL"] == "30.00" and out[0]["DEPTH_NOTE"] == "."
