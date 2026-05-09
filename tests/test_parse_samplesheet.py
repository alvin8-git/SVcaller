import json, subprocess, textwrap
from pathlib import Path
import pytest

SCRIPT = Path(__file__).parent.parent / "bin" / "parse_samplesheet.py"

def run(csv_text: str, tmp_path) -> subprocess.CompletedProcess:
    p = tmp_path / "sheet.csv"
    p.write_text(textwrap.dedent(csv_text))
    return subprocess.run(["python3", str(SCRIPT), str(p)],
                          capture_output=True, text=True)

def test_rejects_missing_sample(tmp_path):
    r = run("sample,bam\n,/nonexistent.bam\n", tmp_path)
    assert r.returncode == 1
    assert "required" in r.stderr

def test_rejects_bam_and_fastq(tmp_path):
    r = run("sample,fastq_1,fastq_2,bam\nS1,/a.fq.gz,/b.fq.gz,/c.bam\n", tmp_path)
    assert r.returncode == 1
    assert "not both" in r.stderr

def test_rejects_missing_fastq2(tmp_path):
    r = run("sample,fastq_1,fastq_2\nS1,/a.fq.gz,\n", tmp_path)
    assert r.returncode == 1
    assert "both" in r.stderr.lower()

def test_accepts_bam_row(tmp_path):
    bam = tmp_path / "test.bam"
    bam.touch()
    r = run(f"sample,bam\nS1,{bam}\n", tmp_path)
    assert r.returncode == 0
    d = json.loads(r.stdout.strip())
    assert d["id"] == "S1"
    assert d["input_type"] == "bam"
