"""provenance.json writer: schema, 'not reported' fallback, and never-fail guarantees."""
import json
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import write_provenance as wp  # noqa: E402

KEYS = {"pipeline", "svcaller_commit", "generated_utc", "sample", "reference", "tools", "images"}
STATE = ("v1-3-gabc", "abc123-dirty")


def _q(image, cmd):
    if "delly" in image:
        raise RuntimeError("boom")          # query raises
    if "annotsv" in image:
        return None                          # tool cannot report
    return "9.9.9"


def test_schema_and_not_reported(tmp_path):
    imgs = {"MANTA_CALL": "q/manta:1", "DELLY_CALL_SVTYPE": "q/delly:1", "ANNOTSV": "q/annotsv:1"}
    p = wp.write(str(tmp_path), "S1", "/ref/x.fa", imgs, querier=_q, state=STATE)
    doc = json.load(open(p))
    assert p == str(tmp_path / "S1" / "provenance.json")
    assert KEYS <= set(doc)
    assert doc["svcaller_commit"] == "abc123-dirty" and doc["sample"] == "S1"
    assert doc["reference"] == {"path": "/ref/x.fa"}
    assert doc["tools"]["manta"] == "9.9.9" and doc["images"]["manta"] == "q/manta:1"
    assert doc["tools"]["delly"] == "not reported"
    assert doc["tools"]["annotsv"] == "not reported"
    assert doc["tools"]["gatk"] == "not reported" and doc["images"]["gatk"] == "not reported"
    assert "abc123" not in "".join(doc["tools"].values())  # repo commit never a tool version


def test_main_exit_zero_when_docker_missing(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", "")
    assert wp.main(["--outdir", str(tmp_path), "--sample", "S2", "--ref-fasta", "/r.fa"]) == 0
    doc = json.load(open(tmp_path / "S2" / "provenance.json"))
    assert KEYS <= set(doc) and set(doc["tools"].values()) == {"not reported"}


def test_unwritable_outdir_does_not_raise(tmp_path):
    f = tmp_path / "afile"
    f.write_text("x")
    assert wp.write(str(f / "sub"), "S3", "/r.fa", {}, querier=_q, state=STATE) is None
    assert wp.main(["--outdir", str(f / "sub"), "--sample", "S3"]) == 0
    assert wp.main([]) == 0  # bad args
