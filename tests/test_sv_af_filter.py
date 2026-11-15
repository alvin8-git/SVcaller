"""bin/sv_af_filter.py: the population-AF column used to drop a common SV is the one of its own type.

AnnotSV fills B_gain/B_loss/B_ins/B_inv for every SV (a benign region that completely contains
it), so reading B_gain_AFmax for a deletion drops deletions that merely sit inside a common
gain region and never drops a common deletion.  Regression for audit item 5 / item 4.
Also pins main.nf's provenance call (GString elements broke ProcessBuilder).
"""
import csv
import os
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO, "bin", "sv_af_filter.py")
COLS = ["SV_type", "B_gain_AFmax", "B_loss_AFmax", "B_ins_AFmax", "B_inv_AFmax", "tag"]


def run(tmp_path, rows):
    src, dst = tmp_path / "in.tsv", tmp_path / "out.tsv"
    with open(src, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(COLS)
        w.writerows(rows)
    subprocess.run([sys.executable, SCRIPT, str(src), str(dst), "0.01"], check=True)
    return [r["tag"] for r in csv.DictReader(open(dst), delimiter="\t")]


def test_del_uses_loss_not_gain(tmp_path):
    kept = run(tmp_path, [
        ["DEL", "0.5", "", "", "", "del_in_common_gain"],      # gain irrelevant -> kept
        ["DEL", "", "0.3", "", "", "common_del"],              # common loss -> dropped
        ["DEL", "", "0.0099", "", "", "rare_del"],
        ["DEL", "", "", "", "", "no_af"],
    ])
    assert kept == ["del_in_common_gain", "rare_del", "no_af"]


def test_dup_still_uses_gain(tmp_path):
    kept = run(tmp_path, [["DUP", "0.2", "", "", "", "common"], ["DUP", "", "0.9", "", "", "loss_ignored"]])
    assert kept == ["loss_ignored"]


def test_ins_inv_use_own_column_and_untyped_kept(tmp_path):
    kept = run(tmp_path, [
        ["INS", "0.9", "0.9", "0.5", "", "common_ins"], ["INS", "0.9", "0.9", "", "", "ins_ok"],
        ["INV", "", "", "", "0.5", "common_inv"], ["INV", "0.9", "", "", "", "inv_ok"],
        ["TRA", "0.9", "0.9", "0.9", "0.9", "tra_kept"], ["BND", "0.9", "", "", "", "bnd_kept"],
    ])
    assert kept == ["ins_ok", "inv_ok", "tra_kept", "bnd_kept"]


def test_unparsable_af_is_not_a_drop(tmp_path):
    assert run(tmp_path, [["DEL", "", "n/a", "", "", "x"]]) == ["x"]


def test_main_nf_provenance_uses_strings_and_runs_at_start():
    text = open(os.path.join(REPO, "main.nf")).read()
    assert "List<String> cmd" in text and "as List<String>" not in text
    wf = text[text.index("\nworkflow {"):]
    assert wf.index("write_provenance_stub()") < wf.index("SVCALLER(")   # stub written before any process can fail
    assert "write_provenance()" in text[text.index("workflow.onComplete"):]
