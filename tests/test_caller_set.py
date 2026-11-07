"""Default caller set (Manta + Delly; GRIDSS and MELT opt-in) and its provenance.

A two-caller run must produce a SUPP_VEC of width 2, say so in provenance, and a Finding
needs BOTH callers.  Jasmine itself is not run here: it emits one SUPP_VEC bit per line of
vcf_list.txt, so the width is the length of the list caller_set.sh writes.
"""
import gzip
import json
import os
import re
import subprocess
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(REPO, "bin")
sys.path.insert(0, BIN)
import caller_consensus as cc  # noqa: E402


def _vcf(path, n):
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for i in range(n):
            f.write(f"chr1\t{1000 + i}\tx{i}\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL\n")
    return os.path.basename(str(path))


def _run_set(tmp_path, specs):
    r = subprocess.run([os.path.join(BIN, "caller_set.sh"), "S1"] + specs,
                       cwd=tmp_path, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    lst = (tmp_path / "vcf_list.txt").read_text().split()
    prov = json.loads((tmp_path / "S1.sv_callers.json").read_text())
    return lst, prov


def test_two_caller_run_has_width_two_and_provenance(tmp_path):
    m, d, g, s, me = (_vcf(tmp_path / n, k) for n, k in
                      [("m.vcf", 3), ("d.vcf", 3), ("g.vcf", 0), ("s.vcf", 0), ("me.vcf", 0)])
    lst, prov = _run_set(tmp_path, [f"Manta={m}:always", f"Delly={d}:always",
                                    f"GRIDSS={g}:off", f"Scramble={s}:off", f"MELT={me}:off"])
    assert lst == [m, d]                                # Jasmine => SUPP_VEC width 2
    assert prov["supp_vec_width"] == 2
    assert prov["callers_in_supp_vec_order"] == ["Manta", "Delly"]
    assert prov["callers_ran"] == ["Manta", "Delly"]
    assert prov["callers_not_run"] == ["GRIDSS", "Scramble", "MELT"]


def test_optin_callers_take_their_bit_in_order(tmp_path):
    f = {n: _vcf(tmp_path / f"{n}.vcf", 2) for n in ("m", "d", "g", "s", "me")}
    lst, prov = _run_set(tmp_path, [f"Manta={f['m']}:always", f"Delly={f['d']}:always",
                                    f"GRIDSS={f['g']}:always", f"Scramble={f['s']}:ifcalls",
                                    f"MELT={f['me']}:ifcalls"])
    assert prov["callers_in_supp_vec_order"] == ["Manta", "Delly", "GRIDSS", "Scramble", "MELT"]
    assert prov["supp_vec_width"] == 5 and len(lst) == 5


def test_ran_but_no_calls_gets_no_bit_but_is_reported(tmp_path):
    m, d, s = (_vcf(tmp_path / n, k) for n, k in [("m.vcf", 2), ("d.vcf", 2), ("s.vcf", 0)])
    _, prov = _run_set(tmp_path, [f"Manta={m}:always", f"Delly={d}:always", f"Scramble={s}:ifcalls"])
    assert prov["supp_vec_width"] == 2
    assert prov["callers_ran"] == ["Manta", "Delly", "Scramble"]
    assert prov["callers_ran_no_calls"] == ["Scramble"]


def test_consensus_needs_both_manta_and_delly():
    c = ["Manta", "Delly"]
    assert cc.classify(c, "SUPP=2;SUPP_VEC=11") == ("Manta,Delly", "PASS")
    assert cc.classify(c, "SUPP=1;SUPP_VEC=10") == ("Manta", "SINGLE_CALLER")
    assert cc.classify(c, "SUPP=1;SUPP_VEC=01") == ("Delly", "SINGLE_CALLER")
    # width mismatch or missing vector is NA, never a silent PASS
    assert cc.classify(c, "SUPP_VEC=110")[1] == "NA"
    assert cc.classify(c, "SVTYPE=DEL")[1] == "NA"
    assert cc.classify(None, "SUPP_VEC=11")[1] == "NA"
    # wider sets keep the historical >=2 rule
    assert cc.classify(["Manta", "Delly", "GRIDSS"], "SUPP_VEC=101")[1] == "PASS"
    assert cc.classify(["Manta", "Delly", "GRIDSS"], "SUPP_VEC=001")[1] == "SINGLE_CALLER"


def test_filtered_table_cli_adds_columns_and_keeps_rows(tmp_path):
    vcf = tmp_path / "s.sv_merged.vcf.gz"
    with gzip.open(vcf, "wt") as f:
        f.write("##fileformat=VCFv4.2\n##SVcaller_callers=Manta,Delly\n#CHROM\tPOS\n")
    tsv = tmp_path / "in.tsv"
    tsv.write_text("SV_chrom\tINFO\nchr1\tSUPP=2;SUPP_VEC=11\nchr2\tSUPP=1;SUPP_VEC=10\n")
    out = tmp_path / "out.tsv"
    cc.main(str(tsv), str(vcf), str(out))
    rows = [l.split("\t") for l in out.read_text().splitlines()]
    assert rows[0] == ["SV_chrom", "INFO"] + cc.NEW_COLS
    assert rows[1][2:] == ["Manta,Delly", "2", "Manta,Delly", "PASS"]
    assert rows[2][2:] == ["Manta,Delly", "2", "Manta", "SINGLE_CALLER"]


def test_legacy_vcf_without_header_is_marked_legacy(tmp_path):
    vcf = tmp_path / "old.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\n")
    tsv = tmp_path / "in.tsv"
    tsv.write_text("SV_chrom\tINFO\nchr1\tSUPP=2;SUPP_VEC=1100\n")
    cc.main(str(tsv), str(vcf), str(tmp_path / "o.tsv"))
    row = (tmp_path / "o.tsv").read_text().splitlines()[1].split("\t")
    assert row[2:] == ["legacy", ".", ".", "NA"]


def test_html_report_uses_recorded_callers():
    import html_report as hr
    hr._RAN_CALLERS = ["Manta", "Delly"]
    try:
        assert hr._parse_supp_vec("SUPP=2;SUPP_VEC=11")["callers"] == "Manta/Delly"
        assert hr._sv_callers_text() == "Manta + Delly"
    finally:
        hr._RAN_CALLERS = None


def test_config_defaults_are_two_callers_with_optin_flags():
    conf = open(os.path.join(REPO, "nextflow.config")).read()
    for k, v in (("skip_gridss", "true"), ("skip_melt", "true"), ("skip_scramble", "true"),
                 ("run_gridss", "false"), ("run_melt", "false"), ("run_scramble", "false")):
        assert re.search(rf"^\s*{k}\s*=\s*{v}\b", conf, re.M), f"{k} default must be {v}"
    merge = open(os.path.join(REPO, "modules", "jasmine", "merge.nf")).read()
    assert "##SVcaller_callers=" in merge and "caller_set.sh" in merge
    # the effective-rule expression must exist for all three opt-in callers
    for c in ("gridss", "melt", "scramble"):
        assert f"params.run_{c} || !params.skip_{c}" in merge
