"""Guard: the two STR callers in sv_calling.nf run on the SAME BAM.

WHY. EXPANSIONHUNTER used to take the raw `ch_bam` while STRling took
`ch_validated_bam`. For a full-hg38 BAM input the raw BAM carries 3366 contigs and
reads whose mate is on an alt contig; the canonical-filtered BAM drops them. So
the two STR callers genotyped different read sets and fed a mismatched STR
consensus, and EH alone skipped the VALIDATE_REF_BAM pre-flight. This keeps them
consistent, and consistent with every SV caller.
"""
import os
import re

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SV_NF = os.path.join(REPO, "subworkflows", "sv_calling.nf")


def _first_arg(src, process):
    m = re.search(process + r"\s*\(\s*(\w+)", src)
    assert m, f"{process} call not found in sv_calling.nf"
    return m.group(1)


@pytest.fixture
def src():
    return open(SV_NF).read()


def test_expansionhunter_uses_the_validated_bam(src):
    got = _first_arg(src, "EXPANSIONHUNTER")
    assert got == "ch_validated_bam", (
        f"EXPANSIONHUNTER runs on {got}; it must run on ch_validated_bam so it sees "
        "the same canonical-filtered, validated BAM as STRling and the SV callers.")


def test_both_str_callers_see_the_same_bam(src):
    assert _first_arg(src, "EXPANSIONHUNTER") == _first_arg(src, "STRLING_CALL"), (
        "EXPANSIONHUNTER and STRLING_CALL are both STR callers feeding one STR "
        "consensus; they must genotype the same BAM.")


def test_sv_calling_has_no_raw_bam_reference(src):
    """FILTER_CHROMS + VALIDATE_REF_BAM were lifted to svcaller.nf. SV_CALLING now
    receives an already-validated BAM, so `ch_bam` must not appear in its code at
    all (comments aside); every caller uses ch_validated_bam."""
    code = "\n".join(re.sub(r"//.*$", "", ln) for ln in src.splitlines())
    assert "ch_bam" not in code, \
        "sv_calling.nf still references ch_bam; it should take ch_validated_bam only"


# --- top-level BAM routing (svcaller.nf) --------------------------------------

SVCALLER_NF = os.path.join(REPO, "workflows", "svcaller.nf")


@pytest.fixture
def top():
    return open(SVCALLER_NF).read()


@pytest.mark.parametrize("sub", ["SV_CALLING", "CNV_CALLING", "SMN_CALLING", "CNV_TRAITS"])
def test_non_alpha_modules_run_on_the_validated_bam(top, sub):
    """The filter is lifted here so CNV/SMN/traits see the same validated BAM as
    the SV callers, not the raw full-hg38 BAM (whose alt contigs bias CNVpytor)."""
    assert _first_arg(top, sub) == "ch_validated_bam", (
        f"{sub} runs on {_first_arg(top, sub)}; it must run on ch_validated_bam")


def test_alpha_globin_deliberately_runs_on_the_raw_bam(top):
    """Alpha-globin is the ONE module kept on the raw ch_bam, on purpose: its GIAB
    depth baselines were calibrated on the raw BAM, so querying the filtered BAM
    against them would drift score = ratio / baseline. Moving it needs the baselines
    re-derived on filtered GIAB BAMs first. If this ever changes, that re-derivation
    must land in the same change."""
    assert _first_arg(top, "ALPHA_GLOBIN") == "ch_bam", (
        "ALPHA_GLOBIN must stay on the raw ch_bam until its baselines are re-derived "
        "on filtered BAMs; see validation/giab_alpha_baseline.tsv")
