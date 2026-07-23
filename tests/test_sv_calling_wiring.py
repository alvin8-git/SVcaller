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


def test_raw_ch_bam_only_feeds_the_filter_branch(src):
    """After the FILTER_CHROMS branch, no caller should take the raw ch_bam.
    ch_bam is legitimate only in the `.branch { needs_filter ... }` that decides
    what to filter; everything downstream runs on ch_validated_bam."""
    # strip comments so a `// ... ch_bam ...` note does not trip the check
    code = "\n".join(re.sub(r"//.*$", "", ln) for ln in src.splitlines())
    # every ch_bam occurrence must be the take decl or the branch, never CALLER(ch_bam
    offenders = re.findall(r"([A-Z][A-Z0-9_]*)\s*\(\s*ch_bam\b", code)
    assert not offenders, (
        f"these run on the raw ch_bam instead of ch_validated_bam: {offenders}")
