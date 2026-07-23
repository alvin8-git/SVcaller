"""The documented SV ensemble size must match what the merge actually merges.

WHY. For months the README, CLAUDE.md and VERSION.md advertised a "6-caller
ensemble (… + SvABA)". It was never true: `modules/jasmine/merge.nf` decompresses
`vcfs[0..4]` and never references `vcfs[5]`, so SvABA's VCF was staged into the
merge task directory and ignored. Every published benchmark, including the run16
F1 baseline, was a 5-caller result labelled as 6.

That is not a cosmetic error. It propagated: OmniGen's generated condition list
told users their structural variants came from "Manta/DELLY/GRIDSS/MELT/SvABA",
naming a caller that contributed nothing and omitting one that did.

These tests are bidirectional on purpose:
  - a doc claiming 6 fails while the code merges 5;
  - the code starting to merge a 6th fails while the docs still say 5.
Either way somebody has to reconcile them, which is the point.
"""
import os
import re

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MERGE_NF = os.path.join(REPO, "modules", "jasmine", "merge.nf")

# Files that discuss the historical error and must be allowed to quote it.
# Keep this list SHORT — it is an exemption from a correctness check.
HISTORY_OK = {
    "docs/reorg-plan-2026-07-21.md",   # the doc that first identified the error
    "docs/CHANGES.md",                 # dated record of the fix
    "tests/test_ensemble_caller_count.py",
    "CLAUDE.md",                       # states the correction explicitly, in context
}

WRONG_COUNT = re.compile(r"\b6[- ]caller|\bsix[- ]caller|\b6 callers\b", re.I)


def _merge_src():
    assert os.path.exists(MERGE_NF), f"{MERGE_NF} missing"
    return open(MERGE_NF).read()


def merged_vcf_indices():
    """Which vcfs[N] the merge module actually touches."""
    return sorted({int(m) for m in re.findall(r"vcfs\[(\d+)\]", _merge_src())})


def test_merge_handles_exactly_five_vcfs():
    idx = merged_vcf_indices()
    assert idx == [0, 1, 2, 3, 4], (
        f"merge.nf references vcfs{idx}. If a 6th caller was wired in, update "
        "README.md, CLAUDE.md, VERSION.md and OmniGen's condition list, then "
        "change this test deliberately — do not just widen it.")


def test_svaba_is_not_merged():
    """The specific claim. vcfs[5] is SvABA in subworkflows/sv_calling.nf."""
    assert "vcfs[5]" not in _merge_src(), (
        "merge.nf now references vcfs[5] (SvABA). If that is intentional, the "
        "ensemble is 6 callers: update the docs and flip skip_svaba back to false, "
        "because it currently defaults to true precisely BECAUSE it was not merged.")


def test_skip_svaba_default_matches_the_wiring():
    """A caller that is not merged must not run by default, and vice versa.

    These two facts drifted apart for months: skip_svaba defaulted to false while
    the merge ignored its output, so every run paid ~4 h/sample at 16 pinned CPUs
    for nothing.
    """
    conf = open(os.path.join(REPO, "nextflow.config")).read()
    m = re.search(r"^\s*skip_svaba\s*=\s*(true|false)", conf, re.M)
    assert m, "skip_svaba not found in nextflow.config"
    skipped = m.group(1) == "true"
    merged = "vcfs[5]" in _merge_src()
    assert skipped != merged, (
        f"skip_svaba={'true' if skipped else 'false'} but SvABA is "
        f"{'merged' if merged else 'NOT merged'} — running a caller whose output "
        "is discarded, or discarding a caller you paid to run.")


def _docs():
    out = []
    for root, dirs, files in os.walk(REPO):
        dirs[:] = [d for d in dirs
                   if not d.startswith(".") and d not in
                   ("work", "results", "__pycache__", "node_modules")
                   and not d.startswith(("work_", "results_"))]
        for f in files:
            if f.endswith(".md"):
                out.append(os.path.join(root, f))
    return sorted(out)


def test_no_doc_claims_six_callers():
    offenders = []
    for path in _docs():
        rel = os.path.relpath(path, REPO)
        if rel in HISTORY_OK:
            continue
        for i, line in enumerate(open(path, errors="replace"), 1):
            if WRONG_COUNT.search(line):
                offenders.append(f"{rel}:{i}: {line.strip()[:100]}")
    assert not offenders, (
        "these claim a 6-caller ensemble; the merge takes 5 (SvABA is never "
        "merged):\n  " + "\n  ".join(offenders))


def test_report_does_not_name_svaba_while_it_is_unmerged():
    """The user-facing HTML report must not credit SvABA while merge.nf ignores it.

    This is not hypothetical: the report template's caller footer read
    "Manta + Delly + GRIDSS + Scramble + MELT + SvABA + STRling" long after SvABA
    was known to contribute nothing, because the doc-scan above only walks .md
    files and never the template or the report builder. A clinician reading the
    report would credit a caller that produced no calls.

    Tied to the merge wiring so it inverts correctly the day SvABA is really
    wired in (vcfs[5]) — then the report SHOULD name it and this test flips.
    """
    if "vcfs[5]" in _merge_src():
        return  # SvABA re-enabled: naming it in the report is now correct
    targets = [
        os.path.join(REPO, "assets", "report_template.html"),
        os.path.join(REPO, "bin", "html_report.py"),
    ]
    named = []
    for p in targets:
        src = open(p, errors="replace").read()
        # "not merged" next to the name is the deliberate historical note, allowed.
        if re.search(r"svaba", src, re.I) and "not merged" not in src.lower():
            named.append(os.path.relpath(p, REPO))
    assert not named, (
        f"{named} name SvABA but merge.nf does not merge it (vcfs[5] absent). "
        "Remove it from the caller list, or wire vcfs[5] in first.")


@pytest.mark.parametrize("doc", ["README.md", "CLAUDE.md"])
def test_key_docs_state_the_real_count(doc):
    """Silence is not enough — the primary docs must say 5 explicitly, since a
    consumer like OmniGen reads them to describe the pipeline to end users."""
    src = open(os.path.join(REPO, doc), errors="replace").read()
    assert re.search(r"\b5[- ]caller|\bfive[- ]caller", src, re.I), \
        f"{doc} does not state the 5-caller ensemble anywhere"
