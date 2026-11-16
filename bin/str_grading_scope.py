#!/usr/bin/env python3
"""Which repeat-expansion calls the HTML report may grade, and how an unresolved allele is worded.

The report must never contradict OmniGen's, so these rules are the same ones (owner decisions 2026-11-15) and read the same
curation tables, carried as copies with a comparison test (tests/test_str_curation_scope.py):

  assets/str_retired_loci.tsv      repeat conditions with no curated gene-condition link: listed, never graded
  assets/str_count_convention.tsv  narrow-gap loci are graded only with a cited counting convention (new catalogue)
  assets/str_not_tested_legacy.tsv loci whose OLD-catalogue coordinates hold no repeat in GRCh38: not tested, never normal

An allele ExpansionHunter could only bound (str_summary sizing "unresolved", shown ">=N") is "unresolved, at least N repeats",
never intermediate or pathogenic; ">=0" has no bound at all and reads "not measurable".

A call whose only intermediate or pathogenic-range allele eh_sizing marked `low_support` (a short allele with under 25% of the
locus's spanning reads, owner decision 2026-11-16) reads "not confirmed (low read support)", never a flag; the supported alleles are graded on their own (html_report.parse_str_loci).
"""
import os
import re

ASSETS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "assets")
RETIRED_FILE = os.path.join(ASSETS, "str_retired_loci.tsv")
CONVENTION_FILE = os.path.join(ASSETS, "str_count_convention.tsv")
NOT_TESTED_FILE = os.path.join(ASSETS, "str_not_tested_legacy.tsv")

RETIRED_REASON = "no curated gene-condition link (ClinGen, GenCC, PanelApp)"
CONVENTION_REASON = "counting convention not confirmed by a cited source"
NOT_MEASURABLE = "not measurable"
# Owner decision 2026-11-15 (policy threshold, not derived): a locus whose normal maximum and pathogenic minimum differ by this many
# units or fewer is graded only when str_count_convention.tsv cites how its units are counted (same constant as OmniGen str_scope).
NARROW_GAP_MAX = 2


def _table(path, need):
    lines = [l for l in open(path).read().splitlines() if l.strip() and not l.startswith("#")]
    hdr = lines[0].split("\t") if lines else []
    if len(lines) < 2 or not set(need) <= set(hdr):
        raise ValueError(f"{path} is empty or lacks columns {list(need)}")
    return [dict(zip(hdr, l.split("\t"))) for l in lines[1:]]


def retired(path=None) -> dict:
    """{locus: reason}: audit verdict RETIRE. A missing or malformed table raises: silence would grade the locus."""
    return {r["gene"]: RETIRED_REASON for r in _table(path or RETIRED_FILE, ["gene", "verdict"]) if r["verdict"] == "RETIRE"}


def convention_loci(path=None) -> set:
    """Loci with a cited counting convention (row needs a convention, a citation and a date)."""
    rows = _table(path or CONVENTION_FILE, ["locus", "convention", "citation", "verified_date"])
    bad = [r["locus"] for r in rows if not all((r[k] or "").strip() for k in ("convention", "citation", "verified_date"))]
    if bad:
        raise ValueError(f"{path or CONVENTION_FILE}: incomplete row for {bad}")
    return {r["locus"] for r in rows}


def narrow_gap(ranges: dict) -> set:
    """Loci whose normal maximum and pathogenic minimum are NARROW_GAP_MAX or fewer apart. ranges = {locus: row with normal_max,
    path_min}, computed from the table in use, never a gene list; a locus with no normal maximum has no gap to measure."""
    out = set()
    for k, v in ranges.items():
        nm, pm = str(v.get("normal_max", "")).strip(), str(v.get("path_min", "")).strip()
        if nm.isdigit() and pm.isdigit() and int(pm) - int(nm) <= NARROW_GAP_MAX:
            out.add(k)
    return out


def not_tested_legacy(path=None) -> dict:
    """{locus: reason} for loci of the OLD catalogue whose coordinates hold no repeat in GRCh38."""
    out = {}
    for r in _table(path or NOT_TESTED_FILE, ["locus", "coordinate_build"]):
        if r["coordinate_build"] not in ("hg19", "neither"):
            raise ValueError(f"{path or NOT_TESTED_FILE}: bad row {r}")
        where = "they match the hg19 position" if r["coordinate_build"] == "hg19" else "they match neither build"
        out[r["locus"]] = f"the catalog's coordinates hold no repeat in GRCh38 ({where}); not tested"
    return out


def load_scope(ranges: dict, new_catalogue: bool, retired_path=None, convention_path=None, not_tested_path=None) -> dict:
    """Read the three curation tables once (paths default to the carried assets) and compute the narrow-gap set once.
    Pass the result to not_graded_reason for every locus of a report."""
    return {"retired": retired(retired_path), "conv": convention_loci(convention_path), "nt": not_tested_legacy(not_tested_path),
            "narrow": narrow_gap(ranges) if new_catalogue else set(), "new": new_catalogue}


def not_graded_reason(locus: str, scope: dict):
    """Why a locus is listed but not graded, or None when it may be graded. Same order as OmniGen: retired, then (old catalogue)
    not tested, then (new catalogue, narrow gap) convention not confirmed. `scope` = load_scope(...)."""
    if locus in scope["retired"]:
        return scope["retired"][locus]
    if not scope["new"] and locus in scope["nt"]:
        return scope["nt"][locus]
    if scope["new"] and locus in scope["narrow"] and locus not in scope["conv"]:
        return CONVENTION_REASON
    return None


def bound_text(lb) -> str:
    """'at least N repeats', or 'not measurable' when there is no bound (>=0: no spanning or flanking read gave any size)."""
    return NOT_MEASURABLE if not lb else f"at least {lb} repeats"


def unresolved_bound(display: str):
    """Lower bound of an allele shown '>=N' (None when the text is not an unresolved bound)."""
    m = re.fullmatch(r"\s*>=\s*(\d+)\s*", display or "")
    return int(m.group(1)) if m else None


LOW_SUPPORT_STATUS = "NOT CONFIRMED"
LOW_SUPPORT_TEXT = "not confirmed (low read support)"

