#!/usr/bin/env python3
"""Integrate the four alpha-globin evidence channels into the frozen contract TSV.

    channel 1  hba_depth.py     -> <S>.alpha_depth.tsv      per-segment score
    channel 2  THIS FILE                                    allele naming
    channel 3  hba_junction.py  -> <S>.alpha_junction.tsv   breakpoints
    channel 4  hba_sites.py     -> <S>.alpha_sites.tsv      targeted pileup

Output: results/<S>/alpha_globin/<S>.alpha_globin.tsv, exactly the 11 columns in
docs/contracts/alpha_globin_contract.md, validated against
validation/examples/SAMPLE.alpha_globin.tsv.

THIS MODULE MEASURES. IT DOES NOT INTERPRET.
No HbH / Bart's / trait classification, no couple risk, no "carrier"/"clear"
wording, and `interpretation_complete` is structurally always `false` — see
CONTRACT_INTERPRETATION_COMPLETE below. OmniGen owns every clinical statement.

----------------------------------------------------------------------------
FOUR RULES THIS FILE EXISTS TO ENFORCE. Each was a real, costly bug.
----------------------------------------------------------------------------

1. NEVER NAME ONE ALLELE FROM A DEGENERATE GROUP.
   `--SEA|--MED` and `--FIL|--THAI` have IDENTICAL depth signatures; see the
   `depth_distinguishable` column of assets/hba_deletion_alleles.tsv. Picking
   `--SEA` because the sample looks SE Asian is a population inference dressed
   up as a measurement. We emit the group. See `_collapse_group()` for why we
   never collapse it today.

2. THRESHOLD ON `score`, NEVER ON THE RAW RATIO. Enforced upstream in
   hba_depth.py; this file consumes `call`, never `ratio`.

3. A SEGMENT FLAGGED `do_not_average` IS NOT EVIDENCE. INTER_Z_A reads 0.99
   ("intact") in a sample where a --SEA deletion covers half of it. It is
   treated as no observation at all, never as a vote for "intact".

4. ABSENT != EMPTY != POPULATED. An absent channel (Nextflow NO_* sentinel) is
   a legitimate skip and moves its tier into `not_screened`. A channel file that
   is PRESENT but has no header is a crashed caller: raise, never default to a
   reassuring 4-gene "normal" result. Mirrors bin/smn_report.py's SmnInputError.
"""
import argparse
import hashlib
import os
import sys
from pathlib import Path

SCORABLE = ("HBZ", "HBA2", "INTER_A2_A1", "HBA1")

# The contract's column list, in order. Frozen: docs/contracts/alpha_globin_contract.md.
CONTRACT_COLUMNS = [
    "sample", "alpha_genes_called", "alpha_genes_confidence", "deletion_alleles",
    "deletion_evidence", "site_genotypes", "site_panel_version", "genotype",
    "screened", "not_screened", "interpretation_complete",
]

# Not a variable. SVcaller measures; it must never be able to claim it interpreted.
CONTRACT_INTERPRETATION_COMPLETE = "false"

# Tier ids. Kept to exactly the set in the frozen example fixture — OmniGen
# renders these strings, so inventing a new id here is an interface change.
TIER_DELETIONAL = "alpha_deletional"
TIER_SITES = "alpha_targeted_sites"
NEVER_SCREENED = ["beta_globin", "alpha_nondeletional_outside_panel"]

_SENTINELS = {"NO_FILE", "NO_DEPTH", "NO_JUNCTION", "NO_SITES", ""}

# call -> copies present at that segment, out of a normal 2.
CALL_TO_COPIES = {"intact": 2, "het_loss": 1, "hom_loss": 0, "gain": 3}


class AlphaGlobinInputError(RuntimeError):
    """A channel file is PRESENT but unreadable/empty => that caller RAN and
    FAILED. Refuse to emit a contract row rather than render a crashed channel
    as a normal result. This is the SMN incident (an empty smn.tsv shown as
    "0 Carrier findings ... Clear") and it must not be repeated here."""


# --------------------------------------------------------------------------- #
# input parsing
# --------------------------------------------------------------------------- #
def is_absent(path) -> bool:
    """True for a Nextflow NO_* sentinel or a genuinely missing path."""
    if not path:
        return True
    name = Path(str(path)).name
    return name in _SENTINELS or name.startswith("NO_") or not Path(str(path)).exists()


def read_tsv(path, what):
    """Parse a '#'-header channel TSV into a list of dicts.

    Absent -> None (legitimate skip). Present-but-headerless -> raise."""
    if is_absent(path):
        return None
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]
    if not lines:
        raise AlphaGlobinInputError(
            f"{what} file '{path}' is present but empty -- that channel ran and "
            f"failed. Refusing to emit an alpha-globin result that would read as "
            f"a normal finding.")
    header = lines[0].lstrip("#").split("\t")
    return [dict(zip(header, l.split("\t"))) for l in lines[1:]]


def group_order(alleles):
    """allele name -> its row index in assets/hba_deletion_alleles.tsv.

    Group members are joined in FILE order, not alphabetically. Both the asset's
    own `depth_distinguishable` cells (`no:--SEA|--MED`) and every occurrence in
    docs/contracts/alpha_globin_contract.md write `--SEA|--MED` and
    `--FIL|--THAI`; an alphabetical sort emits `--MED|--SEA`, which is the same
    measurement spelled a way the consumer was not told to expect. OmniGen
    renders this string, so the spelling is interface, not cosmetics."""
    return {a["allele"]: i for i, a in enumerate(alleles)}


def parse_alleles(path):
    """assets/hba_deletion_alleles.tsv -> list of allele dicts."""
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip() and not l.startswith("#")]
    header = lines[0].split("\t")
    rows = [dict(zip(header, l.split("\t"))) for l in lines[1:]]
    for r in rows:
        if "d_INTER_A2_A1" not in r:
            raise AlphaGlobinInputError(
                "hba_deletion_alleles.tsv has no d_INTER_A2_A1 column. INTER_A2_A1 "
                "is the ONLY diagnostic segment for -a3.7 (the commonest alpha-thal "
                "allele worldwide); without that column -a3.7 cannot be matched. "
                "Regenerate with bin/make_hba_deletion_alleles.py.")
    return rows


def observed_copies(depth_rows):
    """{segment: copies} from channel 1, for scorable segments only.

    Segments flagged do_not_average / uncalibrated, and any segment the depth
    caller could not score, are OMITTED rather than assumed intact. Omission is
    "no observation"; it must never become a vote."""
    out = {}
    for r in depth_rows or []:
        seg = r.get("segment")
        if seg not in SCORABLE:
            continue                      # INTER_Z_A: do_not_average, rule 3
        c = CALL_TO_COPIES.get(r.get("call"))
        if c is not None:
            out[seg] = c
    return out


def marginal_segments(depth_rows):
    """Segments the depth caller flagged as sitting near a decision boundary."""
    out = []
    for r in depth_rows or []:
        if r.get("segment") in SCORABLE and str(r.get("marginal", "")).lower() in ("1", "true", "yes"):
            out.append(r["segment"])
    return out


# --------------------------------------------------------------------------- #
# channel 2 — allele naming
# --------------------------------------------------------------------------- #
def _deltas(sym):
    """A d_* cell -> the SET of copy changes this haplotype may make there.

    Three symbols, and the difference between the last two is a bug that was
    found by adversarial re-derivation after the naive version shipped:

      0 / -1 / +   exact: {0}, {-1}, {+1}.

      'h'  disrupted/hybrid -> {0}. An EXPECTATION, not a wildcard. An -a3.7
           hybrid gene fuses HBA2's 5' end to HBA1's 3' end, so neither gene
           body vanishes; with mosdepth --mapq 0 the hybrid's reads still pile
           onto both paralogues and HBA2/HBA1 land in the `intact` band. This is
           load-bearing: an earlier version skipped 'h' segments entirely, which
           made -a3.7 a near-wildcard that matched THAL1's --SEA signature
           (2,1,1,1) and named the wrong allele on the only real sample we have.

      '?'  unconstrained -> {-1, 0}. The expectation genuinely straddles a
           decision boundary. -a4.2 is ~4.2 kb while HBA2's gene body is only
           835 bp, so 2-3.4 kb of the 2969 bp INTER_A2_A1 goes with it and that
           segment scores anywhere from ~0.43 to ~0.77 depending on the 5'
           breakpoint. Asserting {0} there (what 'h' would do) misnames a silent
           3-gene -a4.2 carrier as a 2-gene alpha-thal trait. '+1' is excluded:
           a deletion allele cannot gain copies.

    UNVALIDATED: no -a3.7 or -a4.2 sample exists in this project, so both are
    reasoning from the gene model, not measurements. The failure mode is
    deliberately one-sided — when nothing matches we emit NA, and NA is safe
    where a confidently wrong allele is not."""
    s = str(sym).strip()
    if s == "h":
        return frozenset({0})
    if s == "?":
        return frozenset({-1, 0})
    if s == "+":
        return frozenset({+1})
    return frozenset({int(s)})


def _cell(hap, seg):
    """hap['d_<seg>'], with a real error instead of a raw KeyError."""
    key = "d_" + seg
    if key not in hap:
        raise AlphaGlobinInputError(
            f"allele {hap.get('allele', '?')!r} has no {key} column, so segment "
            f"{seg} cannot be matched. Either the segment is not one an allele "
            f"may be defined against (INTER_Z_A is do_not_average and "
            f"deliberately has no column), or the allele table needs "
            f"regenerating with bin/make_hba_deletion_alleles.py.")
    return hap[key]


def expected_copies(hap_a, hap_b, seg):
    """The SET of copy counts consistent with genotype hap_a/hap_b at `seg`.

    Two copies to start; each haplotype (None = wild-type 'aa') contributes one
    of its allowed deltas."""
    out = {2}
    for hap in (hap_a, hap_b):
        if hap is None:
            continue
        out = {t + d for t in out for d in _deltas(_cell(hap, seg))}
    return out


def has_hybrid(hap, seg):
    """True when this haplotype's expectation at `seg` is not a single value."""
    return hap is not None and str(_cell(hap, seg)).strip() in ("h", "?")


def match_genotypes(obs, alleles):
    """Every (hap_a, hap_b) genotype consistent with the observed segment copies.

    Wild-type is None. EVERY observed segment must be consistent — there is no
    nearest-match. A signature matching nothing yields no hits, and the caller
    reports NA rather than rounding onto the closest allele.

    Returns (name_a, name_b, alpha_genes_called, rests_on_an_inexact_expectation)."""
    cands = [None] + list(alleles)
    hits = []
    for i, a in enumerate(cands):
        for b in cands[i:]:
            if any(seen not in expected_copies(a, b, seg)
                   for seg, seen in obs.items()):
                continue
            lost = sum(int(h["alpha_genes_lost"]) for h in (a, b) if h is not None)
            hyb = any(has_hybrid(h, seg) for h in (a, b) for seg in obs)
            hits.append((a["allele"] if a else "aa",
                         b["allele"] if b else "aa",
                         4 - lost, hyb))
    return hits


def _collapse_group(candidates, junction_rows):
    """Would a junction let us name ONE allele out of a degenerate group?

    Today: NO, always. The contract permits collapsing only on a junction read
    or a measured extent that EXCLUDES the alternative, and neither is available:

      * assets/hba_deletion_alleles.tsv carries no breakpoints at all, by design
        — alpha-cluster NAHR breakpoints sit inside near-identical homology
        boxes, are published against varying builds, and are not single-valued
        for a given allele. `approx_size` is documentary and a committed test
        asserts it can never be mistaken for a coordinate.
      * So a measured extent of ~20 kb cannot be said to exclude --MED (~17 kb)
        without comparing against numbers we have deliberately refused to write
        down.

    Collapsing on anything else — population, prior, the supplier's label — is a
    population inference dressed as a measurement, and is forbidden by the brief.
    This function is the single place that would change if a per-allele,
    provenance-carrying breakpoint table ever lands. It is deliberately NOT
    wired to `junction_rows` yet; the parameter documents what the input would
    be."""
    return None


def name_alleles(obs, alleles, junction_rows=None):
    """Channel 2. -> (deletion_alleles, alpha_genes_called, note).

    deletion_alleles is '/'-separated haplotypes; '|' INSIDE a haplotype marks
    an unresolved degenerate group, exactly as the contract specifies."""
    if not obs:
        return "NA", "NA", "no scorable segment was measured"

    hits = match_genotypes(obs, alleles)
    if not hits:
        return ("NA", "NA",
                "observed segment copies match no allele or combination of "
                "alleles in the panel; reported as unresolved rather than "
                "forced onto the nearest allele")

    # PARSIMONY. Some genotypes cancel exactly at segment level: --SEA on one
    # chromosome and anti-3.7 (a triplication) on the other restores every
    # segment to 2 copies, so it is indistinguishable from plain aa/aa. Without
    # this filter EVERY normal sample would be reported as ambiguous between
    # "no deletion" and "a deletion masked by a triplication" — technically true,
    # operationally useless, and a guaranteed source of false alarms.
    #
    # We keep only the genotypes with the fewest non-wild-type haplotypes. What
    # that discards is a genuine, if remote, blind spot: a compensated
    # deletion+triplication carrier reads as normal. It is recorded in
    # docs/CHANGES.md and in the module docstring rather than in every row,
    # because a caveat printed on every sample is a caveat nobody reads.
    fewest = min(sum(1 for n in (h[0], h[1]) if n != "aa") for h in hits)
    hits = [h for h in hits if sum(1 for n in (h[0], h[1]) if n != "aa") == fewest]

    order = group_order(alleles)
    rank = lambda n: (0, order[n]) if n in order else (1, 0)   # noqa: E731
    slot_a = sorted({h[0] for h in hits}, key=rank)
    slot_b = sorted({h[1] for h in hits}, key=rank)
    genes = {h[2] for h in hits}

    collapsed = _collapse_group(slot_a, junction_rows)
    if collapsed:
        slot_a = [collapsed]

    note = ""
    if len(genes) > 1:
        note = ("more than one genotype is consistent with the depth signature "
                "and they imply different alpha-gene counts")
    elif len(slot_a) > 1 or len(slot_b) > 1:
        note = "degenerate group: depth alone cannot separate these alleles"

    haps = ["|".join(slot_a), "|".join(slot_b)]
    # Put any deletion haplotype first so 'aa' trails, matching --SEA|--MED/aa.
    haps.sort(key=lambda h: h == "aa")
    called = "none" if haps == ["aa", "aa"] else "/".join(haps)

    if len(genes) != 1:
        return called, "NA", (note + "; candidate genotypes disagree on the "
                              "alpha-gene count").lstrip("; ")
    n = genes.pop()
    if not 0 <= n <= 6:
        # RESOLVED 2026-07-22: the contract now declares 0-6, so a triplication
        # is reported as the 5 (or 6) genes it actually is. The old 0-4 range
        # assumed alpha variation only ever REMOVES genes; anti-3.7 is the
        # reciprocal product of the -a3.7 NAHR and adds one, so a carrier has 5
        # and a homozygote 6. Emitting NA for a perfectly determined count was
        # reporting a measurement failure that had not happened.
        #
        # This branch now catches only genuinely impossible counts — a bug in
        # the allele table or in expected_copies(), not a real genotype. NA here
        # means "could not be determined", which is what it should always have
        # meant.
        return called, "NA", (
            f"{n} alpha genes implied, which is outside the possible range 0-6; "
            f"this indicates a defect in the allele table or the copy-number "
            f"expectations, not a real genotype")
    return called, str(n), note


# --------------------------------------------------------------------------- #
# evidence / confidence
# --------------------------------------------------------------------------- #
def deletion_evidence(has_depth_call, junction_rows):
    """Which channels supported the call: depth|junction|both|none."""
    has_junction = bool(junction_rows)
    if has_depth_call and has_junction:
        return "both"
    if has_depth_call:
        return "depth"
    if has_junction:
        return "junction"
    return "none"


def confidence(depth_rows, junction_rows, has_depth_call, marginal):
    """high = depth and junction agree · medium = depth only · low = marginal
    or the two channels disagree."""
    if depth_rows is None:
        return "low"
    if marginal:
        return "low"
    if junction_rows is None:
        return "medium"
    if bool(junction_rows) == bool(has_depth_call):
        return "high"
    return "low"


# --------------------------------------------------------------------------- #
# channel 4 rollup
# --------------------------------------------------------------------------- #
def site_genotypes(site_rows):
    """';'-joined GENE:HGVS:zygosity for sites actually found, else 'none'.

    Only rows the site caller marked `present` are listed. A `no_call` row
    (depth too low) is NOT silently dropped into 'none' — 'none' would read as
    "we looked and it is absent", which is exactly the false-negative shape this
    module exists to remove. It downgrades confidence instead."""
    if site_rows is None:
        return "NA", False
    present = [r for r in site_rows if r.get("call") == "present"]
    nocall = any(r.get("call") == "no_call" for r in site_rows)
    if not present:
        return "none", nocall
    return ";".join(
        "{}:{}:{}".format(r.get("gene", "?"), r.get("hgvs_c", "?"),
                          r.get("zygosity", "NA"))
        for r in present), nocall


def panel_version(panel_path):
    """'hba_pathogenic_sites.tsv@<sha1[:7]>' — which panel actually ran.

    "No pathogenic site found" means nothing without knowing what was
    interrogated, so the contract pins the panel's hash."""
    if is_absent(panel_path):
        return "NA"
    h = hashlib.sha1(Path(panel_path).read_bytes()).hexdigest()[:7]
    return "{}@{}".format(Path(panel_path).name, h)


def integrated_genotype(deletion_alleles, site_rows):
    """The integrated genotype string.

    DELIBERATE DIVERGENCE FROM THE CONTRACT'S EXAMPLE. The contract illustrates
    this field as `--SEA/aQSa`, i.e. the site variant written INTO a haplotype.
    We do not do that, because it asserts phase we have not measured: short
    reads at this locus do not tell us which chromosome a site variant sits on,
    and on a deletion background that placement is the whole clinical question.
    We emit the deletion genotype, then the site findings after ' +', so the two
    measurements stay separable and neither implies the other."""
    base = deletion_alleles if deletion_alleles not in ("none", "NA") else "aa/aa"
    if deletion_alleles == "NA":
        base = "NA"
    sites = [r for r in (site_rows or []) if r.get("call") == "present"]
    if not sites:
        return base
    tail = ";".join("{}:{}:{}".format(r.get("gene", "?"), r.get("hgvs_c", "?"),
                                      r.get("zygosity", "NA")) for r in sites)
    return "{} +{}".format(base, tail)


# --------------------------------------------------------------------------- #
def build_row(sample, depth_rows, junction_rows, site_rows, alleles, panel_path):
    obs = observed_copies(depth_rows)
    marginal = marginal_segments(depth_rows)

    if depth_rows is None:
        del_alleles, genes, note = "NA", "NA", "depth channel did not run"
    else:
        del_alleles, genes, note = name_alleles(obs, alleles, junction_rows)

    has_depth_call = del_alleles not in ("none", "NA")
    sites, has_nocall = site_genotypes(site_rows)

    conf = confidence(depth_rows, junction_rows, has_depth_call, marginal)
    if has_nocall and conf == "high":
        conf = "medium"

    screened, not_screened = [], list(NEVER_SCREENED)
    (screened if depth_rows is not None else not_screened).append(TIER_DELETIONAL)
    (screened if site_rows is not None else not_screened).append(TIER_SITES)

    return {
        "sample": sample,
        "alpha_genes_called": genes,
        "alpha_genes_confidence": conf,
        "deletion_alleles": del_alleles,
        "deletion_evidence": deletion_evidence(has_depth_call, junction_rows),
        "site_genotypes": sites,
        "site_panel_version": panel_version(panel_path),
        "genotype": integrated_genotype(del_alleles, site_rows),
        "screened": ",".join(screened) if screened else "none",
        "not_screened": ",".join(not_screened),
        "interpretation_complete": CONTRACT_INTERPRETATION_COMPLETE,
    }, note


def write_contract(out_path, row):
    with open(out_path, "w") as fh:
        fh.write("\t".join(CONTRACT_COLUMNS) + "\n")
        fh.write("\t".join(str(row[c]) for c in CONTRACT_COLUMNS) + "\n")


def main(argv=None):
    here = os.path.dirname(os.path.abspath(__file__))
    default_alleles = os.path.join(here, "..", "assets", "hba_deletion_alleles.tsv")

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sample", required=True)
    ap.add_argument("--depth", help="<S>.alpha_depth.tsv (channel 1)")
    ap.add_argument("--junction", help="<S>.alpha_junction.tsv (channel 3)")
    ap.add_argument("--sites", help="<S>.alpha_sites.tsv (channel 4)")
    ap.add_argument("--alleles", default=default_alleles)
    ap.add_argument("--panel", help="assets/hba_pathogenic_sites.tsv (for the hash)")
    ap.add_argument("--out")
    ap.add_argument("--print-alpha-genes", action="store_true",
                    help="Print ONLY alpha_genes_called (from --depth) and exit. "
                         "Channel 4 needs the alpha-gene count before it can turn "
                         "a VAF into a zygosity, but that count is produced by "
                         "channel 2 -- which lives here. This mode lets the "
                         "subworkflow get it without duplicating the allele "
                         "matcher, and without a two-pass run.")
    a = ap.parse_args(argv)

    if a.print_alpha_genes:
        rows = read_tsv(a.depth, "alpha depth (channel 1)")
        if rows is None:
            print("NA")
            return 0
        _, genes, _ = name_alleles(observed_copies(rows), parse_alleles(a.alleles))
        print(genes)
        return 0

    if not a.out:
        ap.error("--out is required unless --print-alpha-genes is given")

    depth_rows = read_tsv(a.depth, "alpha depth (channel 1)")
    junction_rows = read_tsv(a.junction, "alpha junction (channel 3)")
    site_rows = read_tsv(a.sites, "alpha sites (channel 4)")

    if depth_rows is None and site_rows is None:
        raise AlphaGlobinInputError(
            "neither the depth nor the site channel produced output for "
            f"{a.sample}; there is nothing to report and an empty contract row "
            "would be read downstream as a negative result.")

    row, note = build_row(a.sample, depth_rows, junction_rows, site_rows,
                          parse_alleles(a.alleles), a.panel)
    write_contract(a.out, row)
    print("alpha_globin: {} genes={} alleles={} evidence={} conf={}".format(
        a.sample, row["alpha_genes_called"], row["deletion_alleles"],
        row["deletion_evidence"], row["alpha_genes_confidence"]))
    if note:
        print("  note: " + note)
    return 0


if __name__ == "__main__":
    sys.exit(main())
