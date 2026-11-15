#!/usr/bin/env python3
"""Per-allele sizing class for ExpansionHunter calls -> <sample>.str_summary.tsv.

ExpansionHunter (EH) reads are short (150 bp).  An allele longer than one read can only be
bounded, not measured: EH still prints a point estimate (e.g. ATXN2 22/62, CI 52-102) that rests
on in-repeat and flanking reads alone.  This script classifies every allele and, for those that
cannot be sized, replaces the point estimate with a lower bound.  The EH json is never modified.

  sizing = spanned     >= MIN_SPANNING_READS spanning reads at exactly that repeat count
           flanked     flanking/in-repeat evidence only, allele fits within the read length
           unresolved  estimate longer than the read length with no spanning-read support,
                       or a CI wider than MAX_CI_WIDTH repeats

Hemizygous loci: when the sample is male (SampleParameters.Sex in the EH json), a chrX locus has ONE allele.
EH run with --sex male already reports one; if a diploid call is found (older runs that defaulted to female)
the larger allele is kept and allele2 is printed as '.'.

Compound (multi-block) loci: the catalog (assets/eh_catalog.json, built from STRchive by bin/make_eh_catalog.py) models
e.g. CNBP as (CAGG)*(CAGA)*(CA)* and EH reports one Variant per block, each counted in its OWN motif.  The disease-relevant
block is named in assets/eh_graded_variants.tsv; this script writes ONE row per locus for that block only, so a threshold
is compared in repeat units of the disease motif (CNBP: CCTG units, RFC1: AAGGG units), never the whole tract divided by
one motif.  The other blocks stay in <sample>.str_profile.json.  Column `variant_id` names the variant the row describes.

Count offsets: assets/eh_count_offsets.tsv lists loci where EH's repeat count is known to differ from the unit count the STRchive
range is written in (PABPN1: EH counts 11 for the reference's 10 GCN codons, owner decision 2026-11-15).  The offset is added to
genotype, ci, alleleN_est and the displays; eh_raw_genotype keeps EH's value and count_offset the correction.  Sizing (spanned /
flanked / unresolved) is decided on EH's raw counts because spanning-read counts are keyed by them.

Usage: eh_sizing.py <sample.str_profile.json> <out.tsv> [--sex-json <sample.sex.json>] [--graded <eh_graded_variants.tsv>] [--offsets <eh_count_offsets.tsv>]
  --sex-json  adds the inferred sex, its depth ratios and method (bin/infer_sex.py) to the header comment
  --graded    graded-variant table; default assets/eh_graded_variants.tsv beside this script's repo.  Without it every
              variant is reported, which is refused for a profile that has a multi-block locus.
"""
import json
import os
import re
import sys

# ---- constants (single place to tune) -----------------------------------------------------
MIN_SPANNING_READS = 2    # spanning reads at exactly the allele's size for "spanned"
MAX_CI_WIDTH = 30         # CI wider than this many repeats => "unresolved"
DEFAULT_READ_LEN = 150    # used only if the json carries no ReadLength
# --------------------------------------------------------------------------------------------

COLUMNS = ["locus", "repeat_unit", "genotype", "ci", "read_length",
           "allele1_est", "allele1_ci", "allele1_sizing", "allele1_display",
           "allele2_est", "allele2_ci", "allele2_sizing", "allele2_display",
           "genotype_display", "conf", "spanning_reads", "flanking_reads", "inrepeat_reads",
           "variant_id", "disease_motif",      # these two appended 2026-11: parsers that read by header name are unaffected
           "eh_raw_genotype", "count_offset"]  # appended 2026-11-15: EH's own genotype and the correction applied to it
DEFAULT_OFFSETS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "assets", "eh_count_offsets.tsv")
DEFAULT_GRADED = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "assets", "eh_graded_variants.tsv")

HEADER_NOTE = [
    "# str_summary.tsv: ExpansionHunter calls with per-allele sizing (bin/eh_sizing.py). EH json untouched.",
    f"# sizing: spanned (>= {MIN_SPANNING_READS} spanning reads at that size) | flanked (flanking reads only, fits in read) | "
    f"unresolved (estimate > read length with no spanning read, or CI width > {MAX_CI_WIDTH}).",
    "# unresolved alleles are shown as '>=<CI lower bound>' (alleleN_display / genotype_display), never a point estimate;",
    "# conf=unresolved when any allele is unresolved, else high|low (CI width of the larger allele <= max(6, 15%) = high).",
    "# alleles longer than the read length cannot be sized from short reads; they are reported as a lower bound.",
    "# count_offset: a cited per-locus correction (assets/eh_count_offsets.tsv) is ADDED to every EH count before grading; genotype, ci,",
    "# alleleN_est and *_display carry the corrected count, eh_raw_genotype keeps EH's own; spanning/flanking/inrepeat read counts stay in EH units.",
    "# one row per locus: the graded variant (variant_id) in repeat units of disease_motif (assets/eh_graded_variants.tsv).",
]


def load_graded(path):
    """{locus: {variant_id, disease_motif}} from assets/eh_graded_variants.tsv; a missing or malformed file raises."""
    rows = [l.rstrip("\n").split("\t") for l in open(path) if l.strip() and not l.startswith("#")]
    if len(rows) < 2 or not {"locus", "variant_id", "disease_motif"} <= set(rows[0]):
        raise SystemExit(f"ERROR: {path} is empty or lacks locus/variant_id/disease_motif columns")
    return {d["locus"]: d for d in (dict(zip(rows[0], r)) for r in rows[1:])}


def load_offsets(path):
    """{(locus, variant_id): int offset} from assets/eh_count_offsets.tsv; a missing or malformed file raises."""
    rows = [l.rstrip("\n").split("\t") for l in open(path) if l.strip() and not l.startswith("#")]
    if not rows or not {"locus", "variant_id", "offset", "evidence"} <= set(rows[0]):
        raise SystemExit(f"ERROR: {path} lacks locus/variant_id/offset/evidence columns")
    out = {}
    for r in rows[1:]:
        d = dict(zip(rows[0], r))
        if not d.get("evidence", "").strip():
            raise SystemExit(f"ERROR: {path}: offset for {d.get('locus')} has no cited evidence")
        out[(d["locus"], d["variant_id"])] = int(d["offset"])
    return out


def parse_counts(s):
    """'(14, 4), (15, 1)' -> {14: 4, 15: 1}; '()' -> {}."""
    return {int(a): int(b) for a, b in re.findall(r"\((\d+),\s*(\d+)\)", s or "")}


def parse_ci(s, n):
    """'11-14/14-42' -> [(11,14),(14,42)]; falls back to (None,None) per allele."""
    out = []
    segs = (s or "").split("/")
    for i in range(n):
        try:
            lo, hi = segs[i].split("-")
            out.append((int(lo), int(hi)))
        except (IndexError, ValueError):
            out.append((None, None))
    return out


def size_allele(est, ci, motif_len, read_len, spanning):
    """Return 'spanned' | 'flanked' | 'unresolved' for one allele."""
    if spanning.get(est, 0) >= MIN_SPANNING_READS:
        return "spanned"
    if ci[0] is not None and (ci[1] - ci[0]) > MAX_CI_WIDTH:
        return "unresolved"
    if est * motif_len > read_len:
        return "unresolved"
    return "flanked"


def classify_variant(vv, read_len, hemizygous=False):
    """Rows' worth of fields for one EH Variants entry; None if no genotype."""
    gt = vv.get("Genotype")
    if not gt:
        return None
    try:
        ests = [int(x) for x in gt.split("/")]
    except ValueError:
        return None
    if hemizygous and len(ests) == 2:                 # male chrX: one allele; keep the larger
        keep = max(range(2), key=lambda i: ests[i])
        ci_all = (vv.get("GenotypeConfidenceInterval", "") or "").split("/")
        vv = dict(vv, GenotypeConfidenceInterval=ci_all[keep] if len(ci_all) == 2 else "")
        ests = [ests[keep]]
    motif_len = len(vv.get("RepeatUnit", "")) or 1
    cis = parse_ci(vv.get("GenotypeConfidenceInterval", ""), len(ests))
    span = parse_counts(vv.get("CountsOfSpanningReads"))
    al = []
    for e, ci in zip(ests, cis):
        cls = size_allele(e, ci, motif_len, read_len, span)
        if cls == "unresolved":
            disp = f">={ci[0] if ci[0] is not None else e}"
        else:
            disp = str(e)
        al.append(dict(est=e, ci=ci, sizing=cls, display=disp))
    mi = max(range(len(al)), key=lambda i: al[i]["est"])
    if any(a["sizing"] == "unresolved" for a in al):
        conf = "unresolved"
    else:
        lo, hi = al[mi]["ci"]
        w = (hi - lo) if lo is not None else 0
        conf = "high" if w <= max(6, 0.15 * al[mi]["est"]) else "low"
    return dict(alleles=al, conf=conf, genotype_display="/".join(a["display"] for a in al))


def apply_offset(r, offset):
    """Shift every allele's count (and CI, and a '>=N' display) by `offset`; sizing was decided on EH's raw counts
    (spanning reads are keyed by the raw count).  A shifted count below 0 is refused."""
    if not offset:
        return r
    for a in r["alleles"]:
        if a["est"] + offset < 0:
            raise SystemExit(f"ERROR: count offset {offset} makes allele {a['est']} negative")
        a["est"] += offset
        a["ci"] = tuple(None if x is None else x + offset for x in a["ci"])
        if a["sizing"] == "unresolved":
            a["display"] = ">=" + str((a["ci"][0] if a["ci"][0] is not None else a["est"]))
        else:
            a["display"] = str(a["est"])
    r["genotype_display"] = "/".join(a["display"] for a in r["alleles"])
    return r


def summarise(profile, graded=None, offsets=None):
    """Yield one dict per locus (column name -> value).  graded = load_graded() result: only the graded variant of a
    listed locus is reported, and a profile that lacks it is refused (stale profile or catalog)."""
    male = str(profile.get("SampleParameters", {}).get("Sex", "")).lower() == "male"
    for loc, ld in profile["LocusResults"].items():
        read_len = int(ld.get("ReadLength") or DEFAULT_READ_LEN)
        variants = ld.get("Variants", {})
        g = (graded or {}).get(loc)
        if g and g["variant_id"] not in variants:
            raise SystemExit(f"ERROR: locus {loc}: graded variant {g['variant_id']} is not in the profile "
                             f"(variants {sorted(variants)}); the profile was made with a different catalog than "
                             "assets/eh_graded_variants.tsv")
        if g is None and len(variants) > 1:
            raise SystemExit(f"ERROR: locus {loc} has {len(variants)} variants and is not in the graded-variant table; "
                             "which block is the disease repeat cannot be guessed")
        for vid, vv in variants.items():
            if g and vid != g["variant_id"]:
                continue
            hemi = male and vv.get("ReferenceRegion", "").split(":")[0] in ("chrX", "X")
            r = classify_variant(vv, read_len, hemi)
            if r is None:
                continue
            raw_gt = "/".join(str(a["est"]) for a in r["alleles"])
            off = (offsets or {}).get((loc, vid), 0)
            r = apply_offset(r, off)
            row = {"locus": loc, "repeat_unit": vv.get("RepeatUnit", ""),
                   "genotype": "/".join(str(a["est"]) for a in r["alleles"]),
                   "ci": ("/".join("%s-%s" % a["ci"] for a in r["alleles"]) if (hemi or (off and all(x is not None for a in r["alleles"] for x in a["ci"])))
                          else vv.get("GenotypeConfidenceInterval", "")),
                   "read_length": read_len, "genotype_display": r["genotype_display"], "conf": r["conf"],
                   "spanning_reads": vv.get("CountsOfSpanningReads", ""),
                   "flanking_reads": vv.get("CountsOfFlankingReads", ""),
                   "inrepeat_reads": vv.get("CountsOfInrepeatReads", ""),
                   "variant_id": vid, "disease_motif": g["disease_motif"] if g else vv.get("RepeatUnit", ""),
                   "eh_raw_genotype": raw_gt, "count_offset": off}
            for i in (1, 2):
                if i <= len(r["alleles"]):
                    a = r["alleles"][i - 1]
                    lo, hi = a["ci"]
                    row.update({f"allele{i}_est": a["est"],
                                f"allele{i}_ci": f"{lo}-{hi}" if lo is not None else ".",
                                f"allele{i}_sizing": a["sizing"], f"allele{i}_display": a["display"]})
                else:
                    row.update({f"allele{i}_{k}": "." for k in ("est", "ci", "sizing", "display")})
            yield row


def main(argv):
    sex_json, graded_path, offsets_path = None, None, None
    argv = list(argv)
    for flag in ("--sex-json", "--graded", "--offsets"):
        if flag in argv:
            i = argv.index(flag)
            if i + 1 >= len(argv):
                sys.exit(__doc__)
            val = argv[i + 1]
            del argv[i:i + 2]
            if flag == "--sex-json":
                sex_json = val
            elif flag == "--offsets":
                offsets_path = val
            else:
                graded_path = val
    if len(argv) != 3:
        sys.exit(__doc__)
    if graded_path is None and os.path.exists(DEFAULT_GRADED):
        graded_path = DEFAULT_GRADED
    graded = load_graded(graded_path) if graded_path else None
    if offsets_path is None and os.path.exists(DEFAULT_OFFSETS):
        offsets_path = DEFAULT_OFFSETS
    offsets = load_offsets(offsets_path) if offsets_path else {}
    with open(argv[1]) as fh:
        profile = json.load(fh)
    if not profile.get("LocusResults"):
        sys.exit(f"ERROR: {argv[1]} has no LocusResults; refusing to write a summary")
    rows = list(summarise(profile, graded, offsets))
    with open(argv[2], "w") as out:
        out.write("\n".join(HEADER_NOTE) + "\n")
        eh_sex = profile.get("SampleParameters", {}).get("Sex", "unknown")
        out.write(f"# sex passed to ExpansionHunter (json SampleParameters.Sex): {eh_sex}; "
                  "male chrX loci are hemizygous (allele2 = '.').\n")
        if sex_json:
            sx = json.load(open(sex_json))
            out.write(f"# sex inferred: {sx['inferred_sex']} (chrX/auto {sx['chrX_ratio']}, chrY/auto {sx['chrY_ratio']}); "
                      f"method: {sx['method']}\n")
            if sx.get("sex_flag"):
                out.write(f"# FLAG: {sx['sex_flag']}\n")
        out.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            out.write("\t".join(str(r[c]) for c in COLUMNS) + "\n")


if __name__ == "__main__":
    main(sys.argv)
