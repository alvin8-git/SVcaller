#!/usr/bin/env python3
"""Compare ExpansionHunter repeat sizes with GIAB HG002 assembly-derived truth.

Usage: eh_vs_q100.py <str_profile.json> <truth.vcf.gz> <benchmark.bed> <out.tsv>

Truth allele length per haplotype = reference repeat count (ReferenceRegion length / motif length)
plus the net indel length inside the region on that haplotype / motif length, rounded.  Only
indels count (SNVs do not change length); a deletion/insertion is clipped to the part inside the
region.  A locus whose region is not wholly inside the benchmark bed has "no truth".
Reference repeat count: EH's own REF from <prefix>.str.vcf.gz when present, else region length / motif.
Agreement: both | one | none  (an EH allele agrees when within +-TOL repeats of a truth allele,
best pairing) | no_truth.  EH alleles use the point estimate EH prints, so an "unresolved"
allele is judged on its (unreliable) estimate; its sizing class is shown alongside.
Also reports truth_normal (truth max allele < path_min in assets/str_disease_ranges.tsv, i.e. not in
the pathogenic range; n/a when the locus has no path_min).
"""
import csv
import importlib.util
import itertools
import json
import os
import sys

import pysam

TOL = 1
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
CATALOG = os.path.join(REPO, "assets", "eh_catalog.json")
RANGES = os.path.join(REPO, "assets", "str_disease_ranges.tsv")

_spec = importlib.util.spec_from_file_location("eh_sizing", os.path.join(REPO, "bin", "eh_sizing.py"))
eh_sizing = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(eh_sizing)


def load_bed(path):
    iv = {}
    with open(path) as fh:
        for ln in fh:
            f = ln.split()
            if len(f) >= 3:
                iv.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    return iv


def in_bed(iv, chrom, start0, end):
    return any(s <= start0 and end <= e for s, e in iv.get(chrom, []))


def truth_alleles(vcf, chrom, start1, end1, motif_len, ref_units=None):
    """Return (list of per-haplotype repeat counts, list of per-haplotype bp deltas)."""
    r0, r1 = start1 - 1, end1                      # 0-based half-open region
    delta = [0, 0]
    for r in vcf.fetch(chrom, max(0, r0 - 1), r1 + 1):
        if r.alts is None:
            continue
        gt = r.samples[0]["GT"]
        haps = [(0, gt[0]), (1, gt[0])] if len(gt) == 1 else list(enumerate(gt[:2]))  # hemizygous X/Y
        for hap, g in haps:
            if not g:
                continue
            alt, ref = r.alts[g - 1], r.ref
            if len(alt) == len(ref) or alt.startswith("<"):
                continue
            if len(alt) > len(ref):                 # insertion after anchor base at r.start
                if r0 <= r.start < r1:
                    delta[hap] += len(alt) - len(ref)
            else:                                   # deletion of ref[1+len(alt)-1:], starts at r.start+len(alt)
                d0, d1 = r.start + len(alt), r.stop
                delta[hap] -= max(0, min(d1, r1) - max(d0, r0))
    if ref_units is None:
        ref_units = (r1 - r0) / motif_len
    return [round(ref_units + d / motif_len) for d in delta], delta


def eh_vcf_refs(profile_path):
    """EH's own reference repeat count (INFO REF) per locus from <prefix>.str.vcf.gz next to the json.
    ReferenceRegion/motif over-counts interrupted loci (e.g. ATXN1 CAT interruptions); EH's REF is the
    pure-repeat count EH itself uses, so truth and EH are on the same scale. {} if no VCF is found."""
    p = profile_path.replace(".str_profile.json", ".str.vcf.gz")
    refs = {}
    if p != profile_path and os.path.exists(p):
        for r in pysam.VariantFile(p):
            if "REF" in r.info and "VARID" in r.info:
                refs[r.info["VARID"]] = int(r.info["REF"])
    return refs


def agreement(eh, truth):
    best = 0
    for perm in itertools.permutations(truth, len(eh)) if len(eh) <= len(truth) else [tuple(truth)]:
        best = max(best, sum(abs(a - b) <= TOL for a, b in zip(eh, perm)))
    n = len(eh)
    return "both" if best == n else ("one" if best >= 1 else "none")


def main(argv):
    if len(argv) != 5:
        sys.exit(__doc__)
    prof, vcf_path, bed_path, out = argv[1:]
    profile = json.load(open(prof))
    cat = {c["LocusId"]: c for c in json.load(open(CATALOG))}
    pmin = {}   # first pathogenic repeat count per locus
    with open(RANGES) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                pmin[row["locus_id"]] = int(row["path_min"])
            except ValueError:
                pass
    refs = eh_vcf_refs(prof)
    bed, vcf = load_bed(bed_path), pysam.VariantFile(vcf_path)
    summ = {r["locus"]: r for r in eh_sizing.summarise(profile)}
    cols = ["locus", "eh_genotype", "eh_ci", "eh_sizing", "eh_display", "truth_genotype",
            "truth_delta_bp", "agreement", "unresolved", "truth_normal"]
    rows = []
    for loc, c in cat.items():
        s = summ.get(loc)
        if s is None:
            rows.append({"locus": loc, "eh_genotype": "not_called", "agreement": "no_eh_call"})
            continue
        chrom, rng = c["ReferenceRegion"].split(":")
        a, b = (int(x) for x in rng.split("-"))
        motif = len(s["repeat_unit"]) or 1
        eh = [int(x) for x in s["genotype"].split("/")]
        sizing = "/".join(s[f"allele{i}_sizing"] for i in (1, 2) if s[f"allele{i}_sizing"] != ".")
        row = {"locus": loc, "eh_genotype": s["genotype"], "eh_ci": s["ci"], "eh_sizing": sizing,
               "eh_display": s["genotype_display"], "unresolved": "yes" if "unresolved" in sizing else "no"}
        if in_bed(bed, chrom, a - 1, b):
            t, d = truth_alleles(vcf, chrom, a, b, motif, refs.get(loc))
            row.update(truth_genotype="/".join(map(str, sorted(t))), truth_delta_bp="/".join(map(str, d)),
                       agreement=agreement(sorted(eh), sorted(t)),
                       truth_normal=("yes" if max(t) < pmin[loc] else "no") if loc in pmin else "n/a")
        else:
            row.update(truth_genotype="no_truth", truth_delta_bp=".", agreement="no_truth", truth_normal="n/a")
        rows.append(row)
    with open(out, "w") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", restval=".")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main(sys.argv)
