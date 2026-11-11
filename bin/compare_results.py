#!/usr/bin/env python3
"""Compare the final outputs of two SVcaller result directories; exit 1 on any difference.

Used to prove that a pipeline change (scheduling, sharding, resources) did not change a
single call. Text outputs are compared byte for byte after decompression; JSON by parsed
value; VCFs on header (minus volatile lines: file date, tool command lines) and records.

Usage: compare_results.py <results_dir_A> <results_dir_B> <sample>
  <results_dir_*> is the --outdir of a run (it contains <sample>/ or the files directly).
"""
import gzip
import json
import os
import re
import sys

VOLATILE = re.compile(r"^##(fileDate|reference=|bcftools_\w+(Version|Command)|.*[Cc]ommand|source=.*bcftools)")

# (relative path under <outdir>/<sample>, kind)
OUTPUTS = [
    ("{s}.filtered.tsv", "text"),
    ("{s}.smn.tsv", "text"),
    ("{s}.str_profile.json", "json"),
    ("{s}.str_summary.tsv", "text"),
    ("{s}.str.vcf.gz", "vcf"),
    ("{s}.sv_merged.vcf.gz", "vcf"),
    ("{s}.sv_callers.json", "json"),
    ("{s}.cnv_consensus.bed", "text"),
    ("{s}.trait_depth.regions.bed.gz", "text"),
    ("cnv_traits/{s}.amy1.tsv", "text"),
    ("cnv_traits/{s}.gst_null.tsv", "text"),
    ("cnv_traits/{s}.lpa_kiv2.tsv", "text"),
    ("bloodgroup/{s}.rh_status.tsv", "text"),
    ("alpha_globin/{s}.alpha_globin.tsv", "text"),
    ("alpha_globin/{s}.alpha_depth.tsv", "text"),
    ("alpha_globin/{s}.alpha_junction.tsv", "text"),
    ("alpha_globin/{s}.alpha_sites.tsv", "text"),
    ("qc/{s}.flagstat.txt", "text"),
    ("qc/{s}.mosdepth.summary.txt", "text"),
    ("qc/{s}.insert_size_metrics.txt", "text_nohdr"),
]


def read(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", errors="replace") as fh:
        return fh.read().splitlines()


def locate(root, sample, rel):
    for base in (os.path.join(root, sample), root):
        p = os.path.join(base, rel.format(s=sample))
        if os.path.exists(p):
            return p
    return None


def first_diff(a, b):
    for i, (x, y) in enumerate(zip(a, b)):
        if x != y:
            return f"line {i + 1}: {x[:120]!r} != {y[:120]!r}"
    return f"length {len(a)} != {len(b)}"


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    ra, rb, sample = sys.argv[1:]
    bad = 0
    for rel, kind in OUTPUTS:
        pa, pb = locate(ra, sample, rel), locate(rb, sample, rel)
        name = rel.format(s=sample)
        if pa is None and pb is None:
            print(f"SKIP       {name} (absent in both)")
            continue
        if pa is None or pb is None:
            print(f"DIFFERENT  {name}: present only in {'B' if pa is None else 'A'}")
            bad += 1
            continue
        if kind == "json":
            same = json.load(open(pa)) == json.load(open(pb))
            detail = "parsed JSON differs"
        else:
            la, lb = read(pa), read(pb)
            if kind == "vcf":
                la = [l for l in la if not VOLATILE.match(l)]
                lb = [l for l in lb if not VOLATILE.match(l)]
            if kind == "text_nohdr":      # Picard writes a date + command header
                la = [l for l in la if not l.startswith("# ") and not l.startswith("## ")]
                lb = [l for l in lb if not l.startswith("# ") and not l.startswith("## ")]
            same = la == lb
            detail = first_diff(la, lb) if not same else ""
        print(f"{'IDENTICAL ' if same else 'DIFFERENT '} {name}" + ("" if same else f": {detail}"))
        bad += not same
    print("RESULT:", "IDENTICAL" if not bad else f"{bad} output(s) differ")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
