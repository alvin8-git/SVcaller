#!/usr/bin/env python3
"""AMY1 (salivary amylase) copy number by targeted normalized read depth.

AMY1 is a high-copy tandem array (0-20+ copies) that the genome-wide SV/CNV
ensemble cannot resolve, so this is a depth-only estimate: reads over the AMY1
cluster window are summed with mosdepth (--mapq 0, multi-mapping reads retained)
and divided by the diploid control baseline.

Contract: results/<S>/cnv_traits/amy1.tsv
    columns: sample  AMY1_copies  method
"""
import argparse

import cnv_traits_common as c

METHOD = "read-depth-normalized"


def call_amy1(depths):
    try:
        baseline = c.control_baseline(depths)
    except ValueError:
        baseline = None
    region_mean = depths.get("AMY1_CLUSTER")
    if region_mean is None or not baseline:
        return "NA"
    copies = c.estimate_copies(region_mean, baseline)
    return int(round(copies)) if copies is not None else "NA"


def main():
    ap = argparse.ArgumentParser(description="AMY1 copy number")
    ap.add_argument("--depth", required=True, help="mosdepth --by regions.bed.gz")
    ap.add_argument("--cnv-bed", default=None, help="unused; accepted for a uniform interface")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    depths = c.read_region_depths(args.depth)
    amy1_copies = call_amy1(depths)

    c.write_tsv(args.out,
                ["sample", "AMY1_copies", "method"],
                [args.sample, amy1_copies, METHOD])
    print("wrote {}: AMY1_copies={} method={}".format(args.out, amy1_copies, METHOD))


if __name__ == "__main__":
    main()
