#!/usr/bin/env python3
"""LPA KIV-2 tandem-repeat copy number (Lp(a) / cardiovascular risk).

The LPA kringle-IV type-2 (KIV-2) VNTR varies from ~5 to ~50 copies per allele.
All copies are near-identical, so mosdepth over the single reference KIV-2 window
(--mapq 0, multi-mapping reads retained) accumulates coverage proportional to the
total copy number; dividing by the diploid single-copy control baseline yields the
total (diploid) KIV-2 copy estimate.

Contract: results/<S>/cnv_traits/lpa_kiv2.tsv
    columns: sample  KIV2_copies  method

NOTE: the absolute scaling depends on how many KIV-2 units the reference window
already spans; the estimate must be calibrated against a truth sample before
clinical interpretation (tracked in docs/omnigen-additions-plan.md).
"""
import argparse

import cnv_traits_common as c

METHOD = "read-depth-ratio"


def call_kiv2(depths):
    try:
        baseline = c.control_baseline(depths)
    except ValueError:
        baseline = None
    region_mean = depths.get("LPA_KIV2")
    if region_mean is None or not baseline:
        return "NA"
    copies = c.estimate_copies(region_mean, baseline)
    return int(round(copies)) if copies is not None else "NA"


def main():
    ap = argparse.ArgumentParser(description="LPA KIV-2 repeat copy number")
    ap.add_argument("--depth", required=True, help="mosdepth --by regions.bed.gz")
    ap.add_argument("--cnv-bed", default=None, help="unused; accepted for a uniform interface")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    depths = c.read_region_depths(args.depth)
    kiv2_copies = call_kiv2(depths)

    c.write_tsv(args.out,
                ["sample", "KIV2_copies", "method"],
                [args.sample, kiv2_copies, METHOD])
    print("wrote {}: KIV2_copies={} method={}".format(args.out, kiv2_copies, METHOD))


if __name__ == "__main__":
    main()
