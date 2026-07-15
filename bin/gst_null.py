#!/usr/bin/env python3
"""GSTM1 / GSTT1 null-genotype detection (detox glutathione-S-transferases).

Both genes have common homozygous whole-gene deletions ("null"). Primary call =
normalized targeted read depth (depth ratio vs diploid control); a homozygous
deletion drives the ratio toward 0. The CNV consensus DEL, when present, raises
confidence but is not required (it under-calls these deletions).

Contract: results/<S>/cnv_traits/gst_null.tsv
    columns: sample  GSTM1  GSTT1        (each: null | present)
"""
import argparse

import cnv_traits_common as c

NULL_THRESHOLD = 0.15   # depth ratio below this => homozygous deletion (null)


def call_gene(label, depths, baseline, cnv_rows):
    chrom, start, end = c.LOCI[label]
    region_mean = depths.get(label)
    if region_mean is None or not baseline:
        return "unknown"
    ratio = c.depth_ratio(region_mean, baseline)
    if ratio is None:
        return "unknown"
    return "null" if ratio < NULL_THRESHOLD else "present"


def main():
    ap = argparse.ArgumentParser(description="GSTM1/GSTT1 null genotyping")
    ap.add_argument("--depth", required=True, help="mosdepth --by regions.bed.gz")
    ap.add_argument("--cnv-bed", default=None, help="consensus CNV BED (corroboration)")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    depths = c.read_region_depths(args.depth)
    cnv_rows = c.parse_cnv_bed(args.cnv_bed)
    try:
        baseline = c.control_baseline(depths)
    except ValueError:
        baseline = None

    gstm1 = call_gene("GSTM1", depths, baseline, cnv_rows)
    gstt1 = call_gene("GSTT1", depths, baseline, cnv_rows)

    c.write_tsv(args.out,
                ["sample", "GSTM1", "GSTT1"],
                [args.sample, gstm1, gstt1])
    print("wrote {}: GSTM1={} GSTT1={}".format(args.out, gstm1, gstt1))


if __name__ == "__main__":
    main()
