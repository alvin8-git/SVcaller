#!/usr/bin/env python3
"""Rh factor (Rh+/-) from RHD copy number.

The common Rh-negative haplotype is a homozygous whole-gene deletion of RHD.
Primary call = normalized targeted read depth over the RHD locus; the CNV
consensus is used only to corroborate a called deletion.

Contract: results/<S>/bloodgroup/rh_status.tsv
    columns: sample  RHD_copies  Rh_status  confidence

Note: RHD and its paralog RHCE are near-identical; residual RHCE cross-mapping can
lift a true 0-copy toward ~0.3, so RHD_copies < 0.5 is treated as a deletion.
"""
import argparse

import cnv_traits_common as c

DEL_THRESHOLD = 0.5   # normalized copies below this => RHD deleted (Rh-negative)


def call_rh(depths, cnv_rows):
    """Return (RHD_copies:int|'NA', Rh_status, confidence)."""
    chrom, start, end = c.LOCI["RHD"]
    try:
        baseline = c.control_baseline(depths)
    except ValueError:
        baseline = None

    region_mean = depths.get("RHD")
    copies = None
    if region_mean is not None and baseline:
        copies = c.estimate_copies(region_mean, baseline)

    del_hit = bool(c.consensus_overlaps(cnv_rows, chrom, start, end, svtype="DEL"))

    if copies is None:
        return "NA", "unknown", "LOW"

    rhd_copies = int(round(copies))
    if copies < DEL_THRESHOLD:
        rh_status = "neg"
        confidence = "HIGH" if del_hit else "MEDIUM"
    else:
        rh_status = "pos"
        # Clear diploid presence is high confidence; a single ambiguous copy is medium.
        confidence = "HIGH" if copies >= 1.5 else "MEDIUM"
    return rhd_copies, rh_status, confidence


def main():
    ap = argparse.ArgumentParser(description="Rh factor / RHD copy number")
    ap.add_argument("--depth", required=True, help="mosdepth --by regions.bed.gz")
    ap.add_argument("--cnv-bed", default=None, help="consensus CNV BED (corroboration)")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    depths = c.read_region_depths(args.depth)
    cnv_rows = c.parse_cnv_bed(args.cnv_bed)
    rhd_copies, rh_status, confidence = call_rh(depths, cnv_rows)

    c.write_tsv(args.out,
                ["sample", "RHD_copies", "Rh_status", "confidence"],
                [args.sample, rhd_copies, rh_status, confidence])
    print("wrote {}: RHD_copies={} Rh_status={} confidence={}".format(
        args.out, rhd_copies, rh_status, confidence))


if __name__ == "__main__":
    main()
