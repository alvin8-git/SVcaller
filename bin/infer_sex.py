#!/usr/bin/env python3
"""Infer genetic sex from read depth in a few fixed windows (samtools coverage; seconds on a 30x BAM).

Same rule and thresholds as OmniGen prototype/vcf_qc.py:
  chrY(X-degenerate window)/autosome < 0.10 and chrX/autosome >= 0.75   -> female
  chrY/autosome >= 0.25 and chrX/autosome < 0.75                        -> male
  anything else                                                         -> unknown
Unknown is run through ExpansionHunter as female (its default) with `sex_flag` set; it is never silent.

Usage: infer_sex.py <bam> <sample_id> <out.json> [<out.sex.txt>]
  out.json  {sample, inferred_sex, eh_sex, sex_flag, chrX_ratio, chrY_ratio, depths, windows, method}
  out.txt   one word, the value passed to `ExpansionHunter --sex` (male|female)
"""
import json
import subprocess
import sys

# ---- constants -------------------------------------------------------------------------------
WIN_Y = (6500000, 10500000)                       # chrY X-degenerate region (GRCh38)
WIN_X = (40000000, 44000000)                      # chrX non-PAR
WIN_AUTO = [("1", 20000000, 24000000), ("2", 20000000, 24000000), ("3", 20000000, 24000000)]
Y_FEMALE_MAX, X_FEMALE_MIN = 0.10, 0.75
Y_MALE_MIN, X_MALE_MAX = 0.25, 0.75
MIN_AUTOSOME_DEPTH = 5.0                          # below this the ratios are noise -> unknown
# ----------------------------------------------------------------------------------------------
METHOD = ("samtools coverage mean depth: chrY %d-%d, chrX %d-%d vs mean of chr1/2/3 20-24 Mb; "
          "female if Y/auto<%.2f and X/auto>=%.2f, male if Y/auto>=%.2f and X/auto<%.2f, else unknown"
          % (WIN_Y + WIN_X + (Y_FEMALE_MAX, X_FEMALE_MIN, Y_MALE_MIN, X_MALE_MAX)))


def decide(y_ratio, x_ratio):
    if y_ratio < Y_FEMALE_MAX and x_ratio >= X_FEMALE_MIN:
        return "female"
    if y_ratio >= Y_MALE_MIN and x_ratio < X_MALE_MAX:
        return "male"
    return "unknown"


def mean_depth(bam, chrom, a, b, prefix):
    out = subprocess.run(["samtools", "coverage", "-r", f"{prefix}{chrom}:{a}-{b}", bam],
                         capture_output=True, text=True, check=True).stdout.strip().splitlines()
    return float(out[-1].split("\t")[6])


def infer(bam):
    hdr = subprocess.run(["samtools", "view", "-H", bam], capture_output=True, text=True, check=True).stdout
    prefix = "chr" if "SN:chr1\t" in hdr else ""
    dy = mean_depth(bam, "Y", *WIN_Y, prefix)
    dx = mean_depth(bam, "X", *WIN_X, prefix)
    da = sum(mean_depth(bam, c, a, b, prefix) for c, a, b in WIN_AUTO) / len(WIN_AUTO)
    if da < MIN_AUTOSOME_DEPTH:
        return dict(inferred_sex="unknown", chrX_ratio=None, chrY_ratio=None,
                    depths=dict(chrX=dx, chrY=dy, autosome=da))
    yr, xr = dy / da, dx / da
    return dict(inferred_sex=decide(yr, xr), chrX_ratio=round(xr, 3), chrY_ratio=round(yr, 3),
                depths=dict(chrX=round(dx, 2), chrY=round(dy, 2), autosome=round(da, 2)))


def main(argv):
    if len(argv) not in (4, 5):
        sys.exit(__doc__)
    bam, sample, out = argv[1:4]
    r = infer(bam)
    r.update(sample=sample, eh_sex="male" if r["inferred_sex"] == "male" else "female",
             sex_flag=("sex undetermined from depth; ExpansionHunter run as female (chrX loci diploid)"
                       if r["inferred_sex"] == "unknown" else ""),
             method=METHOD)
    with open(out, "w") as fh:
        json.dump(r, fh, indent=1)
    if len(argv) == 5:
        with open(argv[4], "w") as fh:
            fh.write(r["eh_sex"] + "\n")


if __name__ == "__main__":
    main(sys.argv)
