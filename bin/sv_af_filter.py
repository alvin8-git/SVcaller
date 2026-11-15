#!/usr/bin/env python3
"""Drop common SVs from an AnnotSV TSV by the population-AF column that matches the SV type.

AnnotSV (3.4.x, AnnotSV-benignsv.tcl, proc benignSVannotation) fills B_gain_*, B_loss_*,
B_ins_* and B_inv_* for EVERY SV: each one is the max AF of benign records of that kind
(AF >= -benignAF, default 0.01) that COMPLETELY contain the SV (bedtools -f 1.0). So
B_gain_AFmax on a deletion only says "a common gain region contains this deletion", which is
not evidence the deletion is common. The column that applies is the one of the SV's own
type: DEL->B_loss, DUP->B_gain, INS->B_ins, INV->B_inv.  Types with no such column
(TRA/BND/...) have no population AF here and are kept.

usage: sv_af_filter.py <annotated.tsv> <out.tsv> <af_threshold>
"""
import csv
import sys

AF_COLUMN_BY_TYPE = {
    "DEL": "B_loss_AFmax",
    "DUP": "B_gain_AFmax",
    "INS": "B_ins_AFmax",
    "INV": "B_inv_AFmax",
}


def row_af(row):
    col = AF_COLUMN_BY_TYPE.get((row.get("SV_type") or "").upper())
    if not col:
        return 0.0
    try:
        return float((row.get(col) or "0").strip() or "0")
    except ValueError:
        return 0.0


def main(argv):
    src, dst, thr = argv[1], argv[2], float(argv[3])
    with open(src) as fh, open(dst, "w") as out:
        reader = csv.DictReader(fh, delimiter="\t")
        writer = csv.DictWriter(out, fieldnames=reader.fieldnames or [], delimiter="\t")
        writer.writeheader()
        for row in reader:
            if row_af(row) < thr:
                writer.writerow(row)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
