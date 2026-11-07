#!/usr/bin/env python3
"""Append caller-provenance columns to the AF-filtered AnnotSV table.  No row is dropped.

Columns appended (after the AnnotSV ones; DEPTH_* are appended after these by depth_check.py):
  CALLERS_RAN        callers in SUPP_VEC bit order, e.g. "Manta,Delly"  (from the merged VCF header
                     ##SVcaller_callers=; "legacy" when the VCF predates it)
  SUPP_VEC_WIDTH     number of SUPP_VEC bits = len(CALLERS_RAN)  ("." for legacy)
  SUPP_CALLERS       callers whose SUPP_VEC bit is 1 for this row, e.g. "Manta,Delly"
  CALLER_CONSENSUS   PASS | SINGLE_CALLER | NA
      only Manta+Delly in the vector : PASS iff BOTH bits are 1
      any other caller set           : PASS iff >=2 bits are 1 (the historical SUPP>=2 rule)
      no SUPP_VEC / legacy header    : NA
SINGLE_CALLER rows stay in the table; they are flagged, not removed.

Usage: caller_consensus.py <in.tsv> <sv_merged.vcf[.gz]> <out.tsv>
"""
import csv
import gzip
import re
import sys

NEW_COLS = ["CALLERS_RAN", "SUPP_VEC_WIDTH", "SUPP_CALLERS", "CALLER_CONSENSUS"]


def callers_from_vcf(path):
    """Caller names in SUPP_VEC order from ##SVcaller_callers=, or None (legacy VCF)."""
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                return None
            if line.startswith("##SVcaller_callers="):
                v = line.rstrip("\n").split("=", 1)[1].strip()
                return [x for x in v.split(",") if x]
    return None


def classify(callers, info):
    """-> (supp_callers_str, consensus) for one row's INFO string."""
    m = re.search(r"(?:^|;)SUPP_VEC=([01]+)", info or "")
    if not callers or not m or len(m.group(1)) != len(callers):
        return ".", "NA"
    bits = m.group(1)
    sup = [c for c, b in zip(callers, bits) if b == "1"]
    if set(callers) == {"Manta", "Delly"}:
        ok = len(sup) == 2
    else:
        ok = len(sup) >= 2
    return ",".join(sup) or ".", "PASS" if ok else "SINGLE_CALLER"


def main(tsv, vcf, out):
    callers = callers_from_vcf(vcf)
    ran = ",".join(callers) if callers else "legacy"
    width = str(len(callers)) if callers else "."
    with open(tsv) as fh, open(out, "w", newline="") as o:
        rd = csv.DictReader(fh, delimiter="\t")
        fields = list(rd.fieldnames or [])
        wr = csv.DictWriter(o, fieldnames=fields + [c for c in NEW_COLS if c not in fields],
                            delimiter="\t", lineterminator="\n")
        wr.writeheader()
        for row in rd:
            sc, cons = classify(callers, row.get("INFO", ""))
            row.update(CALLERS_RAN=ran, SUPP_VEC_WIDTH=width, SUPP_CALLERS=sc, CALLER_CONSENSUS=cons)
            wr.writerow(row)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
