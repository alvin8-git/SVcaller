#!/usr/bin/env python3
"""Merge CNVpytor and GATK gCNV output into a consensus BED file.

Usage: cnv_consensus.py --cnvpytor <tsv> --gatk <tsv> --sample <id> --out <bed>
"""
import argparse, csv, sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional


@dataclass
class CNVSegment:
    chrom: str
    start: int
    end: int
    cn: int
    svtype: str   # DEL or DUP
    caller: str
    quality: Optional[float] = None


def reciprocal_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> float:
    overlap = max(0, min(a_end, b_end) - max(a_start, b_start))
    len_a, len_b = a_end - a_start, b_end - b_start
    if len_a == 0 or len_b == 0:
        return 0.0
    return overlap / min(len_a, len_b)


def load_cnvpytor(path: str) -> List[CNVSegment]:
    segs = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            svtype_raw = parts[0].lower()
            region = parts[1]  # e.g. chr1:1000-2000
            if ":" not in region or "-" not in region:
                continue
            chrom, coords = region.split(":")
            start, end = map(int, coords.split("-"))
            cn_raw = float(parts[3]) if len(parts) > 3 else 2.0
            cn = round(cn_raw)
            if cn == 2:          # copy-neutral — not a CNV (mirror load_gatk's cn==2 skip)
                continue
            svtype = "DEL" if svtype_raw == "deletion" or cn < 2 else "DUP"
            segs.append(CNVSegment(chrom, start, end, cn, svtype, "CNVpytor"))
    return segs


def load_gatk(path: str) -> List[CNVSegment]:
    """Load GATK segments from the converted TSV (CONTIG/START/END/CALL_COPY_NUMBER/QUALITY).

    Fails loud, not silent: if every data row fails to parse (e.g. the raw .seg was
    passed instead of the converted .tsv, so CALL_COPY_NUMBER is absent), warn to
    stderr rather than returning an empty list that looks like a true negative.
    """
    segs = []
    n_data = 0
    n_parse_fail = 0
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            n_data += 1
            try:
                chrom = row["CONTIG"]
                start = int(row["START"])
                end   = int(row["END"])
                cn    = int(row["CALL_COPY_NUMBER"])
                qual  = float(row.get("QUALITY", 0))
            except (KeyError, ValueError):
                n_parse_fail += 1
                continue
            if cn == 2:          # diploid / neutral call — not a CNV
                continue
            svtype = "DEL" if cn < 2 else "DUP"
            segs.append(CNVSegment(chrom, start, end, cn, svtype, "GATK_gCNV", qual))
    if n_data > 0 and n_parse_fail == n_data:
        sys.stderr.write(
            f"WARNING: cnv_consensus parsed 0/{n_data} GATK rows from {path} "
            "(missing CALL_COPY_NUMBER column — was the raw .seg passed instead of "
            "the converted .tsv?). CNV consensus will be empty.\n")
    return segs


def merge(cnvpytor: List[CNVSegment], gatk: List[CNVSegment],
          min_reciprocal: float = 0.5, gatk_qual_threshold: float = 30.0,
          cnvpytor_only_min_bp: int = 1_000_000) -> List[dict]:
    """Three confidence tiers:
      BOTH/HIGH        — CNVpytor + GATK agree (reciprocal overlap >= min_reciprocal)
      GATK_only/MEDIUM — GATK call passing the quality threshold, no CNVpytor match
      CNVpytor_only/LOW — large (>= cnvpytor_only_min_bp) CNVpytor call, no GATK match.
                          Size-gated so sub-clinical read-depth noise doesn't flood the report.
    """
    results = []
    matched_gatk = set()
    matched_cnvpytor = set()

    for j, a in enumerate(cnvpytor):
        best_match = None
        best_overlap = 0.0
        for i, b in enumerate(gatk):
            if a.chrom != b.chrom or a.svtype != b.svtype:
                continue
            ovl = reciprocal_overlap(a.start, a.end, b.start, b.end)
            if ovl >= min_reciprocal and ovl > best_overlap:
                best_overlap = ovl
                best_match = (i, b)
        if best_match:
            idx, b = best_match
            matched_gatk.add(idx)
            matched_cnvpytor.add(j)
            results.append({
                "chrom": a.chrom, "start": a.start, "end": a.end,
                "cn": b.cn, "svtype": a.svtype,
                "caller_support": "BOTH", "confidence": "HIGH",
                "quality": b.quality if b.quality is not None else "."
            })

    for i, b in enumerate(gatk):
        if i in matched_gatk:
            continue
        if (b.quality or 0) >= gatk_qual_threshold:
            results.append({
                "chrom": b.chrom, "start": b.start, "end": b.end,
                "cn": b.cn, "svtype": b.svtype,
                "caller_support": "GATK_only", "confidence": "MEDIUM",
                "quality": b.quality
            })

    for j, a in enumerate(cnvpytor):
        if j in matched_cnvpytor:
            continue
        if (a.end - a.start) < cnvpytor_only_min_bp:
            continue
        results.append({
            "chrom": a.chrom, "start": a.start, "end": a.end,
            "cn": a.cn, "svtype": a.svtype,
            "caller_support": "CNVpytor_only", "confidence": "LOW",
            "quality": a.quality if a.quality is not None else "."
        })

    results.sort(key=lambda r: (r["chrom"], r["start"]))
    return results


def main():
    parser = argparse.ArgumentParser(description="Merge CNVpytor and GATK gCNV calls")
    parser.add_argument("--cnvpytor", required=True)
    parser.add_argument("--gatk",     required=True)
    parser.add_argument("--sample",   required=True)
    parser.add_argument("--out",      required=True)
    args = parser.parse_args()

    cnvpytor_segs = load_cnvpytor(args.cnvpytor)
    gatk_segs     = load_gatk(args.gatk)
    consensus     = merge(cnvpytor_segs, gatk_segs)

    with open(args.out, "w") as fh:
        fh.write("#chrom\tstart\tend\tcn\tsvtype\tcaller_support\tconfidence\tquality\tsample\n")
        if not consensus:
            # Self-documenting empty result: a clinician/auditor can tell an empty
            # sheet apart from a parse failure. Comment line — not a data row.
            fh.write(f"# 0 consensus CNV segments "
                     f"(CNVpytor input={len(cnvpytor_segs)}, GATK input={len(gatk_segs)})\n")
        for r in consensus:
            fh.write(f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r['cn']}\t"
                     f"{r['svtype']}\t{r['caller_support']}\t{r['confidence']}\t"
                     f"{r['quality']}\t{args.sample}\n")

    print(f"Written {len(consensus)} consensus CNV segments to {args.out}")


if __name__ == "__main__":
    main()
