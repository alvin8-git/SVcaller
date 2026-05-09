#!/usr/bin/env python3
"""Validate and emit samplesheet rows as JSON lines for Nextflow."""
import csv, json, sys
from pathlib import Path

def validate(row: dict) -> dict:
    sid = row.get("sample", "").strip()
    if not sid:
        raise ValueError("'sample' column is required and must not be empty")
    fq1 = row.get("fastq_1", "").strip()
    fq2 = row.get("fastq_2", "").strip()
    bam = row.get("bam", "").strip()
    if bam and (fq1 or fq2):
        raise ValueError(f"Sample {sid}: provide fastq_1/fastq_2 OR bam, not both")
    if (fq1 and not fq2) or (fq2 and not fq1):
        raise ValueError(f"Sample {sid}: fastq_1 and fastq_2 must both be provided")
    if not bam and not fq1:
        raise ValueError(f"Sample {sid}: must provide either bam or fastq_1/fastq_2")
    entry = {"id": sid, "single_end": False}
    if bam:
        if not Path(bam).exists():
            raise FileNotFoundError(f"BAM not found: {bam}")
        entry["bam"] = bam
        entry["input_type"] = "bam"
    else:
        for f in [fq1, fq2]:
            if not Path(f).exists():
                raise FileNotFoundError(f"FASTQ not found: {f}")
        entry["fastq_1"] = fq1
        entry["fastq_2"] = fq2
        entry["input_type"] = "fastq"
    return entry

def main():
    if len(sys.argv) < 2:
        print("Usage: parse_samplesheet.py <samplesheet.csv>", file=sys.stderr)
        sys.exit(1)
    path = sys.argv[1]
    if not Path(path).exists():
        print(f"ERROR: Samplesheet not found: {path}", file=sys.stderr)
        sys.exit(1)
    with open(path) as fh:
        for row in csv.DictReader(fh):
            try:
                print(json.dumps(validate(row)))
            except (ValueError, FileNotFoundError) as e:
                print(f"ERROR: {e}", file=sys.stderr)
                sys.exit(1)

if __name__ == "__main__":
    main()
