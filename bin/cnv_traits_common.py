#!/usr/bin/env python3
"""Shared helpers for the CNV / blood-group copy-number trait interpreters.

Consumed by: rh_status.py, amy1_cn.py, gst_null.py, lpa_kiv2.py.

Primary signal for every trait is a targeted, normalized read-depth estimate
computed from a mosdepth `--by cnv_trait_regions.bed` output. The CNV consensus
BED (chrom,start,end,cn,svtype,caller_support,confidence,quality,sample) is used
only as a corroborating signal — empirically it under-calls the homozygous
whole-gene deletions (RHD/GSTM1/GSTT1) these traits need, so depth leads.

GRCh38 / hg38 (chr-prefixed) coordinates. Verified within the contig bounds of
hg38.canonical.fa: chr1=248,956,422; chr6=170,805,979; chr22=50,818,468.
Coordinates should be biologically re-confirmed against the exact RefSeq/GENCODE
annotation for the reference in params.ref_fasta before production use.
"""
import gzip
import statistics
from pathlib import Path

# label -> (chrom, start, end)
LOCI = {
    "RHD":          ("chr1",  25272393,  25330445),   # RHD gene, 1p36.11 (minus strand)
    "AMY1_CLUSTER": ("chr1",  103655000, 103760000),  # AMY1A/B/C salivary amylase array
    "GSTM1":        ("chr1",  109687814, 109693020),  # GSTM1, 1p13.3
    "GSTT1":        ("chr22", 24376133,  24384680),   # GSTT1, 22q11.23
    "LPA_KIV2":     ("chr6",  160605000, 160650000),  # LPA KIV-2 VNTR block within LPA
}

_SENTINELS = {"NO_FILE", "NO_DEPTH", "NO_CNV", ""}


def _open(path):
    p = str(path)
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def read_region_depths(path):
    """Parse a mosdepth `--by regions.bed(.gz)` file into {label: mean_depth}.

    mosdepth region output columns are: chrom start end [name] mean_depth.
    When the input BED carried a 4th name column, that name is preserved and used
    as the label; otherwise the region is keyed by 'chrom:start-end'.
    """
    depths = {}
    with _open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5:
                label, mean = f[3], f[4]
            elif len(f) == 4:
                label, mean = "{}:{}-{}".format(f[0], f[1], f[2]), f[3]
            else:
                continue
            try:
                depths[label] = float(mean)
            except ValueError:
                continue
    return depths


def control_baseline(depths):
    """Median mean-depth across the CTRL_* control regions.

    Control regions are diploid, copy-number-stable single-locus windows, so their
    median depth represents the depth of a normal 2-copy (2n) region. Median (not
    mean) is used so a stray non-stable control cannot skew the baseline.
    """
    ctrl = [v for k, v in depths.items() if k.upper().startswith("CTRL")]
    if not ctrl:
        raise ValueError("no CTRL_* control regions found in depth file")
    return statistics.median(ctrl)


def depth_ratio(region_mean, baseline):
    """region_depth / diploid-baseline. ~1.0 for a normal 2-copy region, ~0 for a
    homozygous deletion, >1 for amplification. Returns None if baseline invalid."""
    if baseline is None or baseline <= 0:
        return None
    return region_mean / baseline


def estimate_copies(region_mean, baseline, ploidy=2):
    """Normalized diploid copy number = ploidy * region_depth / diploid-baseline.

    baseline is the depth of a 2-copy control region, so a region at the same depth
    yields `ploidy` (=2) copies. Returns None if the baseline is invalid.
    """
    r = depth_ratio(region_mean, baseline)
    return None if r is None else ploidy * r


# --------------------------------------------------------------------------- #
# CNV consensus corroboration
# --------------------------------------------------------------------------- #
def parse_cnv_bed(path):
    """Parse the consensus CNV BED. Returns [] for missing / sentinel inputs."""
    rows = []
    if not path:
        return rows
    p = str(path)
    if Path(p).name in _SENTINELS or not Path(p).exists():
        return rows
    with open(p) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            try:
                rows.append({
                    "chrom": f[0],
                    "start": int(f[1]),
                    "end": int(f[2]),
                    "cn": f[3],
                    "svtype": f[4],
                    "caller_support": f[5] if len(f) > 5 else "",
                    "confidence": f[6] if len(f) > 6 else "",
                })
            except ValueError:
                continue
    return rows


def consensus_overlaps(cnv_rows, chrom, start, end, svtype=None):
    """Consensus rows overlapping [chrom:start-end], optionally filtered by svtype."""
    hits = []
    for r in cnv_rows:
        if r["chrom"] != chrom:
            continue
        if r["end"] < start or r["start"] > end:
            continue
        if svtype and r["svtype"].upper() != svtype.upper():
            continue
        hits.append(r)
    return hits


def write_tsv(out_path, header_cols, data_row):
    """Write a one-record contract TSV: '#'+tab-joined header, then the data row."""
    with open(out_path, "w") as fh:
        fh.write("#" + "\t".join(header_cols) + "\n")
        fh.write("\t".join(str(x) for x in data_row) + "\n")
