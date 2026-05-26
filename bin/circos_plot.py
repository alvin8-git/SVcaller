#!/usr/bin/env python3
"""Generate a Circos plot from SV VCF and CNV BED using pycirclize.

Usage:
  circos_plot.py --sv-vcf merged.vcf.gz --cnv-bed cnv.bed \
                 --cytobands GRCh38_cytobands.txt \
                 --sample SAMPLE_ID --out circos.svg
"""
import argparse, re, gzip
from pathlib import Path
from typing import List, Tuple, Dict, Optional


CHROM_ORDER = [f"chr{i}" for i in list(range(1, 23)) + ["X", "Y"]]

SV_COLOURS = {
    "DEL": "#1F77B4",
    "DUP": "#D62728",
    "INV": "#9467BD",
    "BND": "#FF7F0E",
    "TRA": "#FF7F0E",
    "INS": "#2CA02C",
}


def sv_colour(svtype: str) -> str:
    return SV_COLOURS.get(svtype.upper(), "#7F7F7F")


def parse_cnv_bed(path: str) -> Tuple[List[dict], List[dict]]:
    gains, losses = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            rec = {"chrom": parts[0], "start": int(parts[1]), "end": int(parts[2]),
                   "cn": int(parts[3]), "svtype": parts[4]}
            if rec["svtype"] == "DUP":
                gains.append(rec)
            elif rec["svtype"] == "DEL":
                losses.append(rec)
    return gains, losses


def parse_str_vcf(path: Optional[str]) -> List[dict]:
    """Parse ExpansionHunter VCF; return positions of all genotyped STR loci."""
    loci = []
    if not path or path in ("NO_STR", ""):
        return loci
    opener = gzip.open if str(path).endswith(".gz") else open
    try:
        with opener(path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 8:
                    continue
                chrom, pos = parts[0], int(parts[1])
                if chrom not in CHROM_ORDER:
                    continue
                loci.append({"chrom": chrom, "pos": pos})
    except Exception:
        pass
    return loci


def parse_sv_vcf_links(path: str,
                       min_svlen_intra: int = 50_000,
                       max_links: int = 100) -> List[dict]:
    """Return SV links filtered for clinical Circos display.

    Filtering rationale (ACMG/ClinGen germline SV practice):
    - BND/TRA: all inter-chromosomal rearrangements (always significant if confirmed)
    - DEL/DUP/INV: >= min_svlen_intra (default 50 kb; smaller are typically benign)
    - INS/MEI: excluded (< 1 kb insertions are unreadable in Circos; table is better)
    - Multi-caller only (SUPP ones >= 2) to reduce false-positive visual noise
    - Capped at max_links total, prioritised by SVLEN descending
    """
    all_links = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            chrom1, pos1, info = parts[0], int(parts[1]), parts[7]
            svtype_m = re.search(r"SVTYPE=([^;]+)", info)
            if not svtype_m:
                continue
            svtype = svtype_m.group(1).upper()

            # Skip MEI insertions — too small for meaningful Circos display
            if svtype == "INS":
                continue

            # Multi-caller filter: require SUPP_VEC with >= 2 supporting callers
            supp_m = re.search(r"SUPP_VEC=([01]+)", info)
            if supp_m:
                ones = supp_m.group(1).count("1")
                if ones < 2:
                    continue

            svlen_m = re.search(r"SVLEN=(-?\d+)", info)
            svlen = abs(int(svlen_m.group(1))) if svlen_m else 0
            end_m = re.search(r"END=(\d+)", info)

            if svtype in ("BND", "TRA"):
                alt = parts[4]
                mate_m = re.search(r"[\[\]]([^:[\]]+):(\d+)[\[\]]", alt)
                if not mate_m:
                    continue
                chrom2, pos2 = mate_m.group(1), int(mate_m.group(2))
                if chrom2 == chrom1:
                    continue  # intra-chromosomal BND: treat as intra, apply size filter
            else:
                if svlen < min_svlen_intra:
                    continue
                chrom2 = chrom1
                pos2 = int(end_m.group(1)) if end_m else pos1 + svlen

            if chrom1 not in CHROM_ORDER or chrom2 not in CHROM_ORDER:
                continue
            all_links.append({
                "chrom1": chrom1, "pos1": pos1,
                "chrom2": chrom2, "pos2": pos2,
                "svtype": svtype, "svlen": svlen,
                "colour": sv_colour(svtype),
            })

    # Prioritise by SVLEN descending; BND/TRA (svlen=0) get lowest priority by size
    # but keep all BND/TRA up to the cap
    bnd = [l for l in all_links if l["svtype"] in ("BND", "TRA")]
    intra = sorted([l for l in all_links if l["svtype"] not in ("BND", "TRA")],
                   key=lambda x: x["svlen"], reverse=True)
    # Fill cap: BND first (up to half), then largest intra-chromosomal SVs
    max_bnd = min(len(bnd), max_links // 2)
    links = bnd[:max_bnd] + intra[:max_links - max_bnd]
    return links


def load_chrom_sizes(cytobands_path: str) -> Dict[str, int]:
    sizes = {}
    with open(cytobands_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            chrom, end = parts[0], int(parts[2])
            sizes[chrom] = max(sizes.get(chrom, 0), end)
    return {c: sizes[c] for c in CHROM_ORDER if c in sizes}


def _chrom_colour(chrom: str) -> str:
    idx = CHROM_ORDER.index(chrom) if chrom in CHROM_ORDER else 0
    palette = [
        "#1f77b4","#aec7e8","#ffbb78","#2ca02c","#98df8a","#d62728","#ff9896",
        "#9467bd","#c5b0d5","#8c564b","#c49c94","#e377c2","#f7b6d2","#7f7f7f",
        "#c7c7c7","#bcbd22","#dbdb8d","#17becf","#9edae5","#393b79","#5254a3",
        "#6b6ecf","#9c9ede","#637939",
    ]
    return palette[idx % len(palette)]


def make_circos(sv_vcf: str, cnv_bed: str, cytobands: str,
                sample_id: str, out_svg: str, out_png: str,
                str_vcf: Optional[str] = None) -> None:
    from pycirclize import Circos
    import matplotlib.pyplot as plt

    chrom_sizes = load_chrom_sizes(cytobands)
    gains, losses = parse_cnv_bed(cnv_bed)
    links = parse_sv_vcf_links(sv_vcf)  # filtered: BND/TRA + DEL/DUP/INV>=50kb, multi-caller, cap 100
    str_loci = parse_str_vcf(str_vcf)
    print(f"Circos: {len(links)} SV links selected for display")

    circos = Circos(chrom_sizes, space=1.5)
    circos.text(f"SVcaller\n{sample_id}", size=10, r=15)

    # Ring 1: chromosome ideograms
    for sector in circos.sectors:
        track = sector.add_track((95, 100))
        track.axis(fc=_chrom_colour(sector.name))
        track.text(sector.name.replace("chr", ""), size=6, color="white")

    # Ring 2: CNV gains (red)
    for sector in circos.sectors:
        track = sector.add_track((80, 93), r_pad_ratio=0.1)
        track.axis()
        sector_gains = [(g["start"], g["end"], g["cn"] - 2)
                        for g in gains if g["chrom"] == sector.name]
        for start, end, height in sector_gains:
            if height > 0:
                track.rect(start, end, fc="#D62728", alpha=0.7)

    # Ring 3: CNV losses (blue)
    for sector in circos.sectors:
        track = sector.add_track((67, 80), r_pad_ratio=0.1)
        track.axis()
        sector_losses = [(l["start"], l["end"]) for l in losses if l["chrom"] == sector.name]
        for start, end in sector_losses:
            track.rect(start, end, fc="#1F77B4", alpha=0.7)

    # Ring 4: STR expansion markers (brown dots per locus)
    for sector in circos.sectors:
        str_track = sector.add_track((55, 65), r_pad_ratio=0.1)
        str_track.axis()
        for locus in str_loci:
            if locus["chrom"] == sector.name:
                chrom_len = chrom_sizes.get(sector.name, 1)
                marker_width = max(500_000, chrom_len // 500)
                str_track.rect(locus["pos"], locus["pos"] + marker_width,
                               fc="#8C564B", alpha=0.9)

    # Ring 5: SMN locus highlight on chr5 (gold)
    chr5_sector = next((s for s in circos.sectors if s.name == "chr5"), None)
    if chr5_sector:
        smn_track = chr5_sector.add_track((44, 53))
        smn_track.rect(70_924_941, 70_953_015, fc="#FFBF00", alpha=0.9)

    for link in links:
        try:
            circos.link(
                (link["chrom1"], link["pos1"], link["pos1"] + 1),
                (link["chrom2"], link["pos2"], link["pos2"] + 1),
                color=link["colour"], alpha=0.4, lw=0.5,
            )
        except Exception:
            continue

    fig = circos.plotfig(figsize=(12, 12))
    fig.savefig(out_svg)
    fig.savefig(out_png, dpi=150)  # 150 dpi sufficient for HTML; 1200 dpi produced 60 MB PNG
    plt.close(fig)
    print(f"Circos plot saved: {out_svg}, {out_png}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sv-vcf",    required=True)
    parser.add_argument("--cnv-bed",   required=True)
    parser.add_argument("--cytobands", required=True)
    parser.add_argument("--sample",    required=True)
    parser.add_argument("--out",       required=True, help="Output SVG path")
    parser.add_argument("--str-vcf",   default=None,  help="ExpansionHunter VCF for STR ring")
    args = parser.parse_args()
    out_png = args.out.replace(".svg", ".png")
    make_circos(args.sv_vcf, args.cnv_bed, args.cytobands,
                args.sample, args.out, out_png, str_vcf=args.str_vcf)


if __name__ == "__main__":
    main()
