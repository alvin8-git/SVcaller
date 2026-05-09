#!/usr/bin/env python3
"""Generate a Circos plot from SV VCF and CNV BED using pycirclize.

Usage:
  circos_plot.py --sv-vcf merged.vcf.gz --cnv-bed cnv.bed \
                 --cytobands GRCh38_cytobands.txt \
                 --sample SAMPLE_ID --out circos.svg
"""
import argparse, re, gzip
from pathlib import Path
from typing import List, Tuple, Dict


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


def parse_sv_vcf_links(path: str) -> List[dict]:
    links = []
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
            end_m = re.search(r"END=(\d+)", info)
            if svtype == "BND":
                alt = parts[4]
                mate_m = re.search(r"[\[\]]([^:[\]]+):(\d+)[\[\]]", alt)
                if not mate_m:
                    continue
                chrom2, pos2 = mate_m.group(1), int(mate_m.group(2))
            else:
                chrom2 = chrom1
                pos2 = int(end_m.group(1)) if end_m else pos1 + 1000
            if chrom1 not in CHROM_ORDER or chrom2 not in CHROM_ORDER:
                continue
            links.append({
                "chrom1": chrom1, "pos1": pos1,
                "chrom2": chrom2, "pos2": pos2,
                "svtype": svtype, "colour": sv_colour(svtype),
            })
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
                sample_id: str, out_svg: str, out_png: str) -> None:
    from pycirclize import Circos
    import matplotlib.pyplot as plt

    chrom_sizes = load_chrom_sizes(cytobands)
    gains, losses = parse_cnv_bed(cnv_bed)
    links = parse_sv_vcf_links(sv_vcf)

    circos = Circos(chrom_sizes, space=1.5)
    circos.text(f"SVcaller\n{sample_id}", size=10, r=15)

    for sector in circos.sectors:
        track = sector.add_track((95, 100))
        track.axis(fc=_chrom_colour(sector.name))
        track.text(sector.name.replace("chr", ""), size=6, color="white")

    for sector in circos.sectors:
        track = sector.add_track((80, 93), r_pad_ratio=0.1)
        track.axis()
        sector_gains = [(g["start"], g["end"], g["cn"] - 2)
                        for g in gains if g["chrom"] == sector.name]
        for start, end, height in sector_gains:
            if height > 0:
                track.rect(start, end, fc="#D62728", alpha=0.7)

    for sector in circos.sectors:
        track = sector.add_track((67, 80), r_pad_ratio=0.1)
        track.axis()
        sector_losses = [(l["start"], l["end"]) for l in losses if l["chrom"] == sector.name]
        for start, end in sector_losses:
            track.rect(start, end, fc="#1F77B4", alpha=0.7)

    chr5_sector = next((s for s in circos.sectors if s.name == "chr5"), None)
    if chr5_sector:
        smn_track = chr5_sector.add_track((60, 66))
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
    fig.savefig(out_svg, dpi=150)
    fig.savefig(out_png, dpi=300)
    plt.close(fig)
    print(f"Circos plot saved: {out_svg}, {out_png}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sv-vcf",    required=True)
    parser.add_argument("--cnv-bed",   required=True)
    parser.add_argument("--cytobands", required=True)
    parser.add_argument("--sample",    required=True)
    parser.add_argument("--out",       required=True, help="Output SVG path")
    args = parser.parse_args()
    out_png = args.out.replace(".svg", ".png")
    make_circos(args.sv_vcf, args.cnv_bed, args.cytobands,
                args.sample, args.out, out_png)


if __name__ == "__main__":
    main()
