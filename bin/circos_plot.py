#!/usr/bin/env python3
"""Generate a Circos plot from SV VCF, CNV BED, depth BED, and AnnotSV TSV using pycirclize.

Ring layout (radius):
  95-100  chromosome ideograms
  64-95   coverage depth dot plot (log2 ratio vs median; -2=CN0, 0=CN2, +1=CN4)
  54-64   STR expansion loci
  44-54   AnnotSV gene loci (top 30 by ranking score) + SMN marker
  34-44   ACMG class dots (class 3/4/5, cap 50)
  0-34    SV links centre (34% of diameter — ≥ 1/3)
"""
import argparse, re, gzip, csv, statistics, math
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

ACMG_COLOURS = {
    "5": "#D62728",
    "4": "#FF7F0E",
    "3": "#7F7F7F",
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


def parse_depth_bed(path: Optional[str]) -> List[dict]:
    """Parse mosdepth --by 50000 regions BED (chrom, start, end, mean_depth)."""
    windows = []
    if not path or path in ("NO_FILE", ""):
        return windows
    opener = gzip.open if str(path).endswith(".gz") else open
    try:
        with opener(path, "rt") as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 4:
                    continue
                chrom = parts[0]
                if chrom not in CHROM_ORDER:
                    continue
                try:
                    windows.append({
                        "chrom": chrom,
                        "start": int(parts[1]),
                        "end":   int(parts[2]),
                        "depth": float(parts[3]),
                    })
                except (ValueError, IndexError):
                    continue
    except Exception:
        pass
    return windows


def parse_annotsv_tsv(path: Optional[str],
                      max_genes: int = 30,
                      max_acmg: int = 50) -> Tuple[List[dict], List[dict]]:
    """Parse raw AnnotSV TSV for gene loci (Ring B) and ACMG class dots (Ring C).

    Uses Annotation_mode=='full' rows (one per SV, not per transcript).
    """
    gene_rows, acmg_rows = [], []
    if not path or path in ("NO_FILE", ""):
        return gene_rows, acmg_rows
    try:
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                if row.get("Annotation_mode", "") != "full":
                    continue
                chrom = row.get("SV_chrom", "")
                if not chrom.startswith("chr"):
                    chrom = "chr" + chrom
                if chrom not in CHROM_ORDER:
                    continue
                try:
                    start = int(row.get("SV_start", 0))
                    end   = int(row.get("SV_end", start))
                except ValueError:
                    continue
                gene    = (row.get("Gene_name", "") or "").split(";")[0].strip()
                svtype  = row.get("SV_type", "")
                try:
                    score = float(row.get("AnnotSV_ranking_score", "0") or "0")
                except ValueError:
                    score = 0.0
                acmg_class = str(row.get("ACMG_class", "")).strip()

                gene_rows.append({
                    "chrom": chrom, "start": start, "end": end,
                    "gene": gene, "svtype": svtype, "score": score,
                })
                if acmg_class in ("3", "4", "5"):
                    acmg_rows.append({
                        "chrom": chrom, "start": start, "end": end,
                        "gene": gene, "acmg_class": acmg_class, "score": score,
                    })
    except Exception:
        pass

    gene_rows.sort(key=lambda x: x["score"], reverse=True)
    acmg_rows.sort(key=lambda x: (x["acmg_class"], x["score"]), reverse=True)
    return gene_rows[:max_genes], acmg_rows[:max_acmg]


def parse_sv_vcf_links(path: str,
                       min_svlen_intra: int = 50_000,
                       max_links: int = 100) -> List[dict]:
    """Return filtered SV links for clinical Circos display."""
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
            if svtype == "INS":
                continue
            supp_m = re.search(r"SUPP_VEC=([01]+)", info)
            if supp_m and supp_m.group(1).count("1") < 2:
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
                    continue
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

    bnd   = [l for l in all_links if l["svtype"] in ("BND", "TRA")]
    intra = sorted([l for l in all_links if l["svtype"] not in ("BND", "TRA")],
                   key=lambda x: x["svlen"], reverse=True)
    max_bnd = min(len(bnd), max_links // 2)
    links = bnd[:max_bnd] + intra[:max_links - max_bnd]
    print(f"Circos: {len(links)} SV links selected for display")
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


def _marker_width(chrom_len: int) -> int:
    return max(500_000, chrom_len // 300)


def make_circos(sv_vcf: str, cnv_bed: str, cytobands: str,
                sample_id: str, out_svg: str, out_png: str,
                str_vcf: Optional[str] = None,
                depth_bed: Optional[str] = None,
                annotsv_tsv: Optional[str] = None) -> None:
    from pycirclize import Circos
    import matplotlib.pyplot as plt

    chrom_sizes  = load_chrom_sizes(cytobands)
    gains, losses = parse_cnv_bed(cnv_bed)
    links        = parse_sv_vcf_links(sv_vcf)
    str_loci     = parse_str_vcf(str_vcf)
    depth_wins   = parse_depth_bed(depth_bed)
    gene_rows, acmg_rows = parse_annotsv_tsv(annotsv_tsv)

    # Global median depth for fold-change normalization
    global_median = 1.0
    if depth_wins:
        depths = [w["depth"] for w in depth_wins if w["depth"] > 0]
        if depths:
            global_median = statistics.median(depths)

    print(f"Circos: depth windows={len(depth_wins)}, median={global_median:.1f}x, "
          f"gene_loci={len(gene_rows)}, acmg_dots={len(acmg_rows)}")

    circos = Circos(chrom_sizes, space=1.5)
    circos.text(f"SVcaller\n{sample_id}", size=10, r=17)

    # --- Ring 1: chromosome ideograms (95-100) ---
    for sector in circos.sectors:
        t = sector.add_track((95, 100))
        t.axis(fc=_chrom_colour(sector.name))
        t.text(sector.name.replace("chr", ""), size=6, color="white")

    # --- Ring 2: Coverage depth dot plot (64-95) ---
    # Y-axis = log2(window_depth / genome_median), normalised to [0,1]:
    #   0.0  → log2 = -2  (CN=0, homozygous deletion)
    #   0.25 → log2 = -1  (CN=1, hemizygous loss)
    #   0.50 → log2 =  0  (CN=2, normal diploid)  ← dashed reference line
    #   0.65 → log2 ≈+0.58 (CN=3, +1 copy gain)
    #   0.75 → log2 = +1  (CN=4, +2 copy gain)
    #   1.0  → log2 = +2  (high-level amplification)
    for sector in circos.sectors:
        t = sector.add_track((64, 95), r_pad_ratio=0.05)
        t.axis()
        clen = chrom_sizes.get(sector.name, 1)
        # Dashed reference line at CN=2 (log2=0 → y=0.5)
        try:
            t.line([0, clen], [0.5, 0.5], color="#888888", lw=0.8, alpha=0.5)
        except Exception:
            pass
        # Group windows by colour to avoid per-point colour lists
        normal_x, normal_y = [], []
        gain_x,   gain_y   = [], []
        loss_x,   loss_y   = [], []
        for w in depth_wins:
            if w["chrom"] != sector.name:
                continue
            mid = (w["start"] + w["end"]) // 2
            lr  = math.log2(max(w["depth"], 0.01) / global_median) if global_median > 0 else 0.0
            lr  = max(-2.0, min(2.0, lr))
            y_n = (lr + 2.0) / 4.0     # [-2, +2] → [0, 1]
            if lr > 0.3:                # CN ≥ 2.5 (gain)
                gain_x.append(mid);   gain_y.append(y_n)
            elif lr < -0.5:             # CN ≤ 1.4 (loss)
                loss_x.append(mid);   loss_y.append(y_n)
            else:
                normal_x.append(mid); normal_y.append(y_n)
        # Subsample normal dots (every 10th) to keep SVG size manageable.
        # Gain/loss dots always shown at full resolution — they are the signal.
        if normal_x:
            t.scatter(normal_x[::10], normal_y[::10], color="#AAAAAA", s=0.4, alpha=0.5)
        if gain_x:
            t.scatter(gain_x,   gain_y,   color="#D62728", s=0.8, alpha=0.85)
        if loss_x:
            t.scatter(loss_x,   loss_y,   color="#1F77B4", s=0.8, alpha=0.85)

    # --- Ring 3: STR loci (54-64) ---
    for sector in circos.sectors:
        t = sector.add_track((54, 64), r_pad_ratio=0.1)
        t.axis()
        clen = chrom_sizes.get(sector.name, 1)
        for locus in str_loci:
            if locus["chrom"] == sector.name:
                w = max(500_000, clen // 500)
                t.rect(locus["pos"], locus["pos"] + w, fc="#8C564B", alpha=0.9)

    # --- Ring 4: AnnotSV gene loci (44-54) + SMN locus (chr5, gold) ---
    top5_genes = {r["gene"] for r in gene_rows[:5] if r["gene"]}
    for sector in circos.sectors:
        t = sector.add_track((44, 54), r_pad_ratio=0.1)
        t.axis()
        clen = chrom_sizes.get(sector.name, 1)
        if sector.name == "chr5":
            t.rect(70_924_941, 70_953_015, fc="#FFBF00", alpha=0.9)
        for row in gene_rows:
            if row["chrom"] != sector.name:
                continue
            mid = (row["start"] + row["end"]) // 2
            w = _marker_width(clen)
            t.rect(mid, mid + w, fc=sv_colour(row["svtype"]), alpha=0.8)
            if row["gene"] in top5_genes:
                try:
                    t.text(row["gene"], mid, size=5, color="black")
                except Exception:
                    pass

    # --- Ring 5: ACMG class dots (34-44) ---
    top3_acmg_genes = set(list({r["gene"] for r in acmg_rows
                                 if r["acmg_class"] in ("4", "5") and r["gene"]})[:3])
    for sector in circos.sectors:
        t = sector.add_track((34, 44), r_pad_ratio=0.1)
        t.axis()
        clen = chrom_sizes.get(sector.name, 1)
        for row in acmg_rows:
            if row["chrom"] != sector.name:
                continue
            mid = (row["start"] + row["end"]) // 2
            w = _marker_width(clen)
            fc    = ACMG_COLOURS.get(row["acmg_class"], "#7F7F7F")
            alpha = 0.9 if row["acmg_class"] in ("4", "5") else 0.55
            t.rect(mid, mid + w, fc=fc, alpha=alpha)
            if row["gene"] in top3_acmg_genes:
                try:
                    t.text(row["gene"], mid, size=5, color="black")
                except Exception:
                    pass

    # --- SV links (center, r<=27) ---
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

    # --- Legend ---
    import matplotlib.patches as mpatches
    import matplotlib.lines as mlines

    def _hdr(label):
        return mpatches.Patch(fc="none", ec="none", label=label)

    handles = [
        _hdr("── Rings (outer → inner) ──"),
        mpatches.Patch(fc="#888888", alpha=0.8,  label="Chromosomes"),
        mpatches.Patch(fc="#AAAAAA", alpha=0.7,  label="Coverage depth (50 kb dots, log₂ ratio)"),
        mpatches.Patch(fc="#8C564B", alpha=0.9,  label="STR expansion loci"),
        mpatches.Patch(fc="#888888", alpha=0.5,  label="Gene loci (top 30 by AnnotSV score)"),
        mpatches.Patch(fc="#7F7F7F", alpha=0.6,  label="ACMG class dots"),
        mlines.Line2D([], [], color="#FF7F0E", lw=1.5, alpha=0.6, label="SV links (centre)"),
        _hdr("── Coverage depth (log₂ ratio) ──"),
        mpatches.Patch(fc="#D62728", alpha=0.85, label="Gain  log₂ > +0.3  (CN ≥ 3)"),
        mpatches.Patch(fc="#AAAAAA", alpha=0.6,  label="Normal  log₂ ≈ 0  (CN = 2)"),
        mpatches.Patch(fc="#1F77B4", alpha=0.85, label="Loss  log₂ < −0.5  (CN ≤ 1)"),
        _hdr("── SV Types ──"),
        mpatches.Patch(fc="#1F77B4", label="DEL"),
        mpatches.Patch(fc="#D62728", label="DUP"),
        mpatches.Patch(fc="#9467BD", label="INV"),
        mpatches.Patch(fc="#FF7F0E", label="BND / TRA"),
        mpatches.Patch(fc="#2CA02C", label="INS"),
        _hdr("── ACMG Classification ──"),
        mpatches.Patch(fc="#D62728", label="Class 5 — Pathogenic"),
        mpatches.Patch(fc="#FF7F0E", label="Class 4 — Likely Pathogenic"),
        mpatches.Patch(fc="#7F7F7F", label="Class 3 — VUS"),
    ]
    fig.legend(
        handles=handles,
        loc="lower right",
        fontsize=6,
        framealpha=0.92,
        frameon=True,
        edgecolor="#AAAAAA",
        borderpad=0.8,
        handlelength=1.0,
        handleheight=0.9,
        labelspacing=0.35,
        bbox_to_anchor=(0.99, 0.01),
        bbox_transform=fig.transFigure,
    )

    fig.savefig(out_svg)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"Circos plot saved: {out_svg}, {out_png}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sv-vcf",      required=True)
    parser.add_argument("--cnv-bed",     required=True)
    parser.add_argument("--cytobands",   required=True)
    parser.add_argument("--sample",      required=True)
    parser.add_argument("--out",         required=True, help="Output SVG path")
    parser.add_argument("--str-vcf",     default=None,  help="ExpansionHunter VCF for STR ring")
    parser.add_argument("--depth-bed",   default=None,  help="mosdepth 50kb regions.bed.gz")
    parser.add_argument("--annotsv-tsv", default=None,  help="Raw AnnotSV TSV for gene/ACMG rings")
    args = parser.parse_args()
    out_png = args.out.replace(".svg", ".png")
    make_circos(args.sv_vcf, args.cnv_bed, args.cytobands,
                args.sample, args.out, out_png,
                str_vcf=args.str_vcf,
                depth_bed=args.depth_bed,
                annotsv_tsv=args.annotsv_tsv)


if __name__ == "__main__":
    main()
