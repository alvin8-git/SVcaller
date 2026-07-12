#!/usr/bin/env python3
"""Generate a Circos plot from SV VCF, CNV BED, depth BED, and AnnotSV TSV using pycirclize.

Ring layout (radius, outer → inner) — "genome fingerprint", links-first:
  96-100  chromosome ideograms
  88-95   copy-ratio HEATMAP (1 Mb bins, diverging blue→white→red = loss→neutral→gain)
  81-86   CNV consensus blocks (DUP=red, DEL=blue), min-width ticks
  74-79   STR loci barcode (EH=brown, STRling novel=orange), full-height ticks
  67-72   clinical SV ring: top SVs, coloured by ACMG class (5=red,4=orange,3=grey),
          SMN locus gold; gene/ACMG merged into one ring (was two)
  59-64   CNV-trait loci (RHD/AMY1/GSTM1/GSTT1/LPA KIV-2), labelled gene + call;
          present/normal=grey, deletion/null=blue, high-CN=red (optional; skipped
          when no trait TSVs are supplied)
  0-57    SV links centre — large canvas; interchromosomal BND/TRA emphasised

The plot is a gestalt overview, not a positional record: exact coordinates,
copy number, and confidence live in the HTML report sections and the .xlsx.
Coverage is 1 Mb-binned (re-binned from mosdepth 50 kb) for a readable fingerprint.
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

# Fixed genomic windows for the CNV-trait ring (GRCh38, chr-prefixed). These are
# reference coordinates only — the CALL for each locus (Rh status, present/null,
# copy number) is read at run time from the trait contract TSVs, never hardcoded.
TRAIT_LOCI = {
    "RHD":      ("chr1",  25272393,  25330445),
    "AMY1":     ("chr1",  103571000, 103760000),
    "GSTM1":    ("chr1",  109687814, 109693020),
    "GSTT1":    ("chr22", 24376133,  24384680),
    "LPA_KIV2": ("chr6",  160605000, 160650000),
}
# Colour buckets for the trait ring: deletion/null (loss) = blue, high copy = red,
# present/normal = grey.
TRAIT_COLOURS = {"del": "#1F77B4", "amp": "#D62728", "norm": "#7F7F7F"}


def sv_colour(svtype: str) -> str:
    return SV_COLOURS.get(svtype.upper(), "#7F7F7F")


def _read_trait_tsv(path: Optional[str]) -> Optional[dict]:
    """Read a one-record trait contract TSV ('#col1\\tcol2...' header + one data
    row) into {col: value}. Returns None for missing / sentinel inputs."""
    if not path or path in ("NO_FILE", ""):
        return None
    try:
        with open(path) as fh:
            lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
        if len(lines) < 2:
            return None
        cols = lines[0].lstrip("#").split("\t")
        vals = lines[1].split("\t")
        return dict(zip(cols, vals))
    except Exception:
        return None


def parse_cnv_traits(rh_path: Optional[str], amy1_path: Optional[str],
                     gst_path: Optional[str], lpa_path: Optional[str]) -> List[dict]:
    """Build the CNV-trait ring markers from the four trait contract TSVs.

    Returns [{chrom, start, end, label, bucket}], one per locus that has a real
    call. Values come straight from the TSVs (no hardcoded calls); an absent or
    sentinel TSV simply contributes no marker (ring skipped if none present).
    """
    def locus(name, label, bucket):
        c, s, e = TRAIT_LOCI[name]
        return {"chrom": c, "start": s, "end": e, "label": label, "bucket": bucket}

    def _as_int(v):
        try:
            return int(v)
        except (TypeError, ValueError):
            return None

    markers = []
    rh = _read_trait_tsv(rh_path)
    if rh and rh.get("Rh_status") and rh["Rh_status"] not in ("unknown", "NA"):
        st = rh["Rh_status"]; cn = rh.get("RHD_copies", "?")
        if st == "neg":
            markers.append(locus("RHD", "RHD Rh− (del)", "del"))
        else:
            markers.append(locus("RHD", "RHD Rh+ (CN{})".format(cn), "norm"))

    amy1 = _read_trait_tsv(amy1_path)
    if amy1 and amy1.get("AMY1_copies") not in (None, "NA"):
        cn = amy1["AMY1_copies"]; cni = _as_int(cn)
        bucket = "amp" if (cni is not None and cni > 2) else \
                 ("del" if (cni is not None and cni < 2) else "norm")
        markers.append(locus("AMY1", "AMY1 CN{}".format(cn), bucket))

    gst = _read_trait_tsv(gst_path)
    if gst:
        for gene in ("GSTM1", "GSTT1"):
            v = gst.get(gene)
            if v and v != "unknown":
                markers.append(locus(gene, "{} {}".format(gene, v),
                                     "del" if v == "null" else "norm"))

    lpa = _read_trait_tsv(lpa_path)
    if lpa and lpa.get("KIV2_copies") not in (None, "NA"):
        cn = lpa["KIV2_copies"]
        # LPA KIV-2 is a normal tandem VNTR (typ. ~5-40 copies); label the count,
        # keep it neutral (grey) rather than flagging it as an amplification.
        markers.append(locus("LPA_KIV2", "LPA KIV-2 CN{}".format(cn), "norm"))

    return markers


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
                    "acmg_class": acmg_class,
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


def parse_strling_tsv(path: Optional[str], min_allele2: float = 100.0,
                      min_prob: float = 0.90) -> List[dict]:
    """Parse STRling TSV for high-confidence novel expansion loci (circos STR ring)."""
    loci = []
    if not path or path in ("NO_STRLING", ""):
        return loci
    try:
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                chrom = row.get("#chrom", row.get("chrom", ""))
                if chrom not in CHROM_ORDER:
                    continue
                try:
                    a2   = float(row.get("allele2_est", 0) or 0)
                    prob = float(row.get("prob_expansion", 0) or 0)
                    pos  = int(row.get("left", 0) or 0)
                except (ValueError, TypeError):
                    continue
                if a2 < min_allele2 and prob < min_prob:
                    continue
                loci.append({"chrom": chrom, "pos": pos})
    except Exception:
        pass
    return loci


def parse_sv_vcf_links(path: str,
                       min_svlen_intra: int = 50_000,
                       max_links: int = 150) -> List[dict]:
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
            end_m = re.search(r"END=(\d+)", info)
            # Prefer SVLEN; fall back to END-POS for intrachromosomal SVs that carry
            # END but no SVLEN (some callers omit SVLEN for INV/DEL/DUP). Without this
            # fallback such SVs size to 0 and are dropped by the min_svlen_intra gate.
            if svlen_m:
                svlen = abs(int(svlen_m.group(1)))
            elif end_m:
                svlen = abs(int(end_m.group(1)) - pos1)
            else:
                svlen = 0

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

    # Interchromosomal BND/TRA are the sample's visual signature — give them the
    # lion's share of the cap (up to 110), fill the remainder with the largest
    # intrachromosomal links. Cross-genome links are always drawn richly.
    bnd   = [l for l in all_links if l["svtype"] in ("BND", "TRA")]
    intra = sorted([l for l in all_links if l["svtype"] not in ("BND", "TRA")],
                   key=lambda x: x["svlen"], reverse=True)
    max_bnd = min(len(bnd), 110)
    links = bnd[:max_bnd] + intra[:max_links - max_bnd]
    print(f"Circos: {len(links)} SV links selected ({max_bnd} interchromosomal)")
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
    # Visible tick floor: a point feature spans at least this many bp so it is not
    # sub-pixel on a 250 Mb chromosome. ~2 Mb floor, or chrom/300 for big chroms.
    return max(2_000_000, chrom_len // 300)


def _rebin_log2(depth_wins: List[dict], chrom: str, global_median: float,
                bin_bp: int = 1_000_000) -> List[Tuple[int, int, float]]:
    """Aggregate 50 kb depth windows into bin_bp (default 1 Mb) bins.

    Returns [(start, end, log2_ratio)] using the median depth per bin. Coarser bins
    turn the coverage haze into a readable copy-ratio fingerprint (uniform genome =
    flat band; CNV = a coloured block that pops out).
    """
    from collections import defaultdict
    bins = defaultdict(list)
    for w in depth_wins:
        if w["chrom"] != chrom or w["depth"] <= 0:
            continue
        bins[w["start"] // bin_bp].append(w["depth"])
    out = []
    for b in sorted(bins):
        med = statistics.median(bins[b])
        lr = math.log2(max(med, 0.01) / global_median) if global_median > 0 else 0.0
        out.append((b * bin_bp, (b + 1) * bin_bp, lr))
    return out


def make_circos(sv_vcf: str, cnv_bed: Optional[str], cytobands: str,
                sample_id: str, out_svg: str, out_png: str,
                str_vcf: Optional[str] = None,
                depth_bed: Optional[str] = None,
                annotsv_tsv: Optional[str] = None,
                strling_tsv: Optional[str] = None,
                rh_status_tsv: Optional[str] = None,
                amy1_tsv: Optional[str] = None,
                gst_null_tsv: Optional[str] = None,
                lpa_kiv2_tsv: Optional[str] = None) -> None:
    from pycirclize import Circos
    import matplotlib.pyplot as plt
    import datetime

    chrom_sizes  = load_chrom_sizes(cytobands)
    gains, losses = parse_cnv_bed(cnv_bed) if cnv_bed and cnv_bed != "NO_FILE" else ([], [])
    links        = parse_sv_vcf_links(sv_vcf)
    str_loci     = parse_str_vcf(str_vcf)       # EH catalog loci
    strling_loci = parse_strling_tsv(strling_tsv)  # high-confidence STRling novel
    depth_wins   = parse_depth_bed(depth_bed)
    gene_rows, acmg_rows = parse_annotsv_tsv(annotsv_tsv)
    trait_markers = parse_cnv_traits(rh_status_tsv, amy1_tsv, gst_null_tsv, lpa_kiv2_tsv)

    # Global median depth for fold-change normalization
    global_median = 1.0
    if depth_wins:
        depths = [w["depth"] for w in depth_wins if w["depth"] > 0]
        if depths:
            global_median = statistics.median(depths)

    print(f"Circos: depth windows={len(depth_wins)}, median={global_median:.1f}x, "
          f"cnv_gains={len(gains)}, cnv_losses={len(losses)}, "
          f"str_eh={len(str_loci)}, str_novel={len(strling_loci)}, "
          f"gene_loci={len(gene_rows)}, acmg_dots={len(acmg_rows)}, "
          f"cnv_traits={len(trait_markers)}")

    circos = Circos(chrom_sizes, space=1.5)
    circos.text(f"SVcaller\n{sample_id}", size=10, r=25)

    # --- Ring 1: chromosome ideograms (95-100) ---
    for sector in circos.sectors:
        t = sector.add_track((95, 100))
        t.axis(fc=_chrom_colour(sector.name))
        t.text(sector.name.replace("chr", ""), size=6, color="white")

    # --- Ring 2: copy-ratio HEATMAP (88-95) ---
    # 1 Mb bins (re-binned from mosdepth 50 kb), diverging colormap:
    #   blue = loss (log2 < 0), white = neutral (log2 ≈ 0), red = gain (log2 > 0).
    # vmin/vmax = ±1.0 log2 ≈ CN1..CN4, the clinically meaningful copy-ratio band.
    # A uniform genome reads as a flat near-white band; a CNV is a coloured block.
    import matplotlib.colors as _mcolors
    _cmap = plt.get_cmap("RdBu_r")
    _norm = _mcolors.TwoSlopeNorm(vmin=-1.0, vcenter=0.0, vmax=1.0)
    for sector in circos.sectors:
        t = sector.add_track((88, 95))
        t.axis(fc="#F7F7F7", ec="#CCCCCC", lw=0.3)
        clen = chrom_sizes.get(sector.name, 1)
        for (s, e, lr) in _rebin_log2(depth_wins, sector.name, global_median):
            e = min(e, clen - 1)
            if e <= s:
                continue
            lr_c = max(-1.0, min(1.0, lr))
            try:
                t.rect(s, e, fc=_mcolors.to_hex(_cmap(_norm(lr_c))), ec="none", alpha=0.95)
            except Exception:
                continue

    # --- Ring 3: CNV consensus blocks (81-86) — DUP=red, DEL=blue ---
    # Only substantial CNVs (>= 1 Mb) are drawn as discrete blocks; the heatmap ring
    # above already carries fine-scale copy-ratio, and every call (incl. small ones)
    # is in the CNV report + .xlsx. Keeps the fingerprint clean instead of a speckle.
    _CNV_MIN_BP = 1_000_000
    for sector in circos.sectors:
        t = sector.add_track((81, 86), r_pad_ratio=0.05)
        t.axis(fc="#FCFCFC", ec="#CCCCCC", lw=0.3)
        clen  = chrom_sizes.get(sector.name, 1)
        floor = _marker_width(clen)
        for rec, col in ([(r, "#D62728") for r in gains]
                         + [(r, "#1F77B4") for r in losses]):
            if rec["chrom"] != sector.name:
                continue
            if (rec["end"] - rec["start"]) < _CNV_MIN_BP:
                continue
            s = rec["start"]
            e = min(max(rec["end"], s + floor), clen - 1)
            if e <= s:
                continue
            t.rect(s, e, fc=col, ec="none", alpha=0.9)

    # --- Ring 4: STR loci barcode (74-79) — EH=brown, STRling novel=orange ---
    for sector in circos.sectors:
        t = sector.add_track((74, 79), r_pad_ratio=0.1)
        t.axis(fc="#FCFCFC", ec="#CCCCCC", lw=0.3)
        clen  = chrom_sizes.get(sector.name, 1)
        w_str = _marker_width(clen)
        for locus in str_loci:
            if locus["chrom"] == sector.name:
                end = min(locus["pos"] + w_str, clen - 1)
                t.rect(locus["pos"], end, fc="#8C564B", ec="none", alpha=0.95)
        for locus in strling_loci:
            if locus["chrom"] == sector.name:
                end = min(locus["pos"] + w_str, clen - 1)
                t.rect(locus["pos"], end, fc="#FF7F0E", ec="none", alpha=0.9)

    # --- Ring 5: clinical SV ring (67-72) — gene + ACMG merged ---
    # One ring (was two): each top SV is a min-width tick at its locus, coloured by
    # ACMG class (5=red, 4=orange, 3=grey) when classified, else by SV type. SMN gold.
    # acmg_rows is unused now — class is carried on gene_rows. Labels: top-5 by score
    # plus every class 4/5 gene (so pathogenic SVs are always named on the plot).
    label_genes = {r["gene"] for r in gene_rows[:5] if r["gene"]}
    label_genes |= {r["gene"] for r in gene_rows
                    if r.get("acmg_class") in ("4", "5") and r["gene"]}
    for sector in circos.sectors:
        t = sector.add_track((67, 72), r_pad_ratio=0.1)
        t.axis(fc="#FCFCFC", ec="#CCCCCC", lw=0.3)
        clen = chrom_sizes.get(sector.name, 1)
        w = _marker_width(clen)
        if sector.name == "chr5":
            t.rect(70_924_941, min(max(70_953_015, 70_924_941 + w), clen - 1),
                   fc="#FFBF00", ec="none", alpha=0.95)   # SMN1 locus
        for row in gene_rows:
            if row["chrom"] != sector.name:
                continue
            mid = (row["start"] + row["end"]) // 2
            end = min(mid + w, clen - 1)
            cls = row.get("acmg_class", "")
            if cls in ("3", "4", "5"):
                fc    = ACMG_COLOURS.get(cls, "#7F7F7F")
                alpha = 0.95 if cls in ("4", "5") else 0.6
            else:
                fc, alpha = sv_colour(row["svtype"]), 0.7
            if end > mid:
                t.rect(mid, end, fc=fc, ec="none", alpha=alpha)
            if row["gene"] in label_genes:
                try:
                    t.text(row["gene"], mid, size=5, color="black")
                except Exception:
                    pass

    # --- Ring 6: CNV-trait loci (59-64) — RHD/AMY1/GSTM1/GSTT1/LPA KIV-2 ---
    # Optional ring: plotted only when trait TSVs were supplied (skipped for
    # BAM-less/trait-less samples). Each locus is a labelled marker coloured by its
    # call read from the contract TSVs: deletion/null=blue, high-CN=red, else grey.
    if trait_markers:
        for sector in circos.sectors:
            t = sector.add_track((59, 64), r_pad_ratio=0.1)
            t.axis(fc="#FCFCFC", ec="#CCCCCC", lw=0.3)
            clen = chrom_sizes.get(sector.name, 1)
            w = _marker_width(clen)
            for m in trait_markers:
                if m["chrom"] != sector.name:
                    continue
                mid = (m["start"] + m["end"]) // 2
                end = min(mid + w, clen - 1)
                if end <= mid:
                    continue
                t.rect(mid, end, fc=TRAIT_COLOURS.get(m["bucket"], "#7F7F7F"),
                       ec="none", alpha=0.95)
                try:
                    t.text(m["label"], mid, size=5, color="black")
                except Exception:
                    pass

    # --- SV links (centre, r < ~57) — interchromosomal BND/TRA emphasised ---
    # The innermost ring now sits at r=67, so links own a much larger central canvas
    # than before (was r<46). Cross-genome translocations are the sample signature:
    # draw them thicker and more opaque; intrachromosomal links stay thin and faint.
    for link in links:
        interchrom = link["chrom1"] != link["chrom2"]
        lw    = 1.1 if interchrom else 0.4
        alpha = 0.55 if interchrom else 0.28
        try:
            circos.link(
                (link["chrom1"], link["pos1"], link["pos1"] + 1),
                (link["chrom2"], link["pos2"], link["pos2"] + 1),
                color=link["colour"], alpha=alpha, lw=lw,
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
        mpatches.Patch(fc="#D62728", alpha=0.95, label="Copy-ratio heatmap: gain (red)"),
        mpatches.Patch(fc="#1F77B4", alpha=0.95, label="Copy-ratio heatmap: loss (blue), 1 Mb bins"),
        mpatches.Patch(fc="#D62728", alpha=0.9,  label="CNV block: gain / DUP (red)"),
        mpatches.Patch(fc="#1F77B4", alpha=0.9,  label="CNV block: loss / DEL (blue)"),
        mpatches.Patch(fc="#8C564B", alpha=0.95, label="STR: EH catalog locus (brown)"),
        mpatches.Patch(fc="#FF7F0E", alpha=0.9,  label="STR: STRling novel (orange)"),
        mpatches.Patch(fc="#7F7F7F", alpha=0.7,  label="Clinical SV ring (ACMG-coloured)"),
        mpatches.Patch(fc="#FFBF00", alpha=0.95, label="SMN1 locus (gold)"),
        mpatches.Patch(fc="#7F7F7F", alpha=0.95, label="CNV-trait locus: present / normal (grey)"),
        mpatches.Patch(fc="#1F77B4", alpha=0.95, label="CNV-trait locus: deletion / null (blue)"),
        mlines.Line2D([], [], color="#FF7F0E", lw=1.5, alpha=0.6, label="SV links — interchromosomal"),
        _hdr("── Copy-ratio heatmap (log₂) ──"),
        mpatches.Patch(fc="#D62728", alpha=0.95, label="Gain  log₂ > 0  (red)"),
        mpatches.Patch(fc="#F7F7F7", alpha=1.0,  label="Neutral  log₂ ≈ 0  (white)"),
        mpatches.Patch(fc="#1F77B4", alpha=0.95, label="Loss  log₂ < 0  (blue); ±1 ≈ CN1..CN4"),
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
    # Drop the CNV-trait legend rows when the trait ring was not drawn, so
    # trait-less samples (e.g. BAM-less COLO829) keep a clean legend.
    if not trait_markers:
        handles = [h for h in handles
                   if not getattr(h, "get_label", lambda: "")().startswith("CNV-trait locus")]
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

    # Provenance: visible render-date footer + embedded SVG metadata, so a
    # carried-forward SVG is distinguishable from a fresh render.
    today = datetime.date.today().isoformat()
    fig.text(0.01, 0.01, f"SVcaller circos · {sample_id} · rendered {today}",
             size=6, color="#999999", ha="left", va="bottom")

    fig.savefig(out_svg, metadata={
        "Title": f"{sample_id} genome circos",
        "Creator": "SVcaller bin/circos_plot.py",
        "Date": today,
    })
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"Circos plot saved: {out_svg}, {out_png} (rendered {today})")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sv-vcf",      required=True)
    parser.add_argument("--cnv-bed",     default=None,  help="CNV consensus BED (optional)")
    parser.add_argument("--cytobands",   required=True)
    parser.add_argument("--sample",      required=True)
    parser.add_argument("--out",         required=True, help="Output SVG path")
    parser.add_argument("--str-vcf",     default=None,  help="ExpansionHunter VCF for STR ring")
    parser.add_argument("--depth-bed",   default=None,  help="mosdepth 50kb regions.bed.gz")
    parser.add_argument("--annotsv-tsv", default=None,  help="Raw AnnotSV TSV for gene/ACMG rings")
    parser.add_argument("--strling-tsv", default=None,  dest="strling_tsv",
                        help="STRling genotype TSV for novel STR ring")
    # CNV-trait ring inputs (optional): the four blood-group/copy-number trait
    # contract TSVs. Any subset may be given; the ring is skipped if none are.
    parser.add_argument("--rh-status", default=None, dest="rh_status_tsv",
                        help="rh_status.tsv (RHD/Rh) for the CNV-trait ring")
    parser.add_argument("--amy1",      default=None, dest="amy1_tsv",
                        help="amy1.tsv (AMY1 copy number) for the CNV-trait ring")
    parser.add_argument("--gst-null",  default=None, dest="gst_null_tsv",
                        help="gst_null.tsv (GSTM1/GSTT1) for the CNV-trait ring")
    parser.add_argument("--lpa-kiv2",  default=None, dest="lpa_kiv2_tsv",
                        help="lpa_kiv2.tsv (LPA KIV-2) for the CNV-trait ring")
    args = parser.parse_args()
    out_png = args.out.replace(".svg", ".png")
    make_circos(args.sv_vcf, args.cnv_bed, args.cytobands,
                args.sample, args.out, out_png,
                str_vcf=args.str_vcf,
                depth_bed=args.depth_bed,
                annotsv_tsv=args.annotsv_tsv,
                strling_tsv=args.strling_tsv,
                rh_status_tsv=args.rh_status_tsv,
                amy1_tsv=args.amy1_tsv,
                gst_null_tsv=args.gst_null_tsv,
                lpa_kiv2_tsv=args.lpa_kiv2_tsv)


if __name__ == "__main__":
    main()
