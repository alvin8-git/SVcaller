#!/usr/bin/env python3
"""Build per-sample SVcaller HTML report using Jinja2."""
import argparse, csv, json
from datetime import date
from pathlib import Path
from jinja2 import Environment, FileSystemLoader


TEMPLATE_DIR = Path(__file__).parent.parent / "assets"


def parse_sv_summary(sv_tsv_path: str) -> list:
    counts: dict = {}
    high: dict = {}
    try:
        with open(sv_tsv_path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                st = row.get("SV_type", row.get("SVTYPE", "UNK")).upper()
                supp = row.get("SUPPORT", row.get("caller_support", ""))
                counts[st] = counts.get(st, 0) + 1
                if "BOTH" in supp or "HIGH" in supp:
                    high[st] = high.get(st, 0) + 1
    except (FileNotFoundError, KeyError):
        pass
    return [{"svtype": k, "total": v, "high": high.get(k, 0)}
            for k, v in sorted(counts.items())]


def parse_top_svs(sv_tsv_path: str, n: int = 20) -> list:
    rows = []
    try:
        with open(sv_tsv_path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                acmg = row.get("AnnotSV_ranking_score", row.get("Ranking", ""))
                try:
                    acmg_score = float(acmg)
                except (ValueError, TypeError):
                    acmg_score = 0
                if acmg_score >= 0.9:
                    chrom = row.get("SV_chrom", row.get("Chr", ""))
                    start = row.get("SV_start", row.get("Start", ""))
                    end   = row.get("SV_end",   row.get("End",   ""))
                    size  = abs(int(end) - int(start)) if start and end else 0
                    rows.append({
                        "chrom": chrom, "start": start, "end": end,
                        "svtype": row.get("SV_type", ""),
                        "size": _fmt_size(size),
                        "gene": row.get("Gene_name", row.get("Gene", "")),
                        "acmg": acmg,
                        "omim": row.get("OMIM_morbid", row.get("OMIM", "")),
                    })
    except (FileNotFoundError, KeyError):
        pass
    return rows[:n]


def _fmt_size(bp: int) -> str:
    if bp >= 1_000_000:
        return f"{bp/1_000_000:.1f} Mb"
    if bp >= 1_000:
        return f"{bp/1_000:.1f} kb"
    return f"{bp} bp"


def parse_qc(coverage_path: str, metrics_path: str) -> dict:
    qc = {"mean_depth": "N/A", "dup_rate": "N/A", "mapped_pct": "N/A"}
    _parse_mosdepth(coverage_path, qc)
    _parse_picard(metrics_path, qc)
    return qc


def _parse_mosdepth(summary_path: str, qc: dict) -> None:
    if not summary_path or summary_path in ("NO_FILE", "null"):
        return
    try:
        with open(summary_path) as fh:
            for line in fh:
                if line.startswith("total\t"):
                    parts = line.strip().split("\t")
                    if len(parts) >= 4:
                        qc["mean_depth"] = f"{float(parts[3]):.1f}"
                    break
    except (FileNotFoundError, ValueError):
        pass


def _parse_picard(metrics_path: str, qc: dict) -> None:
    if not metrics_path or metrics_path in ("NO_FILE", "null"):
        return
    try:
        with open(metrics_path) as fh:
            lines = fh.readlines()
        for i, line in enumerate(lines):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if "PERCENT_DUPLICATION" in stripped:
                header = stripped.split("\t")
                data = lines[i + 1].strip().split("\t") if i + 1 < len(lines) else []
                d = dict(zip(header, data))
                pct = d.get("PERCENT_DUPLICATION", "")
                if pct:
                    qc["dup_rate"] = f"{float(pct) * 100:.2f}"
                break
    except (FileNotFoundError, ValueError, IndexError):
        pass


def parse_str_loci(str_vcf_path: str) -> list:
    if not str_vcf_path or str_vcf_path in ("NO_FILE", "NO_STR", "null"):
        return []
    loci = []
    try:
        header_fields: list = []
        with open(str_vcf_path) as fh:
            for line in fh:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    header_fields = line.strip().lstrip("#").split("\t")
                    continue
                parts = line.strip().split("\t")
                if not parts or len(parts) < 8:
                    continue
                d = dict(zip(header_fields, parts))
                locus_id = d.get("ID", ".")
                chrom = d.get("CHROM", "")
                pos   = d.get("POS", "")
                info  = d.get("INFO", "")
                fmt   = d.get("FORMAT", "")
                sample_val = parts[9] if len(parts) > 9 else ""
                fmt_d = dict(zip(fmt.split(":"), sample_val.split(":")))
                repcn   = fmt_d.get("REPCN", ".")
                repunit = next(
                    (f.split("=", 1)[1] for f in info.split(";") if f.startswith("REPUNIT=")),
                    "."
                )
                loci.append({
                    "locus": locus_id, "chrom": chrom, "pos": pos,
                    "repunit": repunit, "repcn": repcn,
                })
    except FileNotFoundError:
        pass
    return loci


def parse_benchmark(json_path: str) -> list:
    """Parse Truvari summary.json — flat dict with top-level precision/recall/f1."""
    try:
        with open(json_path) as fh:
            data = json.load(fh)
        if "precision" in data:
            return [{"svtype": "Overall",
                     "precision": data.get("precision", 0),
                     "recall":    data.get("recall",    0),
                     "f1":        data.get("f1",        0)}]
        # Fallback: per-svtype format
        return [{"svtype": k,
                 "precision": v.get("precision", 0),
                 "recall": v.get("recall", 0),
                 "f1": v.get("f1", 0)}
                for k, v in data.items() if isinstance(v, dict)]
    except (FileNotFoundError, json.JSONDecodeError):
        return []


def parse_benchmark_sizebin(json_path: str) -> list:
    """Parse per-size-bin Truvari sizebin JSON produced by TRUVARI_BENCH."""
    try:
        with open(json_path) as fh:
            data = json.load(fh)
        return [{"bin": k,
                 "precision": v.get("precision", 0),
                 "recall":    v.get("recall",    0),
                 "f1":        v.get("f1",        0)}
                for k, v in data.items()]
    except (FileNotFoundError, json.JSONDecodeError):
        return []


def render_report(sample_id: str, smn_html_path: str, cnv_bed_path: str,
                  sv_tsv_path: str, circos_svg_path: str, out_path: str,
                  pipeline_version: str = "1.0.0",
                  benchmark_json: str = None,
                  sizebin_json: str = None,
                  coverage_path: str = None,
                  metrics_path: str = None,
                  str_vcf_path: str = None) -> None:
    env = Environment(loader=FileSystemLoader(str(TEMPLATE_DIR)))
    template = env.get_template("report_template.html")

    smn_html = Path(smn_html_path).read_text()
    circos_svg_inline = Path(circos_svg_path).read_text()
    sv_summary = parse_sv_summary(sv_tsv_path)
    top_svs    = parse_top_svs(sv_tsv_path)
    benchmark      = parse_benchmark(benchmark_json) if benchmark_json else None
    benchmark_bins = parse_benchmark_sizebin(sizebin_json) if sizebin_json else None
    qc             = parse_qc(coverage_path or "", metrics_path or "")
    str_loci       = parse_str_loci(str_vcf_path or "")

    html = template.render(
        sample_id=sample_id,
        pipeline_version=pipeline_version,
        run_date=date.today().isoformat(),
        qc=qc,
        sv_summary=sv_summary,
        top_svs=top_svs,
        smn_html=smn_html,
        circos_svg_inline=circos_svg_inline,
        benchmark=benchmark,
        benchmark_bins=benchmark_bins,
        str_loci=str_loci,
    )
    Path(out_path).write_text(html)
    print(f"HTML report written to {out_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample",           required=True)
    parser.add_argument("--smn-html",         required=True)
    parser.add_argument("--cnv-bed",          required=True)
    parser.add_argument("--sv-tsv",           required=True)
    parser.add_argument("--circos-svg",       required=True)
    parser.add_argument("--out",              required=True)
    parser.add_argument("--pipeline-version", default="1.0.0")
    parser.add_argument("--benchmark",        default=None)
    parser.add_argument("--sizebin",          default=None, help="per-size-bin Truvari JSON")
    parser.add_argument("--coverage",         default=None, help="mosdepth summary.txt")
    parser.add_argument("--metrics",          default=None, help="Picard MarkDup metrics")
    parser.add_argument("--str-vcf",          default=None, help="ExpansionHunter VCF")
    args = parser.parse_args()
    render_report(
        sample_id=args.sample,
        smn_html_path=args.smn_html,
        cnv_bed_path=args.cnv_bed,
        sv_tsv_path=args.sv_tsv,
        circos_svg_path=args.circos_svg,
        out_path=args.out,
        pipeline_version=args.pipeline_version,
        benchmark_json=args.benchmark,
        sizebin_json=args.sizebin,
        coverage_path=args.coverage,
        metrics_path=args.metrics,
        str_vcf_path=args.str_vcf,
    )


if __name__ == "__main__":
    main()
