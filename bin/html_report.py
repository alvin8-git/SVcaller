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


def parse_qc_stub() -> dict:
    return {"mean_depth": "N/A", "dup_rate": "N/A", "mapped_pct": "N/A"}


def parse_benchmark(json_path: str) -> list:
    try:
        with open(json_path) as fh:
            data = json.load(fh)
        rows = []
        for svtype, metrics in data.items():
            rows.append({
                "svtype": svtype,
                "precision": metrics.get("precision", 0),
                "recall":    metrics.get("recall",    0),
                "f1":        metrics.get("f1",        0),
            })
        return rows
    except (FileNotFoundError, json.JSONDecodeError):
        return []


def render_report(sample_id: str, smn_html_path: str, cnv_bed_path: str,
                  sv_tsv_path: str, circos_svg_path: str, out_path: str,
                  pipeline_version: str = "1.0.0",
                  benchmark_json: str = None) -> None:
    env = Environment(loader=FileSystemLoader(str(TEMPLATE_DIR)))
    template = env.get_template("report_template.html")

    smn_html = Path(smn_html_path).read_text()
    circos_svg_inline = Path(circos_svg_path).read_text()
    sv_summary = parse_sv_summary(sv_tsv_path)
    top_svs    = parse_top_svs(sv_tsv_path)
    benchmark  = parse_benchmark(benchmark_json) if benchmark_json else None
    qc         = parse_qc_stub()

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
    )


if __name__ == "__main__":
    main()
