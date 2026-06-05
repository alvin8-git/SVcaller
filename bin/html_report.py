#!/usr/bin/env python3
"""Build per-sample SVcaller HTML report using Jinja2."""
import argparse, csv, gzip, json, re
from datetime import date
from pathlib import Path
from jinja2 import Environment, FileSystemLoader


TEMPLATE_DIR = Path(__file__).parent.parent / "assets"
STR_RANGES_PATH = TEMPLATE_DIR / "str_disease_ranges.tsv"

# Canonical autosomes + sex chrs for per-chr QC table
_CANONICAL = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]
_LOW_COV_THRESHOLD = 15.0
_MAX_ARTIFACT_SIZE = 50_000_000  # chromosome-spanning calls — always artifacts
_LARGE_SV_SUPP_MIN = 2           # require multi-caller support for SVs >1 Mb


_SD_PAIR_POS_WINDOW = 50_000   # bp — start position proximity for SD pair detection
_SD_PAIR_SIZE_TOL   = 0.20     # fractional size difference tolerance

_GNOMAD_AF_THRESHOLD = 0.01   # SVs with population AF > 1% are common variation, not P/LP


def _dedup_pathogenic(rows: list) -> list:
    """Deduplicate P/LP calls: collapse same-(gene,svtype) groups, then remove
    reciprocal DUP+DEL pairs at the same locus (segmental duplication artifacts).

    SD artifact signature: same gene, one DUP + one DEL, starts within 50 kb,
    sizes within 20%. Both callers being fooled by the same SD copy-number
    ambiguity produces concordant but opposite-sign calls — not a real event.
    """
    from collections import defaultdict

    # Step 1: collapse same-(gene,svtype) groups → best-supported representative
    groups: dict = defaultdict(list)
    for r in rows:
        groups[(r["gene"], r["svtype"])].append(r)
    deduped = []
    for group in groups.values():
        group.sort(key=lambda r: (
            -(int(r["supp"]) if str(r["supp"]).isdigit() else 0),
            -(int(r["acmg_class"]) if str(r["acmg_class"]).isdigit() else 0),
        ))
        best = group[0]
        if len(group) > 1:
            best = dict(best)
            best["collapsed"] = len(group)
        deduped.append(best)

    # Step 2: remove reciprocal DUP+DEL pairs at the same locus.
    # Match by CHROMOSOME + POSITION, not gene name, because paired calls at the
    # same locus often get different primary gene annotations from AnnotSV.
    dups_all = [r for r in deduped if r["svtype"] == "DUP"]
    dels_all = [r for r in deduped if r["svtype"] == "DEL"]

    sd_artifact_ids: set = set()
    for dup in dups_all:
        for del_ in dels_all:
            if dup.get("chrom") != del_.get("chrom"):
                continue
            try:
                pos_diff  = abs(int(dup["start"]) - int(del_["start"]))
                dup_size  = dup.get("size_bp", 0)
                del_size  = del_.get("size_bp", 0)
                size_diff = abs(dup_size - del_size) / max(dup_size, del_size, 1)
            except (ValueError, TypeError):
                continue
            if pos_diff <= _SD_PAIR_POS_WINDOW and size_diff <= _SD_PAIR_SIZE_TOL:
                sd_artifact_ids.add(id(dup))
                sd_artifact_ids.add(id(del_))

    deduped = [r for r in deduped if id(r) not in sd_artifact_ids]
    deduped.sort(key=lambda r: (-int(r["acmg_class"]), r["gene"], r["svtype"]))
    return deduped


def _parse_supp_vec(info_str: str) -> dict:
    """Extract SUPP, SUPP_VEC and caller breakdown from Jasmine INFO field."""
    info_d = dict(
        (f.split("=", 1) if "=" in f else (f, ""))
        for f in info_str.split(";")
    )
    supp     = info_d.get("SUPP", "?")
    supp_vec = info_d.get("SUPP_VEC", "")
    callers  = []
    names    = ["Manta", "Delly", "GRIDSS", "Scramble", "MELT"]
    for i, bit in enumerate(supp_vec):
        if bit == "1" and i < len(names):
            callers.append(names[i])
    return {"supp": supp, "supp_vec": supp_vec, "callers": "/".join(callers) or "?"}


def _load_str_ranges() -> dict:
    """Load curated STR disease/repeat-range metadata keyed by locus_id."""
    ranges = {}
    try:
        with open(STR_RANGES_PATH) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                lid = row.get("locus_id", "").strip()
                if lid:
                    ranges[lid] = row
    except FileNotFoundError:
        pass
    return ranges


# ---------------------------------------------------------------------------
# SV parsing
# ---------------------------------------------------------------------------

def parse_sv_summary(sv_tsv_path: str) -> list:
    """Count SVs per type; flag ACMG class 4/5 (Path/LP)."""
    counts: dict = {}
    pathogenic: dict = {}
    try:
        with open(sv_tsv_path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                if row.get("Annotation_mode", "") != "full":
                    continue
                st = row.get("SV_type", row.get("SVTYPE", "UNK")).upper()
                counts[st] = counts.get(st, 0) + 1
                if str(row.get("ACMG_class", "")).strip() in ("4", "5"):
                    pathogenic[st] = pathogenic.get(st, 0) + 1
    except (FileNotFoundError, KeyError):
        pass
    return [{"svtype": k, "total": v, "high": pathogenic.get(k, 0)}
            for k, v in sorted(counts.items())]


def parse_sv_pathogenic(sv_tsv_path: str) -> list:
    """Return Class 4/5 SVs ≤10 Mb with caller support for the findings highlight."""
    rows = []
    try:
        with open(sv_tsv_path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                if row.get("Annotation_mode", "") != "full":
                    continue
                cls = str(row.get("ACMG_class", "")).strip()
                if cls not in ("4", "5"):
                    continue
                try:
                    size_bp = abs(int(row.get("SV_end", 0) or 0) - int(row.get("SV_start", 0) or 0))
                except (ValueError, TypeError):
                    size_bp = 0
                gene = row.get("Gene_name", row.get("Gene", ""))
                primary_gene = gene.split(";")[0] if gene else "—"
                supp = _parse_supp_vec(row.get("INFO", ""))
                supp_n = int(supp["supp"]) if str(supp["supp"]).isdigit() else 0
                if size_bp > _MAX_ARTIFACT_SIZE:
                    continue   # chromosome-spanning: always artifact
                if size_bp > 1_000_000 and supp_n < _LARGE_SV_SUPP_MIN:
                    continue   # large SV with single-caller support: likely artifact

                # P2 — gnomAD-SV population AF filter (hard): common variants are not P/LP
                svt = row.get("SV_type", "").upper()
                if svt in ("DUP",):
                    af_raw = row.get("B_gain_AFmax", "") or ""
                elif svt == "DEL":
                    af_raw = row.get("B_loss_AFmax", "") or ""
                elif svt == "INS":
                    af_raw = row.get("B_ins_AFmax", "") or ""
                elif svt == "INV":
                    af_raw = row.get("B_inv_AFmax", "") or ""
                else:
                    af_raw = ""
                try:
                    pop_af = float(af_raw) if af_raw and af_raw not in (".", "") else 0.0
                except ValueError:
                    pop_af = 0.0
                if pop_af > _GNOMAD_AF_THRESHOLD:
                    continue   # common variant (gnomAD-SV AF > 1%): not P/LP

                # P1 — SD boundary annotation (soft flag): both breakpoints in segdup
                sd_left  = (row.get("SegDup_left",  "") or "").strip()
                sd_right = (row.get("SegDup_right", "") or "").strip()
                both_sd  = bool(sd_left and sd_right
                                and sd_left  not in (".", "")
                                and sd_right not in (".", ""))

                rows.append({
                    "svtype":         row.get("SV_type", ""),
                    "chrom":          row.get("SV_chrom", ""),
                    "start":          row.get("SV_start", ""),
                    "end":            row.get("SV_end", ""),
                    "gene":           primary_gene,
                    "omim":           row.get("OMIM_morbid", ""),
                    "acmg_class":     cls,
                    "size":           _fmt_size(size_bp),
                    "size_bp":        size_bp,
                    "callers":        supp["callers"],
                    "supp":           supp["supp"],
                    "collapsed":      0,
                    "pop_af":         f"{pop_af:.4f}" if pop_af > 0 else "",
                    "sd_boundary":    both_sd,
                })
    except (FileNotFoundError, KeyError, ValueError):
        pass
    return _dedup_pathogenic(rows)


def parse_top_svs(sv_tsv_path: str, n: int = 20) -> list:
    """Top N annotated SVs by AnnotSV score; skip >10 Mb spanning artifacts."""
    rows = []
    try:
        with open(sv_tsv_path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                if row.get("Annotation_mode", "") != "full":
                    continue
                start = row.get("SV_start", "") or ""
                end   = row.get("SV_end",   "") or ""
                try:
                    size_bp = abs(int(end) - int(start))
                except (ValueError, TypeError):
                    size_bp = 0
                if size_bp > 10_000_000:
                    continue   # skip chromosome-spanning annotation artifacts
                acmg_score_raw = row.get("AnnotSV_ranking_score", row.get("Ranking", ""))
                try:
                    acmg_score = float(acmg_score_raw)
                except (ValueError, TypeError):
                    acmg_score = 0.0
                gene_full    = row.get("Gene_name", row.get("Gene", ""))
                primary_gene = gene_full.split(";")[0] if gene_full else "—"
                acmg_class   = str(row.get("ACMG_class", "")).strip()
                gnomad_af    = row.get("GD_AF", row.get("gnomAD_SV_AF", row.get("GNOMAD_AF", "")))
                inh_mode     = row.get("OMIM_inheritance", row.get("Gene_inheritance", ""))
                supp = _parse_supp_vec(row.get("INFO", ""))
                rows.append({
                    "chrom": row.get("SV_chrom", row.get("Chr", "")),
                    "start": start, "end": end,
                    "svtype": row.get("SV_type", ""),
                    "size":   _fmt_size(size_bp),
                    "gene":   primary_gene,
                    "gene_full": gene_full,
                    "acmg":  acmg_score_raw,
                    "acmg_class": acmg_class,
                    "omim":  row.get("OMIM_morbid", row.get("OMIM", "")),
                    "gnomad_af": gnomad_af,
                    "inheritance": inh_mode,
                    "callers": supp["callers"],
                    "supp":    supp["supp"],
                    "_score": acmg_score,
                })
    except (FileNotFoundError, KeyError):
        pass
    rows.sort(key=lambda x: x["_score"], reverse=True)
    for r in rows:
        del r["_score"]
    return rows[:n]


def _fmt_size(bp: int) -> str:
    if bp >= 1_000_000:
        return f"{bp/1_000_000:.1f} Mb"
    if bp >= 1_000:
        return f"{bp/1_000:.1f} kb"
    return f"{bp} bp"


# ---------------------------------------------------------------------------
# CNV parsing
# ---------------------------------------------------------------------------

def parse_cnv_summary(cnv_bed_path: str) -> dict:
    """Count CNV calls by confidence tier from cnv_consensus.bed."""
    summary = {"total": 0, "del": 0, "dup": 0, "both": 0, "high": 0, "medium": 0}
    try:
        with open(cnv_bed_path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 7:
                    continue
                svtype     = parts[4].upper() if len(parts) > 4 else ""
                confidence = parts[6].upper() if len(parts) > 6 else ""
                summary["total"] += 1
                if svtype == "DEL":
                    summary["del"] += 1
                elif svtype == "DUP":
                    summary["dup"] += 1
                if confidence == "BOTH":
                    summary["both"] += 1
                elif confidence in ("HIGH", "BOTH"):
                    summary["high"] += 1
                elif confidence == "MEDIUM":
                    summary["medium"] += 1
    except (FileNotFoundError, IndexError):
        pass
    return summary


# ---------------------------------------------------------------------------
# QC parsing
# ---------------------------------------------------------------------------

def parse_qc(coverage_path: str, metrics_path: str, flagstat_path: str = "",
             insert_size_path: str = "") -> dict:
    qc = {
        "mean_depth": "N/A", "dup_rate": "N/A",
        "mapped_pct": "N/A", "properly_paired_pct": "N/A",
        "total_reads": "N/A", "supplementary_pct": "N/A",
        "median_insert_size": "N/A", "mean_insert_size": "N/A",
        "insert_size_std": "N/A", "read_pairs": "N/A",
        "chr_coverage": [],    # list of {chrom, mean, flag}
        "low_cov_chrs": [],    # chrs below threshold
        "qc_pass": True,
        "qc_warnings": [],
    }
    _parse_mosdepth(coverage_path, qc)
    _parse_picard(metrics_path, qc)
    _parse_flagstat(flagstat_path, qc)
    _parse_insert_size(insert_size_path, qc)
    _compute_qc_status(qc)
    return qc


def _parse_mosdepth(summary_path: str, qc: dict) -> None:
    if not summary_path or summary_path in ("NO_FILE", "null"):
        return
    chr_cov = []
    try:
        with open(summary_path) as fh:
            for line in fh:
                if line.startswith("total\t"):
                    parts = line.strip().split("\t")
                    if len(parts) >= 4:
                        qc["mean_depth"] = f"{float(parts[3]):.1f}"
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 4 or parts[0] == "chrom":
                    continue
                chrom = parts[0]
                if chrom not in _CANONICAL:
                    continue
                try:
                    mean = float(parts[3])
                except ValueError:
                    continue
                flag = mean < _LOW_COV_THRESHOLD and chrom not in ("chrY",)
                chr_cov.append({"chrom": chrom, "mean": f"{mean:.1f}", "flag": flag})
        qc["chr_coverage"] = chr_cov
        qc["low_cov_chrs"] = [c["chrom"] for c in chr_cov if c["flag"]]
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


def _parse_insert_size(insert_size_path: str, qc: dict) -> None:
    if not insert_size_path or insert_size_path in ("NO_FILE", "null"):
        return
    try:
        with open(insert_size_path) as fh:
            lines = fh.readlines()
        for i, line in enumerate(lines):
            if "MEDIAN_INSERT_SIZE" in line:
                header = line.strip().split("\t")
                data = lines[i + 1].strip().split("\t") if i + 1 < len(lines) else []
                d = dict(zip(header, data))
                if d.get("MEDIAN_INSERT_SIZE"):
                    qc["median_insert_size"] = d["MEDIAN_INSERT_SIZE"]
                if d.get("MEAN_INSERT_SIZE"):
                    qc["mean_insert_size"] = f"{float(d['MEAN_INSERT_SIZE']):.0f}"
                if d.get("STANDARD_DEVIATION"):
                    qc["insert_size_std"] = f"{float(d['STANDARD_DEVIATION']):.0f}"
                if d.get("READ_PAIRS"):
                    n = int(d["READ_PAIRS"])
                    qc["read_pairs"] = f"{n/1_000_000:.1f}M"
                break
    except (FileNotFoundError, ValueError, IndexError):
        pass


def _parse_flagstat(flagstat_path: str, qc: dict) -> None:
    if not flagstat_path or flagstat_path in ("NO_FILE", "null"):
        return
    try:
        primary = 0
        primary_dup = 0
        supplementary = 0
        properly_paired = 0
        with open(flagstat_path) as fh:
            for line in fh:
                if " mapped (" in line and "primary" not in line:
                    m = re.search(r'\(([0-9.]+)%', line)
                    if m:
                        qc["mapped_pct"] = m.group(1)
                m = re.match(r'(\d+) \+ \d+ primary duplicates', line)
                if m:
                    primary_dup = int(m.group(1))
                m = re.match(r'(\d+) \+ \d+ primary\s*$', line)
                if m:
                    primary = int(m.group(1))
                m = re.match(r'(\d+) \+ \d+ supplementary', line)
                if m:
                    supplementary = int(m.group(1))
                m = re.match(r'(\d+) \+ \d+ properly paired', line)
                if m:
                    properly_paired = int(m.group(1))
        if primary > 0:
            qc["dup_rate"] = f"{primary_dup / primary * 100:.2f}"
            qc["total_reads"] = f"{primary/1_000_000:.1f}M"
            if properly_paired:
                qc["properly_paired_pct"] = f"{properly_paired / primary * 100:.1f}"
            total_with_supp = primary + supplementary
            if supplementary and total_with_supp > 0:
                qc["supplementary_pct"] = f"{supplementary / total_with_supp * 100:.2f}"
    except (FileNotFoundError, ValueError):
        pass


def _compute_qc_status(qc: dict) -> None:
    warnings = []
    try:
        depth = float(qc.get("mean_depth", 0))
        if depth < 20:
            warnings.append(f"Low coverage: {depth:.1f}x (threshold 20x)")
    except (ValueError, TypeError):
        pass
    try:
        dup = float(qc.get("dup_rate", 0))
        if dup > 20:
            warnings.append(f"High duplicate rate: {dup:.1f}% (threshold 20%)")
    except (ValueError, TypeError):
        pass
    if qc.get("low_cov_chrs"):
        warnings.append(f"Low coverage chromosomes: {', '.join(qc['low_cov_chrs'])}")
    try:
        mapped = float(qc.get("mapped_pct", 100))
        if mapped < 90:
            warnings.append(f"Low mapping rate: {mapped:.1f}% (threshold 90%)")
    except (ValueError, TypeError):
        pass
    qc["qc_warnings"] = warnings
    qc["qc_pass"] = len(warnings) == 0


# ---------------------------------------------------------------------------
# STR parsing
# ---------------------------------------------------------------------------

def parse_str_loci(str_vcf_path: str, str_ranges: dict = None) -> list:
    """Parse ExpansionHunter VCF; join with clinical range metadata."""
    if not str_vcf_path or str_vcf_path in ("NO_FILE", "NO_STR", "null"):
        return []
    if str_ranges is None:
        str_ranges = {}
    loci = []
    try:
        header_fields: list = []
        opener = gzip.open if str_vcf_path.endswith(".gz") else open
        with opener(str_vcf_path, "rt") as fh:
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
                chrom = d.get("CHROM", "")
                pos   = d.get("POS", "")
                info  = d.get("INFO", "")
                fmt   = d.get("FORMAT", "")
                sample_val = parts[9] if len(parts) > 9 else ""
                info_d = dict(
                    (f.split("=", 1) if "=" in f else (f, ""))
                    for f in info.split(";")
                )
                # VARID holds the locus name in EH output (e.g. VARID=FMR1)
                locus_id = info_d.get("VARID") or info_d.get("REPID") or d.get("ID", ".")
                repunit  = info_d.get("RU") or info_d.get("REPUNIT") or "."
                ref_cn   = info_d.get("REF", ".")
                fmt_d = dict(zip(fmt.split(":"), sample_val.split(":")))
                repcn     = fmt_d.get("REPCN", ".")
                repci     = fmt_d.get("REPCI", ".")
                adsp      = fmt_d.get("ADSP", ".")  # spanning reads
                adfl      = fmt_d.get("ADFL", ".")  # flanking reads
                gt        = fmt_d.get("GT", ".")
                # Join with clinical ranges
                meta = str_ranges.get(locus_id, {})
                status = _str_status(repcn, meta)
                loci.append({
                    "locus": locus_id,
                    "chrom": chrom, "pos": pos,
                    "repunit": repunit,
                    "ref_cn": ref_cn,
                    "repcn": repcn,
                    "repci": repci,
                    "spanning": adsp,
                    "flanking": adfl,
                    "gt": gt,
                    "disease":  meta.get("disease", ""),
                    "inheritance": meta.get("inheritance", ""),
                    "normal_max":  meta.get("normal_max", ""),
                    "path_min":    meta.get("path_min", ""),
                    "notes":       meta.get("notes", ""),
                    "status": status,
                })
    except FileNotFoundError:
        pass
    return loci


def _str_status(repcn: str, meta: dict) -> str:
    """Classify STR call as NORMAL / PREMUTATION / EXPANDED / UNKNOWN."""
    if not meta:
        return "UNKNOWN"
    try:
        # repcn may be allele1/allele2 — take max allele
        alleles = [int(a) for a in repcn.split("/") if a.isdigit()]
        if not alleles:
            return "UNKNOWN"
        max_cn = max(alleles)
        path_min = meta.get("path_min", "")
        premut_max = meta.get("premut_max", "")
        normal_max = meta.get("normal_max", "")
        path_min_i = int(path_min) if str(path_min).isdigit() else None
        normal_max_i = int(normal_max) if str(normal_max).isdigit() else None
        premut_max_i = int(premut_max) if str(premut_max).isdigit() else None
        if path_min_i and max_cn >= path_min_i:
            return "EXPANDED"
        if premut_max_i and normal_max_i and normal_max_i < max_cn <= premut_max_i:
            return "PREMUTATION"
        if normal_max_i and max_cn <= normal_max_i:
            return "NORMAL"
        return "INTERMEDIATE"
    except (ValueError, TypeError):
        return "UNKNOWN"


# ---------------------------------------------------------------------------
# Benchmark parsing
# ---------------------------------------------------------------------------

def parse_benchmark(json_path: str) -> list:
    """Parse Truvari summary.json — returns list with precision/recall/f1/tp/fp/fn."""
    try:
        with open(json_path) as fh:
            data = json.load(fh)
        if "precision" in data:
            return [{"svtype": "Overall",
                     "precision": data.get("precision") or 0,
                     "recall":    data.get("recall")    or 0,
                     "f1":        data.get("f1")        or 0,
                     "tp":        data.get("TP-comp")   or 0,
                     "fp":        data.get("FP")        or 0,
                     "fn":        data.get("FN")        or 0,
                     "base_cnt":  data.get("base cnt")  or 0,
                     "comp_cnt":  data.get("comp cnt")  or 0}]
        return [{"svtype": k,
                 "precision": v.get("precision") or 0,
                 "recall":    v.get("recall")    or 0,
                 "f1":        v.get("f1")        or 0,
                 "tp": 0, "fp": 0, "fn": 0, "base_cnt": 0, "comp_cnt": 0}
                for k, v in data.items() if isinstance(v, dict)]
    except (FileNotFoundError, json.JSONDecodeError):
        return []


def parse_benchmark_sizebin(json_path: str) -> list:
    """Parse per-size-bin Truvari sizebin JSON."""
    try:
        with open(json_path) as fh:
            data = json.load(fh)
        return [{"bin": k,
                 "precision": v.get("precision") or 0,
                 "recall":    v.get("recall")    or 0,
                 "f1":        v.get("f1")        or 0}
                for k, v in data.items()]
    except (FileNotFoundError, json.JSONDecodeError):
        return []


# ---------------------------------------------------------------------------
# Render
# ---------------------------------------------------------------------------

def render_report(sample_id: str, smn_html_path: str, cnv_bed_path: str,
                  sv_tsv_path: str, circos_svg_path: str, out_path: str,
                  pipeline_version: str = "1.0.0",
                  benchmark_json: str = None,
                  sizebin_json: str = None,
                  benchmark_v5q_json: str = None,
                  sizebin_v5q_json: str = None,
                  coverage_path: str = None,
                  metrics_path: str = None,
                  flagstat_path: str = None,
                  insert_size_path: str = None,
                  str_vcf_path: str = None) -> None:
    env = Environment(loader=FileSystemLoader(str(TEMPLATE_DIR)))
    template = env.get_template("report_template.html")

    str_ranges = _load_str_ranges()

    smn_html         = Path(smn_html_path).read_text()
    circos_svg_inline = Path(circos_svg_path).read_text()
    sv_summary       = parse_sv_summary(sv_tsv_path)
    sv_pathogenic    = parse_sv_pathogenic(sv_tsv_path)
    top_svs          = parse_top_svs(sv_tsv_path)
    cnv_summary      = parse_cnv_summary(cnv_bed_path)
    benchmark          = parse_benchmark(benchmark_json)           if benchmark_json      else None
    benchmark_bins     = parse_benchmark_sizebin(sizebin_json)     if sizebin_json        else None
    benchmark_v5q      = parse_benchmark(benchmark_v5q_json)       if benchmark_v5q_json  else None
    benchmark_bins_v5q = parse_benchmark_sizebin(sizebin_v5q_json) if sizebin_v5q_json    else None
    qc               = parse_qc(coverage_path or "", metrics_path or "",
                                flagstat_path or "", insert_size_path or "")
    str_loci         = parse_str_loci(str_vcf_path or "", str_ranges)

    html = template.render(
        sample_id=sample_id,
        pipeline_version=pipeline_version,
        run_date=date.today().isoformat(),
        qc=qc,
        sv_summary=sv_summary,
        sv_pathogenic=sv_pathogenic,
        top_svs=top_svs,
        cnv_summary=cnv_summary,
        smn_html=smn_html,
        circos_svg_inline=circos_svg_inline,
        benchmark=benchmark,
        benchmark_bins=benchmark_bins,
        benchmark_v5q=benchmark_v5q,
        benchmark_bins_v5q=benchmark_bins_v5q,
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
    parser.add_argument("--sizebin",          default=None)
    parser.add_argument("--benchmark-v5q",    default=None, dest="benchmark_v5q")
    parser.add_argument("--sizebin-v5q",      default=None, dest="sizebin_v5q")
    parser.add_argument("--coverage",         default=None)
    parser.add_argument("--metrics",          default=None)
    parser.add_argument("--flagstat",         default=None)
    parser.add_argument("--insert-size",      default=None, dest="insert_size")
    parser.add_argument("--str-vcf",          default=None)
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
        benchmark_v5q_json=args.benchmark_v5q,
        sizebin_v5q_json=args.sizebin_v5q,
        coverage_path=args.coverage,
        metrics_path=args.metrics,
        flagstat_path=args.flagstat,
        insert_size_path=args.insert_size,
        str_vcf_path=args.str_vcf,
    )


if __name__ == "__main__":
    main()
