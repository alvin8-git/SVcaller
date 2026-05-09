#!/usr/bin/env python3
"""SMN1/SMN2 copy number parsing, classification, and HTML report section generator."""
import json, argparse
from pathlib import Path
from typing import Optional


def classify_sma(smn1_cn: int, smn2_cn: int) -> dict:
    """Classify SMA status from SMN1 copy number."""
    if smn1_cn == 0:
        status = "Affected"
        badge_class = "badge-danger"
    elif smn1_cn == 1:
        status = "Carrier"
        badge_class = "badge-warning"
    else:
        status = "Normal"
        badge_class = "badge-success"
    return {
        "status": status,
        "badge_class": badge_class,
        "smn1_cn": smn1_cn,
        "smn2_cn": smn2_cn,
        "interpretation": _interpretation(smn1_cn, smn2_cn),
    }


def detect_two_plus_zero(smn1_cn: int, smn1_allele1: int, smn1_allele2: int) -> bool:
    """Detect 2+0 haplotype: appears as CN=2 but one allele carries 0 copies."""
    return smn1_cn == 2 and (smn1_allele1 == 0 or smn1_allele2 == 0)


def _interpretation(smn1_cn: int, smn2_cn: int) -> str:
    if smn1_cn == 0:
        sma_type = {1: "Type I (severe)", 2: "Type II/III", 3: "Type III", 4: "Type IV (mild)"}.get(smn2_cn, "severity uncertain")
        return (f"Homozygous SMN1 deletion. Consistent with SMA. "
                f"SMN2 copy number = {smn2_cn}: predicts {sma_type}.")
    if smn1_cn == 1:
        return f"SMA carrier (1 functional SMN1 copy). SMN2 CN = {smn2_cn}."
    return f"Normal SMN1 copy number ({smn1_cn}). SMN2 CN = {smn2_cn}."


def parse_smn_tsv(tsv_path: str) -> dict:
    """Parse SMNCopyNumberCaller TSV output and return structured dict."""
    result = {
        "smn1_cn": None, "smn2_cn": None,
        "smn1_allele1": None, "smn1_allele2": None,
        "confidence": "UNKNOWN",
    }
    with open(tsv_path) as fh:
        lines = [l.strip() for l in fh if l.strip() and not l.startswith("#")]
    if not lines:
        return result
    header = lines[0].split("\t")
    values = lines[1].split("\t") if len(lines) > 1 else []
    d = dict(zip(header, values))
    try:
        result["smn1_cn"]      = int(d.get("SMN1_CN", d.get("smn1", 2)))
        result["smn2_cn"]      = int(d.get("SMN2_CN", d.get("smn2", 2)))
        result["smn1_allele1"] = int(d.get("SMN1_allele1", result["smn1_cn"]))
        result["smn1_allele2"] = int(d.get("SMN1_allele2", 0))
        result["confidence"]   = d.get("Confidence", "HIGH")
    except (ValueError, KeyError):
        pass
    return result


def render_html_section(sample_id: str, tsv_path: str) -> str:
    """Return an HTML string for the SMN section of the per-sample report."""
    parsed   = parse_smn_tsv(tsv_path)
    smn1     = parsed["smn1_cn"] or 2
    smn2     = parsed["smn2_cn"] or 2
    a1       = parsed["smn1_allele1"] or smn1
    a2       = parsed["smn1_allele2"] or 0
    two_zero = detect_two_plus_zero(smn1, a1, a2)
    cls_info = classify_sma(smn1, smn2)
    warn     = '<div class="alert alert-warning">&#x26A0; 2+0 haplotype detected: sample appears CN=2 but may be an SMA carrier.</div>' if two_zero else ""
    return f"""
<div class="card mb-3">
  <div class="card-header"><h5>SMN1/SMN2 Copy Number — {sample_id}</h5></div>
  <div class="card-body">
    {warn}
    <table class="table table-sm">
      <tr><th>Gene</th><th>Copy Number</th></tr>
      <tr><td>SMN1</td><td><strong>{smn1}</strong></td></tr>
      <tr><td>SMN2</td><td><strong>{smn2}</strong></td></tr>
    </table>
    <p><span class="badge {cls_info['badge_class']}">{cls_info['status']}</span>
       &nbsp; {cls_info['interpretation']}</p>
    <small class="text-muted">Confidence: {parsed['confidence']}</small>
  </div>
</div>"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tsv",    required=True, help="SMNCopyNumberCaller TSV output")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--out",    required=True, help="Output HTML snippet path")
    args = parser.parse_args()
    html = render_html_section(args.sample, args.tsv)
    Path(args.out).write_text(html)
    print(f"SMN HTML section written to {args.out}")


if __name__ == "__main__":
    main()
