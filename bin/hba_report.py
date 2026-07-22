#!/usr/bin/env python3
"""Alpha-globin HTML card for the per-sample report. Mirrors bin/smn_report.py.

DELIBERATELY THIN AND FACTUAL. It renders the measurement and the scope of the
measurement — nothing else. It must never:

  * use the word "thalassemia" in a result line,
  * say "carrier", "clear", "negative", "normal" or "low risk",
  * classify HbH / Bart's / trait, or state any couple-level risk.

OmniGen owns every clinical statement; this module owns the numbers. The
`not_screened` declaration is rendered INLINE, not as a footnote, because a
footnote is what let an empty SMN result render as "0 Carrier findings ...
Clear" in production.
"""
import argparse
import html
from pathlib import Path

# The tiers OmniGen renders. Human-readable text for each id it may carry.
TIER_TEXT = {
    "beta_globin":
        "Beta-globin (HBB) was not examined by this module at all. Most "
        "beta-globin alleles are point mutations and are read from a germline "
        "VCF elsewhere.",
    "alpha_nondeletional_outside_panel":
        "Only the named positions in the pinned site panel were examined. Any "
        "alpha-globin sequence change outside that list was not looked for.",
    "alpha_deletional":
        "The deletion channel did not run for this sample, so alpha-globin "
        "copy number was not measured.",
    "alpha_targeted_sites":
        "The targeted-site channel did not run for this sample, so no panel "
        "position was examined.",
}


class AlphaReportInputError(RuntimeError):
    """The contract TSV exists but carries no data row => the module ran and
    failed. Never render that as an absence of findings."""


def parse_contract(path):
    """Read <S>.alpha_globin.tsv. Absent -> None. Present-but-empty -> raise."""
    name = Path(str(path)).name
    if name.startswith("NO_") or not Path(str(path)).exists():
        return None
    lines = [l.rstrip("\n") for l in Path(str(path)).read_text().splitlines()
             if l.strip()]
    if len(lines) < 2:
        raise AlphaReportInputError(
            f"alpha-globin contract '{path}' has no data row -- the module ran "
            f"and produced no result. Refusing to render that as 'nothing "
            f"found', which is how an empty SMN artifact once became a clean "
            f"bill of health.")
    return dict(zip(lines[0].lstrip("#").split("\t"), lines[1].split("\t")))


def _e(v):
    return html.escape(str(v if v not in (None, "") else "NA"))


def _rows(d):
    fields = [
        ("Alpha genes measured", d.get("alpha_genes_called"),
         "Functional alpha-globin gene count, 0-4. NA where the measurement "
         "does not resolve to a single count."),
        ("Measurement confidence", d.get("alpha_genes_confidence"),
         "high = depth and junction evidence agree; medium = depth only; "
         "low = a segment sat near a decision boundary, or the channels "
         "disagreed."),
        ("Deletion alleles", d.get("deletion_alleles"),
         "'/'-separated haplotypes. A '|' INSIDE a haplotype means depth cannot "
         "separate those alleles and the group is reported as-is; picking one "
         "would be a population inference, not a measurement."),
        ("Supporting evidence", d.get("deletion_evidence"), ""),
        ("Panel positions found", d.get("site_genotypes"),
         "'none' means every position in the pinned panel was examined and no "
         "variant was seen there. It says nothing about positions not on the "
         "panel."),
        ("Site panel version", d.get("site_panel_version"),
         "Which panel file actually ran. 'No site found' is uninterpretable "
         "without it."),
        ("Integrated genotype", d.get("genotype"),
         "Deletion genotype, then any panel-site findings after '+'. The two "
         "are kept separate because short reads do not establish which "
         "chromosome a site variant sits on."),
    ]
    out = []
    for label, value, note in fields:
        hint = f'<br><small class="text-muted">{_e(note)}</small>' if note else ""
        out.append(f"<tr><th style=\"width:32%\">{_e(label)}</th>"
                   f"<td><strong>{_e(value)}</strong>{hint}</td></tr>")
    return "\n      ".join(out)


def _scope(d):
    not_screened = [t for t in (d.get("not_screened") or "").split(",") if t]
    screened = [t for t in (d.get("screened") or "").split(",") if t]
    items = "\n        ".join(
        f"<li><code>{_e(t)}</code> &mdash; {_e(TIER_TEXT.get(t, 'not examined'))}</li>"
        for t in not_screened) or "<li>none declared</li>"
    return f"""
    <div class="alert alert-warning">
      <strong>What this did NOT examine.</strong> Their absence from the results
      above is <strong>not</strong> a negative result &mdash; nothing was tested,
      so nothing can be ruled out.
      <ul>
        {items}
      </ul>
      <small>Examined: <code>{_e(','.join(screened) or 'none')}</code>.
      This module reports measurements only; it does not interpret them and
      <code>interpretation_complete</code> is
      <code>{_e(d.get('interpretation_complete', 'false'))}</code>.</small>
    </div>"""


def render_html_section(sample_id, contract_path):
    d = parse_contract(contract_path)
    if d is None:
        return (f'<div class="card mb-3"><div class="card-header">'
                f'<h5>Alpha-globin (HBA1/HBA2) &mdash; {html.escape(sample_id)}</h5>'
                f'</div><div class="card-body"><div class="alert alert-secondary">'
                f'The alpha-globin module did not run for this sample. Nothing '
                f'was measured at this locus, so nothing can be ruled out.'
                f'</div></div></div>')
    return f"""
<div class="card mb-3">
  <div class="card-header"><h5>Alpha-globin (HBA1/HBA2) &mdash; {html.escape(sample_id)}</h5></div>
  <div class="card-body">
    {_scope(d)}
    <table class="table table-sm">
      {_rows(d)}
    </table>
  </div>
</div>"""


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", required=True, help="<S>.alpha_globin.tsv")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    Path(a.out).write_text(render_html_section(a.sample, a.tsv))
    print(f"alpha-globin HTML section written to {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
