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
import re
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


# The pinned point-mutation panel. Loaded only to NAME site hits ("HBA2:c.377"
# -> "Hb Quong Sze") and to list what the scan covered, instead of dumping raw
# coordinates a reader cannot interpret. Source of truth is
# assets/hba_pathogenic_sites.tsv; if it cannot be found the card degrades to the
# raw site string rather than failing.
def _load_panel():
    for base in (Path(__file__).resolve().parent.parent / "assets",
                 Path("/usr/local/assets")):
        p = base / "hba_pathogenic_sites.tsv"
        if not p.exists():
            continue
        by_site, ordered = {}, []
        for line in p.read_text().splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 3 or f[0] == "gene":
                continue
            gene, allele, hgvs = f[0], f[1], f[2]
            by_site[(gene, hgvs)] = allele
            base_name = re.sub(r"\s*\(HBA1\)\s*$", "", allele)  # de-dup the paralogue
            if base_name not in ordered:
                ordered.append(base_name)
        return by_site, ordered
    return {}, []


def _dosage(d):
    n = (d.get("alpha_genes_called") or "").strip()
    return "not resolved to a single count" if n in ("", "NA") else f"{_e(n)} of 4 genes"


_EVIDENCE = {"both": "read-depth + junction read", "depth": "read-depth",
             "junction": "junction read", "none": "&mdash;", "": "&mdash;"}

_ALPHA_AA = ("aa", "αα")


def _describe_deletion(alleles):
    """Plain description of the deletion genotype, keeping the FULL group string.
    Never collapses a '|' group to one allele — that would invent precision."""
    alleles = (alleles or "").strip()
    haps = alleles.split("/")
    dels = [h for h in haps if h and h != "none" and h not in _ALPHA_AA]
    if not dels or alleles in ("", "none"):
        return "None detected by depth."
    parts = []
    for h in dels:
        two_gene = "--" in h
        kind = "&alpha;0 (two-gene)" if two_gene else "&alpha;+ (one-gene)"
        grp = " &mdash; subtype indistinguishable by depth" if "|" in h else ""
        parts.append(f"{kind} deletion{grp}")
    zyg = ("heterozygous" if len(dels) == 1 and len(dels) < len(haps)
           else "homozygous" if len(dels) == 2 and len(set(dels)) == 1
           else "biallelic")
    return f"{'; '.join(parts)}, {zyg}. Full call: <code>{_e(alleles)}</code>."


def _describe_sites(genotypes, by_site, ordered):
    n = len(by_site) or "the pinned"   # genomic positions screened (Adana is 2)
    g = (genotypes or "").strip()
    if g in ("", "none"):
        return f"{n} sites screened; none detected."
    out = []
    for entry in g.split(","):
        bits = entry.strip().split(":")
        if len(bits) >= 3:
            gene, hgvs, zyg = bits[0], bits[1], bits[2]
            name = by_site.get((gene, hgvs))
            label = f"{name} ({gene} {hgvs})" if name else f"{gene} {hgvs}"
            zyg = {"het": "heterozygous", "hom": "homozygous"}.get(zyg, zyg)
            out.append(f"{_e(label)}, {_e(zyg)}")
        elif entry.strip():
            out.append(_e(entry.strip()))
    return "; ".join(out) + "."


def _genotype_note(d):
    """One measurement-level clause after the genotype — never a clinical call."""
    alleles = (d.get("deletion_alleles") or "").strip()
    haps = alleles.split("/")
    dels = [h for h in haps if h and h != "none" and h not in _ALPHA_AA]
    if not dels or alleles in ("", "none"):
        return " &mdash; no large deletion detected by depth"
    if any(h in _ALPHA_AA for h in haps):
        return " &mdash; one chromosome carries a deletion, the other is intact for large deletions"
    return " &mdash; both chromosomes carry a deletion"


def render_html_section(sample_id, contract_path):
    d = parse_contract(contract_path)
    if d is None:
        return (f'<div class="card mb-3"><div class="card-header">'
                f'<h5>Alpha-globin (HBA1/HBA2) &mdash; {html.escape(sample_id)}</h5>'
                f'</div><div class="card-body"><div class="alert alert-secondary">'
                f'The alpha-globin module did not run for this sample. Nothing '
                f'was measured at this locus, so nothing can be ruled out.'
                f'</div></div></div>')

    by_site, ordered = _load_panel()
    conf = _e(d.get("alpha_genes_confidence") or "NA")
    site_list = ", ".join(ordered) if ordered else "a pinned set of positions"
    not_screened = [t for t in (d.get("not_screened") or "").split(",") if t]
    scope_items = "; ".join(TIER_TEXT.get(t, t) for t in not_screened) or "none declared"

    # Headline (dosage + genotype) mirrors the SMN card: the answer first, in
    # plain terms. Everything technical — the panel, its version, what was NOT
    # looked at — moves to a muted footer. The footer still states, explicitly,
    # that an untested locus is not a cleared one; that guard is non-negotiable
    # even in the compact layout (an empty SMN artifact once read as "Clear").
    return f"""
<div class="card mb-3">
  <div class="card-header"><h5>Alpha-globin (HBA1/HBA2) &mdash; {html.escape(sample_id)}</h5></div>
  <div class="card-body">
    <p class="mb-1" style="font-size:1.08rem"><strong>&alpha;-globin dosage: {_dosage(d)}</strong>
       <span class="badge bg-secondary align-middle">confidence: {conf}</span></p>
    <p class="mb-3">Genotype: <strong>{_e(d.get('genotype') or 'NA')}</strong><span class="text-muted">{_genotype_note(d)}</span></p>
    <table class="table table-sm mb-3">
      <tr><th style="width:34%">Deletion</th><td>{_describe_deletion(d.get('deletion_alleles'))}</td></tr>
      <tr><th>Supporting evidence</th><td>{_EVIDENCE.get((d.get('deletion_evidence') or '').strip(), _e(d.get('deletion_evidence')))}</td></tr>
      <tr><th>Point-mutation scan</th><td>{_describe_sites(d.get('site_genotypes'), by_site, ordered)}</td></tr>
    </table>
    <div class="small text-muted border-top pt-2" data-role="alpha-scope-footer">
      <div><strong>Screened:</strong> &alpha;-gene dosage and named point mutations ({_e(site_list)}).</div>
      <div><strong>Not examined:</strong> {_e(scope_items)}
        Absence above is not the same as tested-and-absent &mdash; these were not
        examined here, so nothing about them can be ruled out.</div>
      <div>Measurements only; no clinical interpretation is made here
        (<code>interpretation_complete={_e(d.get('interpretation_complete', 'false'))}</code>).
        Panel <code>{_e(d.get('site_panel_version') or 'NA')}</code>.</div>
    </div>
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
