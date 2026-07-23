"""The alpha-globin HTML card.

The card is where a false negative would actually be READ, so most of these
tests are about what it must NOT say. SVcaller measures; OmniGen interprets.
"""
import os
import re
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))

import hba_report as hr  # noqa: E402

COLS = ["sample", "alpha_genes_called", "alpha_genes_confidence",
        "deletion_alleles", "deletion_evidence", "site_genotypes",
        "site_panel_version", "genotype", "screened", "not_screened",
        "interpretation_complete"]

THAL1_ROW = ["THAL1", "2", "high", "--SEA|--MED/aa", "both", "none",
             "hba_pathogenic_sites.tsv@dfd6ccf", "--SEA|--MED/aa",
             "alpha_deletional,alpha_targeted_sites",
             "beta_globin,alpha_nondeletional_outside_panel", "false"]


def _write(tmp_path, row=THAL1_ROW, name="S.alpha_globin.tsv"):
    p = tmp_path / name
    p.write_text("\t".join(COLS) + "\n" + "\t".join(row) + "\n")
    return str(p)


def _text(html):
    """Visible text, tags stripped and whitespace collapsed — what a reader
    actually sees. Collapsing matters: inline <strong> inside a sentence leaves
    double spaces that would defeat a naive substring check."""
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html)).strip()


def _results_text(html):
    """Visible text of the RESULTS TABLE only.

    The scope block deliberately contains the words "not a negative result" —
    that is the warning doing its job. A blanket ban on those words across the
    whole card would forbid the warning itself, so the dangerous-wording checks
    apply to the results region, where a reader looks for the answer."""
    return _text(html[html.index("<table"):])


def test_renders_the_measurements(tmp_path):
    html = hr.render_html_section("THAL1", _write(tmp_path))
    body = _text(html)
    assert "--SEA|--MED/aa" in body
    assert "hba_pathogenic_sites.tsv@dfd6ccf" in body


def test_scope_is_present_and_explicit_in_the_footer(tmp_path):
    """The card was made compact (2026-07-23): the scope moved from an alert box
    above the table to a muted footer below it. What must NOT change is that the
    scope is explicit and states that an untested locus is not a cleared one — a
    footnote-that-reads-as-'Clear' is exactly the SMN failure this guards against.

    Note the phrasing is 'not the same as tested-and-absent', not 'not a negative
    result': once the footer sits BELOW the table, the word 'negative' would fall
    inside the results region that test_a_normal_sample_is_not_described_as_negative
    forbids. Same meaning, no banned word."""
    html = hr.render_html_section("THAL1", _write(tmp_path))
    assert "alpha-scope-footer" in html, "the scope footer must be present"
    body = _text(html)
    assert "beta" in body.lower() and "HBB" in body
    assert "Not examined" in body
    assert "nothing about them can be ruled out" in body
    assert "not the same as tested-and-absent" in body


def test_every_not_screened_tier_gets_readable_text(tmp_path):
    """A bare tier id in the report is a tier nobody reads — the footer renders
    each not_screened tier's human sentence, not its slug."""
    html = hr.render_html_section("THAL1", _write(tmp_path))
    for tier in ("beta_globin", "alpha_nondeletional_outside_panel"):
        assert hr.TIER_TEXT[tier][:30] in _text(html)


def test_card_never_interprets(tmp_path):
    """No HbH / Bart's / trait classification, no couple risk, and never the
    unqualified word 'thalassemia' in a result line."""
    html = hr.render_html_section("THAL1", _write(tmp_path))
    body = _text(html).lower()
    for token in ("thalassemia", "thalassaemia", "carrier", "hbh", "bart",
                  "affected", "clear", "low risk", "healthy", "diagnos"):
        assert token not in body, f"{token!r} must not appear in the card"


def test_a_normal_sample_is_not_described_as_negative(tmp_path):
    """The single most dangerous string this module could emit."""
    row = list(THAL1_ROW)
    row[1], row[3], row[5], row[7] = "4", "none", "none", "aa/aa"
    html = hr.render_html_section("S", _write(tmp_path, row))
    body = _results_text(html).lower()
    for token in ("negative", "no findings", "nothing found", "normal",
                  "not a carrier", "clear"):
        assert token not in body, f"{token!r} must not appear for a 4-gene sample"


def test_empty_contract_raises_rather_than_rendering_nothing_found(tmp_path):
    """A header-only contract means the integrator ran and produced no row."""
    p = tmp_path / "S.alpha_globin.tsv"
    p.write_text("\t".join(COLS) + "\n")
    with pytest.raises(hr.AlphaReportInputError):
        hr.render_html_section("S", str(p))

    p.write_text("")
    with pytest.raises(hr.AlphaReportInputError):
        hr.render_html_section("S", str(p))


def test_absent_contract_says_nothing_was_measured(tmp_path):
    """Absent (module skipped) is a legitimate state, but it must read as
    'nothing was tested', never as an absence of findings."""
    html = hr.render_html_section("S", str(tmp_path / "NO_FILE"))
    body = _text(html).lower()
    assert "did not run" in body
    assert "nothing can be ruled out" in body
    for token in ("negative", "normal", "clear"):
        assert token not in body


def test_degenerate_group_is_shown_as_the_group(tmp_path):
    """OmniGen must render the group as-is; so must we. Displaying only the
    first member would re-invent the precision the group exists to deny."""
    html = hr.render_html_section("THAL1", _write(tmp_path))
    body = _text(html)
    assert "--SEA|--MED/aa" in body
    assert not re.search(r"(?<![|\-])--SEA/aa", body.replace("--SEA|--MED/aa", ""))


def test_site_hits_are_named_from_the_panel(tmp_path):
    """A raw 'HBA2:c.377:het' is uninterpretable; the card names it from the
    pinned panel. This is the whole reason the panel is loaded."""
    row = list(THAL1_ROW)
    row[5] = "HBA2:c.377:het"     # Hb Quong Sze in assets/hba_pathogenic_sites.tsv
    body = _text(hr.render_html_section("THAL2", _write(tmp_path, row)))
    assert "Hb Quong Sze" in body, "panel-driven site naming did not fire"
    assert "heterozygous" in body
    # and the footer must still list what the scan covered, by name
    assert "Constant Spring" in body


def test_panel_load_is_optional(tmp_path, monkeypatch):
    """If the panel file cannot be found the card must still render (degraded to
    the raw site string), never crash — reports must not depend on asset layout."""
    monkeypatch.setattr(hr, "_load_panel", lambda: ({}, []))
    row = list(THAL1_ROW); row[5] = "HBA2:c.377:het"
    body = _text(hr.render_html_section("S", _write(tmp_path, row)))
    assert "HBA2 c.377" in body and "heterozygous" in body


def test_values_are_html_escaped(tmp_path):
    row = list(THAL1_ROW)
    row[3] = "<script>alert(1)</script>"
    html = hr.render_html_section("S", _write(tmp_path, row))
    assert "<script>" not in html
    assert "&lt;script&gt;" in html


def test_interpretation_complete_is_surfaced(tmp_path):
    html = hr.render_html_section("THAL1", _write(tmp_path))
    assert "interpretation_complete" in html
    assert "false" in _text(html)


def test_cli(tmp_path):
    import subprocess
    out = tmp_path / "card.html"
    r = subprocess.run(
        [sys.executable, os.path.join(REPO, "bin", "hba_report.py"),
         "--tsv", _write(tmp_path), "--sample", "THAL1", "--out", str(out)],
        capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "--SEA|--MED/aa" in out.read_text()
