#!/usr/bin/env python3
"""Regenerate assets/str_disease_ranges.tsv from the pinned STRchive release (assets/sources/, see PROVENANCE.md).

Single cited source for repeat-expansion thresholds (owner decision 2026-11-15).  Counts are in repeat units of the GRADED
block of each locus (assets/eh_graded_variants.tsv, written by bin/make_eh_catalog.py): normal_max = STRchive benign_max,
intermediate_min/max = STRchive intermediate_min/max, path_min = STRchive pathogenic_min.  Nothing is typed by hand.
A locus for which STRchive gives no benign_max / intermediate range leaves that cell BLANK and says so in `notes`; a blank is
"STRchive gives none", never "no limit".

Column meanings (this table belongs to the 2026-11 STRchive catalogue; assets/str_disease_ranges.legacy.tsv belongs to the old one):
  normal_max            largest benign count (STRchive benign_max)
  intermediate_min/max  STRchive's intermediate range CLAMPED to lie strictly between normal_max and path_min
  premut_min/premut_max premutation range, ONLY for loci where STRchive's own text calls its intermediate alleles a premutation
                        (PREMUTATION below, quote pinned by tests/test_str_ranges.py); blank = STRchive defines none, so an
                        intermediate count stays INTERMEDIATE (e.g. TBP: reduced penetrance, not a premutation)
  path_min              smallest pathogenic count (STRchive pathogenic_min)
  repeat_unit           the motif EH COUNTS for the graded block (reference orientation, from assets/eh_catalog.json)
  disease_motif         STRchive's pathogenic motif in gene orientation (documentation only; same unit count)
  inheritance           XD/XR are written XLD/XLR as in the previous table

Usage: make_str_ranges.py [--check]
"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOCI = os.path.join(ROOT, "assets", "sources", "STRchive-loci.json")
GRADED = os.path.join(ROOT, "assets", "eh_graded_variants.tsv")
OUT = os.path.join(ROOT, "assets", "str_disease_ranges.tsv")
RELEASE = "STRchive v2.28.0 (commit e964fa88df3b2189ee478cbd6eaad15bb28e3a0a)"
CATALOG = os.path.join(ROOT, "assets", "eh_catalog.json")
COLS = ["locus_id", "gene", "disease", "inheritance", "repeat_unit", "normal_max", "premut_max", "path_min", "notes",
        "premut_min", "intermediate_min", "intermediate_max", "graded_variant", "disease_motif", "strchive_id", "source"]
INHERITANCE = {"XD": "XLD", "XR": "XLR"}      # previous table's codes
# loci where STRchive's `details` calls intermediate alleles a premutation: locus -> (premut_min, premut_max, evidence).
# None = take STRchive's (clamped) intermediate bound.  Each evidence sentence is in STRchive-loci.json `details`.
PREMUTATION = {
    "FMR1": (55, None, "gray zone 45-54 'may expand into the premutation range'; premutation = intermediate above the gray zone"),
    "AFF2": (None, None, "'Intermediate alleles correspond to a premutation' [pmid:23914978]"),
    "ATXN1": (None, None, "'intermediate alleles are considered premutations' [genereviews:NBK1184]"),
    "ATXN7": (None, 33, "'premutations (28-33) to reduced penetrance (34-36)' [genereviews:NBK1256]"),
    "CACNA1A": (None, None, "'intermediate range (19-20 motifs) can be associated with a premutation' [genereviews:NBK1140]"),
    "CNBP": (30, 74, "'~30-~54 repeats are considered premutations, ~55-74 premutations with possible reduced penetrance' [genereviews:NBK1466]"),
    "DMPK": (None, None, "'Intermediate alleles of 35-49 repeats are associated with premutation' [genereviews:NBK1165]"),
    "FXN": (None, None, "'Intermediate alleles (34-55) are associated with premutations' [genereviews:NBK1281]"),
    "HTT": (None, None, "'27-35 motifs are unstable/premutations' [pmid:39572770]"),
}


def counted_motif(entry, variant_id):
    """Motif of the graded block as the catalog models it (what EH counts)."""
    import re
    motifs = re.findall(r"\(([A-Z]+)\)\*", entry["LocusStructure"])
    ids = entry.get("VariantId") or [entry["LocusId"]]
    return motifs[ids.index(variant_id)]


def build():
    loci = {x["id"]: x for x in json.load(open(LOCI))}
    cat = {e["LocusId"]: e for e in json.load(open(CATALOG))}
    rows = []
    for line in open(GRADED):
        if line.startswith("#") or line.startswith("locus\t"):
            continue
        g = dict(zip(["locus", "variant_id", "strchive_id", "disease_motif"], line.rstrip("\n").split("\t")))
        x = loci[g["strchive_id"]]
        blank = lambda v: "" if v is None else str(v)
        notes = []
        imin, imax = x["intermediate_min"], x["intermediate_max"]
        if imin is not None:                            # clamp to (normal_max, path_min), exclusive
            if x["benign_max"] is not None:
                imin = max(imin, x["benign_max"] + 1)
            if x["pathogenic_min"] is not None:
                imax = min(imax, x["pathogenic_min"] - 1)
        if x["benign_max"] is None:
            notes.append("STRchive gives no benign/normal maximum")
        if x["intermediate_min"] is None:
            notes.append("STRchive gives no intermediate range")
        else:
            if (imin, imax) != (x["intermediate_min"], x["intermediate_max"]):
                notes.append(f"STRchive intermediate {x['intermediate_min']}-{x['intermediate_max']} overlaps its "
                             f"normal/pathogenic bound; clamped to {imin}-{imax}")
        pmin = pmax = ""
        if g["locus"] in PREMUTATION:
            a, b, why = PREMUTATION[g["locus"]]
            pmin = "" if a is None else str(a)
            pmax = str(imax if b is None else b)
            notes.append(f"premutation range per STRchive text: {why}")
        if g["variant_id"] != g["locus"]:
            notes.append(f"compound locus: counts are repeats of the {g['variant_id']} block ({g['disease_motif']})")
        refs = x.get("references") or []
        src = f"{RELEASE} {g['strchive_id']}; https://strchive.org/loci/{g['strchive_id']}/; {len(refs)} refs, first: " + \
              ",".join(refs[:3])
        rows.append({"locus_id": g["locus"], "gene": x["gene"], "disease": x["disease"], "inheritance": "/".join(INHERITANCE.get(c, c) for c in x["inheritance"]),
                     "repeat_unit": counted_motif(cat[g["locus"]], g["variant_id"]),
                     "normal_max": blank(x["benign_max"]), "premut_max": pmax, "path_min": blank(x["pathogenic_min"]),
                     "notes": "; ".join(notes), "premut_min": pmin, "intermediate_min": blank(imin),
                     "intermediate_max": blank(imax), "graded_variant": g["variant_id"],
                     "disease_motif": (x["pathogenic_motif_gene_orientation"] or x["reference_motif_reference_orientation"])[0],
                     "strchive_id": g["strchive_id"], "source": src})
    return rows


def render(rows):
    return "\t".join(COLS) + "\n" + "".join("\t".join(r[c] for c in COLS) + "\n" for r in rows)


def main(argv):
    txt = render(build())
    if "--check" in argv:
        sys.exit(0 if open(OUT).read() == txt else "ERROR: committed ranges differ from STRchive; rerun bin/make_str_ranges.py")
    open(OUT, "w").write(txt)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main(sys.argv)
