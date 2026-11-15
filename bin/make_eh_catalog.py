#!/usr/bin/env python3
"""Regenerate assets/eh_catalog.json and assets/eh_graded_variants.tsv from the pinned STRchive hg38 catalogs.

Nothing in the output is typed by hand: structure, coordinates and VariantIds come from
assets/sources/STRchive-disease-loci.hg38.{expansionhunter,stranger}.json (see assets/sources/PROVENANCE.md);
the only local input is assets/eh_locus_map.tsv (OmniGen LocusId -> STRchive id).

eh_catalog.json    ExpansionHunter's variant catalog.  LocusId = OmniGen's name.  A multi-block (compound) locus
                   carries one Variant per block with STRchive's VariantIds.
eh_graded_variants.tsv   per locus, WHICH variant is the disease-relevant block (STRchive PathologicRegion), its motif
                   in disease orientation (DisplayRU) and STRchive's NormalMax / PathologicMin for that block.  bin/eh_sizing.py
                   reads it so thresholds are compared in repeat units of that block, never the whole tract.

Usage: make_eh_catalog.py [--check]   (--check exits 1 if the committed outputs differ from what the sources produce)
"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "assets", "sources")
EH_SRC = os.path.join(SRC, "STRchive-disease-loci.hg38.expansionhunter.json")
ST_SRC = os.path.join(SRC, "STRchive-disease-loci.hg38.stranger.json")
ILL_SRC = os.path.join(SRC, "Illumina-ExpansionHunter.grch38.variant_catalog.json")
MAP = os.path.join(ROOT, "assets", "eh_locus_map.tsv")
OUT_CAT = os.path.join(ROOT, "assets", "eh_catalog.json")
OUT_GRADED = os.path.join(ROOT, "assets", "eh_graded_variants.tsv")
GRADED_COLS = ["locus", "variant_id", "strchive_id", "disease_motif", "n_variants",
               "strchive_normal_max", "strchive_path_min", "inheritance", "region", "all_variant_ids"]


def read_map():
    out = []
    for line in open(MAP):
        line = line.rstrip("\n")
        if not line.strip() or line.startswith("#") or line.startswith("locus\t"):
            continue
        f = line.split("\t")
        out.append((f[0], f[1], f[2] if len(f) > 2 else ""))
    return out


def build():
    eh = {e["LocusId"]: e for e in json.load(open(EH_SRC))}
    st = {e["LocusId"]: e for e in json.load(open(ST_SRC))}
    catalog, graded = [], []
    ill = {e["LocusId"]: e for e in json.load(open(ILL_SRC))}
    for locus, sid, region_source in read_map():
        if sid not in eh or sid not in st:
            raise SystemExit(f"ERROR: {locus} -> {sid} is not in the pinned STRchive catalogs")
        e, s = eh[sid], st[sid]
        if (e["LocusStructure"], e["ReferenceRegion"]) != (s["LocusStructure"], s["ReferenceRegion"]):
            raise SystemExit(f"ERROR: STRchive's EH and Stranger catalogs disagree on {sid}")
        regions = e["ReferenceRegion"] if isinstance(e["ReferenceRegion"], list) else [e["ReferenceRegion"]]
        entry = {"LocusId": locus, "LocusStructure": e["LocusStructure"], "ReferenceRegion": e["ReferenceRegion"],
                 "VariantType": e["VariantType"]}
        if region_source == "illumina":                  # owner-approved override (assets/eh_locus_map.tsv)
            il = ill[locus]
            reg = il["ReferenceRegion"]
            entry["LocusStructure"] = il["LocusStructure"]
            entry["ReferenceRegion"] = reg if reg.startswith("chr") else "chr" + reg
            regions = [entry["ReferenceRegion"]]
        elif region_source:
            raise SystemExit(f"ERROR: {locus}: unknown region_source {region_source!r}")
        if len(regions) > 1:
            vids = s["VariantId"]
            entry["VariantId"] = vids
            vid = vids[regions.index(s["PathologicRegion"])]
        else:
            vid = locus                                 # EH names a single-variant locus's variant after the locus
        catalog.append(entry)
        regs = entry["ReferenceRegion"] if isinstance(entry["ReferenceRegion"], list) else [entry["ReferenceRegion"]]
        vids_all = entry.get("VariantId") or [locus]
        graded.append({"region": regs[vids_all.index(vid)], "all_variant_ids": ",".join(vids_all),"locus": locus, "variant_id": vid, "strchive_id": sid, "disease_motif": s["DisplayRU"],
                       "n_variants": len(regions), "strchive_normal_max": s["NormalMax"],
                       "strchive_path_min": s["PathologicMin"], "inheritance": "/".join(s["InheritanceMode"])})
    return catalog, graded


def render_catalog(catalog):
    return "[\n" + ",\n".join("  " + json.dumps(c) for c in catalog) + "\n]\n"


def render_graded(graded):
    head = ("# variant EH reports for each locus that is graded; derived from STRchive by bin/make_eh_catalog.py, do not edit\n"
            "# disease_motif = STRchive DisplayRU (disease orientation; EH reports the reference-orientation motif, same unit count)\n"
            "# region = ReferenceRegion of the graded block (an EH VCF record has POS=start, END=end of it); all_variant_ids = every variant EH reports for the locus\n")
    return head + "\t".join(GRADED_COLS) + "\n" + "".join("\t".join(str(g[c]) for c in GRADED_COLS) + "\n" for g in graded)


def main(argv):
    catalog, graded = build()
    cat_txt, gr_txt = render_catalog(catalog), render_graded(graded)
    if "--check" in argv:
        same = open(OUT_CAT).read() == cat_txt and open(OUT_GRADED).read() == gr_txt
        sys.exit(0 if same else "ERROR: committed catalog differs from the STRchive sources; rerun bin/make_eh_catalog.py")
    open(OUT_CAT, "w").write(cat_txt)
    open(OUT_GRADED, "w").write(gr_txt)
    print(f"wrote {OUT_CAT} ({len(catalog)} loci) and {OUT_GRADED}")


if __name__ == "__main__":
    main(sys.argv)
