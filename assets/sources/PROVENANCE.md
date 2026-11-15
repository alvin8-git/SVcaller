# STRchive source catalogs (pinned)

| file | what |
|---|---|
| `STRchive-disease-loci.hg38.expansionhunter.json` | STRchive's ExpansionHunter catalog, GRCh38 |
| `STRchive-disease-loci.hg38.stranger.json` | STRchive's Stranger catalog, GRCh38 (adds VariantId, PathologicRegion, NormalMax, PathologicMin, DisplayRU) |

| `STRchive-loci.json` | STRchive locus database (benign/intermediate/pathogenic ranges, disease, references); same release, `data/STRchive-loci.json` |
| `Illumina-ExpansionHunter.grch38.variant_catalog.json` | Illumina/ExpansionHunter `variant_catalog/grch38/variant_catalog.json`, commit `7ff9aa5335fdb5b3124cc0f365a64bb7a185803a`; used ONLY for AR and AFF2 regions (owner decision 2026-11-15, see `assets/eh_locus_map.tsv`); its coordinates have no `chr` prefix |

- Upstream: https://github.com/dashnowlab/STRchive, directory `data/catalogs/`
- Release: **v2.28.0** (2026-10-08), commit `e964fa88df3b2189ee478cbd6eaad15bb28e3a0a`
- Files are byte-for-byte the upstream files at those commits (sha256 in the tests).
- `python3 bin/make_str_ranges.py` regenerates `assets/str_disease_ranges.tsv`.
- ExpansionHunter catalog syntax (multi-block loci, `VariantId`, 0-based half-open `ReferenceRegion`):
  https://github.com/Illumina/ExpansionHunter/blob/master/docs/04_VariantCatalogFiles.md
- Regenerate `assets/eh_catalog.json` and `assets/eh_graded_variants.tsv` with `python3 bin/make_eh_catalog.py`.
- To update STRchive: replace the two json files from a newer tag, update this file and the sha256 in the test, rerun the script.
