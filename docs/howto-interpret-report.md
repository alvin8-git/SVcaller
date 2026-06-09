# How to Interpret the SVcaller HTML Report

Open `{outdir}/{sample}/{sample}.report.html` in any modern browser. The report has seven sections, described below.

## Prerequisites

- A completed pipeline run producing `{sample}.report.html`
- Basic familiarity with structural variant (SV) terminology: deletion (DEL), duplication (DUP), inversion (INV), translocation (BND/TRA), mobile element insertion (MEI/INS)
- For clinical interpretation: access to OMIM and ClinVar for variant follow-up

## Sections

### 1. Alignment QC

| Field | Source | What to look for |
|-------|--------|-----------------|
| Mean coverage | mosdepth | Should be ≥30×. Values below 25× reduce SV sensitivity, especially for Scramble MEI calls. |
| Duplicate rate | Picard MarkDuplicates | <10% is normal for WGS. >20% suggests library quality issues. |
| Mapped reads | samtools flagstat | >95% is expected. Lower values indicate sample quality problems or reference mismatch. |

### 2. Structural Variant Summary

Table of SV counts by type (DEL, DUP, INV, BND, INS) with the count of class 4/5 (Pathogenic/Likely Pathogenic) calls highlighted.

**What the numbers mean:**
- A normal 30× WGS sample typically has 8,000-15,000 SV calls after merging all callers
- After gnomAD SV frequency filtering (`filtered.tsv`), expect 5,000-12,000 calls
- High INV and TRA counts (>5,000 each) from a single caller (especially GRIDSS) are often artefacts — filter by `SUPP_VEC ≥ 2` to see multi-caller-supported calls only

### 3. Genome-wide Circos Plot

The Circos plot shows the full genome at a glance. Rings from outside to inside:

| Ring | What it shows | What to look for |
|------|--------------|-----------------|
| Chromosomes | Chr1-22+X+Y | Orientation reference |
| Coverage depth | Log₂ depth ratio vs. median (50 kb windows) | Red blocks = copy gain region; blue blocks = copy loss. For a diploid sample, this should be uniformly grey. |
| STR loci | ExpansionHunter repeat expansions | Brown bars at known loci (HTT, FMR1, C9orf72, etc.) |
| Gene loci | Top 30 SVs by AnnotSV ranking score | Coloured bars by SV type; gold bar = SMN locus on chr5 |
| ACMG class dots | Class 5 (red), 4 (orange), 3 (grey) | Red/orange dots indicate potentially pathogenic SVs |
| SV links | Multi-caller SVs ≥50 kb | Coloured arcs connecting breakpoints. Dense chr regions or inter-chromosomal links warrant review. |

**Normal appearance:** Uniform grey depth ring, sparse SV links, no red/orange ACMG dots.

**Abnormal patterns to investigate:**
- Large red/blue block in depth ring → large CNV; cross-reference with `cnv_consensus.bed`
- Red/orange ACMG dot in a known disease locus → check **Top Annotated SVs** table
- Dense SV links from one chromosome → possible chromothripsis or complex rearrangement

### 4. SMN Copy Number

SMN1 and SMN2 copy number results from SMNCopyNumberCaller. Normal is SMN1=2, SMN2=2.

| Finding | Interpretation |
|---------|---------------|
| SMN1=0, SMN2=4 | SMA-affected (homozygous deletion of SMN1) |
| SMN1=1, SMN2=3 | SMA carrier |
| SMN1=2, SMN2=2 | Normal |
| SMN1≥3 | Confirm with orthogonal method — high CN calls have higher uncertainty |

### 5. STR Expansion Loci

Table of ExpansionHunter results for the catalogued repeat loci. Check `REPCN` (repeat count) against established pathogenic thresholds:

- **HTT (Huntington's):** >36 repeats is pathogenic
- **FMR1 (Fragile X):** >200 repeats (full mutation), 55-200 (premutation)
- **C9orf72 (ALS/FTD):** >30 repeats is a risk factor

Repeat counts above threshold should be confirmed with PCR or long-read sequencing.

### 6. Top Annotated SVs

Top 20 SVs sorted by AnnotSV ranking score (higher = more likely pathogenic). Columns:

| Column | Meaning |
|--------|---------|
| Chr / Start / End | Genomic coordinates |
| Type | SV type (DEL, DUP, INV, BND, INS) |
| Size | Formatted size (bp, kb, Mb) |
| Gene | First overlapping gene (from AnnotSV) |
| Score | AnnotSV ranking score |
| Class | ACMG class (1=benign → 5=pathogenic) |
| OMIM | OMIM morbid gene status |

**Class 4/5 SVs should always be reviewed.** Check the full `filtered.tsv` for additional context (multiple transcripts, population frequency, loss-of-function status).

### 7. GIAB Benchmark (validation runs only)

Appears only when `--giab_truth` was provided. Shows Truvari precision, recall, and F1 for the overall callset and four size bins (50-300 bp, 300 bp-1 kb, 1-10 kb, >10 kb).

F1 score interpretation:
- F1 ≥ 0.5: acceptable for research use
- F1 ≥ 0.4: typical for multi-caller WGS pipelines on HG002
- F1 < 0.3: investigate caller failures or reference mismatch

## Verification

The report was generated correctly if:
- File size is 1.5-3 MB (larger suggests SVG rendering issue)
- All 7 sections are present and non-empty
- Circos plot renders inline (not a broken image)
- Coverage depth is reported as a number (not "N/A")

## Troubleshooting

**Mean coverage or mapped-reads percentage shows N/A**
Check that `{sample}.flagstat.txt` was published to the output directory. The Alignment QC section shows mean coverage (mosdepth), duplicate rate (Picard), and mapped-reads percentage (samtools flagstat). If flagstat is missing, re-run the pipeline; flagstat is produced during the PREPROCESS step alongside MarkDuplicates.

**SV Summary and Top Annotated SVs sections are empty**
AnnotSV failed or produced 0-byte output. Check `--annotsv_db` path — it must point to the parent of `Annotations_Human/`, not to `Annotations_Human/` itself.

**Circos plot is not visible / broken**
The inline SVG may have failed to render. Check `{sample}.circos.png` directly in the output directory as a fallback.

## Related

- [Architecture reference — M6/M7 Report](reference-architecture.md#m6m7-report-subworkflowsreportnf)
- [How to run the GIAB validation](howto-run-validation.md)
- [Parameter reference](reference-parameters.md)
