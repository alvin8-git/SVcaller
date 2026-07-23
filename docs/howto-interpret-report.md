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

### 8. Blood Group & Copy-Number Traits

Populated from the `CNV_TRAITS` subworkflow (`results/<S>/bloodgroup/*.rh_status.tsv`
and `results/<S>/cnv_traits/*.{amy1,gst_null,lpa_kiv2}.tsv`). Each call is a
**targeted, normalized read-depth** estimate: mosdepth (`--mapq 0`, multi-mapping
reads kept) over the trait window divided by the median depth of the copy-number-stable
`CTRL_*` control windows (the diploid 2n baseline). The CNV consensus BED is used
only to *corroborate* a called deletion — it under-calls these paralogous
whole-gene deletions, so depth leads.

| Trait | Reported as | Read this way |
|-------|-------------|---------------|
| Rh(D) blood group | pos / neg + RHD copies + confidence | RHD CN≈2 ⇒ Rh(D) positive; homozygous RHD deletion (CN→0) ⇒ Rh(D) negative. RHD depth runs slightly high from RHCE paralog cross-mapping. |
| AMY1 | copies | Salivary-amylase array size. **Proportional-only until calibrated** — the window averages the array against flanking single-copy sequence, so absolute copies are approximate. |
| GSTM1 / GSTT1 | present / null | `null` = homozygous whole-gene deletion. **A ratio near 0.5 is a heterozygous deletion (1 copy) and reports as `present`** — the current call is binary and does not distinguish het-del from 2-copy. |
| Lp(a) — LPA KIV-2 | repeat copies | Cardiovascular-risk VNTR. **Proportional-only until calibrated** — the constant depends on how many KIV-2 units the reference window already spans. |

**Worked example — HG002 (validated 2026-07-12).** Control baseline 32.31X (matches
the pipeline's ~31.6X autosomal mean). Rh(D) **positive** (RHD 36.4X, CN 2; consensus
BED has no RHD deletion — the two channels agree). GSTT1 **present** (32.5X, ratio 1.00).
GSTM1 reads **present** at ratio 0.46 (~1 copy): GIAB v5.0q truth confirms a heterozygous
18.4 kb GSTM1 deletion (`GT 0|1`), i.e. one functional copy — so "present" is correct and
this is *not* a homozygous null. AMY1=2 and KIV-2=5 are the uncalibrated depth estimates.
See `docs/omnigen-additions-plan.md` for the full validation table.

**Worked example — HG001 / NA12878 (validated 2026-07-12).** Control baseline 30.81X
(matches the pipeline's ~29.9X autosomal mean). Rh(D) **positive** (RHD 19.0X → CN 1,
confidence MEDIUM; consensus BED has no RHD deletion). GSTT1 **present** (31.4X, ratio 1.02).
**GSTM1 caveat (settled by a paralog-aware PSV test).** NA12878 is often cited as *the*
canonical homozygous GSTM1-null, but this BAM reads GSTM1 at 16.8X, **ratio 0.545 (~1 copy)**.
Depth alone can't tell a het deletion from a homozygous null, because GSTM1's near-identical
paralog **GSTM2** can cross-map in; so we checked paralog-specific variants (fixed
GSTM1-vs-GSTM2 differences). The result is unambiguous **heterozygous deletion, not null**:
99% of window reads (611/617) are MAPQ≥20 (uniquely GSTM1, not ambiguous cross-mappers); at
39 GSTM1-vs-GSTM2 discriminating sites the reads are **433 GSTM1-allele : 0 GSTM2-allele**;
and the 11 `1/1` calls sit off the PSVs and don't match GSTM2 bases (real GSTM1 hemizygous
variants). So `present` is the correct call **for the data**, and the expected homozygous-null
is a documented **mismatch** with this alignment — the numbers were not forced. (`gst_null.py`
needs no paralog fix here: cross-mapping is negligible, so a true null would read ~0 depth and
fire `null` correctly.) AMY1=2, KIV-2=6 uncalibrated. Neither HG002 nor HG001 exercises the
`null` branch — a genuine homozygous-null sample is still needed to validate it end-to-end.

### 9. Alpha-globin (HBA1/HBA2)

A succinct card from the M8 alpha-globin measurement (`bin/hba_report.py`, wired into
`subworkflows/report.nf`). It reports measurements only and makes no clinical
interpretation. The layout mirrors the SMN card: a headline, a three-row body, and a
muted footer.

| Field | What it shows |
|-------|---------------|
| Alpha-gene dosage | Called alpha-gene count as "N of 4 genes", or "not resolved to a single count" when depth is ambiguous. |
| Genotype | The measured genotype string (e.g. `aa/aa`), with a one-clause measurement note (no large deletion detected, one chromosome carries a deletion, or both do). |
| Deletion | The deletion genotype, keeping the FULL group string. When depth cannot tell two alleles apart, it is reported as a GROUP (e.g. `--SEA|--MED`), never collapsed to one allele. |
| Supporting evidence | What backed the deletion call: read-depth, junction read, or both. |
| Point-mutation scan | Named hits from the pinned site panel (e.g. `HBA2:c.377` is named "Hb Quong Sze"). Only the panel positions were scanned. |

The footer lists what was screened (alpha-gene dosage and the named point mutations) and,
explicitly, what was NOT examined. Absence in this card is not the same as
tested-and-absent: anything outside the panel was not looked at, so nothing about it is
ruled out. The footer also records the panel version and that no clinical interpretation
is made here.

## Verification

The report was generated correctly if:
- File size is 1.5-3 MB (larger suggests SVG rendering issue)
- All 7 sections are present and non-empty
- Circos plot renders inline (not a broken image)
- Coverage depth is reported as a number (not "N/A")

## Troubleshooting

**Mean coverage or mapped-reads percentage shows N/A**
The Alignment QC section shows mean coverage (mosdepth), duplicate rate (Picard), and mapped-reads percentage (samtools flagstat). All three are wired. mosdepth summary, flagstat, and insert-size now publish to `results/<sample>/qc/`, so the QC section is reproducible after a work-dir cleanup. Look for `results/<sample>/qc/<sample>.flagstat.txt`. The mapped-reads percentage is the true alignment rate: flagstat runs on the input BAM before FILTER_CHROMS, so it is not a post-filter 100%. N/A now only appears if the `qc/` directory is missing (an interrupted run, or the sample never reached the flagstat step); re-run the pipeline to regenerate it.

**SV Summary and Top Annotated SVs sections are empty**
AnnotSV failed or produced 0-byte output. Check `--annotsv_db` path — it must point to the parent of `Annotations_Human/`, not to `Annotations_Human/` itself.

**Circos plot is not visible / broken**
The inline SVG may have failed to render. Check `{sample}.circos.png` directly in the output directory as a fallback.

## Related

- [Architecture reference — M6/M7 Report](reference-architecture.md#m6m7-report-subworkflowsreportnf)
- [How to run the GIAB validation](howto-run-validation.md)
- [Parameter reference](reference-parameters.md)
