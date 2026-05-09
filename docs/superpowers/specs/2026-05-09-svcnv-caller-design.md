# SVcaller: Structural Variant & Copy Number Variant Pipeline — Design Spec

**Date:** 2026-05-09  
**Status:** Approved  
**Reference genome:** GRCh38/hg38  
**Workflow engine:** Nextflow DSL2 v25.10.4  
**Containers:** Docker  

---

## 1. Overview

SVcaller is a modular Nextflow DSL2 pipeline for detecting structural variants (SVs) and copy number variants (CNVs) from human whole-genome sequencing data (Illumina PE150, 30x minimum coverage). It includes a dedicated SMN1/SMN2 copy number module for SMA clinical classification, a Circos-based genome-wide visualization, and a GIAB benchmarking framework.

The pipeline operates in **single-sample mode** — each sample is processed independently. A separate one-time workflow builds a Panel of Normals (PoN) from GIAB reference samples for improved CNV sensitivity.

---

## 2. Goals

- Detect all major SV classes: DEL, INS, INV, DUP, BND/TRA, and STR expansions
- Call genome-wide CNVs using an ensemble of depth-based callers + GIAB PoN
- Accurately determine SMN1 and SMN2 copy numbers from 30x WGS, including 2+0 haplotype detection
- Validate against GIAB truth sets (HG001–HG007) using truvari
- Produce a per-sample HTML report with embedded Circos plot, CNV tracks, SMN copy number, and QC metrics
- Be architecturally extensible to Nanopore long-read input without modifying core modules

## 3. Non-Goals

- Somatic (tumor/normal) SV calling
- Exome or targeted panel input
- SNP/indel calling (handled by separate pipeline)
- Trio joint calling (each trio member is validated independently)
- Nanopore implementation in v1 (architecture supports it; code deferred)

---

## 4. Architecture

### 4.1 Module Map

```
Input: FASTQ ──┐
               ├──→ [M1: PREPROCESS]
Input: BAM ────┘     BWA-MEM2 align (if FASTQ)
                      Picard MarkDuplicates
                      mosdepth coverage QC (fail < 30x)
                      samtools flagstat + FastQC
                              │
             ┌────────────────┼─────────────────┐
             ↓                ↓                  ↓
  [M2: SV CALLING]  [M3: CNV CALLING]   [M4: SMN1/SMN2]
  Manta v1.6         CNVpytor v1.3       SMNCopyNumberCaller v1.1
  DELLY v1.2.6       GATK gCNV v4.5      (per-sample, independent)
  GRIDSS v2.13       (cohort PoN mode)
  ExpansionHunter v5
  → JASMINE merge    → consensus CNV
             │                │                  ↓
             └────────────────┘           SMN1/2 copy numbers
                      ↓                   SMA classification
             [M5: ANNOTATION]              2+0 haplotype flag
              AnnotSV v3.4
              gnomAD-SV v2.1 freq filter
                      ↓
          ┌───────────┴──────────────┐
          ↓                          ↓
  [M6: VISUALIZATION]       [M7: REPORTING]
  pycirclize Circos SVG       Jinja2 HTML report
  CNV depth plots             MultiQC QC summary
  SMN barplot                 truvari GIAB benchmark
```

### 4.2 Parallel Execution

Modules M2, M3, and M4 execute in parallel on the same preprocessed BAM. Nextflow channels fork at the end of M1 and merge into M5.

### 4.3 PoN Build Workflow (one-time, offline)

A separate `pon_build.nf` workflow runs once on HG001–HG007:
- GATK `PreprocessIntervals` → `CollectReadCounts` (per sample) → `CreateReadCountPanelOfNormals`
- Output: `giab_pon.hdf5` — reused by all future sample runs via `--pon` parameter

---

## 5. Tool Selection

### M1 — Pre-processing

| Step | Tool | Version | Notes |
|------|------|---------|-------|
| Alignment | BWA-MEM2 | 2.2.1 | 3× faster than BWA; best PE150 accuracy |
| Sort/index | samtools | 1.20 | Standard |
| Deduplication | Picard MarkDuplicates | 3.2 | Required for depth-based CNV callers |
| Coverage QC | mosdepth | 0.3.8 | Pipeline halts with error if mean depth < 30x |
| Raw read QC | FastQC + MultiQC | latest | Standard metrics |

### M2 — SV Calling Ensemble

| Tool | Version | SV Types | Rationale |
|------|---------|---------|-----------|
| Manta | 1.6.0 | DEL, INS, INV, DUP, BND | Fast, well-validated, Illumina-maintained |
| DELLY | 1.2.6 | DEL, INS, INV, DUP, TRA | Best for tandem dups; used for joint genotyping |
| GRIDSS | 2.13.2 | All + complex | Assembly-based; highest sensitivity for BNDs and translocations |
| ExpansionHunter | 5.0.0 | STR expansions | Known loci catalog: FMR1, HTT, ATXN1/3, AR, RFC1, C9orf72, 60+ |
| JASMINE | 1.1.5 | Merge | SUPP≥2 filter for consensus; handles complex overlapping calls |

**Merge logic:** Each caller produces a per-sample VCF. JASMINE merges all three structural callers (Manta + DELLY + GRIDSS) into a consensus VCF. Calls supported by ≥2 callers are flagged HIGH confidence. ExpansionHunter output is kept separate as a supplementary STR VCF.

### M3 — CNV Calling

| Tool | Version | Method | Notes |
|------|---------|--------|-------|
| CNVpytor | 1.3.1 | Read-depth bins | Effective for large CNVs (>5 kb); no PoN required |
| GATK gCNV | 4.5.0.0 | Denoising + Bayesian HMM | High sensitivity for small CNVs with GIAB PoN |

**Consensus CNV:** Segments reported by both callers, OR single-caller GATK gCNV calls with quality score ≥ Q30. Final output: BED file with `CN`, `TYPE` (gain/loss), `CALLER_SUPPORT` fields.

### M4 — SMN1/SMN2

**Tool:** SMNCopyNumberCaller v1.1 (Illumina, open-source, Python)

**Mechanism:**
- Evaluates depth ratios and allele frequencies at 16 informative SNP positions in the SMN1/SMN2 paralogs (chr5:70,000,000–72,000,000, GRCh38)
- Key discriminant: `c.840C` (SMN1) vs `c.840T` (SMN2)
- Models copy number from depth + allele ratio → integer CN per gene

**Outputs:**
- `SMN1_CN` and `SMN2_CN` (integer, 0–6+)
- SMA classification: Affected (SMN1=0) / Carrier (SMN1=1) / Normal (SMN1≥2)
- **2+0 haplotype flag:** detects samples with 2 copies on one chromosome and 0 on the other — critical for accurate carrier detection as these appear as CN=2 (normal) without phasing

**Validation:** Run on GIAB trio (HG002 proband + HG003 father + HG004 mother) with known SMN1/2 truth table, and on additional GIAB samples with published SMN copy number data.

### M5 — Annotation

| Tool | Version | Annotations |
|------|---------|------------|
| AnnotSV | 3.4.2 | Gene overlap, OMIM disease, ClinVar SVs, DGV frequency, ACMG scoring, regulatory regions |
| gnomAD-SV | 2.1 | Population AF filter — flag SVs with AF > 1% in gnomAD as likely benign |

### M6 — Visualization

**Tool:** pycirclize v1.7 (Python, actively maintained)

**Circos plot ring layout (outermost → center):**

| Ring | Content | Colour |
|------|---------|--------|
| 1 | Chromosome ideograms + GRCh38 cytobands | Standard ISCN |
| 2 | CNV gains histogram | Red `#D62728` |
| 3 | CNV losses histogram | Blue `#1F77B4` |
| 4 | STR expansion markers (dot per locus) | Brown `#8C564B` |
| 5 | SMN locus highlight (chr5q13.2) with CN label | Gold `#FFBF00` |
| Center | SV links | Translocation: orange `#FF7F0E` / Inversion: purple `#9467BD` / Duplication arc: red / Tandem: green `#2CA02C` |

**Output:** `circos.svg` (embedded in HTML report) + `circos.png` (1200 dpi).  
Colour palette is colourblind-accessible (tested against deuteranopia/protanopia simulation).

### M7 — Reporting and Benchmarking

| Component | Tool | Output |
|-----------|------|--------|
| Per-sample HTML report | Jinja2 + Python | Embedded Circos SVG, CNV table, SMN section, QC summary |
| Alignment QC | MultiQC | Aggregated `multiqc_report.html` |
| GIAB benchmarking | truvari v4 | Precision, recall, F1 per SV type (DEL/INS/INV/DUP/BND); per-size-bin metrics |

**HTML report sections:**
1. Sample metadata + pipeline version
2. Alignment QC (coverage, duplication rate, insert size)
3. SV summary table (counts by type, confidence tier)
4. CNV genome plot (chromosome-level view)
5. Circos plot (embedded SVG)
6. SMN1/SMN2 copy number barchart + SMA classification badge
7. STR expansion results (flagged loci only)
8. Top annotated SVs (ACMG class 4/5, OMIM disease genes)
9. GIAB benchmark metrics (if truth set provided)

---

## 6. Repository Structure

```
SVcaller/
├── main.nf                         # Pipeline entry point, parameter parsing
├── nextflow.config                  # Docker profiles, resource limits
├── conf/
│   ├── base.config                  # CPU/mem defaults per module
│   └── docker.config                # Container image registry mappings
├── workflows/
│   ├── svcaller.nf                  # Top-level workflow orchestration
│   └── pon_build.nf                 # One-time GIAB PoN build workflow
├── subworkflows/
│   ├── preprocess.nf                # M1
│   ├── sv_calling.nf                # M2
│   ├── cnv_calling.nf               # M3
│   ├── smn_calling.nf               # M4
│   ├── annotate.nf                  # M5
│   └── report.nf                    # M6-M7
├── modules/
│   ├── bwamem2/align.nf
│   ├── picard/markduplicates.nf
│   ├── mosdepth/coverage.nf
│   ├── manta/call.nf
│   ├── delly/call.nf
│   ├── gridss/call.nf
│   ├── expansionhunter/call.nf
│   ├── jasmine/merge.nf
│   ├── cnvpytor/call.nf
│   ├── gatk/gcnv_call.nf
│   ├── gatk/gcnv_pon.nf
│   ├── smn_caller/call.nf
│   ├── annotsv/annotate.nf
│   └── truvari/bench.nf
├── bin/
│   ├── circos_plot.py               # pycirclize Circos generator
│   ├── cnv_consensus.py             # CNVpytor + gCNV merge logic
│   ├── smn_report.py                # SMN HTML section + 2+0 haplotype logic
│   └── html_report.py               # Jinja2 report builder
├── assets/
│   ├── report_template.html         # Jinja2 HTML template
│   ├── GRCh38_cytobands.txt         # UCSC cytobands for Circos ideogram
│   └── eh_catalog.json              # ExpansionHunter variant catalog
├── validation/
│   ├── giab_benchmark.sh            # truvari wrapper script for HG001–HG007
│   └── smn_truth_table.tsv          # Known SMN1/2 copy numbers for GIAB samples
├── docs/
│   └── superpowers/specs/
│       └── 2026-05-09-svcnv-caller-design.md   # This file
├── .gitignore
└── results/                         # Gitignored; per-sample output
    └── {sample_id}/
        ├── sv_merged.vcf.gz
        ├── sv_merged.vcf.gz.tbi
        ├── cnv_consensus.bed
        ├── str_expansions.vcf.gz
        ├── smn_copy_numbers.tsv
        ├── annotated_sv.tsv
        ├── circos.svg
        ├── circos.png
        └── report.html              # PRIMARY DELIVERABLE
```

---

## 7. Pipeline Parameters

```
Required:
  --input        Path to FASTQ (R1,R2 CSV samplesheet) or BAM samplesheet
  --genome       GRCh38 (only supported value in v1)
  --ref_fasta    Path to GRCh38 FASTA

Optional:
  --pon          Path to GATK gCNV PoN HDF5 (built from GIAB HG001-7)
  --outdir       Output directory (default: results/)
  --min_depth    Minimum mean depth to proceed (default: 30)
  --eh_catalog   ExpansionHunter catalog JSON (default: assets/eh_catalog.json)
  --giab_truth   Path to GIAB SV truth VCF for truvari benchmarking
  --max_cpus     Max CPUs per process (default: 64)
  --max_memory   Max memory (default: 120.GB)
```

---

## 8. Validation Plan

### 8.1 GIAB Samples

| Sample | Population | Special value |
|--------|-----------|---------------|
| HG001 (NA12878) | CEU | Most-benchmarked sample; extensive SV truth |
| HG002 | AJ (Ashkenazi) | Best SV truth set (GIAB v0.6); primary benchmark |
| HG003 | AJ | HG002's father; trio validation |
| HG004 | AJ | HG002's mother; trio validation |
| HG005 | Han Chinese | Ethnic diversity |
| HG006 | Han Chinese | HG005's father |
| HG007 | Han Chinese | HG005's mother |

All 7 samples processed independently. HG002 is the primary benchmark sample (GIAB SV v0.6 truth set).

### 8.2 SV Benchmarking (truvari)

Metrics reported per SV type and per size bin (50–300 bp, 300 bp–1 kb, 1–10 kb, >10 kb):
- Precision, Recall, F1
- TP, FP, FN counts

### 8.3 SMN Validation

| Sample | Known SMN1 CN | Known SMN2 CN |
|--------|--------------|--------------|
| HG002  | 2            | 1            |
| HG003  | 2            | 2            |
| HG004  | 2            | 1            |
| (additional GIAB samples TBD from published data) | | |

Trio with known SMN1/2 CNVs (user-provided) will be added to `validation/smn_truth_table.tsv`.

---

## 9. Nanopore Extension (v2 Architecture)

The pipeline is designed so long-read support requires **no changes to M3–M7**. Extension points:

- `subworkflows/preprocess_ont.nf` — dorado basecalling + minimap2 alignment
- `subworkflows/sv_calling_ont.nf` — Sniffles2 v2 + CuteSV v2 (replace Manta/DELLY/GRIDSS)
- `main.nf` branches on `--sequencing_type [illumina|ont]`
- JASMINE merge, AnnotSV, pycirclize, and HTML report are shared unchanged

Long-read CNV (e.g., HiFiCNV or Spectre) plugs into M3 as an alternative subworkflow.

---

## 10. Docker Container Strategy

All tools run in pre-built biocontainer Docker images. Custom images are built only for:
- `SMNCopyNumberCaller` (no official biocontainer) — Dockerfile in `modules/smn_caller/`
- `bin/` Python scripts (pycirclize, Jinja2, pandas) — single `svcaller-utils` image

Container images pinned by SHA256 digest in `conf/docker.config` for full reproducibility.

---

## 11. Resource Profile (default, 64 CPU / 125 GB RAM server)

| Module | CPUs | Memory | Expected runtime (30x WGS) |
|--------|------|--------|---------------------------|
| M1 (BWA-MEM2) | 32 | 16 GB | ~45 min |
| M2 Manta | 16 | 8 GB | ~20 min |
| M2 DELLY | 4 | 8 GB | ~60 min |
| M2 GRIDSS | 8 | 32 GB | ~90 min |
| M2 ExpansionHunter | 4 | 4 GB | ~15 min |
| M3 CNVpytor | 4 | 8 GB | ~20 min |
| M3 GATK gCNV | 8 | 16 GB | ~30 min |
| M4 SMNCopyNumberCaller | 2 | 4 GB | ~10 min |
| M5–M7 | 4 | 8 GB | ~15 min |
| **Total wall-clock (parallelised)** | — | — | **~2.5 hours** |

---

## 12. Technology Summary

| Layer | Technology |
|-------|-----------|
| Workflow | Nextflow DSL2 v25.10.4 |
| Containers | Docker v28.3.2 |
| Alignment | BWA-MEM2 v2.2.1 |
| SV callers | Manta 1.6, DELLY 1.2.6, GRIDSS 2.13.2, ExpansionHunter 5.0 |
| SV merge | JASMINE 1.1.5 |
| CNV callers | CNVpytor 1.3.1, GATK gCNV 4.5 |
| SMN | SMNCopyNumberCaller 1.1 |
| Annotation | AnnotSV 3.4.2, gnomAD-SV 2.1 |
| Visualization | pycirclize 1.7, matplotlib |
| Reporting | Jinja2, MultiQC |
| Benchmarking | truvari 4.x |
| Reference | GRCh38/hg38 (Ensembl release 112) |
