# SVcaller

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A525.10.4-23aa62.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/docker-enabled-2496ED.svg?logo=docker)](https://www.docker.com/)
[![Python](https://img.shields.io/badge/python-3.11-3776AB.svg?logo=python)](https://www.python.org/)
[![Tests](https://img.shields.io/badge/tests-25%20passing-brightgreen.svg)](#testing)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**A production-ready Nextflow DSL2 pipeline for calling structural variants (SVs), copy number variants (CNVs), and SMN1/SMN2 copy numbers from Illumina paired-end WGS data (PE150, ≥30×).**

Accepts FASTQ or pre-aligned BAM input. Produces a per-sample HTML report with an embedded Circos plot, annotated SV/CNV tables, SMN1/SMN2 classification, and optional GIAB benchmark metrics.

---

## Features

- **Ensemble SV calling** — Manta + DELLY + GRIDSS merged with JASMINE (min 2/3 callers)
- **STR genotyping** — ExpansionHunter with full Illumina repeat catalog
- **Dual CNV calling** — CNVpytor + GATK gCNV with consensus merging
- **SMN1/SMN2 copy number** — SMNCopyNumberCaller with 2+0 haplotype detection
- **Clinical annotation** — AnnotSV 3.4 with gnomAD-SV allele frequency filtering
- **Genome-wide visualization** — pycirclize Circos plot (SVs, CNV gains/losses, SMN locus)
- **Self-contained HTML report** — Bootstrap 5, embedded SVG, optional truvari GIAB benchmark
- **Reproducible** — every tool pinned in a Docker container; no conda required

---

## Pipeline Overview

```
FASTQ or BAM
     │
     ▼
┌─────────────────────────────────────────────────────────────────┐
│  M1 · Pre-processing                                            │
│  BWA-MEM2 align → samtools sort → Picard MarkDup → mosdepth QC │
└────────────────────────────┬────────────────────────────────────┘
                             │ BAM
          ┌──────────────────┼──────────────────┐
          ▼                  ▼                  ▼
┌─────────────────┐ ┌───────────────┐ ┌──────────────────┐
│ M2 · SV Calling │ │ M3 · CNV Call │ │ M4 · SMN Calling │
│ Manta + DELLY   │ │ CNVpytor      │ │ SMNCopyNumber    │
│ GRIDSS + EH     │ │ GATK gCNV     │ │ Caller v1.1      │
│ JASMINE merge   │ │ consensus BED │ │                  │
└────────┬────────┘ └───────┬───────┘ └───────┬──────────┘
         │                  │                 │
         ▼                  │                 │
┌────────────────┐          │                 │
│ M5 · Annotate  │          │                 │
│ AnnotSV 3.4    │          │                 │
│ gnomAD-SV AF   │          │                 │
│ filter <1%     │          │                 │
└────────┬───────┘          │                 │
         └──────────────────┴─────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────┐
│  M6/M7 · Visualization & Reporting          │
│  pycirclize Circos plot + Jinja2 HTML report │
│  (optional) truvari GIAB benchmark          │
└─────────────────────────────────────────────┘
```

---

## Requirements

| Dependency | Version | Notes |
|---|---|---|
| [Nextflow](https://nextflow.io) | ≥ 25.10.4 | `curl -s get.nextflow.io \| bash` |
| Docker | ≥ 24 | All tools run in containers |
| Python | ≥ 3.11 | For local test/report scripts |
| Java | ≥ 17 | Required by Nextflow |

**System:** 16+ CPU cores, 64+ GB RAM, 500 GB disk recommended for 30× WGS.

---

## Quick Start

```bash
# 1. Clone
git clone https://github.com/alvin8-git/SVcaller.git
cd SVcaller

# 2. Install Nextflow (if not already installed)
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

# 3. Download reference data (~60 min first time)
bash validation/download_refs.sh

# 4. Run on a sample
nextflow run main.nf \
    -profile docker \
    --input /path/to/samplesheet.csv \
    --ref_fasta /data/alvin/ref/GRCh38/GRCh38.fasta \
    --outdir results
```

---

## Input Samplesheet

A CSV file with one row per sample. Provide either `fastq_1`/`fastq_2` **or** `bam` — not both.

```csv
sample,fastq_1,fastq_2,bam
HG002,/data/HG002_R1.fastq.gz,/data/HG002_R2.fastq.gz,
HG003,,,/data/HG003.GRCh38.bam
```

| Column | Required | Description |
|--------|----------|-------------|
| `sample` | Yes | Unique sample ID (no spaces) |
| `fastq_1` | FASTQ mode | Absolute path to R1 FASTQ.gz |
| `fastq_2` | FASTQ mode | Absolute path to R2 FASTQ.gz |
| `bam` | BAM mode | Absolute path to sorted BAM |

---

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Path to samplesheet CSV |
| `--ref_fasta` | required | GRCh38 reference FASTA |
| `--outdir` | `results` | Output directory |
| `--min_depth` | `30` | Minimum mean coverage (fails pipeline if below) |
| `--pon` | null | GATK gCNV Panel of Normals HDF5 (see [PoN Build](#panel-of-normals)) |
| `--intervals` | null | Preprocessed intervals BED for GATK gCNV |
| `--annotsv_db` | null | AnnotSV annotation directory |
| `--eh_catalog` | `assets/eh_catalog.json` | ExpansionHunter variant catalog |
| `--giab_truth` | null | GIAB SV truth VCF.gz for truvari benchmarking |
| `--max_cpus` | `64` | Max CPUs per process |
| `--max_memory` | `120.GB` | Max memory per process |

---

## Running the Pipeline

### FASTQ input

```bash
nextflow run main.nf \
    -profile docker \
    --input samplesheet.csv \
    --ref_fasta /data/alvin/ref/GRCh38/GRCh38.fasta \
    --pon /data/alvin/SVcaller/pon/giab_cnv_pon.hdf5 \
    --annotsv_db /data/alvin/ref/annotsv \
    --outdir results \
    --max_cpus 32 \
    --max_memory '64.GB'
```

### BAM input (skip alignment)

```bash
nextflow run main.nf \
    -profile docker \
    --input bam_samplesheet.csv \
    --ref_fasta /data/alvin/ref/GRCh38/GRCh38.fasta \
    --outdir results
```

### With GIAB benchmarking

```bash
nextflow run main.nf \
    -profile docker \
    --input samplesheet.csv \
    --ref_fasta /data/alvin/ref/GRCh38/GRCh38.fasta \
    --giab_truth /data/alvin/ref/GIAB/HG002_SV_v0.6.vcf.gz \
    --outdir results
```

### Resume a failed run

```bash
nextflow run main.nf -profile docker --input samplesheet.csv \
    --ref_fasta /data/alvin/ref/GRCh38/GRCh38.fasta -resume
```

---

## Panel of Normals

GATK gCNV requires a Panel of Normals (PoN) built from ≥10 normal samples. Build once from GIAB HG001–HG007:

```bash
# 1. Edit validation/giab_samplesheet.csv with your BAM paths
# 2. Build PoN
nextflow run workflows/pon_build.nf \
    -profile docker \
    --input validation/giab_samplesheet.csv \
    --ref_fasta /data/alvin/ref/GRCh38/GRCh38.fasta \
    --outdir /data/alvin/SVcaller/pon

# PoN output: /data/alvin/SVcaller/pon/giab_cnv_pon.hdf5
```

Pass it to the main pipeline with `--pon /data/alvin/SVcaller/pon/giab_cnv_pon.hdf5`.

---

## Output

```
results/
└── <sample_id>/
    ├── <sample_id>.report.html          # Self-contained HTML report
    ├── <sample_id>.sv_merged.vcf.gz     # Ensemble SV calls (Manta+DELLY+GRIDSS)
    ├── <sample_id>.sv_merged.vcf.gz.tbi
    ├── <sample_id>.str.vcf.gz           # STR calls (ExpansionHunter)
    ├── <sample_id>.cnv_consensus.bed    # Consensus CNV calls
    ├── <sample_id>.smn.tsv              # SMN1/SMN2 copy numbers
    ├── <sample_id>.annotated.tsv        # AnnotSV-annotated SVs
    ├── <sample_id>.circos.svg           # Genome-wide Circos plot
    ├── <sample_id>.circos.png
    └── <sample_id>.truvari/             # GIAB benchmark (if --giab_truth set)
        └── summary.json
```

The HTML report includes:
- Alignment QC (coverage, duplicate rate)
- SV summary table by type
- Embedded Circos plot
- SMN1/SMN2 copy number with SMA classification
- Top annotated SVs (ACMG pathogenicity score ≥ 0.9)
- GIAB benchmark precision/recall/F1 (optional)

---

## GIAB Benchmarking

Validate SV calls against GIAB truth sets:

```bash
bash validation/giab_benchmark.sh \
    HG002 \
    results/HG002/HG002.sv_merged.vcf.gz
```

**Current benchmark on HG002 (4 callers — Manta + Delly + GRIDSS + Scramble, run69):**

| Metric | Result |
|--------|--------|
| Precision | 0.648 |
| Recall | 0.231 |
| F1 | **0.341** |

Recall is limited by MEI insertions (dominated by L1/ALU not caught by Scramble at default settings). See [Design decisions](docs/explanation-design.md#scramble-mei-canonical-svlen-estimates) for context.

---

## Documentation

| Document | What it covers |
|----------|---------------|
| [Getting started tutorial](docs/tutorial-getting-started.md) | First run from scratch to open report |
| [How to run the GIAB validation](docs/howto-run-validation.md) | Benchmarking against GIAB HG002 truth set |
| [How to build a Panel of Normals](docs/howto-build-pon.md) | GATK gCNV PON construction |
| [How to interpret the HTML report](docs/howto-interpret-report.md) | Clinical interpretation of all report sections and Circos rings |
| [Parameter reference](docs/reference-parameters.md) | All CLI flags, samplesheet format, output files, resource labels |
| [Architecture reference](docs/reference-architecture.md) | Module-by-module technical description with I/O and design notes |
| [Design decisions](docs/explanation-design.md) | Why the pipeline is built the way it is (channel patterns, sentinel files, PON choices) |

---

## Testing

```bash
# Install test dependencies
pip install pytest pycirclize matplotlib jinja2

# Run all 25 unit tests
pytest tests/ -v
```

Tests cover: samplesheet validation, CNV consensus merging, SMN classification, Circos plot parsing, HTML report rendering.

---

## Docker Containers

All tools run in pinned Docker containers. No local tool installation required beyond Docker itself.

| Module | Container |
|--------|-----------|
| BWA-MEM2 | `quay.io/biocontainers/bwa-mem2:2.2.1` |
| Picard MarkDup | `quay.io/biocontainers/picard:3.2.0` |
| mosdepth | `quay.io/biocontainers/mosdepth:0.3.8` |
| Manta | `quay.io/biocontainers/manta:1.6.0` |
| DELLY | `quay.io/biocontainers/delly:1.2.6` |
| GRIDSS | `gridss/gridss:2.13.2` |
| ExpansionHunter | `quay.io/biocontainers/expansionhunter:5.0.0` |
| JASMINE | `quay.io/biocontainers/jasminesv:1.1.5` |
| CNVpytor | `quay.io/biocontainers/cnvpytor:1.3.1` |
| GATK | `broadinstitute/gatk:4.5.0.0` |
| SMNCopyNumberCaller | `svcaller/smncopynum:1.1` |
| AnnotSV | `quay.io/biocontainers/annotsv:3.4.2` |
| truvari | `quay.io/biocontainers/truvari:4.2.2` |

---

## Genome Build

GRCh38 / hg38 only. `chr`-prefixed chromosome names required.

---

## License

MIT — see [LICENSE](LICENSE).

---

## Citation

If you use SVcaller in your research, please cite the individual tools (Manta, DELLY, GRIDSS, CNVpytor, GATK, SMNCopyNumberCaller, AnnotSV, ExpansionHunter) as appropriate for your analysis.
