# SVcaller Parameter Reference

Complete reference for all CLI parameters, samplesheet columns, and output files.

## Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `--input` | Path | Samplesheet CSV. See [Samplesheet format](#samplesheet-format). |
| `--ref_fasta` | Path | GRCh38 reference FASTA. Index files (`.fai`, `.0123`, `.bwt.2bit.64`) are inferred from the same path prefix. For FASTQ inputs use `hg38.canonical.fa` (chr1-22+X+Y+M only) to skip FILTER_CHROMS. |
| `--eh_catalog` | Path | ExpansionHunter locus catalog JSON. Default catalog at `assets/eh_catalog.json`. |

## Optional Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--pon` | Path | null | GATK gCNV Panel of Normals HDF5. Required for CNV calling via GATK. Without it, only CNVpytor runs. |
| `--intervals` | Path | null | Target capture BED (WES) or autosomal BED (WGS). Required for GATK CNV. Use `/data/alvin/ref/GRCh38/wgs_autosomal.bed` for WGS. |
| `--annotsv_db` | Path | null | AnnotSV annotations directory. Point to the parent of `Annotations_Human/`. Without it, ANNOTSV emits an empty stub TSV and annotation sections in the report are blank. |
| `--giab_truth` | Path | null | GIAB truth VCF.gz (enables Truvari benchmarking in the report). |
| `--min_depth` | Integer | 30 | Minimum mean coverage (mosdepth). Pipeline halts if the sample is below this threshold. |
| `--outdir` | Path | `results` | Output directory. Per-sample subdirectories are created automatically. |
| `--utils_container` | String | `svcaller/utils:1.1` | Docker image for Python bin/ scripts. |
| `--skip_gridss` | Boolean | false | Skip GRIDSS (saves 4-6 h and 60 GB RAM). Manta + Delly + Scramble still run. Use for SMN validation runs or resource-constrained environments. |
| `--tiered_gridss` | Boolean | false | Run GRIDSS only on Manta residual regions (smaller input → ~2.5 h, 40 GB). Reduces sensitivity for complex rearrangements in non-Manta regions. |

## Samplesheet Format

CSV with header row. Each row is one sample. Provide either a FASTQ pair **or** a BAM — not both.

```
sample,fastq_1,fastq_2,bam
HG002,/path/HG002_R1.fq.gz,/path/HG002_R2.fq.gz,
HG003,,,/path/HG003.bam
```

| Column | Required | Description |
|--------|----------|-------------|
| `sample` | Always | Sample identifier. Used as prefix for all output files and subdirectory name. |
| `fastq_1` | If no BAM | Path to R1 FASTQ (gzipped). Pipeline aligns with BWA-MEM2. |
| `fastq_2` | If no BAM | Path to R2 FASTQ (gzipped). |
| `bam` | If no FASTQ | Path to pre-aligned, sorted, duplicate-marked BAM. Must have a `.bai` index at the same path. |

**FASTQ vs BAM behaviour:**
- FASTQ samples aligned to `hg38.canonical.fa` skip the FILTER_CHROMS step (saves ~25 min/sample). The canonical reference contains only chr1-22+X+Y+M.
- BAM samples always run FILTER_CHROMS regardless of reference, because pre-aligned BAMs may contain ALT/decoy contigs.

## Output Files

All outputs land in `{outdir}/{sample}/`:

| File | Description |
|------|-------------|
| `{sample}.report.html` | Per-sample clinical HTML report. Contains QC metrics, SV summary, Circos plot, top annotated SVs, SMN copy number, STR loci, and optionally GIAB benchmark results. |
| `{sample}.sv_merged.vcf.gz` | Merged SV calls from Manta + Delly + GRIDSS + Scramble (Jasmine merge). Tabix-indexed. |
| `{sample}.sv_merged.vcf.gz.tbi` | Tabix index for the merged VCF. |
| `{sample}.filtered.tsv` | AnnotSV-annotated SV table after gnomAD SV frequency filter. Tab-separated, with ACMG classification and gene annotations. |
| `{sample}.cnv_consensus.bed` | Consensus CNV calls from CNVpytor + GATK gCNV. Columns: chrom, start, end, cn, svtype, caller_support, confidence, quality, sample. |
| `{sample}.smn.tsv` | SMN1/SMN2 copy number table from SMNCopyNumberCaller. |
| `{sample}.circos.png` | Genome-wide Circos plot (PNG, 150 DPI). |

MultiQC HTML is written to `{outdir}/multiqc_report.html`.

## Resource Labels

Processes are labelled by resource tier. Retry on OOM (exit codes 137, 143, 104, 134, 139) doubles memory automatically.

| Label | CPUs | Memory | Typical processes |
|-------|------|--------|-------------------|
| `process_single` | 1 | 2 GB | Report assembly, samplesheet parsing |
| `process_low` | 2 | 4 GB | mosdepth, Scramble, FastQC |
| `process_medium` | 4 | 8 GB | Manta, CNVpytor, SMNCopyNumberCaller, AnnotSV |
| `process_high` | 8 | 16 GB | BWA-MEM2 alignment, Delly, GATK gCNV |
| `process_gridss` | 20 | 60 GB | GRIDSS (full BAM mode) |

## Related

- [Architecture reference](reference-architecture.md) — module-by-module technical description
- [How to run the GIAB validation](howto-run-validation.md)
- [How to build a Panel of Normals](howto-build-pon.md)
