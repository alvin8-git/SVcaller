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
| `--intervals` | Path | null | Target capture BED (WES) or autosomal BED (WGS). Required for GATK CNV. Use `/path/to/ref/GRCh38/wgs_autosomal.bed` for WGS. |
| `--annotsv_db` | Path | null | AnnotSV annotations directory. Point to the parent of `Annotations_Human/`. Without it, ANNOTSV emits an empty stub TSV and annotation sections in the report are blank. |
| `--giab_truth` | Path | null | GIAB truth VCF.gz (enables Truvari benchmarking in the report). |
| `--min_depth` | Integer | 25 | Minimum mean coverage (mosdepth). Pipeline halts if the sample is below this threshold. |
| `--outdir` | Path | `results` | Output directory. Per-sample subdirectories are created automatically. |
| `--tmp_dir` | Path | `$TMPDIR` or `/tmp` | Host temp directory bind-mounted into Docker containers as `/tmp` and exported as `TMPDIR`. Defaults to the invoking shell's `$TMPDIR`, falling back to `/tmp`, so no host path is baked in. Override on machines with a dedicated scratch volume (e.g. `--tmp_dir /data/scratch/tmp`). |
| `--utils_container` | String | `svcaller/utils:1.3` | Docker image for Python bin/ scripts. |
| `--auto_cleanup` | Boolean | false | Delete the run's `-work-dir` automatically on **successful** completion (via the `workflow.onComplete` hook). Off by default so intermediates survive for `-resume`. When false, the run prints a tip to run `bash bin/nf-cleanup.sh <sampleId>` instead. See [Storage & Cache Management](#storage--cache-management). |
| `--max_cpus` | Integer | 64 | Upper bound on CPUs any single process may request (caps the `task.attempt` scaling). Lower it on small machines to prevent over-subscription. |
| `--max_memory` | String | `120.GB` | Upper bound on memory any single process may request. Lower it on small machines so retries don't request more RAM than exists. |
| `--max_time` | String | `240.h` | Upper bound on wall-clock time any single process may request. |
| `--skip_gridss` | Boolean | false | Skip GRIDSS (saves 4-6 h and 60 GB RAM). Manta + Delly + Scramble still run. Use for SMN validation runs or resource-constrained environments. |
| `--skip_melt` | Boolean | false | Skip MELT MEI calling (saves ~2 h). Use when the `svcaller/melt:2.2.2` container is unavailable or when mobile element insertions are not clinically relevant. |
| `--tiered_gridss` | Boolean | false | Run GRIDSS only on Manta residual regions (smaller input → ~2.5 h, 40 GB). Reduces sensitivity for complex rearrangements in non-Manta regions. |
| `--skip_svaba` | Boolean | true | Off by default. SvABA is staged but not part of the ensemble: `JASMINE_MERGE` lists only `vcfs[0..4]` and never references SvABA's VCF, so even with `--skip_svaba false` its calls do not reach the merge (only the runtime changes). Re-enabling means wiring `vcfs[5]` into `vcf_list.txt` first. Opting in also requires the classic BWA index (`--bwa_index`); `main.nf` fails loud if it is absent. |
| `--melt_min_reads` | Integer | 3 | Passed as MELT's `-sr` flag: minimum split-read support to keep a call. MELT's own `-sr` default is `0` (no filtering), so `3` is **more** stringent than stock MELT. Set `0` to maximise recall. |
| `--sv_pon` | Path | null | GIAB 7-sample SV Panel of Normals BED. SVs matching population variants in the PON are flagged as `COMMON_SV` in the report. Use `pon/sv_pon/giab_sv_pon.bed`. |
| `--giab_truth_v5q` | Path | null | GIAB v5.0q truth VCF.gz. Enables a second Truvari benchmark pass alongside `--giab_truth` (T2TQ100-V1.0). |
| `--melt_refs` | Path | null | Path to MELT `me_refs/` directory. Auto-detected from the container at `/opt/melt/me_refs` if unset. Only needed when running MELT outside Docker. |

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
- FASTQ samples aligned to `hg38.canonical.fa` skip the FILTER_CHROMS step (saves ~70 min/sample on a ~30× BAM — measured; see CLAUDE.md "FILTER_CHROMS is slow"). The canonical reference contains only chr1-22+X+Y+M.
- BAM samples always run FILTER_CHROMS regardless of reference, because pre-aligned BAMs may contain ALT/decoy contigs.

## Output Files

All outputs land in `{outdir}/{sample}/`:

| File | Description |
|------|-------------|
| `{sample}.report.html` | Per-sample clinical HTML report. Contains QC metrics, SV summary, Circos plot, top annotated SVs, SMN copy number, STR loci, and optionally GIAB benchmark results. |
| `{sample}.sv_merged.vcf.gz` | Merged SV calls from the 5-caller ensemble: Manta + Delly + GRIDSS (core) plus Scramble + MELT (added only when each has calls), via Jasmine merge. Tabix-indexed. |
| `{sample}.sv_merged.vcf.gz.tbi` | Tabix index for the merged VCF. |
| `{sample}.filtered.tsv` | AnnotSV-annotated SV table after gnomAD SV frequency filter. Tab-separated, with ACMG classification and gene annotations. |
| `{sample}.cnv_consensus.bed` | Consensus CNV calls from CNVpytor + GATK gCNV. Columns: chrom, start, end, cn, svtype, caller_support, confidence, quality, sample. |
| `{sample}.smn.tsv` | SMN1/SMN2 copy number table from SMNCopyNumberCaller. |
| `{sample}.circos.svg` | Genome-wide Circos plot (SVG, embedded inline in the HTML report). |
| `{sample}.circos.png` | Genome-wide Circos plot (PNG, 150 DPI, fallback if SVG rendering fails). |
| `{sample}.variants.xlsx` | Excel workbook with four sheets: SVs, CNVs, STRs, SMN. Suitable for clinical reporting workflows. |

QC metrics (mosdepth summary, flagstat, insert-size) publish to `{outdir}/{sample}/qc/`, so the report's QC section is reproducible after a work-dir cleanup.

MultiQC HTML is written to `{outdir}/multiqc_report.html`.

## Resource Labels

Processes are labelled by resource tier in `conf/base.config`. Values shown are the **first-attempt** request; on a retry triggered by an OOM exit code (137, 143, 104, 134, 139) the memory (and most CPU) requests scale by `task.attempt`. Every request is clamped by `--max_cpus` / `--max_memory` / `--max_time`.

| Label | CPUs | Memory (attempt 1) | Per-retry scaling | Typical processes |
|-------|------|--------------------|-------------------|-------------------|
| `process_single` | 1 | 6 GB | memory ×attempt | Report assembly, samplesheet parsing, CNV consensus, Circos |
| `process_low` | 4 | 12 GB | cpus + memory ×attempt | mosdepth, Scramble, FastQC |
| `process_medium` | 8 | 36 GB | cpus + memory ×attempt | CNVpytor, SMNCopyNumberCaller, AnnotSV |
| `process_high` | 16 | 32 GB | cpus + memory ×attempt | GATK gCNV |
| `process_gridss` | 16 | 32 GB | memory only (32→64→96 GB) | GRIDSS (full BAM mode) |

Several processes pin resources by name (overriding the tier): `BWAMEM2_ALIGN` (16 CPU / 24 GB), `MANTA_CALL` (16 CPU / 32 GB), `PICARD_MARKDUP` (4 CPU / 16 GB), `GATK_COLLECT_COUNTS` (4 CPU / 16 GB), `SAMTOOLS_FILTER_CHROMS` (8 CPU / 8 GB), `DELLY_CALL_SVTYPE` (2 CPU — carries the `process_medium` label for memory/time but `delly call` is effectively single-threaded and parallelises across 5 SV-type tasks, so an 8-CPU reservation only oversubscribes the scheduler and serialises the fan-out; the 2-CPU pin is a reservation-only change and does not alter Delly's output). `GRIDSS_CALL` uses `maxRetries = 3` (4 attempts: 32 → 64 → 96 GB) because it is the most likely caller to silently OOM. `SVABA_CALL` fixes CPUs at 16 so retries don't claim extra slots (SvABA is memory-bound, not CPU-bound).

For a low-overhead workstation profile that caps these tiers so the suite runs without resource exhaustion, see [Local vs cluster profiles](#local-vs-cluster-profiles) below.

## Storage & Cache Management

Long pipeline runs accumulate large intermediates in the `-work-dir`. SVcaller hardens storage three ways:

**1. Per-sample work directories.** Always pass one `-work-dir` per sample or batch — `work_<sampleId>` for single samples, `work_<batchName>` for batches. Never share a `work/` across run types: it causes Nextflow session-lock conflicts, blocks per-sample cleanup, and accumulates unreclaimable intermediates.

**2. Automatic cleanup on success.** Pass `--auto_cleanup true` to delete the `-work-dir` when (and only when) the run completes successfully. This trades away `-resume` for the next run, so use it only for one-shot production runs you won't re-enter.

**3. Post-run cleanup script.** When `--auto_cleanup` is off (the default), run `bash bin/nf-cleanup.sh <sampleId>` after results are published. It verifies the sample's outputs exist under `--outdir`, removes `work_<sampleId>`, and prunes orphaned `.nextflow/cache` sessions (skipping any still holding a lock).

**4. `storeDir` caches survive `nextflow clean`.** Three processes write to a persistent `storeDir` so re-runs and new samples skip expensive recomputation. These caches are keyed on inputs (reference + intervals) and are safe to share across runs against the same reference — do **not** delete them between samples:

| Process | storeDir path | What it caches |
|---------|---------------|----------------|
| `SAMTOOLS_FILTER_CHROMS` | `${outdir}/.cache/filter_chroms` | Chrom-filtered BAMs (BAM inputs) |
| `GRIDSS_SETUP` | `${outdir}/cache/gridss_ref` | GRIDSS reference setup |
| `GATK_PREPROCESS_INTERVALS` | `${outdir}/cache/gatk_preprocess` | Binned interval_list (stable per ref+BED) |

## Environment Variables

For multi-user compute and to keep temp/container state off small root partitions, export these before launching:

| Variable | Purpose |
|----------|---------|
| `NXF_ANSI_LOG=false` | **Required for background/nohup runs.** Without it Nextflow's ANSI renderer deadlocks all JVM threads when there is no TTY. |
| `TMPDIR` | Redirect temp files off the (small) root partition. The pipeline already sets `TMPDIR=/path/to/tmp` in `nextflow.config`; override it for your site. |
| `NXF_HOME` | Per-user Nextflow home (plugins, assets). Set a distinct path per user on shared machines to avoid plugin-cache contention. |
| `NXF_SINGULARITY_CACHEDIR` | Per-user image cache for Singularity/Apptainer deployments. Set distinct paths per user so concurrent runs don't collide while building/pulling the same image. |
| `-cache <dir>` (CLI flag) | Give each concurrent sample run its own `.nextflow/cache` to prevent session-lock conflicts when running several samples in parallel from the same checkout. |

## Local vs cluster profiles

Select a profile with `-profile <name>` (combinable, comma-separated):

| Profile | File | Use when |
|---------|------|----------|
| `docker` | `conf/docker.config` | **Default for all runs.** Every process runs in its verified quay.io biocontainer. |
| `local` | `conf/local.config` | Conda-based execution, no Docker. Points processes at the `svcaller` conda env and puts `bin/` on PATH. Use on a workstation without a container runtime. |
| `test` | `conf/test.config` | Minimal smoke-test resources. |

All profiles inherit the tiers in `conf/base.config`. On a resource-constrained workstation, cap the tiers without editing `base.config` by lowering the bounds at launch, e.g. `--max_cpus 8 --max_memory 32.GB`; combine with `--skip_gridss true` (saves ~60 GB RAM peak) so a single sample fits on a laptop-class machine. For an enterprise cluster, add an executor block (`executor.name = 'slurm'`, queue, `scratch = true`) in a site config and pass it with `-c site.config`.

## Related

- [Architecture reference](reference-architecture.md) — module-by-module technical description
- [How to run the GIAB validation](howto-run-validation.md)
- [How to build a Panel of Normals](howto-build-pon.md)
