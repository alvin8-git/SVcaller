# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Skill routing

When the user's request matches an available skill, invoke it via the Skill tool. When in doubt, invoke the skill.

Key routing rules:
- Product ideas/brainstorming → invoke /office-hours
- Strategy/scope → invoke /plan-ceo-review
- Architecture → invoke /plan-eng-review
- Design system/plan review → invoke /design-consultation or /plan-design-review
- Full review pipeline → invoke /autoplan
- Bugs/errors → invoke /investigate
- QA/testing site behavior → invoke /qa or /qa-only
- Code review/diff check → invoke /review
- Visual polish → invoke /design-review
- Ship/deploy/PR → invoke /ship or /land-and-deploy
- Save progress → invoke /context-save
- Resume context → invoke /context-restore
- Author a backlog-ready spec/issue → invoke /spec

## What This Pipeline Does

SVcaller is a Nextflow DSL2 pipeline for human WGS structural variant (SV), copy-number variant (CNV), short tandem repeat (STR), and SMN1/SMN2 copy-number calling. It targets GRCh38 at ≥30× coverage and produces per-sample HTML reports.

## Running the Pipeline

**Use `-profile docker` for all runs.** All quay.io biocontainer tags have been verified and fixed.

**Always use `NXF_ANSI_LOG=false` for background runs.** Without it, Nextflow's ANSI renderer deadlocks all JVM threads when there's no TTY (nohup/background). After launching, verify cached/submitted task lines appear in the log within 3 minutes.

**Work directory convention: one `-work-dir` per sample or batch run.** Never reuse a shared `work/` across different run types — it causes session lock conflicts, prevents per-sample cleanup, and accumulates 3+ TB of unreclaimable intermediates. Use `work_<sampleId>` for single-sample runs and `work_<batchName>` for multi-sample batches. Clean a sample's work dir once its results are published: `rm -rf work_<sampleId>`.

**Path variables — set once per environment.** The commands below use `${REF}`, `${PROJ}`, `${TMP}` so they stay portable and copy-paste runnable:
```bash
export REF=/path/to/ref          # reference-data root: GRCh38, GIAB, AnnotSV
export PROJ="$(pwd)"             # this checkout (absolute path)
export TMP=/path/to/tmp          # scratch + logs; keep off the small root partition
```

```bash
# SV/CNV validation — HG002 only; Truvari benchmark against GIAB SV truth
# Use hg38.canonical.fa (chr1-22+X+Y+M only) for FASTQ inputs — skips 25-min FILTER_CHROMS step
NXF_ANSI_LOG=false nohup nextflow run main.nf -profile docker \
  --input validation/validation_samplesheet.csv \
  --ref_fasta ${REF}/GRCh38/hg38.canonical.fa \
  --intervals ${REF}/GRCh38/wgs_autosomal.bed \
  --pon ${PROJ}/pon/pon/giab_cnv_pon.hdf5 \
  --giab_truth ${REF}/GIAB/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz \
  --giab_truth_v5q ${REF}/GIAB/HG002_GRCh38_v5.0q_stvar.vcf.gz \
  --eh_catalog assets/eh_catalog.json \
  --annotsv_db ${REF}/annotsv/Annotations_Human \
  --sv_pon ${PROJ}/pon/sv_pon/giab_sv_pon.bed \
  --outdir ${PROJ}/results \
  -work-dir ${PROJ}/work_HG002 \
  -resume > ${TMP}/main_runN.log 2>&1 &

# SMN validation — SMA trio only; no Truvari (no SV truth for clinical samples)
# --skip_gridss saves 4-6 h; only SMN1/2 CN and CNV results matter here
# Must use hg38.canonical.fa: FILTER_CHROMS strips alt contigs from BAM; if ref still has them Manta fails
NXF_ANSI_LOG=false nohup nextflow run main.nf -profile docker \
  --input validation/smn_validation_samplesheet.csv \
  --ref_fasta ${REF}/GRCh38/hg38.canonical.fa \
  --intervals ${REF}/GRCh38/wgs_autosomal.bed \
  --pon ${PROJ}/pon/pon/giab_cnv_pon.hdf5 \
  --eh_catalog assets/eh_catalog.json \
  --sv_pon ${PROJ}/pon/sv_pon/giab_sv_pon.bed \
  --skip_gridss true \
  --outdir ${PROJ}/results_smn \
  -work-dir ${PROJ}/work_smn \
  -resume > ${TMP}/smn_runN.log 2>&1 &

# GIAB PON sample reports — HG001, HG003-HG007 (HG002 done via validation run)
# --skip_melt: MELT adds 12 h with no SV truth benefit for these report-only samples
# SV calling not cached from pon_build.nf — expect 8-12 h for 6 samples
NXF_ANSI_LOG=false nohup nextflow run main.nf -profile docker \
  --input validation/giab_reports_samplesheet.csv \
  --ref_fasta ${REF}/GRCh38/hg38.canonical.fa \
  --intervals ${REF}/GRCh38/wgs_autosomal.bed \
  --pon ${PROJ}/pon/pon/giab_cnv_pon.hdf5 \
  --eh_catalog assets/eh_catalog.json \
  --annotsv_db ${REF}/annotsv/Annotations_Human \
  --sv_pon ${PROJ}/pon/sv_pon/giab_sv_pon.bed \
  --skip_gridss true \
  --skip_melt true \
  --outdir ${PROJ}/results_giab \
  -work-dir ${PROJ}/work_giab_reports \
  -resume > ${TMP}/giab_reports_runN.log 2>&1 &

# Clinical sample (single sample, new patient)
# Replace SAMPLEID with the actual sample name — never reuse another sample's work dir
NXF_ANSI_LOG=false nohup nextflow run main.nf -profile docker \
  --input /path/to/SAMPLEID_samplesheet.csv \
  --ref_fasta ${REF}/GRCh38/hg38.canonical.fa \
  --intervals ${REF}/GRCh38/wgs_autosomal.bed \
  --pon ${PROJ}/pon/pon/giab_cnv_pon.hdf5 \
  --eh_catalog assets/eh_catalog.json \
  --annotsv_db ${REF}/annotsv/Annotations_Human \
  --sv_pon ${PROJ}/pon/sv_pon/giab_sv_pon.bed \
  --outdir ${PROJ}/results_SAMPLEID \
  -work-dir ${PROJ}/work_SAMPLEID \
  > ${TMP}/SAMPLEID_run1.log 2>&1 &

# PON build (already complete — only re-run if GIAB BAMs change)
NXF_ANSI_LOG=false nohup nextflow run workflows/pon_build.nf -profile docker \
  --input validation/giab_samplesheet.csv \
  --ref_fasta ${REF}/GRCh38/hg38.fa \
  --intervals ${REF}/GRCh38/wgs_autosomal.bed \
  --outdir ${PROJ}/pon \
  -work-dir ${PROJ}/work_pon \
  -resume > ${TMP}/pon_run.log 2>&1 &

# Check pipeline progress (update log path to latest run)
tail -20 ${TMP}/smn_runN.log

# Clean up a completed sample's work dir (safe once results are published)
rm -rf ${PROJ}/work_SAMPLEID
```

**PON location:** `${PROJ}/pon/pon/giab_cnv_pon.hdf5` (446 MB, built from HG001-HG007)

**Canonical reference:** `${REF}/GRCh38/hg38.canonical.fa` (chr1-22+X+Y+M, 1.49 GB). Use for FASTQ inputs so aligned BAMs contain only canonical chromosomes. The pipeline automatically skips FILTER_CHROMS for FASTQ-derived BAMs (saves ~25 min/sample). BWA-MEM2 index at same path prefix. BAM inputs always run FILTER_CHROMS regardless of reference used.

## Python Tests

```bash
# Run all tests
pytest tests/

# Run a single test file
pytest tests/test_cnv_consensus.py

# Run a specific test
pytest tests/test_cnv_consensus.py::test_reciprocal_overlap
```

The Python scripts in `bin/` are tested with pytest but execute inside `svcaller/utils:1.2` Docker container during the pipeline run.

## Samplesheet Format

CSV with header: `sample,fastq_1,fastq_2,bam`. Each row provides either a FASTQ pair or a pre-aligned BAM (leave the other columns blank).

```
sample,fastq_1,fastq_2,bam
HG002,/path/HG002_R1.fq.gz,/path/HG002_R2.fq.gz,
HG003,,,/path/HG003.bam
```

## Architecture

```
main.nf                          # Entry: parse samplesheet, set up channels, call SVCALLER
└── workflows/svcaller.nf        # Top-level orchestration
    ├── subworkflows/preprocess.nf   # M1: BWA-MEM2 align → SAMTOOLS_SORT → Picard MarkDup → Mosdepth QC
    ├── subworkflows/sv_calling.nf   # M2: Manta + Delly + GRIDSS + Scramble + MELT (parallel) → Jasmine merge; ExpansionHunter (STRs)
    ├── subworkflows/cnv_calling.nf  # M3: CNVpytor + GATK gCNV → cnv_consensus.py (reciprocal overlap merge)
    ├── subworkflows/smn_calling.nf  # M4: SMNCopyNumberCaller
    ├── subworkflows/annotate.nf     # M5: AnnotSV
    └── subworkflows/report.nf       # M6/M7: pycirclize Circos SVG + optional Truvari bench → HTML report
```

**Key design points:**
- M2, M3, M4 run in parallel on the same BAM channel after preprocessing.
- M3 case mode runs `GATK_PREPROCESS_INTERVALS` (bin-length 1000) before `CollectReadCounts` — must match PON build intervals.
- SV merge (Jasmine) requires Manta+Delly+GRIDSS to succeed; Scramble and MELT are optional (only added to vcf_list.txt if they have calls).
- `TRA_CONSENSUS` (`bin/tra_consensus.py`) runs after Jasmine: Jasmine does not co-cluster interchromosomal breakends, so every translocation lands at SUPP=1 and the SUPP≥2 gate erased all of them (0/26 on COLO829). This step re-clusters TRA by breakend proximity (`--tra_window`, default 1000 bp) across callers and recomputes SUPP/SUPP_VEC → 24/26 recall on COLO829. It (not JASMINE_MERGE) publishes the final `sv_merged.vcf.gz`. GRIDSS still drops inter-chr pairs in `gridss_convert_bnd.py`, so it contributes no TRA vote yet (the 2 missed are Delly-only).
- CNV consensus (`bin/cnv_consensus.py`) uses reciprocal overlap ≥0.5 to call `BOTH`/`HIGH` calls; GATK-only calls with quality ≥30 are included as `MEDIUM`.
- Mosdepth halts the pipeline if coverage < `params.min_depth` (default 30).
- Optional inputs (PoN, intervals, AnnotSV DB, GIAB truth) use `Channel.value(file("NO_PON"))` / `Channel.empty()` sentinel patterns. ANNOTSV emits a stub empty TSV when `--annotsv_db` is not provided.
- REPORT workflow joins 9 channels into BUILD_HTML_REPORT: sv_tsv, cnv_bed, smn_tsv, circos_svg, benchmark_json, sizebin_json, coverage_summary, picard_metrics, str_vcf. Each optional channel uses `remainder: true` join + `?: file("NO_FILE")` fallback.
- Truvari runs overall + 4 size bins (50–300 bp, 300 bp–1 kb, 1–10 kb, >10 kb) — both JSONs wired to HTML report.

## Known Issues / Environment Notes

- **PON built without GC correction**: `CreateReadCountPanelOfNormals` omits `--annotated-intervals` because the GRCh38 `.dict` file has alphabetical chromosome order while BAM headers use numeric order — GATK dict comparison fails. Acceptable for WGS at 30x.
- **Nextflow channel exhaustion**: shared reference files (FASTA, FAI, dict, bwt_index) must use `Channel.value()` not `Channel.fromPath()` — queue channels are consumed after the first subworkflow and all subsequent subworkflows receive nothing. Fixed in `main.nf`, `pon_build.nf`, and `cnv_calling.nf`.
- **Meta-map consistency for exact joins**: REPORT's `ch_report_in` uses exact `.join()` on the full `meta` map, so any channel whose meta differs by even one key is silently dropped (no report produced). `ch_final_bam` adds `needs_chr_filter` (preprocess.nf) which propagates to coverage/flagstat/insert_size/SV/CNV/SMN, but `PICARD_MARKDUP.out.metrics` is produced before the tag — it must be re-tagged (`needs_chr_filter: false`) to match, or every FASTQ-derived sample's report silently vanishes. When adding meta keys mid-pipeline, tag ALL sibling channels or normalize joins to `meta.id`.
- **BWA-MEM2 index path**: The bwt_index is staged as a directory in the work dir. BWA-MEM2 must be called with `${bwt_index}/${fasta}` (not just `${fasta}`) so it finds index files inside the staged directory rather than the work dir root.
- **samtools not in bwa-mem2 container**: `quay.io/biocontainers/bwa-mem2:2.2.1--he70b90d_8` does not include samtools. Sorting is done in a separate `SAMTOOLS_SORT` process using `quay.io/biocontainers/samtools:1.23.1--ha83d96e_0`.
- **annotsv not in conda env**: runs via Docker (`quay.io/biocontainers/annotsv:3.4.6--py313hdfd78af_0`) or gracefully skipped (emits empty TSV header) when `--annotsv_db` not provided.
- **AnnotSV `-annotationsDir` path**: pass `\$(dirname ${annotsv_db})` (parent of `Annotations_Human`), not `${annotsv_db}` itself — AnnotSV appends `Annotations_Human` internally.
- **Delly BCF format**: Delly 1.2.6 outputs BCF binary even with `.vcf` extension. `DELLY_MERGE` uses `bcftools concat | bcftools sort` in `broadinstitute/gatk:4.5.0.0` container (has bcftools 1.13).
- **Jasmine unsorted output**: Jasmine does not sort its merged VCF. `JASMINE_MERGE` runs `sort -k1,1 -k2,2n` after Jasmine before `bgzip | tabix`.
- **svcaller/utils:1.1**: rebuilt from `Dockerfile.utils` to add `COPY assets/ /usr/local/assets/` (report template) and fix STR VCF gzip reading and null Truvari precision/recall values.
- **samtools flagstat**: wired and working — `mapped_pct` and `dup_rate` both parsed from `HG002.flagstat.txt` and shown in HTML QC section.
- **MELT container**: `svcaller/melt:2.2.2` — must be built locally from `MELTv2.2.2.tar.gz` (MELT.jar requires registration at melt.igs.umaryland.edu; not in bioconda/biocontainers). Build: `docker build -f Dockerfile.melt -t svcaller/melt:2.2.2 .`. Requires bowtie2 in PATH (included in container). ME types run alphabetically: ALU→HERVK→LINE1→SVA (~2h total at 30×).
- **MELT INFO headers**: MELT outputs many INFO fields (DIFF/LP/RP/RA/PRIOR/SR/MEINFO etc.) that Jasmine drops from the merged VCF header, causing bcftools/Truvari to fatal-exit. Fixed: `call.nf` strips INFO to SVTYPE/MEITYPE/SVLEN/END; `merge.nf` injects `##INFO=<ID=MEITYPE,...>` before `#CHROM`.
- **MELT -n argument**: takes gene annotation BED12 (`Hg38.genes.bed` inside container at `/opt/melt/add_bed_files/Hg38/Hg38.genes.bed`), NOT prior insertion sites. This argument is mandatory in v2.2.2.

## Python Scripts (`bin/`)

All run inside `svcaller/utils:1.1`. Each is a standalone CLI tool:

| Script | Purpose |
|---|---|
| `cnv_consensus.py` | Merge CNVpytor TSV + GATK gCNV SEG → consensus BED |
| `html_report.py` | Assemble per-sample HTML report from all sections |
| `smn_report.py` | Generate SMN1/SMN2 HTML section |
| `circos_plot.py` | Generate SVG circos plot via pycirclize |
| `parse_samplesheet.py` | Samplesheet parsing utility |

## Module Conventions

Each module under `modules/<tool>/` follows: `input` tuple → `script` block → `output` tuple with named `emit`. Resource labels (`process_single`, `process_low`, `process_medium`, `process_high`, `process_gridss`) map to CPU/memory tiers in `conf/base.config` (single 1c/6G, low 4c/12G, medium 8c/36G, high 16c/32G, gridss 16c/32G — memory scales ×`task.attempt` on retry, clamped by `--max_*`). Retry on OOM exit codes (137, 143, 104, 134, 139) is automatic. `GRIDSS_CALL` has `maxRetries=3` (32→64→96 GB); `SVABA_CALL` pins CPUs at 16 (memory-bound).

**storeDir caches** (survive `nextflow clean`, reused across samples against the same reference — do not delete between runs): `SAMTOOLS_FILTER_CHROMS` → `${outdir}/.cache/filter_chroms`, `GRIDSS_SETUP` → `${outdir}/cache/gridss_ref`, `GATK_PREPROCESS_INTERVALS` → `${outdir}/cache/gatk_preprocess`. **Pre-flight:** `VALIDATE_REF_BAM` (start of M2) fails fast if the reference has contigs the BAM lacks — catches the canonical-vs-full hg38 mismatch in seconds instead of a silent Manta crash ~4 h in.

**Storage/optimization:** MOSDEPTH runs `--no-per-base` (skips the unused ~4.5 GB/sample per-base BED; only `regions.bed.gz` is consumed). Reclaim orphaned work dirs from failed/superseded runs with `bash bin/nf-cleanup.sh --reclaim` (dry-run; `--force` to delete) — keeps the latest successful run for `-resume`, refuses while a run is active. Scatter-gather is safe only for depth/per-locus tools (Delly already internal; GATK counts/CNVpytor/EH shardable); Manta/GRIDSS/SvABA need whole-genome context — do not chunk. Do NOT `docker rmi`/`prune -a` the locally-built `svcaller/{melt:2.2.2,utils,smncopynum}` images (not re-pullable; `prune -f` dangling-only is safe).

## Current F1 Baseline (2026-06-08, run16: 6-caller Manta+Delly+GRIDSS+Scramble+MELT+SvABA, GRIDSS BND→SV fix)

| Benchmark | F1 | Precision | Recall | TP-base | FP | comp cnt |
|---|---|---|---|---|---|---|
| GIAB T2TQ100-V1.0 | 0.378 | 0.733 | 0.255 | 7571 | 2604 | 9739 |
| GIAB v5.0q | 0.383 | 0.738 | 0.259 | 7286 | 2449 | 9333 |

GRIDSS contributes 3095 SVs (14% of 21,735 merged): 2713 DEL, 246 DUP, 136 INV.
Next target: improve recall (currently ~25%) — low-QUAL Manta/Delly rescue, or soft GRIDSS QUAL floor.

## Validation

```bash
# Download GIAB truth files to ${REF}/GIAB/
bash validation/download_refs.sh

# Standalone Truvari benchmark (runs Docker internally)
bash validation/giab_benchmark.sh HG002 /path/to/HG002.sv_merged.vcf.gz

# GIAB sample sheet for pipeline-integrated benchmarking
validation/giab_samplesheet.csv
```

## Key Parameters

| Param | Default | Notes |
|---|---|---|
| `--input` | required | Samplesheet CSV |
| `--ref_fasta` | required | GRCh38 FASTA (index .fai and .0123 inferred from path) |
| `--pon` | null | GATK gCNV Panel of Normals HDF5 |
| `--intervals` | null | Target capture BED |
| `--annotsv_db` | null | AnnotSV database directory — use `${REF}/annotsv/Annotations_Human` |
| `--giab_truth` | null | GIAB T2TQ100-V1.0 truth VCF.gz (enables Truvari T2T benchmark) |
| `--giab_truth_v5q` | null | GIAB v5.0q truth VCF.gz (enables Truvari v5q benchmark) |
| `--skip_melt` | false | Skip MELT MEI calling (saves ~2h; use when MELT container unavailable) |
| `--melt_refs` | null | Path to MELT me_refs dir (auto-detected from container if unset) |
| `--tra_window` | 1000 | Breakend-proximity window (bp) for cross-caller translocation consensus (TRA_CONSENSUS) |
| `--min_depth` | 25 | Mosdepth coverage threshold |
| `--outdir` | results | Output directory |
| `--auto_cleanup` | false | Delete `-work-dir` on successful completion (drops `-resume` cache; one-shot runs only). Otherwise run `bash bin/nf-cleanup.sh <sampleId>` post-run. |
| `--max_cpus` / `--max_memory` / `--max_time` | 64 / 120.GB / 240.h | Per-process resource caps. Lower on small machines (e.g. `--max_memory 32.GB`) to avoid over-subscription. |
| `--utils_container` | `svcaller/utils:1.2` | Container for Python bin/ scripts |

# context-mode — MANDATORY routing rules

You have context-mode MCP tools available. These rules are NOT optional — they protect your context window from flooding. A single unrouted command can dump 56 KB into context and waste the entire session.

## BLOCKED commands — do NOT attempt these

### curl / wget — BLOCKED
Any Bash command containing `curl` or `wget` is intercepted and replaced with an error message. Do NOT retry.
Instead use:
- `ctx_fetch_and_index(url, source)` to fetch and index web pages
- `ctx_execute(language: "javascript", code: "const r = await fetch(...)")` to run HTTP calls in sandbox

### Inline HTTP — BLOCKED
Any Bash command containing `fetch('http`, `requests.get(`, `requests.post(`, `http.get(`, or `http.request(` is intercepted and replaced with an error message. Do NOT retry with Bash.
Instead use:
- `ctx_execute(language, code)` to run HTTP calls in sandbox — only stdout enters context

### WebFetch — BLOCKED
WebFetch calls are denied entirely. The URL is extracted and you are told to use `ctx_fetch_and_index` instead.
Instead use:
- `ctx_fetch_and_index(url, source)` then `ctx_search(queries)` to query the indexed content

## REDIRECTED tools — use sandbox equivalents

### Bash (>20 lines output)
Bash is ONLY for: `git`, `mkdir`, `rm`, `mv`, `cd`, `ls`, `npm install`, `pip install`, and other short-output commands.
For everything else, use:
- `ctx_batch_execute(commands, queries)` — run multiple commands + search in ONE call
- `ctx_execute(language: "shell", code: "...")` — run in sandbox, only stdout enters context

### Read (for analysis)
If you are reading a file to **Edit** it → Read is correct (Edit needs content in context).
If you are reading to **analyze, explore, or summarize** → use `ctx_execute_file(path, language, code)` instead. Only your printed summary enters context. The raw file content stays in the sandbox.

### Grep (large results)
Grep results can flood context. Use `ctx_execute(language: "shell", code: "grep ...")` to run searches in sandbox. Only your printed summary enters context.

## Tool selection hierarchy

1. **GATHER**: `ctx_batch_execute(commands, queries)` — Primary tool. Runs all commands, auto-indexes output, returns search results. ONE call replaces 30+ individual calls.
2. **FOLLOW-UP**: `ctx_search(queries: ["q1", "q2", ...])` — Query indexed content. Pass ALL questions as array in ONE call.
3. **PROCESSING**: `ctx_execute(language, code)` | `ctx_execute_file(path, language, code)` — Sandbox execution. Only stdout enters context.
4. **WEB**: `ctx_fetch_and_index(url, source)` then `ctx_search(queries)` — Fetch, chunk, index, query. Raw HTML never enters context.
5. **INDEX**: `ctx_index(content, source)` — Store content in FTS5 knowledge base for later search.

## Subagent routing

When spawning subagents (Agent/Task tool), the routing block is automatically injected into their prompt. Bash-type subagents are upgraded to general-purpose so they have access to MCP tools. You do NOT need to manually instruct subagents about context-mode.

## Output constraints

- Keep responses under 500 words.
- Write artifacts (code, configs, PRDs) to FILES — never return them as inline text. Return only: file path + 1-line description.
- When indexing content, use descriptive source labels so others can `ctx_search(source: "label")` later.

## ctx commands

| Command | Action |
|---------|--------|
| `ctx stats` | Call the `ctx_stats` MCP tool and display the full output verbatim |
| `ctx doctor` | Call the `ctx_doctor` MCP tool, run the returned shell command, display as checklist |
| `ctx upgrade` | Call the `ctx_upgrade` MCP tool, run the returned shell command, display as checklist |
