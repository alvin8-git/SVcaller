# SVcaller Architecture Reference

Technical reference for every module, process, and data flow in the pipeline.

## Top-Level Flow

```
main.nf
└── SVCALLER (workflows/svcaller.nf)
    ├── PREPROCESS      (subworkflows/preprocess.nf)   ← M1
    ├── SV_CALLING      (subworkflows/sv_calling.nf)   ← M2  ┐
    ├── CNV_CALLING     (subworkflows/cnv_calling.nf)  ← M3  ├─ parallel on same BAM
    ├── SMN_CALLING     (subworkflows/smn_calling.nf)  ← M4  ┘
    ├── ANNOTATE        (subworkflows/annotate.nf)     ← M5
    └── REPORT          (subworkflows/report.nf)       ← M6/M7
```

M2, M3, and M4 run in parallel on the preprocessed BAM channel. M5 and M6/M7 run after M2.

## M1: Preprocessing (`subworkflows/preprocess.nf`)

**Purpose:** Convert raw input (FASTQ or BAM) to a sorted, duplicate-marked BAM with QC metrics.

### Processes

| Process | Module | Tool | Input | Output |
|---------|--------|------|-------|--------|
| `BWAMEM2_ALIGN` | `modules/bwamem2/align.nf` | BWA-MEM2 2.2.1 | FASTQ pair, reference FASTA + index | Unsorted BAM |
| `SAMTOOLS_SORT` | `modules/samtools/sort.nf` | samtools 1.23.1 | Unsorted BAM | Sorted BAM |
| `FILTER_CHROMS` | `modules/samtools/filter_chroms.nf` | samtools | Sorted BAM | BAM with canonical chroms only |
| `PICARD_MARKDUP` | `modules/picard/markduplicates.nf` | Picard | Sorted/filtered BAM | Duplicate-marked BAM + metrics |
| `MOSDEPTH` | `modules/mosdepth/coverage.nf` | mosdepth | Final BAM | Summary TXT + 50kb regions BED |
| `FASTQC` | `modules/fastqc/qc.nf` | FastQC | FASTQ | QC reports for MultiQC |
| `PICARD_INSERT_SIZE` | `modules/picard/insert_size.nf` | Picard | Final BAM | Insert size metrics |
| `SAMTOOLS_FLAGSTAT` | `modules/samtools/flagstat.nf` | samtools | Final BAM | Flagstat TXT |

**FILTER_CHROMS** is skipped for FASTQ samples aligned to `hg38.canonical.fa`. BAM inputs always run FILTER_CHROMS. The filter removes both reads and `@SQ` header lines for non-canonical chromosomes — keeping non-canonical `@SQ` entries even without corresponding reads causes Manta to silently produce 0 assembly output (scanner generates candidates, assembly phase crashes internally with no error message).

**MOSDEPTH** uses `--by 50000` to produce 50 kb window coverage used in the Circos depth ring. It also checks mean depth against `--min_depth` and halts the pipeline if coverage is insufficient.

**Outputs emitted:** `ch_final_bam`, `ch_markdup_metrics`, `ch_coverage` (summary TXT), `ch_depth_bed` (50 kb BED), `ch_flagstat`, `ch_insert_size`, `ch_multiqc_files`.

## M2: SV Calling (`subworkflows/sv_calling.nf`)

**Purpose:** Call structural variants using complementary callers and merge into a single VCF.

### Pre-flight: `VALIDATE_REF_BAM`

Before any caller runs, `VALIDATE_REF_BAM` (`modules/samtools/validate_ref_bam.nf`) checks that the BAM's chromosome set is a superset of the reference FASTA's. If the reference carries contigs absent from the BAM, it fails fast with a hint to use `hg38.canonical.fa` (not `hg38.fa`) for BAM inputs. This catches the mismatch at second-zero instead of letting Manta crash silently ~4 hours into the run.

### SV Callers

| Caller | Module | Detects | Container |
|--------|--------|---------|-----------|
| Manta | `modules/manta/call.nf` | DEL, DUP, INV, BND, INS | quay.io/biocontainers/manta:1.6.0 |
| Delly | `modules/delly/call.nf` + `merge.nf` | DEL, DUP, INV, TRA, INS | quay.io/biocontainers/delly:1.2.6 |
| GRIDSS | `modules/gridss/call.nf` | BND (precise breakpoints) | gridss/gridss:2.13.2 |
| Scramble | `modules/scramble/call.nf` | MEI (ALU, L1, SVA insertions) | quay.io/biocontainers/scramble:1.0.2 |
| MELT | `modules/melt/call.nf` | MEI (ALU, HERVK, LINE1, SVA) | svcaller/melt:2.2.2 (local build) |
| SvABA | `modules/svaba/call.nf` | DEL, DUP, INV, INS (local assembly) | quay.io/biocontainers/svaba:1.2.0 |

**SvABA** performs local assembly and calls `bwa_idx_load_from_disk` internally, so it
requires the **classic** BWA index (`ref_fasta.{amb,ann,bwt,pac,sa}`) staged next to the
reference — a *different* format from the bwa-mem2 alignment index (`.0123`/`.bwt.2bit.64`),
not interchangeable. `SVABA_CALL` declares `path bwa_index`; the index is threaded from
`main.nf` (`ch_bwa_index`, prefix `--bwa_index`, default `ref_fasta`) through
`subworkflows/sv_calling.nf` so Nextflow symlinks all five files into the task work dir.
Unless `--skip_svaba` is set, `main.nf` fails loud if the index is absent, pointing the
operator at `bwa index <ref>`. Historically these files were never staged and a
`2>&1 || true` masked the resulting crash, so SvABA silently produced nothing (see
[CHANGES.md](CHANGES.md), 2026-07-15).

**Delly** outputs BCF binary even with `.vcf` extension. `DELLY_MERGE` uses `bcftools concat | bcftools sort` (inside the GATK container which includes bcftools 1.13) to convert and sort.

**GRIDSS** can be skipped with `--skip_gridss true` or run in tiered mode (`--tiered_gridss true`) where it only processes Manta residual regions (smaller input, faster runtime). GRIDSS BND pairs are converted to DEL/DUP/INV calls by `bin/gridss_convert_bnd.py` before Jasmine merge. `GRIDSS_CALL` carries `maxRetries = 3` (4 attempts: 32 → 64 → 96 GB) because it is the caller most prone to silent OOM. `GRIDSS_SETUP` writes its reference setup to a `storeDir` cache (see [storeDir caches](#storedir-caches)) so it runs once per reference, not once per sample.

**Scramble** calls mobile element insertions (MEI). Each MEI type gets a canonical SVLEN assigned at the Jasmine merge step: ALU=300 bp, L1=6000 bp, SVA=1500 bp.

**MELT** (Mobile Element Locator Tool) calls four ME types: ALU, HERVK, LINE1, SVA. Runs alphabetically (~2 h at 30×). Requires `svcaller/melt:2.2.2` built locally from `MELTv2.2.2.tar.gz` (MELT.jar requires registration; not in bioconda/biocontainers). Skip with `--skip_melt true`. MELT INFO headers (DIFF/LP/RP/RA/PRIOR/SR/MEINFO etc.) are stripped to SVTYPE/MEITYPE/SVLEN/END before Jasmine to prevent header mismatch crashes.

**ExpansionHunter** runs in parallel on the BAM to call short tandem repeat (STR) expansions. Uses `assets/eh_catalog.json`. Output: `ch_str_vcf`.

### Merge: Jasmine

`JASMINE_MERGE` (modules/jasmine/merge.nf) takes all caller VCFs and produces a single merged VCF with:
- `SUPP_VEC` field: bitmask indicating which callers support each call (e.g., `1110` = Manta+Delly+GRIDSS, `10000` = MELT only)
- `SUPP`: count of supporting callers

Jasmine output is not sorted. The module sorts it with `sort -k1,1 -k2,2n` before `bgzip | tabix`.

**Manta, Delly, and GRIDSS must all succeed** for Jasmine to run (inner join). Scramble and MELT are optional — they are added to `vcf_list.txt` only if they produced calls. If a mandatory caller fails, the sample fails.

**Outputs emitted:** `ch_sv_vcf`, `ch_sv_tbi`.

Also in M2: **ExpansionHunter** runs in parallel on the BAM to call short tandem repeat (STR) expansions. Output: `ch_str_vcf`.

## M3: CNV Calling (`subworkflows/cnv_calling.nf`)

**Purpose:** Call copy number variants using read-depth methods and merge into a consensus BED.

### CNV Callers

| Caller | Module | Method | Requires |
|--------|--------|--------|---------|
| CNVpytor | `modules/cnvpytor/call.nf` | Read depth, population statistics | Nothing extra |
| GATK gCNV | `modules/gatk/gcnv_call.nf` | GC-corrected read depth vs. PON | `--pon`, `--intervals` |

**GATK gCNV pipeline steps:**
1. `GATK_PREPROCESS_INTERVALS` — bins the genome into 1000 bp intervals
2. `GATK_COLLECT_READ_COUNTS` — counts reads per interval
3. `GATK_DENOISE_READ_COUNTS` — denoises against the PON
4. `GATK_MODEL_SEGMENTS` — segments denoised copy ratios
5. `GATK_CALL_COPY_RATIO_SEGMENTS` — calls CNV segments

**Consensus merge** (`bin/cnv_consensus.py`): Overlapping calls from both callers are merged using reciprocal overlap ≥0.5. Confidence levels:
- `BOTH` — both callers agree with ≥50% reciprocal overlap
- `HIGH` — GATK-only with quality ≥30
- `MEDIUM` — GATK-only with quality <30 or CNVpytor-only with strong signal
- `LOW` — CNVpytor-only with weak signal

**Output emitted:** `ch_cnv_bed`.

## M4: SMN Calling (`subworkflows/smn_calling.nf`)

**Purpose:** Determine SMN1 and SMN2 copy numbers (SMA diagnosis support).

Tool: SMNCopyNumberCaller. Runs from the full pre-filtered BAM (before canonical chrom filtering) because SMN loci span the chr5 boundary in some reference builds.

**Output emitted:** `ch_smn_tsv` — TSV with per-sample SMN1/SMN2 CN.

## M5: Annotation (`subworkflows/annotate.nf`)

**Purpose:** Annotate merged SV calls with gene overlaps, population frequencies, and ACMG classification.

### Processes

| Process | Tool | Input | Output |
|---------|------|-------|--------|
| `ANNOTSV` | AnnotSV 3.4.6 | `sv_merged.vcf.gz` | TSV with gene/OMIM/ACMG_class per SV |
| `GNOMAD_SV_FILTER` | Python (bin/) | AnnotSV TSV | Frequency-filtered TSV |

**AnnotSV** runs with `-annotationMode both` to produce:
- `full` rows: one per SV (used for counting and ranking)
- `split` rows: one per transcript overlap (used for gene-level details)

**AnnotSV subdirectory bug workaround:** AnnotSV writes output to a date-stamped subdirectory (`YYYYMMDD_AnnotSV/`) despite being given an absolute `-outputFile` path. The module uses a `find -maxdepth 2` fallback to locate and move the output file.

**Two TSVs are emitted:**
- `ch_annotsv_tsv` — raw AnnotSV output (all SVs, passed to Circos for gene loci + ACMG rings)
- `ANNOTATE.out.tsv` / `ch_sv_tsv` — gnomAD-filtered TSV (passed to HTML report for clinical tables)

## M6/M7: Report (`subworkflows/report.nf`)

**Purpose:** Assemble the per-sample HTML report, Circos plot, optional Truvari benchmark, and MultiQC report.

### Processes

| Process | Container | Input | Output |
|---------|-----------|-------|--------|
| `CIRCOS_PLOT` | svcaller/utils:1.3 | sv_vcf, cnv_bed, str_vcf, depth_bed, annotsv_tsv, cytobands | SVG + PNG |
| `TRUVARI_BENCH` | quay.io/.../truvari | sv_vcf, truth_vcf, truth_bed | summary.json + sizebin.json |
| `MULTIQC` | quay.io/.../multiqc | QC files from M1 | HTML report |
| `BUILD_HTML_REPORT` | svcaller/utils:1.3 | All above + coverage, metrics, smn_tsv | .report.html + .variants.xlsx |

### Circos Plot Ring Layout (outer → inner)

Remodeled into a "genome fingerprint" — links-first, balanced rings (the plot is a gestalt overview; exact coordinates/CN/confidence live in the report sections + `.xlsx`):

| Ring | Radius | Content |
|------|--------|---------|
| Chromosomes | 95–100 | Ideograms chr1–22+X+Y |
| Copy-ratio heatmap | 88–95 | Log₂ depth ratio in **1 Mb bins** (re-binned from mosdepth 50 kb), diverging RdBu_r (blue=loss, white=neutral, red=gain); ±1 ≈ CN1..CN4. Uniform genome reads near-white; CNVs pop as colour blocks |
| CNV blocks | 81–86 | CNV consensus calls **≥1 Mb** (DUP=red, DEL=blue); smaller calls are in the heatmap + xlsx |
| STR loci | 74–79 | ExpansionHunter (brown) + STRling novel (orange), full-height barcode ticks |
| Clinical SV | 67–72 | Top 30 SVs by AnnotSV score, **coloured by ACMG class** (5=red, 4=orange, 3=grey, else by SV type); gold = SMN locus. Merges the former gene-loci + ACMG rings |
| SV links | 0–65 | DEL/DUP/INV ≥50 kb and BND/TRA, **all requiring multi-caller support (SUPP_VEC≥2)** — same consensus as the callset/xlsx, so the plot never shows SVs the pipeline filtered out. Cap 150, up to 110 interchromosomal, drawn thicker/opaque as the sample signature (sparse for normal genomes; a real web only for samples with consensus-supported rearrangements). Point features use a ~2 Mb min-width floor. Note: GIAB truth is INS/DEL only, so interchromosomal recall is unmeasured by GIAB |

### Channel Join Pattern

The `REPORT` workflow joins 9 mandatory channels with `join()` (inner join) and 3 optional channels with `join(remainder: true)` + `?: file("NO_FILE")` fallback. Any meta-map mismatch on mandatory channels silently drops the sample — all channels must carry identical meta maps.

**Output:** `html = BUILD_HTML_REPORT.out.html`, `multiqc = MULTIQC.out.html`.

## Python Scripts (`bin/`)

Scripts run inside `svcaller/utils:1.3`. The `export PATH=${projectDir}/bin:$PATH` override in process scripts ensures the host's `bin/` takes precedence over container-baked versions.

| Script | CLI entrypoint | Key logic |
|--------|---------------|-----------|
| `cnv_consensus.py` | `cnv_consensus.py --cnvpytor ... --gatk ...` | Reciprocal overlap merge; confidence scoring |
| `html_report.py` | `html_report.py --sample ... --sv-tsv ...` | Jinja2 template rendering; parses mosdepth/Picard/flagstat/AnnotSV |
| `smn_report.py` | `smn_report.py --tsv ... --sample ...` | Generates HTML section for SMN copy number |
| `circos_plot.py` | `circos_plot.py --sv-vcf ... --depth-bed ...` | pycirclize-based ring plot; log₂ depth scatter; AnnotSV gene/ACMG rings |
| `parse_samplesheet.py` | Internal | Validates CSV; emits per-row channel entries |

## Configuration Files

| File | Purpose |
|------|---------|
| `conf/base.config` | Resource tiers (CPU/memory per label), per-name overrides, OOM retry logic. `GRIDSS_CALL` `maxRetries = 3`; `SVABA_CALL` CPUs fixed at 16. |
| `conf/docker.config` | Docker run options (`--user`, `--volume`) |
| `conf/local.config` | Conda-based execution profile (no Docker) — points processes at the `svcaller` conda env and prepends `bin/` to PATH |
| `conf/test.config` | Minimal smoke-test resources |
| `nextflow.config` | Parameter defaults, profile definitions, manifest, `env { TMPDIR }`, the `auto_cleanup` `workflow.onComplete` hook, and `--max_cpus`/`--max_memory`/`--max_time` caps |

## storeDir caches

Three processes persist outputs to a `storeDir` so they survive `nextflow clean` and are reused across runs and samples against the same reference. The caches are keyed on their inputs (reference + intervals); deleting them only forces recomputation.

| Process | storeDir path | Cached artifact |
|---------|---------------|-----------------|
| `SAMTOOLS_FILTER_CHROMS` | `${outdir}/.cache/filter_chroms` | Chrom-filtered BAM (BAM inputs) |
| `GRIDSS_SETUP` | `${outdir}/cache/gridss_ref` | GRIDSS reference setup |
| `GATK_PREPROCESS_INTERVALS` | `${outdir}/cache/gatk_preprocess` | Binned `interval_list` |

See [Storage & Cache Management](reference-parameters.md#storage--cache-management) for cleanup workflow and `--auto_cleanup`.

## Related

- [Parameter reference](reference-parameters.md) — all CLI flags and samplesheet columns
- [Design decisions explained](explanation-design.md) — why the architecture is this way
- [How to interpret the HTML report](howto-interpret-report.md)
