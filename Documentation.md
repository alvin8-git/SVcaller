# SVcaller Pipeline Documentation

**Version:** 1.0.0  
**Target genome:** GRCh38  
**Minimum coverage:** ≥30×  
**Execution engine:** Nextflow DSL2 with Docker (`-profile docker`)

---

## Table of Contents

1. [Why SVcaller? (Battlecard)](#1-why-svcaller-battlecard)
2. [Overview](#2-overview)
3. [M1: Preprocessing](#3-m1-preprocessing)
4. [PON Build (pre-run, one-time)](#4-pon-build-pre-run-one-time)
5. [M2: SV Calling](#5-m2-sv-calling)
6. [M3: CNV Calling](#6-m3-cnv-calling)
7. [M4: SMN Copy Number Calling](#7-m4-smn-copy-number-calling)
8. [M5: Annotation](#8-m5-annotation)
9. [M6/M7: Reporting](#9-m6m7-reporting)
10. [Known Limitations](#10-known-limitations)
11. [Alternative Pipeline Strategies](#11-alternative-pipeline-strategies)
12. [Key Parameters Reference](#12-key-parameters-reference)

---

## 1. Why SVcaller? (Battlecard)

### Strengths and weaknesses at a glance

| | SVcaller | nf-core/sarek | DRAGEN (Illumina) | Manual (per-caller) |
|---|---|---|---|---|
| **SV calling** | 3-caller ensemble (Manta + Delly + GRIDSS) | Manta only (default) | Proprietary SVs | User-defined |
| **CNV calling** | Dual method (CNVpytor + GATK gCNV PON) | GATK gCNV or Control-FREEC | DRAGEN CNV | User-defined |
| **STR genotyping** | ExpansionHunter (catalog-based) | — | ExpansionHunter | User-defined |
| **SMN1/SMN2 CN** | SMNCopyNumberCaller (built-in) | — | — | Separate tool, manual |
| **SV annotation** | AnnotSV (clinical-grade, gnomAD filter) | VEP | — | Separate tool |
| **GIAB benchmarking** | Built-in Truvari (5 size bins) | — | — | Manual |
| **Per-sample report** | HTML (SVs, CNVs, STRs, SMN, QC, circos) | MultiQC only | PDF summary | Manual |
| **License** | Open source | Open source | Commercial | Open source |
| **Hardware** | Local/HPC, Docker | Local/HPC/Cloud, Docker/Singularity | Illumina hardware only | User-defined |
| **Setup effort** | PON build required (one-time) | Reference files only | Minimal | High |
| **Short reads** | ✅ | ✅ | ✅ | ✅ |
| **Long reads** | ❌ | ❌ (separate wf) | ✅ | ✅ |
| **Somatic mode** | ❌ | ✅ | ✅ | User-defined |
| **Min. coverage** | 30× | 30× | 20× | Caller-dependent |

### Pros

| Strength | Detail |
|---|---|
| **Higher SV sensitivity** | Three independent callers (split-read, paired-end, assembly) catch different SV classes; Jasmine inner-join merge reduces false positives while retaining multi-evidence calls |
| **SMN1/SMN2 built-in** | Most pipelines omit SMN copy number entirely; SVcaller integrates it for SMA diagnostics without a separate workflow |
| **Clinically oriented annotation** | AnnotSV adds OMIM, ClinVar, DGV, ACMG pathogenicity fields directly to the SV VCF — essential for diagnostic interpretation |
| **Integrated benchmarking** | Truvari against GIAB truth sets with five size-bin breakdowns (50–300 bp, 300 bp–1 kb, 1–10 kb, >10 kb, overall) is automatic when `--giab_truth` is provided |
| **Unified HTML report** | SVs, CNVs, STRs, SMN CN, QC metrics, and a genome-wide circos plot in one report per sample — no manual assembly required |
| **Reproducible containers** | All tools pinned to verified quay.io biocontainer tags; `-profile docker` is fully tested end-to-end |
| **Dual CNV methods** | CNVpytor (read-depth HMM) + GATK gCNV (PON-denoised) with reciprocal-overlap consensus — complementary false-positive profiles |
| **Nextflow -resume** | Any failed step resumes from the last checkpoint; large WGS runs do not restart from scratch |

### Cons / when NOT to use SVcaller

| Limitation | Implication |
|---|---|
| **Short-read WGS only** | Not suitable for Oxford Nanopore or PacBio long reads; use Sniffles2 + PBSV + Paraphase instead |
| **No somatic/tumour-normal mode** | Designed for germline calling; paired tumour/normal analysis is not supported |
| **PON required for GATK CNV** | A panel of normals (≥7 samples recommended) must be built before the first case run — adds one-time setup overhead |
| **≥30× coverage required** | Mosdepth gates the pipeline; lower-coverage samples (e.g. 10× WGS or WES) will fail the depth check |
| **GRIDSS memory-intensive** | GRIDSS assembly requires 32 GB RAM minimum; unsuitable for machines with <32 GB |
| **No CNV for WES/panel** | CNVpytor and GATK gCNV are tuned for WGS read depth; targeted panels need CNVKit or GATK4 ExomeDepth instead |
| **GC bias correction absent** | PON was built without `--annotated-intervals` (dict ordering mismatch); may reduce CNV sensitivity in extreme GC regions |
| **Single-sample only** | Family/trio analysis (de novo filtering, inheritance phasing) is not implemented |

### Decision guide

```
Need germline SV + CNV + SMN from 30× WGS?  →  SVcaller
Need somatic calling or tumour/normal?        →  nf-core/sarek or DRAGEN
Have long-read data (ONT / PacBio)?          →  Sniffles2 + Paraphase (manual) or DRAGEN LR
Need maximum throughput, commercial support?  →  DRAGEN
Have <30× coverage or WES data?              →  nf-core/sarek or CNVKit
```

---

## 2. Overview

SVcaller is a Nextflow DSL2 pipeline for comprehensive short-read whole-genome sequencing (WGS) variant calling focused on structural variants (SVs), copy-number variants (CNVs), short tandem repeats (STRs), and SMN1/SMN2 copy number — all from a single GRCh38-aligned dataset at ≥30× coverage. The pipeline accepts either paired-end FASTQ files or pre-aligned BAM files per sample (mixed within the same run), performs alignment and duplicate marking, then fans out into three parallel variant-calling arms (M2: SVs + STRs, M3: CNVs, M4: SMN). SV calls from three independent callers (Manta, Delly, GRIDSS) are merged using Jasmine with an inner-join requirement — every sample must successfully complete all three callers before the merged VCF is emitted, providing fail-fast behavior on caller errors. CNV calls from two complementary methods (CNVpytor read-depth and GATK gCNV PON-based) are reconciled by a custom reciprocal-overlap merge script. Merged SV calls are annotated with AnnotSV and optionally benchmarked against GIAB truth sets using Truvari across five size bins. A per-sample HTML report integrating all variant types, QC metrics, a genome-wide circos plot, and optional benchmark statistics is produced for each sample.

**Pipeline architecture:**

```
main.nf
└── workflows/svcaller.nf
    ├── subworkflows/preprocess.nf      M1: FastQC → BWA-MEM2 → SAMTOOLS_SORT → Picard MarkDup → Mosdepth
    ├── subworkflows/sv_calling.nf      M2: Manta + Delly + GRIDSS (parallel) → Jasmine; ExpansionHunter
    ├── subworkflows/cnv_calling.nf     M3: CNVpytor + GATK gCNV → cnv_consensus.py
    ├── subworkflows/smn_calling.nf     M4: SMNCopyNumberCaller
    ├── subworkflows/annotate.nf        M5: AnnotSV → gnomAD SV filter
    └── subworkflows/report.nf          M6/M7: pycirclize + Truvari + MultiQC → HTML report
```

**Key design choices:**

- M2, M3, and M4 run in parallel on the same deduplicated BAM channel, maximising throughput on multi-core systems.
- Reference files (FASTA, FAI, dict, intervals) are propagated as `Channel.value()` — not queue channels — so all samples in a batch share them without channel exhaustion.
- Optional inputs (PON, intervals, AnnotSV DB, GIAB truth) use sentinel file patterns (`NO_PON`, `NO_INTERVALS`, `NO_ANNOTSV`) rather than conditional branching, keeping the workflow topology uniform.
- All containers are verified quay.io biocontainer tags; the custom Python scripts run inside `svcaller/utils:1.0`.

---

## 3. M1: Preprocessing

**Subworkflow:** `subworkflows/preprocess.nf`

The preprocessing subworkflow converts raw reads to a clean, duplicate-marked, coverage-QC-gated BAM. It accepts both FASTQ pairs and pre-aligned BAMs in the same samplesheet; FASTQ samples undergo full alignment while BAM samples bypass alignment and feed directly into duplicate marking.

### 3.1 FastQC

| Aspect | Detail |
|---|---|
| Module | `modules/fastqc/qc.nf` |
| Input | Paired FASTQ files (FASTQ-mode samples only) |
| Output | Per-sample ZIP archives collected into MultiQC |
| Purpose | Per-base quality score distributions, GC content, adapter contamination, duplicate level estimates from raw reads |

FastQC runs only on FASTQ-input samples; BAM-supplied samples skip this step since their raw reads are not available. FastQC output is collected for MultiQC aggregation at the reporting stage.

**Why FastQC:** It is the de facto standard for Illumina QC, has zero dependencies beyond Java, and its output format is natively parsed by MultiQC. Alternatives (fastp, Trimmomatic QC mode) are heavier and add trimming side effects; this pipeline intentionally delegates trimming decisions to the user upstream.

### 5.2 BWA-MEM2

| Aspect | Detail |
|---|---|
| Module | `modules/bwamem2/align.nf` |
| Input | Paired FASTQs, reference FASTA, FAI, BWA-MEM2 index directory |
| Output | Unsorted BAM with read-group tags |
| Resource label | `process_high` (maximum CPU/memory tier) |
| Key parameters | `-t ${task.cpus}` (all available CPUs); read-group `@RG\tID\tSM\tPL:ILLUMINA\tLB` set per sample |

BWA-MEM2 is a drop-in replacement for BWA-MEM that uses SIMD (AVX-512/AVX2/SSE4.1) vectorised Smith-Waterman for roughly 2× alignment speed with identical output. The read-group line is mandatory for downstream GATK tools and Picard.

**Why BWA-MEM2 over alternatives:**

| Tool | Approach | Relative speed | Notes |
|---|---|---|---|
| BWA-MEM2 | Seed-extend, FM-index | ~2× BWA-MEM | Best short-read accuracy for SVs; reads spanning breakpoints align correctly |
| BWA-MEM | Seed-extend, FM-index | 1× baseline | Still widely used; identical accuracy, slower |
| Bowtie2 | Seed-extend, BWT | Similar | Optimised for ChIP-seq/ATAC-seq; lower sensitivity for reads spanning large deletions |
| DRAGEN | Hardware-accelerated | 5–10× | FPGA-based; expensive hardware, closed-source, not portable |
| Minimap2 | Minimiser-based | Faster for long reads | Designed for long reads; suboptimal for 150 bp Illumina, lower split-read sensitivity |

**Known limitation:** The bwa-mem2 Docker container does not include samtools, requiring a separate SAMTOOLS_SORT step (see §2.3).

### 5.3 SAMTOOLS_SORT

| Aspect | Detail |
|---|---|
| Module | `modules/samtools/sort.nf` |
| Input | Unsorted BAM from BWA-MEM2 |
| Output | Coordinate-sorted BAM + BAI index |
| Resource label | `process_medium` |
| Key parameters | `-@ ${task.cpus}`, `-m 2G` per thread |

This step is a separate process because samtools is not bundled in the bwa-mem2 biocontainer image. It coordinate-sorts and indexes the BAM, making it compatible with Picard MarkDuplicates (which requires coordinate order).

### 5.4 Picard MarkDuplicates

| Aspect | Detail |
|---|---|
| Module | `modules/picard/markduplicates.nf` |
| Input | Coordinate-sorted BAM + BAI (from SAMTOOLS_SORT or pre-supplied) |
| Output | Duplicate-marked BAM + BAI + metrics TSV |
| Purpose | Flags PCR and optical duplicates to prevent false-positive variant calls |

FASTQ-derived BAMs (post-sort) and pre-supplied BAMs are mixed (`ch_all_bam = SAMTOOLS_SORT.out.bam.join(...).mix(ch_bam_in)`) before entering MarkDup, ensuring uniform treatment. Picard metrics (duplication rate) are surfaced in the HTML report QC section.

**Why Picard over sambamba markdup:** Picard is the GATK-ecosystem standard; its metrics file format is directly parsed by MultiQC. sambamba is faster but produces a different metrics format, and optical duplicate detection requires the Picard coordinate-based model for patterned flowcells (NovaSeq).

### 5.5 Mosdepth

| Aspect | Detail |
|---|---|
| Module | `modules/mosdepth/coverage.nf` |
| Input | Deduplicated BAM + BAI, `min_depth` threshold |
| Output | Summary TSV with mean coverage per chromosome |
| Threshold | Pipeline halts if mean autosomal coverage < `params.min_depth` (default: 30) |

Mosdepth is the fastest BAM-based depth calculator available for short reads. The pipeline uses it as a QC gate: if a sample fails the 30× threshold, the pipeline exits immediately rather than wasting compute on low-confidence variant calls. The summary file is passed to the HTML report.

**Why Mosdepth over samtools depth:** mosdepth produces a compact summary in seconds; `samtools depth` emits per-base output which is hundreds of times larger and slower to parse for a WGS BAM.

**Limitation:** `samtools flagstat` is not wired into the pipeline, so the mapping rate percentage shows "N/A" in the HTML QC section. Coverage depth (mosdepth) and duplication rate (Picard) are reported.

---

## 4. PON Build (pre-run, one-time)

**Workflow:** `workflows/pon_build.nf`  
**Completed PON:** `/data/alvin/SVcaller/pon/pon/giab_cnv_pon.hdf5` (446 MB, built from GIAB HG001–HG007)

The Panel of Normals (PON) is built once from a set of well-characterised normal genomes and reused across all case samples. It captures systematic read-depth biases arising from GC content, mappability, and library preparation artefacts that are constant across samples processed with the same protocol.

### 3.1 PreprocessIntervals

Bins the input BED intervals into fixed-width 1000 bp bins using GATK `PreprocessIntervals` (`--bin-length 1000 --interval-merging-rule OVERLAPPING_ONLY`). This discretisation is required before read counting and must match exactly between PON build and case-mode calling (both use the same preprocessed interval list).

### 5.2 AnnotateIntervals

Runs GATK `AnnotateIntervals` to compute GC content and mappability scores per bin. The annotation TSV is produced but **not passed to `CreateReadCountPanelOfNormals`** — see the GC-correction limitation below.

A `sed` post-processing step strips `M5:` and `UR:` fields from `@SQ` header lines in the annotation TSV to avoid a GATK dictionary-comparison error caused by mismatched metadata between the FASTA `.dict` file and BAM headers.

### 5.3 CollectReadCounts

Runs GATK `CollectReadCounts` per sample in parallel, counting reads overlapping each preprocessed bin and writing a per-sample HDF5 file. The preprocessed interval list is shared via `.first()` (converting the single-item queue channel to a value channel) so all 7 GIAB samples can each receive it without channel exhaustion.

### 5.4 CreateReadCountPanelOfNormals

Aggregates all 7 per-sample HDF5 count files with GATK `CreateReadCountPanelOfNormals` into a single PON HDF5. The PON learns the mean and variance of read depth in each bin across the cohort and is used in case mode to denoise individual sample read counts.

**Why GIAB HG001–HG007 as normals:** These are the seven NIST Genome in a Bottle reference samples with high-confidence whole-genome truth sets, high-quality deep WGS data, and broad ethnic diversity (European, Ashkenazi Jewish, Chinese, Han Chinese × 2). They are the best-characterised human genomic references available and represent typical short-read library and coverage profiles.

**GC correction omission:** `CreateReadCountPanelOfNormals` supports `--annotated-intervals` to enable GC-bias correction during denoising. This flag is intentionally omitted in this pipeline because the GRCh38 `.dict` file uses alphabetical chromosome ordering (chr1, chr10, chr11, ...) while the BAM headers (and embedded HDF5 dictionary) use numeric ordering (chr1, chr2, chr3, ...). GATK's sequence dictionary comparison fails when the orderings differ, causing the PON build to abort. At ≥30× WGS coverage, GC-bias effects on copy-number calling are modest, making this an acceptable trade-off.

**Alternative PON approaches:**

| Tool | Method | Notes |
|---|---|---|
| CNVKit reference | Median read-depth per target bin | Simpler; no GC correction by default; better for targeted panels |
| WisecondorX | z-score normalisation across bins | Designed for low-pass WGS; overkill at 30× |
| DRAGEN CNV PON | Proprietary normalisation | Requires DRAGEN hardware |

---

## 5. M2: SV Calling

**Subworkflow:** `subworkflows/sv_calling.nf`

Three SV callers run in parallel on the deduplicated BAM, covering complementary algorithmic approaches and SV types. ExpansionHunter runs in parallel for STR genotyping. All three structural callers must complete successfully before Jasmine merge (inner join, fail-fast semantics).

### 4.1 Manta

| Aspect | Detail |
|---|---|
| Module | `modules/manta/call.nf` |
| SV types detected | DEL, DUP, INV, INS, BND (translocation breakends) |
| Algorithm | Split-read + paired-end read evidence; local assembly at candidate breakpoints |
| Strengths | High sensitivity for deletions and inversions ≥50 bp; fast (graph-based candidate generation); low false-positive rate at high coverage |
| Limitations | Lower sensitivity for insertions >500 bp; requires ≥20× coverage; does not call tandem repeats |

Manta (Illumina) uses a two-stage approach: candidate SV generation from anomalous read pairs and soft-clipped reads, followed by local de novo assembly at each candidate. It was chosen because it provides the best balance of sensitivity, specificity, and speed among short-read SV callers, and is the most widely validated tool in clinical and research WGS pipelines (e.g., nf-core/sarek).

**Alternatives:** LUMPY (now superseded), SVaba (slower assembly-based), PBSV (PacBio long-read only), Sniffles2 (Oxford Nanopore only).

### 5.2 Delly

| Aspect | Detail |
|---|---|
| Module | `modules/delly/call.nf` |
| SV types detected | DEL, DUP, INV, BND, INS |
| Algorithm | Paired-end read orientation and insert-size analysis; split-read refinement |
| Strengths | Strong paired-end signal for translocations and inter-chromosomal events; complements Manta's assembly evidence |
| Limitations | Lower sensitivity for insertions; slower than Manta on WGS |
| Output format | Five per-type VCFs merged with shell + bgzip + tabix (bcftools is not available in the Delly biocontainer) |

Delly was chosen to provide independent paired-end evidence that can confirm or refute Manta calls, particularly for inter-chromosomal translocations (BND/TRA) where Manta is less sensitive.

**Implementation note:** The Delly biocontainer (`quay.io/biocontainers/delly:1.2.6`) includes `bgzip` and `tabix` but not `bcftools`. The module calls `delly call -o VCF` for each of the five SV types (DEL INS INV DUP BND) separately, then merges the VCF bodies with a shell `grep` + `sort -k1,1V -k2,2n` pipeline before compressing with bgzip and indexing with tabix.

### 5.3 GRIDSS

| Aspect | Detail |
|---|---|
| Module | `modules/gridss/call.nf` |
| SV types detected | All SV types including small insertions and complex rearrangements |
| Algorithm | Assembly-based; builds a de Bruijn graph over soft-clipped and discordant reads; highly sensitive for sequence-resolved breakpoints |
| Strengths | Best sensitivity for complex SVs and mobile element insertions; sequence-level resolution of breakpoints |
| Limitations | High computational requirements (resource label `process_gridss` — the highest tier); higher false-positive rate; slower than Manta/Delly |
| Resource label | `process_gridss` (dedicated highest tier with OOM retry) |

GRIDSS provides the assembly-based layer that split-read callers miss, particularly for insertions and complex rearrangements.

**GRIDSS_SETUP pre-build optimisation:** A `GRIDSS_SETUP` process (`modules/gridss/setup.nf`) runs `gridss --steps setupreference` once per pipeline invocation to build the BWA index (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`) for the reference FASTA. The index is cached in `storeDir` (`${params.outdir}/cache/gridss_ref`) so it survives across runs. The five index files are staged into each `GRIDSS_CALL` task; because GRIDSS detects the pre-existing `.bwt` file, it skips the ~40-minute per-sample BWA index rebuild, saving approximately 40 × N sample minutes of critical-path time. The `.gridsscache` and `.img` files used in older GRIDSS versions are **not** produced by `setupreference` in GRIDSS 2.13.2 — they are built on demand during each `GRIDSS_CALL` invocation (fast, ~10 s).

**Comparison of M2 callers:**

| Tool | Approach | Speed | Best for |
|---|---|---|---|
| Manta | Split-read + local assembly | Fast | DEL, INV, BND — high precision |
| Delly | Paired-end + split-read | Medium | BND, inter-chromosomal events |
| GRIDSS | Full assembly (de Bruijn graph) | Slow | Complex SVs, insertions, mobile elements |

### 5.4 Jasmine Merge

| Aspect | Detail |
|---|---|
| Module | `modules/jasmine/merge.nf` |
| Input | List of three per-caller VCFs per sample (Manta, Delly, GRIDSS) |
| Output | Merged, deduplicated VCF with caller-support annotations |
| Merge criterion | `min_support=2` (SV must be supported by at least 2 of 3 callers) |
| Join semantics | Inner join — all three callers must succeed for a sample to reach this step |

Jasmine (Jackpot-Aware SV Merger) clusters SVs from multiple callers by breakpoint proximity and reciprocal overlap, then generates a consensus call with caller-support annotations. Setting `min_support=2` retains calls seen by at least two callers, reducing single-caller false positives while preserving sensitivity.

The inner-join channel pattern (`MANTA_CALL.out.vcf.join(DELLY_CALL.out.vcf).join(GRIDSS_CALL.out.vcf)`) means that if any single caller fails for a sample, that sample drops out of the merge channel and the pipeline fails fast — preventing silent partial calls from propagating to downstream steps.

### 5.5 ExpansionHunter (STRs)

| Aspect | Detail |
|---|---|
| Module | `modules/expansionhunter/call.nf` |
| Input | BAM + BAI, reference FASTA, repeat catalog JSON (`assets/eh_catalog.json`) |
| Output | Per-sample VCF with genotyped repeat allele sizes |
| Purpose | Targeted genotyping of known short tandem repeat loci associated with repeat expansion diseases |

ExpansionHunter genotypes only loci defined in the repeat catalog (a curated JSON specifying repeat unit, reference coordinates, and expected allele size range for each disease locus). It does not perform de novo STR discovery. The STR VCF is displayed as ring 4 of the circos plot and as a dedicated section in the HTML report.

STR results are independent of the SV merge — ExpansionHunter runs in parallel and its output is passed directly to the report without merging with SV calls.

---

## 6. M3: CNV Calling

**Subworkflow:** `subworkflows/cnv_calling.nf`

Two complementary CNV callers run in parallel. Their outputs are merged by a custom Python script using reciprocal overlap, producing a consensus BED with confidence labels.

### 5.1 CNVpytor

| Aspect | Detail |
|---|---|
| Module | `modules/cnvpytor/call.nf` |
| Input | BAM + reference FASTA |
| Output | TSV of called CNV segments |
| Method | Read-depth analysis in variable-size bins; Gaussian HMM segmentation; GC and mappability correction |
| SV types | DEL, DUP |

CNVpytor (successor to CNVnator) computes read depth across the genome in adaptive bins, corrects for GC content and mappability biases, and then applies a hidden Markov model to segment the depth profile into integer copy-number states. It works without a panel of normals, making it suitable as a standalone caller and as a complement to the PON-based GATK approach.

### 5.2 GATK gCNV (case mode)

| Aspect | Detail |
|---|---|
| Module | `modules/gatk/gcnv_call.nf` |
| Processes used | `PreprocessIntervals` (bin-length 1000) → `CollectReadCounts` → `DenoiseReadCounts` → `ModelSegments` → `CallCopyRatioSegments` |
| Input | BAM + PON HDF5 + preprocessed interval list |
| Output | SEG file with per-segment copy-number calls and quality scores |
| Method | PON-based denoising followed by CBS (circular binary segmentation) |

GATK gCNV case mode first denoises the sample's read counts by projecting out the systematic noise components learned from the PON, then runs CBS to detect segments of altered copy number. Quality scores per segment reflect the statistical confidence of the call.

The interval list used at case calling time must match the one used during PON construction (both use `--bin-length 1000 --interval-merging-rule OVERLAPPING_ONLY`). This is enforced by routing the output of `GATK_PREPROCESS_INTERVALS` as input to both `CollectReadCounts` (case mode) and `GATK_GCNV_CALL`.

### 5.3 CNV Consensus Merge

**Script:** `bin/cnv_consensus.py`

The consensus script takes one TSV from CNVpytor and one SEG file from GATK gCNV and produces a single BED file per sample. Its merge logic:

1. For each CNVpytor segment, find the best-overlapping GATK segment on the same chromosome with the same SV type (DEL or DUP).
2. If the reciprocal overlap (overlap / min(len_A, len_B)) is ≥0.5, emit a `BOTH / HIGH` consensus call using the GATK copy number (higher confidence from PON denoising).
3. CNVpytor segments with no GATK match are emitted as `CNVpytor / LOW`.
4. Unmatched GATK segments with quality ≥30 are emitted as `GATK_gCNV / MEDIUM`.

**Output BED columns:** chrom, start, end, copy_number, svtype, caller_support, confidence, quality

**Why two callers:**

| Aspect | CNVpytor | GATK gCNV |
|---|---|---|
| Requires PON | No | Yes |
| GC correction | Yes (built-in) | Yes (if annotated intervals used) |
| Sensitivity for small CNVs | Moderate | Higher (1 kb bins) |
| False positive control | Moderate | Better (PON denoising) |
| Copy number resolution | Continuous depth | Integer call + quality score |

The two methods are complementary: CNVpytor provides PON-independent evidence and GC correction; GATK gCNV provides statistically rigorous PON-based denoising. Their overlap (reciprocal overlap ≥0.5, same type) defines high-confidence calls.

**Alternatives:**

| Tool | Method | Notes |
|---|---|---|
| CNVKit | Binned depth + PON | Strong for targeted panels; WGS support is less mature |
| WisecondorX | Within-sample z-scores | Designed for low-pass (1–5×); not appropriate at 30× |
| DRAGEN CNV | Proprietary | Requires DRAGEN hardware; closed source |
| Lumpy (read-depth component) | Read-depth + split-read | Superseded; no longer maintained |

---

## 7. M4: SMN Copy Number Calling

**Subworkflow:** `subworkflows/smn_calling.nf`  
**Module:** `modules/smn_caller/call.nf`  
**Tool:** SMNCopyNumberCaller

### 6.1 The SMN1/SMN2 Problem

SMN1 and SMN2 are highly homologous paralogs located 500 kb apart on chromosome 5q13.2. They share >99% sequence identity across most of their length, differing at only a handful of positions (most critically c.840C in SMN1, which is c.840T in SMN2). Standard short-read aligners cannot reliably distinguish reads originating from one paralog versus the other, making copy-number calling with standard tools (GATK gCNV, CNVpytor) unreliable at this locus.

SMN1 copy number is the molecular basis of spinal muscular atrophy (SMA): two functional SMN1 copies are normal, one copy indicates a carrier, and zero copies (with variable SMN2 CN) indicates SMA with severity inversely correlated with SMN2 copy number.

### 6.2 SMNCopyNumberCaller

SMNCopyNumberCaller is a purpose-built tool developed by Illumina that:

1. Collects read depth and paralog-discriminating variant (PDV) ratios at SMN1/SMN2-specific positions from the BAM.
2. Uses a hidden Markov model trained on WGS data to jointly estimate SMN1 and SMN2 copy numbers and per-allele breakdowns.
3. Detects the 2+0 haplotype (two copies of SMN1 on one chromosome, zero on the other — total CN=2 but carrier status) which would be missed by simple copy-number counting.
4. Reports confidence levels per call.

**Output TSV fields parsed by `bin/smn_report.py`:** `SMN1_CN`, `SMN2_CN`, `SMN1_allele1`, `SMN1_allele2`, `Confidence`

**Classification logic (from `smn_report.py`):**

| SMN1 CN | SMA status | Badge |
|---|---|---|
| 0 | Affected | Red — SMA severity predicted from SMN2 CN |
| 1 | Carrier | Yellow |
| ≥2 | Normal | Green |
| 2 (2+0 haplotype) | Carrier (hidden) | Yellow — flagged with warning |

**Validation truth values for GIAB SMA samples:**

| Sample | SMN1 CN | SMN2 CN | Exon 7 (SMN1/SMN2) | Exon 8 (SMN1/SMN2) | Clinical status |
|---|---|---|---|---|---|
| SMAPB | 0 | 3 | 0 / 3 | 0 / 3 | Affected (SMA) — homozygous SMN1 deletion |
| SMAM | 1 | 5 | 1 / 5 | 1 / 5 | Carrier (SMA mother) — 5× SMN2 |
| SMAD | 1 | 1 | 1 / 1 | 1 / 1 | Carrier (SMA father) — 1× SMN2 |

Note: SMAM and SMAD labels were verified from clinical records in May 2026 (previous documentation had them transposed). The distinguishing feature is SMN2 copy number: SMAM has 5× SMN2 (consistent with maternal carrier of severe SMA family), SMAD has only 1× SMN2.

### 6.3 Limitations and Alternatives

- **Input requirement:** Requires standard short-read WGS BAM aligned to GRCh38 including the SMN region (must not be masked). Works poorly with minimap2-aligned BAMs due to differences in multi-mapping handling at the paralog region.
- **Coverage requirement:** Performance degrades below 30× at the SMN locus.
- **Long-read superiority:** Long-read platforms (PacBio HiFi, Oxford Nanopore) can phase SMN1 and SMN2 alleles directly using tools like paraphase (Pacific Biosciences), which gives near-perfect accuracy including 2+0 detection. For clinical SMA testing, long-read sequencing is increasingly preferred.
- **No de novo pathogenic variant detection:** SMNCopyNumberCaller reports copy number only; it does not call the c.840C>T or other SMN2-to-SMN1 conversion variants.

---

## 8. M5: Annotation

**Subworkflow:** `subworkflows/annotate.nf`

### 7.1 AnnotSV

| Aspect | Detail |
|---|---|
| Module | `modules/annotsv/annotate.nf` |
| Input | Merged SV VCF from Jasmine |
| Output | TSV with one row per SV per overlapping gene annotation |
| Database | User-supplied `--annotsv_db` directory (AnnotSV 3.x format) |
| Container | `quay.io/biocontainers/annotsv:3.4.6--py313hdfd78af_0` |

AnnotSV annotates each SV with:

- **Gene overlaps:** full vs. partial overlap with RefSeq/Ensembl genes; split/full annotation modes
- **gnomAD SV population frequencies:** allele frequency in gnomAD-SV v2.1 (12,920 unrelated genomes)
- **ClinVar pathogenic SVs:** known pathogenic structural variants
- **DGV (Database of Genomic Variants):** known benign population SVs
- **OMIM:** disease associations for overlapping genes
- **ACMG/ClinGen SV classification scores:** semi-automated pathogenicity scoring

If `--annotsv_db` is not provided, ANNOTSV emits a stub empty TSV (with correct column headers) so that downstream steps do not fail.

**Why AnnotSV over alternatives:**

| Tool | Focus | Notes |
|---|---|---|
| AnnotSV | Clinical-grade SV annotation; ACMG scoring | Most comprehensive SV-specific annotation; actively maintained |
| VEP (SV mode) | SNV-centric; SV support is limited | Does not compute SV-specific population frequencies or ACMG class |
| SVAnnotation | Lightweight SV overlapper | Less comprehensive databases; no ACMG scoring |
| AnnotSV + ClassifyCNV | AnnotSV extended | ClassifyCNV adds ACMG CNV classification — not yet integrated here |

### 7.2 gnomAD SV Filter

After AnnotSV, a `GNOMAD_SV_FILTER` process filters the TSV to remove SVs with gnomAD-SV allele frequency > 0.01 (1%). This removes common benign variants from the clinical output. The 1% threshold is the standard population-genetics cutoff for rare variant analysis and aligns with ACMG rare variant guidelines. The filtered TSV is emitted as the final annotation output.

---

## 9. M6/M7: Reporting

**Subworkflow:** `subworkflows/report.nf`

The reporting stage assembles all per-sample results into a self-contained HTML report. It runs after M5 (annotation) but is independent of PON and SMN calling latency because Nextflow joins all channels before `BUILD_HTML_REPORT`.

### 8.1 Circos Plot (pycirclize)

**Script:** `bin/circos_plot.py`  
**Module:** `modules/pycirclize/plot.nf`

Generates a genome-wide SVG circos plot with five concentric rings:

| Ring | Content | Colour |
|---|---|---|
| 1 (outermost) | Chromosome ideograms (chr1–22, X, Y) with labels | Per-chromosome palette |
| 2 | CNV gains (DUP) from consensus BED | Red (#D62728) |
| 3 | CNV losses (DEL) from consensus BED | Blue (#1F77B4) |
| 4 | STR expansion loci from ExpansionHunter VCF | Brown (#8C564B) |
| 5 (chr5 only) | SMN locus highlight (chr5:70,924,941–70,953,015) | Gold (#FFBF00) |

Structural variant breakpoint links are drawn as arcs between genomic coordinates, coloured by SV type (DEL blue, DUP red, INV purple, BND/TRA orange, INS green). The SVG is embedded inline in the HTML report (no external file dependency).

### 8.2 MultiQC

**Module:** `modules/multiqc/report.nf`

MultiQC aggregates FastQC ZIP archives, Picard duplicate metrics, and Mosdepth summary files from all samples into a single multi-sample QC report. This provides a batch-level view of coverage uniformity, duplication rates, and read quality across the cohort.

### 8.3 Truvari Benchmarking (optional)

**Module:** `modules/truvari/bench.nf`  
**Activated by:** `--giab_truth <vcf.gz>`

When a GIAB truth VCF is provided, Truvari benchmarks the Jasmine-merged SV VCF against the truth set. It produces two JSON outputs:

1. **Overall benchmark JSON:** Precision, recall, and F1 across all SV types (parsed by `parse_benchmark()` in `html_report.py`).
2. **Size-bin benchmark JSON:** Precision, recall, and F1 broken down across four size categories:
   - 50–300 bp (small SVs — most challenging)
   - 300 bp–1 kb
   - 1–10 kb
   - >10 kb (large SVs — typically highest recall)

Both JSONs are wired to the HTML report, presenting separate benchmark tables for overall and size-stratified performance.

### 8.4 HTML Report Assembly

**Script:** `bin/html_report.py`  
**Process:** `BUILD_HTML_REPORT` in `report.nf`

The report joins nine input channels per sample using Nextflow's `remainder: true` join with `?: file("NO_FILE")` fallback for optional inputs:

| Input | Source | HTML section |
|---|---|---|
| `sv_tsv` | AnnotSV filtered TSV | SV summary table + top SVs |
| `cnv_bed` | CNV consensus BED | CNV section |
| `smn_tsv` | SMNCopyNumberCaller TSV | SMN section (rendered by `smn_report.py`) |
| `circos_svg` | pycirclize SVG | Genome-wide circos (inline SVG) |
| `benchmark_json` | Truvari overall JSON | Benchmark table (if truth provided) |
| `sizebin_json` | Truvari size-bin JSON | Size-bin benchmark table |
| `coverage_summary` | Mosdepth summary | QC section: depth per chromosome |
| `picard_metrics` | Picard MarkDup metrics | QC section: duplication rate |
| `str_vcf` | ExpansionHunter VCF | STR loci section |

The `smn_report.py` script runs first within the `BUILD_HTML_REPORT` process to generate an HTML fragment for the SMN section, which is then embedded into the full report by `html_report.py` using a Jinja2 template.

The final per-sample `<sample>.report.html` is a self-contained file with all content (including the circos SVG) inlined — no external assets required for viewing.

---

## 10. Known Limitations

### 9.1 SMA Samples Not Yet in Validation Samplesheet

The clinical SMA trio (SMAPB, SMAM, SMAD) have FASTQs located in `ValidationBAM/SMA_BAM/` but are not yet added to `validation/validation_samplesheet.csv`. Once added, they will be aligned through BWA-MEM2 (preferred over external minimap2 BAMs for SMNCopyNumberCaller accuracy — minimap2 multi-mapping handling at the SMN paralog region differs from BWA-MEM2 and may reduce calling accuracy). See TODO §Next Steps item 1.

### 9.2 PON Built Without GC Correction

`CreateReadCountPanelOfNormals` omits `--annotated-intervals` (which would enable GC-bias correction during denoising) because the GRCh38 `.dict` file uses alphabetical chromosome ordering while BAM headers use numeric ordering. GATK's sequence dictionary comparison fails on this mismatch. At ≥30× WGS, GC-bias effects on CNV calls are moderate but non-zero; GC correction would improve sensitivity for CNVs in high-GC regions (e.g., pericentromeric repeats).

**Workaround options:** (1) Re-sort the `.dict` file to numeric chromosome order to match BAM headers; (2) use a FASTA/dict generated by the same tool that created the BAM headers.

### 9.3 Nextflow Channel Exhaustion

In Nextflow DSL2, a `Channel.fromPath()` channel (queue channel) is consumed after its first use. Reference files shared across multiple subworkflows (FASTA, FAI, dict, intervals) must be declared as `Channel.value()` (broadcast value channels) in `main.nf` so that PREPROCESS, SV_CALLING, CNV_CALLING, and SMN_CALLING can each receive the same reference file without the channel being exhausted after the first sample or subworkflow consumes it. This pattern is applied in `main.nf` and `pon_build.nf`. Forgetting this — using `Channel.fromPath()` for shared references — causes silent processing of only the first sample in a multi-sample run.

### 9.4 AnnotSV Not in Conda Environment

AnnotSV requires complex Perl dependencies and is not available in the conda environment. It runs exclusively via Docker (`quay.io/biocontainers/annotsv:3.4.6--py313hdfd78af_0`). When `--annotsv_db` is not provided, the ANNOTSV process emits an empty TSV (correct column headers, no data rows) so downstream processes do not fail.

### 9.5 samtools flagstat Not Wired

The mapping rate percentage (`mapped_pct`) shows "N/A" in the HTML report QC section because `samtools flagstat` is not run as a pipeline step. Coverage depth comes from mosdepth and duplication rate from Picard. Adding a SAMTOOLS_FLAGSTAT module and wiring its output to the report would resolve this.

### 9.6 Jasmine Inner-Join Fail-Fast

If any one of Manta, Delly, or GRIDSS fails for a sample, that sample's SV merge channel is silently empty (Nextflow inner join semantics) and the sample will not appear in any downstream SV output. This is intentional (prevents partial calls from propagating) but means that transient caller failures require investigation of the work directory rather than an explicit error message.

### 9.7 GRIDSS Memory Demands

GRIDSS is the most memory-intensive process and uses the dedicated `process_gridss` resource label. Automatic retry on OOM exit codes (137, 143, 104, 134, 139) is configured in `conf/base.config`, but very large SVs or highly repetitive regions can cause GRIDSS to exceed even the retry memory allocation on machines with limited RAM.

---

## 11. Alternative Pipeline Strategies

### 10.1 DRAGEN End-to-End

Illumina DRAGEN provides alignment, duplicate marking, SV calling (Manta-like), CNV calling, and STR calling in a single hardware-accelerated pipeline on FPGA hardware. Advantages: fastest available; single vendor support. Disadvantages: requires expensive DRAGEN hardware or cloud-based DRAGEN licensing; closed-source algorithms; limited customisability; no equivalent to the SMNCopyNumberCaller integration. DRAGEN is appropriate when throughput is the primary concern and institutional hardware is available.

### 10.2 nf-core/sarek

The nf-core/sarek pipeline (https://nf-co.re/sarek) is a community-maintained Nextflow DSL2 pipeline for germline and somatic variant calling. It includes BWA-MEM2 alignment, Picard/GATK preprocessing, and SV calling with Manta and Tiddit. Advantages: well-tested, nf-core community support, CI/CD, schema validation. Disadvantages: does not integrate CNVpytor, SMNCopyNumberCaller, ExpansionHunter, or AnnotSV in the same workflow; SV benchmarking and HTML reporting are not included. SVcaller was built to add these capabilities in an integrated, report-centric workflow.

### 10.3 Long-Read Approaches (ONT + Sniffles2; PacBio HiFi + PBSV/paraphase)

| Aspect | Short-read (this pipeline) | Long-read (ONT/HiFi) |
|---|---|---|
| SV detection | ≥50 bp deletions/duplications best; insertions harder | All SV types including mobile element insertions; better for inversions and complex rearrangements |
| STR genotyping | ExpansionHunter (catalog-based) | TRGT, Straglr (de novo expansion discovery possible) |
| SMN1/SMN2 | SMNCopyNumberCaller (PDV-based) | paraphase, NanoVar (direct phasing of alleles — near-perfect accuracy) |
| Coverage required | ≥30× short-read | 15–30× long-read (lower for structural variants) |
| Cost | Lower (short-read sequencing) | Higher (long-read sequencing) |
| Repeat regions | Poor (reads shorter than repeat) | Excellent (reads span entire repeats) |

For clinical SMN1 testing, rare disease SV diagnosis, or centromeric/telomeric regions, long-read sequencing with Sniffles2 (ONT) or PBSV + paraphase (PacBio HiFi) is increasingly preferred. SVcaller remains the appropriate choice for high-throughput short-read WGS cohorts where cost and throughput are priorities.

### 10.4 When to Use This Pipeline

SVcaller is appropriate when:
- Short-read Illumina WGS at ≥30× is the sequencing modality.
- Integrated SV + CNV + STR + SMN calling from a single BAM is required.
- Per-sample HTML reports with circos visualisation are needed.
- Optional GIAB-based benchmarking with size-stratified metrics is useful.
- A pre-built GATK gCNV PON from GIAB normals is available (or one has been built from institutional controls).

---

## 12. Key Parameters Reference

| Parameter | Default | Description |
|---|---|---|
| `--input` | required | Samplesheet CSV (`sample,fastq_1,fastq_2,bam`) |
| `--ref_fasta` | required | GRCh38 FASTA; `.fai` and `.dict` inferred from path |
| `--pon` | null | GATK gCNV PON HDF5 file; sentinel `NO_PON` used if absent |
| `--intervals` | null | Target BED (WGS autosomal intervals); required for GATK gCNV |
| `--annotsv_db` | null | AnnotSV database directory; annotation skipped if absent |
| `--giab_truth` | null | GIAB truth VCF.gz; enables Truvari benchmarking |
| `--eh_catalog` | required | ExpansionHunter repeat catalog JSON |
| `--min_depth` | 30 | Minimum mean coverage; pipeline halts below this threshold |
| `--outdir` | `results` | Output directory |
| `--utils_container` | `svcaller/utils:1.0` | Docker image for Python bin/ scripts |
| `-profile docker` | — | Use Docker for all containers (recommended; all tags verified) |
| `-resume` | — | Resume from last successful checkpoint |

**Resource labels (defined in `conf/base.config`):**

| Label | Tier | Typical use |
|---|---|---|
| `process_single` | Minimal CPU/memory | Python scripts, consensus merge |
| `process_low` | Low | FastQC, light processing |
| `process_medium` | Medium | Samtools sort, GATK CollectReadCounts |
| `process_high` | High | BWA-MEM2, CreatePON |
| `process_gridss` | Maximum | GRIDSS (dedicated highest tier) |

All resource labels support automatic retry on OOM exit codes: 137, 143, 104, 134, 139.

---

*Documentation generated from source code at `/data/alvin/SVcaller/` (updated 2026-05-23: SMN truth table corrected, DELLY VCF output, GRIDSS pre-build, Known Limitations revised).*

---

## Glossary of Tools

### Pipeline tools (used in this pipeline)

| Tool | Role in pipeline | GitHub |
|---|---|---|
| **Nextflow** | Workflow orchestration (DSL2) | [nextflow-io/nextflow](https://github.com/nextflow-io/nextflow) |
| **FastQC** | Raw read quality control | [s-andrews/FastQC](https://github.com/s-andrews/FastQC) |
| **BWA-MEM2** | Short-read alignment (GRCh38) | [bwa-mem2/bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) |
| **SAMtools** | BAM sorting and indexing | [samtools/samtools](https://github.com/samtools/samtools) |
| **Picard MarkDuplicates** | PCR duplicate marking | [broadinstitute/picard](https://github.com/broadinstitute/picard) |
| **Mosdepth** | Coverage QC and depth gating | [brentp/mosdepth](https://github.com/brentp/mosdepth) |
| **Manta** | SV calling (split-read + paired-end) | [Illumina/manta](https://github.com/Illumina/manta) |
| **Delly** | SV calling (paired-end + split-read) | [dellytools/delly](https://github.com/dellytools/delly) |
| **GRIDSS** | SV calling (assembly-based breakend) | [PapenfussLab/gridss](https://github.com/PapenfussLab/gridss) |
| **Jasmine** | Multi-caller SV merging | [mkirsche/Jasmine](https://github.com/mkirsche/Jasmine) |
| **ExpansionHunter** | Short tandem repeat (STR) genotyping | [Illumina/ExpansionHunter](https://github.com/Illumina/ExpansionHunter) |
| **CNVpytor** | Read-depth CNV calling | [abyzovlab/CNVpytor](https://github.com/abyzovlab/CNVpytor) |
| **GATK** | gCNV PON build + case calling | [broadinstitute/gatk](https://github.com/broadinstitute/gatk) |
| **SMNCopyNumberCaller** | SMN1/SMN2 paralog copy number | [Illumina/SMNCopyNumberCaller](https://github.com/Illumina/SMNCopyNumberCaller) |
| **AnnotSV** | Clinical SV annotation | [lgmgeo/AnnotSV](https://github.com/lgmgeo/AnnotSV) |
| **Truvari** | SV benchmarking against truth sets | [ACEnglish/truvari](https://github.com/ACEnglish/truvari) |
| **pycirclize** | Circos-style genome visualization | [moshi4/pyCirclize](https://github.com/moshi4/pyCirclize) |
| **MultiQC** | Aggregate QC report | [MultiQC/MultiQC](https://github.com/MultiQC/MultiQC) |

### Alternative tools (referenced in comparisons)

| Tool | Category | GitHub |
|---|---|---|
| **BWA-MEM** | Short-read aligner (predecessor to BWA-MEM2) | [lh3/bwa](https://github.com/lh3/bwa) |
| **Bowtie2** | Short-read aligner | [BenLangmead/bowtie2](https://github.com/BenLangmead/bowtie2) |
| **minimap2** | Long-read / splice-aware aligner | [lh3/minimap2](https://github.com/lh3/minimap2) |
| **LUMPY** | SV caller (probabilistic framework) | [arq5x/lumpy-sv](https://github.com/arq5x/lumpy-sv) |
| **SVaba** | SV caller (local assembly) | [walaj/svaba](https://github.com/walaj/svaba) |
| **PBSV** | SV caller for PacBio long reads | [PacificBiosciences/pbsv](https://github.com/PacificBiosciences/pbsv) |
| **Sniffles2** | SV caller for Oxford Nanopore reads | [fritzsedlazeck/Sniffles](https://github.com/fritzsedlazeck/Sniffles) |
| **CNVKit** | Read-depth CNV calling (capture / WGS) | [etal/cnvkit](https://github.com/etal/cnvkit) |
| **WisecondorX** | PON-based CNV calling (low-coverage) | [CenterForMedicalGeneticsGhent/WisecondorX](https://github.com/CenterForMedicalGeneticsGhent/WisecondorX) |
| **Paraphase** | SMN1/SMN2 phased copy number (PacBio) | [PacificBiosciences/paraphase](https://github.com/PacificBiosciences/paraphase) |
| **nf-core/sarek** | End-to-end Nextflow WGS pipeline | [nf-core/sarek](https://github.com/nf-core/sarek) |
