# SVcaller Design Decisions

Why the pipeline is built the way it is. Each section covers a non-obvious choice, the problem it solves, and what was traded off.

## Channel.value() for shared reference files

**The problem.** Nextflow queue channels are single-consumption: once a process reads from a channel, it is exhausted. When multiple subworkflows (M2, M3, M4, M5) all need the reference FASTA and its index, using `Channel.fromPath()` means only the first subworkflow gets the file. Subsequent subworkflows receive an empty channel and silently produce no output.

**The approach.** All shared reference files (FASTA, FAI, dict, BWA-MEM2 index directory, cytobands, EH catalog) are published as `Channel.value()` in `main.nf`. Value channels are broadcast: every consumer gets a copy.

**Trade-off.** Value channels cannot be filtered or transformed after creation. Any per-sample branching on reference content (e.g., a different reference per sample) requires a different pattern.

---

## Sentinel files for optional inputs

**The problem.** Nextflow process `input:` declarations require every slot to be filled. Optional inputs (PON, AnnotSV DB, truth VCF, STR VCF) may or may not exist depending on the run. Conditional process definitions would require duplicating the entire process block.

**The approach.** Optional inputs use sentinel file paths: `file("NO_PON")`, `file("NO_FILE")`, `file("NO_STR")`. The process script checks `param_name.name != "NO_FILE"` before constructing CLI arguments. This keeps topology uniform — the DAG shape is identical whether or not optional inputs are provided.

**Trade-off.** Sentinel strings must be checked consistently in every process and Python script that handles optional inputs. A mismatch (e.g., checking `"NO_FILE"` but receiving `"null"`) silently skips a section. See `_parse_mosdepth`, `_parse_flagstat` etc. in `bin/html_report.py` for the canonical check pattern.

---

## Inner join on Jasmine (fail-fast on any caller)

**The problem.** If one SV caller produces no output (crash, empty VCF, zero variants), Jasmine can still run with the remaining callers. The merged VCF would then have incorrect SUPP_VEC values — a call appearing in only Manta and Delly would show `SUPP_VEC=110` when GRIDSS simply did not run, not because GRIDSS rejected the call.

**The approach.** The four caller output channels are joined with an inner join before Jasmine. If any caller channel produces no tuple for a sample, the join produces nothing, and Jasmine never runs for that sample. The pipeline fails visibly rather than silently producing corrupt SUPP_VEC data.

**Trade-off.** A single caller crash fails the entire sample. For clinical pipelines, a known-bad SUPP_VEC is worse than a failed run — operators can re-run, but a clinician cannot detect a silently wrong support vector.

---

## No GC correction in the GATK gCNV PON

**The problem.** GATK `CreateReadCountPanelOfNormals` accepts `--annotated-intervals` for GC-corrected denoising. This requires a `PreprocessIntervals` output whose `.dict` header lists chromosomes in the same order as the BAM headers. The GRCh38 `.dict` created by `samtools dict` sorts chromosomes alphabetically (`chr1, chr10, chr11...`), while BAM headers from typical aligners use numeric order (`chr1, chr2, chr3...`). GATK's dict comparison fails with a sequence dictionary mismatch error.

**The approach.** The PON was built without `--annotated-intervals`. GC correction is omitted. At ≥30× WGS coverage, GC bias is small relative to true CNV signal, and the PON built from 7 GIAB samples provides adequate noise modelling.

**Trade-off.** CNV calls in high-GC or low-GC regions (centromeres, telomeres, GC-rich genes) will have higher false-positive rates than a GC-corrected model. For targeted panels or low-coverage WGS this trade-off would be unacceptable. At 30× WGS it is acceptable.

---

## Canonical reference to skip FILTER_CHROMS

**The problem.** Standard GRCh38 FASTA files contain ALT contigs, decoy sequences, HLA variants, and unplaced scaffolds. Aligners place reads on these contigs, and SV callers emit calls on `chrUn_*`, `HLA-*`, and `*_alt` contigs. Downstream tools (Jasmine, Truvari, GRIDSS) have inconsistent handling of non-canonical chromosomes and may crash or produce unexpected output.

**The approach.** FILTER_CHROMS (`modules/samtools/filter_chroms.nf`) strips non-canonical contigs from BAM files. For FASTQ inputs, aligning directly to `hg38.canonical.fa` (chr1-22+X+Y+M only) means the aligned BAM never contains non-canonical contigs, and FILTER_CHROMS can be skipped entirely (saves ~70 min/sample on a ~30× BAM — measured; see CLAUDE.md "FILTER_CHROMS is slow"). Pre-aligned BAM inputs always run FILTER_CHROMS because their alignment reference is unknown.

**Trade-off.** Variants on unplaced scaffolds and ALT contigs are silently dropped. For germline WGS of canonical chromosomes this is acceptable. For studies targeting HLA diversity or structural variation on ALT contigs, this pipeline is not appropriate.

---

## Canonical @SQ headers in FILTER_CHROMS output BAM

**The problem.** A subtle but fatal failure mode: FILTER_CHROMS kept all 3,366 `@SQ` header lines from the original hg38 BAM even though it filtered the reads down to canonical chromosomes. Manta's scanner found 526K candidate SV edges from the read alignments, but Manta's assembly phase silently produced 0 variants. The root cause: when the BAM header references chromosomes (`@SQ SN:chrUn_*`, `@SQ SN:HLA-*`, etc.) that have no corresponding reads, Manta's assembly internally crashes and emits empty output — no error, no warning, just zero calls.

**The approach.** The awk `@SQ` filter in `filter_chroms.nf` was extended from `if (c in order)` to `if (c in order && c in can)`. The `can` array contains exactly the 25 canonical chromosomes (chr1-22+X+Y+M). Only matching `@SQ` lines are emitted into the output BAM header. All 3,341 non-canonical `@SQ` entries are dropped, even though their corresponding reads were already absent.

**Why this was non-obvious.** The filtered BAM looked healthy: 628M reads, correct MAPQ distribution, normal insert size distribution. The Manta scanner statistics (3.7M anomalous pairs, 1.7M split reads, 526K candidate edges) confirmed valid SV evidence was present. The failure was entirely in the assembly phase, which is not separately logged. The `@SQ` count mismatch was only found by comparing `samtools view -H` output before and after filtering.

**Trade-off.** None: this is a pure fix. BAMs with non-canonical @SQ lines and no corresponding reads carry no information in those headers; dropping them cannot affect variant calling.

---

## Pre-flight reference/BAM validation (`VALIDATE_REF_BAM`)

**The problem.** The `@SQ`-mismatch failure above (and similar reference/BAM inconsistencies) surfaces *late* — Manta is the first caller to crash, ~4 hours into a run, and it crashes silently with zero output. An operator discovers the bad input only after burning a half-day of compute, then has to diagnose a non-erroring failure.

**The approach.** `VALIDATE_REF_BAM` (`modules/samtools/validate_ref_bam.nf`) runs as a pre-flight gate at the very start of M2, before any caller. It compares the chromosome set of the input BAM against the reference FASTA and fails fast — in seconds, with a clear error — if the reference carries contigs the BAM lacks, hinting to use `hg38.canonical.fa` rather than `hg38.fa` for BAM inputs.

**Trade-off.** A few extra seconds at the front of every run, and the check is a superset relationship (it does not catch every conceivable mismatch). In exchange, the most common and most expensive misconfiguration is caught at second zero instead of hour four. For a clinical pipeline where a silent 4-hour failure is the worst outcome, this is strongly positive.

---

## Per-sample work directories

**The problem.** Nextflow stores intermediate files in `work/`. If multiple samples share a work directory and one sample fails, re-running with `-resume` can create session lock conflicts. More critically, a shared `work/` accumulates intermediates across all samples and run types, making targeted cleanup impossible without also deleting cache for other samples.

**The approach.** Each sample run uses its own work directory: `work_<sampleId>` for single samples, `work_<batchName>` for multi-sample batches. Once a sample's results are published to `{outdir}/{sample}/`, its work directory can be deleted safely (`rm -rf work_SAMPLE`) without affecting any other sample's cache.

**Trade-off.** Managing multiple work directories requires more explicit bookkeeping — you cannot `-resume` across work directories. The convention `work_<sampleId>` makes it unambiguous which work directory belongs to which sample.

---

## Scramble MEI canonical SVLEN estimates

**The problem.** Scramble reports mobile element insertions but does not always report the full inserted sequence length. Truvari compares SVs using a ±30% size similarity threshold. An L1 insertion truncated to 1000 bp in the sample would fail the size match against a truth-set L1 with SVLEN=6000, even if it is the same insertion.

**The approach.** At the Jasmine merge step, Scramble VCFs are post-processed with awk to assign canonical SVLEN values: ALU=300 bp, L1=6000 bp, SVA=1500 bp. These represent full-length MEI sizes, not the observed insertion lengths.

**Trade-off.** This inflates SVLEN for truncated L1 insertions, which are common. A truncated L1 at ~1000 bp with canonical SVLEN=6000 may fail Truvari's ±30% size gate against a truth L1 reported at 2000 bp. This is a known source of false negatives in the Truvari benchmark. Using actual MEINFO coordinates for SVLEN estimation is a planned improvement.

---

## Depth scatter subsampling (1-in-10 normal dots)

**The problem.** mosdepth at 50 kb resolution produces 61,775 windows across the canonical genome. Rendering each window as an SVG element produces a 6.8 MB SVG, which inflates the HTML report to ~7 MB and is slow to render in browsers.

**The approach.** Normal-coverage windows (log₂ ratio within [-0.5, +0.3]) are subsampled at 1-in-10 before rendering. Gain and loss windows (the clinical signal) are always rendered at full density.

**Trade-off.** Fine-scale coverage variation in diploid regions is less visible. For CNV detection the relevant signal is the outlier windows, which are unaffected. The SVG drops from 6.8 MB to ~1 MB; the HTML report from ~7 MB to ~2 MB.

---

## Inlined Bootstrap CSS for air-gapped reports

**The problem.** The clinical HTML report styled itself by linking the Bootstrap CSS from a CDN. Clinical and diagnostic networks are frequently air-gapped — no outbound internet. On those machines the `<link>` silently fails and the report renders as unstyled HTML: tables lose their grid, the layout collapses, and a document meant for clinicians looks broken.

**The approach.** `bin/html_report.py` inlines Bootstrap directly into a `<style>` block in the generated HTML. It resolves the CSS in priority order: the bundled `assets/bootstrap.min.css` baked into the `svcaller/utils` container first, a CDN fetch as a fallback for connected build environments, and an empty string as a last resort so the report still assembles. The result is a single self-contained `.report.html` with no external dependencies.

**Trade-off.** Each report carries ~25 KB of inlined CSS instead of sharing one cached CDN file across reports, and the bundled Bootstrap is pinned to one version rather than tracking the CDN. For single-file clinical artifacts that must render identically on a disconnected workstation years later, self-containment outweighs the duplication.

---

## A failed caller must never publish an output (no empty placeholders)

**The problem.** Several modules ended their command blocks with an `|| touch <output>`
fallback, or wrote `echo '{}' > <output>.json`. The intent was benign — satisfy the
Nextflow output declaration so the pipeline keeps moving when a tool's output filename
varies across versions. The effect was not. When a caller genuinely **failed**, the
process still exited 0 and published a **zero-byte file**. Nextflow saw its declared
output, marked the task successful, and the run completed.

On disk, a crashed stage and a negative result are then byte-for-byte indistinguishable.
Downstream, OmniGen gated its reads on `os.path.exists()` — `True` for a 0-byte file —
and rendered a **complete-looking consumer report** reading *"0 Carrier findings,
0 Medication flags, Clear"*. A crashed SMN caller became a clean bill of health for a
human being. (Full write-up: [CHANGES.md](CHANGES.md).)

**The approach.** A caller that fails, fails. No placeholder is written, the process exits
non-zero, and `errorStrategy = 'finish'` (`conf/base.config`) stops the run with an
actionable error — the real stderr, exit code and work dir. The three states are kept
strictly distinct:

| State | On disk | How it is reached |
|---|---|---|
| Real negative result | Well-formed artifact, header, zero data rows | Caller ran successfully and found nothing |
| Intentionally skipped | Well-formed header-only VCF from a `*_STUB` process | Operator passed `--skip_melt` / `--skip_scramble` / `--skip_svaba` / `--skip_gridss` |
| Failure | **Nothing** | Caller crashed, or its dependencies are missing |

The key design point is the middle row: *optionality is expressed in the workflow, not by
writing a fake file in the module.* The `--skip_*` flags and their `*_STUB` processes
already existed — which is precisely why the in-module "tool not found → emit empty VCF,
exit 0" paths were pure failure-masking and could be deleted outright. "MELT is not
installed" is a misconfiguration to surface, not a genome with no mobile-element
insertions.

**Trade-off.** Runs that previously "succeeded" on a broken or missing caller will now
fail. That is the point: they were not succeeding, they were silently reporting nothing
found. The cost is that a genuinely optional stage must be skipped explicitly rather than
allowed to quietly no-op. For a pipeline whose output a clinician reads, an error you must
act on beats a report you cannot trust.

`tests/test_no_empty_placeholders.py` statically guards the patterns against reintroduction.

---

## Related

- [Architecture reference](reference-architecture.md) — module-by-module technical description
- [Parameter reference](reference-parameters.md) — all CLI flags
- [Changes](CHANGES.md) — behavioural changes and the incidents behind them
