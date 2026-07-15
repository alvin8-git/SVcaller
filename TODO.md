# SVcaller TODO

Tracks implementation status against the design spec (`docs/superpowers/specs/2026-05-09-svcnv-caller-design.md`).

---

## M1 — Preprocessing

- [x] BWA-MEM2 alignment (`modules/bwamem2/align.nf`)
- [x] Picard MarkDuplicates (`modules/picard/markduplicates.nf`)
- [x] mosdepth coverage QC with min-depth halt (`modules/mosdepth/coverage.nf`)
- [x] FastQC raw read QC (`modules/fastqc/qc.nf`)
- [x] samtools flagstat — mapped_pct wired to HTML QC section (`modules/samtools/flagstat.nf`) — 2026-05-16

## M2 — SV Calling

- [x] Manta (`modules/manta/call.nf`)
- [x] DELLY (`modules/delly/call.nf`)
- [x] GRIDSS (`modules/gridss/call.nf`) + BND→DEL/DUP/INV converter (`bin/gridss_convert_bnd.py`)
- [x] Scramble MEI caller (`modules/scramble/call.nf`)
- [x] MELT MEI caller (`modules/melt/call.nf`) — local container build required
- [x] SvABA local-assembly caller (`modules/svaba/call.nf`) — 2026-07-15 fixed latent staging bug: SvABA had NEVER produced a variant because its classic BWA index (`.amb/.ann/.bwt/.pac/.sa`) was never staged next to `ref_fasta` (masked by a since-removed `|| true`). Now declares `path bwa_index`, threaded via `ch_bwa_index` with `--bwa_index` param + fail-loud check. NOTE: HG002/HG003 must be re-run to actually gain SvABA calls.
- [x] STRling genome-wide STR caller (`modules/strling/call.nf`)
- [x] ExpansionHunter — 32 disease STR loci (`modules/expansionhunter/call.nf`)
- [x] JASMINE merge (min_support=2) (`modules/jasmine/merge.nf`)
- [x] SV PON annotation — GIAB 7-sample BED (`modules/annotsv/sv_pon_annotate.nf`)
- [x] sv_tbi emitted for Truvari benchmarking

## M3 — CNV Calling

- [x] CNVpytor (`modules/cnvpytor/call.nf`)
- [x] GATK gCNV case mode (`modules/gatk/gcnv_call.nf`)
- [x] Consensus CNV BED (`bin/cnv_consensus.py`)
- [x] GATK PreprocessIntervals in case mode (`subworkflows/cnv_calling.nf`)
- [x] GATK PreprocessIntervals in PON build (`workflows/pon_build.nf`)

## M4 — SMN1/SMN2

- [x] SMNCopyNumberCaller v1.1 (`modules/smn_caller/call.nf`)
- [x] SMN HTML section (`bin/smn_report.py`)
- [x] 2+0 haplotype detection — verified correct in `bin/smn_report.py`
- [x] SMN truth table with GIAB + SMA trio (`validation/smn_truth_table.tsv`)

## M5 — Annotation

- [x] AnnotSV v3.4 (`modules/annotsv/annotate.nf`)
- [x] gnomAD-SV AF > 1% filter (`modules/annotsv/annotate.nf::GNOMAD_SV_FILTER`)
- [ ] AnnotSV database download — requires `--annotsv_db` parameter; not automated

## M6 — Visualization

- [x] pycirclize Circos SVG + PNG (`bin/circos_plot.py`, `modules/pycirclize/plot.nf`)
- [x] Ring 1: chromosome ideograms
- [x] Ring 2: CNV gains (red `#D62728`)
- [x] Ring 3: CNV losses (blue `#1F77B4`)
- [x] Ring 4: STR expansion markers (brown `#8C564B`)
- [x] Ring 5: SMN locus highlight chr5q13.2 (gold `#FFBF00`)
- [x] Center: SV links (type-coloured)
- [x] PNG at 1200 dpi
- [x] circos.png published to outdir

## M7 — Reporting & Benchmarking

- [x] Jinja2 HTML per-sample report (`bin/html_report.py`)
- [x] MultiQC aggregated QC report (`modules/multiqc/report.nf`)
- [x] Truvari GIAB benchmarking (HG002 SV v0.6) — wired via `--giab_truth`
- [x] HTML report section 2: Alignment QC — mosdepth + Picard + flagstat mapped_pct
- [x] HTML report section 7: STR expansion loci — ExpansionHunter VCF parsing
- [x] HTML report section 8: 3-tier clinical SV classification (Tier 1=ACMG SF v3.2, Tier 2=OMIM morbid, Tier 3=all; top 10 in HTML + full XLS export)
- [x] Truvari per-size-bin metrics — 4 bins wired to HTML report

## M8 — Blood Group & Copy-Number Traits (OmniGen)

- [x] Targeted depth process `TRAIT_DEPTH` (`modules/traits/depth.nf`) + region/control BED (`assets/cnv_trait_regions.bed`)
- [x] Interpreters (`bin/rh_status.py`, `amy1_cn.py`, `gst_null.py`, `lpa_kiv2.py`) sharing `bin/cnv_traits_common.py`
- [x] `CNV_TRAITS` subworkflow + report card in `bin/html_report.py` / `subworkflows/report.nf`
- [x] **Real HG002 depth run + validation — 2026-07-12.** Ran targeted depth (`samtools bedcov -Q 0`, mosdepth-equivalent) on the 30X GRCh38 HG002 BAM; ran all four interpreters. Artifacts: `results/HG002/HG002.trait_depth.regions.bed.gz`, `results/HG002/HG002.bam_stats.json`, `results/HG002/bloodgroup/HG002.rh_status.tsv`, `results/HG002/cnv_traits/HG002.{amy1,gst_null,lpa_kiv2}.tsv`. Control baseline 32.31X ≈ pipeline autosomal ~31.6X (PASS). Rh(D) **positive** (RHD CN 2; consensus BED has 0 RHD DELs — agree). GSTT1 **present** (ratio 1.00). GSTM1 **present** (ratio 0.46 ≈ 1 copy): GIAB v5.0q truth shows the 18.4 kb GSTM1 deletion is **heterozygous** (`GT 0|1`), so HG002 is het-del not homozygous-null — depth and truth agree; see `docs/omnigen-additions-plan.md`. AMY1=2 / KIV-2=5 reported but uncalibrated. Demo report `docs/demo/HG002_report.html` regenerated with real values.
- [x] **Real HG001 depth run + validation — 2026-07-12.** Ran targeted depth (`samtools bedcov -Q 0`, mosdepth-equivalent) on the 30X GRCh38 HG001/NA12878 BAM; ran all four interpreters. Artifacts: `results/HG001/HG001.trait_depth.regions.bed.gz`, `results/HG001/HG001.bam_stats.json`, `results/HG001/bloodgroup/HG001.rh_status.tsv`, `results/HG001/cnv_traits/HG001.{amy1,gst_null,lpa_kiv2}.tsv`. Control baseline 30.81X ≈ pipeline autosomal ~29.9X (PASS). Rh(D) **positive** (RHD depth 19.0X → CN 1, confidence MEDIUM; RHCE cross-mapping — depth drops to 15.8X at MAPQ≥20; consensus BED has 0 RHD DELs — agree). GSTT1 **present** (ratio 1.02). **GSTM1 flagged:** NA12878 is widely cited as the canonical homozygous GSTM1-null anchor, but this BAM reads GSTM1 at **16.8X, ratio 0.545 (~1 copy)**. A **paralog-aware PSV test** (GSTM1 vs its near-identical paralog GSTM2) settles het-deletion vs homozygous-null, which depth alone cannot: (1) **99% of window reads (611/617) are MAPQ≥20** — uniquely placed at GSTM1, not ambiguous GSTM2 cross-mappers (only 6 are MAPQ0); a homozygous null relying on GSTM2 cross-map would render these reads MAPQ0, since GSTM1 exists in the *reference* regardless of the sample. (2) At **39 fixed GSTM1-vs-GSTM2 discriminating positions, 433/433 reads carry the GSTM1-specific allele and 0 the GSTM2 allele** → no detectable cross-mapping. (3) The **11 GSTM1-window `1/1` HaplotypeCaller calls are none on a PSV and none match the GSTM2 base** → genuine GSTM1 hemizygous variants (1 surviving copy ⇒ reads homozygous), not cross-map artifacts. Verdict: **heterozygous GSTM1 deletion (~1 real copy)**; the interpreter's `present` call is correct **for the data**, and this is a **MISMATCH vs the assumed homozygous-null truth** — that premise is not supported by this BAM. Not fabricated. (No paralog-aware fix to `gst_null.py` is needed here: cross-mapping is minimal, so a *true* homozygous GSTM1 null would read ~0 depth and correctly fire `null`; the real gap is the het-del tier, tracked below.) AMY1=2 / KIV-2=6 reported but uncalibrated. Demo report `docs/demo/HG001_report.html` generated with real values.
- [x] **Circos CNV-trait ring + clean HG002/HG001 re-render — 2026-07-12.** Added an optional labelled trait ring (r=59-64) to `bin/circos_plot.py` (`--rh-status/--amy1/--gst-null/--lpa-kiv2`), plotting the 5 loci with the call read from the contract TSVs (present/normal=grey, deletion/null=blue, high-CN=red); ring is skipped when no trait TSVs are supplied (verified: COLO829 renders with the ring omitted). Wired the four trait channels into `modules/pycirclize/plot.nf` + `subworkflows/report.nf`. Re-rendered HG002 and HG001 circos cleanly (fresh 2026-07-12 provenance in SVG metadata + on-plot footer) and regenerated both demo reports (`docs/demo/{HG002,HG001}_report.html`) so they inline the new SVG. NOTE (superseded 2026-07-13 — see next item): the HG002 Truvari-benchmark and STR report sections were absent from this rebuild because their `benchmark.json`/STR inputs lived in `work_HG002` (reclaimed during the results/ consolidation); coverage + flagstat QC were reconstructed from the BAM (`samtools bedcov`/`flagstat`) so the QC section stayed intact. Tests: `test_circos_plot.py`, `test_cnv_display.py`, `test_html_report.py`, `test_cnv_traits.py` — 32 passed.
- [x] **HG002 demo report: Truvari + STR sections restored — 2026-07-13.** Re-ran the two stages whose inputs were lost with `work_HG002`, using the exact module commands against the existing `results/HG002/HG002.sv_merged.vcf.gz` (Manta+DELLY, `--skip_gridss`) and the 30X GRCh38 HG002 BAM, then regenerated the report with the real generator (`bin/smn_report.py` + `bin/html_report.py` in `svcaller/utils:1.2`). No pipeline logic changed.
  - **Truvari v4.3.1** (`quay.io/biocontainers/truvari:4.3.1`, module flags: `--passonly --pick multi --pctseq 0 --typeignore --sizemin 50`), both truth sets:
    - **T2TQ100-V1.0** (`/data/alvin/ref/GIAB/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz`): TP-base 7,571 / TP-comp 7,135 · FP 2,605 · FN 22,120 · **precision 0.733 · recall 0.255 · F1 0.378** (base 29,691 / comp 9,740).
    - **v5.0q** (`/data/alvin/ref/GIAB/HG002_GRCh38_v5.0q_stvar.vcf.gz`): TP-base 7,286 / TP-comp 6,884 · FP 2,450 · FN 20,837 · **precision 0.738 · recall 0.259 · F1 0.383** (base 28,123 / comp 9,334). Recall is short-read/`--skip_gridss` limited (truth is assembly-based and INS-heavy) — these are the real numbers for this call set, not a regression.
  - **ExpansionHunter 5.0.0** on the 32-locus `assets/eh_catalog.json`: 32 loci genotyped. Report calls 2 EXPANDED (ATXN2 22/62 ⚠INREPEAT; NIPA1 4/18) and 2 INTERMEDIATE (CNBP 13/51 ⚠INREPEAT; RFC1 12/33 ⚠INREPEAT) — all short-read INREPEAT estimates, long-read confirmation flagged in the report.
  - Artifacts now in `results/HG002/`: `HG002.str.vcf.gz{,.tbi}`, `HG002.str_profile.json`, `benchmark/HG002.{T2T,v5q}.truvari_{summary,sizebin}.json`. `docs/demo/HG002_report.html` regenerated — a strict superset of the 2026-07-12 report: adds the STR card + both Validation Benchmark cards while keeping the fresh circos with the 5-locus CNV-trait ring, the blood-group/CNV-trait card, and coverage/flagstat QC. Picard insert-size/duplicate metrics remain N/A (not re-run). Tests: 32 passed.
- [ ] Calibrate AMY1 and LPA KIV-2 absolute scaling against a truth sample (proportional-only today)
- [ ] Add heterozygous-deletion tier (ratio ~0.35–0.65) to `bin/gst_null.py` so 1-copy GSTM1/GSTT1 states report as het-del instead of collapsing to "present"

---

## Infrastructure

- [x] **Results directory consolidation — 2026-07-12.** Merged all `results_*/` run dirs into a single per-sample `results/` tree and a single `work/` scratch dir; de-duplicated the shared PON; reclaimed ~772 GB of disposable `work_*` scratch + regenerable cache. See [ORGANIZATION.md](ORGANIZATION.md).
- [x] Nextflow DSL2 pipeline (`main.nf`, `workflows/svcaller.nf`)
- [x] Docker profile (`conf/docker.config`) — all quay.io biocontainer tags verified
- [x] PON build workflow (`workflows/pon_build.nf`) — complete; HDF5 at `pon/pon/giab_cnv_pon.hdf5`
- [x] GIAB samplesheet (`validation/giab_samplesheet.csv`)
- [x] SV validation samplesheet — HG002 only (`validation/validation_samplesheet.csv`); Truvari runs against GIAB truth
- [x] SMN validation samplesheet — SMA trio only (`validation/smn_validation_samplesheet.csv`); run with `--skip_gridss true`, no `--giab_truth`
- [x] WGS intervals BED (`/data/alvin/ref/GRCh38/wgs_autosomal.bed`)
- [x] `svcaller/utils:1.2` Docker image — built (openpyxl 3.1.5 + report template)
- [x] `svcaller/melt:2.2.2` Docker image — built locally from MELTv2.2.2.tar.gz
- [x] `svcaller/smncopynum:1.1` Docker image — built
- [x] GIAB PON sample report samplesheet (`validation/giab_reports_samplesheet.csv`) — HG001, HG003-HG007
- [x] TMPDIR fix — containers bind-mount `/data/alvin/tmp` as `/tmp` to avoid root fs overflow

## Known Bugs / Blockers

- [ ] **AnnotSV DB** — requires `--annotsv_db /data/alvin/ref/annotsv/Annotations_Human`; annotation step emits empty stub TSV without it.
- [ ] **samtools flagstat not in MultiQC** — flagstat file mixed into `ch_multiqc_files` but MultiQC may need explicit module config to parse it.

## Recently Fixed

- [x] **SMN rename bug** — SMNCopyNumberCaller outputs `{sample}.tsv` (not `{sample}_smn.tsv`); touch fallback was always triggered → fixed rename order in `modules/smn_caller/call.nf` to check `${meta.id}.tsv` before fallback — 2026-05-23
- [x] **SMAD/SMAM truth table swap** — labels were transposed in `validation/smn_truth_table.tsv` and `Documentation.md`; corrected from clinical records: SMAM=SMN2×5, SMAD=SMN2×1 — 2026-05-23
- [x] **DELLY bcftools missing** — bcftools not in `quay.io/biocontainers/delly:1.2.6` container; rewrote `modules/delly/call.nf` to emit VCF directly merged with shell + bgzip + tabix — 2026-05-23
- [x] **GRIDSS_SETUP flag** — `--setupworkingdir` was invalid; corrected to `--steps setupreference`; removed `.img`/`.gridsscache` from outputs (not produced by setupreference in GRIDSS 2.13.2) — 2026-05-23
- [x] **GRIDSS BND→SV conversion** — GRIDSS outputs BND pairs; `bin/gridss_convert_bnd.py` converts them to typed DEL/DUP/INV before JASMINE merge; GRIDSS now contributes 3,095 SVs (14% of 21,735 merged) — 2026-06-08
- [x] **MELT INFO header stripping** — MELT INFO fields (DIFF/LP/RP/RA/PRIOR/SR/MEINFO) dropped by JASMINE causing bcftools/Truvari fatal exit; fixed by stripping to SVTYPE/MEITYPE/SVLEN/END in `modules/melt/call.nf` + injecting MEITYPE header in `modules/jasmine/merge.nf` — 2026-06-08
- [x] **SV_PON annotation** — `--sv_pon` param was never wired into run commands; SV_PON_ANNOTATE silently skipped → 0 PON hits on all previous runs; fixed by adding `--sv_pon pon/sv_pon/giab_sv_pon.bed` to all run commands; 609/7818 HG002 SVs now correctly flagged (7.8%) — 2026-06-08
- [x] **3-tier clinical HTML report** — replaced flat "top SVs" table with Tier 1 (ACMG SF v3.2), Tier 2 (OMIM morbid, top 10), Tier 3 (all, top 10 + XLS download); HTML reduced from 7 MB → 2 MB — 2026-06-08
- [x] **XLS export** — 4-sheet openpyxl workbook (SV/CNV/STR/SMN) with full untruncated tables; download button in HTML report header — 2026-06-08
- [x] **svcaller/utils:1.2** — added openpyxl==3.1.5 and `import sys` fix to `bin/html_report.py`; rebuilt container — 2026-06-08

## Performance

### Implemented
- [x] Conditional sort skip in `SAMTOOLS_FILTER_CHROMS` — skips ~3h re-sort when BAM chr order already matches FAI order
- [x] `GRIDSS_SETUP` process with `storeDir` — pre-builds BWA index via `--steps setupreference`, staged into all GRIDSS_CALL tasks to skip per-sample rebuild (~40 min × N samples saved)
- [x] `--skip_gridss` flag — Manta+DELLY only (min_support=1); saves ~4–5h critical-path time when throughput > sensitivity
- [x] **DELLY parallel SV types** — `DELLY_CALL_SVTYPE` fans out DEL/INS/INV/DUP/BND as 5 concurrent Nextflow processes; `DELLY_MERGE` collects via `groupTuple(size:5)` — ~5× DELLY speedup (2026-05-23)
- [x] **Tiered GRIDSS** (`--tiered_gridss`) — Manta runs first; `MANTA_RESIDUAL_REGIONS` extracts non-PASS loci ±1 kb; `SAMTOOLS_SUBSET --fetch-pairs` creates a region-subset BAM; GRIDSS runs only on that subset — estimated wall time reduction from ~4–6 h to ~30–60 min per sample (2026-05-23)

### Future (lower priority)
- [ ] Increase GRIDSS thread count above 16 (test 24–32; diminishing returns expected)
- [ ] Publish filtered BAMs to `outdir` — allow re-supply as BAM input to skip FILTER_CHROMS on re-runs
- [ ] Pre-built canonical-only reference FASTA — eliminate FILTER_CHROMS entirely for BAMs aligned to it
- [ ] Benchmark Lumpy/Smoove as faster GRIDSS replacement
- [ ] Chromosome-scatter GRIDSS — split BAM by chr1-22+X+Y, run 24 parallel GRIDSS processes, merge VCFs; saves ~24× on large servers but requires BND/TRA post-merge deduplication

## Validation

### Implemented
- [x] GIAB SV v0.6 Truvari benchmark (HG002, ~12K SVs; deletion-biased)
- [x] `validation/download_refs.sh` extended to attempt GIAB SV v1.0 download with FTP probe + graceful fallback

### Pending
- [ ] **GIAB SV v1.0** — multi-platform HiFi+ONT+short-read benchmark (~75K SVs, insertion-inclusive); verify URL at NCBI FTP `NIST_SV_v1.0/` and update `download_refs.sh`; pass `--giab_truth ${GIAB_DIR}/HG002_SV_v1.0.vcf.gz`; expect lower recall (harder benchmark) — that is expected and a more honest measure
- [ ] **CMRG benchmark** — 273 clinically medically relevant genes (CYP450, BRCA, PMS2, etc.); run Truvari with CMRG BED as `--includebed`; wire second benchmark JSON to HTML report as "Clinical genes" section; directly relevant for rare disease diagnostics use case

## Robustness / Output Integrity

### Implemented
- [x] **No empty placeholder outputs on failure** (2026-07-13) — removed every `|| touch <output>` / `echo '{}' > <output>` fallback that let a failed caller exit 0 and publish a zero-byte artifact. Motivated by a production incident: an empty `smn.tsv` from a crashed `SMN_CALLER` was read by OmniGen (which gated on `os.path.exists()`, true for a 0-byte file) and rendered as a complete consumer report reading "0 Carrier findings, 0 Medication flags, Clear" — a crashed caller shown as a clean bill of health. Sites fixed: `modules/smn_caller/call.nf`, `modules/annotsv/annotate.nf`, `modules/expansionhunter/call.nf`, `modules/svaba/call.nf`, `modules/scramble/call.nf`, `modules/melt/call.nf`. Legitimate empty results (header + zero rows) and explicit `--skip_*`/`*_STUB` paths are preserved. Guarded by `tests/test_no_empty_placeholders.py`. See `docs/CHANGES.md`.

### Pending
- [x] **Audit `bin/html_report.py` / `bin/smn_report.py` for `os.path.exists()`-style silent-failure gating** (2026-07-15) — audited both generators and their helpers. Two genuine findings: `smn_report.parse_smn_tsv` defaulted a 0-byte/header-only SMN TSV to "Normal (CN=2)", and `html_report.render_report` read its required inputs (`--smn-html`, `--circos-svg`, `--sv-tsv`, `--cnv-bed`) with no non-empty check, rendering a 0-byte crash placeholder as "no findings". Added fail-loud guards (`SmnInputError`; `UpstreamEmptyError` + `_is_absent`/`_require_nonempty`) that raise on a present-but-empty required input while preserving the three-state distinction (absent `NO_*` sentinel = skip, header-only = legit negative, empty = failure). Regression tests added; suite 67 → 77 passed. Healthy HG002 path verified unchanged. See `docs/CHANGES.md`. (The `validation/*.sh` `|| true` sweep remains open.)
- [ ] **Assert non-empty published artifacts in an end-to-end smoke run** — the current guard is static (grep-based). A cheap `-stub-run` or small-BAM run asserting every file in `results/<sample>/` is non-zero-byte would close the loop. Not done here: no pipeline run was performed for this change.

## Next Steps (Priority Order)

1. Improve recall (~25%) — explore low-QUAL Manta/Delly rescue or soft GRIDSS QUAL floor
2. Complete CMRG benchmark — 273 clinically relevant genes with Truvari `--includebed`
3. Generate GIAB PON sample reports (HG001, HG003-HG007) once giab_reports_run1 completes
4. Run `nextflow clean -but last -f` to recover disk space from work dir after active runs complete
