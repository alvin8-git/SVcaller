# SVcaller TODO

Tracks implementation status against the design spec (`docs/superpowers/specs/2026-05-09-svcnv-caller-design.md`).

---

## M1 — Preprocessing

- [x] BWA-MEM2 alignment (`modules/bwamem2/align.nf`)
- [x] Picard MarkDuplicates (`modules/picard/markduplicates.nf`)
- [x] mosdepth coverage QC with min-depth halt (`modules/mosdepth/coverage.nf`)
- [x] FastQC raw read QC (`modules/fastqc/qc.nf`)
- [ ] samtools flagstat output (spec mentions it; not yet wired)

## M2 — SV Calling

- [x] Manta (`modules/manta/call.nf`)
- [x] DELLY (`modules/delly/call.nf`)
- [x] GRIDSS (`modules/gridss/call.nf`)
- [x] ExpansionHunter — STR loci (`modules/expansionhunter/call.nf`)
- [x] JASMINE merge (min_support=2) (`modules/jasmine/merge.nf`)
- [x] sv_tbi emitted for Truvari benchmarking

## M3 — CNV Calling

- [x] CNVpytor (`modules/cnvpytor/call.nf`)
- [x] GATK gCNV case mode (`modules/gatk/gcnv_call.nf`)
- [x] Consensus CNV BED (`bin/cnv_consensus.py`)
- [x] GATK PreprocessIntervals in case mode (`subworkflows/cnv_calling.nf`) — added 2026-05-12
- [x] GATK PreprocessIntervals in PON build (`workflows/pon_build.nf`) — added 2026-05-12

## M4 — SMN1/SMN2

- [x] SMNCopyNumberCaller v1.1 (`modules/smn_caller/call.nf`)
- [x] SMN HTML section (`bin/smn_report.py`)
- [x] 2+0 haplotype detection — verified correct in `bin/smn_report.py` 2026-05-12
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
- [x] Ring 4: STR expansion markers (brown `#8C564B`) — added 2026-05-12
- [x] Ring 5: SMN locus highlight chr5q13.2 (gold `#FFBF00`)
- [x] Center: SV links (type-coloured)
- [x] PNG at 1200 dpi (per spec)
- [x] circos.png published to outdir — publishDir added to `modules/pycirclize/plot.nf` 2026-05-12

## M7 — Reporting & Benchmarking

- [x] Jinja2 HTML per-sample report (`bin/html_report.py`)
- [x] MultiQC aggregated QC report (`modules/multiqc/report.nf`) — added 2026-05-12
- [x] Truvari GIAB benchmarking (HG002 SV v0.6) — wired via `--giab_truth`
- [x] HTML report section 2: Alignment QC — mosdepth + Picard metrics parsed in `bin/html_report.py` 2026-05-12
- [x] HTML report section 7: STR expansion loci — ExpansionHunter VCF parsing + template section added 2026-05-12
- [x] HTML report section 8: Top annotated SVs ACMG 4/5 — `parse_top_svs()` filters AnnotSV_ranking_score ≥ 0.9
- [x] Truvari per-size-bin metrics — 4 bins added to TRUVARI_BENCH, sizebin JSON wired to HTML report 2026-05-12

---

## Infrastructure

- [x] Nextflow DSL2 pipeline (`main.nf`, `workflows/svcaller.nf`)
- [x] Docker profile (`conf/docker.config`) — image tags need verification against registry
- [x] Local conda profile (`conf/local.config`, `environment.yml`) — conda env building
- [x] PON build workflow (`workflows/pon_build.nf`) — entry workflow added, PreprocessIntervals added
- [x] GIAB samplesheet (`validation/giab_samplesheet.csv`) — paths corrected
- [x] Validation samplesheet HG002 + SMA trio (`validation/validation_samplesheet.csv`)
- [x] WGS intervals BED (`/data/alvin/ref/GRCh38/wgs_autosomal.bed`) — generated
- [x] `svcaller/utils:1.0` Docker image — built
- [x] `svcaller/smncopynum:1.1` Docker image — built

## Known Bugs / Blockers

- [ ] **Conda env** — `svcaller` env still failing due to pywfa conflict; build5 running with cnvpytor moved to pip (`/data/alvin/tmp/conda_svcaller_build5.log`).
- [ ] **PON build re-run** — PON build run2 in progress (`/data/alvin/tmp/pon_build_run2.log`). All 7 GIAB BAM paths now correct, PreprocessIntervals added, using `-profile docker`.
- [ ] **Docker image tags** — quay.io biocontainer tags in `conf/docker.config` return "manifest unknown". Use `local` profile (conda) for main pipeline; docker works for PON build (GATK only).
- [ ] **AnnotSV DB** — not downloaded; annotation step emits empty stub TSV without `--annotsv_db`. AnnotSV not in conda env (use Docker or skip).
- [ ] **samtools flagstat** — `mapped_pct` shows "N/A" in HTML QC section.

## Next Steps (Priority Order)

1. Wait for conda env build to complete (`/data/alvin/tmp/conda_svcaller_build3.log`), verify with `conda activate svcaller`
2. Re-run PON build (giab_samplesheet paths now correct, PreprocessIntervals added):
   `nextflow run workflows/pon_build.nf --input validation/giab_samplesheet.csv --ref_fasta /data/alvin/ref/GRCh38/GRCh38.fasta --intervals /data/alvin/ref/GRCh38/wgs_autosomal.bed --outdir /data/alvin/SVcaller/pon -profile local`
3. Re-run main validation pipeline with PON: `nextflow run main.nf -profile local --input validation/validation_samplesheet.csv --pon /data/alvin/SVcaller/pon/gcnv_pon.hdf5 ...`
4. Add per-size-bin Truvari metrics (50–300 bp, 300 bp–1 kb, 1–10 kb, >10 kb)
5. Download AnnotSV database and set `--annotsv_db` parameter
6. Add samtools flagstat output (spec item, not yet wired in preprocess subworkflow)
