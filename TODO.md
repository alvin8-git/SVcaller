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
- [x] GRIDSS (`modules/gridss/call.nf`)
- [x] ExpansionHunter — STR loci (`modules/expansionhunter/call.nf`)
- [x] JASMINE merge (min_support=2) (`modules/jasmine/merge.nf`)
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
- [x] HTML report section 8: Top annotated SVs ACMG 4/5
- [x] Truvari per-size-bin metrics — 4 bins wired to HTML report

---

## Infrastructure

- [x] Nextflow DSL2 pipeline (`main.nf`, `workflows/svcaller.nf`)
- [x] Docker profile (`conf/docker.config`) — all quay.io biocontainer tags verified
- [x] PON build workflow (`workflows/pon_build.nf`) — complete; HDF5 at `pon/pon/giab_cnv_pon.hdf5`
- [x] GIAB samplesheet (`validation/giab_samplesheet.csv`)
- [x] Validation samplesheet HG002 (`validation/validation_samplesheet.csv`)
- [x] WGS intervals BED (`/data/alvin/ref/GRCh38/wgs_autosomal.bed`)
- [x] `svcaller/utils:1.0` Docker image — built
- [x] `svcaller/smncopynum:1.1` Docker image — built
- [x] TMPDIR fix — containers bind-mount `/data/alvin/tmp` as `/tmp` to avoid root fs overflow

## Known Bugs / Blockers

- [ ] **AnnotSV DB** — not downloaded; annotation step emits empty stub TSV without `--annotsv_db`.
- [ ] **samtools flagstat not in MultiQC** — flagstat file mixed into `ch_multiqc_files` but MultiQC may need explicit module config to parse it.
- [ ] **SMA BAMs** — SMAD/SMAM/SMAPB alignments done externally; awaiting BAMs to add back to samplesheet.

## Next Steps (Priority Order)

1. Add SMA samples back to `validation/validation_samplesheet.csv` once external BAMs arrive
2. Download AnnotSV database: set `--annotsv_db` parameter for full annotation
3. Verify flagstat appears correctly in HTML QC section after first HG002 run completes
