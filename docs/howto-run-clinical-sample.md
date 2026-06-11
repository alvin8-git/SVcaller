# How to Run a Clinical Sample

Step-by-step guide for running a new patient sample through the SVcaller pipeline.

## Prerequisites

- Pre-aligned BAM file with a `.bai` index at the same path, **or** paired FASTQ files (R1/R2, gzipped)
- Reference FASTA at `/data/alvin/ref/GRCh38/hg38.canonical.fa` (canonical chromosomes, chr1-22+X+Y+M)
- PON at `/data/alvin/SVcaller/pon/pon/giab_cnv_pon.hdf5`
- AnnotSV database at `/data/alvin/ref/annotsv/Annotations_Human`
- At least 500 GB free disk space for a 30× WGS sample

## Step 1: Create the samplesheet

Create a CSV at a convenient location. Replace `SAMPLEID` with the actual sample identifier (no spaces, no special characters).

**BAM input:**
```
sample,fastq_1,fastq_2,bam
SAMPLEID,,,/path/to/SAMPLEID.bam
```

**FASTQ input:**
```
sample,fastq_1,fastq_2,bam
SAMPLEID,/path/to/SAMPLEID_R1.fq.gz,/path/to/SAMPLEID_R2.fq.gz,
```

Save as `/path/to/SAMPLEID_samplesheet.csv`.

## Step 2: Launch the pipeline

```bash
NXF_ANSI_LOG=false nohup nextflow run /data/alvin/SVcaller/main.nf \
  -profile docker \
  --input       /path/to/SAMPLEID_samplesheet.csv \
  --ref_fasta   /data/alvin/ref/GRCh38/hg38.canonical.fa \
  --intervals   /data/alvin/ref/GRCh38/wgs_autosomal.bed \
  --pon         /data/alvin/SVcaller/pon/pon/giab_cnv_pon.hdf5 \
  --eh_catalog  /data/alvin/SVcaller/assets/eh_catalog.json \
  --annotsv_db  /data/alvin/ref/annotsv/Annotations_Human \
  --sv_pon      /data/alvin/SVcaller/pon/sv_pon/giab_sv_pon.bed \
  --outdir      /data/alvin/SVcaller/results_SAMPLEID \
  -work-dir     /data/alvin/SVcaller/work_SAMPLEID \
  > /data/alvin/tmp/SAMPLEID_run1.log 2>&1 &
echo "PID: $!"
```

Key points:
- Use a dedicated work directory per sample (`work_SAMPLEID`). Never share work dirs across samples.
- Use `hg38.canonical.fa` for FASTQ inputs — the pipeline skips the 25-minute FILTER_CHROMS step.
- BAM inputs always run FILTER_CHROMS regardless of reference.
- The pipeline runs MELT by default (~2 h). Add `--skip_melt true` if the `svcaller/melt:2.2.2` container is unavailable.

## Step 3: Monitor progress

```bash
tail -f /data/alvin/tmp/SAMPLEID_run1.log
```

Expected timeline for a 30× WGS sample with BAM input:

| Step | Time |
|------|------|
| FILTER_CHROMS | ~25 min |
| Manta | ~45 min |
| Delly | ~90 min |
| GRIDSS | ~4-6 h |
| MELT | ~2 h |
| CNV (GATK + CNVpytor) | ~60 min |
| SMN caller | ~20 min |
| AnnotSV | ~30 min |
| Report | ~5 min |
| **Total (all parallel)** | **~8-10 h** |

Use `--skip_gridss true` to reduce runtime to ~3-4 h (reduces SV sensitivity).

## Step 4: Verify outputs

After the pipeline completes, check:

```bash
ls /data/alvin/SVcaller/results_SAMPLEID/SAMPLEID/
```

Expected files:
- `SAMPLEID.report.html` — clinical HTML report
- `SAMPLEID.variants.xlsx` — Excel workbook (SVs, CNVs, STRs, SMN sheets)
- `SAMPLEID.sv_merged.vcf.gz` — merged SV calls
- `SAMPLEID.filtered.tsv` — AnnotSV-annotated, frequency-filtered SVs
- `SAMPLEID.cnv_consensus.bed` — CNV calls
- `SAMPLEID.smn.tsv` — SMN1/SMN2 copy numbers
- `SAMPLEID.circos.svg` + `SAMPLEID.circos.png` — genome Circos plot

Open `SAMPLEID.report.html` in a browser to review the clinical findings.

## Step 5: Clean up

A 30× WGS run leaves large intermediates in the work directory. Once results are confirmed complete, reclaim that space. Three options, in order of preference:

**Option A — `nf-cleanup.sh` (recommended).** Verifies the sample's outputs exist under `--outdir`, then removes the work dir and prunes orphaned `.nextflow/cache` sessions (skipping any still locked):

```bash
bash /data/alvin/SVcaller/bin/nf-cleanup.sh SAMPLEID
```

**Option B — `--auto_cleanup` at launch.** Add `--auto_cleanup true` to the Step 2 command to delete the work dir automatically on **successful** completion. Only use this for one-shot runs — it removes the `-resume` cache, so a failed or re-run sample starts from scratch.

**Option C — manual.** Safe because all results live in `results_SAMPLEID/`, which is untouched:

```bash
rm -rf /data/alvin/SVcaller/work_SAMPLEID
```

Do **not** delete the `storeDir` caches under `results_SAMPLEID/cache/` and `results_SAMPLEID/.cache/` if you plan to run more samples against the same reference — they let later runs skip the GRIDSS reference setup, interval binning, and chrom filtering.

## Troubleshooting

**Pipeline exits at MOSDEPTH with "Coverage below minimum threshold"**
The sample has <30× mean coverage. Lower the threshold with `--min_depth 20` (reduces sensitivity) or investigate the BAM for alignment issues.

**Manta produces 0 variants despite finding candidates**
Check that the BAM was filtered correctly: `samtools view -H SAMPLEID.filtered.bam | grep "^@SQ" | wc -l` should return exactly 25. If it returns 3366+, the BAM was aligned to a non-canonical reference and the `@SQ` headers contain alt-contig entries. Re-run FILTER_CHROMS manually or check the reference used.

**GRIDSS OOM (exit 137)**
`GRIDSS_CALL` auto-retries up to 3 times, climbing 32 → 64 → 96 GB. If it still OOMs (or the machine can't supply 96 GB), use `--skip_gridss true` and rely on Manta+Delly+Scramble+MELT for SV calling. On a small machine, also lower the global cap, e.g. `--max_memory 32.GB`, so no process requests more RAM than exists.

**AnnotSV sections empty in report**
The `--annotsv_db` path must point to the parent directory of `Annotations_Human/`, not to `Annotations_Human/` itself. AnnotSV appends the folder name internally.

## Related

- [Parameter reference](reference-parameters.md) — all CLI flags
- [How to run BAM-input samples](howto-run-bam-inputs.md) — additional guidance for pre-aligned BAMs
- [How to interpret the HTML report](howto-interpret-report.md)
- [Design decisions explained](explanation-design.md)
