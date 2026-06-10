# How to Run Pre-Aligned BAM Inputs

Guidance specific to running the pipeline on pre-aligned BAMs (clinical samples, SMA validation cohort, externally aligned WGS).

## When to use this guide

Use this guide when your sample is a BAM file rather than paired FASTQs. Common cases:
- Clinical samples delivered as aligned BAMs from a sequencing provider
- SMA validation cohort (`SMAM`, `SMAD`, `SMAPB`)
- Samples from external studies where re-alignment is not needed

## Key differences from FASTQ inputs

| Aspect | FASTQ input | BAM input |
|--------|------------|-----------|
| FILTER_CHROMS | Skipped (use `hg38.canonical.fa`) | Always runs (~25 min) |
| Reference FASTA | `hg38.canonical.fa` recommended | `hg38.canonical.fa` required (see note below) |
| Alignment step | BWA-MEM2 + SAMTOOLS_SORT | Skipped |
| @SQ header risk | None (canonical ref produces canonical headers) | Present — must be corrected by FILTER_CHROMS |

> **Note:** `hg38.canonical.fa` is required for BAM inputs, not optional. FILTER_CHROMS strips alt contigs from the BAM's reads and `@SQ` headers, leaving only chr1-22+X+Y+M. If the reference still contains alt contigs (as `hg38.fa` does), Manta fails with "BAM/CRAM file is missing a chromosome found in the reference fasta file." Using the canonical reference ensures both the BAM and the reference agree on exactly 25 chromosomes.

## Why BAM inputs always run FILTER_CHROMS

Pre-aligned BAMs may have been aligned to the full GRCh38 reference containing 3,366 chromosomes (canonical + ALT + decoy + HLA + unplaced scaffolds). Even if reads on non-canonical chromosomes are absent, the BAM header retains all `@SQ` lines. Manta's assembly phase silently crashes when the header references chromosomes without reads, producing 0 SV calls despite the scanner finding hundreds of thousands of candidates.

FILTER_CHROMS removes both the non-canonical reads and the non-canonical `@SQ` header entries, leaving exactly 25 canonical chromosomes (chr1-22+X+Y+M) in both header and reads.

## Running multiple BAM samples in parallel

Each sample must have its own work directory to prevent session lock conflicts.

```bash
# Create samplesheets for each sample first:
# validation/smn_SMAM_samplesheet.csv
# validation/smn_SMAD_samplesheet.csv
# validation/smn_SMAPB_samplesheet.csv

for SAMPLE in SMAM SMAD SMAPB; do
  NXF_ANSI_LOG=false nohup nextflow run /data/alvin/SVcaller/main.nf \
    -profile docker \
    --input       /data/alvin/SVcaller/validation/smn_${SAMPLE}_samplesheet.csv \
    --ref_fasta   /data/alvin/ref/GRCh38/hg38.canonical.fa \
    --intervals   /data/alvin/ref/GRCh38/wgs_autosomal.bed \
    --pon         /data/alvin/SVcaller/pon/pon/giab_cnv_pon.hdf5 \
    --eh_catalog  /data/alvin/SVcaller/assets/eh_catalog.json \
    --sv_pon      /data/alvin/SVcaller/pon/sv_pon/giab_sv_pon.bed \
    --skip_gridss true \
    --outdir      /data/alvin/SVcaller/results_smn \
    -work-dir     /data/alvin/SVcaller/work_smn/${SAMPLE} \
    > /data/alvin/tmp/smn_${SAMPLE}_run1.log 2>&1 &
  echo "${SAMPLE} PID: $!"
done
```

Or use the pre-built script: `bash validation/run_smn_parallel.sh`

## Verifying the filtered BAM

After FILTER_CHROMS completes, verify the output BAM has correct canonical headers:

```bash
# Should return exactly 25
samtools view -H work_SAMPLE/*/SAMTOOLS_FILTER_CHROMS/*/SAMPLE.filtered.bam \
  | grep "^@SQ" | wc -l

# Should return chr1..chr22, chrX, chrY, chrM (in that order)
samtools view -H work_SAMPLE/*/SAMTOOLS_FILTER_CHROMS/*/SAMPLE.filtered.bam \
  | grep "^@SQ" | awk '{print $2}' | sed 's/SN://'
```

A count of 3366 (or any value >25) indicates FILTER_CHROMS did not apply the canonical @SQ filter correctly, and Manta will produce 0 variants.

## SMN-specific runs (SMA validation cohort)

For the SMA trio (SMAM = affected child, SMAD = father carrier, SMAPB = mother carrier), GRIDSS can be skipped since the goal is SMN copy number, not SV benchmarking:

```bash
--skip_gridss true   # saves 4-6 h; Manta+Delly+Scramble+MELT still run
```

The AnnotSV database is not required for SMA samples since they are not clinical unknowns:

```bash
# Omit --annotsv_db — report will have blank annotation sections, which is acceptable
```

Expected SMN1/SMN2 copy numbers for the SMA trio:

| Sample | Relationship | Expected SMN1 | Expected SMN2 |
|--------|-------------|---------------|---------------|
| SMAM   | Affected child | 0 | 4 |
| SMAD   | Father (carrier) | 1 | 3 |
| SMAPB  | Mother (carrier) | 1 | 3 |

## BAM requirements

The input BAM must be:
- Sorted by coordinate (not name-sorted)
- Duplicate-marked (Picard MarkDuplicates or equivalent)
- Indexed (`.bai` at the same path — not `.bam.bai` in a different directory)
- Aligned to GRCh38 (hg38) — hg19/GRCh37 is not supported

Check index existence before running:
```bash
ls -lh /path/to/SAMPLE.bam /path/to/SAMPLE.bam.bai
```

## Related

- [How to run a clinical sample](howto-run-clinical-sample.md) — general workflow with all options
- [Design decisions — FILTER_CHROMS canonical @SQ headers](explanation-design.md#canonical-sq-headers-in-filter_chroms-output-bam)
- [Parameter reference](reference-parameters.md)
