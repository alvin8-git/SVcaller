# How to Build a Panel of Normals

Build a GATK gCNV Panel of Normals (PON) from a set of normal WGS samples. The PON is used to denoise read-depth signals and improve CNV specificity.

## Prerequisites

- Docker running on the host
- Nextflow installed
- A set of normal WGS BAMs — minimum 5, ideally 10-30 samples with similar library preparation
- Reference FASTA and autosomal intervals BED
- At least 8 GB RAM per sample (GATK CollectReadCounts is memory-intensive)

**The built PON for GIAB samples is already at:**
`/data/alvin/SVcaller/pon/pon/giab_cnv_pon.hdf5` (446 MB, built from HG001-HG007)

Only re-build if you are validating against a different cohort or adding more normals.

## Steps

1. Create a samplesheet with your normal BAMs. Each row is one normal sample:

   ```
   sample,fastq_1,fastq_2,bam
   Normal01,,,/path/to/Normal01.bam
   Normal02,,,/path/to/Normal02.bam
   Normal03,,,/path/to/Normal03.bam
   ```

   Save as `validation/my_pon_samplesheet.csv`.

2. Run the PON build workflow:

   ```bash
   nextflow run workflows/pon_build.nf -profile docker \
     --input validation/my_pon_samplesheet.csv \
     --ref_fasta /data/alvin/ref/GRCh38/hg38.fa \
     --intervals /data/alvin/ref/GRCh38/wgs_autosomal.bed \
     --outdir /data/alvin/SVcaller/pon_new \
     -work-dir /data/alvin/SVcaller/work \
     -resume \
     > /data/alvin/tmp/pon_build.log 2>&1
   ```

3. Wait for completion. Runtime is roughly 30-60 min per sample for CollectReadCounts, plus a final PON creation step. A 7-sample PON takes ~4-5 h.

## Verification

```bash
ls -lh /data/alvin/SVcaller/pon_new/pon/
# Should show giab_cnv_pon.hdf5 (typically 200-600 MB depending on sample count)
```

Use the new PON in the main pipeline:

```bash
nextflow run main.nf -profile docker \
  --input validation/validation_samplesheet.csv \
  --pon /data/alvin/SVcaller/pon_new/pon/giab_cnv_pon.hdf5 \
  ...
```

## Troubleshooting

**GATK PreprocessIntervals fails with "sequence dictionary mismatch"**
This happens when the reference `.dict` file uses alphabetical chromosome order (`chr1, chr10...`) but BAM headers use numeric order (`chr1, chr2...`). The PON build workflow omits `--annotated-intervals` to avoid this. If you see this error, check that you are not passing `--annotated-intervals` manually.

**PON file is very small (<10 MB)**
CollectReadCounts likely produced empty or near-empty outputs. Check that:
- The BAMs are sorted and indexed
- The intervals BED covers the same chromosomes as the BAMs
- Coverage is ≥10× — very low coverage BAMs produce too few counts for reliable denoising

**GATK gCNV calls look noisy after using the new PON**
Add more samples to the PON (aim for 20+). Fewer than 5 samples gives poor noise characterisation. Samples must have similar library prep and sequencing depth.

## Related

- [How to run the GIAB validation](howto-run-validation.md)
- [Architecture reference — M3 CNV calling](reference-architecture.md#m3-cnv-calling-subworkflowscnv_callingnf)
- [Design decisions — No GC correction in the PON](explanation-design.md#no-gc-correction-in-the-gatk-gcnv-pon)
