# How to Run the GIAB Validation Benchmark

Run the full SV/CNV pipeline on HG002 and benchmark against the GIAB truth set to measure precision, recall, and F1 score.

## Prerequisites

- Docker running on the host
- Nextflow installed (`nextflow -version` should show ≥23.x)
- GIAB truth files downloaded to `/path/to/ref/GIAB/`:
  - `HG002.bwa.sortdup.bqsr.bam` and `.bai`
  - `GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz` and `.tbi`
- PON built at `/path/to/SVcaller/pon/pon/giab_cnv_pon.hdf5`
- AnnotSV database at `/path/to/ref/annotsv/Annotations_Human/`
- Canonical reference FASTA at `/path/to/ref/GRCh38/hg38.canonical.fa`

## Steps

1. Confirm the validation samplesheet points to the HG002 BAM:

   ```
   sample,fastq_1,fastq_2,bam
   HG002,,,/path/to/ref/GIAB/HG002.bwa.sortdup.bqsr.bam
   ```

   File is at `validation/validation_samplesheet.csv`.

2. Run the pipeline with `-profile docker`:

   ```bash
   nextflow run main.nf -profile docker \
     --input validation/validation_samplesheet.csv \
     --ref_fasta /path/to/ref/GRCh38/hg38.canonical.fa \
     --intervals /path/to/ref/GRCh38/wgs_autosomal.bed \
     --pon /path/to/SVcaller/pon/pon/giab_cnv_pon.hdf5 \
     --giab_truth /path/to/ref/GIAB/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz \
     --eh_catalog assets/eh_catalog.json \
     --annotsv_db /path/to/ref/annotsv/Annotations_Human \
     --outdir results \
     -work-dir work \
     -resume \
     > /path/to/tmp/run.log 2>&1
   ```

   Monitor progress in another terminal: `tail -f /path/to/tmp/run.log`

3. Wait for completion. Typical runtimes:
   - M1 Preprocessing (incl. FILTER_CHROMS for BAM input): ~25-35 min
   - M2 SV Calling (Manta + Delly + GRIDSS full mode): ~6-7 h
   - M2 with `--skip_gridss true`: ~45 min
   - Full pipeline: ~8-9 h (with GRIDSS), ~2 h (without)

## Verification

Check the benchmark results appear in the HTML report:

```bash
ls -lh results/HG002/HG002.report.html
# Should be ~1-3 MB

# Extract F1 from the Truvari JSON directly:
cat work/*/*/summary.json | python3 -c "import json,sys; d=json.load(sys.stdin); print(f\"F1={d['f1']:.3f} P={d['precision']:.3f} R={d['recall']:.3f}\")"
```

Open `results/HG002/HG002.report.html` in a browser. The **GIAB Benchmark** section shows overall F1 and per-size-bin metrics (50-300 bp, 300 bp-1 kb, 1-10 kb, >10 kb).

## Troubleshooting

**Pipeline halts at MOSDEPTH with "mean depth below threshold"**
The HG002 BAM must be ≥30× mean coverage. Run `mosdepth --no-abbrev HG002 HG002.bam` and check `HG002.mosdepth.summary.txt`. If coverage is lower, pass `--min_depth 20` to lower the threshold.

**GRIDSS runs out of memory (exit code 137)**
GRIDSS requires 60 GB RAM. Use `--skip_gridss true` to skip it or `--tiered_gridss true` to run it on Manta residual regions only (~40 GB).

**BUILD_HTML_REPORT receives zero inputs (report section blank)**
This usually means a channel meta-map mismatch. All channels joining into BUILD_HTML_REPORT must carry identical meta maps. Check that the BAM stub in `preprocess.nf` uses `[*:meta, needs_chr_filter: true]`.

**Truvari section missing from report**
Verify `--giab_truth` was passed and the `.tbi` index file exists at the same path as the VCF.

## Related

- [How to build a Panel of Normals](howto-build-pon.md)
- [How to interpret the HTML report](howto-interpret-report.md)
- [Parameter reference](reference-parameters.md)
