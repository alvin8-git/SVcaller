# Tutorial: Your First SVcaller Run

You'll run the SVcaller pipeline on a real WGS sample (HG002, a GIAB reference standard) and open the resulting HTML report. By the end you'll understand the full pipeline flow and know what each output file contains.

## What you'll need

- Linux machine with Docker installed and running
- Nextflow: `curl -s https://get.nextflow.io | bash` or `conda install nextflow`
- 60 GB free disk space (GRIDSS intermediate files are large)
- 16 GB RAM minimum; 64 GB recommended if running GRIDSS
- HG002 BAM and the GIAB truth set (see prerequisites below)

**Quick check — confirm Docker and Nextflow work:**

```bash
docker run hello-world
nextflow -version   # expect ≥ 23.x
```

## Step 1: Clone the repository and verify assets

```bash
git clone https://github.com/alvin8-git/SVcaller.git
cd SVcaller
ls assets/
# Should show: eh_catalog.json  GRCh38_cytobands.txt  report_template.html  schema_input.json
```

## Step 2: Set up reference files

You need three things on disk before running:

| File | Path | Notes |
|------|------|-------|
| GRCh38 canonical FASTA | `/data/alvin/ref/GRCh38/hg38.canonical.fa` | chr1-22+X+Y+M only; BWA-MEM2 index at same prefix |
| Autosomal intervals BED | `/data/alvin/ref/GRCh38/wgs_autosomal.bed` | Used by GATK CNV |
| GATK gCNV PON | `/data/alvin/SVcaller/pon/pon/giab_cnv_pon.hdf5` | Built from 7 GIAB normals |

**If these already exist, skip to Step 3.** If not, see [How to build a Panel of Normals](howto-build-pon.md) for the PON, and download the GRCh38 reference from NCBI/Ensembl.

## Step 3: Create a samplesheet

Create a CSV pointing to your BAM (or FASTQ pair):

```bash
cat > my_samplesheet.csv << 'EOF'
sample,fastq_1,fastq_2,bam
HG002,,,/data/alvin/ref/GIAB/HG002.bwa.sortdup.bqsr.bam
EOF
```

You should see something in your browser within seconds after the pipeline finishes.

## Step 4: Run the pipeline (skip GRIDSS to save time)

For a first run, skip GRIDSS (saves 5-6 h and 60 GB RAM). You still get Manta, Delly, and Scramble calls.

```bash
nextflow run main.nf -profile docker \
  --input my_samplesheet.csv \
  --ref_fasta /data/alvin/ref/GRCh38/hg38.canonical.fa \
  --intervals /data/alvin/ref/GRCh38/wgs_autosomal.bed \
  --pon /data/alvin/SVcaller/pon/pon/giab_cnv_pon.hdf5 \
  --eh_catalog assets/eh_catalog.json \
  --annotsv_db /data/alvin/ref/annotsv/Annotations_Human \
  --outdir results \
  -work-dir work \
  --skip_gridss true \
  -resume \
  > /data/alvin/tmp/run.log 2>&1 &

# Watch progress:
tail -f /data/alvin/tmp/run.log
```

You should see Nextflow start printing process status lines:

```
executor >  local (1)
[xx/xxxxxx] SVCALLER:PREPROCESS:BWAMEM2_ALIGN (HG002) | 0 of 1
...
```

## Step 5: Check progress

After ~25 min (FILTER_CHROMS for BAM input), then another ~30-45 min for the SV callers, the pipeline completes. Look for:

```
[xx/xxxxxx] SVC…:BUILD_HTML_REPORT (HG002) | 1 of 1 ✔
Completed at: ...
Duration   : ...
CPU hours  : ...
```

Any line showing `FAILED` means a process crashed — check `work/<hash>/.command.err` for the error.

## Step 6: Open the report

```bash
ls -lh results/HG002/
# HG002.report.html   ~2 MB
# HG002.sv_merged.vcf.gz
# HG002.cnv_consensus.bed
# HG002.smn.tsv
# HG002.filtered.tsv
# HG002.circos.png

xdg-open results/HG002/HG002.report.html   # Linux
# or: open results/HG002/HG002.report.html  # macOS
```

The report opens with:
- **Alignment QC** — coverage, duplicate rate
- **SV Summary** — counts by type with ACMG class 4/5 highlights
- **Circos plot** — genome-wide overview (see the coverage ring, SV links, ACMG dots)
- **SMN copy number** — SMN1/SMN2 CN status
- **Top Annotated SVs** — ranked by AnnotSV score

## What you built

You ran a complete germline WGS SV/CNV pipeline that:
- Called SVs with 3 callers (Manta, Delly, Scramble) and merged them into a single VCF
- Called CNVs with CNVpytor and GATK gCNV
- Determined SMN1/SMN2 copy number
- Annotated all SVs with gene overlaps, population frequency, and ACMG classification
- Assembled everything into a clinical HTML report with a genome-wide Circos plot

**Next steps:**
- Add `--giab_truth` to enable Truvari benchmarking against GIAB truth
- Add GRIDSS for precise BND breakpoints: remove `--skip_gridss true`
- See [How to interpret the HTML report](howto-interpret-report.md) for clinical interpretation guidance
- See [Parameter reference](reference-parameters.md) for all available options
