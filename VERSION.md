# SVcaller Version History

## v1.0.0 (in development — 2026-05-12)

**Initial implementation.** All core modules implemented per design spec.

### Implemented
- M1: BWA-MEM2 → Picard MarkDup → mosdepth → FastQC
- M2: Manta + DELLY + GRIDSS + ExpansionHunter → JASMINE merge
- M3: CNVpytor + GATK gCNV → consensus CNV BED
- M4: SMNCopyNumberCaller v1.1
- M5: AnnotSV + gnomAD-SV AF filter
- M6: pycirclize Circos (5 rings: ideogram, CNV gain/loss, STR, SMN locus; SV links)
- M7: Jinja2 HTML report + MultiQC + Truvari GIAB benchmarking
- PON build workflow (GATK gCNV Panel of Normals from GIAB HG001–HG007)
- Validation samplesheet: HG002 (GIAB) + SMA trio (SMAPB/SMAD/SMAM)

### Containers
| Tool | Image |
|------|-------|
| All pipeline tools | `svcaller` conda env (`environment.yml`) |
| Custom utils | `svcaller/utils:1.0` |
| SMNCopyNumberCaller | `svcaller/smncopynum:1.1` |

### Reference Data
| File | Path |
|------|------|
| GRCh38 FASTA | `/data/alvin/ref/GRCh38/GRCh38.fasta` |
| WGS intervals BED | `/data/alvin/ref/GRCh38/wgs_autosomal.bed` |
| GIAB BAMs HG001–HG007 | `/data/alvin/ref/GIAB/` |
| GIAB SV truth (HG002) | `/data/alvin/ref/GIAB/HG002_SV_v0.6.vcf.gz` |

### Additional Fixes (2026-05-12 session 2)
- M3: PreprocessIntervals added to case CNV workflow (`subworkflows/cnv_calling.nf`)
- M7: HTML report section 2 (Alignment QC) — mosdepth + Picard metrics parsed
- M7: HTML report section 7 (STR expansion loci) — ExpansionHunter VCF parsing + template
- M7: HTML report section 8 (Top annotated SVs) — verified ACMG score threshold
- M7: Truvari per-size-bin metrics (4 bins) wired into HTML report
- M7: Truvari `parse_benchmark()` bug fixed (was not parsing flat JSON format)
- M6: `circos.png` publishDir added to `modules/pycirclize/plot.nf`
- M5: AnnotSV graceful skip when `--annotsv_db` not provided
- M4: SMN 2+0 haplotype detection verified correct
- `environment.yml`: added fastqc, multiqc; fixed python=3.9 for pywfa compatibility
- giab_samplesheet.csv BAM paths corrected; PON re-run pending

### Known Issues
- PON build not yet re-run (giab_samplesheet paths fixed, awaiting conda env completion)
- Docker biocontainer tags invalid; use `-profile local` (conda)
- Conda env build requires python=3.9 (pywfa 0.5.1 not available for python≥3.10)
- AnnotSV excluded from conda env (run via Docker or provide `--annotsv_db` separately)
- samtools flagstat not yet wired (mapped_pct shows "N/A" in QC section)

---

## Planned: v1.1.0

- Nanopore long-read support (Sniffles2 + CuteSV, minimap2 alignment)
- Full HTML report sections 2, 7, 8
- Per-size-bin Truvari metrics
- AnnotSV database auto-download script
