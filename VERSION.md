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

---

## v1.1.0 (in development — 2026-06-08)

**Major feature release.** 6-caller ensemble, 3-tier clinical report, SV PON, XLS export.

### New Callers
- Scramble MEI caller (L1/ALU/SVA mobile element insertions)
- MELT MEI caller (ALU/HERVK/LINE1/SVA; local build from MELTv2.2.2.tar.gz)
- SvABA local-assembly caller (DEL/INV/BND from local read assembly)
- STRling genome-wide STR scanner (complements ExpansionHunter's 32-locus catalog)
- GRIDSS BND→SV converter (`bin/gridss_convert_bnd.py`): BND pairs auto-converted to typed DEL/DUP/INV before JASMINE merge

### Annotation
- SV Panel of Normals (`pon/sv_pon/giab_sv_pon.bed`, 10,576 sites from 7 GIAB samples)
- SegDup boundary badge
- ENCODE blacklist badge
- gnomAD-SV AF badge (soft annotation, not hard filter)

### Clinical Report
- 3-tier SV classification: Tier 1 (ACMG SF v3.2), Tier 2 (OMIM morbid, top 10), Tier 3 (all, top 10 + XLS)
- XLS export — 4-sheet openpyxl workbook (SV / CNV / STR / SMN), full untruncated tables
- HTML report size reduced from ~7 MB to ~2 MB

### STR Report
- ExpansionHunter catalog expanded 1 → 32 disease loci
- INREPEAT / INTERMEDIATE / NORMAL status badge per locus
- Novel STRling candidates suppressed from clinical table

### Performance
- Canonical reference (`hg38.canonical.fa`) skips FILTER_CHROMS for FASTQ-derived BAMs (~25 min saved)

### Containers
| Tool | Image |
|------|-------|
| Python utils / report | `svcaller/utils:1.3` |
| MELT | `svcaller/melt:2.2.2` |
| SMNCopyNumberCaller | `svcaller/smncopynum:1.1` |

### GIAB Benchmark (HG002, run16)
| Benchmark | Precision | Recall | F1 |
|-----------|-----------|--------|-----|
| T2TQ100-V1.0 | 0.733 | 0.255 | 0.378 |
| v5.0q | 0.738 | 0.259 | 0.383 |

### Planned: v1.2.0

- Nanopore long-read support (Sniffles2 + CuteSV, minimap2 alignment)
- CMRG benchmark — 273 clinically relevant genes
- Improved recall: low-QUAL Manta/Delly rescue, soft GRIDSS QUAL floor
