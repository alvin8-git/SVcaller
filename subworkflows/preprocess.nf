include { BWAMEM2_ALIGN      } from '../modules/bwamem2/align'
include { SAMTOOLS_SORT      } from '../modules/samtools/sort'
include { SAMTOOLS_FLAGSTAT  } from '../modules/samtools/flagstat'
include { PICARD_MARKDUP     } from '../modules/picard/markduplicates'
include { PICARD_INSERT_SIZE } from '../modules/picard/insert_size'
include { MOSDEPTH           } from '../modules/mosdepth/coverage'
include { FASTQC             } from '../modules/fastqc/qc'

workflow PREPROCESS {
    take:
    ch_samplesheet  // [ meta, fq1|null, fq2|null, bam|null ]
    ch_fasta        // path
    ch_fai          // path
    ch_bwt_index    // path (directory)
    min_depth       // integer

    main:
    // Split into FASTQ and BAM channels
    ch_fastq = ch_samplesheet
        .filter { meta, fq1, fq2, bam -> fq1 != null }
        .map    { meta, fq1, fq2, bam -> [meta, [fq1, fq2]] }

    ch_bam_in = ch_samplesheet
        .filter { meta, fq1, fq2, bam -> bam != null }
        .map    { meta, fq1, fq2, bam -> [meta, bam, file("${bam}.bai")] }

    // Raw read QC (FASTQ samples only)
    FASTQC(ch_fastq)

    // Align FASTQs
    BWAMEM2_ALIGN(ch_fastq, ch_fasta, ch_fai, ch_bwt_index)

    // Sort aligned BAM (samtools is not in the bwa-mem2 container)
    SAMTOOLS_SORT(BWAMEM2_ALIGN.out.bam)

    // MarkDup only on FASTQ-derived BAMs; pre-supplied BAMs are already dup-marked
    PICARD_MARKDUP(SAMTOOLS_SORT.out.bam.join(SAMTOOLS_SORT.out.bai))

    // Tag each BAM so sv_calling.nf can skip FILTER_CHROMS for FASTQ-derived BAMs.
    // FASTQ-derived BAMs are aligned to hg38.canonical.fa → no alt contigs → filter not needed.
    // Pre-supplied BAMs may be aligned to full hg38 (with alt contigs) → filter needed.
    ch_final_bam = PICARD_MARKDUP.out.bam
        .join(PICARD_MARKDUP.out.bai)
        .map { meta, bam, bai -> [[*:meta, needs_chr_filter: false], bam, bai] }
        .mix(
            ch_bam_in.map { meta, bam, bai -> [[*:meta, needs_chr_filter: true], bam, bai] }
        )

    // Coverage QC — halts pipeline if < min_depth
    MOSDEPTH(ch_final_bam, min_depth)

    // Mapping rate QC
    SAMTOOLS_FLAGSTAT(ch_final_bam)

    // Insert size distribution QC
    PICARD_INSERT_SIZE(ch_final_bam)

    // Stub metrics for pre-supplied BAMs so downstream report join has an entry
    ch_markdup_metrics = PICARD_MARKDUP.out.metrics
        .mix(ch_bam_in.map { meta, bam, bai -> [[*:meta, needs_chr_filter: true], file("NO_METRICS")] })

    emit:
    bam           = ch_final_bam
    coverage      = MOSDEPTH.out.summary
    regions_bed   = MOSDEPTH.out.regions_bed
    metrics       = ch_markdup_metrics
    insert_size   = PICARD_INSERT_SIZE.out.metrics
    flagstat      = SAMTOOLS_FLAGSTAT.out.flagstat
    fastqc_zip    = FASTQC.out.zip
}
