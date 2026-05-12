include { BWAMEM2_ALIGN    } from '../modules/bwamem2/align'
include { PICARD_MARKDUP   } from '../modules/picard/markduplicates'
include { MOSDEPTH         } from '../modules/mosdepth/coverage'
include { FASTQC           } from '../modules/fastqc/qc'

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

    // Merge aligned + pre-supplied BAMs into MarkDup
    ch_all_bam = BWAMEM2_ALIGN.out.bam
        .join(BWAMEM2_ALIGN.out.bai)
        .mix(ch_bam_in)

    PICARD_MARKDUP(ch_all_bam)

    ch_final_bam = PICARD_MARKDUP.out.bam
        .join(PICARD_MARKDUP.out.bai)

    // Coverage QC — halts pipeline if < min_depth
    MOSDEPTH(ch_final_bam, min_depth)

    emit:
    bam        = ch_final_bam
    coverage   = MOSDEPTH.out.summary
    metrics    = PICARD_MARKDUP.out.metrics
    fastqc_zip = FASTQC.out.zip
}
