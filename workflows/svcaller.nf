include { PREPROCESS   } from '../subworkflows/preprocess'
include { SV_CALLING   } from '../subworkflows/sv_calling'
include { CNV_CALLING  } from '../subworkflows/cnv_calling'
include { SMN_CALLING  } from '../subworkflows/smn_calling'
include { ANNOTATE     } from '../subworkflows/annotate'
include { REPORT       } from '../subworkflows/report'

workflow SVCALLER {
    take:
    ch_input      // parsed samplesheet channel
    ch_fasta
    ch_fai
    ch_bwt_index
    ch_dict
    ch_pon
    ch_intervals
    ch_annotsv_db
    ch_cytobands
    ch_eh_catalog

    main:
    // M1: Preprocess
    PREPROCESS(ch_input, ch_fasta, ch_fai, ch_bwt_index, params.min_depth)

    ch_bam = PREPROCESS.out.bam

    // M2 + M3 + M4: run in parallel on same BAM
    SV_CALLING(ch_bam, ch_fasta, ch_fai, ch_eh_catalog)
    CNV_CALLING(ch_bam, ch_fasta, ch_fai, ch_dict, ch_pon, ch_intervals)
    SMN_CALLING(ch_bam, ch_fasta, ch_fai)

    // M5: Annotate SVs
    ANNOTATE(SV_CALLING.out.sv_vcf, ch_annotsv_db)

    // Optional truvari truth channels
    ch_truth     = params.giab_truth
        ? Channel.fromPath(params.giab_truth, checkIfExists: true)
        : Channel.empty()
    ch_truth_tbi = params.giab_truth
        ? Channel.fromPath("${params.giab_truth}.tbi", checkIfExists: true)
        : Channel.empty()
    ch_truth_bed = params.giab_truth
        ? Channel.fromPath("${params.giab_truth}".replaceAll(/\.vcf\.gz$/, '.bed'), checkIfExists: true)
        : Channel.empty()

    // Collect QC files for MultiQC aggregation
    ch_multiqc_files = Channel.empty()
        .mix(PREPROCESS.out.fastqc_zip.map { meta, zips -> zips }.flatten())
        .mix(PREPROCESS.out.metrics.map { meta, m -> m })
        .mix(PREPROCESS.out.coverage.map { meta, s -> s })
        .collect()

    // M6 + M7: Visualize and report
    REPORT(
        ANNOTATE.out.tsv,
        CNV_CALLING.out.cnv_bed,
        SMN_CALLING.out.tsv,
        SV_CALLING.out.sv_vcf,
        SV_CALLING.out.sv_tbi,
        SV_CALLING.out.str_vcf,
        ch_cytobands,
        ch_truth,
        ch_truth_tbi,
        ch_truth_bed,
        ch_multiqc_files,
        PREPROCESS.out.coverage,
        PREPROCESS.out.metrics,
    )

    emit:
    sv_vcf   = SV_CALLING.out.sv_vcf
    str_vcf  = SV_CALLING.out.str_vcf
    cnv_bed  = CNV_CALLING.out.cnv_bed
    smn_tsv  = SMN_CALLING.out.tsv
    html     = REPORT.out.html
    multiqc  = REPORT.out.multiqc
}
