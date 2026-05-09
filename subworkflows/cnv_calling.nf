include { CNVPYTOR_CALL  } from '../modules/cnvpytor/call'
include { GATK_GCNV_CALL } from '../modules/gatk/gcnv_call'

process CNV_CONSENSUS {
    tag "${meta.id}"
    label 'process_single'
    container params.utils_container ?: 'svcaller/utils:1.0'

    input:
    tuple val(meta), path(cnvpytor_tsv), path(gatk_tsv)

    output:
    tuple val(meta), path("${meta.id}.cnv_consensus.bed"), emit: bed

    script:
    """
    cnv_consensus.py \\
        --cnvpytor ${cnvpytor_tsv} \\
        --gatk     ${gatk_tsv} \\
        --sample   ${meta.id} \\
        --out      ${meta.id}.cnv_consensus.bed
    """
}

workflow CNV_CALLING {
    take:
    ch_bam       // [ meta, bam, bai ]
    ch_fasta     // path
    ch_fai       // path
    ch_dict      // path (.dict)
    ch_pon       // path (PoN HDF5, may be null)
    ch_intervals // path

    main:
    CNVPYTOR_CALL(ch_bam, ch_fasta)
    GATK_GCNV_CALL(ch_bam, ch_fasta, ch_fai, ch_dict, ch_pon, ch_intervals)

    ch_for_consensus = CNVPYTOR_CALL.out.tsv
        .join(GATK_GCNV_CALL.out.seg.map { meta, seg -> [meta, seg] })

    CNV_CONSENSUS(ch_for_consensus)

    emit:
    cnv_bed = CNV_CONSENSUS.out.bed
}
