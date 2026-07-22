include { TRAIT_DEPTH } from '../modules/traits/depth'

// One targeted read-depth pass (TRAIT_DEPTH) feeds four small interpreter
// processes, each emitting a stable per-sample contract file consumed by OmniGen.
// The CNV consensus BED is passed in only as a corroborating signal.

process RH_STATUS {
    tag "${meta.id}"
    label 'process_single'
    container params.utils_container ?: 'svcaller/utils:1.3'
    publishDir "${params.outdir}/${meta.id}/bloodgroup", mode: 'copy', pattern: "*.rh_status.tsv"

    input:
    tuple val(meta), path(depth_bed), path(cnv_bed)

    output:
    tuple val(meta), path("${meta.id}.rh_status.tsv"), emit: tsv

    script:
    """
    export PATH=${projectDir}/bin:\$PATH
    rh_status.py \\
        --depth   ${depth_bed} \\
        --cnv-bed ${cnv_bed} \\
        --sample  ${meta.id} \\
        --out     ${meta.id}.rh_status.tsv
    """
}

process AMY1_CN {
    tag "${meta.id}"
    label 'process_single'
    container params.utils_container ?: 'svcaller/utils:1.3'
    publishDir "${params.outdir}/${meta.id}/cnv_traits", mode: 'copy', pattern: "*.amy1.tsv"

    input:
    tuple val(meta), path(depth_bed), path(cnv_bed)

    output:
    tuple val(meta), path("${meta.id}.amy1.tsv"), emit: tsv

    script:
    """
    export PATH=${projectDir}/bin:\$PATH
    amy1_cn.py \\
        --depth  ${depth_bed} \\
        --sample ${meta.id} \\
        --out    ${meta.id}.amy1.tsv
    """
}

process GST_NULL {
    tag "${meta.id}"
    label 'process_single'
    container params.utils_container ?: 'svcaller/utils:1.3'
    publishDir "${params.outdir}/${meta.id}/cnv_traits", mode: 'copy', pattern: "*.gst_null.tsv"

    input:
    tuple val(meta), path(depth_bed), path(cnv_bed)

    output:
    tuple val(meta), path("${meta.id}.gst_null.tsv"), emit: tsv

    script:
    """
    export PATH=${projectDir}/bin:\$PATH
    gst_null.py \\
        --depth   ${depth_bed} \\
        --cnv-bed ${cnv_bed} \\
        --sample  ${meta.id} \\
        --out     ${meta.id}.gst_null.tsv
    """
}

process LPA_KIV2 {
    tag "${meta.id}"
    label 'process_single'
    container params.utils_container ?: 'svcaller/utils:1.3'
    publishDir "${params.outdir}/${meta.id}/cnv_traits", mode: 'copy', pattern: "*.lpa_kiv2.tsv"

    input:
    tuple val(meta), path(depth_bed), path(cnv_bed)

    output:
    tuple val(meta), path("${meta.id}.lpa_kiv2.tsv"), emit: tsv

    script:
    """
    export PATH=${projectDir}/bin:\$PATH
    lpa_kiv2.py \\
        --depth  ${depth_bed} \\
        --sample ${meta.id} \\
        --out    ${meta.id}.lpa_kiv2.tsv
    """
}

workflow CNV_TRAITS {
    take:
    ch_bam           // [ meta, bam, bai ]
    ch_trait_regions // path  assets/cnv_trait_regions.bed
    ch_cnv_bed       // [ meta, cnv_consensus.bed ]  from CNV_CALLING.out.cnv_bed

    main:
    TRAIT_DEPTH(ch_bam, ch_trait_regions)

    // [ meta, depth_bed, cnv_bed ] — cnv_bed optional (NO_FILE sentinel keeps the
    // interpreters running for samples without a consensus BED).
    ch_in = TRAIT_DEPTH.out.depth
        .join(ch_cnv_bed, remainder: true)
        .filter { it[1] != null }   // drop cnv_bed-only remainders (depth absent)
        .map { meta, depth, cnv -> [meta, depth, cnv ?: file("NO_FILE")] }

    RH_STATUS(ch_in)
    AMY1_CN(ch_in)
    GST_NULL(ch_in)
    LPA_KIV2(ch_in)

    emit:
    rh   = RH_STATUS.out.tsv
    amy1 = AMY1_CN.out.tsv
    gst  = GST_NULL.out.tsv
    lpa  = LPA_KIV2.out.tsv
}
