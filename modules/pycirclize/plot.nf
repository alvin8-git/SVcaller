process CIRCOS_PLOT {
    tag "${meta.id}"
    label 'process_single'
    container params.utils_container ?: 'svcaller/utils:1.3'
    publishDir "${params.outdir}/${meta.id}", mode: 'copy', pattern: "*.circos.png"

    input:
    tuple val(meta), path(sv_vcf), path(cnv_bed), path(str_vcf), path(depth_bed), path(annotsv_tsv), path(strling_tsv), path(rh_status_tsv), path(amy1_tsv), path(gst_null_tsv), path(lpa_kiv2_tsv)
    path cytobands

    output:
    tuple val(meta), path("${meta.id}.circos.svg"), emit: svg
    tuple val(meta), path("${meta.id}.circos.png"), emit: png

    script:
    def cnv_arg      = cnv_bed.name     != "NO_FILE"     ? "--cnv-bed     ${cnv_bed}"     : ""
    def str_arg      = str_vcf.name     != "NO_STR"      ? "--str-vcf     ${str_vcf}"     : ""
    def depth_arg    = depth_bed.name   != "NO_DEPTH"    ? "--depth-bed   ${depth_bed}"   : ""
    def annotsv_arg  = annotsv_tsv.name != "NO_ANNOTSV"  ? "--annotsv-tsv ${annotsv_tsv}" : ""
    def strling_arg  = strling_tsv.name != "NO_STRLING"  ? "--strling-tsv ${strling_tsv}" : ""
    // CNV-trait ring inputs (optional): the four trait contract TSVs. Absent → ring skipped.
    def rh_arg       = rh_status_tsv.name != "NO_FILE"   ? "--rh-status   ${rh_status_tsv}" : ""
    def amy1_arg     = amy1_tsv.name      != "NO_FILE"   ? "--amy1        ${amy1_tsv}"      : ""
    def gst_arg      = gst_null_tsv.name  != "NO_FILE"   ? "--gst-null    ${gst_null_tsv}"  : ""
    def lpa_arg      = lpa_kiv2_tsv.name  != "NO_FILE"   ? "--lpa-kiv2    ${lpa_kiv2_tsv}"  : ""
    """
    export PATH=${projectDir}/bin:\$PATH
    circos_plot.py \\
        --sv-vcf    ${sv_vcf} \\
        --cytobands ${cytobands} \\
        --sample    ${meta.id} \\
        ${cnv_arg} \\
        ${str_arg} \\
        ${depth_arg} \\
        ${annotsv_arg} \\
        ${strling_arg} \\
        ${rh_arg} \\
        ${amy1_arg} \\
        ${gst_arg} \\
        ${lpa_arg} \\
        --out       ${meta.id}.circos.svg
    """
}
