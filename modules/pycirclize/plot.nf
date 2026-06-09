process CIRCOS_PLOT {
    tag "${meta.id}"
    label 'process_single'
    container 'svcaller/utils:1.0'
    publishDir "${params.outdir}/${meta.id}", mode: 'copy', pattern: "*.circos.png"

    input:
    tuple val(meta), path(sv_vcf), path(cnv_bed), path(str_vcf), path(depth_bed), path(annotsv_tsv), path(strling_tsv)
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
    // v10: CNV ring (58-61), STR ring EH+STRling (54-57), depth (62-92)
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
        --out       ${meta.id}.circos.svg
    """
}
