process CIRCOS_PLOT {
    tag "${meta.id}"
    label 'process_single'
    container 'svcaller/utils:1.0'
    publishDir "${params.outdir}/${meta.id}", mode: 'copy', pattern: "*.circos.png"

    input:
    tuple val(meta), path(sv_vcf), path(cnv_bed), path(str_vcf), path(depth_bed), path(annotsv_tsv)
    path cytobands

    output:
    tuple val(meta), path("${meta.id}.circos.svg"), emit: svg
    tuple val(meta), path("${meta.id}.circos.png"), emit: png

    script:
    def str_arg     = str_vcf.name    != "NO_STR"  ? "--str-vcf    ${str_vcf}"    : ""
    def depth_arg   = depth_bed.name  != "NO_FILE" ? "--depth-bed  ${depth_bed}"  : ""
    def annotsv_arg = annotsv_tsv.name != "NO_FILE" ? "--annotsv-tsv ${annotsv_tsv}" : ""
    // v6: remove CNV rings; expand depth to log2 scatter (63-93); update legend
    """
    export PATH=${projectDir}/bin:\$PATH
    circos_plot.py \\
        --sv-vcf    ${sv_vcf} \\
        --cnv-bed   ${cnv_bed} \\
        --cytobands ${cytobands} \\
        --sample    ${meta.id} \\
        ${str_arg} \\
        ${depth_arg} \\
        ${annotsv_arg} \\
        --out       ${meta.id}.circos.svg
    """
}
