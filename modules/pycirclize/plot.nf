process CIRCOS_PLOT {
    tag "${meta.id}"
    label 'process_single'
    container 'svcaller/utils:1.0'

    input:
    tuple val(meta), path(sv_vcf), path(cnv_bed)
    path cytobands

    output:
    tuple val(meta), path("${meta.id}.circos.svg"), emit: svg
    tuple val(meta), path("${meta.id}.circos.png"), emit: png

    script:
    """
    circos_plot.py \\
        --sv-vcf    ${sv_vcf} \\
        --cnv-bed   ${cnv_bed} \\
        --cytobands ${cytobands} \\
        --sample    ${meta.id} \\
        --out       ${meta.id}.circos.svg
    """
}
