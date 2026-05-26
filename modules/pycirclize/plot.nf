process CIRCOS_PLOT {
    tag "${meta.id}"
    label 'process_single'
    container 'svcaller/utils:1.0'
    publishDir "${params.outdir}/${meta.id}", mode: 'copy', pattern: "*.circos.png"

    input:
    tuple val(meta), path(sv_vcf), path(cnv_bed), path(str_vcf)
    path cytobands

    output:
    tuple val(meta), path("${meta.id}.circos.svg"), emit: svg
    tuple val(meta), path("${meta.id}.circos.png"), emit: png

    script:
    def str_arg = str_vcf.name != "NO_STR" ? "--str-vcf ${str_vcf}" : ""
    // v2: filtered links (BND/TRA + DEL/DUP/INV>=50kb, multi-caller, cap 100); 150dpi PNG
    """
    export PATH=${projectDir}/bin:\$PATH
    circos_plot.py \\
        --sv-vcf    ${sv_vcf} \\
        --cnv-bed   ${cnv_bed} \\
        --cytobands ${cytobands} \\
        --sample    ${meta.id} \\
        ${str_arg} \\
        --out       ${meta.id}.circos.svg
    """
}
