process GRIDSS_CONVERT_BND {
    tag "${meta.id}"
    label 'process_single'
    container params.utils_container ?: 'svcaller/utils:1.3'

    input:
    tuple val(meta), path(bnd_vcf), path(bnd_tbi)

    output:
    tuple val(meta), path("${meta.id}.gridss.simple.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.gridss.simple.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                                            emit: versions

    script:
    """
    export PATH=${projectDir}/bin:\$PATH

    gridss_convert_bnd.py \\
        ${bnd_vcf} \\
        --out ${meta.id}.gridss.simple.vcf \\
        --min-qual-del  500 \\
        --min-qual-dup 1000 \\
        --min-qual-inv 1000 \\
        --min-svlen 50

    bgzip ${meta.id}.gridss.simple.vcf
    tabix -p vcf ${meta.id}.gridss.simple.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss_convert_bnd: 1.0
    END_VERSIONS
    """
}
