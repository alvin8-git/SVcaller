process MELT_STUB {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.melt.vcf.gz"), emit: vcf

    script:
    """
    printf "##fileformat=VCFv4.1\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t${meta.id}\\n" | gzip > ${meta.id}.melt.vcf.gz
    """
}
