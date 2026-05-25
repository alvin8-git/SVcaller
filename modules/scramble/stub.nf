process SCRAMBLE_STUB {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.scramble.stub.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.scramble.stub.vcf.gz.tbi"), emit: tbi

    script:
    """
    printf '##fileformat=VCFv4.2\\n##FILTER=<ID=PASS,Description="All filters passed">\\n##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\\n##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length">\\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t${meta.id}\\n' \\
        | bgzip > ${meta.id}.scramble.stub.vcf.gz
    tabix -p vcf ${meta.id}.scramble.stub.vcf.gz
    """
}
