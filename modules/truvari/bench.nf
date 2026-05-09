process TRUVARI_BENCH {
    tag "${meta.id}"
    label 'process_low'
    container 'quay.io/biocontainers/truvari:4.2.2--pyhdfd78af_0'

    input:
    tuple val(meta), path(query_vcf), path(query_tbi)
    path truth_vcf
    path truth_tbi
    path truth_bed

    output:
    tuple val(meta), path("${meta.id}.truvari/summary.json"), emit: summary
    tuple val(meta), path("${meta.id}.truvari/"),             emit: dir

    script:
    """
    truvari bench \\
        -b ${truth_vcf} \\
        -c ${query_vcf} \\
        --includebed ${truth_bed} \\
        -o ${meta.id}.truvari \\
        --passonly \\
        --pick multi \\
        --sizemin 50
    """
}
