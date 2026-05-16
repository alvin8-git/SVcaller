process PICARD_INSERT_SIZE {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.insert_size_metrics.txt"), emit: metrics
    path "versions.yml",                                          emit: versions

    script:
    def avail_mem = task.memory ? "${(task.memory.giga * 0.8).intValue()}g" : "4g"
    """
    picard -Xmx${avail_mem} CollectInsertSizeMetrics \\
        -I ${bam} \\
        -O ${meta.id}.insert_size_metrics.txt \\
        -H ${meta.id}.insert_size_histogram.pdf \\
        --MINIMUM_PCT 0.05

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CollectInsertSizeMetrics --version 2>&1 | grep -oE '[0-9]+\\.[0-9.]+' | head -1)
    END_VERSIONS
    """
}
