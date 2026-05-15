process SAMTOOLS_SORT {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"),     emit: bam
    tuple val(meta), path("${meta.id}.sorted.bam.bai"), emit: bai
    path "versions.yml",                                 emit: versions

    script:
    """
    samtools sort -@ ${task.cpus} -m 2G \\
        -o ${meta.id}.sorted.bam ${bam}

    samtools index -@ ${task.cpus} ${meta.id}.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
