process SAMTOOLS_SUBSET {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai), path(bed)

    output:
    tuple val(meta), path("${meta.id}.residual.bam"), path("${meta.id}.residual.bam.bai"), emit: bam
    path "versions.yml",                                                                     emit: versions

    script:
    """
    # --fetch-pairs: include complete read pairs even when mate is outside BED regions
    # This is critical for GRIDSS which relies on discordant pair orientation
    samtools view \\
        -L ${bed} \\
        --fetch-pairs \\
        -@ ${task.cpus} \\
        -b \\
        -o ${meta.id}.residual.bam \\
        ${bam}
    samtools index -@ ${task.cpus} ${meta.id}.residual.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | awk '{print \$2}')
    END_VERSIONS
    """
}
