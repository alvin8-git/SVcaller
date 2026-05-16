process SAMTOOLS_FILTER_CHROMS {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)
    path fai

    output:
    tuple val(meta), path("${meta.id}.filtered.bam"), path("${meta.id}.filtered.bam.bai"), emit: bam

    script:
    """
    REF_CHROMS=\$(awk '{print \$1}' ${fai} | tr '\\n' ' ')
    samtools view -h -@ ${task.cpus} -b -o ${meta.id}.filtered.bam ${bam} \$REF_CHROMS
    samtools index ${meta.id}.filtered.bam
    """
}
