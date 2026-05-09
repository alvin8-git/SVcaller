process BWAMEM2_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(reads)
    path fasta
    path fai
    path bwt_index  // directory containing bwa-mem2 index files

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"),     emit: bam
    tuple val(meta), path("${meta.id}.sorted.bam.bai"), emit: bai
    path "versions.yml",                                 emit: versions

    script:
    def rg = "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:${meta.id}"
    """
    bwa-mem2 mem \\
        -t ${task.cpus} \\
        -R "${rg}" \\
        ${fasta} \\
        ${reads} \\
        | samtools sort -@ ${task.cpus} -m 2G \\
        -o ${meta.id}.sorted.bam

    samtools index -@ ${task.cpus} ${meta.id}.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>/dev/null | head -1 | tr -d '\\n')
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
