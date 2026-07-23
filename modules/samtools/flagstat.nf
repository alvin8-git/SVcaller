process SAMTOOLS_FLAGSTAT {
    tag "${meta.id}"
    label 'process_single'
    // Publish so the report's QC section is reproducible after the work dir is
    // cleaned (Rule 5). Consumed from work/ only used to silently lose QC.
    publishDir "${params.outdir}/${meta.id}/qc", mode: 'copy', pattern: "*.flagstat.txt"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.flagstat.txt"), emit: flagstat
    path "versions.yml",                               emit: versions

    script:
    """
    samtools flagstat -@ ${task.cpus} ${bam} > ${meta.id}.flagstat.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
