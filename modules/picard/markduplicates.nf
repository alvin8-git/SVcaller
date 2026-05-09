process PICARD_MARKDUP {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.markdup.bam"),     emit: bam
    tuple val(meta), path("${meta.id}.markdup.bam.bai"), emit: bai
    tuple val(meta), path("${meta.id}.dup_metrics.txt"), emit: metrics
    path "versions.yml",                                  emit: versions

    script:
    """
    mkdir -p tmp
    picard MarkDuplicates \\
        -Xmx${(task.memory.toGiga() * 0.85).intValue()}g \\
        INPUT=${bam} \\
        OUTPUT=${meta.id}.markdup.bam \\
        METRICS_FILE=${meta.id}.dup_metrics.txt \\
        REMOVE_DUPLICATES=false \\
        VALIDATION_STRINGENCY=LENIENT \\
        CREATE_INDEX=true \\
        TMP_DIR=./tmp

    mv ${meta.id}.markdup.bai ${meta.id}.markdup.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard MarkDuplicates --version 2>&1 | grep -o 'Version:.*' | sed 's/Version://')
    END_VERSIONS
    """
}
