process FASTQC {
    tag "${meta.id}"
    label 'process_low'
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"),  emit: zip
    path  "versions.yml",            emit: versions

    script:
    """
    fastqc --threads ${task.cpus} ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version 2>&1 | grep FastQC | awk '{print \$NF}')
    END_VERSIONS
    """
}
