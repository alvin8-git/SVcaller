process MULTIQC {
    label 'process_single'
    container 'quay.io/biocontainers/multiqc:1.34--pyhdfd78af_0'

    input:
    path multiqc_files, stageAs: "?/*"

    output:
    path "multiqc_report.html", emit: html
    path "multiqc_data/",       emit: data
    path "versions.yml",        emit: versions

    script:
    """
    multiqc --force .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version 2>&1 | awk '{print \$NF}')
    END_VERSIONS
    """
}
