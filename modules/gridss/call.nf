process GRIDSS_CALL {
    tag "${meta.id}"
    label 'process_gridss'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path("${meta.id}.gridss.sv.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.gridss.sv.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                                        emit: versions

    script:
    def heap = (task.memory.toGiga() * 0.85).intValue()
    """
    gridss \\
        --reference ${fasta} \\
        --output ${meta.id}.gridss.sv.vcf.gz \\
        --workingdir ./gridss_work \\
        --threads ${task.cpus} \\
        --jvmheap ${heap}g \\
        --picardoptions VALIDATION_STRINGENCY=LENIENT \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: \$(gridss --version 2>&1 | grep -oP '(?<=GRIDSS v)[^ ]+' | head -1)
    END_VERSIONS
    """
}
