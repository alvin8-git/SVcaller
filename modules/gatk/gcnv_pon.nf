process GATK_PREPROCESS_INTERVALS {
    label 'process_single'
    container 'broadinstitute/gatk:4.5.0.0'
    publishDir "${params.outdir}/pon", mode: 'copy', pattern: "*.interval_list"

    input:
    path fasta
    path fai
    path dict
    path intervals

    output:
    path "preprocessed.interval_list", emit: preprocessed
    path "versions.yml",                emit: versions

    script:
    def heap = (task.memory.toGiga() * 0.85).intValue()
    """
    gatk --java-options "-Xmx${heap}g" PreprocessIntervals \\
        -R ${fasta} \\
        -L ${intervals} \\
        --bin-length 1000 \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        -O preprocessed.interval_list

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep -oP '(?<=GATK v)[0-9.]+' | head -1)
    END_VERSIONS
    """
}

process GATK_COLLECT_COUNTS {
    tag "${meta.id}"
    label 'process_medium'
    container 'broadinstitute/gatk:4.5.0.0'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai
    path fasta_dict
    path intervals

    output:
    tuple val(meta), path("${meta.id}.counts.hdf5"), emit: hdf5
    path "versions.yml",                              emit: versions

    script:
    def heap = (task.memory.toGiga() * 0.85).intValue()
    """
    gatk --java-options "-Xmx${heap}g" CollectReadCounts \\
        -I ${bam} \\
        -L ${intervals} \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        -R ${fasta} \\
        -O ${meta.id}.counts.hdf5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep -oP '(?<=GATK v)[0-9.]+' | head -1)
    END_VERSIONS
    """
}

process GATK_CREATE_PON {
    label 'process_high'
    container 'broadinstitute/gatk:4.5.0.0'
    publishDir "${params.outdir}/pon", mode: 'copy'

    input:
    path hdf5_files   // list of all sample HDF5 count files
    path annotated_intervals

    output:
    path "giab_cnv_pon.hdf5", emit: pon
    path "versions.yml",       emit: versions

    script:
    def heap = (task.memory.toGiga() * 0.85).intValue()
    def inputs = hdf5_files.collect { "-I $it" }.join(" \\\n        ")
    """
    gatk --java-options "-Xmx${heap}g" CreateReadCountPanelOfNormals \\
        ${inputs} \\
        --annotated-intervals ${annotated_intervals} \\
        --output giab_cnv_pon.hdf5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep -oP '(?<=GATK v)[0-9.]+' | head -1)
    END_VERSIONS
    """
}
