process EXPANSIONHUNTER {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai
    path catalog

    output:
    tuple val(meta), path("${meta.id}.str.vcf.gz"),       emit: vcf
    tuple val(meta), path("${meta.id}.str.vcf.gz.tbi"),   emit: tbi
    tuple val(meta), path("${meta.id}.str_profile.json"), emit: json
    path "versions.yml",                                   emit: versions

    script:
    """
    ExpansionHunter \\
        --reads ${bam} \\
        --reference ${fasta} \\
        --variant-catalog ${catalog} \\
        --output-prefix ${meta.id}.str \\
        --threads ${task.cpus}

    bgzip ${meta.id}.str.vcf
    tabix -p vcf ${meta.id}.str.vcf.gz

    # Normalise output json filename across EH versions
    if [ -f ${meta.id}.str.json ]; then
        cp ${meta.id}.str.json ${meta.id}.str_profile.json
    else
        echo '{}' > ${meta.id}.str_profile.json
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        expansionhunter: \$(ExpansionHunter --version 2>&1 | head -1 | awk '{print \$NF}')
    END_VERSIONS
    """
}
