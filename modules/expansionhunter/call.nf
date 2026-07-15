process EXPANSIONHUNTER {
    tag "${meta.id}"
    label 'process_low'

    // Publish the repeat-expansion profile to a stable results path so downstream
    // consumers (e.g. OmniGen) read it from results/ instead of the ephemeral
    // Nextflow work/ hash dir, which is garbage-collected and changes every run.
    publishDir "${params.outdir}/${meta.id}", mode: 'copy', pattern: "*.str_profile.json"

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

    # Normalise output json filename across EH versions.
    # An empty '{}' profile is NOT a valid "no repeat expansions" result -- it is what a
    # broken ExpansionHunter run looks like. This file is published to results/ and read
    # by OmniGen, which previously could not tell the two apart and rendered a crashed
    # caller as a clean report. If EH exited 0 without a profile, fail loudly.
    if [ -f ${meta.id}.str.json ]; then
        cp ${meta.id}.str.json ${meta.id}.str_profile.json
    else
        echo "ERROR: ExpansionHunter exited 0 but wrote no ${meta.id}.str.json profile." >&2
        echo "Refusing to publish a fake empty '{}' STR profile: downstream consumers" >&2
        echo "would read it as 'no repeat expansions detected'." >&2
        ls -la >&2
        exit 1
    fi

    if [ ! -s ${meta.id}.str_profile.json ]; then
        echo "ERROR: ${meta.id}.str_profile.json is empty." >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        expansionhunter: \$(ExpansionHunter --version 2>&1 | head -1 | awk '{print \$NF}')
    END_VERSIONS
    """
}
