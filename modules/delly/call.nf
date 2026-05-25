// DELLY_CALL_SVTYPE runs one SV type per process invocation.
// sv_calling.nf fans out 5 types in parallel then collects into DELLY_MERGE.
process DELLY_CALL_SVTYPE {
    tag "${meta.id}:${svtype}"
    label 'process_medium'
    maxForks 4

    input:
    tuple val(meta), path(bam), path(bai), val(svtype)
    path fasta
    path fai

    output:
    tuple val(meta), path("${meta.id}.delly.${svtype}.vcf"), emit: vcf
    path "versions.yml",                                       emit: versions

    script:
    """
    delly call \\
        -t ${svtype} \\
        -g ${fasta} \\
        -o ${meta.id}.delly.${svtype}.vcf \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$(delly --version 2>&1 | grep "DELLY" | head -1 | awk '{print \$2}')
    END_VERSIONS
    """
}
