process CNVPYTOR_CALL {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta

    output:
    tuple val(meta), path("${meta.id}.cnvpytor.tsv"), emit: tsv
    path "versions.yml",                               emit: versions

    script:
    """
    cnvpytor -root ${meta.id}.pytor -rd ${bam}
    cnvpytor -root ${meta.id}.pytor -his 1000 10000 100000
    cnvpytor -root ${meta.id}.pytor -partition 1000 10000 100000
    cnvpytor -root ${meta.id}.pytor -call 1000  > ${meta.id}.cnvpytor_1kb.tsv
    cnvpytor -root ${meta.id}.pytor -call 10000 > ${meta.id}.cnvpytor_10kb.tsv

    # Merge 1kb and 10kb calls
    cat ${meta.id}.cnvpytor_1kb.tsv ${meta.id}.cnvpytor_10kb.tsv \\
        | sort -k2,2 -k3,3n \\
        > ${meta.id}.cnvpytor.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvpytor: \$(cnvpytor --version 2>&1 | head -1 | awk '{print \$NF}')
    END_VERSIONS
    """
}
