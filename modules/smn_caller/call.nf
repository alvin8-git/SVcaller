process SMN_CALLER {
    tag "${meta.id}"
    label 'process_low'
    container 'svcaller/smncopynum:1.1'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path("${meta.id}.smn.tsv"),         emit: tsv
    tuple val(meta), path("${meta.id}.smn_detail.json"), emit: json
    path "versions.yml",                                  emit: versions

    script:
    """
    # smn_caller.py takes a manifest (one BAM path per line)
    echo "\${PWD}/${bam}" > manifest.txt

    python /opt/smn/smn_caller.py \\
        --manifest manifest.txt \\
        --genome 38 \\
        --outDir . \\
        --prefix ${meta.id} \\
        --threads ${task.cpus} \\
        --reference ${fasta}

    # Rename output files to standard names
    mv ${meta.id}_smn.tsv ${meta.id}.smn.tsv 2>/dev/null || \\
        mv smn_result*.tsv ${meta.id}.smn.tsv 2>/dev/null || \\
        touch ${meta.id}.smn.tsv

    mv ${meta.id}_smn_detail.json ${meta.id}.smn_detail.json 2>/dev/null || \\
        echo '{}' > ${meta.id}.smn_detail.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        smncopynum: 1.1
    END_VERSIONS
    """
}
