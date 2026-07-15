process SMN_CALLER {
    tag "${meta.id}"
    label 'process_low'
    container 'svcaller/smncopynum:1.1'
    publishDir "${params.outdir}/${meta.id}", mode: 'copy', pattern: "*.smn.tsv"

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

    # Rename output files to standard names (tool emits {prefix}.tsv / {prefix}.json).
    # NEVER fall back to 'touch': an empty .smn.tsv is indistinguishable from a real
    # "no SMN findings" result downstream, and consumers (e.g. OmniGen) rendered a
    # crashed SMN caller as a clean bill of health. If the caller produced nothing,
    # this process MUST fail.
    mv ${meta.id}_smn.tsv ${meta.id}.smn.tsv 2>/dev/null || \\
        mv ${meta.id}.tsv   ${meta.id}.smn.tsv 2>/dev/null || \\
        mv smn_result*.tsv  ${meta.id}.smn.tsv 2>/dev/null || true

    if [ ! -s "${meta.id}.smn.tsv" ]; then
        echo "ERROR: SMN_CALLER produced no non-empty TSV for sample '${meta.id}'." >&2
        echo "smn_caller.py exited 0 but no result table was found (or it was empty)." >&2
        echo "Refusing to publish a zero-byte placeholder. Work dir contents:" >&2
        ls -la >&2
        exit 1
    fi

    # A single-BAM manifest must yield exactly one result row. A header-only table
    # means the caller did not genotype the sample -> that is a failure, not a result.
    # Count non-blank lines: a valid table is header + >=1 sample row = >=2 lines.
    n_lines=\$(grep -c '[^[:space:]]' ${meta.id}.smn.tsv || true)
    if [ "\${n_lines:-0}" -lt 2 ]; then
        echo "ERROR: SMN_CALLER TSV for '${meta.id}' has no sample row (header only)." >&2
        echo "The caller did not genotype this sample. Refusing to publish an empty result." >&2
        cat ${meta.id}.smn.tsv >&2
        exit 1
    fi

    mv ${meta.id}_smn_detail.json ${meta.id}.smn_detail.json 2>/dev/null || \\
        mv ${meta.id}.json          ${meta.id}.smn_detail.json 2>/dev/null || true

    if [ ! -s "${meta.id}.smn_detail.json" ]; then
        echo "ERROR: SMN_CALLER produced no detail JSON for sample '${meta.id}'." >&2
        echo "smn_caller.py emits {prefix}.tsv and {prefix}.json together; a missing JSON" >&2
        echo "alongside a present TSV means the run is incomplete. Not writing a fake '{}'." >&2
        ls -la >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        smncopynum: 1.1
    END_VERSIONS
    """
}
