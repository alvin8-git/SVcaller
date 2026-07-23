process MOSDEPTH {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/${meta.id}/qc", mode: 'copy', pattern: "*.mosdepth.summary.txt"

    input:
    tuple val(meta), path(bam), path(bai)
    val   min_depth

    output:
    tuple val(meta), path("${meta.id}.mosdepth.summary.txt"),  emit: summary
    tuple val(meta), path("${meta.id}.regions.bed.gz"),        emit: regions_bed
    path "versions.yml",                                        emit: versions

    script:
    """
    # --no-per-base: skip the ~4.5 GB/sample per-base.bed.gz. It is never declared as
    # an output nor read downstream (only summary.txt + regions.bed.gz are used), so
    # generating it only wastes disk and time (per-base is mosdepth's slowest stage).
    mosdepth \\
        --threads ${task.cpus} \\
        --no-per-base \\
        --quantize 0:5:30:500: \\
        --by 50000 \\
        ${meta.id} \\
        ${bam}

    # Fail pipeline if mean depth below threshold.
    # awk /^total/{exit} stops at first match — avoids multi-line MEAN_DEPTH when --by is used.
    MEAN_DEPTH=\$(awk '/^total/{print \$4; exit}' ${meta.id}.mosdepth.summary.txt)
    if awk -v d="\$MEAN_DEPTH" -v m="${min_depth}" 'BEGIN{exit (d >= m) ? 0 : 1}'; then
        echo "PASS: mean depth \$MEAN_DEPTH >= ${min_depth}x"
    else
        echo "ERROR: mean depth \$MEAN_DEPTH < ${min_depth}x (required). Aborting." >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/mosdepth //')
    END_VERSIONS
    """
}
