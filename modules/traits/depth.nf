process TRAIT_DEPTH {
    tag "${meta.id}"
    label 'process_low'
    // container assigned in conf/docker.config: mosdepth 0.3.14 biocontainer

    input:
    tuple val(meta), path(bam), path(bai)
    path  regions_bed          // assets/cnv_trait_regions.bed (trait + CTRL_* control regions)

    output:
    tuple val(meta), path("${meta.id}.trait_depth.regions.bed.gz"), emit: depth
    path "versions.yml",                                            emit: versions

    script:
    """
    # --mapq 0 (do NOT MAPQ-filter): AMY1 array copies, LPA KIV-2 repeat units and
    # paralogous RHD/GST reads are multi-mapping; filtering destroys the copy-number
    # signal. Normalizing against the CTRL_* control regions (also --mapq 0) cancels
    # the multi-mapping bias to first order.
    mosdepth \\
        --threads ${task.cpus} \\
        --no-per-base \\
        --mapq 0 \\
        --by ${regions_bed} \\
        ${meta.id}.trait_depth \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/mosdepth //')
    END_VERSIONS
    """
}
