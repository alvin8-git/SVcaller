// Depth floor for the SMN1/SMN2 result (bin/smn_depth_floor.py). Below params.smn_min_depth
// mean coverage the caller's copy numbers are replaced by an explicit no-call
// (Info=NoCall:BelowValidatedDepth); the caller's own row is kept as <id>.smn.raw.tsv.
// Mean depth is the value MOSDEPTH already measured; nothing is re-measured here.
// The input is staged under in/ so --out never writes through the staged symlink into SMN_CALLER's cache.
// Container: svcaller/utils (python3 only).
process SMN_DEPTH_FLOOR {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/${meta.id}", mode: 'copy', pattern: "*.smn.tsv"
    publishDir "${params.outdir}/${meta.id}", mode: 'copy', pattern: "*.smn.raw.tsv"

    input:
    tuple val(meta), path(smn_tsv, stageAs: 'in/*'), path(mosdepth_summary)
    val   min_depth

    output:
    tuple val(meta), path("${meta.id}.smn.tsv"), emit: tsv
    path "${meta.id}.smn.raw.tsv", optional: true, emit: raw

    script:
    """
    smn_depth_floor.py in/${smn_tsv.name} \\
        --summary ${mosdepth_summary} \\
        --min-depth ${min_depth} \\
        --out ${meta.id}.smn.tsv \\
        --raw ${meta.id}.smn.raw.tsv
    """
}
