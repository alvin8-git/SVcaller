// Read-depth consistency annotation for DEL/DUP calls (bin/depth_check.py).
// Last step before <sample>.filtered.tsv is published. Annotation only: no row is
// dropped, seven columns (DEPTH_IN, DEPTH_FLANK, DEPTH_GLOBAL, DEPTH_RATIO, DEPTH_RATIO_GLOBAL, DEPTH_SUPPORT, DEPTH_NOTE) are appended.
// Container: default svcaller/utils (python3 + samtools); no new dependency.
process DEPTH_CHECK {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/${meta.id}", mode: 'copy', pattern: "*.filtered.tsv"

    input:
    tuple val(meta), path(gnomad_tsv), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.filtered.tsv"), emit: tsv

    script:
    """
    depth_check.py ${gnomad_tsv} ${bam} ${meta.id}.filtered.tsv --threads ${task.cpus}
    """
}
