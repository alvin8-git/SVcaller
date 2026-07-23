process SAMTOOLS_FILTER_CHROMS {
    tag "${meta.id}"
    label 'process_medium'
    // Cache filtered BAM across sessions — avoids re-writing 130 GB on every new run.
    // Stored in outdir/.cache/filter_chroms; survives nextflow clean; delete manually when input BAM changes.
    storeDir "${params.outdir}/.cache/filter_chroms"

    input:
    tuple val(meta), path(bam), path(bai)
    path fai

    output:
    tuple val(meta), path("${meta.id}.filtered.bam"), path("${meta.id}.filtered.bam.bai"), emit: bam

    script:
    // The filter logic lives in bin/filter_chroms.sh (auto-staged onto PATH by
    // Nextflow). It fans the per-read work over the 24 canonical chromosomes with
    // xargs -P${task.cpus}, then samtools cat, instead of a single-threaded awk
    // that pinned one core for ~70 min. Output is content-identical (same reads,
    // same header, canonical-sorted); see the script header and tests/test_filter_chroms.py.
    // Only canonical @SQ are emitted (Manta assembles 0 variants when the header
    // keeps all 3366 hg38 alt-contig @SQ), and reads whose mate is on a non-canonical
    // contig are dropped (Manta FATAL_ERROR on mate_tid=-1).
    """
    filter_chroms.sh ${bam} ${fai} ${task.cpus} ${meta.id}.filtered.bam
    """
}
