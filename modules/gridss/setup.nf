process GRIDSS_SETUP {
    tag "genome_index"
    label 'process_gridss'
    storeDir "${params.outdir}/cache/gridss_ref"

    input:
    path fasta
    path fai

    output:
    path "${fasta}.amb", emit: amb
    path "${fasta}.ann", emit: ann
    path "${fasta}.bwt", emit: bwt
    path "${fasta}.pac", emit: pac
    path "${fasta}.sa",  emit: sa

    script:
    def heap = (task.memory.toGiga() * 0.85).intValue()
    """
    gridss \\
        --reference ${fasta} \\
        --steps setupreference \\
        --threads ${task.cpus} \\
        --jvmheap ${heap}g \\
        --workingdir ./ref_cache
    """
}
