process STRLING_CALL {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/strling:0.5.2--py39h3e45d22_2'

    input:
    tuple val(meta), path(bam), path(bai)
    path ref_fasta
    path ref_fai

    output:
    tuple val(meta), path("${meta.id}.strling.tsv"), emit: tsv

    script:
    """
    # Index reference for STRling (produces ref.str)
    strling index ${ref_fasta}

    # Extract STR-spanning reads per sample
    strling extract \\
        -f ${ref_fasta} \\
        -b ${ref_fasta}.str \\
        ${bam} \\
        ${meta.id}

    # Call STR expansions genome-wide
    strling call \\
        -f ${ref_fasta} \\
        ${meta.id}.txt \\
        ${bam}

    # Rename genotype output
    mv ${meta.id}-genotype.txt ${meta.id}.strling.tsv || \\
        printf "#chrom\\tleft\\tright\\trepeatunit\\tallele1_est\\tallele2_est\\tspanning\\tflanking\\tsum_str_counts\\tprob_expansion\\n" > ${meta.id}.strling.tsv
    """
}
