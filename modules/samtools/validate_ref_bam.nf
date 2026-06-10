process VALIDATE_REF_BAM {
    tag "$meta.id"
    label 'process_single'

    container 'quay.io/biocontainers/samtools:1.23.1--ha83d96e_0'

    input:
    tuple val(meta), path(bam), path(bai)
    path fai

    output:
    tuple val(meta), path(bam), path(bai), emit: bam

    script:
    """
    bam_chroms=\$(samtools view -H ${bam} | awk '/^@SQ/{gsub("SN:","",$2); print \$2}' | sort)
    ref_chroms=\$(cut -f1 ${fai} | sort)
    missing=\$(comm -23 <(echo "\$ref_chroms") <(echo "\$bam_chroms"))
    if [ -n "\$missing" ]; then
        echo "ERROR: reference has chromosomes not present in BAM ${bam}:" >&2
        echo "\$missing" >&2
        echo "" >&2
        echo "Hint: Use hg38.canonical.fa (not hg38.fa) for BAM inputs." >&2
        echo "      FILTER_CHROMS strips alt contigs from the BAM; the reference must match." >&2
        exit 1
    fi
    echo "VALIDATE_REF_BAM: ${meta.id} OK (ref and BAM chromosomes consistent)"
    """
}
