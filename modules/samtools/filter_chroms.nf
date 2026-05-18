process SAMTOOLS_FILTER_CHROMS {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)
    path fai

    output:
    tuple val(meta), path("${meta.id}.filtered.bam"), path("${meta.id}.filtered.bam.bai"), emit: bam

    script:
    """
    REF_CHROMS=\$(awk '{print \$1}' ${fai} | tr '\\n' ' ')

    # Filter reads to ref chroms AND strip non-ref @SQ lines from header in one pass
    samtools view -h -@ ${task.cpus} ${bam} \$REF_CHROMS | \\
        awk 'NR==FNR {chroms[\$1]=1; next}
             /^@SQ/ {
                 for (i=1; i<=NF; i++) {
                     if (\$i ~ /^SN:/) {
                         chrom = \$i; sub(/^SN:/, "", chrom)
                         if (!chroms[chrom]) next
                     }
                 }
                 print; next
             }
             {print}' ${fai} - | \\
        samtools view -@ ${task.cpus} -b -o ${meta.id}.filtered.bam
    samtools index ${meta.id}.filtered.bam
    """
}
