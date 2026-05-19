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

    # Filter reads to ref chroms, reorder @SQ lines to match FAI order, re-sort reads
    # to match new header — required when input BAM uses different chrom ordering than ref
    samtools view -h -@ ${task.cpus} ${bam} \$REF_CHROMS | \\
        awk 'NR==FNR { order[\$1]=NR; total=NR; next }
             /^@HD/ { hd=\$0; next }
             /^@SQ/ {
                 for (i=1; i<=NF; i++) {
                     if (\$i ~ /^SN:/) {
                         c=\$i; sub(/^SN:/, "", c)
                         if (c in order) sq[order[c]]=\$0
                     }
                 }
                 next
             }
             /^@/ { other[++on]=\$0; next }
             {
                 if (!header_done) {
                     header_done=1
                     if (hd) print hd
                     for (i=1; i<=total; i++) if (i in sq) print sq[i]
                     for (i=1; i<=on; i++) print other[i]
                 }
                 print
             }' ${fai} - | \\
        samtools view -@ ${task.cpus} -b | \\
        samtools sort -@ ${task.cpus} -o ${meta.id}.filtered.bam
    samtools index ${meta.id}.filtered.bam
    """
}
