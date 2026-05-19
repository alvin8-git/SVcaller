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
    # Restrict to canonical chromosomes only — alt/decoy contigs cause ambiguous
    # read mapping that inflates false CNV/SV calls and create header ordering mismatches
    CANONICAL="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM"

    # Get canonical chroms in FAI order (reference may use alphabetical or numeric ordering)
    REF_CHROMS=\$(awk 'BEGIN{n=split("chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM",c); for(i=1;i<=n;i++) can[c[i]]=1} \$1 in can {print \$1}' ${fai} | tr '\\n' ' ')

    # Filter reads to canonical chroms, output @SQ in FAI order, re-sort to match header
    samtools view -h -@ ${task.cpus} ${bam} \$REF_CHROMS | \\
        awk 'BEGIN{n=split("chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM",c); for(i=1;i<=n;i++) can[c[i]]=1}
             NR==FNR { if (\$1 in can) order[\$1]=++idx; next }
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
                     for (i=1; i<=idx; i++) if (i in sq) print sq[i]
                     for (i=1; i<=on; i++) print other[i]
                 }
                 print
             }' ${fai} - | \\
        samtools view -@ ${task.cpus} -b | \\
        samtools sort -@ ${task.cpus} -o ${meta.id}.filtered.bam
    samtools index ${meta.id}.filtered.bam
    """
}
