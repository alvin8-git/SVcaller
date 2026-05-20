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
    # Reads restricted to canonical chromosomes (no alt/decoy noise for CNV/SV calling).
    # ALL @SQ lines kept in FAI order so the sequence dictionary matches the reference
    # exactly — GRIDSS/Picard require dictionary size parity with the reference.
    # Reads are re-sorted to match the new (FAI-ordered) header, fixing order mismatches
    # between BAMs aligned to different hg38 builds (e.g. numeric vs alphabetical order).
    CANONICAL="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM"

    # Get canonical chroms in FAI order for samtools view region filter
    CANONICAL_ORDERED=\$(awk 'BEGIN{n=split("chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM",c); for(i=1;i<=n;i++) can[c[i]]=1} \$1 in can {print \$1}' ${fai} | tr '\\n' ' ')

    samtools view -h -@ ${task.cpus} ${bam} \$CANONICAL_ORDERED | \\
        awk 'BEGIN{
                 OFS="\t"
                 n=split("chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM",c)
                 for(i=1;i<=n;i++) can[c[i]]=1
             }
             NR==FNR { order[\$1]=NR; total=NR; next }
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
                 if (\$7 != "=" && \$7 != "*" && !(\$7 in can)) { \$7="*"; \$8=0 }
                 if (\$7 == "*" && \$8+0 != 0) { \$8=0 }
                 print
             }' ${fai} - | \\
        samtools view -@ ${task.cpus} -b | \\
        samtools sort -@ ${task.cpus} -o ${meta.id}.filtered.bam
    samtools index ${meta.id}.filtered.bam
    """
}
