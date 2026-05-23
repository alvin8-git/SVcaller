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
    # Re-sort only when BAM chromosome order differs from FAI order (saves ~3h when already matching).
    CANONICAL_ORDERED=\$(awk 'BEGIN{n=split("chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM",c); for(i=1;i<=n;i++) can[c[i]]=1} \$1 in can {print \$1}' ${fai} | tr '\\n' ' ')

    # Detect whether BAM canonical chr order already matches FAI order
    BAM_CANON=\$(samtools view -H ${bam} | \\
        awk '/^@SQ/{for(i=1;i<=NF;i++) if(\$i ~ /^SN:/){sub(/SN:/,"",\$i); if(\$i ~ "^chr([0-9]+|X|Y|M)\$") print \$i}}' | \\
        tr '\\n' ',')
    FAI_CANON=\$(awk '\$1 ~ "^chr([0-9]+|X|Y|M)\$" {print \$1}' ${fai} | tr '\\n' ',')

    # Shared filter pipeline: view → reheader/fixup awk → compress
    filter_pipe() {
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
                 for (fi=12; fi<=NF; fi++) {
                     if (\$fi ~ /^SA:Z:/) {
                         n_sa=split(substr(\$fi,6),saparts,";")
                         new_sa=""
                         for (sj=1; sj<=n_sa; sj++) {
                             if (saparts[sj]=="") continue
                             split(saparts[sj],saf,",")
                             if (saf[1] in can) new_sa=new_sa saf[1]","saf[2]","saf[3]","saf[4]","saf[5]","saf[6]";"
                         }
                         if (new_sa!="") \$fi="SA:Z:"new_sa
                         else { for (fk=fi; fk<NF; fk++) \$fk=\$(fk+1); NF--; fi-- }
                         break
                     }
                 }
                 print
             }' ${fai} - | \\
        samtools view -@ ${task.cpus} -b
    }

    if [ "\$BAM_CANON" = "\$FAI_CANON" ]; then
        filter_pipe > ${meta.id}.filtered.bam
    else
        filter_pipe | samtools sort -@ ${task.cpus} -o ${meta.id}.filtered.bam
    fi
    samtools index ${meta.id}.filtered.bam
    """
}
