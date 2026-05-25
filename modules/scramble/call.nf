process SCRAMBLE_CALL {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path("${meta.id}.scramble.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.scramble.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                                       emit: versions

    script:
    """
    # Step 1: find split-read clusters
    cluster_identifier ${bam} > clusters.txt

    # Locate MEI consensus sequences bundled with the container
    MEI_REFS=\$(find /usr/local/share/scramble /opt/conda/share/scramble 2>/dev/null \
        -name 'MEI_consensus_seqs.fa' | head -1)
    [ -z "\$MEI_REFS" ] && MEI_REFS=/usr/local/share/scramble/MEI_consensus_seqs.fa

    # Step 2: call MEI (|| true: exits 1 when no clusters found)
    scramble \\
        --cluster-file clusters.txt \\
        --mei-refs "\$MEI_REFS" \\
        --ref ${fasta} \\
        --sample ${meta.id} \\
        --out-name ${meta.id}.scramble \\
        --eval-meis || true

    if [ -f "${meta.id}.scramble.MEI.vcf" ] && grep -qv '^#' "${meta.id}.scramble.MEI.vcf" 2>/dev/null; then
        grep '^#' ${meta.id}.scramble.MEI.vcf > ${meta.id}.scramble.sorted.vcf
        grep -v '^#' ${meta.id}.scramble.MEI.vcf | sort -k1,1 -k2,2n >> ${meta.id}.scramble.sorted.vcf
        bgzip ${meta.id}.scramble.sorted.vcf
        mv ${meta.id}.scramble.sorted.vcf.gz ${meta.id}.scramble.vcf.gz
    else
        printf '##fileformat=VCFv4.2\\n##FILTER=<ID=PASS,Description="All filters passed">\\n##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\\n##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length">\\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t${meta.id}\\n' \\
            | bgzip > ${meta.id}.scramble.vcf.gz
    fi
    tabix -p vcf ${meta.id}.scramble.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scramble: \$(cluster_identifier --version 2>&1 | grep -oP '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1 || echo "unknown")
    END_VERSIONS
    """
}
