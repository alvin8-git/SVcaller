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

    # Step 2: call MEI — scramble.sh hardcodes --install-dir and --mei-refs internally
    scramble.sh \\
        --cluster-file clusters.txt \\
        --ref ${fasta} \\
        --out-name ${meta.id}.scramble \\
        --eval-meis || true

    if [ -f "${meta.id}.scramble.vcf" ] && grep -qv '^#' "${meta.id}.scramble.vcf" 2>/dev/null; then
        grep '^#' ${meta.id}.scramble.vcf > ${meta.id}.scramble.sorted.vcf
        grep -v '^#' ${meta.id}.scramble.vcf | sort -k1,1 -k2,2n >> ${meta.id}.scramble.sorted.vcf
        bgzip ${meta.id}.scramble.sorted.vcf
        mv ${meta.id}.scramble.sorted.vcf.gz ${meta.id}.scramble.vcf.gz
    else
        printf '##fileformat=VCFv4.2\\n##FILTER=<ID=PASS,Description="All filters passed">\\n##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\\n##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length">\\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t${meta.id}\\n' \\
            | bgzip > ${meta.id}.scramble.vcf.gz
    fi
    tabix -p vcf ${meta.id}.scramble.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scramble: \$(cluster_identifier --version 2>&1 | grep -o '[0-9]*\\.[0-9]*\\.[0-9]*' | head -1 || echo "unknown")
    END_VERSIONS
    """
}
