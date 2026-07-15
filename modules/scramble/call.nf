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

    # Step 2: call MEI — scramble.sh hardcodes --install-dir and --mei-refs internally.
    # SCRAMble.R does setwd(INSTALL.DIR) before reading the cluster file, so absolute path required.
    # No '|| true': a SCRAMble crash previously fell through to the empty-VCF branch
    # below and was indistinguishable downstream from a real "no MEIs detected" result.
    # To run without SCRAMble on purpose, use --skip_scramble (SCRAMBLE_STUB).
    scramble.sh \\
        --cluster-file \$PWD/clusters.txt \\
        --ref \$PWD/${fasta} \\
        --out-name \$PWD/${meta.id}.scramble \\
        --eval-meis

    # SCRAMble exited 0. A VCF with only a header (or no MEI rows) is a VALID empty
    # result -- zero mobile-element insertions is a real finding -- so the else branch
    # below still emits a well-formed header-only VCF (never a zero-byte file).
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
