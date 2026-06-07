process SVABA_CALL {
    tag "${meta.id}"
    label 'process_high'
    container 'quay.io/biocontainers/svaba:1.2.0--h69ac913_1'

    input:
    tuple val(meta), path(bam), path(bai)
    path ref_fasta
    path ref_fai

    output:
    tuple val(meta), path("${meta.id}.svaba.vcf.gz"), emit: vcf

    script:
    """
    # Germline single-sample SvABA run
    svaba run \\
        -t ${bam} \\
        -G ${ref_fasta} \\
        -a ${meta.id} \\
        -p ${task.cpus} \\
        --germline \\
        2>&1 || true

    # Normalise: merge sv + indel VCFs, strip to canonical chrs, add SVLEN where absent
    (
        # Header from sv.vcf (most complete), with sample GT field
        grep "^##" ${meta.id}.svaba.sv.vcf 2>/dev/null || true
        printf '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n'
        printf '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">\n'
        printf '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n'
        grep "^#CHROM" ${meta.id}.svaba.sv.vcf 2>/dev/null || printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${meta.id}\n"

        for vcf in ${meta.id}.svaba.sv.vcf ${meta.id}.svaba.indel.vcf; do
            [ -f "\$vcf" ] || continue
            awk 'BEGIN{OFS="\t"}
                /^#/{next}
                {
                    # keep canonical chromosomes only
                    if (\$1 !~ /^chr([0-9XYM]|1[0-9]|2[0-2])\$/) next
                    # skip BND (translocations) — Jasmine handles these but adds complexity
                    if (\$5 ~ /^(\\.|\\/|[ACGTN]+\\[|[ACGTN]+\\]|\\[[^\\]]+\\]|\\][^\\[]+\\[)/) next
                    svtype = ""; svlen = 0; end = \$2
                    n = split(\$8, info, ";")
                    for (i=1;i<=n;i++) {
                        if (info[i] ~ /^SVTYPE=/) svtype = substr(info[i],8)
                        if (info[i] ~ /^SVLEN=/)  svlen  = substr(info[i],7)+0
                        if (info[i] ~ /^END=/)    end    = substr(info[i],5)+0
                    }
                    if (svtype == "") next
                    if (svlen == 0 && end > \$2) svlen = (svtype == "DEL") ? -(end-\$2) : (end-\$2)
                    # Rebuild INFO with essential tags
                    print \$1, \$2, \$3, \$4, \$5, \$6, "PASS",
                          "SVTYPE="svtype";SVLEN="svlen";END="end,
                          \$9, \$10
                }' "\$vcf"
        done
    ) | sort -k1,1 -k2,2n | bgzip > ${meta.id}.svaba.vcf.gz

    tabix -p vcf ${meta.id}.svaba.vcf.gz || true

    # Emit empty VCF stub if SvABA produced nothing
    if ! tabix -H ${meta.id}.svaba.vcf.gz 2>/dev/null | grep -q "^#CHROM"; then
        printf "##fileformat=VCFv4.1\n##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV type\">\n##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"SV length\">\n##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${meta.id}\n" \\
            | bgzip > ${meta.id}.svaba.vcf.gz
        tabix -p vcf ${meta.id}.svaba.vcf.gz || true
    fi
    """
}

process SVABA_STUB {
    tag "${meta.id}"
    label 'process_single'
    container 'quay.io/biocontainers/htslib:1.21--h566b1c6_0'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.svaba.vcf.gz"), emit: vcf

    script:
    """
    printf "##fileformat=VCFv4.1\n##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV type\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${meta.id}\n" \\
        | bgzip > ${meta.id}.svaba.vcf.gz
    tabix -p vcf ${meta.id}.svaba.vcf.gz
    """
}
