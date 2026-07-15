process SVABA_CALL {
    tag "${meta.id}"
    label 'process_high'
    container 'quay.io/biocontainers/svaba:1.2.0--h69ac913_1'

    input:
    tuple val(meta), path(bam), path(bai)
    path ref_fasta
    path ref_fai
    // CLASSIC bwa index (.amb/.ann/.bwt/.pac/.sa). SvABA calls bwa_idx_load_from_disk,
    // which loads these from disk next to ${ref_fasta}. They MUST be staged here or SvABA
    // dies with "[E::bwa_idx_load_from_disk] fail to locate the index files". This is a
    // different format from the bwa-mem2 alignment index and is not interchangeable.
    // Historically these were never declared, so Nextflow never symlinked them and SvABA
    // failed on every run — masked by a since-removed `2>&1 || true`.
    path bwa_index

    output:
    tuple val(meta), path("${meta.id}.svaba.vcf.gz"), emit: vcf

    script:
    """
    # Germline single-sample SvABA run.
    # No '|| true' here: if SvABA crashes, this process must fail. Swallowing the exit
    # code let a crashed caller fall through to the empty-VCF stub below and render
    # downstream as a legitimate "0 SVs found".
    # To run without SvABA on purpose, use --skip_svaba (SVABA_STUB).
    svaba run \\
        -t ${bam} \\
        -G ${ref_fasta} \\
        -a ${meta.id} \\
        -p ${task.cpus} \\
        --germline

    # SvABA can exit 0 without writing its SV VCF if it died mid-run; treat that as failure.
    if [ ! -f "${meta.id}.svaba.sv.vcf" ]; then
        echo "ERROR: svaba run exited 0 but produced no ${meta.id}.svaba.sv.vcf." >&2
        echo "Refusing to emit an empty VCF that would look like '0 SVs found'." >&2
        ls -la >&2
        exit 1
    fi

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

    # A header-only VCF (SvABA ran fine, found nothing) is a VALID empty result and is
    # kept as such. What is no longer tolerated is reaching this point via a crash.
    tabix -p vcf ${meta.id}.svaba.vcf.gz

    # Guard the published artifact: it must at least carry a #CHROM header line.
    if ! tabix -H ${meta.id}.svaba.vcf.gz 2>/dev/null | grep -q "^#CHROM"; then
        echo "ERROR: normalised ${meta.id}.svaba.vcf.gz has no #CHROM header." >&2
        echo "Refusing to publish a malformed/empty VCF." >&2
        exit 1
    fi
    rm -f ${meta.id}.svaba.sv.vcf ${meta.id}.svaba.indel.vcf \
          ${meta.id}.svaba.bps.txt.gz ${meta.id}.svaba.discordant.txt.gz \
          ${meta.id}.svaba.contigs.bam ${meta.id}.svaba.alignments.txt.gz \
          ${meta.id}.svaba.log 2>/dev/null || true
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
