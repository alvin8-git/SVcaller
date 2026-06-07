process SV_PON_ANNOTATE {
    tag "${meta.id}"
    label 'process_single'
    container 'quay.io/biocontainers/annotsv:3.4.6--py313hdfd78af_0'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path pon_bed   // recurrent SV sites (chrom, start, end)

    output:
    tuple val(meta), path("${meta.id}.sv_pon.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.sv_pon.vcf.gz.tbi"), emit: tbi

    script:
    """
    # Normalize BED: swap inverted coords (start > end) before bedtools
    awk 'BEGIN{OFS="\\t"} NF>=3 && \$2+0>\$3+0 {t=\$2;\$2=\$3;\$3=t} {print}' ${pon_bed} > pon_norm.bed

    # Get IDs of SVs overlapping PON sites (reciprocal 50% overlap)
    bedtools intersect -a ${vcf} -b pon_norm.bed -u -f 0.5 -r -wa \\
        | grep -v "^#" | cut -f3 > pon_hit_ids.txt || true

    # Annotate matching SVs with SV_PON flag in INFO
    # Use sprintf("%c",34) to embed double-quote in awk without shell escaping issues
    zcat ${vcf} | awk -v hitfile="pon_hit_ids.txt" 'BEGIN {
            OFS = "\\t"
            while ((getline line < hitfile) > 0) hits[line] = 1
            close(hitfile)
            q = sprintf("%c", 34)
            infohdr = "##INFO=<ID=SV_PON,Number=0,Type=Flag,Description=" q "SV overlaps recurrent site in GIAB multi-sample PON (>=2 of HG001/HG003-HG007)" q ">"
        }
        /^##/ { print; next }
        /^#CHROM/ { print infohdr; print; next }
        { if (\$3 in hits) \$8 = "SV_PON;" \$8; print }
    ' | bgzip > ${meta.id}.sv_pon.vcf.gz

    tabix -p vcf ${meta.id}.sv_pon.vcf.gz
    """
}
