process SV_PON_ANNOTATE {
    tag "${meta.id}"
    label 'process_single'
    container 'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path pon_bed   // recurrent SV sites (chrom, start, end)

    output:
    tuple val(meta), path("${meta.id}.sv_pon.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.sv_pon.vcf.gz.tbi"), emit: tbi

    script:
    """
    # Get IDs of SVs overlapping PON sites (reciprocal 50% overlap)
    bedtools intersect -a ${vcf} -b ${pon_bed} -u -f 0.5 -r -wa \\
        | grep -v "^#" | cut -f3 > pon_hit_ids.txt || true

    # Annotate matching SVs with SV_PON=1 in INFO
    zcat ${vcf} | awk -v hitfile="pon_hit_ids.txt" '
        BEGIN {
            while ((getline line < hitfile) > 0) hits[line] = 1
            close(hitfile)
        }
        /^##/ { print; next }
        /^#CHROM/ {
            print "##INFO=<ID=SV_PON,Number=0,Type=Flag,Description=\"SV overlaps recurrent site in GIAB multi-sample PON (present in >=2 of HG001/HG003-HG007)\">"
            print; next
        }
        {
            if (\$3 in hits) {
                sub(/INFO=/, "INFO=SV_PON;")
            }
            print
        }
    ' | bgzip > ${meta.id}.sv_pon.vcf.gz

    tabix -p vcf ${meta.id}.sv_pon.vcf.gz
    """
}
