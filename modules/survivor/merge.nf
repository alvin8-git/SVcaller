process SURVIVOR_MERGE {
    tag "sv_pon"
    label 'process_low'
    container 'quay.io/biocontainers/survivor:1.0.7--hd03093a_2'

    input:
    path vcf_list   // text file with one VCF path per line
    val  max_dist   // max breakpoint distance (bp) to merge
    val  min_callers // min number of samples a site must appear in
    val  output_name

    output:
    path "${output_name}.vcf", emit: vcf

    script:
    """
    SURVIVOR merge ${vcf_list} ${max_dist} ${min_callers} 0 0 0 50 ${output_name}.vcf
    """
}
