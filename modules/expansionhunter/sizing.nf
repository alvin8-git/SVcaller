// Per-allele sizing class for ExpansionHunter calls (bin/eh_sizing.py).
// Reads <id>.str_profile.json (left untouched) and writes <id>.str_summary.tsv with the
// alleleN_sizing / alleleN_display / conf columns. Container: svcaller/utils (python3).
process EH_SIZING {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/${meta.id}", mode: 'copy', pattern: "*.str_summary.tsv"

    input:
    tuple val(meta), path(json), path(sex_json)

    output:
    tuple val(meta), path("${meta.id}.str_summary.tsv"), emit: tsv

    script:
    """
    eh_sizing.py ${json} ${meta.id}.str_summary.tsv --sex-json ${sex_json}
    """
}
