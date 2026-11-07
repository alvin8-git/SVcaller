// Sex from read depth (bin/infer_sex.py, same rule as OmniGen vcf_qc.py). Feeds ExpansionHunter --sex
// (EH defaults to female, which genotypes chrX loci as diploid in a male) and the str_summary.tsv header.
// Container: svcaller/utils (python3 + samtools).
process INFER_SEX {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/${meta.id}", mode: 'copy', pattern: "*.sex.json"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.sex.json"), path("${meta.id}.sex.txt"), emit: sex

    script:
    """
    infer_sex.py ${bam} ${meta.id} ${meta.id}.sex.json ${meta.id}.sex.txt
    """
}
