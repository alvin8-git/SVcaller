process ANNOTSV {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/annotsv:3.4.6--py313hdfd78af_0'

    input:
    tuple val(meta), path(sv_vcf)
    path annotsv_db    // path to AnnotSV annotation directory

    output:
    tuple val(meta), path("${meta.id}.annotated.tsv"), emit: tsv
    path "versions.yml",                                emit: versions

    script:
    def skip_annotsv = annotsv_db.name == 'NO_ANNOTSV'
    if (skip_annotsv)
    """
    printf 'SV_chrom\tSV_start\tSV_end\tSV_type\tAnnotSV_ranking_score\tGene_name\tOMIM_morbid\tB_gain_AFmax\tB_loss_AFmax\n' \
        > ${meta.id}.annotated.tsv
    printf '"${task.process}":\n    annotsv: skipped (no --annotsv_db)\n' > versions.yml
    """
    else
    """
    AnnotSV \\
        -SVinputFile ${sv_vcf} \\
        -annotationsDir \$(dirname ${annotsv_db}) \\
        -genomeBuild GRCh38 \\
        -outputFile ${meta.id}.annotated \\
        -SVminSize 50 \\
        -tx ENSEMBL \\
        -annotationMode both

    [ -f "${meta.id}.annotated.tsv" ] || mv ${meta.id}.annotated.tsv.gz ${meta.id}.annotated.tsv 2>/dev/null || touch ${meta.id}.annotated.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotsv: \$(AnnotSV -help 2>&1 | grep "AnnotSV" | head -1 | awk '{print \$NF}')
    END_VERSIONS
    """
}

process GNOMAD_SV_FILTER {
    tag "${meta.id}"
    label 'process_single'
    container 'quay.io/biocontainers/annotsv:3.4.6--py313hdfd78af_0'
    publishDir "${params.outdir}/${meta.id}", mode: 'copy', pattern: "*.filtered.tsv"

    input:
    tuple val(meta), path(annotated_tsv)
    val af_threshold  // default 0.01

    output:
    tuple val(meta), path("${meta.id}.filtered.tsv"), emit: tsv

    script:
    """
    python3 -c "
import sys, csv
with open('${annotated_tsv}') as fh:
    reader = csv.DictReader(fh, delimiter='\\t')
    fieldnames = reader.fieldnames or []
    with open('${meta.id}.filtered.tsv', 'w') as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter='\\t')
        writer.writeheader()
        for row in reader:
            pop_af_str = row.get('B_gain_AFmax', row.get('B_loss_AFmax', '0')) or '0'
            try:
                pop_af = float(pop_af_str)
            except ValueError:
                pop_af = 0.0
            if pop_af < ${af_threshold}:
                writer.writerow(row)
"
    """
}
