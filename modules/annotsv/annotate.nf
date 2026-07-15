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
        -outputFile \$PWD/${meta.id}.annotated \\
        -SVminSize 50 \\
        -tx ENSEMBL \\
        -annotationMode both

    # AnnotSV may write to a date-stamped subdir even with an absolute -outputFile path
    if [ ! -s "${meta.id}.annotated.tsv" ]; then
        f=\$(find . -mindepth 2 -maxdepth 3 -name "${meta.id}.annotated.tsv" | head -1)
        if [ -n "\$f" ]; then
            mv "\$f" ./${meta.id}.annotated.tsv
        fi
    fi

    # Still nothing? Distinguish a LEGITIMATE empty result from a masked failure.
    # AnnotSV can exit 0 and emit no file when there is genuinely nothing to annotate
    # (no SV records >= -SVminSize 50 in the input). That is a real, valid empty result
    # and must be published as a header-only TSV -- never as a zero-byte file.
    # But if the input DID contain SV records and AnnotSV still produced nothing, the
    # annotation silently failed and this process must fail. Publishing an empty
    # placeholder here previously let a crashed stage render downstream as "no findings".
    if [ ! -s "${meta.id}.annotated.tsv" ]; then
        n_in=\$(zcat -f ${sv_vcf} | grep -vc '^#' || true)
        if [ "\${n_in:-0}" -eq 0 ]; then
            echo "NOTE: input ${sv_vcf} contains 0 SV records; emitting header-only annotation." >&2
            printf 'SV_chrom\tSV_start\tSV_end\tSV_type\tAnnotSV_ranking_score\tGene_name\tOMIM_morbid\tB_gain_AFmax\tB_loss_AFmax\n' \\
                > ${meta.id}.annotated.tsv
        else
            echo "ERROR: AnnotSV produced no output for '${meta.id}' but the input VCF" >&2
            echo "(${sv_vcf}) contains \${n_in} SV record(s). The annotation step failed silently." >&2
            echo "Refusing to publish an empty placeholder TSV. Work dir contents:" >&2
            ls -laR . >&2
            exit 1
        fi
    fi

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
