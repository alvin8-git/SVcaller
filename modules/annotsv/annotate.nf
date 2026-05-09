process ANNOTSV {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/annotsv:3.4.2--pl5321hdfd78af_0'

    input:
    tuple val(meta), path(sv_vcf)
    path annotsv_db    // path to AnnotSV annotation directory

    output:
    tuple val(meta), path("${meta.id}.annotated.tsv"), emit: tsv
    path "versions.yml",                                emit: versions

    script:
    """
    AnnotSV \\
        -SVinputFile ${sv_vcf} \\
        -annotationsDir ${annotsv_db} \\
        -genome GRCh38 \\
        -outputFile ${meta.id}.annotated \\
        -SVminSize 50 \\
        -tx ENSEMBL \\
        -annotationMode both

    # Ensure output file is present (AnnotSV may add .tsv extension)
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
    container 'quay.io/biocontainers/annotsv:3.4.2--pl5321hdfd78af_0'

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
