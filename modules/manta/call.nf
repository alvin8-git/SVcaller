process MANTA_CALL {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path("${meta.id}.manta.sv.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.manta.sv.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                                       emit: versions

    script:
    """
    configManta.py \\
        --bam ${bam} \\
        --referenceFasta ${fasta} \\
        --runDir manta_run

    python manta_run/runWorkflow.py \\
        -j ${task.cpus} \\
        -g ${task.memory.toGiga()}

    cp manta_run/results/variants/diploidSV.vcf.gz     ${meta.id}.manta.sv.vcf.gz
    cp manta_run/results/variants/diploidSV.vcf.gz.tbi ${meta.id}.manta.sv.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        manta: \$(configManta.py --version 2>&1 | head -1)
    END_VERSIONS
    """
}
