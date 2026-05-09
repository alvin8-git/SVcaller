process DELLY_CALL {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path("${meta.id}.delly.sv.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.delly.sv.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                                       emit: versions

    script:
    """
    # Call all SV types
    for SVTYPE in DEL INS INV DUP TRA; do
        delly call \\
            -t \${SVTYPE} \\
            -g ${fasta} \\
            -o ${meta.id}.delly.\${SVTYPE}.bcf \\
            ${bam}
    done

    # Merge all SV types and convert to VCF
    bcftools concat -a ${meta.id}.delly.*.bcf \\
        | bcftools sort -O z -o ${meta.id}.delly.sv.vcf.gz
    bcftools index -t ${meta.id}.delly.sv.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$(delly --version 2>&1 | grep "DELLY" | head -1 | awk '{print \$2}')
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
