process DELLY_MERGE {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(vcfs)

    output:
    tuple val(meta), path("${meta.id}.delly.sv.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.delly.sv.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                                       emit: versions

    script:
    """
    # Delly outputs BCF (BGZF binary) even with .vcf extension — bcftools concat handles
    # both BCF and VCF uniformly; grep-based text merge silently corrupts BCF input.
    bcftools concat -a --allow-overlaps \\
        ${meta.id}.delly.DEL.vcf \\
        ${meta.id}.delly.INS.vcf \\
        ${meta.id}.delly.INV.vcf \\
        ${meta.id}.delly.DUP.vcf \\
        ${meta.id}.delly.BND.vcf | \\
    bcftools sort -O z -o ${meta.id}.delly.sv.vcf.gz
    tabix -p vcf ${meta.id}.delly.sv.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """
}
