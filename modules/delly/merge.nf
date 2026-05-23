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
    { grep "^#" ${meta.id}.delly.DEL.vcf
      for sv in DEL INS INV DUP BND; do
          grep -v "^#" ${meta.id}.delly.\${sv}.vcf
      done | sort -k1,1V -k2,2n
    } | bgzip -c > ${meta.id}.delly.sv.vcf.gz
    tabix -p vcf ${meta.id}.delly.sv.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bgzip: \$(bgzip --version 2>&1 | head -1 | awk '{print \$NF}')
        tabix: \$(tabix --version 2>&1 | head -1 | awk '{print \$NF}')
    END_VERSIONS
    """
}
