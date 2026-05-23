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
    # Call each SV type to VCF (bgzip/tabix available; bcftools not in this container)
    for SVTYPE in DEL INS INV DUP BND; do
        delly call \\
            -t \${SVTYPE} \\
            -g ${fasta} \\
            -o ${meta.id}.delly.\${SVTYPE}.vcf \\
            ${bam}
    done

    # Merge: header from first file, data from all types sorted by chrom+pos
    { grep "^#" ${meta.id}.delly.DEL.vcf
      for sv in DEL INS INV DUP BND; do
          grep -v "^#" ${meta.id}.delly.\${sv}.vcf
      done | sort -k1,1V -k2,2n
    } | bgzip -c > ${meta.id}.delly.sv.vcf.gz
    tabix -p vcf ${meta.id}.delly.sv.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$(delly --version 2>&1 | grep "DELLY" | head -1 | awk '{print \$2}')
        bgzip: \$(bgzip --version 2>&1 | head -1 | awk '{print \$NF}')
    END_VERSIONS
    """
}
