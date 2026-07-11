process TRA_CONSENSUS {
    tag "${meta.id}"
    label 'process_single'
    container "${params.utils_container}"
    publishDir "${params.outdir}/${meta.id}", mode: 'copy', pattern: "*.sv_merged.vcf.gz*"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.id}.sv_merged.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.sv_merged.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                                        emit: versions

    script:
    """
    # Re-cluster translocations across callers and recompute SUPP/SUPP_VEC.
    # Jasmine leaves every cross-caller TRA at SUPP=1 (it doesn't co-cluster
    # interchromosomal breakends), so the downstream SUPP>=2 gate erased all real
    # translocations. This rebuilds genuine multi-caller support for TRA.
    export PATH=${projectDir}/bin:\$PATH

    tra_consensus.py ${vcf} --out ${meta.id}.tra.vcf --window ${params.tra_window}

    # Jasmine/this step do not coordinate-sort; tabix needs contiguous chrom blocks.
    grep '^#' ${meta.id}.tra.vcf > ${meta.id}.sv_merged.sorted.vcf
    grep -v '^#' ${meta.id}.tra.vcf | sort -T . -k1,1 -k2,2n >> ${meta.id}.sv_merged.sorted.vcf
    bgzip ${meta.id}.sv_merged.sorted.vcf
    mv ${meta.id}.sv_merged.sorted.vcf.gz ${meta.id}.sv_merged.vcf.gz
    tabix -p vcf ${meta.id}.sv_merged.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tra_consensus: 1.0
    END_VERSIONS
    """
}
