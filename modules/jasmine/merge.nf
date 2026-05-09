process JASMINE_MERGE {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(vcfs)   // list of 3 VCF.gz files [manta, delly, gridss]
    path fasta
    path fai

    output:
    tuple val(meta), path("${meta.id}.sv_merged.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.sv_merged.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                                        emit: versions

    script:
    def vcf_list_str = vcfs.collect { "\"${it.baseName}\"" }.join(' ')
    """
    # Decompress each input VCF explicitly
    bgzip -d -c ${vcfs[0]} > ${vcfs[0].baseName}
    bgzip -d -c ${vcfs[1]} > ${vcfs[1].baseName}
    bgzip -d -c ${vcfs[2]} > ${vcfs[2].baseName}

    # Build vcf_list.txt from known filenames (avoids ls glob ambiguity)
    printf '%s\\n' ${vcfs[0].baseName} ${vcfs[1].baseName} ${vcfs[2].baseName} > vcf_list.txt

    jasmine \\
        file_list=vcf_list.txt \\
        out_file=${meta.id}.sv_merged.vcf \\
        genome_file=${fasta} \\
        min_support=2 \\
        --dup_to_ins \\
        --normalize_type \\
        --ignore_strand

    bgzip ${meta.id}.sv_merged.vcf
    tabix -p vcf ${meta.id}.sv_merged.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jasmine: \$(jasmine --version 2>&1 | head -1 | awk '{print \$NF}')
    END_VERSIONS
    """
}
