process JASMINE_MERGE {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/${meta.id}", mode: 'copy', pattern: "*.sv_merged.vcf.gz*"

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
    # Canonical-chromosome filter for all callers.
    # GRIDSS INFO fields are stripped to essentials (SVTYPE, MATEID, END, SVLEN, CIPOS, CIEND):
    # the full GRIDSS INFO contains multi-KB assembly sequences that cause Jasmine's merge to
    # produce binary-corrupt lines, crashing DuplicationsToInsertions post-processing.
    # DUPs are converted to INS for Manta pre-merge; --dup_to_ins removed from Jasmine to
    # avoid the post-merge DuplicationsToInsertions crash on complex GRIDSS-merged records.
    zcat ${vcfs[0]} | awk '
        BEGIN{OFS="\\t"}
        /^#/{print;next}
        \$1~/^chr([0-9]+|X|Y|M)\$/{
            if(\$5=="<DUP>"){\$5="<INS>"}
            print
        }' > ${vcfs[0].baseName}
    zcat ${vcfs[1]} | awk '/^#/{print;next} \$1~/^chr([0-9]+|X|Y|M)\$/{print}' > ${vcfs[1].baseName}
    zcat ${vcfs[2]} | awk '
        BEGIN{OFS="\\t"}
        /^#/{print;next}
        \$1!~/^chr([0-9]+|X|Y|M)\$/{next}
        \$7 != "PASS" {next}
        {
            n=split(\$8,info,";"); new_info=""
            for(i=1;i<=n;i++){
                if(info[i]~/^(SVTYPE|MATEID|END|SVLEN|CIPOS|CIEND|HOMLEN|HOMSEQ|INSLEN|IMPRECISE|EVENT)=/ || info[i]=="IMPRECISE")
                    new_info=(new_info=="")?info[i]:new_info";"info[i]
            }
            \$8=(new_info=="")?".":new_info; print
        }' > ${vcfs[2].baseName}

    # Build vcf_list.txt from known filenames (avoids ls glob ambiguity)
    printf '%s\\n' ${vcfs[0].baseName} ${vcfs[1].baseName} ${vcfs[2].baseName} > vcf_list.txt

    jasmine \\
        file_list=vcf_list.txt \\
        out_file=${meta.id}.sv_merged.vcf \\
        genome_file=${fasta} \\
        min_support=${params.skip_gridss ? 1 : 2} \\
        --normalize_type \\
        --ignore_strand

    # Jasmine output is not coordinate-sorted; tabix requires contiguous chromosome blocks
    grep '^#' ${meta.id}.sv_merged.vcf > ${meta.id}.sv_merged.sorted.vcf
    grep -v '^#' ${meta.id}.sv_merged.vcf | sort -k1,1 -k2,2n >> ${meta.id}.sv_merged.sorted.vcf
    bgzip ${meta.id}.sv_merged.sorted.vcf
    mv ${meta.id}.sv_merged.sorted.vcf.gz ${meta.id}.sv_merged.vcf.gz
    tabix -p vcf ${meta.id}.sv_merged.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jasmine: \$(jasmine --version 2>&1 | head -1 | awk '{print \$NF}')
    END_VERSIONS
    """
}
