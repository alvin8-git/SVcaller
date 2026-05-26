process JASMINE_MERGE {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/${meta.id}", mode: 'copy', pattern: "*.sv_merged.vcf.gz*"

    input:
    tuple val(meta), path(vcfs)   // list of VCF.gz files: [manta, delly, gridss, scramble]
    path fasta
    path fai

    output:
    tuple val(meta), path("${meta.id}.sv_merged.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.sv_merged.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                                        emit: versions

    script:
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
    zcat ${vcfs[1]} | awk '
        BEGIN{OFS="\\t"; keep="^(SVTYPE|END|SVLEN|CIPOS|CIEND|HOMLEN|HOMSEQ|INSSEQ)\$"}
        /^#/{print;next}
        \$1~/^chr([0-9]+|X|Y|M)\$/{
            n=split(\$8,info,";"); new_info=""
            for(i=1;i<=n;i++){
                key=info[i]; if(index(key,"=")) key=substr(key,1,index(key,"=")-1)
                if(key ~ keep || info[i]=="IMPRECISE" || info[i]=="PRECISE")
                    new_info=(new_info=="")?info[i]:new_info";"info[i]
            }
            \$8=(new_info=="")?".":new_info; print
        }' > ${vcfs[1].baseName}
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
    # Scramble MEI: canonical-chromosome filter; add FORMAT/GT (Scramble VCF has none);
    # add SVTYPE=INS and canonical SVLEN (SCRAMble.R declares these in header but omits
    # them from data lines). Jasmine needs both to merge INS records correctly.
    zcat ${vcfs[3]} | awk '
        BEGIN{OFS="\\t"}
        /^##FORMAT/{next}
        /^#CHROM/{
            print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\\"Genotype\\">"
            print \$0 "\\tFORMAT\\t${meta.id}"
            next
        }
        /^#/{print;next}
        \$1~/^chr([0-9]+|X|Y|M)\$/{
            if(\$7 != "PASS" && \$7 != ".") next
            if(\$6+0 < 70) next
            svlen=300
            if(\$5 ~ /L1/) svlen=6000
            else if(\$5 ~ /SVA/) svlen=1500
            \$8="SVTYPE=INS;SVLEN=" svlen
            print \$0 "\\tGT\\t0/1"
        }' > ${vcfs[3].baseName}

    # Build vcf_list.txt. If Scramble produced no calls (stub or empty sample),
    # use 3 files so Jasmine SUPP_VEC stays 3-char (keeps GRIDSS TRA filter correct).
    scramble_cnt=\$(grep -cv '^#' ${vcfs[3].baseName} || true)
    if [ "\$scramble_cnt" -gt 0 ]; then
        printf '%s\\n' ${vcfs[0].baseName} ${vcfs[1].baseName} ${vcfs[2].baseName} ${vcfs[3].baseName} > vcf_list.txt
    else
        printf '%s\\n' ${vcfs[0].baseName} ${vcfs[1].baseName} ${vcfs[2].baseName} > vcf_list.txt
    fi

    jasmine \\
        file_list=vcf_list.txt \\
        out_file=${meta.id}.sv_merged.vcf \\
        genome_file=${fasta} \\
        min_support=1 \\
        --normalize_type \\
        --ignore_strand

    # Post-merge filters (applied in one pass):
    # 1. GRIDSS-only TRA: single-caller TRA where position 3 of SUPP_VEC is "1" (GRIDSS).
    #    Works for both 3-caller ("001") and 4-caller ("0010") SUPP_VEC.
    #    These are 100K+ GRIDSS BND noise records with no support from other callers.
    # 2. Tiered min_support for large SVs: single-caller SVs with |SVLEN| > 10 kb.
    #    Precision for large single-caller calls is very low (0.19); removing them improves
    #    F1 without hurting recall much (Scramble MEI are typically < 6 kb).
    # 3. FORMAT → GT only: at min_support=1, Delly-only records enter the merged VCF with
    #    Delly-specific FORMAT tags not declared in Manta's header, crashing bcftools/Truvari.
    awk 'BEGIN{OFS="\\t"}
        /^##FORMAT/{if(/ID=GT,/)print; next}
        /^#/{print;next}
        {
            supp=""; svlen=0; is_tra=0
            n=split(\$8,info,";")
            for(i=1;i<=n;i++){
                if(index(info[i],"SUPP_VEC=")==1) supp=substr(info[i],10)
                if(index(info[i],"SVTYPE=")==1){
                    t=substr(info[i],8)
                    if(t=="TRA"||t=="BND") is_tra=1
                }
                if(index(info[i],"SVLEN=")==1){ svlen=substr(info[i],7)+0; if(svlen<0) svlen=-svlen }
            }
            if(supp!=""){
                ones=0; for(j=1;j<=length(supp);j++) if(substr(supp,j,1)=="1") ones++
                if(is_tra && ones==1 && substr(supp,3,1)=="1") next
                if(ones==1 && svlen>10000) next
            }
            m=split(\$9,f,":"); gt_i=1
            for(i=1;i<=m;i++){if(f[i]=="GT"){gt_i=i;break}}
            split(\$10,s,":"); \$9="GT"; \$10=s[gt_i]
            print
        }' ${meta.id}.sv_merged.vcf > ${meta.id}.sv_merged_filt.vcf
    mv ${meta.id}.sv_merged_filt.vcf ${meta.id}.sv_merged.vcf

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
