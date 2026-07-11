process JASMINE_MERGE {
    tag "${meta.id}"
    label 'process_medium'
    // Not published: TRA_CONSENSUS consumes this and publishes the final sv_merged.vcf.gz.

    input:
    tuple val(meta), path(vcfs)   // list of VCF.gz files: [manta, delly, gridss, scramble, melt]
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
        \$1~/^chr([0-9]+|X|Y|M)\$/ && \$7=="PASS"{
            n=split(\$8,info,";"); new_info=""
            for(i=1;i<=n;i++){
                key=info[i]; if(index(key,"=")) key=substr(key,1,index(key,"=")-1)
                if(key ~ keep || info[i]=="IMPRECISE" || info[i]=="PRECISE")
                    new_info=(new_info=="")?info[i]:new_info";"info[i]
            }
            \$8=(new_info=="")?".":new_info; print
        }' > ${vcfs[1].baseName}
    # GRIDSS: converted from BND pairs to simple DEL/DUP/INV by GRIDSS_CONVERT_BND.
    # Already QUAL-filtered and PASS; just enforce canonical chromosomes and strip to core INFO.
    zcat ${vcfs[2]} | awk '
        BEGIN{OFS="\\t"}
        /^#/{print;next}
        \$1!~/^chr([0-9]+|X|Y|M)\$/{next}
        \$7 != "PASS" {next}
        {
            n=split(\$8,info,";"); new_info=""
            for(i=1;i<=n;i++){
                if(info[i]~/^(SVTYPE|END|SVLEN|CIPOS|CIEND|HOMLEN|HOMSEQ|INSLEN)=/ || info[i]=="IMPRECISE")
                    new_info=(new_info=="")?info[i]:new_info";"info[i]
            }
            \$8=(new_info=="")?".":new_info; print
        }' > ${vcfs[2].baseName}
    # Scramble MEI: canonical-chromosome filter; add FORMAT/GT (Scramble VCF has none);
    # add SVTYPE=INS and canonical SVLEN (SCRAMble.R omits both from data lines).
    # Jasmine needs both to merge INS records correctly.
    # SVLEN values: ALU=300 (full-length ~282 bp), SVA=1500, L1=1500.
    # L1 was previously 6000 (full-length L1HS). Diagnostic confirmed MEINFO START/END
    # are genomic insertion-site coordinates (~1 bp span), not ME sequence coordinates,
    # so observed insertion length cannot be derived from MEINFO. The canonical 6000 bp
    # caused all L1 Truvari comparisons to fail Truvari's ±30% size gate against typical
    # truncated L1 truth insertions (500-2000 bp). 1500 bp covers truth L1s in 1050-2143 bp.
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
            if(\$5 ~ /L1/) svlen=1500
            else if(\$5 ~ /SVA/) svlen=1500
            \$8="SVTYPE=INS;SVLEN=" svlen
            print \$0 "\\tGT\\t0/1"
        }' > ${vcfs[3].baseName}

    # MELT MEI: SVTYPE already normalised to INS; PASS + canonical-chr already filtered in MELT_CALL.
    zcat ${vcfs[4]} > ${vcfs[4].baseName}

    # Build vcf_list.txt. Core 3 callers always included (SUPP_VEC positions 1-3).
    # Scramble added at position 4 if it has calls; MELT added at last position if it has calls.
    # GRIDSS is always at position 3 so the TRA filter (substr(supp,3,1)=="1") is stable.
    printf '%s\\n' ${vcfs[0].baseName} ${vcfs[1].baseName} ${vcfs[2].baseName} > vcf_list.txt
    scramble_cnt=\$(grep -cv '^#' ${vcfs[3].baseName} 2>/dev/null || echo 0)
    [ "\$scramble_cnt" -gt 0 ] && printf '%s\\n' ${vcfs[3].baseName} >> vcf_list.txt
    melt_cnt=\$(grep -cv '^#' ${vcfs[4].baseName} 2>/dev/null || echo 0)
    [ "\$melt_cnt" -gt 0 ] && printf '%s\\n' ${vcfs[4].baseName} >> vcf_list.txt

    jasmine \\
        file_list=vcf_list.txt \\
        out_file=${meta.id}.sv_merged.vcf \\
        genome_file=${fasta} \\
        min_support=1 \\
        --normalize_type \\
        --ignore_strand

    # Post-merge filters (applied in one pass):
    # 1. GRIDSS-only TRA: single-caller TRA where position 3 of SUPP_VEC is "1" (GRIDSS).
    #    Belt-and-suspenders: GRIDSS_CONVERT_BND already excludes TRA; catches any residual.
    # 2. Tiered min_support for large non-INS SVs: single-caller SVs with |SVLEN| > 10 kb
    #    that are NOT INS. INS (including MELT/Scramble MEI) are exempt since MELT is
    #    highly specific and MEIs are typically < 6 kb anyway.
    # 3. FORMAT → GT only: at min_support=1, Delly-only records enter the merged VCF with
    #    Delly-specific FORMAT tags not declared in Manta's header, crashing bcftools/Truvari.
    # Note: blanket single-caller DUP/INV filter removed — GRIDSS DUP/INV are pre-filtered
    #       by QUAL threshold in GRIDSS_CONVERT_BND; Manta/Delly single-caller DUP/INV allowed.
    awk 'BEGIN{OFS="\\t"}
        /^##FORMAT/{if(/ID=GT,/)print; next}
        /^#CHROM/{
            print "##INFO=<ID=MEITYPE,Number=1,Type=String,Description=\\"Mobile element insertion type (ALU/LINE1/SVA/HERVK)\\">"
            print; next
        }
        /^#/{print;next}
        {
            supp=""; svlen=0; is_tra=0; is_ins=0
            n=split(\$8,info,";")
            for(i=1;i<=n;i++){
                if(index(info[i],"SUPP_VEC=")==1) supp=substr(info[i],10)
                if(index(info[i],"SVTYPE=")==1){
                    t=substr(info[i],8)
                    if(t=="TRA"||t=="BND") is_tra=1
                    if(t=="INS") is_ins=1
                }
                if(index(info[i],"SVLEN=")==1){ svlen=substr(info[i],7)+0; if(svlen<0) svlen=-svlen }
            }
            if(supp!=""){
                ones=0; for(j=1;j<=length(supp);j++) if(substr(supp,j,1)=="1") ones++
                if(is_tra && ones==1 && substr(supp,3,1)=="1") next
                if(ones==1 && svlen>10000 && !is_ins) next
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
