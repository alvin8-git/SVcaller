// v2: container now includes bowtie2 (required by MELT Single for discordant read extraction)
process MELT_CALL {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)
    path ref_fasta
    path ref_fai

    output:
    tuple val(meta), path("${meta.id}.melt.vcf.gz"), emit: vcf

    script:
    def mem_gb = Math.max(task.memory.toGiga().toInteger() - 2, 4)
    def refs_override = params.melt_refs ? "MELT_REFS=\"${params.melt_refs}\"" : ""
    """
    ${refs_override}
    bowtie2 --version 2>&1 | head -1 || true

    # Auto-detect MELT installation if not provided via --melt_refs
    if [ -z "\${MELT_REFS:-}" ]; then
        melt_jar=\$(find /usr/share /opt/conda/share /opt -name "MELT.jar" 2>/dev/null | head -1 || true)
        if [ -z "\$melt_jar" ]; then
            echo "WARNING: MELT.jar not found; emitting empty VCF" >&2
            printf "##fileformat=VCFv4.1\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t${meta.id}\\n" | gzip > ${meta.id}.melt.vcf.gz
            exit 0
        fi
        MELT_DIR=\$(dirname "\$melt_jar")
        MELT_REFS="\${MELT_DIR}/me_refs/Hg38"
        MELT_BEDS="\${MELT_DIR}/add_bed_files/1KGP_Hg38"
    else
        MELT_BEDS="\${MELT_REFS}/../add_bed_files/1KGP_Hg38"
        melt_jar=\$(find /usr/share /opt/conda/share /opt -name "MELT.jar" 2>/dev/null | head -1 || true)
    fi

    # Check ME refs exist
    n_refs=\$(ls "\${MELT_REFS}"/*_MELT.zip 2>/dev/null | wc -l || echo 0)
    if [ "\$n_refs" -eq 0 ]; then
        echo "WARNING: no *_MELT.zip found in \${MELT_REFS}; emitting empty VCF" >&2
        printf "##fileformat=VCFv4.1\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t${meta.id}\\n" | gzip > ${meta.id}.melt.vcf.gz
        exit 0
    fi

    # Run MELT Single for each ME type (ALU, LINE1/LINE, SVA, HERVK)
    mkdir -p melt_tmp
    for zip_file in "\${MELT_REFS}"/*_MELT.zip; do
        me=\$(basename "\$zip_file" _MELT.zip)
        # -n takes gene annotation BED12 (not prior sites); Hg38.genes.bed is standard BED12
        gene_bed="\${MELT_DIR}/add_bed_files/Hg38/Hg38.genes.bed"
        n_arg=""
        [ -f "\$gene_bed" ] && n_arg="-n \$gene_bed"
        java -Xmx${mem_gb}g -jar "\$melt_jar" Single \\
            -bamfile ${bam} \\
            -h       ${ref_fasta} \\
            -t       "\$zip_file" \\
            \${n_arg} \\
            -w       melt_tmp/\${me} \\
            -c       ${params.min_depth} \\
            -reads   ${params.melt_min_reads} \\
            2>&1 || true
    done

    # Merge per-type VCFs; filter to canonical chr + PASS; normalise SVTYPE → INS
    first=""
    for vcf in melt_tmp/*/*.final_comp.vcf; do
        [ -f "\$vcf" ] || continue
        first="\$vcf"; break
    done

    if [ -z "\$first" ]; then
        printf "##fileformat=VCFv4.1\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t${meta.id}\\n" | gzip > ${meta.id}.melt.vcf.gz
        exit 0
    fi

    (
        grep "^#" "\$first"
        for vcf in melt_tmp/*/*.final_comp.vcf; do
            [ -f "\$vcf" ] || continue
            grep -v "^#" "\$vcf"
        done
    ) | awk '
        BEGIN{OFS="\\t"}
        /^##fileformat/{
            print
            print "##INFO=<ID=MEITYPE,Number=1,Type=String,Description=\\"Mobile element insertion type (ALU/LINE1/SVA/HERVK)\\">"
            next
        }
        /^#/{print;next}
        \$1!~/^chr([0-9]+|X|Y|M)\$/{next}
        \$7!="PASS" && \$7!="."{next}
        {
            # Strip INFO to essentials; normalise SVTYPE → INS; store ME type as MEITYPE
            svtype="INS"; meitype=""; svlen=""; end_=""
            n=split(\$8,info,";")
            for(i=1;i<=n;i++){
                f=info[i]
                if(f~/^SVTYPE=/){
                    t=substr(f,8)
                    if(t~/^(ALU|LINE1|LINE|SVA|HERVK)\$/) { meitype=t; svtype="INS" }
                    else svtype=t
                }
                else if(f~/^SVLEN=/) svlen=f
                else if(f~/^END=/)   end_=f
            }
            new_info="SVTYPE=" svtype
            if(meitype!="") new_info=new_info ";MEITYPE=" meitype
            if(svlen!="")   new_info=new_info ";" svlen
            if(end_!="")    new_info=new_info ";" end_
            \$8=new_info
            # Add FORMAT/GT if absent
            if(\$9=="" || \$9=="."){
                \$9="GT"; \$10="0/1"
            }
            print
        }' | sort -k1,1 -k2,2n | gzip > ${meta.id}.melt.vcf.gz
    """
}
