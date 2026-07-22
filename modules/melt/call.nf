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

    # Auto-detect MELT installation if not provided via --melt_refs.
    # A missing MELT install is a MISCONFIGURATION, not a result. Emitting an empty VCF
    # here made "MELT was never installed" look identical to "this genome has zero mobile
    # element insertions". If MELT is intentionally unavailable, run with --skip_melt,
    # which routes to MELT_STUB and records the stage as explicitly skipped.
    if [ -z "\${MELT_REFS:-}" ]; then
        melt_jar=\$(find /usr/share /opt/conda/share /opt -name "MELT.jar" 2>/dev/null | head -1 || true)
        if [ -z "\$melt_jar" ]; then
            echo "ERROR: MELT.jar not found on this system and --melt_refs was not set." >&2
            echo "MELT cannot run. Refusing to emit an empty VCF that would be read as" >&2
            echo "'no mobile element insertions found'." >&2
            echo "Fix: install MELT / set --melt_refs, or run with --skip_melt to skip MELT explicitly." >&2
            exit 1
        fi
        MELT_DIR=\$(dirname "\$melt_jar")
        MELT_REFS="\${MELT_DIR}/me_refs/Hg38"
        MELT_BEDS="\${MELT_DIR}/add_bed_files/1KGP_Hg38"
    else
        MELT_BEDS="\${MELT_REFS}/../add_bed_files/1KGP_Hg38"
        melt_jar=\$(find /usr/share /opt/conda/share /opt -name "MELT.jar" 2>/dev/null | head -1 || true)
        if [ -z "\$melt_jar" ]; then
            echo "ERROR: --melt_refs was set to '\${MELT_REFS}' but no MELT.jar was found." >&2
            echo "Fix: install MELT, or run with --skip_melt to skip MELT explicitly." >&2
            exit 1
        fi
        MELT_DIR=\$(dirname "\$melt_jar")
    fi

    # Check ME refs exist -- again a misconfiguration, not an empty result.
    n_refs=\$(ls "\${MELT_REFS}"/*_MELT.zip 2>/dev/null | wc -l || echo 0)
    if [ "\$n_refs" -eq 0 ]; then
        echo "ERROR: no *_MELT.zip mobile-element references found in '\${MELT_REFS}'." >&2
        echo "MELT cannot detect any ME type without them. Refusing to emit an empty VCF." >&2
        echo "Fix: point --melt_refs at a valid me_refs/Hg38 dir, or run with --skip_melt." >&2
        exit 1
    fi

    # Run MELT Single for each ME type (ALU, LINE1/LINE, SVA, HERVK).
    # A non-zero exit from MELT is a crash, not "zero insertions of this type" -- the
    # previous '|| true' let every ME type die and still produce a confident empty VCF.
    failed_types=""
    mkdir -p melt_tmp
    for zip_file in "\${MELT_REFS}"/*_MELT.zip; do
        me=\$(basename "\$zip_file" _MELT.zip)
        # -n takes gene annotation BED12 (not prior sites); Hg38.genes.bed is standard BED12
        gene_bed="\${MELT_DIR}/add_bed_files/Hg38/Hg38.genes.bed"
        n_arg=""
        [ -f "\$gene_bed" ] && n_arg="-n \$gene_bed"
        if ! java -Xmx${mem_gb}g -jar "\$melt_jar" Single \\
            -bamfile ${bam} \\
            -h       ${ref_fasta} \\
            -t       "\$zip_file" \\
            \${n_arg} \\
            -w       melt_tmp/\${me} \\
            -c       ${params.min_depth} \\
            -sr      ${params.melt_min_reads} \\
            2>&1; then
            failed_types="\${failed_types} \${me}"
        fi
    done

    if [ -n "\$failed_types" ]; then
        echo "ERROR: MELT Single failed for ME type(s):\${failed_types}" >&2
        echo "A crashed MELT run must not be published as 'no mobile element insertions'." >&2
        echo "Fix the MELT failure above, or run with --skip_melt to skip MELT explicitly." >&2
        exit 1
    fi

    # Merge per-type VCFs; filter to canonical chr + PASS; normalise SVTYPE → INS
    first=""
    for vcf in melt_tmp/*/*.final_comp.vcf; do
        [ -f "\$vcf" ] || continue
        first="\$vcf"; break
    done

    # Every MELT type exited 0, so each should have written a *.final_comp.vcf (with zero
    # rows if it found nothing). No VCF at all means MELT did not actually run.
    if [ -z "\$first" ]; then
        echo "ERROR: MELT exited 0 for all ME types but wrote no *.final_comp.vcf." >&2
        echo "Refusing to emit an empty VCF that would be read as 'no MEIs found'." >&2
        find melt_tmp -type f >&2 || true
        exit 1
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
    rm -rf melt_tmp
    """
}
