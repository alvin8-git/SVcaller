// v3: MELT Single runs ONE mobile-element type per task. sv_calling.nf fans out
// the 4 types (ALU, HERVK, LINE1, SVA) in parallel then collects into MELT_MERGE.
// Previously one task looped the 4 types serially (~1.2 h). Fanning out drops the
// wall clock to the slowest single type (~30 min) at the cost of 4 concurrent
// java heaps, bounded by maxForks 4 and per-task memory (see the memory note below).
process MELT_CALL {
    tag "${meta.id}:${metype}"
    label 'process_medium'
    maxForks 4

    input:
    tuple val(meta), path(bam), path(bai), val(metype)
    path ref_fasta
    path ref_fai

    output:
    tuple val(meta), path("${meta.id}.${metype}.final_comp.vcf"), emit: vcf

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

    # Check THIS type's ME ref exists -- a misconfiguration, not an empty result.
    # (The old task globbed all *_MELT.zip; per-type we validate only the one we run.)
    zip_file="\${MELT_REFS}/${metype}_MELT.zip"
    if [ ! -f "\$zip_file" ]; then
        echo "ERROR: mobile-element reference '\$zip_file' not found." >&2
        echo "MELT cannot detect ${metype} without it. Refusing to emit an empty VCF." >&2
        echo "Fix: point --melt_refs at a valid me_refs/Hg38 dir, or run with --skip_melt." >&2
        exit 1
    fi

    # Run MELT Single for the single ME type ${metype}.
    # A non-zero exit from MELT is a crash, not "zero insertions of this type" -- no
    # '|| true' here; a failed type must fail the task so the merge is never fed a
    # confident empty VCF.
    mkdir -p melt_tmp/${metype}
    # -n takes gene annotation BED12 (not prior sites); Hg38.genes.bed is standard BED12
    gene_bed="\${MELT_DIR}/add_bed_files/Hg38/Hg38.genes.bed"
    n_arg=""
    [ -f "\$gene_bed" ] && n_arg="-n \$gene_bed"
    if ! java -Xmx${mem_gb}g -jar "\$melt_jar" Single \\
        -bamfile ${bam} \\
        -h       ${ref_fasta} \\
        -t       "\$zip_file" \\
        \${n_arg} \\
        -w       melt_tmp/${metype} \\
        -c       ${params.min_depth} \\
        -sr      ${params.melt_min_reads} \\
        2>&1; then
        echo "ERROR: MELT Single failed for ME type ${metype}." >&2
        echo "A crashed MELT run must not be published as 'no mobile element insertions'." >&2
        echo "Fix the MELT failure above, or run with --skip_melt to skip MELT explicitly." >&2
        exit 1
    fi

    # MELT exited 0, so it should have written a *.final_comp.vcf (zero rows if it found
    # nothing). No VCF at all means MELT did not actually run for this type.
    comp=\$(ls melt_tmp/${metype}/*.final_comp.vcf 2>/dev/null | head -1 || true)
    if [ -z "\$comp" ]; then
        echo "ERROR: MELT exited 0 for ${metype} but wrote no *.final_comp.vcf." >&2
        echo "Refusing to emit an empty VCF that would be read as 'no MEIs found'." >&2
        find melt_tmp/${metype} -type f >&2 || true
        exit 1
    fi

    # Hand the raw per-type VCF to MELT_MERGE untouched; all filtering/normalisation
    # (canonical chr + PASS, INFO strip, MEITYPE header) happens once in the merge.
    cp "\$comp" ${meta.id}.${metype}.final_comp.vcf
    rm -rf melt_tmp
    """
}
