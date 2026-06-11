process GATK_GCNV_CALL {
    tag "${meta.id}"
    label 'process_medium'
    container 'broadinstitute/gatk:4.5.0.0'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai
    path fasta_dict
    path pon_hdf5      // Panel of Normals from pon_build.nf
    path intervals     // preprocessed intervals BED

    output:
    tuple val(meta), path("${meta.id}.gatk_cnv.seg"),     emit: seg
    tuple val(meta), path("${meta.id}.gatk_cnv.tsv"),     emit: tsv
    tuple val(meta), path("${meta.id}.gatk_cnv.vcf.gz"),  emit: vcf
    tuple val(meta), path("${meta.id}.gatk_cnv.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                                   emit: versions

    script:
    def heap = (task.memory.toGiga() * 0.85).intValue()
    def pon_arg = pon_hdf5.name != 'NO_PON' ? "--count-panel-of-normals ${pon_hdf5}" : ""
    """
    # Collect read counts
    gatk --java-options "-Xmx${heap}g" CollectReadCounts \\
        -I ${bam} \\
        -L ${intervals} \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        -R ${fasta} \\
        -O ${meta.id}.counts.hdf5

    # Denoise against PoN
    gatk --java-options "-Xmx${heap}g" DenoiseReadCounts \\
        -I ${meta.id}.counts.hdf5 \\
        ${pon_arg} \\
        --standardized-copy-ratios ${meta.id}.standardizedCR.tsv \\
        --denoised-copy-ratios ${meta.id}.denoisedCR.tsv

    # Model segments
    gatk --java-options "-Xmx${heap}g" ModelSegments \\
        --denoised-copy-ratios ${meta.id}.denoisedCR.tsv \\
        --output . \\
        --output-prefix ${meta.id}

    # Call CNV segments
    gatk --java-options "-Xmx${heap}g" CallCopyRatioSegments \\
        --input ${meta.id}.cr.seg \\
        --output ${meta.id}.gatk_cnv.seg

    # Convert to simple TSV for cnv_consensus.py.
    # CallCopyRatioSegments columns: CONTIG START END NUM_POINTS_COPY_RATIO MEAN_LOG2_COPY_RATIO CALL
    # The categorical call (+/-/0) is column 6, NOT column 5 (which is MEAN_LOG2_COPY_RATIO).
    grep -v "^@" ${meta.id}.gatk_cnv.seg \\
        | awk 'NR==1{print "CONTIG\tSTART\tEND\tCALL_COPY_NUMBER\tQUALITY"; next}
               {cn=(\$6=="+")?3:(\$6=="-")?1:2; print \$1"\t"\$2"\t"\$3"\t"cn"\t50"}' \\
        > ${meta.id}.gatk_cnv.tsv

    # Produce compressed stub VCF for reporting
    printf "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" \\
        | bgzip > ${meta.id}.gatk_cnv.vcf.gz
    tabix -p vcf ${meta.id}.gatk_cnv.vcf.gz

    # Clean up intermediates
    rm -f ${meta.id}.counts.hdf5 ${meta.id}.standardizedCR.tsv ${meta.id}.denoisedCR.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep -oP '(?<=GATK v)[0-9.]+' | head -1)
    END_VERSIONS
    """
}
