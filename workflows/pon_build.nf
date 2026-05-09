include { GATK_COLLECT_COUNTS } from '../modules/gatk/gcnv_pon'
include { GATK_CREATE_PON     } from '../modules/gatk/gcnv_pon'

/*
 * One-time workflow: build GATK gCNV Panel of Normals from GIAB samples HG001-HG007.
 * Run BEFORE the main svcaller pipeline.
 *
 * Usage:
 *   nextflow run workflows/pon_build.nf \
 *     --input giab_bam_samplesheet.csv \
 *     --ref_fasta /data/alvin/ref/GRCh38/GRCh38.fasta \
 *     --outdir /data/alvin/SVcaller/pon \
 *     -profile docker
 */

process GATK_ANNOTATE_INTERVALS {
    label 'process_single'
    container 'broadinstitute/gatk:4.5.0.0'

    input:
    path fasta
    path fai
    path dict
    path intervals

    output:
    path "annotated_intervals.tsv", emit: annotated
    path "versions.yml",             emit: versions

    script:
    def heap = (task.memory.toGiga() * 0.85).intValue()
    """
    gatk --java-options "-Xmx${heap}g" AnnotateIntervals \\
        -R ${fasta} \\
        -L ${intervals} \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        -O annotated_intervals.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep -oP '(?<=GATK v)[0-9.]+' | head -1)
    END_VERSIONS
    """
}

workflow PON_BUILD {
    take:
    ch_bam        // [ meta, bam, bai ]
    ch_fasta
    ch_fai
    ch_dict
    ch_intervals

    main:
    // Annotate intervals once (compute GC content, mappability)
    GATK_ANNOTATE_INTERVALS(ch_fasta, ch_fai, ch_dict, ch_intervals)

    GATK_COLLECT_COUNTS(ch_bam, ch_fasta, ch_fai, ch_dict, ch_intervals)

    ch_all_hdf5 = GATK_COLLECT_COUNTS.out.hdf5.map { meta, h -> h }.collect()

    GATK_CREATE_PON(ch_all_hdf5, GATK_ANNOTATE_INTERVALS.out.annotated)

    emit:
    pon = GATK_CREATE_PON.out.pon
}
