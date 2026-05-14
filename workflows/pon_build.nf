include { GATK_PREPROCESS_INTERVALS } from '../modules/gatk/gcnv_pon'
include { GATK_COLLECT_COUNTS       } from '../modules/gatk/gcnv_pon'
include { GATK_CREATE_PON           } from '../modules/gatk/gcnv_pon'

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

    # Strip M5/UR from @SQ header lines so the dict matches the minimal dict
    # embedded in CollectReadCounts HDF5 files (derived from BAM headers)
    sed -i '/^@SQ/{ s/\tM5:[^\t]*//g; s/\tUR:[^\t]*//g; }' annotated_intervals.tsv

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
    // Bin intervals for read counting (required for GATK gCNV cohort mode)
    GATK_PREPROCESS_INTERVALS(ch_fasta, ch_fai, ch_dict, ch_intervals)

    // Annotate preprocessed intervals with GC content and mappability
    GATK_ANNOTATE_INTERVALS(ch_fasta, ch_fai, ch_dict, GATK_PREPROCESS_INTERVALS.out.preprocessed)

    // Collect read counts per sample using preprocessed (binned) intervals
    // .first() converts the single-item queue channel to a value channel so all 7 samples share it
    GATK_COLLECT_COUNTS(ch_bam, ch_fasta, ch_fai, ch_dict, GATK_PREPROCESS_INTERVALS.out.preprocessed.first())

    ch_all_hdf5 = GATK_COLLECT_COUNTS.out.hdf5.map { meta, h -> h }.collect()

    GATK_CREATE_PON(ch_all_hdf5)

    emit:
    pon = GATK_CREATE_PON.out.pon
}

workflow {
    if (!params.input)     error "ERROR: --input is required"
    if (!params.ref_fasta) error "ERROR: --ref_fasta is required"
    if (!params.intervals) error "ERROR: --intervals is required"

    ch_bam = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample]
            def bam  = file(row.bam, checkIfExists: true)
            def bai  = file("${row.bam}.bai", checkIfExists: true)
            [meta, bam, bai]
        }

    // Value channels so each sample can reuse the same ref files without exhausting the channel
    ch_fasta     = Channel.value(file(params.ref_fasta, checkIfExists: true))
    ch_fai       = Channel.value(file("${params.ref_fasta}.fai", checkIfExists: true))
    ch_dict      = Channel.value(file("${params.ref_fasta}".replaceAll(/\.fa(sta)?$/, ".dict")))
    ch_intervals = Channel.value(file(params.intervals, checkIfExists: true))

    PON_BUILD(ch_bam, ch_fasta, ch_fai, ch_dict, ch_intervals)
}
