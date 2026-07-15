#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SVCALLER } from './workflows/svcaller'

def validate_params() {
    if (!params.input)     error "ERROR: --input is required"
    if (!params.ref_fasta) error "ERROR: --ref_fasta is required"
}

workflow {
    validate_params()

    // Parse samplesheet → channel of [meta, fq1|null, fq2|null, bam|null]
    ch_input = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample]
            def fq1  = row.fastq_1 ? file(row.fastq_1, checkIfExists: true) : null
            def fq2  = row.fastq_2 ? file(row.fastq_2, checkIfExists: true) : null
            def bam  = row.bam     ? file(row.bam,     checkIfExists: true) : null
            [meta, fq1, fq2, bam]
        }

    // Value channels so all subworkflows (PREPROCESS, SV_CALLING, CNV_CALLING, SMN_CALLING)
    // can each receive the same reference files without queue-channel exhaustion
    ch_fasta     = Channel.value(file(params.ref_fasta, checkIfExists: true))
    ch_fai       = Channel.value(file("${params.ref_fasta}.fai", checkIfExists: true))
    ch_dict      = Channel.value(file("${params.ref_fasta}".replaceAll(/\.fa(sta)?$/, ".dict"),
                                      checkIfExists: false))
    ch_bwt_index = Channel.value(file(params.ref_fasta).parent)

    // CLASSIC bwa index for SvABA. SvABA calls bwa_idx_load_from_disk internally, which
    // needs hg38.canonical.fa.{amb,ann,bwt,pac,sa} co-located with ref_fasta in the task
    // work dir. This is a DIFFERENT format from the bwa-mem2 index (.0123/.bwt.2bit.64)
    // used for alignment, so ch_bwt_index does NOT satisfy SvABA. Historically these files
    // were never staged and `2>&1 || true` masked the crash, so SvABA silently produced
    // nothing. Stage them explicitly (see modules/svaba/call.nf) — only required when SvABA
    // actually runs, so operators using --skip_svaba need not build the classic index.
    def bwa_prefix = params.bwa_index ?: params.ref_fasta
    ch_bwa_index = Channel.value(file("NO_BWA_INDEX"))
    if (!params.skip_svaba) {
        def bwa_exts    = ['amb', 'ann', 'bwt', 'pac', 'sa']
        def bwa_missing = bwa_exts.findAll { !file("${bwa_prefix}.${it}").exists() }
        if (bwa_missing) {
            error """ERROR: classic BWA index not found for ${bwa_prefix} """ +
                  """(missing: ${bwa_missing.collect { "${bwa_prefix}.${it}" }.join(', ')}).
SvABA needs the CLASSIC bwa index (.amb/.ann/.bwt/.pac/.sa) alongside the reference; the
bwa-mem2 alignment index (.0123/.bwt.2bit.64) is a different format and will NOT work.
Build it with:   bwa index ${bwa_prefix}
or run without SvABA:   --skip_svaba"""
        }
        ch_bwa_index = Channel.value(bwa_exts.collect {
            file("${bwa_prefix}.${it}", checkIfExists: true)
        })
    }

    ch_pon       = params.pon
                    ? Channel.value(file(params.pon, checkIfExists: true))
                    : Channel.value(file("NO_PON"))
    ch_intervals = params.intervals
                    ? Channel.value(file(params.intervals, checkIfExists: true))
                    : Channel.value(file("NO_INTERVALS"))
    ch_annotsv   = params.annotsv_db
                    ? Channel.value(file(params.annotsv_db, checkIfExists: true))
                    : Channel.value(file("NO_ANNOTSV"))
    // Value channels — shared reference files reused across every sample. A queue
    // channel (fromPath) holds one item and is consumed on the FIRST sample, so a
    // multi-sample process (e.g. CIRCOS_PLOT, EXPANSIONHUNTER) only runs once and the
    // rest of the samples silently vanish downstream. Must be Channel.value.
    ch_cytobands = Channel.value(file("${projectDir}/assets/GRCh38_cytobands.txt"))
    ch_catalog   = Channel.value(file(params.eh_catalog, checkIfExists: true))
    ch_trait_regions = Channel.value(file("${projectDir}/assets/cnv_trait_regions.bed",
                                          checkIfExists: true))

    SVCALLER(
        ch_input, ch_fasta, ch_fai, ch_bwt_index, ch_bwa_index,
        ch_dict, ch_pon, ch_intervals, ch_annotsv,
        ch_cytobands, ch_catalog, ch_trait_regions,
    )
}

// Event handlers must live in the main script, not nextflow.config — in the
// config `workflow` is a ConfigObject and `onComplete` is rejected at parse time.
workflow.onComplete {
    if (workflow.success && params.auto_cleanup) {
        def workDir = new File(workflow.workDir.toString())
        if (workDir.exists()) {
            log.info "auto_cleanup: removing work dir ${workDir}"
            workDir.deleteDir()
        }
    }
    if (workflow.success) {
        log.info "Pipeline complete. Results: ${params.outdir}"
        if (!params.auto_cleanup) {
            log.info "Tip: run 'bash bin/nf-cleanup.sh <sampleId>' to remove intermediates."
        }
    }
}
