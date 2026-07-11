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
    ch_pon       = params.pon
                    ? Channel.value(file(params.pon, checkIfExists: true))
                    : Channel.value(file("NO_PON"))
    ch_intervals = params.intervals
                    ? Channel.value(file(params.intervals, checkIfExists: true))
                    : Channel.value(file("NO_INTERVALS"))
    ch_annotsv   = params.annotsv_db
                    ? Channel.value(file(params.annotsv_db, checkIfExists: true))
                    : Channel.value(file("NO_ANNOTSV"))
    ch_cytobands = Channel.fromPath("${projectDir}/assets/GRCh38_cytobands.txt",
                                     checkIfExists: false)
    ch_catalog   = Channel.fromPath(params.eh_catalog, checkIfExists: true)
    ch_trait_regions = Channel.value(file("${projectDir}/assets/cnv_trait_regions.bed",
                                          checkIfExists: true))

    SVCALLER(
        ch_input, ch_fasta, ch_fai, ch_bwt_index,
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
