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

    ch_fasta     = Channel.fromPath(params.ref_fasta, checkIfExists: true)
    ch_fai       = Channel.fromPath("${params.ref_fasta}.fai", checkIfExists: true)
    ch_dict      = Channel.fromPath("${params.ref_fasta}".replaceAll(/\.fa(sta)?$/, ".dict"),
                                    checkIfExists: false)
    ch_bwt_index = Channel.fromPath("${params.ref_fasta}.0123", checkIfExists: false)
                          .map { file(it.parent) }
    ch_pon       = params.pon
                    ? Channel.fromPath(params.pon, checkIfExists: true)
                    : Channel.value(file("NO_PON"))
    ch_intervals = params.intervals
                    ? Channel.fromPath(params.intervals, checkIfExists: true)
                    : Channel.value(file("NO_INTERVALS"))
    ch_annotsv   = params.annotsv_db
                    ? Channel.fromPath(params.annotsv_db, checkIfExists: true)
                    : Channel.value(file("NO_ANNOTSV"))
    ch_cytobands = Channel.fromPath("${projectDir}/assets/GRCh38_cytobands.txt",
                                     checkIfExists: false)
    ch_catalog   = Channel.fromPath(params.eh_catalog, checkIfExists: true)

    SVCALLER(
        ch_input, ch_fasta, ch_fai, ch_bwt_index,
        ch_dict, ch_pon, ch_intervals, ch_annotsv,
        ch_cytobands, ch_catalog,
    )
}
