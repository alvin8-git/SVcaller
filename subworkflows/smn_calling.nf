include { SMN_CALLER      } from '../modules/smn_caller/call'
include { SMN_DEPTH_FLOOR } from '../modules/smn_caller/depth_floor'

workflow SMN_CALLING {
    take:
    ch_bam       // [ meta, bam, bai ]
    ch_fasta     // path
    ch_fai       // path
    ch_coverage  // [ meta, <id>.mosdepth.summary.txt ]  (PREPROCESS.out.coverage; mean depth already measured)

    main:
    SMN_CALLER(ch_bam, ch_fasta, ch_fai)

    // Join on sample id only: the coverage meta carries needs_chr_filter, which may differ.
    ch_floor_in = SMN_CALLER.out.tsv
        .map { meta, tsv -> [meta.id, meta, tsv] }
        .join(ch_coverage.map { meta, summary -> [meta.id, summary] }, failOnMismatch: true)
        .map { id, meta, tsv, summary -> [meta, tsv, summary] }
    SMN_DEPTH_FLOOR(ch_floor_in, params.smn_min_depth)

    emit:
    tsv  = SMN_DEPTH_FLOOR.out.tsv
    json = SMN_CALLER.out.json
}
