include { SMN_CALLER } from '../modules/smn_caller/call'

workflow SMN_CALLING {
    take:
    ch_bam    // [ meta, bam, bai ]
    ch_fasta  // path
    ch_fai    // path

    main:
    SMN_CALLER(ch_bam, ch_fasta, ch_fai)

    emit:
    tsv  = SMN_CALLER.out.tsv
    json = SMN_CALLER.out.json
}
