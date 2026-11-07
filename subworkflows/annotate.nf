include { ANNOTSV           } from '../modules/annotsv/annotate'
include { GNOMAD_SV_FILTER  } from '../modules/annotsv/annotate'
include { DEPTH_CHECK       } from '../modules/depth/depth_check'

workflow ANNOTATE {
    take:
    ch_sv_vcf      // [ meta, sv_vcf.gz ]
    ch_annotsv_db  // path to AnnotSV db directory
    ch_bam         // [ meta, bam, bai ] -- read-depth consistency check

    main:
    ANNOTSV(ch_sv_vcf, ch_annotsv_db)
    GNOMAD_SV_FILTER(ANNOTSV.out.tsv.join(ch_sv_vcf), 0.01)
    // Last step before <id>.filtered.tsv is written: append DEPTH_* columns (no rows dropped).
    DEPTH_CHECK(GNOMAD_SV_FILTER.out.tsv.join(ch_bam))

    emit:
    tsv           = DEPTH_CHECK.out.tsv
    annotated_tsv = ANNOTSV.out.tsv
}
