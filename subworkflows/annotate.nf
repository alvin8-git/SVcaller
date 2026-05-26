include { ANNOTSV           } from '../modules/annotsv/annotate'
include { GNOMAD_SV_FILTER  } from '../modules/annotsv/annotate'

workflow ANNOTATE {
    take:
    ch_sv_vcf      // [ meta, sv_vcf.gz ]
    ch_annotsv_db  // path to AnnotSV db directory

    main:
    ANNOTSV(ch_sv_vcf, ch_annotsv_db)
    GNOMAD_SV_FILTER(ANNOTSV.out.tsv, 0.01)

    emit:
    tsv           = GNOMAD_SV_FILTER.out.tsv
    annotated_tsv = ANNOTSV.out.tsv
}
