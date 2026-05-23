include { MANTA_CALL             } from '../modules/manta/call'
include { DELLY_CALL             } from '../modules/delly/call'
include { GRIDSS_CALL            } from '../modules/gridss/call'
include { GRIDSS_SETUP           } from '../modules/gridss/setup'
include { GRIDSS_STUB            } from '../modules/gridss/stub'
include { EXPANSIONHUNTER        } from '../modules/expansionhunter/call'
include { JASMINE_MERGE          } from '../modules/jasmine/merge'
include { SAMTOOLS_FILTER_CHROMS } from '../modules/samtools/filter_chroms'

workflow SV_CALLING {
    take:
    ch_bam        // [ meta, bam, bai ]
    ch_fasta      // path
    ch_fai        // path
    ch_eh_catalog // path

    main:
    // Filter BAM to reference chromosomes (Manta, DELLY, GRIDSS all require BAM/ref parity)
    SAMTOOLS_FILTER_CHROMS(ch_bam, ch_fai)
    ch_filtered_bam = SAMTOOLS_FILTER_CHROMS.out.bam

    // Run Manta and DELLY in parallel on filtered BAM
    MANTA_CALL(ch_filtered_bam, ch_fasta, ch_fai)
    DELLY_CALL(ch_filtered_bam, ch_fasta, ch_fai)
    EXPANSIONHUNTER(ch_bam, ch_fasta, ch_fai, ch_eh_catalog)

    // Pre-build GRIDSS reference index once (storeDir caches across runs); skip per-sample rebuild
    if (!params.skip_gridss) {
        GRIDSS_SETUP(ch_fasta, ch_fai)
        GRIDSS_CALL(
            ch_filtered_bam, ch_fasta, ch_fai,
            GRIDSS_SETUP.out.amb, GRIDSS_SETUP.out.ann,
            GRIDSS_SETUP.out.bwt, GRIDSS_SETUP.out.pac,
            GRIDSS_SETUP.out.sa
        )
        ch_gridss_vcf = GRIDSS_CALL.out.vcf
    } else {
        GRIDSS_STUB(ch_filtered_bam)
        ch_gridss_vcf = GRIDSS_STUB.out.vcf
    }

    // Collect VCFs per sample and merge with JASMINE
    // min_support=2 (3 callers) or 1 (2 callers when skip_gridss=true)
    ch_to_merge = MANTA_CALL.out.vcf
        .join(DELLY_CALL.out.vcf)
        .join(ch_gridss_vcf)
        .map { meta, manta_vcf, delly_vcf, gridss_vcf ->
            [meta, [manta_vcf, delly_vcf, gridss_vcf]]
        }

    JASMINE_MERGE(ch_to_merge, ch_fasta, ch_fai)

    emit:
    sv_vcf  = JASMINE_MERGE.out.vcf
    sv_tbi  = JASMINE_MERGE.out.tbi
    str_vcf = EXPANSIONHUNTER.out.vcf
}
