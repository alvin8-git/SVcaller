include { MANTA_CALL       } from '../modules/manta/call'
include { DELLY_CALL       } from '../modules/delly/call'
include { GRIDSS_CALL      } from '../modules/gridss/call'
include { EXPANSIONHUNTER  } from '../modules/expansionhunter/call'
include { JASMINE_MERGE    } from '../modules/jasmine/merge'

workflow SV_CALLING {
    take:
    ch_bam        // [ meta, bam, bai ]
    ch_fasta      // path
    ch_fai        // path
    ch_eh_catalog // path

    main:
    // Run 3 structural callers in parallel
    MANTA_CALL(ch_bam, ch_fasta, ch_fai)
    DELLY_CALL(ch_bam, ch_fasta, ch_fai)
    GRIDSS_CALL(ch_bam, ch_fasta, ch_fai)
    EXPANSIONHUNTER(ch_bam, ch_fasta, ch_fai, ch_eh_catalog)

    // Collect 3 structural VCFs per sample and merge with JASMINE (min_support=2)
    ch_to_merge = MANTA_CALL.out.vcf
        .join(DELLY_CALL.out.vcf)
        .join(GRIDSS_CALL.out.vcf)
        .map { meta, manta_vcf, delly_vcf, gridss_vcf ->
            [meta, [manta_vcf, delly_vcf, gridss_vcf]]
        }

    // Inner join: sample must complete all 3 callers to reach merge step (fail-fast on caller error)
    JASMINE_MERGE(ch_to_merge, ch_fasta, ch_fai)

    emit:
    sv_vcf  = JASMINE_MERGE.out.vcf
    str_vcf = EXPANSIONHUNTER.out.vcf
}
