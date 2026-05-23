include { MANTA_CALL             } from '../modules/manta/call'
include { MANTA_RESIDUAL_REGIONS } from '../modules/manta/residual_regions'
include { DELLY_CALL_SVTYPE      } from '../modules/delly/call'
include { DELLY_MERGE            } from '../modules/delly/merge'
include { GRIDSS_CALL            } from '../modules/gridss/call'
include { GRIDSS_SETUP           } from '../modules/gridss/setup'
include { GRIDSS_STUB            } from '../modules/gridss/stub'
include { SAMTOOLS_SUBSET        } from '../modules/samtools/subset'
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

    // DELLY: fan out 5 SV types in parallel; collect groupTuple(size:5) → merge
    ch_delly_input = ch_filtered_bam.combine(Channel.from(['DEL', 'INS', 'INV', 'DUP', 'BND']))
    DELLY_CALL_SVTYPE(ch_delly_input, ch_fasta, ch_fai)
    DELLY_MERGE(DELLY_CALL_SVTYPE.out.vcf.groupTuple(size: 5))

    EXPANSIONHUNTER(ch_bam, ch_fasta, ch_fai, ch_eh_catalog)

    // Manta always runs first; its output feeds tiered GRIDSS when enabled
    MANTA_CALL(ch_filtered_bam, ch_fasta, ch_fai)

    if (!params.skip_gridss) {
        GRIDSS_SETUP(ch_fasta, ch_fai)

        if (params.tiered_gridss) {
            // Tiered mode: GRIDSS runs only on regions where Manta was not PASS.
            // Trades GRIDSS parallelism with Manta for a much smaller input BAM
            // (~5–20% of reads in most 30× WGS samples), cutting GRIDSS wall time
            // from ~4–6 h to ~30–60 min per sample.
            MANTA_RESIDUAL_REGIONS(MANTA_CALL.out.vcf, ch_fai)
            SAMTOOLS_SUBSET(
                ch_filtered_bam.join(MANTA_RESIDUAL_REGIONS.out.bed)
            )
            GRIDSS_CALL(
                SAMTOOLS_SUBSET.out.bam, ch_fasta, ch_fai,
                GRIDSS_SETUP.out.amb, GRIDSS_SETUP.out.ann,
                GRIDSS_SETUP.out.bwt, GRIDSS_SETUP.out.pac,
                GRIDSS_SETUP.out.sa
            )
        } else {
            // Standard mode: GRIDSS runs on full filtered BAM in parallel with Manta
            GRIDSS_CALL(
                ch_filtered_bam, ch_fasta, ch_fai,
                GRIDSS_SETUP.out.amb, GRIDSS_SETUP.out.ann,
                GRIDSS_SETUP.out.bwt, GRIDSS_SETUP.out.pac,
                GRIDSS_SETUP.out.sa
            )
        }
        ch_gridss_vcf = GRIDSS_CALL.out.vcf
    } else {
        GRIDSS_STUB(ch_filtered_bam)
        ch_gridss_vcf = GRIDSS_STUB.out.vcf
    }

    // Collect VCFs per sample and merge with JASMINE
    // min_support=2 (3 callers) or 1 (2 callers when skip_gridss=true)
    ch_to_merge = MANTA_CALL.out.vcf
        .join(DELLY_MERGE.out.vcf)
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
