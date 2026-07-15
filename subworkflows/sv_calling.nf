include { MANTA_CALL             } from '../modules/manta/call'
include { MANTA_RESIDUAL_REGIONS } from '../modules/manta/residual_regions'
include { DELLY_CALL_SVTYPE      } from '../modules/delly/call'
include { DELLY_MERGE            } from '../modules/delly/merge'
include { GRIDSS_CALL            } from '../modules/gridss/call'
include { GRIDSS_CONVERT_BND     } from '../modules/gridss/convert_bnd'
include { GRIDSS_SETUP           } from '../modules/gridss/setup'
include { GRIDSS_STUB            } from '../modules/gridss/stub'
include { SAMTOOLS_SUBSET        } from '../modules/samtools/subset'
include { EXPANSIONHUNTER        } from '../modules/expansionhunter/call'
include { JASMINE_MERGE          } from '../modules/jasmine/merge'
include { TRA_CONSENSUS          } from '../modules/jasmine/tra_consensus'
include { SCRAMBLE_CALL          } from '../modules/scramble/call'
include { SCRAMBLE_STUB          } from '../modules/scramble/stub'
include { MELT_CALL              } from '../modules/melt/call'
include { MELT_STUB              } from '../modules/melt/stub'
include { SVABA_CALL             } from '../modules/svaba/call'
include { SVABA_STUB             } from '../modules/svaba/call'
include { STRLING_CALL           } from '../modules/strling/call'
include { SAMTOOLS_FILTER_CHROMS } from '../modules/samtools/filter_chroms'
include { VALIDATE_REF_BAM      } from '../modules/samtools/validate_ref_bam'

workflow SV_CALLING {
    take:
    ch_bam        // [ meta, bam, bai ]
    ch_fasta      // path
    ch_fai        // path
    ch_eh_catalog // path
    ch_bwa_index  // classic bwa index files (.amb/.ann/.bwt/.pac/.sa) for SvABA

    main:
    // FILTER_CHROMS: skip for FASTQ-derived BAMs (aligned to hg38.canonical.fa — no alt contigs).
    ch_bam.branch {
        needs_filter: it[0].get('needs_chr_filter', true)
        canonical:    true
    }.set { ch_bam_branched }

    SAMTOOLS_FILTER_CHROMS(ch_bam_branched.needs_filter, ch_fai)

    ch_filtered_bam = SAMTOOLS_FILTER_CHROMS.out.bam
        .mix(ch_bam_branched.canonical)

    // Pre-flight: verify BAM and reference chromosomes are consistent before any caller runs.
    // Catches hg38.fa-vs-canonical-BAM mismatch at second 0 instead of hour 4 (Manta crash).
    VALIDATE_REF_BAM(ch_filtered_bam, ch_fai)
    ch_validated_bam = VALIDATE_REF_BAM.out.bam

    // DELLY: fan out 5 SV types in parallel; collect groupTuple(size:5) → merge
    ch_delly_input = ch_validated_bam.combine(Channel.from(['DEL', 'INS', 'INV', 'DUP', 'BND']))
    DELLY_CALL_SVTYPE(ch_delly_input, ch_fasta, ch_fai)
    DELLY_MERGE(DELLY_CALL_SVTYPE.out.vcf.groupTuple(size: 5))

    EXPANSIONHUNTER(ch_bam, ch_fasta, ch_fai, ch_eh_catalog)

    // P6: STRling genome-wide STR expansion detection (parallel with EH)
    if (!params.skip_strling) {
        STRLING_CALL(ch_validated_bam, ch_fasta, ch_fai)
        ch_strling_tsv = STRLING_CALL.out.tsv
    } else {
        ch_strling_tsv = Channel.empty()
    }

    // Manta always runs first; its output feeds tiered GRIDSS when enabled
    MANTA_CALL(ch_validated_bam, ch_fasta, ch_fai)

    if (!params.skip_gridss) {
        GRIDSS_SETUP(ch_fasta, ch_fai)

        if (params.tiered_gridss) {
            MANTA_RESIDUAL_REGIONS(MANTA_CALL.out.vcf, ch_fai)
            SAMTOOLS_SUBSET(
                ch_validated_bam.join(MANTA_RESIDUAL_REGIONS.out.bed)
            )
            GRIDSS_CALL(
                SAMTOOLS_SUBSET.out.bam, ch_fasta, ch_fai,
                GRIDSS_SETUP.out.amb, GRIDSS_SETUP.out.ann,
                GRIDSS_SETUP.out.bwt, GRIDSS_SETUP.out.pac,
                GRIDSS_SETUP.out.sa
            )
        } else {
            GRIDSS_CALL(
                ch_validated_bam, ch_fasta, ch_fai,
                GRIDSS_SETUP.out.amb, GRIDSS_SETUP.out.ann,
                GRIDSS_SETUP.out.bwt, GRIDSS_SETUP.out.pac,
                GRIDSS_SETUP.out.sa
            )
        }
        // Convert BND pairs → simple DEL/DUP/INV so Jasmine can merge them correctly.
        // Raw GRIDSS output is all BND; without conversion GRIDSS contributes 0 SVs.
        GRIDSS_CONVERT_BND(
            GRIDSS_CALL.out.vcf.join(GRIDSS_CALL.out.tbi)
        )
        ch_gridss_vcf = GRIDSS_CONVERT_BND.out.vcf
    } else {
        GRIDSS_STUB(ch_validated_bam)
        ch_gridss_vcf = GRIDSS_STUB.out.vcf
    }

    if (!params.skip_scramble) {
        SCRAMBLE_CALL(ch_validated_bam, ch_fasta, ch_fai)
        ch_scramble_vcf = SCRAMBLE_CALL.out.vcf
    } else {
        SCRAMBLE_STUB(ch_validated_bam)
        ch_scramble_vcf = SCRAMBLE_STUB.out.vcf
    }

    if (!params.skip_melt) {
        MELT_CALL(ch_validated_bam, ch_fasta, ch_fai)
        ch_melt_vcf = MELT_CALL.out.vcf
    } else {
        MELT_STUB(ch_validated_bam)
        ch_melt_vcf = MELT_STUB.out.vcf
    }

    // P7: SvABA local-assembly caller; 6th position in SUPP_VEC
    // SUPP_VEC positions: Manta[0] Delly[1] GRIDSS[2] Scramble[3] MELT[4] SvABA[5]
    if (!params.skip_svaba) {
        // ch_bwa_index stages the CLASSIC bwa index next to ref_fasta so SvABA's
        // bwa_idx_load_from_disk can find it. Without it SvABA dies on every run.
        SVABA_CALL(ch_validated_bam, ch_fasta, ch_fai, ch_bwa_index)
        ch_svaba_vcf = SVABA_CALL.out.vcf
    } else {
        SVABA_STUB(ch_validated_bam)
        ch_svaba_vcf = SVABA_STUB.out.vcf
    }

    // Collect VCFs per sample and merge with JASMINE (min_support=1)
    // File order: [manta, delly, gridss, scramble, melt, svaba] — determines SUPP_VEC bit positions
    ch_to_merge = MANTA_CALL.out.vcf
        .join(DELLY_MERGE.out.vcf)
        .join(ch_gridss_vcf)
        .join(ch_scramble_vcf)
        .join(ch_melt_vcf)
        .join(ch_svaba_vcf)
        .map { meta, manta_vcf, delly_vcf, gridss_vcf, scramble_vcf, melt_vcf, svaba_vcf ->
            [meta, [manta_vcf, delly_vcf, gridss_vcf, scramble_vcf, melt_vcf, svaba_vcf]]
        }

    JASMINE_MERGE(ch_to_merge, ch_fasta, ch_fai)

    // Rebuild cross-caller support for translocations (Jasmine leaves all TRA at
    // SUPP=1). Produces the published final sv_merged.vcf.gz.
    TRA_CONSENSUS(JASMINE_MERGE.out.vcf.join(JASMINE_MERGE.out.tbi))

    emit:
    sv_vcf      = TRA_CONSENSUS.out.vcf
    sv_tbi      = TRA_CONSENSUS.out.tbi
    str_vcf     = EXPANSIONHUNTER.out.vcf
    strling_tsv = ch_strling_tsv
}
