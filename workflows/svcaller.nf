include { PREPROCESS      } from '../subworkflows/preprocess'
include { SV_CALLING      } from '../subworkflows/sv_calling'
include { CNV_CALLING     } from '../subworkflows/cnv_calling'
include { CNV_TRAITS      } from '../subworkflows/cnv_traits'
include { SMN_CALLING     } from '../subworkflows/smn_calling'
include { ALPHA_GLOBIN    } from '../subworkflows/alpha_globin'
include { ANNOTATE        } from '../subworkflows/annotate'
include { REPORT          } from '../subworkflows/report'
include { SV_PON_ANNOTATE } from '../modules/sv_pon/annotate'
include { SAMTOOLS_FILTER_CHROMS } from '../modules/samtools/filter_chroms'
include { VALIDATE_REF_BAM       } from '../modules/samtools/validate_ref_bam'

workflow SVCALLER {
    take:
    ch_input      // parsed samplesheet channel
    ch_fasta
    ch_fai
    ch_bwt_index
    ch_bwa_index  // classic bwa index files for SvABA (.amb/.ann/.bwt/.pac/.sa)
    ch_dict
    ch_pon
    ch_intervals
    ch_annotsv_db
    ch_cytobands
    ch_eh_catalog
    ch_trait_regions
    ch_hba_segments
    ch_hba_panel
    ch_hba_alleles

    main:
    // M1: Preprocess
    PREPROCESS(ch_input, ch_fasta, ch_fai, ch_bwt_index, params.min_depth)

    ch_bam = PREPROCESS.out.bam

    // Canonical-chromosome filter + ref/BAM validation, lifted here from inside
    // SV_CALLING so CNV, SMN, and CNV-traits run on the SAME validated BAM as the
    // SV callers instead of the raw one. For a full-hg38 BAM input the raw BAM
    // carries 3366 contigs; CNVpytor's genome-wide depth baseline is then biased by
    // low-coverage alt/decoy contigs (canonical bins read high). FASTQ-derived BAMs
    // are aligned to hg38.canonical.fa (no alt contigs) and skip the filter, so this
    // is a no-op for them. It runs once (storeDir-cached) and the SV branch dominates
    // the critical path, so nothing downstream starts later.
    ch_bam.branch {
        needs_filter: it[0].get('needs_chr_filter', true)
        canonical:    true
    }.set { ch_bam_branched }
    SAMTOOLS_FILTER_CHROMS(ch_bam_branched.needs_filter, ch_fai)
    ch_validated_bam = SAMTOOLS_FILTER_CHROMS.out.bam
        .mix(ch_bam_branched.canonical)
    VALIDATE_REF_BAM(ch_validated_bam, ch_fai)
    ch_validated_bam = VALIDATE_REF_BAM.out.bam

    // M2 + M3 + M4: run in parallel on the validated BAM
    SV_CALLING(ch_validated_bam, ch_fasta, ch_fai, ch_eh_catalog, ch_bwa_index)
    CNV_CALLING(ch_validated_bam, ch_fasta, ch_fai, ch_dict, ch_pon, ch_intervals)
    SMN_CALLING(ch_validated_bam, ch_fasta, ch_fai)

    // Copy-number / blood-group traits: targeted normalized read depth + consensus
    // corroboration → per-sample OmniGen contract files (Rh/RHD, AMY1, GST-null, LPA KIV-2)
    CNV_TRAITS(ch_validated_bam, ch_trait_regions, CNV_CALLING.out.cnv_bed)

    // M8: alpha-globin (HBA1/HBA2). Runs on the RAW ch_bam, DELIBERATELY, not the
    // validated BAM the other modules now use. Its GIAB depth baselines
    // (validation/giab_alpha_baseline.tsv, assets/hba_segments.bed col 5) were
    // calibrated on the raw BAM; querying the filtered BAM against raw-derived
    // baselines would drift the score = ratio / baseline silently. Moving alpha to
    // the validated BAM requires re-deriving those baselines on filtered GIAB BAMs
    // first (a separate, validated change). It MEASURES only: alpha-gene dosage,
    // deletion alleles, targeted site genotypes, and the scope screened. It does not
    // interpret; OmniGen discovers the contract at
    // results/<S>/alpha_globin/<S>.alpha_globin.tsv.
    if (!params.skip_alpha_globin) {
        ALPHA_GLOBIN(ch_bam, ch_hba_segments, ch_trait_regions, ch_hba_panel,
                     ch_hba_alleles, ch_fasta, ch_fai)
        ch_alpha_globin = ALPHA_GLOBIN.out.tsv
    } else {
        ch_alpha_globin = Channel.empty()
    }

    // P3: Optional GIAB SV PON annotation — annotates Jasmine VCF with SV_PON=1
    // before AnnotSV so the flag is preserved in the TSV INFO column.
    if (params.sv_pon) {
        ch_pon_bed = Channel.value(file(params.sv_pon, checkIfExists: true))
        SV_PON_ANNOTATE(
            SV_CALLING.out.sv_vcf.join(SV_CALLING.out.sv_tbi),
            ch_pon_bed
        )
        ch_sv_for_annotate = SV_PON_ANNOTATE.out.vcf
        ch_sv_tbi_for_report = SV_PON_ANNOTATE.out.tbi
    } else {
        ch_sv_for_annotate   = SV_CALLING.out.sv_vcf
        ch_sv_tbi_for_report = SV_CALLING.out.sv_tbi
    }

    // M5: Annotate SVs
    ANNOTATE(ch_sv_for_annotate, ch_annotsv_db)

    // Optional truvari truth channels
    ch_truth     = params.giab_truth
        ? Channel.fromPath(params.giab_truth, checkIfExists: true)
        : Channel.empty()
    ch_truth_tbi = params.giab_truth
        ? Channel.fromPath("${params.giab_truth}.tbi", checkIfExists: true)
        : Channel.empty()
    ch_truth_bed = params.giab_truth
        ? Channel.fromPath("${params.giab_truth}".replaceAll(/\.vcf\.gz$/, '.bed'), checkIfExists: true)
        : Channel.empty()

    // Collect QC files for MultiQC aggregation
    ch_multiqc_files = Channel.empty()
        .mix(PREPROCESS.out.fastqc_zip.map { meta, zips -> zips }.flatten())
        .mix(PREPROCESS.out.metrics.map { meta, m -> m })
        .mix(PREPROCESS.out.coverage.map { meta, s -> s })
        .mix(PREPROCESS.out.flagstat.map { meta, f -> f })
        .mix(PREPROCESS.out.insert_size.map { meta, i -> i })
        .collect()

    // M6 + M7: Visualize and report
    REPORT(
        ANNOTATE.out.tsv,
        CNV_CALLING.out.cnv_bed,
        SMN_CALLING.out.tsv,
        ch_sv_for_annotate,
        ch_sv_tbi_for_report,
        SV_CALLING.out.str_vcf,
        ch_cytobands,
        ch_truth,
        ch_truth_tbi,
        ch_truth_bed,
        ch_multiqc_files,
        PREPROCESS.out.coverage,
        PREPROCESS.out.metrics,
        PREPROCESS.out.flagstat,
        PREPROCESS.out.insert_size,
        PREPROCESS.out.regions_bed,
        ANNOTATE.out.annotated_tsv,
        SV_CALLING.out.strling_tsv,
        CNV_TRAITS.out.rh,
        CNV_TRAITS.out.amy1,
        CNV_TRAITS.out.gst,
        CNV_TRAITS.out.lpa,
        ch_alpha_globin,
    )

    emit:
    alpha_globin = ch_alpha_globin
    sv_vcf   = SV_CALLING.out.sv_vcf
    str_vcf  = SV_CALLING.out.str_vcf
    cnv_bed  = CNV_CALLING.out.cnv_bed
    smn_tsv  = SMN_CALLING.out.tsv
    html     = REPORT.out.html
    multiqc  = REPORT.out.multiqc
}
