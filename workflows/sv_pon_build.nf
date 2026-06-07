/**
 * sv_pon_build.nf — Build a multi-sample SV Panel of Normals from GIAB samples.
 *
 * Runs Manta on each BAM (HG001, HG003-HG007; HG002 excluded to avoid circularity),
 * then merges all per-sample VCFs with SURVIVOR. Sites present in ≥2 samples are
 * written to a BED file for use as a recurrent-SV blacklist.
 *
 * Usage:
 *   nextflow run workflows/sv_pon_build.nf -profile docker \
 *     --input validation/pon_sv_samplesheet.csv \
 *     --ref_fasta /data/alvin/ref/GRCh38/hg38.canonical.fa \
 *     --outdir pon/sv_pon
 */

include { MANTA_CALL             } from '../modules/manta/call'
include { SAMTOOLS_FILTER_CHROMS } from '../modules/samtools/filter_chroms'
include { SURVIVOR_MERGE         } from '../modules/survivor/merge'

process SURVIVOR_VCF_TO_BED {
    label 'process_single'
    publishDir "${params.outdir}", mode: 'copy'
    container 'quay.io/biocontainers/survivor:1.0.7--hd03093a_2'

    input:
    path vcf

    output:
    path "giab_sv_pon.bed", emit: bed

    script:
    """
    awk '!/^#/ {
        split(\$8, info, ";")
        svtype = ""
        endval = \$2 + 1000
        for (i in info) {
            if (info[i] ~ /^SVTYPE=/) svtype = substr(info[i], 8)
            if (info[i] ~ /^END=/)    endval = substr(info[i], 5) + 0
        }
        if (svtype == "TRA" || svtype == "BND") next
        if (endval - \$2 > 10000000) endval = \$2 + 10000000
        start = (\$2 - 500 < 0) ? 0 : \$2 - 500
        print \$1 "\\t" start "\\t" (endval + 500)
    }' ${vcf} | sort -k1,1 -k2,2n > giab_sv_pon.bed
    """
}

params.input   = null
params.outdir  = 'pon/sv_pon'
params.min_sv_callers = 2   // min samples a site must appear in to enter the PON

workflow {
    ch_fasta = Channel.value(file(params.ref_fasta, checkIfExists: true))
    ch_fai   = Channel.value(file("${params.ref_fasta}.fai", checkIfExists: true))

    // Parse samplesheet: sample,bam columns (same format as giab_samplesheet.csv)
    ch_bam = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample, needs_chr_filter: true]
            def bam  = file(row.bam, checkIfExists: true)
            def bai  = file("${row.bam}.bai", checkIfExists: true)
            [meta, bam, bai]
        }

    // Filter to canonical chromosomes
    SAMTOOLS_FILTER_CHROMS(ch_bam, ch_fai)

    // Run Manta (fast; ~1h per sample) for PON sites
    MANTA_CALL(SAMTOOLS_FILTER_CHROMS.out.bam, ch_fasta, ch_fai)

    // Collect all per-sample Manta VCFs and merge across samples
    MANTA_CALL.out.vcf
        .collect { meta, vcf -> vcf.toString() }
        .map { paths ->
            def listFile = file("${workDir}/sv_pon_vcf_list.txt")
            listFile.text = paths.join("\n") + "\n"
            listFile
        }
        .set { ch_vcf_list }

    SURVIVOR_MERGE(ch_vcf_list, 1000, params.min_sv_callers, "giab_sv_pon")

    // Convert SURVIVOR VCF → BED (chrom, start-500, end+500) for bedtools intersect
    SURVIVOR_VCF_TO_BED(SURVIVOR_MERGE.out.vcf)

    SURVIVOR_VCF_TO_BED.out.bed.view { bed ->
        log.info "SV PON BED: ${bed} — use with --sv_pon ${params.outdir}/giab_sv_pon.bed"
    }
}
