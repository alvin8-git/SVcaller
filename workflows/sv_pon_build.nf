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
    SURVIVOR_MERGE.out.vcf
        .collectFile(storeDir: "${params.outdir}") { vcf ->
            def bed = file("${params.outdir}/giab_sv_pon.bed")
            def lines = []
            vcf.eachLine { line ->
                if (line.startsWith("#")) return
                def fields = line.split("\t")
                def chrom = fields[0]
                def start = [0L, (fields[1] as long) - 500L].max()
                def end_pos = fields[7].split(";").find { it.startsWith("END=") }
                def end = end_pos ? ((end_pos.split("=")[1] as long) + 500L) : (start + 1000L)
                lines << "${chrom}\t${start}\t${end}\n"
            }
            ["giab_sv_pon.bed", lines.join("")]
        }

    SURVIVOR_MERGE.out.vcf.view { vcf ->
        log.info "SV PON VCF: ${vcf} — copy giab_sv_pon.bed from ${params.outdir} to use with --sv_pon"
    }
}
