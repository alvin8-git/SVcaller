include { CIRCOS_PLOT   } from '../modules/pycirclize/plot'
include { TRUVARI_BENCH } from '../modules/truvari/bench'
include { MULTIQC       } from '../modules/multiqc/report'

process BUILD_HTML_REPORT {
    tag "${meta.id}"
    label 'process_single'
    container 'svcaller/utils:1.0'

    input:
    tuple val(meta), path(sv_tsv), path(cnv_bed), path(smn_tsv),
                     path(circos_svg), path(benchmark_json), path(sizebin_json),
                     path(coverage_summary), path(picard_metrics), path(str_vcf),
                     path(flagstat_txt), path(insert_size_metrics)

    output:
    tuple val(meta), path("${meta.id}.report.html"), emit: html

    script:
    def bench_arg    = benchmark_json.name   != "NO_FILE" ? "--benchmark ${benchmark_json}"   : ""
    def sizebin_arg  = sizebin_json.name     != "NO_FILE" ? "--sizebin   ${sizebin_json}"     : ""
    def cov_arg      = coverage_summary.name != "NO_FILE" ? "--coverage  ${coverage_summary}" : ""
    def met_arg      = picard_metrics.name   != "NO_FILE" ? "--metrics   ${picard_metrics}"   : ""
    def str_arg      = str_vcf.name          != "NO_STR"  ? "--str-vcf   ${str_vcf}"          : ""
    def flagstat_arg     = flagstat_txt.name        != "NO_FILE" ? "--flagstat    ${flagstat_txt}"        : ""
    def insert_size_arg  = insert_size_metrics.name != "NO_FILE" ? "--insert-size ${insert_size_metrics}" : ""
    """
    smn_report.py \\
        --tsv    ${smn_tsv} \\
        --sample ${meta.id} \\
        --out    ${meta.id}.smn_section.html

    html_report.py \\
        --sample           ${meta.id} \\
        --smn-html         ${meta.id}.smn_section.html \\
        --cnv-bed          ${cnv_bed} \\
        --sv-tsv           ${sv_tsv} \\
        --circos-svg       ${circos_svg} \\
        --out              ${meta.id}.report.html \\
        --pipeline-version ${workflow.manifest.version} \\
        ${bench_arg} \\
        ${sizebin_arg} \\
        ${cov_arg} \\
        ${met_arg} \\
        ${flagstat_arg} \\
        ${insert_size_arg} \\
        ${str_arg}
    """
}

workflow REPORT {
    take:
    ch_sv_tsv        // [ meta, tsv ]
    ch_cnv_bed       // [ meta, bed ]
    ch_smn_tsv       // [ meta, tsv ]
    ch_sv_vcf        // [ meta, vcf.gz ] for Circos
    ch_sv_tbi        // [ meta, vcf.gz.tbi ]
    ch_str_vcf       // [ meta, vcf ] ExpansionHunter STR for Circos ring 4 + HTML section 7
    ch_cytobands     // path
    ch_truth_vcf     // optional path (Channel.empty() if no truth)
    ch_truth_tbi     // optional path (Channel.empty() if no truth)
    ch_truth_bed     // optional path (Channel.empty() if no truth)
    ch_multiqc_files // collected QC files for MultiQC (may be empty)
    ch_coverage      // [ meta, mosdepth_summary ] for HTML section 2
    ch_metrics       // [ meta, picard_metrics ]   for HTML section 2
    ch_flagstat      // [ meta, flagstat_txt ]      for HTML section 2
    ch_insert_size   // [ meta, insert_size_metrics ] for HTML section 2

    main:
    ch_circos_in = ch_sv_vcf
        .join(ch_cnv_bed)
        .join(ch_str_vcf, remainder: true)
        .map { meta, sv, cnv, str -> [meta, sv, cnv, str ?: file("NO_STR")] }
    CIRCOS_PLOT(ch_circos_in, ch_cytobands)

    MULTIQC(ch_multiqc_files.ifEmpty([]))

    ch_bench   = Channel.empty()
    ch_sizebin = Channel.empty()
    if (params.giab_truth) {
        ch_sv_for_bench = ch_sv_vcf
            .join(ch_sv_tbi)
            .map { meta, vcf, tbi -> [meta, vcf, tbi] }
        TRUVARI_BENCH(ch_sv_for_bench, ch_truth_vcf, ch_truth_tbi, ch_truth_bed)
        ch_bench   = TRUVARI_BENCH.out.summary
        ch_sizebin = TRUVARI_BENCH.out.sizebin
    }

    ch_report_in = ch_sv_tsv
        .join(ch_cnv_bed)
        .join(ch_smn_tsv)
        .join(CIRCOS_PLOT.out.svg)
        .join(ch_bench.ifEmpty { [[:], file("NO_FILE")] }, remainder: true)
        .map { meta, sv, cnv, smn, svg, bench ->
            [meta, sv, cnv, smn, svg, bench ?: file("NO_FILE")]
        }
        .join(ch_sizebin.ifEmpty { [[:], file("NO_FILE")] }, remainder: true)
        .map { meta, sv, cnv, smn, svg, bench, sizebin ->
            [meta, sv, cnv, smn, svg, bench, sizebin ?: file("NO_FILE")]
        }
        .join(ch_coverage.ifEmpty { [[:], file("NO_FILE")] }, remainder: true)
        .map { meta, sv, cnv, smn, svg, bench, sizebin, cov ->
            [meta, sv, cnv, smn, svg, bench, sizebin, cov ?: file("NO_FILE")]
        }
        .join(ch_metrics.ifEmpty { [[:], file("NO_FILE")] }, remainder: true)
        .map { meta, sv, cnv, smn, svg, bench, sizebin, cov, met ->
            [meta, sv, cnv, smn, svg, bench, sizebin, cov, met ?: file("NO_FILE")]
        }
        .join(ch_str_vcf.ifEmpty { [[:], file("NO_STR")] }, remainder: true)
        .map { meta, sv, cnv, smn, svg, bench, sizebin, cov, met, str ->
            [meta, sv, cnv, smn, svg, bench, sizebin, cov, met, str ?: file("NO_STR")]
        }
        .join(ch_flagstat.ifEmpty { [[:], file("NO_FILE")] }, remainder: true)
        .map { meta, sv, cnv, smn, svg, bench, sizebin, cov, met, str, flagstat ->
            [meta, sv, cnv, smn, svg, bench, sizebin, cov, met, str, flagstat ?: file("NO_FILE")]
        }
        .join(ch_insert_size.ifEmpty { [[:], file("NO_FILE")] }, remainder: true)
        .map { meta, sv, cnv, smn, svg, bench, sizebin, cov, met, str, flagstat, insert_size ->
            [meta, sv, cnv, smn, svg, bench, sizebin, cov, met, str, flagstat, insert_size ?: file("NO_FILE")]
        }

    BUILD_HTML_REPORT(ch_report_in)

    emit:
    html    = BUILD_HTML_REPORT.out.html
    multiqc = MULTIQC.out.html.ifEmpty([])
}
