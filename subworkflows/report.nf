include { CIRCOS_PLOT                      } from '../modules/pycirclize/plot'
include { TRUVARI_BENCH                    } from '../modules/truvari/bench'
include { TRUVARI_BENCH as TRUVARI_BENCH_V5Q } from '../modules/truvari/bench'
include { MULTIQC                          } from '../modules/multiqc/report'

process BUILD_HTML_REPORT {
    tag "${meta.id}"
    label 'process_single'
    container 'svcaller/utils:1.1'
    publishDir "${params.outdir}/${meta.id}", mode: 'copy', pattern: "*.report.html"

    input:
    tuple val(meta), path(sv_tsv), path(cnv_bed), path(smn_tsv),
                     path(circos_svg), path(benchmark_json), path(sizebin_json),
                     path(coverage_summary), path(picard_metrics), path(str_vcf),
                     path(flagstat_txt), path(insert_size_metrics),
                     path(benchmark_json_v5q), path(sizebin_json_v5q),
                     path(strling_tsv)

    output:
    tuple val(meta), path("${meta.id}.report.html"), emit: html

    script:
    def bench_arg       = benchmark_json.name      != "NO_BENCH"      ? "--benchmark     ${benchmark_json}"      : ""
    def sizebin_arg     = sizebin_json.name        != "NO_SIZEBIN"    ? "--sizebin       ${sizebin_json}"        : ""
    def cov_arg         = coverage_summary.name    != "NO_COV"        ? "--coverage      ${coverage_summary}"    : ""
    def met_arg         = picard_metrics.name      != "NO_METRICS"    ? "--metrics       ${picard_metrics}"      : ""
    def str_arg         = str_vcf.name             != "NO_STR"        ? "--str-vcf       ${str_vcf}"             : ""
    def flagstat_arg    = flagstat_txt.name        != "NO_FLAGSTAT"   ? "--flagstat      ${flagstat_txt}"        : ""
    def insert_size_arg = insert_size_metrics.name != "NO_INSERT"     ? "--insert-size   ${insert_size_metrics}" : ""
    def bench_v5q_arg   = benchmark_json_v5q.name  != "NO_BENCH_V5Q"  ? "--benchmark-v5q ${benchmark_json_v5q}"  : ""
    def sizebin_v5q_arg = sizebin_json_v5q.name    != "NO_SBN_V5Q"    ? "--sizebin-v5q   ${sizebin_json_v5q}"    : ""
    def strling_arg     = strling_tsv.name         != "NO_STRLING"    ? "--strling-tsv   ${strling_tsv}"         : ""
    """
    export PATH=${projectDir}/bin:\$PATH
    smn_report.py \\
        --tsv    ${smn_tsv} \\
        --sample ${meta.id} \\
        --out    ${meta.id}.smn_section.html

    # v15: NORMAL filter, INREPEAT badge, INTERMEDIATE fix
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
        ${bench_v5q_arg} \\
        ${sizebin_v5q_arg} \\
        ${cov_arg} \\
        ${met_arg} \\
        ${flagstat_arg} \\
        ${insert_size_arg} \\
        ${str_arg} \\
        ${strling_arg}
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
    ch_depth_bed     // [ meta, regions.bed.gz ]    mosdepth 50kb windows for Circos ring A
    ch_annotsv_tsv   // [ meta, annotated.tsv ]     raw AnnotSV output for Circos rings B/C
    ch_strling_tsv   // [ meta, tsv ] STRling genome-wide STR expansions (optional)

    main:
    ch_circos_in = ch_sv_vcf
        .join(ch_cnv_bed)
        .join(ch_str_vcf, remainder: true)
        .filter { it[1] != null }   // drop tuples fired before sv_vcf is ready (str_vcf caches before Jasmine)
        .join(ch_depth_bed, remainder: true)
        .join(ch_annotsv_tsv, remainder: true)
        .filter { it[1] != null }   // drop samples re-introduced by depth/annotsv remainders when cnv_bed absent
        .join(ch_strling_tsv, remainder: true)
        .filter { it[1] != null }
        .map { meta, sv, cnv, str, depth, annotsv, strling ->
            [meta, sv, cnv, str ?: file("NO_STR"), depth ?: file("NO_DEPTH"),
             annotsv ?: file("NO_ANNOTSV"), strling ?: file("NO_STRLING")]
        }
    CIRCOS_PLOT(ch_circos_in, ch_cytobands)

    MULTIQC(ch_multiqc_files.ifEmpty([]))

    ch_bench   = Channel.empty()
    ch_sizebin = Channel.empty()
    if (params.giab_truth) {
        ch_sv_for_bench = ch_sv_vcf
            .join(ch_sv_tbi)
            .map { meta, vcf, tbi -> [meta, vcf, tbi] }
        TRUVARI_BENCH(ch_sv_for_bench, ch_truth_vcf, ch_truth_tbi, ch_truth_bed, "T2T")
        ch_bench   = TRUVARI_BENCH.out.summary
        ch_sizebin = TRUVARI_BENCH.out.sizebin
    }

    ch_bench_v5q   = Channel.empty()
    ch_sizebin_v5q = Channel.empty()
    if (params.giab_truth_v5q) {
        ch_truth_v5q     = Channel.fromPath(params.giab_truth_v5q, checkIfExists: true)
        ch_truth_v5q_tbi = Channel.fromPath("${params.giab_truth_v5q}.tbi", checkIfExists: true)
        ch_truth_v5q_bed = Channel.fromPath("${params.giab_truth_v5q}".replaceAll(/\.vcf\.gz$/, '.bed'), checkIfExists: true)
        ch_sv_for_v5q = ch_sv_vcf
            .join(ch_sv_tbi)
            .map { meta, vcf, tbi -> [meta, vcf, tbi] }
        TRUVARI_BENCH_V5Q(ch_sv_for_v5q, ch_truth_v5q, ch_truth_v5q_tbi, ch_truth_v5q_bed, "v5q")
        ch_bench_v5q   = TRUVARI_BENCH_V5Q.out.summary
        ch_sizebin_v5q = TRUVARI_BENCH_V5Q.out.sizebin
    }

    // Mandatory channels joined first (inner join). Optional channels last (remainder: true).
    // filter { it[1] != null } drops spurious right-side remainders from timing races.
    ch_report_in = ch_sv_tsv
        .join(ch_cnv_bed)
        .join(ch_smn_tsv)
        .join(CIRCOS_PLOT.out.svg)
        .join(ch_coverage)
        .join(ch_metrics)
        .join(ch_str_vcf)
        .join(ch_flagstat)
        .join(ch_insert_size)
        // tuple: [meta, sv, cnv, smn, svg, cov, met, str, flagstat, ins]
        .join(ch_bench, remainder: true)
        .filter   { it[1] != null }
        .map { meta, sv, cnv, smn, svg, cov, met, str, flagstat, ins, bench ->
            [meta, sv, cnv, smn, svg, bench ?: file("NO_BENCH"), cov, met, str, flagstat, ins]
        }
        // tuple: [meta, sv, cnv, smn, svg, bench, cov, met, str, flagstat, ins]
        .join(ch_sizebin, remainder: true)
        .filter   { it[1] != null }
        .map { meta, sv, cnv, smn, svg, bench, cov, met, str, flagstat, ins, sizebin ->
            [meta, sv, cnv, smn, svg, bench, sizebin ?: file("NO_SIZEBIN"), cov, met, str, flagstat, ins]
        }
        // tuple: [meta, sv, cnv, smn, svg, bench, sizebin, cov, met, str, flagstat, ins]
        .join(ch_bench_v5q, remainder: true)
        .filter   { it[1] != null }
        .map { meta, sv, cnv, smn, svg, bench, sizebin, cov, met, str, flagstat, ins, bench_v5q ->
            [meta, sv, cnv, smn, svg, bench, sizebin, cov, met, str, flagstat, ins, bench_v5q ?: file("NO_BENCH_V5Q")]
        }
        // tuple: [meta, ..., bench_v5q]
        .join(ch_sizebin_v5q, remainder: true)
        .filter   { it[1] != null }
        .map { meta, sv, cnv, smn, svg, bench, sizebin, cov, met, str, flagstat, ins, bench_v5q, sizebin_v5q ->
            [meta, sv, cnv, smn, svg, bench, sizebin, cov, met, str, flagstat, ins,
             bench_v5q, sizebin_v5q ?: file("NO_SBN_V5Q")]
        }
        // tuple: [meta, ..., sizebin_v5q]
        .join(ch_strling_tsv, remainder: true)
        .filter   { it[1] != null }
        .map { meta, sv, cnv, smn, svg, bench, sizebin, cov, met, str, flagstat, ins, bench_v5q, sizebin_v5q, strling ->
            [meta, sv, cnv, smn, svg, bench, sizebin, cov, met, str, flagstat, ins,
             bench_v5q, sizebin_v5q, strling ?: file("NO_STRLING")]
        }
        // final: [meta, sv_tsv, cnv_bed, smn_tsv, circos_svg, benchmark_json, sizebin_json,
        //         coverage_summary, picard_metrics, str_vcf, flagstat_txt, insert_size_metrics,
        //         benchmark_json_v5q, sizebin_json_v5q, strling_tsv]

    BUILD_HTML_REPORT(ch_report_in)

    emit:
    html    = BUILD_HTML_REPORT.out.html
    multiqc = MULTIQC.out.html.ifEmpty([])
}
