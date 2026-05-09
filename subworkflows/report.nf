include { CIRCOS_PLOT   } from '../modules/pycirclize/plot'
include { TRUVARI_BENCH } from '../modules/truvari/bench'

process BUILD_HTML_REPORT {
    tag "${meta.id}"
    label 'process_single'
    container 'svcaller/utils:1.0'

    input:
    tuple val(meta), path(sv_tsv), path(cnv_bed), path(smn_tsv),
                     path(circos_svg), path(benchmark_json)

    output:
    tuple val(meta), path("${meta.id}.report.html"), emit: html

    script:
    def bench_arg = benchmark_json.name != "NO_FILE" ? "--benchmark ${benchmark_json}" : ""
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
        ${bench_arg}
    """
}

workflow REPORT {
    take:
    ch_sv_tsv        // [ meta, tsv ]
    ch_cnv_bed       // [ meta, bed ]
    ch_smn_tsv       // [ meta, tsv ]
    ch_sv_vcf        // [ meta, vcf.gz ] for Circos
    ch_cytobands     // path
    ch_truth_vcf     // optional path

    main:
    ch_circos_in = ch_sv_vcf.join(ch_cnv_bed)
    CIRCOS_PLOT(ch_circos_in, ch_cytobands)

    ch_bench = Channel.empty()
    if (params.giab_truth) {
        TRUVARI_BENCH(ch_sv_vcf, ch_truth_vcf)
        ch_bench = TRUVARI_BENCH.out.summary
    }

    ch_report_in = ch_sv_tsv
        .join(ch_cnv_bed)
        .join(ch_smn_tsv)
        .join(CIRCOS_PLOT.out.svg)
        .join(ch_bench.ifEmpty { [[:], file("NO_FILE")] }, remainder: true)
        .map { meta, sv, cnv, smn, svg, bench ->
            [meta, sv, cnv, smn, svg, bench ?: file("NO_FILE")]
        }

    BUILD_HTML_REPORT(ch_report_in)

    emit:
    html = BUILD_HTML_REPORT.out.html
}
