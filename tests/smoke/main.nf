#!/usr/bin/env nextflow
// Smoke test: runs the pipeline's real ALPHA_GLOBIN (M8) and CNV_TRAITS (M3b)
// subworkflows on a ~12 MB GIAB HG002 slice (tests/fixtures/smoke/).
// Launch from the repo root with tests/smoke/run.sh, which also checks the outputs.
// ponytail: only the two targeted-depth subworkflows; alignment, SV/CNV/STR/SMN
// calling need a whole genome and are not covered here.
nextflow.enable.dsl = 2

include { ALPHA_GLOBIN } from '../../subworkflows/alpha_globin'
include { CNV_TRAITS   } from '../../subworkflows/cnv_traits'

workflow {
    def root = file("${projectDir}/../..").normalize()
    def fx   = "${root}/tests/fixtures/smoke"

    ch_bam = Channel.value([[id: 'HG002'],
                            file("${fx}/HG002.slice.bam", checkIfExists: true),
                            file("${fx}/HG002.slice.bam.bai", checkIfExists: true)])
    ch_trait_regions = Channel.value(file("${root}/assets/cnv_trait_regions.bed", checkIfExists: true))

    ALPHA_GLOBIN(ch_bam,
                 Channel.value(file("${root}/assets/hba_segments.bed", checkIfExists: true)),
                 ch_trait_regions,
                 Channel.value(file("${root}/assets/hba_pathogenic_sites.tsv", checkIfExists: true)),
                 Channel.value(file("${root}/assets/hba_deletion_alleles.tsv", checkIfExists: true)),
                 Channel.value(file("${fx}/chr16_1-250000.fa", checkIfExists: true)),
                 Channel.value(file("${fx}/chr16_1-250000.fa.fai", checkIfExists: true)))

    // No CNV consensus BED on a slice: the trait callers take NO_FILE, as in a run without one.
    CNV_TRAITS(ch_bam, ch_trait_regions, Channel.empty())
}
