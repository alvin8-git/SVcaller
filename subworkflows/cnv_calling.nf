// Per-chromosome scatter assessment (optimization #4) — both callers stay MONOLITHIC.
//
// CNVpytor: NOT safely scatterable. The `-his` step fits the diploid RD baseline
// genome-wide and builds one GC-correction table across the whole genome; every
// call in the TSV is a copy number normalized to that global baseline (2.0 =
// diploid, see cnv_consensus.load_cnvpytor col 3). Independent per-chromosome
// roots would each self-normalize to their own median, which is provably wrong on
// the sex chromosomes: a male chrX/chrY is haploid, so a per-chrom run makes its
// own coverage the CN=2 level and calls the whole chromosome neutral. CNVpytor
// also runs `-rd` over the full BAM (no interval restriction), so it sees
// chrX/Y/M. There is no supported way to merge per-chrom .pytor roots before a
// single global `-his`, and parallel tasks writing one root would corrupt it. The
// only real speedup here is CNVpytor's own within-task `-j` jobs flag, which is a
// different change and version-dependent. Left monolithic.
//
// GATK gCNV: NOT safely scatterable against this PON. CollectReadCounts alone is
// interval-shardable, but DenoiseReadCounts requires the sample counts to carry
// the EXACT interval set the PON was built on (giab_cnv_pon.hdf5 is genome-wide,
// bin-length 1000). A single-chromosome counts file fails GATK's interval match,
// and GATK ships no tool to concatenate per-chrom read-count HDF5 back into one
// genome-wide file. So the denoise/model/call stages cannot shard, and the count
// shard cannot be gathered for them. Left monolithic.
//
// tests/test_cnv_scatter_contract.py locks the TSV contract and the gather-by-
// concat semantics either caller would have to preserve, so a future attempt
// fails a test instead of silently emptying a clinical CNV sheet.
include { CNVPYTOR_CALL           } from '../modules/cnvpytor/call'
include { GATK_PREPROCESS_INTERVALS } from '../modules/gatk/gcnv_pon'
include { GATK_GCNV_CALL           } from '../modules/gatk/gcnv_call'

process CNV_CONSENSUS {
    tag "${meta.id}"
    label 'process_single'
    container params.utils_container ?: 'svcaller/utils:1.3'
    publishDir "${params.outdir}/${meta.id}", mode: 'copy', pattern: "*.cnv_consensus.bed"

    input:
    tuple val(meta), path(cnvpytor_tsv), path(gatk_tsv)

    output:
    tuple val(meta), path("${meta.id}.cnv_consensus.bed"), emit: bed

    script:
    """
    cnv_consensus.py \\
        --cnvpytor ${cnvpytor_tsv} \\
        --gatk     ${gatk_tsv} \\
        --sample   ${meta.id} \\
        --out      ${meta.id}.cnv_consensus.bed
    """
}

workflow CNV_CALLING {
    take:
    ch_bam       // [ meta, bam, bai ]
    ch_fasta     // path
    ch_fai       // path
    ch_dict      // path (.dict)
    ch_pon       // path (PoN HDF5, may be null)
    ch_intervals // path

    main:
    CNVPYTOR_CALL(ch_bam, ch_fasta)

    // Preprocess intervals to bin-length 1000 interval_list (must match PON build)
    GATK_PREPROCESS_INTERVALS(ch_fasta, ch_fai, ch_dict, ch_intervals)
    // .first() makes the single interval_list a reusable value channel; without it the
    // per-sample ch_bam pairs against a 1-element queue and GATK_GCNV_CALL silently runs
    // for only the first sample (ch_pon is now a value channel for the same reason).
    GATK_GCNV_CALL(ch_bam, ch_fasta, ch_fai, ch_dict, ch_pon, GATK_PREPROCESS_INTERVALS.out.preprocessed.first())

    // Join the awk-converted GATK TSV (columns CONTIG/START/END/CALL_COPY_NUMBER/QUALITY),
    // NOT the raw .seg — cnv_consensus.py reads CALL_COPY_NUMBER, absent from the .seg.
    ch_for_consensus = CNVPYTOR_CALL.out.tsv
        .join(GATK_GCNV_CALL.out.tsv.map { meta, tsv -> [meta, tsv] })

    CNV_CONSENSUS(ch_for_consensus)

    emit:
    cnv_bed = CNV_CONSENSUS.out.bed
}
