// M8 — alpha-globin (HBA1/HBA2) measurement.
//
// Its OWN subworkflow, not a fifth cnv_trait: two of the four channels are a
// junction search and a targeted pileup, which do not fit TRAIT_DEPTH's
// one-depth-pass-then-interpret shape. It does reuse bin/cnv_traits_common.py
// for normalized depth rather than growing a second copy of that logic, and it
// consumes the same CTRL_* windows from assets/cnv_trait_regions.bed.
//
//   ch1  HBA_DEPTH     mosdepth over the 5 diagnostic segments + CTRL windows
//   ch2  ALPHA_GLOBIN  allele naming, inside the integrator
//   ch3  HBA_JUNCTION  split-read / discordant-pair breakpoint search
//   ch4  HBA_SITES     targeted pileup at the pinned panel positions
//
// SVcaller MEASURES. It does not interpret: no HbH / Bart's / trait
// classification, no couple-level risk. OmniGen owns every clinical statement.

process HBA_DEPTH {
    tag "${meta.id}"
    label 'process_low'
    // container assigned in conf/docker.config: mosdepth biocontainer

    input:
    tuple val(meta), path(bam), path(bai)
    path segments_bed        // assets/hba_segments.bed
    path trait_regions_bed   // assets/cnv_trait_regions.bed (for CTRL_* windows)

    output:
    tuple val(meta), path("${meta.id}.hba_depth.regions.bed.gz"), emit: depth
    path "versions.yml",                                          emit: versions

    script:
    """
    # hba_segments.bed carries baseline/reliability in cols 5-6 and a trailing
    # '# note'. mosdepth wants a plain BED, so project to 4 columns; the extra
    # columns are re-joined by hba_depth.py from the same asset.
    awk 'BEGIN{OFS="\\t"} !/^#/ && NF>=4 {print \$1,\$2,\$3,\$4}' ${segments_bed} > hba_regions.bed
    # CTRL_* only: the trait loci (RHD/AMY1/...) are irrelevant here and each one
    # costs an extra region scan.
    awk 'BEGIN{OFS="\\t"} !/^#/ && \$4 ~ /^CTRL/ {print \$1,\$2,\$3,\$4}' ${trait_regions_bed} >> hba_regions.bed
    sort -k1,1 -k2,2n hba_regions.bed -o hba_regions.bed

    # --mapq 0 for the same reason as TRAIT_DEPTH: HBA1 and HBA2 are ~99%
    # identical, so MAPQ-filtering this locus destroys the copy-number signal
    # outright. The CTRL_* windows are measured the same way, and the baselines
    # in hba_segments.bed col 5 were calibrated against exactly this treatment.
    mosdepth \\
        --threads ${task.cpus} \\
        --no-per-base \\
        --mapq 0 \\
        --by hba_regions.bed \\
        ${meta.id}.hba_depth \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/mosdepth //')
    END_VERSIONS
    """
}

process HBA_DEPTH_CALL {
    tag "${meta.id}"
    label 'process_single'
    container params.utils_container ?: 'svcaller/utils:1.2'
    publishDir "${params.outdir}/${meta.id}/alpha_globin", mode: 'copy', pattern: "*.alpha_depth.tsv"

    input:
    tuple val(meta), path(depth_bed)
    path segments_bed

    output:
    tuple val(meta), path("${meta.id}.alpha_depth.tsv"), emit: tsv

    script:
    """
    export PATH=${projectDir}/bin:\$PATH
    hba_depth.py \\
        --depth    ${depth_bed} \\
        --segments ${segments_bed} \\
        --sample   ${meta.id} \\
        --out      ${meta.id}.alpha_depth.tsv
    """
}

process HBA_JUNCTION {
    tag "${meta.id}"
    label 'process_single'
    container params.utils_container ?: 'svcaller/utils:1.2'
    publishDir "${params.outdir}/${meta.id}/alpha_globin", mode: 'copy', pattern: "*.alpha_junction.tsv"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.alpha_junction.tsv"), emit: tsv

    script:
    """
    export PATH=${projectDir}/bin:\$PATH
    # A header-only file is a VALID NEGATIVE (no junction found) and must not be
    # confused with a crashed channel. The script exits non-zero on real errors;
    # there is deliberately no '|| touch' fallback here.
    hba_junction.py \\
        --bam    ${bam} \\
        --region ${params.hba_region ?: 'chr16:1-250000'} \\
        --sample ${meta.id} \\
        --out    ${meta.id}.alpha_junction.tsv
    """
}

process HBA_SITES {
    tag "${meta.id}"
    label 'process_single'
    container params.utils_container ?: 'svcaller/utils:1.2'
    publishDir "${params.outdir}/${meta.id}/alpha_globin", mode: 'copy', pattern: "*.alpha_sites.tsv"

    input:
    tuple val(meta), path(bam), path(bai), path(depth_tsv)
    path panel
    path alleles
    path fasta
    path fai

    output:
    tuple val(meta), path("${meta.id}.alpha_sites.tsv"), emit: tsv

    script:
    """
    export PATH=${projectDir}/bin:\$PATH
    # Zygosity at these sites is copy-number dependent: on a --SEA background the
    # surviving HBA2 is hemizygous, so a real variant sits near 100% VAF, not
    # 50%. The alpha-gene count is therefore an INPUT here.
    #
    # That count is produced by CHANNEL 2, not channel 1 — the per-segment depth
    # TSV has no gene-count column, only per-segment scores. So we ask
    # alpha_globin.py for it directly rather than re-deriving it here, which
    # would duplicate the allele matcher and let the two drift apart.
    #
    # On any failure this yields NA, and hba_sites.py then degrades to a
    # VAF-only report (zygosity_basis=vaf_only) rather than silently assuming a
    # 2-copy diploid background.
    ALPHA_GENES=\$(alpha_globin.py \\
                      --sample ${meta.id} \\
                      --depth  ${depth_tsv} \\
                      --alleles ${alleles} \\
                      --print-alpha-genes 2>/dev/null || echo NA)
    ALPHA_GENES=\${ALPHA_GENES:-NA}

    hba_sites.py \\
        --bam         ${bam} \\
        --panel       ${panel} \\
        --ref         ${fasta} \\
        --sample      ${meta.id} \\
        --alpha-genes "\${ALPHA_GENES}" \\
        --out         ${meta.id}.alpha_sites.tsv
    """
}

process ALPHA_GLOBIN_INTEGRATE {
    tag "${meta.id}"
    label 'process_single'
    container params.utils_container ?: 'svcaller/utils:1.2'
    publishDir "${params.outdir}/${meta.id}/alpha_globin", mode: 'copy', pattern: "*.alpha_globin.tsv"

    input:
    tuple val(meta), path(depth_tsv), path(junction_tsv), path(sites_tsv)
    path alleles
    path panel

    output:
    tuple val(meta), path("${meta.id}.alpha_globin.tsv"), emit: tsv

    script:
    """
    export PATH=${projectDir}/bin:\$PATH
    alpha_globin.py \\
        --sample   ${meta.id} \\
        --depth    ${depth_tsv} \\
        --junction ${junction_tsv} \\
        --sites    ${sites_tsv} \\
        --alleles  ${alleles} \\
        --panel    ${panel} \\
        --out      ${meta.id}.alpha_globin.tsv
    """
}

workflow ALPHA_GLOBIN {
    take:
    ch_bam             // [ meta, bam, bai ]
    ch_segments        // value  assets/hba_segments.bed
    ch_trait_regions   // value  assets/cnv_trait_regions.bed
    ch_panel           // value  assets/hba_pathogenic_sites.tsv
    ch_alleles         // value  assets/hba_deletion_alleles.tsv
    ch_fasta           // value  reference FASTA
    ch_fai             // value  reference .fai

    main:
    // Every shared reference/asset above MUST arrive as Channel.value(). A queue
    // channel is consumed by the first process that reads it and every later one
    // silently receives nothing — see CLAUDE.md "Nextflow channel exhaustion".

    HBA_DEPTH(ch_bam, ch_segments, ch_trait_regions)
    HBA_DEPTH_CALL(HBA_DEPTH.out.depth, ch_segments)
    HBA_JUNCTION(ch_bam)

    // Channel 4 consumes channel 1's output: VAF is uninterpretable without the
    // alpha-gene count. Joined on the FULL meta map, which is safe here because
    // both sides descend from the same ch_bam without any added key.
    HBA_SITES(ch_bam.join(HBA_DEPTH_CALL.out.tsv), ch_panel, ch_alleles,
              ch_fasta, ch_fai)

    // remainder:true + a NO_* sentinel, never an exact .join(): an exact join on
    // the full meta map silently drops the whole sample when one channel is
    // missing, producing no output and no error.
    ch_in = HBA_DEPTH_CALL.out.tsv
        .join(HBA_JUNCTION.out.tsv, remainder: true)
        .join(HBA_SITES.out.tsv,    remainder: true)
        .filter { it[1] != null }   // depth absent => nothing to integrate
        .map { meta, depth, junction, sites ->
            [meta, depth, junction ?: file("NO_JUNCTION"), sites ?: file("NO_SITES")]
        }

    ALPHA_GLOBIN_INTEGRATE(ch_in, ch_alleles, ch_panel)

    emit:
    tsv      = ALPHA_GLOBIN_INTEGRATE.out.tsv
    depth    = HBA_DEPTH_CALL.out.tsv
    junction = HBA_JUNCTION.out.tsv
    sites    = HBA_SITES.out.tsv
}
