process MANTA_RESIDUAL_REGIONS {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(vcf)
    path fai

    output:
    tuple val(meta), path("${meta.id}.manta_residual.bed"), emit: bed
    path "versions.yml",                                     emit: versions

    script:
    """
    # Extract non-PASS Manta calls; extend ±1 kb; merge overlapping intervals
    zcat ${vcf} | awk -v OFS='\\t' '
        /^#/ { next }
        \$7 == "PASS" { next }
        {
            end = \$2
            n = split(\$8, info, ";")
            for (i = 1; i <= n; i++) {
                if (info[i] ~ /^END=/) end = substr(info[i], 5) + 0
            }
            s = (\$2 - 1000 > 0) ? \$2 - 1000 : 0
            e = end + 1000
            print \$1, s, e
        }
    ' | sort -k1,1V -k2,2n | awk -v OFS='\\t' '
        NR == 1 { chr = \$1; s = \$2; e = \$3; next }
        \$1 == chr && \$2 <= e { if (\$3 > e) e = \$3; next }
        { print chr, s, e; chr = \$1; s = \$2; e = \$3 }
        END { if (NR > 0) print chr, s, e }
    ' > ${meta.id}.manta_residual.bed

    # If all Manta calls were PASS (ideal case), create a trivial region so
    # SAMTOOLS_SUBSET and GRIDSS_CALL processes don't fail on empty input
    if [ ! -s ${meta.id}.manta_residual.bed ]; then
        head -1 ${fai} | awk '{print \$1"\\t0\\t1000"}' > ${meta.id}.manta_residual.bed
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -1)
    END_VERSIONS
    """
}
