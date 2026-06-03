process TRUVARI_BENCH {
    tag "${meta.id}/${truth_label}"
    label 'process_low'
    container 'quay.io/biocontainers/truvari:4.3.1--pyhdfd78af_0'

    input:
    tuple val(meta), path(query_vcf), path(query_tbi)
    path truth_vcf
    path truth_tbi
    path truth_bed
    val  truth_label   // e.g. "T2T" or "v06" — used to produce unique output filenames

    output:
    tuple val(meta), path("${meta.id}.${truth_label}.truvari_summary.json"), emit: summary
    tuple val(meta), path("${meta.id}.${truth_label}.truvari_sizebin.json"), emit: sizebin
    tuple val(meta), path("${meta.id}.${truth_label}.truvari/"),             emit: dir

    script:
    """
    # Overall benchmark
    # --pctseq 0: skip sequence similarity (symbolic alleles have no sequence)
    # --typeignore: allow DUP/INV calls to match truth INS/DEL (T2TQ100-V1.0 has no DUP/INV)
    truvari bench \\
        -b ${truth_vcf} \\
        -c ${query_vcf} \\
        --includebed ${truth_bed} \\
        -o ${meta.id}.${truth_label}.truvari \\
        --passonly \\
        --pick multi \\
        --pctseq 0 \\
        --typeignore \\
        --sizemin 50

    cp ${meta.id}.${truth_label}.truvari/summary.json ${meta.id}.${truth_label}.truvari_summary.json

    # Per-size-bin benchmarks: 50-300 bp, 300 bp-1 kb, 1-10 kb, >10 kb
    truvari bench -b ${truth_vcf} -c ${query_vcf} --includebed ${truth_bed} \\
        -o ${meta.id}.${truth_label}.truvari_50_300  --passonly --pick multi --pctseq 0 --typeignore --sizemin 50   --sizemax 300   || true
    truvari bench -b ${truth_vcf} -c ${query_vcf} --includebed ${truth_bed} \\
        -o ${meta.id}.${truth_label}.truvari_300_1k  --passonly --pick multi --pctseq 0 --typeignore --sizemin 300  --sizemax 1000  || true
    truvari bench -b ${truth_vcf} -c ${query_vcf} --includebed ${truth_bed} \\
        -o ${meta.id}.${truth_label}.truvari_1k_10k  --passonly --pick multi --pctseq 0 --typeignore --sizemin 1000 --sizemax 10000 || true
    truvari bench -b ${truth_vcf} -c ${query_vcf} --includebed ${truth_bed} \\
        -o ${meta.id}.${truth_label}.truvari_gt10k   --passonly --pick multi --pctseq 0 --typeignore --sizemin 10000               || true

    # Merge size-bin summaries into one JSON
    python3 - <<'PYEOF'
import json, os
label = "${truth_label}"
sample = "${meta.id}"
bins = {}
for name, bin_label in [
    (f"{sample}.{label}.truvari_50_300/summary.json",  "50-300 bp"),
    (f"{sample}.{label}.truvari_300_1k/summary.json",  "300 bp-1 kb"),
    (f"{sample}.{label}.truvari_1k_10k/summary.json",  "1-10 kb"),
    (f"{sample}.{label}.truvari_gt10k/summary.json",   ">10 kb"),
]:
    try:
        with open(name) as f:
            d = json.load(f)
        bins[bin_label] = {k: d.get(k, 0) for k in ("precision", "recall", "f1")}
    except Exception:
        bins[bin_label] = {"precision": 0, "recall": 0, "f1": 0}
with open(f"{sample}.{label}.truvari_sizebin.json", "w") as f:
    json.dump(bins, f)
PYEOF
    """
}
