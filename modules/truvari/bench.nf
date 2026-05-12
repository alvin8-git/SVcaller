process TRUVARI_BENCH {
    tag "${meta.id}"
    label 'process_low'
    container 'quay.io/biocontainers/truvari:4.2.2--pyhdfd78af_0'

    input:
    tuple val(meta), path(query_vcf), path(query_tbi)
    path truth_vcf
    path truth_tbi
    path truth_bed

    output:
    tuple val(meta), path("${meta.id}.truvari/summary.json"),  emit: summary
    tuple val(meta), path("${meta.id}.truvari_sizebin.json"),  emit: sizebin
    tuple val(meta), path("${meta.id}.truvari/"),              emit: dir

    script:
    """
    # Overall benchmark
    truvari bench \\
        -b ${truth_vcf} \\
        -c ${query_vcf} \\
        --includebed ${truth_bed} \\
        -o ${meta.id}.truvari \\
        --passonly \\
        --pick multi \\
        --sizemin 50

    # Per-size-bin benchmarks: 50-300 bp, 300 bp-1 kb, 1-10 kb, >10 kb
    truvari bench -b ${truth_vcf} -c ${query_vcf} --includebed ${truth_bed} \\
        -o ${meta.id}.truvari_50_300  --passonly --pick multi --sizemin 50   --sizemax 300   || true
    truvari bench -b ${truth_vcf} -c ${query_vcf} --includebed ${truth_bed} \\
        -o ${meta.id}.truvari_300_1k  --passonly --pick multi --sizemin 300  --sizemax 1000  || true
    truvari bench -b ${truth_vcf} -c ${query_vcf} --includebed ${truth_bed} \\
        -o ${meta.id}.truvari_1k_10k  --passonly --pick multi --sizemin 1000 --sizemax 10000 || true
    truvari bench -b ${truth_vcf} -c ${query_vcf} --includebed ${truth_bed} \\
        -o ${meta.id}.truvari_gt10k   --passonly --pick multi --sizemin 10000               || true

    # Merge size-bin summaries into one JSON
    python3 - <<'PYEOF'
import json, os
bins = {}
for name, label in [
    ("${meta.id}.truvari_50_300/summary.json",  "50-300 bp"),
    ("${meta.id}.truvari_300_1k/summary.json",  "300 bp-1 kb"),
    ("${meta.id}.truvari_1k_10k/summary.json",  "1-10 kb"),
    ("${meta.id}.truvari_gt10k/summary.json",   ">10 kb"),
]:
    try:
        with open(name) as f:
            d = json.load(f)
        bins[label] = {k: d.get(k, 0) for k in ("precision", "recall", "f1")}
    except Exception:
        bins[label] = {"precision": 0, "recall": 0, "f1": 0}
with open("${meta.id}.truvari_sizebin.json", "w") as f:
    json.dump(bins, f)
PYEOF
    """
}
