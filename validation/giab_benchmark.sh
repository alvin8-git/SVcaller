#!/usr/bin/env bash
# Run truvari benchmarking for one GIAB sample.
# Usage: giab_benchmark.sh <sample_id> <query_vcf> [giab_truth_dir]
# Example: giab_benchmark.sh HG002 /data/alvin/tmp/svcaller_HG002/results/HG002/HG002.sv_merged.vcf.gz
set -euo pipefail

SAMPLE="${1:?Usage: $0 SAMPLE_ID QUERY_VCF [GIAB_TRUTH_DIR]}"
QUERY_VCF="${2:?QUERY_VCF required}"
TRUTH_DIR="${3:-/data/alvin/ref/GIAB}"
OUT_DIR="/data/alvin/tmp/truvari_${SAMPLE}"

TRUTH_VCF="${TRUTH_DIR}/${SAMPLE}_SV_v0.6.vcf.gz"
TRUTH_TBI="${TRUTH_DIR}/${SAMPLE}_SV_v0.6.vcf.gz.tbi"
TRUTH_BED="${TRUTH_DIR}/${SAMPLE}_SV_v0.6_HC.bed"

if [ ! -f "${TRUTH_VCF}" ]; then
    echo "ERROR: Truth VCF not found: ${TRUTH_VCF}" >&2
    echo "  Run validation/download_refs.sh first." >&2
    exit 1
fi

if [ ! -f "${QUERY_VCF}" ]; then
    echo "ERROR: Query VCF not found: ${QUERY_VCF}" >&2
    exit 1
fi

mkdir -p "${OUT_DIR}"

QUERY_DIR="$(dirname "$(realpath "${QUERY_VCF}")")"

docker run --rm \
  -v "${QUERY_DIR}:/query:ro" \
  -v "${TRUTH_DIR}:/truth:ro" \
  -v "${OUT_DIR}:/out" \
  quay.io/biocontainers/truvari:4.2.2--pyhdfd78af_0 \
  truvari bench \
    -b "/truth/$(basename "${TRUTH_VCF}")" \
    -c "/query/$(basename "${QUERY_VCF}")" \
    --includebed "/truth/$(basename "${TRUTH_BED}")" \
    -o /out \
    --passonly \
    --pick multi \
    --sizemin 50

echo ""
echo "=== ${SAMPLE} Benchmark Results ==="
if [ -f "${OUT_DIR}/summary.json" ]; then
    python3 -c "
import json, sys
with open('${OUT_DIR}/summary.json') as f:
    d = json.load(f)
print(f'  Precision: {d[\"precision\"]:.4f}')
print(f'  Recall:    {d[\"recall\"]:.4f}')
print(f'  F1:        {d[\"f1\"]:.4f}')
print(f'  TP base:   {d[\"TP-base\"]}  FP: {d[\"FP\"]}  FN: {d[\"FN\"]}')
print(f'  Full results: ${OUT_DIR}/summary.json')
"
fi
