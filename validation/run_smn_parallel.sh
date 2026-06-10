#!/usr/bin/env bash
# Run SMA trio samples in parallel, each with an isolated Nextflow cache.
#
# Why isolated caches (-cache per sample):
#   Nextflow uses a single RocksDB lock file per cache directory. Without
#   isolation, concurrent runs compete for the same lock and all but one fail
#   immediately. Giving each sample its own -cache dir removes the contention
#   entirely — runs are genuinely parallel with no serialization.
#
# Why hg38.canonical.fa (not hg38.fa):
#   FILTER_CHROMS strips alt contigs from BAM reads and @SQ headers, leaving
#   only chr1-22+X+Y+M. If the reference still declares alt chromosomes (as
#   hg38.fa does), Manta fails: "BAM is missing a chromosome found in the
#   reference." The canonical reference matches the filtered BAM exactly.
#
# Usage: bash validation/run_smn_parallel.sh [RUN_LABEL]
#   RUN_LABEL defaults to "run1" — change to "run2", "run3", etc. on reruns.

set -uo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUN_LABEL="${1:-run1}"
REF=/data/alvin/ref/GRCh38/hg38.canonical.fa
INTERVALS=/data/alvin/ref/GRCh38/wgs_autosomal.bed
PON=/data/alvin/SVcaller/pon/pon/giab_cnv_pon.hdf5
EH_CATALOG=${REPO_DIR}/assets/eh_catalog.json
SV_PON=/data/alvin/SVcaller/pon/sv_pon/giab_sv_pon.bed
OUTDIR=${REPO_DIR}/results_smn
LOG_DIR=/data/alvin/tmp

mkdir -p "$OUTDIR"

echo "Launching SMA trio in parallel (run label: ${RUN_LABEL})..."
PIDS=()

for SAMPLE in SMAM SMAD SMAPB; do
  SAMPLESHEET=${REPO_DIR}/validation/smn_${SAMPLE}_samplesheet.csv
  WORK_DIR=${REPO_DIR}/work_smn/${SAMPLE}
  CACHE_DIR=${WORK_DIR}/.nxf_cache
  LOG=${LOG_DIR}/smn_${SAMPLE}_${RUN_LABEL}.log

  mkdir -p "$WORK_DIR" "$CACHE_DIR"
  echo "  Starting ${SAMPLE} — work: ${WORK_DIR}, cache: ${CACHE_DIR}, log: ${LOG}"

  NXF_ANSI_LOG=false nohup nextflow run "${REPO_DIR}/main.nf" \
    -profile docker \
    -cache        "$CACHE_DIR" \
    --input       "$SAMPLESHEET" \
    --ref_fasta   "$REF" \
    --intervals   "$INTERVALS" \
    --pon         "$PON" \
    --eh_catalog  "$EH_CATALOG" \
    --sv_pon      "$SV_PON" \
    --skip_gridss true \
    --outdir      "$OUTDIR" \
    -work-dir     "$WORK_DIR" \
    -resume \
    > "$LOG" 2>&1 &

  PIDS+=($!)
  echo "    PID $!"
done

echo ""
echo "All 3 launched. Monitor with:"
for SAMPLE in SMAM SMAD SMAPB; do
  echo "  tail -f ${LOG_DIR}/smn_${SAMPLE}_${RUN_LABEL}.log"
done
echo ""
echo "Waiting for all 3 to complete..."
FAIL=0
for PID in "${PIDS[@]}"; do
  wait "$PID" || FAIL=$((FAIL+1))
done
echo "Done. Failures: ${FAIL}/3"
