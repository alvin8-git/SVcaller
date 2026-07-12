#!/usr/bin/env bash
# Run GIAB report pipeline for each sample in its own work directory.
# Samples run sequentially to avoid CPU/memory contention.
# Usage: bash validation/run_giab_reports.sh [HG001 HG003 ...]
#   (defaults to all 6 samples if no args given)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
REF=/data/alvin/ref/GRCh38/hg38.canonical.fa
PON=/data/alvin/SVcaller/pon/pon/giab_cnv_pon.hdf5
EH_CATALOG=${REPO_DIR}/assets/eh_catalog.json
ANNOTSV_DB=/data/alvin/ref/annotsv/Annotations_Human
SV_PON=/data/alvin/SVcaller/pon/sv_pon/giab_sv_pon.bed
INTERVALS=/data/alvin/ref/GRCh38/wgs_autosomal.bed
OUTDIR=${REPO_DIR}/results
WORKDIR_BASE=${REPO_DIR}/work
LOG_DIR=/data/alvin/tmp

SAMPLES=("${@:-HG001 HG003 HG004 HG005 HG006 HG007}")
if [ $# -eq 0 ]; then
  SAMPLES=(HG001 HG003 HG004 HG005 HG006 HG007)
fi

mkdir -p "$WORKDIR_BASE" "$OUTDIR"

for SAMPLE in "${SAMPLES[@]}"; do
  SAMPLESHEET=${SCRIPT_DIR}/giab_${SAMPLE}_samplesheet.csv
  if [ ! -f "$SAMPLESHEET" ]; then
    echo "ERROR: samplesheet not found: $SAMPLESHEET" >&2
    exit 1
  fi

  WORK_DIR=${WORKDIR_BASE}/${SAMPLE}
  LOG=${LOG_DIR}/giab_${SAMPLE}_run1.log

  echo ">>> Starting ${SAMPLE} — work dir: ${WORK_DIR}"
  echo "    Log: ${LOG}"

  NXF_ANSI_LOG=false nextflow run "${REPO_DIR}/main.nf" \
    -profile docker \
    --input       "$SAMPLESHEET" \
    --ref_fasta   "$REF" \
    --intervals   "$INTERVALS" \
    --pon         "$PON" \
    --eh_catalog  "$EH_CATALOG" \
    --annotsv_db  "$ANNOTSV_DB" \
    --sv_pon      "$SV_PON" \
    --skip_gridss true \
    --skip_melt   true \
    --outdir      "$OUTDIR" \
    -work-dir     "$WORK_DIR" \
    -resume \
    > "$LOG" 2>&1

  EXIT=$?
  if [ $EXIT -ne 0 ]; then
    echo "ERROR: ${SAMPLE} failed (exit ${EXIT}) — check ${LOG}"
    exit $EXIT
  fi
  echo "    ${SAMPLE} complete."
done

echo "All GIAB samples complete. Results in: ${OUTDIR}"
