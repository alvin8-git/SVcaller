#!/usr/bin/env bash
# Run 3 SMA samples in parallel, each in its own work dir.
# Prerequisites: clean work/ directory first to free ~3.5 TB of disk space.
# Usage: bash validation/run_smn_parallel.sh
#
# Each sample gets work_smn/<SAMPLE> to isolate session locks and cache.
# Logs at /data/alvin/tmp/smn_<SAMPLE>_run1.log

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REF=/data/alvin/ref/GRCh38/hg38.fa
INTERVALS=/data/alvin/ref/GRCh38/wgs_autosomal.bed
PON=/data/alvin/SVcaller/pon/pon/giab_cnv_pon.hdf5
EH_CATALOG=${REPO_DIR}/assets/eh_catalog.json
SV_PON=/data/alvin/SVcaller/pon/sv_pon/giab_sv_pon.bed
OUTDIR=${REPO_DIR}/results_smn
LOG_DIR=/data/alvin/tmp

mkdir -p "$OUTDIR"

echo "Launching 3 SMA samples in parallel..."

for SAMPLE in SMAM SMAD SMAPB; do
  SAMPLESHEET=${REPO_DIR}/validation/smn_${SAMPLE}_samplesheet.csv
  WORK_DIR=${REPO_DIR}/work_smn/${SAMPLE}
  LOG=${LOG_DIR}/smn_${SAMPLE}_run1.log

  mkdir -p "$WORK_DIR"
  echo "  Starting ${SAMPLE} — work: ${WORK_DIR}, log: ${LOG}"

  NXF_ANSI_LOG=false nohup nextflow run "${REPO_DIR}/main.nf" \
    -profile docker \
    --input       "$SAMPLESHEET" \
    --ref_fasta   "$REF" \
    --intervals   "$INTERVALS" \
    --pon         "$PON" \
    --eh_catalog  "$EH_CATALOG" \
    --sv_pon      "$SV_PON" \
    --skip_gridss true \
    --outdir      "$OUTDIR" \
    -work-dir     "$WORK_DIR" \
    > "$LOG" 2>&1 &

  echo "    PID $!"
done

echo ""
echo "All 3 launched. Monitor with:"
echo "  tail -f /data/alvin/tmp/smn_SMAM_run1.log"
echo "  tail -f /data/alvin/tmp/smn_SMAD_run1.log"
echo "  tail -f /data/alvin/tmp/smn_SMAPB_run1.log"
