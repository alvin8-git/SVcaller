#!/usr/bin/env bash
# Post-run cleanup: remove work dir and orphaned Nextflow cache sessions.
# Run after confirming results have been published to --outdir.
#
# Usage: bash bin/nf-cleanup.sh <sample_id>
#   e.g. bash bin/nf-cleanup.sh SMAM
#        bash bin/nf-cleanup.sh HG002
#
# Safety: checks results exist before deleting work dir.

set -euo pipefail

SAMPLE="${1:?Usage: nf-cleanup.sh <sample_id>}"
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK_DIR="${REPO_DIR}/work_${SAMPLE}"

# 1. Verify results published before deleting work dir
RESULTS_FOUND=0
for OUTDIR in "${REPO_DIR}/results" "${REPO_DIR}/results_smn" "${REPO_DIR}/results_giab" "${REPO_DIR}/results_${SAMPLE}"; do
  if [ -d "${OUTDIR}/${SAMPLE}" ] || ls "${OUTDIR}"/*.html 2>/dev/null | grep -q "${SAMPLE}" 2>/dev/null; then
    RESULTS_FOUND=1
    echo "Results confirmed at: ${OUTDIR}"
    break
  fi
done

if [ "${RESULTS_FOUND}" -eq 0 ]; then
  echo "ERROR: No results found for ${SAMPLE} in any results_* directory." >&2
  echo "       Run 'ls results*/*/  to verify outputs before cleaning." >&2
  exit 1
fi

# 2. Remove work dir
if [ -d "${WORK_DIR}" ]; then
  WORK_SIZE=$(du -sh "${WORK_DIR}" 2>/dev/null | cut -f1)
  echo "Removing work dir: ${WORK_DIR} (${WORK_SIZE})"
  rm -rf "${WORK_DIR}"
  echo "Done."
else
  echo "Work dir not found: ${WORK_DIR} (already clean)"
fi

# 3. Prune orphaned .nextflow/cache session dirs (those not locked by any process)
NXF_CACHE="${REPO_DIR}/.nextflow/cache"
if [ -d "${NXF_CACHE}" ]; then
  PRUNED=0
  for DB_DIR in "${NXF_CACHE}"/*/db; do
    LOCK="${DB_DIR}/LOCK"
    [ -f "${LOCK}" ] || continue
    if ! lsof "${LOCK}" >/dev/null 2>&1; then
      SESSION_DIR="$(dirname "${DB_DIR}")"
      echo "Pruning orphan cache: ${SESSION_DIR}"
      rm -rf "${SESSION_DIR}"
      PRUNED=$((PRUNED+1))
    fi
  done
  [ "${PRUNED}" -gt 0 ] && echo "Pruned ${PRUNED} orphan cache session(s)." || echo "No orphan cache sessions found."
fi

# 4. Optional: prune dangling Docker layers
echo ""
echo "Optional: run 'docker image prune -f' to reclaim dangling Docker layers."
