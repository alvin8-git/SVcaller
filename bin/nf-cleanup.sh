#!/usr/bin/env bash
# Post-run cleanup. Two modes:
#
#   bash bin/nf-cleanup.sh <sample_id>
#       Remove work/<sample_id> after confirming results are published, and
#       prune orphaned .nextflow/cache sessions.  (per-sample full delete)
#
#   bash bin/nf-cleanup.sh --reclaim [--force]
#       Reclaim disk from orphaned work dirs left by failed/superseded runs that
#       share a -work-dir (e.g. crashed attempts), keeping only the latest
#       SUCCESSFUL run so -resume still works. Dry-run by default — pass --force
#       to actually delete. Refuses to run while a pipeline is still active.
#
# Safety: per-sample mode checks results exist first; reclaim mode is dry-run
# unless --force, and aborts if any Nextflow session is currently locked.

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Return 0 if any Nextflow cache session is currently held (a run is active).
nf_session_active() {
  local cache="${REPO_DIR}/.nextflow/cache" db lock
  [ -d "${cache}" ] || return 1
  for db in "${cache}"/*/db; do
    lock="${db}/LOCK"
    [ -f "${lock}" ] || continue
    if lsof "${lock}" >/dev/null 2>&1; then return 0; fi
  done
  return 1
}

# ---------------------------------------------------------------------------
# Mode: --reclaim  (surgical orphan reclamation via `nextflow clean`)
# ---------------------------------------------------------------------------
if [ "${1:-}" = "--reclaim" ]; then
  FORCE=0
  [ "${2:-}" = "--force" ] && FORCE=1

  command -v nextflow >/dev/null 2>&1 || { echo "ERROR: nextflow not on PATH." >&2; exit 1; }
  cd "${REPO_DIR}"

  if nf_session_active; then
    echo "ERROR: a Nextflow run is still active (cache LOCK held)." >&2
    echo "       Wait for it to finish — cleaning now could delete live task dirs." >&2
    exit 1
  fi

  # Latest SUCCESSFUL run name = the run name immediately before the first 'OK'
  # status field on the most recent OK line of `nextflow log` (anchoring on the
  # status token is robust to the multi-token timestamp/duration columns).
  KEEP="$(nextflow log 2>/dev/null \
    | awk 'NR>1{for(i=1;i<=NF;i++){if($i=="OK"){print $(i-1); break} else if($i=="ERR") break}}' \
    | tail -1)"

  if [ -z "${KEEP}" ]; then
    echo "ERROR: no successful run found in 'nextflow log'. Nothing reclaimed." >&2
    exit 1
  fi

  echo "Latest successful run to KEEP: ${KEEP}"
  echo "WARNING: 'nextflow clean -but ${KEEP}' removes the work dirs of ALL OTHER"
  echo "         runs in this repo's history (any sample). Review the list below."
  echo ""

  if [ "${FORCE}" -eq 1 ]; then
    echo "Deleting (forced):"
    nextflow clean -but "${KEEP}" -f
    echo "Done. Reclaimed orphaned work dirs; kept run '${KEEP}' for -resume."
  else
    echo "[dry-run] would remove (pass --force to delete):"
    nextflow clean -but "${KEEP}" -n
    echo ""
    echo "Re-run with --force to delete:  bash bin/nf-cleanup.sh --reclaim --force"
  fi
  exit 0
fi

# ---------------------------------------------------------------------------
# Mode: per-sample full delete
# ---------------------------------------------------------------------------
SAMPLE="${1:?Usage: nf-cleanup.sh <sample_id> | --reclaim [--force]}"
WORK_DIR="${REPO_DIR}/work/${SAMPLE}"

# 1. Verify results published before deleting work dir
RESULTS_FOUND=0
for OUTDIR in "${REPO_DIR}/results"; do
  if [ -d "${OUTDIR}/${SAMPLE}" ] || ls "${OUTDIR}"/*.html 2>/dev/null | grep -q "${SAMPLE}" 2>/dev/null; then
    RESULTS_FOUND=1
    echo "Results confirmed at: ${OUTDIR}"
    break
  fi
done

if [ "${RESULTS_FOUND}" -eq 0 ]; then
  echo "ERROR: No results found for ${SAMPLE} in ${REPO_DIR}/results/." >&2
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
