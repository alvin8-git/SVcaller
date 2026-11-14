#!/usr/bin/env bash
# Smoke test: run the real ALPHA_GLOBIN + CNV_TRAITS subworkflows on the committed
# GIAB HG002 slice (tests/fixtures/smoke/, ~12 MB) and check the calls.
# Usage (from anywhere):  bash tests/smoke/run.sh [scratch_dir]
# Needs: Nextflow >= 25.10, Docker, the svcaller/utils:1.3 image (built locally,
# see README "Installation"). About 30 s once the images are present.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="${1:-${REPO}/work_smoke}"   # matches the work_* gitignore rule
mkdir -p "${OUT}"
cd "${REPO}"

NXF_ANSI_LOG=false nextflow run tests/smoke/main.nf -profile docker \
    --outdir "${OUT}/results" --tmp_dir "${OUT}" \
    -work-dir "${OUT}/work" > "${OUT}/smoke.log" 2>&1 \
    || { echo "FAIL: nextflow exited non-zero, see ${OUT}/smoke.log"; exit 1; }

R="${OUT}/results/HG002"
fail=0
# check <file> <awk column> <expected>  (first non-comment data row)
check() {
    local got
    got=$(grep -v '^#' "$1" | awk -F'\t' -v c="$2" 'NR==1 && $1=="sample"{next} {print $c; exit}')
    if [[ "${got}" == "$3" ]]; then echo "ok    $(basename "$1") col$2 = $3"
    else echo "FAIL  $(basename "$1") col$2: expected '$3', got '${got}'"; fail=1; fi
}

# Expected values are those of the full 30x HG002 BAM (results/HG002 and
# docs/demo/HG002_report.html); the slice must reproduce them.
check "${R}/bloodgroup/HG002.rh_status.tsv"      2 2        # RHD copies
check "${R}/bloodgroup/HG002.rh_status.tsv"      3 pos      # Rh(D)
check "${R}/cnv_traits/HG002.gst_null.tsv"       2 present  # GSTM1
check "${R}/cnv_traits/HG002.gst_null.tsv"       3 present  # GSTT1
check "${R}/cnv_traits/HG002.amy1.tsv"           2 2        # AMY1 copies
check "${R}/cnv_traits/HG002.lpa_kiv2.tsv"       2 5        # KIV-2 copies
check "${R}/alpha_globin/HG002.alpha_globin.tsv" 2 4        # alpha genes
check "${R}/alpha_globin/HG002.alpha_globin.tsv" 8 aa/aa    # genotype
check "${R}/alpha_globin/HG002.alpha_globin.tsv" 4 none     # deletion alleles

if (( fail )); then echo "SMOKE TEST FAILED"; exit 1; fi
echo "SMOKE TEST PASSED (outputs in ${R})"
