#!/usr/bin/env bash
# Slice the globin loci out of 1000G 30x CRAMs without downloading them.
#
# Each CRAM is 14-16 GB; the regions below total a few MB per sample. Slicing
# remotely over HTTPS turns a 750 GB download into well under a gigabyte.
#
#   bash validation/fetch_1000g_slices.sh [-o OUTDIR] [-j JOBS] [-n]
#
# Input : validation/hbb_1000g_carriers.tsv  (23 carriers, HbE / CD41-42)
#         validation/hbb_1000g_controls.tsv  (25 matched non-carriers)
# Output: OUTDIR/<sample>.globin.bam(.bai) + OUTDIR/manifest.tsv
#
# WHAT THESE SLICES CAN AND CANNOT DO
#
#   CAN  validate the beta site-lookup logic, and alpha-globin module logic
#        (depth + junction) as a standalone caller.
#   CANNOT run the Nextflow pipeline. MOSDEPTH's genome-wide coverage QC fails
#        params.min_depth on a sliced BAM, and SAMTOOLS_FILTER_CHROMS expects
#        whole chromosomes. Keep 2-3 whole CRAMs for one honest end-to-end run.
#
#   The CTRL_* windows below are NOT optional. Alpha copy-number calling
#   normalizes locus depth against control regions; slice chr16 alone and
#   normalization produces confident nonsense instead of failing loudly.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUTDIR="${TMPDIR:-/data/alvin/tmp}/1000g_globin_slices"
JOBS=4
DRYRUN=0

while getopts "o:j:nh" opt; do
  case $opt in
    o) OUTDIR="$OPTARG" ;;
    j) JOBS="$OPTARG" ;;
    n) DRYRUN=1 ;;
    h) sed -n '2,25p' "$0"; exit 0 ;;
    *) exit 2 ;;
  esac
done

# CRAM cannot be decoded without the exact reference it was encoded against.
# The local hg38.fa MD5s do NOT match the 1000G analysis set, so point htslib
# at the ENA registry instead of a local FASTA.
export REF_PATH="${REF_PATH:-https://www.ebi.ac.uk/ena/cram/md5/%s}"

# chr11 HBB cluster + LCR; chr16 alpha cluster; 8 CTRL windows for depth
# normalization (must match assets/cnv_trait_regions.bed).
REGIONS=(
  chr11:5200000-5300000            # HBB cluster + locus control region
  chr16:1-500000                   # HBA1/HBA2/HBZ cluster
  chr2:100000000-100020000         # CTRL_1
  chr3:100000000-100020000         # CTRL_2
  chr4:100000000-100020000         # CTRL_3
  chr5:100000000-100020000         # CTRL_4
  chr7:100000000-100020000         # CTRL_5
  chr8:100000000-100020000         # CTRL_6
  chr11:100000000-100020000        # CTRL_7
  chr12:100000000-100020000        # CTRL_8
)

command -v samtools >/dev/null || { echo "FATAL: samtools not on PATH" >&2; exit 1; }
SAMTOOLS_VER="$(samtools --version | head -1)"
STAMP="$(date -u +%Y-%m-%dT%H:%M:%SZ)"

mkdir -p "$OUTDIR"
MANIFEST="$OUTDIR/manifest.tsv"
[ -s "$MANIFEST" ] || printf 'sample\tclass\tallele\tpop\treads\tbytes\tfetched_utc\tsource_cram\n' > "$MANIFEST"

# ---- build the work list from both TSVs -------------------------------------
# fields: sample, allele, genotype, pop, cram_url  (see the TSV headers)
worklist() {
  awk -F'\t' -v cls="$2" '!/^#/ && $1!="sample" && NF>=10 {
      print $1"\t"cls"\t"$2"\t"$9"\t"$10
  }' "$1"
}
WORK="$OUTDIR/.worklist"
: > "$WORK"
worklist "$REPO/validation/hbb_1000g_carriers.tsv" carrier >> "$WORK"
worklist "$REPO/validation/hbb_1000g_controls.tsv" control >> "$WORK"

# A sample can appear twice in the carrier file (HG02379 carries both alleles).
# Collapse to one slice per sample, joining the allele labels.
sort -t$'\t' -k1,1 "$WORK" | awk -F'\t' '
  { if ($1==prev) { al=al"+"$3 }
    else { if (prev!="") print prev"\t"cls"\t"al"\t"pop"\t"url;
           prev=$1; cls=$2; al=$3; pop=$4; url=$5 } }
  END { if (prev!="") print prev"\t"cls"\t"al"\t"pop"\t"url }
' > "$WORK.uniq"
mv "$WORK.uniq" "$WORK"

TOTAL=$(wc -l < "$WORK")
echo "samples: $TOTAL   regions: ${#REGIONS[@]}   jobs: $JOBS"
echo "outdir : $OUTDIR"
[ "$DRYRUN" = 1 ] && { echo "--- dry run ---"; cat "$WORK"; exit 0; }

# ---- slice one sample -------------------------------------------------------
slice_one() {
  local sample="$1" cls="$2" allele="$3" pop="$4" url="$5"
  local out="$OUTDIR/${sample}.globin.bam"
  local raw="$out.raw" hdr="$out.hdr"
  # REGIONS is an array in the parent; arrays do not survive into this subshell,
  # so rebuild it from the exported string.
  local -a regions
  read -ra regions <<< "$REGIONS_STR"

  if [ -s "$out" ] && samtools quickcheck "$out" 2>/dev/null; then
    echo "  skip  $sample (already present)"; return 0
  fi
  if [ "$url" = "NOT_FOUND" ] || [ -z "$url" ]; then
    echo "  FAIL  $sample - no CRAM url" >&2; return 1
  fi

  if ! samtools view -b "$url" "${regions[@]}" > "$raw"; then
    echo "  FAIL  $sample - slice failed" >&2; rm -f "$raw"; return 1
  fi
  if ! samtools quickcheck "$raw" 2>/dev/null; then
    echo "  FAIL  $sample - truncated slice" >&2; rm -f "$raw"; return 1
  fi

  # Provenance: without this the slices become mystery BAMs within months.
  samtools view -H "$raw" > "$hdr"
  {
    printf '@CO\tslice_source:%s\n' "$url"
    printf '@CO\tslice_regions:%s\n' "$REGIONS_STR"
    printf '@CO\tslice_fetched_utc:%s\n' "$STAMP"
    printf '@CO\tslice_tool:%s\n' "$SAMTOOLS_VER"
    printf '@CO\tslice_note:%s\n' "PARTIAL GENOME - genome-wide QC (mosdepth min_depth) will fail; not pipeline-runnable"
  } >> "$hdr"

  samtools reheader "$hdr" "$raw" > "$out"
  rm -f "$raw" "$hdr"
  samtools index "$out"

  local reads bytes
  reads=$(samtools view -c "$out")
  bytes=$(stat -c%s "$out")
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$sample" "$cls" "$allele" "$pop" "$reads" "$bytes" "$STAMP" "$url" >> "$MANIFEST"
  echo "  ok    $sample  ${cls}/${allele}  ${reads} reads  $((bytes/1024)) KB"
}
export -f slice_one
export OUTDIR MANIFEST STAMP SAMTOOLS_VER REF_PATH
export REGIONS_STR="${REGIONS[*]}"

# xargs -P gives parallelism without a GNU parallel dependency. Each job is
# network-bound, so oversubscribing cores is fine; be polite to EBI though.
FAILED=0
if ! awk -F'\t' '{printf "%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5}' "$WORK" \
   | xargs -P "$JOBS" -I{} bash -c 'IFS=$'"'"'\t'"'"' read -r s c a p u <<< "{}"; slice_one "$s" "$c" "$a" "$p" "$u"'; then
  FAILED=1
fi

DONE=$(( $(wc -l < "$MANIFEST") - 1 ))
echo
echo "slices present: $DONE / $TOTAL"
echo "manifest      : $MANIFEST"
[ "$DONE" -lt "$TOTAL" ] && echo "NOTE: $((TOTAL-DONE)) sample(s) missing - re-run to retry (completed slices are skipped)."
exit 0
