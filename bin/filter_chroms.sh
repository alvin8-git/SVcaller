#!/usr/bin/env bash
# Filter a BAM to the canonical chromosomes (chr1-22, X, Y, M), in PARALLEL by
# chromosome. Replaces the single-threaded awk that pinned one core for ~70 min
# on a 79 GB BAM (measured 2026-07-22: 19 MB/s, awk at 99.7% of one core while
# the samtools reader/writer sat idle). Fanning the per-read work over the 24
# canonical chromosomes lifts the ceiling to roughly chr1's share of the genome
# (~8%), so ~8-10x, ~7-9 min.
#
# The per-read transforms are IDENTICAL to the previous inline awk, so the output
# is content-equivalent (BGZF block boundaries differ, so the bytes do not; that
# is expected and harmless). The four transforms:
#   1. header: keep @HD, emit only canonical @SQ in FAI order, keep other @ lines.
#   2. strip non-canonical entries from SA:Z supplementary-alignment tags.
#   3. drop reads whose mate (RNEXT) is on a non-canonical contig (not "=" / "*").
#      Manta crashes FATAL_ERROR on mate_tid=-1; clearing RNEXT to "*" also crashes,
#      so the read must be dropped entirely.
#   4. when RNEXT=="*" and PNEXT!=0, zero PNEXT.
#
# Region queries (samtools view <bam> <chr>) return reads coordinate-sorted and
# exclude unmapped reads (RNAME=*), exactly as the old `view -h <bam> <regions>`
# did. Iterating chromosomes in FAI order and `samtools cat`-ing the parts (which
# concatenates BGZF blocks without recompressing) yields a canonical-sorted BAM by
# construction, so the old conditional re-sort is gone.
#
# Usage: filter_chroms.sh <in.bam> <fai> <threads> <out.bam>
set -euo pipefail

BAM="${1:?usage: filter_chroms.sh <in.bam> <fai> <threads> <out.bam>}"
FAI="${2:?fai}"
THREADS="${3:?threads}"
OUT="${4:?out.bam}"

CANON="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM"

# Canonical chromosomes present in the FAI, in FAI order. This is both the emit
# order and the set of region queries. A canonical chrom absent from the FAI is
# simply skipped (matches the old code, which keyed @SQ off FAI order).
mapfile -t CHROMS < <(awk '$1 ~ "^chr([0-9]+|X|Y|M)$" {print $1}' "$FAI")
if [ "${#CHROMS[@]}" -eq 0 ]; then
    echo "ERROR: no canonical chromosomes found in FAI '$FAI'." >&2
    exit 1
fi

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# 1. Corrected header, built once. FAI pass sets canonical order; header pass
#    keeps @HD inline, buffers canonical @SQ into FAI order, keeps other @ lines.
samtools view -H "$BAM" | awk -v canon="$CANON" '
    BEGIN { n=split(canon,c," "); for(i=1;i<=n;i++) can[c[i]]=1 }
    FNR==NR { if($1 in can) order[$1]=++ord; next }        # FAI pass
    /^@HD/  { print; next }
    /^@SQ/  {
        sn=""
        for(i=1;i<=NF;i++) if($i ~ /^SN:/){ sn=$i; sub(/^SN:/,"",sn) }
        if(sn in can && sn in order) sq[order[sn]]=$0
        next
    }
    /^@/    { other[++on]=$0; next }
    END {
        for(i=1;i<=ord;i++) if(i in sq) print sq[i]
        for(i=1;i<=on;i++) print other[i]
    }
' "$FAI" - > "$WORK/header.sam"

# 2. Per-chromosome part: reads for the chrom, transformed, prefixed with the
#    corrected header so `samtools cat` sees identical headers, compressed.
export BAM WORK CANON
filter_one() {
    local chr="$1"
    { cat "$WORK/header.sam"
      samtools view "$BAM" "$chr" | awk -v canon="$CANON" '
        BEGIN { OFS="\t"; n=split(canon,c," "); for(i=1;i<=n;i++) can[c[i]]=1 }
        {
            for (fi=12; fi<=NF; fi++) {
                if ($fi ~ /^SA:Z:/) {
                    n_sa=split(substr($fi,6),saparts,";"); new_sa=""
                    for (sj=1; sj<=n_sa; sj++) {
                        if (saparts[sj]=="") continue
                        split(saparts[sj],saf,",")
                        if (saf[1] in can) new_sa=new_sa saf[1]","saf[2]","saf[3]","saf[4]","saf[5]","saf[6]";"
                    }
                    if (new_sa!="") $fi="SA:Z:"new_sa
                    else { for (fk=fi; fk<NF; fk++) $fk=$(fk+1); NF--; fi-- }
                    break
                }
            }
            if ($7 != "=" && $7 != "*" && !($7 in can)) next
            if ($7 == "*" && $8+0 != 0) $8=0
            print
        }'
    } | samtools view -b -o "$WORK/part.$chr.bam" -
}
export -f filter_one

printf '%s\n' "${CHROMS[@]}" | xargs -P"$THREADS" -I{} bash -c 'filter_one "$1"' _ {}

# 3. Concatenate the parts in canonical order (BGZF concat, no recompress), index.
PARTS=()
for chr in "${CHROMS[@]}"; do PARTS+=("$WORK/part.$chr.bam"); done
samtools cat -o "$OUT" "${PARTS[@]}"
samtools index "$OUT"
