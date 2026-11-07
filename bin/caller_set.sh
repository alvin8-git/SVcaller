#!/bin/sh
# Decide which callers enter the Jasmine merge, in SUPP_VEC bit order, and record it.
#
# usage: caller_set.sh <sample> NAME=FILE:MODE [NAME=FILE:MODE ...]   (order = bit order)
#   MODE  always   caller ran; always gets a SUPP_VEC bit (Manta, Delly, GRIDSS when run)
#         ifcalls  caller ran; gets a bit only if its VCF has >=1 record (Scramble, MELT: legacy rule)
#         off      caller did not run; no bit, listed under callers_not_run
# writes (cwd): vcf_list.txt  callers.txt (one name per line = bit order)  <sample>.sv_callers.json
# POSIX sh + grep only: runs inside the Jasmine container. Tested in tests/test_caller_set.py.
set -eu
sample=$1; shift
: > vcf_list.txt; : > callers.txt
ran=""; notrun=""; nocalls=""
for spec in "$@"; do
    name=${spec%%=*}; rest=${spec#*=}; mode=${rest##*:}; file=${rest%:*}
    case "$mode" in
        off)     notrun="$notrun $name" ;;
        always)  echo "$file" >> vcf_list.txt; echo "$name" >> callers.txt; ran="$ran $name" ;;
        ifcalls) n=$(grep -cv '^#' "$file" 2>/dev/null || true)
                 ran="$ran $name"
                 if [ "${n:-0}" -gt 0 ]; then echo "$file" >> vcf_list.txt; echo "$name" >> callers.txt
                 else nocalls="$nocalls $name"; fi ;;
        *) echo "caller_set.sh: bad mode '$mode' in '$spec'" >&2; exit 2 ;;
    esac
done
jl() { printf '['; sep=""; for x in $1; do printf '%s"%s"' "$sep" "$x"; sep=","; done; printf ']'; }
width=$(grep -c . callers.txt || true)
printf '{"sample":"%s","supp_vec_width":%s,"callers_in_supp_vec_order":%s,"callers_ran":%s,"callers_ran_no_calls":%s,"callers_not_run":%s}\n' \
    "$sample" "$width" "$(jl "$(cat callers.txt | tr '\n' ' ')")" "$(jl "$ran")" "$(jl "$nocalls")" "$(jl "$notrun")" \
    > "${sample}.sv_callers.json"
