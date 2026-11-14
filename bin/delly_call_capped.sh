#!/bin/sh
# Run `delly call` for one SV type; for INV only, bound it by a wall-clock cap.
#
# usage: delly_call_capped.sh <sample> <svtype> <limit_seconds> <fasta> <fai> <bam>
#
# Why: on some samples `delly call -t INV` sits in "Split-read assembly" for 7-8 h (one thread,
# no I/O) while the other types finish in ~80 min. When the cap fires the sample CONTINUES with a
# header-only <sample>.delly.INV.vcf (valid input for DELLY_MERGE) and the explicit record
# <sample>.sv_not_assessed.tsv says delly inversions were NOT ASSESSED. It never reads as
# "no inversions found". ERROR is never NEGATIVE.
#
#   limit_seconds <= 0   no cap (INV behaves like the other types)
#   any other type       never capped, never writes the record file
#   INV                  always writes <sample>.sv_not_assessed.tsv: header only = delly INV
#                        finished (or was uncapped); one row = delly INV hit the cap
#   delly fails before the cap  the script exits with delly's code (a real failure fails the task)
# POSIX sh + awk/zcat/strings/bgzip: runs inside the delly container.
# Tested in tests/test_delly_inv_cap.py.
set -eu
sample=$1; svtype=$2; limit=$3; fasta=$4; fai=$5; bam=$6
out="${sample}.delly.${svtype}.vcf"
rec="${sample}.sv_not_assessed.tsv"

if [ "$svtype" != "INV" ]; then
    exec delly call -t "$svtype" -g "$fasta" -o "$out" "$bam"
fi

printf 'sample\tcaller\tsvtype\treason\tlimit\n' > "$rec"
# Own watchdog instead of coreutils `timeout`: the delly image ships BusyBox timeout, which
# reports 128+signal (143) and never 124, so its exit code cannot tell "cap fired" from any
# other SIGTERM. The sentinel file is created only by the watchdog, just before it kills delly.
rm -f inv_cap_fired
rc=0
if [ "$limit" -gt 0 ]; then
    delly call -t INV -g "$fasta" -o "$out" "$bam" &
    dpid=$!
    # Poll in 1 s steps (no long-lived child to orphan); ends as soon as delly is gone.
    ( i=0
      while [ "$i" -lt "$limit" ] && kill -0 "$dpid" 2>/dev/null; do sleep 1; i=$((i+1)); done
      if kill -0 "$dpid" 2>/dev/null; then
          : > inv_cap_fired                      # sentinel first, before any signal
          kill -TERM "$dpid" 2>/dev/null || true
          g=0
          while [ "$g" -lt 10 ] && kill -0 "$dpid" 2>/dev/null; do sleep 1; g=$((g+1)); done
          kill -KILL "$dpid" 2>/dev/null || true  # TERM ignored or slow: escalate after ~10 s
      fi ) &
    wd=$!
    wait "$dpid" || rc=$?
    # Do not kill the watchdog: if it is mid-escalation it must finish. It exits on its own
    # within ~10 s once delly is gone.
    wait "$wd" 2>/dev/null || true
else
    delly call -t INV -g "$fasta" -o "$out" "$bam" || rc=$?
fi
[ "$rc" -eq 0 ] && exit 0
# Capped only if the watchdog fired AND delly died of its TERM (143) or KILL (137); any other status
# is a real delly failure and fails the task.
if [ ! -f inv_cap_fired ] || { [ "$rc" -ne 143 ] && [ "$rc" -ne 137 ]; }; then exit "$rc"; fi
rm -f inv_cap_fired

echo "WARN: delly call -t INV exceeded the ${limit}s cap for ${sample}; INVERSIONS WERE NOT ASSESSED by delly for this sample (this is NOT 'no inversions found'). Continuing without delly INV calls." >&2
rm -f "$out" "$out.csi"
# VCF sample name exactly as delly 1.2.6 writes it (checked against real delly VCFs): the SM of the
# first @RG header line; with no @RG SM, the BAM file name minus its extension. The BAM header text is
# BGZF (zcat); only real @RG lines are parsed, so an @PG command line mentioning SM: cannot match.
sm=$(zcat "$bam" 2>/dev/null | head -c 16000000 | strings -n 1 | awk -F'\t' '
    { for (i = 1; i <= NF; i++) if ($i == "@RG") { for (j = i + 1; j <= NF; j++) if ($j ~ /^SM:/) { print substr($j, 4); exit } ; exit } }' | head -1 || true)
if [ -z "$sm" ]; then
    sm=$(basename "$bam"); sm=${sm%.*}
fi
{
    printf '##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description="All filters passed">\n'
    printf '##ALT=<ID=INV,Description="Inversion">\n'
    printf '##FILTER=<ID=LowQual,Description="Poor quality and insufficient number of PEs and SRs.">\n'
    printf '##INFO=<ID=CIEND,Number=2,Type=Integer,Description="PE confidence interval around END">\n'
    printf '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="PE confidence interval around POS">\n'
    printf '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">\n'
    printf '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n'
    printf '##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">\n'
    printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    awk '{printf "##contig=<ID=%s,length=%s>\n", $1, $2}' "$fai"
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$sm"
} | bgzip -c > "$out"
printf '%s\tdelly\tINV\ttime_limit_exceeded\t%ss\n' "$sample" "$limit" >> "$rec"
