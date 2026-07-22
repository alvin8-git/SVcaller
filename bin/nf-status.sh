#!/usr/bin/env bash
# One-shot status of a running Nextflow pipeline.
#
# WHY THIS EXISTS
# ---------------
# Every ad-hoc way of checking a run has a failure mode that reports something
# FALSE rather than nothing, which is worse — a wrong status looks like a finding
# and gets acted on. All of these bit this project on 2026-07-22:
#
#   pgrep -f "<pattern>"      matches the CHECKING command's own command line.
#                             Reported "docker build still running" after it had
#                             finished, and "SvABA is running" when it was not.
#   pkill -f "<script>"       same, and it kills its own shell.
#   pgrep -x java             does not self-match, but matches ANY java. A
#                             long-running Cromwell server on this box blocked a
#                             pipeline resume for 17 minutes.
#   grep -c ERROR <log>       the log is appended across runs, so a stale failure
#                             from the previous run reads as a new one.
#   grep -o <pat> | wc -l     counts string matches, not processes. Reported
#                             "3 nextflow processes" when there was 1.
#   $? from `nextflow --help` non-zero because a param was missing, not because
#                             compilation failed.
#
# So this script uses only two sources, both unambiguous:
#   liveness -> kill -0 on a RECORDED pid. Cannot self-match, cannot catch an
#               unrelated process.
#   state    -> per-task .exitcode files in the work dir. Nextflow writes one per
#               task; the content IS the exit status. Filesystem truth.
#
# The run log is read ONLY to map a work-dir hash to a human-readable process
# name. That is a label lookup, not a state query — the distinction that the
# grep-the-log bug got wrong.
#
# Usage:
#   nf-status.sh <pidfile> <workdir> [baseline-file] [runlog]
#
#   baseline-file  optional. A list of already-failed task dirs, written by a
#                  previous invocation. Failures present in it are NOT reported
#                  again — this is what stops stale failures re-alerting forever.
#                  The file is refreshed on each run.
#
# Exit status: 0 always. Read the output; do not branch on $? — see above.
set -uo pipefail

PIDFILE="${1:?usage: nf-status.sh <pidfile> <workdir> [baseline] [runlog]}"
WORKDIR="${2:?usage: nf-status.sh <pidfile> <workdir> [baseline] [runlog]}"
BASELINE="${3:-}"
RUNLOG="${4:-}"

# ---- liveness -------------------------------------------------------------
if [ -r "$PIDFILE" ] && kill -0 "$(cat "$PIDFILE")" 2>/dev/null; then
    echo "alive=yes pid=$(cat "$PIDFILE")"
else
    echo "alive=no"
fi

# ---- per-task state from the filesystem -----------------------------------
done_n=0; fail_n=0
failed_list=$(mktemp)
while IFS= read -r f; do
    [ -z "$f" ] && continue
    # NOT `read -r code < "$f" || continue`. Nextflow writes .exitcode with NO
    # trailing newline ("127" is 3 bytes), and bash `read` returns NON-ZERO at
    # EOF-without-newline even though it HAS set the variable. The obvious
    # guard therefore skips every real task and reports zero of everything —
    # a silent wrong answer, which is the exact failure mode this script exists
    # to prevent. `cat` + strip whitespace is newline-agnostic.
    # For the same reason, never `find ... -exec cat {} +` and grep the result:
    # with no trailing newlines every exit code concatenates onto one line.
    code=$(tr -d '[:space:]' < "$f" 2>/dev/null)
    [ -z "$code" ] && continue
    if [ "$code" = "0" ]; then
        done_n=$((done_n + 1))
    else
        fail_n=$((fail_n + 1))
        printf '%s\t%s\n' "$(dirname "$f")" "$code" >> "$failed_list"
    fi
done < <(find "$WORKDIR" -name .exitcode 2>/dev/null)
sort -o "$failed_list" "$failed_list" 2>/dev/null || true

echo "tasks_completed=$done_n tasks_failed=$fail_n"

# ---- NEW failures only ----------------------------------------------------
# Without a baseline every poll re-reports the same old failures, which is how a
# stale error from a previous run gets announced as a fresh one.
if [ -n "$BASELINE" ]; then
    if [ -f "$BASELINE" ]; then
        new=$(comm -13 "$BASELINE" "$failed_list" 2>/dev/null)
    else
        new=""                      # first call: establish the baseline silently
    fi
    cp "$failed_list" "$BASELINE" 2>/dev/null || true
else
    new=$(cat "$failed_list")
fi

if [ -n "$new" ]; then
    echo "NEW_FAILURES:"
    printf '%s\n' "$new" | while IFS=$'\t' read -r dir code; do
        [ -z "$dir" ] && continue
        # hash -> process name. Label lookup only; state came from .exitcode.
        name=""
        if [ -n "$RUNLOG" ] && [ -r "$RUNLOG" ]; then
            short=$(basename "$(dirname "$dir")")/$(basename "$dir" | cut -c1-6)
            name=$(grep -F "[$short" "$RUNLOG" 2>/dev/null | grep -oP 'process > \K.*' | tail -1)
        fi
        echo "  exit=$code ${name:-<name unresolved>}"
        echo "    dir: $dir"
        if [ -s "$dir/.command.err" ]; then
            sed -n '1,3p' "$dir/.command.err" | sed 's/^/    err: /'
        fi
    done
fi

rm -f "$failed_list"
