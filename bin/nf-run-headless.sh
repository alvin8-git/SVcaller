#!/usr/bin/env bash
# Launch a Nextflow run FULLY DETACHED from the calling shell.
#
#   bin/nf-run-headless.sh <run-name> <nextflow args...>
#
# Why this exists. A run started with plain `nohup ... &` from an interactive
# session is still in that session's process group. When the terminal, the agent
# harness, or the ssh connection goes away, the group can be signalled and the
# run dies hours in — which happened here on 2026-07-22 with the THAL run.
#
# `setsid` puts the process in a NEW SESSION with no controlling terminal, so it
# cannot receive SIGHUP from the parent going away. stdin is closed and both
# streams are redirected, so nothing blocks on a terminal that no longer exists.
#
# It also records the PID, because every other way of finding a run again turns
# out to be wrong — see the failure table in CLAUDE.md and bin/nf-status.sh.
#
#   log  -> ${TMP}/<run-name>.log      (appended; runs are resumable)
#   pid  -> ${TMP}/<run-name>.pid      (read by bin/nf-status.sh)
#
# Check on it with:
#   bin/nf-status.sh ${TMP}/<run-name>.pid <workdir> <baseline> ${TMP}/<run-name>.log
set -uo pipefail

NAME="${1:?usage: nf-run-headless.sh <run-name> <nextflow args...>}"
shift
[ $# -gt 0 ] || { echo "no nextflow arguments given" >&2; exit 2; }

TMP="${TMPDIR:-/data/alvin/tmp}"
mkdir -p "$TMP"
LOG="$TMP/${NAME}.log"
PIDFILE="$TMP/${NAME}.pid"

# Refuse to start a second run over the same state. Two Nextflow processes on one
# work dir fight over the session lock and corrupt the resume cache. Liveness is
# `kill -0` on a recorded pid — never a process-name match, which would catch this
# very script or an unrelated JVM.
if [ -r "$PIDFILE" ] && kill -0 "$(cat "$PIDFILE")" 2>/dev/null; then
    echo "refusing to start: run '$NAME' is already alive (pid $(cat "$PIDFILE"))" >&2
    echo "  stop it first, or use a different run name" >&2
    exit 1
fi

# NXF_ANSI_LOG=false is mandatory without a TTY: the ANSI renderer deadlocks every
# JVM thread when there is no terminal, and the run hangs with no output at all.
setsid env NXF_ANSI_LOG=false nextflow "$@" >>"$LOG" 2>&1 </dev/null &
child=$!
disown "$child" 2>/dev/null || true

# The pid recorded is the `nextflow` launcher, which is what actually holds the
# session lock and what `kill -0` should track.
echo "$child" > "$PIDFILE"
sleep 2
if kill -0 "$child" 2>/dev/null; then
    echo "started '$NAME' pid=$child"
    echo "  log: $LOG"
    echo "  pid: $PIDFILE"
else
    echo "FAILED to start '$NAME' — last log lines:" >&2
    tail -20 "$LOG" >&2
    rm -f "$PIDFILE"
    exit 1
fi
