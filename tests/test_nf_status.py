"""Guard bin/nf-status.sh — the pipeline status reporter.

A status tool that reports something FALSE is worse than one that reports
nothing: a wrong status looks like a finding and gets acted on. On 2026-07-22 a
whole sequence of ad-hoc checks did exactly that (see the header of
bin/nf-status.sh), and the replacement written to fix them shipped with the same
class of bug on its first draft:

    read -r code < "$f" || continue

Nextflow writes `.exitcode` with NO trailing newline — "127" is three bytes.
bash `read` returns NON-ZERO at EOF-without-newline *even though it has set the
variable*, so that guard skipped every real task and reported
`tasks_completed=0 tasks_failed=0` while 15 tasks sat on disk. It was caught only
because the output was compared against an independently-counted ground truth
rather than believed.

So these tests write `.exitcode` files WITHOUT trailing newlines, exactly as
Nextflow does. A version of the script that regresses to `read` will report zero
and fail here.
"""
import os
import shutil
import subprocess

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO, "bin", "nf-status.sh")

pytestmark = pytest.mark.skipif(
    not shutil.which("bash"), reason="bash not available")


def _task(workdir, name, code, err=None):
    """Create a Nextflow-shaped task dir. NO trailing newline on .exitcode —
    that is the real on-disk format and the whole point of this fixture."""
    d = os.path.join(workdir, name[:2], name[2:])
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, ".exitcode"), "w") as fh:
        fh.write(str(code))            # deliberately no "\n"
    if err:
        with open(os.path.join(d, ".command.err"), "w") as fh:
            fh.write(err + "\n")
    return d


def _run(pidfile, workdir, baseline="", runlog=""):
    r = subprocess.run(["bash", SCRIPT, pidfile, workdir, baseline, runlog],
                       capture_output=True, text=True)
    return r.stdout


def _counts(out):
    done = int(out.split("tasks_completed=")[1].split()[0])
    fail = int(out.split("tasks_failed=")[1].split()[0])
    return done, fail


@pytest.fixture
def work(tmp_path):
    w = tmp_path / "work"
    w.mkdir()
    for i in range(3):
        _task(str(w), f"aa{i:04d}bbbb", 0)
    _task(str(w), "ff0001cccc", 127, err="mosdepth: command not found")
    return str(w)


def test_counts_exitcode_files_without_trailing_newline(work, tmp_path):
    """The regression. 3 completed + 1 failed must not read as 0 + 0."""
    dead = tmp_path / "dead.pid"
    dead.write_text("999999")
    done, fail = _counts(_run(str(dead), work))
    assert (done, fail) == (3, 1), (
        f"got completed={done} failed={fail}; if both are 0 the script is using "
        "`read` and silently skipping every task")


def test_dead_pid_reports_not_alive(work, tmp_path):
    dead = tmp_path / "dead.pid"
    dead.write_text("999999")
    assert "alive=no" in _run(str(dead), work)


def test_live_pid_reports_alive(work, tmp_path):
    """kill -0 on a real pid. Uses this process's own pid — no pattern matching,
    which is the property that makes it immune to the self-match failures."""
    live = tmp_path / "live.pid"
    live.write_text(str(os.getpid()))
    assert "alive=yes" in _run(str(live), work)


def test_missing_pidfile_is_not_alive(work, tmp_path):
    assert "alive=no" in _run(str(tmp_path / "nope.pid"), work)


def test_baseline_suppresses_stale_failures(work, tmp_path):
    """The other half of the 2026-07-22 mess: a stale failure from a PREVIOUS run
    re-announced as a new one on every poll."""
    dead = tmp_path / "dead.pid"
    dead.write_text("999999")
    base = str(tmp_path / "baseline")

    first = _run(str(dead), work, baseline=base)
    assert "NEW_FAILURES" not in first, "first call must only establish the baseline"

    second = _run(str(dead), work, baseline=base)
    assert "NEW_FAILURES" not in second, "an unchanged failure set must stay silent"


def test_new_failure_is_reported(work, tmp_path):
    dead = tmp_path / "dead.pid"
    dead.write_text("999999")
    base = str(tmp_path / "baseline")
    _run(str(dead), work, baseline=base)               # prime

    _task(work, "dd9999eeee", 9, err="synthetic crash")
    out = _run(str(dead), work, baseline=base)
    assert "NEW_FAILURES" in out, "a genuinely new failure must surface"
    assert "exit=9" in out
    assert "synthetic crash" in out, "the stderr excerpt makes it diagnosable"


def test_always_exits_zero(work, tmp_path):
    """Callers must read the OUTPUT, not branch on $?. A non-zero exit from an
    unrelated internal command is how `nextflow --help` got misread as a compile
    failure."""
    dead = tmp_path / "dead.pid"
    dead.write_text("999999")
    r = subprocess.run(["bash", SCRIPT, str(dead), work], capture_output=True)
    assert r.returncode == 0


def test_script_does_not_use_pattern_matching_for_state():
    """The failures this script replaces all came from pgrep/pkill -f matching
    the checking command's own command line. It must not reintroduce them."""
    src = open(SCRIPT).read()
    body = src.split("set -uo pipefail", 1)[1]     # skip the explanatory header
    for banned in ("pgrep", "pkill", "killall"):
        assert banned not in body, (
            f"nf-status.sh uses {banned}; liveness must be `kill -0` on a "
            "recorded pid, which cannot match the checking process itself")


def test_missing_workdir_refuses_instead_of_reporting_zero(tmp_path):
    """A mistyped work dir must NOT report `tasks_completed=0 tasks_failed=0`.

    `find` on a nonexistent path prints nothing, so the counters stay at zero and
    the output is indistinguishable from a run that has genuinely done nothing.
    That happened for real: `work_thal` vs `work_THAL` reported a healthy 63-task
    run as zero of everything, and it was believed for a moment. Worse, the
    invocation went on to overwrite the failure baseline with an empty list, so
    the next correct call re-alerted three already-known failures as new.

    Refusing (exit 2, message on stderr) is the only safe answer.
    """
    dead = tmp_path / "dead.pid"
    dead.write_text("999999")
    baseline = tmp_path / "baseline"
    baseline.write_text("preexisting\t1\n")
    r = subprocess.run(
        ["bash", SCRIPT, str(dead), str(tmp_path / "work_typo"), str(baseline)],
        capture_output=True, text=True)
    assert r.returncode == 2, f"expected refusal, got rc={r.returncode}"
    assert "tasks_completed" not in r.stdout, \
        "reported task counts for a work dir that does not exist"
    assert "does not exist" in r.stderr
    # and it must not have clobbered the baseline on its way out
    assert baseline.read_text() == "preexisting\t1\n", \
        "a refused call still overwrote the failure baseline"
