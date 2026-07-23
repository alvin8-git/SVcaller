"""If nextflow.config enables trace/report/timeline, the custom container images
MUST install procps (ps).

WHY. Enabling trace/report/timeline makes Nextflow collect per-task resource
metrics by running `ps` inside each task's container. A container without `ps`
then fails EVERY task with "Command 'ps' required by nextflow to collect task
metrics cannot be found" (exit 1). This bit the SMA validation run: svcaller/utils
and svcaller/smncopynum lacked procps, so HBA_JUNCTION and SMN_CALLER died the
moment trace was turned on. Biocontainers ship ps; only our custom images are our
responsibility.
"""
import os
import re

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

CUSTOM_DOCKERFILES = ["Dockerfile.utils", "Dockerfile.smncopynum"]


def _read(rel):
    return open(os.path.join(REPO, rel)).read()


def _metrics_enabled(cfg):
    return any(
        re.search(rf"\b{scope}\s*{{[^}}]*enabled\s*=\s*true", cfg, re.S)
        for scope in ("trace", "report", "timeline"))


def _installs_procps(dockerfile_src):
    """True only if a NON-comment line names procps. A `# procps ...` comment must
    not count, or the test passes even when procps is dropped from the apt line
    (which is exactly how the first version of this test went vacuous)."""
    for line in dockerfile_src.splitlines():
        if line.strip().startswith("#"):
            continue
        if "procps" in line:
            return True
    return False


def test_metrics_detector_is_not_vacuous():
    """Guard the guard: if _metrics_enabled silently stopped matching, the lockstep
    test below would pass by skipping. Pin it against known blocks."""
    assert _metrics_enabled("trace { enabled = true\n file = 'x' }")
    assert _metrics_enabled("report   { enabled = true; overwrite = true }")
    assert not _metrics_enabled("trace { enabled = false }")
    assert not _metrics_enabled("params { foo = true }")


def test_installs_procps_ignores_comments():
    assert _installs_procps("RUN apt-get install -y procps")
    assert not _installs_procps("# procps provides ps\nRUN apt-get install -y samtools")


def test_metrics_collection_and_procps_stay_in_lockstep():
    cfg = _read("nextflow.config")
    if not _metrics_enabled(cfg):
        # No trace/report/timeline -> Nextflow never runs ps -> procps not required.
        return
    missing = [df for df in CUSTOM_DOCKERFILES if not _installs_procps(_read(df))]
    assert not missing, (
        f"nextflow.config enables trace/report/timeline (which run `ps` in every "
        f"task container) but these custom images do not install procps: {missing}. "
        "Every task in them will fail with \"Command 'ps' required by nextflow\".")
