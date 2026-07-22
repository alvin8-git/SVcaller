"""Every process must actually get a container that has its tools.

WHY THIS EXISTS. The first real end-to-end run of M8 died in seconds with

    .command.sh: line 16: mosdepth: command not found      (exit 127)

`HBA_DEPTH` shells out to mosdepth exactly like its sibling `TRAIT_DEPTH`, and
`subworkflows/alpha_globin.nf:20` even carries the comment "container assigned in
conf/docker.config: mosdepth biocontainer" — copied from TRAIT_DEPTH along with
the pattern. But the process was never added to the `withName:` selector, so it
silently fell through to the default utils container, which has samtools and
tabix and no mosdepth.

No unit test could catch that: the Python was correct, the module was correct,
and the comment asserted the wiring existed. It is a config-to-code binding, and
the only cheap way to check it is structurally.

These tests are pure text analysis — no Docker, no Nextflow, no network.
"""
import os
import re

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NF_DIRS = ("modules", "subworkflows", "workflows")
CONF_DIR = os.path.join(REPO, "conf")

# tool -> the processes allowed to need it are those whose resolved container
# provides it. We cannot run Docker here, so we assert the BINDING instead:
# a process that calls `tool` must resolve to a container selector known to
# carry it.
TOOL_SELECTORS = {
    "mosdepth": "mosdepth",     # the biocontainer image name must mention it
}


def _nf_files():
    out = []
    for d in NF_DIRS:
        for root, _, files in os.walk(os.path.join(REPO, d)):
            out += [os.path.join(root, f) for f in files if f.endswith(".nf")]
    return sorted(out)


def _processes():
    """[(process_name, body, path)] — body runs to the next top-level `process`."""
    found = []
    for path in _nf_files():
        src = open(path).read()
        parts = re.split(r"^process\s+([A-Z0-9_]+)\s*\{", src, flags=re.M)
        # parts = [preamble, name1, body1, name2, body2, ...]
        for i in range(1, len(parts), 2):
            found.append((parts[i], parts[i + 1], path))
    return found


def _selectors():
    """[(regex_alternation, container_string)] from every conf/*.config."""
    sels = []
    for f in sorted(os.listdir(CONF_DIR)):
        if not f.endswith(".config"):
            continue
        for line in open(os.path.join(CONF_DIR, f)):
            m = re.search(r"withName:\s*'([^']+)'\s*\{[^}]*container\s*=\s*'([^']+)'", line)
            if m:
                sels.append((m.group(1), m.group(2)))
    return sels


def _resolve(name):
    """The container a process gets from conf/, or None."""
    for pattern, image in _selectors():
        if re.fullmatch(pattern, name):
            return image
    return None


def test_every_process_has_a_container():
    """Inline `container` directive, or a conf/ selector that matches it."""
    orphans = []
    for name, body, path in _processes():
        if re.search(r"^\s*container\s", body, re.M):
            continue
        if _resolve(name) is None:
            orphans.append(f"{name} ({os.path.relpath(path, REPO)})")
    assert not orphans, (
        "these processes set no container inline and match no conf/ withName "
        "selector, so they fall through to the default image:\n  "
        + "\n  ".join(orphans))


@pytest.mark.parametrize("tool,marker", sorted(TOOL_SELECTORS.items()))
def test_processes_calling_a_tool_get_an_image_that_has_it(tool, marker):
    """A process whose script calls `tool` must resolve to an image named for it.

    This is the exact bug: HBA_DEPTH called mosdepth while resolving to the utils
    image. Matching on the image name is coarse, but it is the strongest check
    available without pulling and running every container.
    """
    wrong = []
    for name, body, path in _processes():
        # A mention in a COMMENT does not count. report.nf describes its channels
        # as "[ meta, mosdepth_summary ]" and "mosdepth 50kb windows" — prose
        # about where the data came from, not an invocation. Strip // comments
        # first, then require the tool to sit at a command position.
        script = re.sub(r"//.*$", "", body, flags=re.M)
        if not re.search(rf"(?:^|[|&;(]|\$\()\s*{tool}(?:\s|\\|$)", script, re.M):
            continue
        inline = re.search(r"^\s*container\s+.*?'([^']+)'", body, re.M)
        image = inline.group(1) if inline else _resolve(name)
        if image is None:
            wrong.append(f"{name}: calls {tool} but has NO container at all")
        elif marker not in image:
            wrong.append(f"{name}: calls {tool} but resolves to {image!r} "
                         f"({os.path.relpath(path, REPO)})")
    assert not wrong, "\n  " + "\n  ".join(wrong)


def test_hba_depth_specifically_resolves_to_mosdepth():
    """The regression itself, pinned by name so the fix cannot be reverted quietly."""
    image = _resolve("HBA_DEPTH")
    assert image is not None, "HBA_DEPTH matches no conf/ selector"
    assert "mosdepth" in image, (
        f"HBA_DEPTH resolves to {image!r}; it shells out to mosdepth and the "
        "utils image does not have it — this is exit 127 on the first real run")
    # and its sibling must still be right
    assert "mosdepth" in (_resolve("TRAIT_DEPTH") or "")
