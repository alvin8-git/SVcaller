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


def _matching_selectors(name):
    """Every conf/ selector that matches this process name."""
    return [(p, img) for p, img in _selectors() if re.fullmatch(p, name)]


def _resolve(name):
    """The container a process gets from conf/, or None.

    Returns the LAST match. An earlier draft returned the first, which is not
    what Nextflow does — later config settings take precedence — so a process
    matching two selectors would have been reported with the wrong image while
    the test passed. Rather than rely on getting that precedence right,
    test_no_process_matches_two_selectors makes the ambiguity itself a failure,
    so which one wins never has to be reasoned about.
    """
    hits = _matching_selectors(name)
    return hits[-1][1] if hits else None


def test_no_process_matches_two_selectors():
    """Overlapping withName selectors are a config smell regardless of which
    Nextflow applies: a human reading conf/ cannot tell which image a process
    gets. Ambiguity here is exactly how HBA_DEPTH would have been mis-resolved."""
    ambiguous = []
    for name, body, _ in _processes():
        if re.search(r"^\s*container\s", body, re.M):
            continue
        hits = _matching_selectors(name)
        if len(hits) > 1:
            ambiguous.append(f"{name}: matched by {[p for p, _ in hits]}")
    assert not ambiguous, "overlapping container selectors:\n  " + "\n  ".join(ambiguous)


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


UTILS_RE = re.compile(r"svcaller/utils:[\d.]+")


def test_utils_tag_is_consistent_everywhere():
    """The utils tag is written in 14 places. They must all agree.

    CLAUDE.md flagged this as a known hazard ("three modules hardcode it"), and
    it is exactly the kind of drift nothing detects until a run pulls a stale
    image. The authority is params.utils_container in nextflow.config; every
    other occurrence is a fallback default and must match it.
    """
    conf = open(os.path.join(REPO, "nextflow.config")).read()
    m = re.search(r"utils_container\s*=\s*'(svcaller/utils:[\d.]+)'", conf)
    assert m, "params.utils_container not found in nextflow.config"
    authority = m.group(1)

    mismatches = []
    for root, _, files in os.walk(REPO):
        if any(p in root for p in ("/work_", "/results", "/.git", "__pycache__")):
            continue
        for f in files:
            if not f.endswith((".nf", ".config")):
                continue
            path = os.path.join(root, f)
            for i, line in enumerate(open(path), 1):
                for tag in UTILS_RE.findall(line):
                    if tag != authority:
                        mismatches.append(
                            f"{os.path.relpath(path, REPO)}:{i} has {tag}, "
                            f"authority is {authority}")
    assert not mismatches, "utils container tag drift:\n  " + "\n  ".join(mismatches)


def test_utils_processes_honour_the_param():
    """A hardcoded tag with no `params.utils_container ?:` fallback silently
    ignores --utils_container. Three processes did exactly that (pycirclize,
    gridss/convert_bnd, report) while the parameter was documented as working."""
    offenders = []
    for root, _, files in os.walk(REPO):
        if any(p in root for p in ("/work_", "/results", "/.git", "__pycache__")):
            continue
        for f in files:
            if not f.endswith(".nf"):
                continue
            path = os.path.join(root, f)
            for i, line in enumerate(open(path), 1):
                if not UTILS_RE.search(line) or "container" not in line:
                    continue
                if "params.utils_container" not in line:
                    offenders.append(f"{os.path.relpath(path, REPO)}:{i}: {line.strip()}")
    assert not offenders, (
        "these hardcode the utils image and ignore --utils_container:\n  "
        + "\n  ".join(offenders))


def test_hba_depth_specifically_resolves_to_mosdepth():
    """The regression itself, pinned by name so the fix cannot be reverted quietly."""
    image = _resolve("HBA_DEPTH")
    assert image is not None, "HBA_DEPTH matches no conf/ selector"
    assert "mosdepth" in image, (
        f"HBA_DEPTH resolves to {image!r}; it shells out to mosdepth and the "
        "utils image does not have it — this is exit 127 on the first real run")
    # and its sibling must still be right
    assert "mosdepth" in (_resolve("TRAIT_DEPTH") or "")
