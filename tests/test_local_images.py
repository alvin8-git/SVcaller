"""Pinned local images must exist before launch.

`svcaller/*` images are built locally (Dockerfile.utils, Dockerfile.smncopynum), never pulled, so a
pinned tag that was never built dies mid-run with "pull access denied" (HG001_md run, 2026-11-07,
after 42 min).  Build commands: `docker build -f Dockerfile.utils -t svcaller/utils:<tag> .`
and `docker build -f Dockerfile.smncopynum -t svcaller/smncopynum:<tag> .`.
MELT is opt-in (`--run_melt`) and built only on demand (Dockerfile.melt), so it is exempt.
"""
import os
import re
import shutil
import subprocess

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OPTIN = ("svcaller/melt",)


def _src(*parts):
    return open(os.path.join(REPO, *parts)).read()


def _pinned():
    tags = set()
    for root, dirs, files in os.walk(REPO):
        dirs[:] = [d for d in dirs if d in ("conf", "modules", "subworkflows", "workflows")
                   or root == REPO and d == "conf"]
        for f in files:
            if f.endswith((".nf", ".config")):
                tags |= set(re.findall(r"svcaller/[a-z]+:[0-9][0-9.]*", open(os.path.join(root, f)).read()))
    tags |= set(re.findall(r"svcaller/[a-z]+:[0-9][0-9.]*", _src("nextflow.config")))
    return {t for t in tags if not t.startswith(OPTIN)}


def test_utils_tag_is_one_value_everywhere():
    utils = {t for t in _pinned() if t.startswith("svcaller/utils:")}
    assert len(utils) == 1, f"svcaller/utils pinned at several tags: {sorted(utils)}"


def test_pinned_images_exist_locally():
    if not shutil.which("docker"):
        pytest.skip("no docker CLI")
    if subprocess.run(["docker", "info"], capture_output=True).returncode != 0:
        pytest.skip("docker daemon not reachable")
    missing = [t for t in sorted(_pinned())
               if subprocess.run(["docker", "image", "inspect", t], capture_output=True).returncode != 0]
    assert not missing, (
        f"pinned image(s) not present locally: {missing}. Build with `docker build -f "
        "Dockerfile.utils -t svcaller/utils:<tag> .` (or Dockerfile.smncopynum) before launching.")
