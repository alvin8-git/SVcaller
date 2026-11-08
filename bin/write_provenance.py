#!/usr/bin/env python3
"""Write <outdir>/<sample>/provenance.json (producing-software versions for OmniGen).

Non-blocking by contract: every failure is a warning on stderr and the exit code is
always 0.  Tool versions come from the tool itself, run inside its own container
(`docker run --rm <image> <cmd>`); if that is impossible the value is "not reported".
The repo's git commit is only ever used for the `svcaller_commit` field, never
for a tool version.  Images are those recorded by the run's conf/docker.config.
"""
import argparse
import datetime
import json
import os
import re
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NR = "not reported"

# tool -> (nextflow process name in docker.config, version command inside the image)
TOOLS = {
    "manta": ("MANTA_CALL", ["configManta.py", "--version"]),
    "delly": ("DELLY_CALL_SVTYPE", ["delly", "--version"]),
    "annotsv": ("ANNOTSV", ["AnnotSV", "--version"]),
    "expansionhunter": ("EXPANSIONHUNTER", ["ExpansionHunter", "--version"]),
    "gatk": ("DELLY_MERGE", ["gatk", "--version"]),
    "samtools": ("SAMTOOLS_SORT", ["samtools", "--version"]),
    "bwa-mem2": ("BWAMEM2_ALIGN", ["bwa-mem2", "version"]),
    "picard": ("PICARD_MARKDUP", ["picard", "MarkDuplicates", "--version"]),
    "mosdepth": ("MOSDEPTH", ["mosdepth", "--version"]),
    "svaba": ("SVABA_CALL", ["svaba", "--version"]),
    "strling": ("STRLING_CALL", ["strling", "--version"]),
    "cnvpytor": ("CNVPYTOR_CALL", ["cnvpytor", "--version"]),
    "jasmine": ("JASMINE_MERGE", ["jasmine", "--version"]),
    "truvari": ("TRUVARI_BENCH", ["truvari", "version"]),
    "gridss": ("GRIDSS_CALL", ["gridss", "--version"]),
    "scramble": ("SCRAMBLE_CALL", ["cluster_identifier", "--version"]),
}
VERSION_RE = re.compile(r"\d+\.\d+(?:\.\d+)*[\w+-]*")


def warn(msg):
    print(f"[provenance] WARNING: {msg}", file=sys.stderr)


def git(args):
    return subprocess.run(["git", "-C", REPO] + args, capture_output=True, text=True,
                          timeout=30, check=True).stdout.strip()


def repo_state():
    """(describe, commit) of the SVcaller checkout; commit gets -dirty when the tree is."""
    try:
        commit = git(["rev-parse", "HEAD"])
        dirty = bool(git(["status", "--porcelain", "--untracked-files=no"]))
        describe = git(["describe", "--tags", "--always"]) or commit[:7]
        if dirty:
            commit += "-dirty"
            describe += "-dirty"
        return describe, commit
    except Exception as e:  # noqa: BLE001
        warn(f"git state unavailable: {e}")
        return NR, NR


def parse_images(config_path):
    """process name -> container image, from `withName: 'A|B' { container = '...' }`."""
    images = {}
    try:
        text = open(config_path).read()
        text = re.sub(r"//[^\n]*", "", text)
        for names, image in re.findall(
                r"withName:\s*'([^']+)'\s*\{[^}]*?container\s*=\s*'([^']+)'", text):
            for n in names.split("|"):
                images[n.strip()] = image
    except Exception as e:  # noqa: BLE001
        warn(f"cannot read {config_path}: {e}")
    return images


def first_version_line(out):
    for line in out.splitlines():
        m = VERSION_RE.search(line)
        if m:
            return m.group(0)
    return None


def query_version(image, cmd, timeout=90):
    """Run the tool in its own (already local) image; None when it cannot report."""
    try:
        have = subprocess.run(["docker", "image", "inspect", image], capture_output=True,
                              timeout=30)
        if have.returncode != 0:
            warn(f"image {image} not present locally; not pulling")
            return None
        r = subprocess.run(["docker", "run", "--rm", image] + cmd, capture_output=True,
                           text=True, timeout=timeout)
        return first_version_line(r.stdout + "\n" + r.stderr)  # some tools exit 1 on --version
    except Exception as e:  # noqa: BLE001
        warn(f"version query failed for {image} {cmd[0]}: {e}")
        return None


def build(sample, ref_fasta, images_map, querier=query_version, state=None):
    describe, commit = state or repo_state()
    tools, images = {}, {}
    for tool, (proc, cmd) in TOOLS.items():
        image = images_map.get(proc)
        images[tool] = image or NR
        v = None
        if image:
            try:
                v = querier(image, cmd)
            except Exception as e:  # noqa: BLE001
                warn(f"{tool}: {e}")
        tools[tool] = v or NR
    return {
        "pipeline": f"SVcaller {describe}",
        "svcaller_commit": commit,
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "sample": sample,
        "reference": {"path": ref_fasta},
        "tools": tools,
        "images": images,
    }


def write(outdir, sample, ref_fasta, images_map, querier=query_version, state=None):
    """Returns the path written, or None.  Never raises."""
    try:
        doc = build(sample, ref_fasta, images_map, querier, state)
        d = os.path.join(outdir, sample)
        os.makedirs(d, exist_ok=True)
        path = os.path.join(d, "provenance.json")
        tmp = path + ".tmp"
        with open(tmp, "w") as f:
            json.dump(doc, f, indent=2)
            f.write("\n")
        os.replace(tmp, path)
        return path
    except Exception as e:  # noqa: BLE001
        warn(f"could not write provenance for {sample} under {outdir}: {e}")
        return None


def main(argv=None):
    try:
        ap = argparse.ArgumentParser()
        ap.add_argument("--outdir", required=True)
        ap.add_argument("--sample", action="append", required=True)
        ap.add_argument("--ref-fasta", default=NR)
        ap.add_argument("--docker-config", default=os.path.join(REPO, "conf", "docker.config"))
        a = ap.parse_args(argv)
        images_map = parse_images(a.docker_config)
        for s in a.sample:
            p = write(a.outdir, s, a.ref_fasta, images_map)
            if p:
                print(f"[provenance] wrote {p}")
    except SystemExit:
        warn("bad arguments")
    except Exception as e:  # noqa: BLE001
        warn(str(e))
    return 0


if __name__ == "__main__":
    sys.exit(main())
