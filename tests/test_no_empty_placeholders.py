"""Regression guard: SVcaller must never publish an empty placeholder on failure.

WHY THIS TEST EXISTS
--------------------
SVcaller processes used to end their command blocks with an `|| touch <output>`
fallback. When a caller FAILED, the process still exited 0 and published a
ZERO-BYTE output file. Downstream, the OmniGen consumer report gated its reads on
`os.path.exists()` -- which is True for a 0-byte file -- so a crashed SMN caller
was rendered as a complete, clean-looking consumer report reading
"0 Carrier findings, 0 Medication flags, Clear".

A crashed caller must never be indistinguishable from a negative result. These
tests statically assert the failure-masking patterns stay dead.

Note the distinction this suite deliberately preserves:
  * a VCF/TSV with a header and zero data rows is a LEGITIMATE empty result
    ("we looked, we found nothing") and is still allowed;
  * a zero-byte file, or a tool crash swallowed by `|| true`, is a MASKED FAILURE
    and is not.
"""

import pathlib
import re

import pytest

REPO = pathlib.Path(__file__).resolve().parent.parent
NF_DIRS = ["modules", "subworkflows", "workflows"]


def _nf_files():
    files = []
    for d in NF_DIRS:
        files.extend(sorted((REPO / d).rglob("*.nf")))
    files.append(REPO / "main.nf")
    return [f for f in files if f.exists()]


def _lines(path):
    """Script lines with comments and stderr diagnostics stripped.

    The modules deliberately *describe* the banned patterns in comments and in their
    error messages ("Not writing a fake '{}'"), so scan executable lines only.
    """
    out = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith("#") or stripped.startswith("//"):
            continue
        if ">&2" in line:  # diagnostic written to stderr, not a published artifact
            continue
        out.append(line)
    return out


def test_nf_files_found():
    """Sanity: the scan below is actually looking at something."""
    assert len(_nf_files()) > 10


def test_no_touch_fallback_anywhere():
    """No process may `touch` an output file as a fallback.

    `|| touch out.tsv` converts a caller crash into a published zero-byte file.
    This is the exact pattern that caused the OmniGen clean-bill-of-health bug
    (modules/smn_caller/call.nf and modules/annotsv/annotate.nf).
    """
    offenders = []
    for f in _nf_files():
        for n, line in enumerate(_lines(f), 1):
            if re.search(r"\|\|\s*\\?\s*touch\b", line):
                offenders.append(f"{f.relative_to(REPO)}:{n}: {line.strip()}")
            # a bare `touch <something>` used to create a declared output
            elif re.match(r"^\s*touch\s+\S", line):
                offenders.append(f"{f.relative_to(REPO)}:{n}: {line.strip()}")
    assert not offenders, (
        "Found `touch` fallbacks that publish empty placeholder outputs:\n"
        + "\n".join(offenders)
    )


def test_no_fake_empty_json_fallback():
    """`echo '{}' > out.json` is the JSON flavour of the same bug."""
    offenders = []
    for f in _nf_files():
        for n, line in enumerate(_lines(f), 1):
            if re.search(r"echo\s+['\"]\{\}['\"]\s*>", line):
                offenders.append(f"{f.relative_to(REPO)}:{n}: {line.strip()}")
    assert not offenders, (
        "Found fake-empty-JSON fallbacks:\n" + "\n".join(offenders)
    )


# Tool invocations whose non-zero exit must propagate. Swallowing these with
# `|| true` let a crashed caller fall through to an "empty result" code path.
CALLER_INVOCATIONS = {
    "modules/svaba/call.nf": "svaba run",
    "modules/scramble/call.nf": "scramble.sh",
}


@pytest.mark.parametrize("relpath,tool", sorted(CALLER_INVOCATIONS.items()))
def test_caller_invocation_not_swallowed(relpath, tool):
    """The caller's own exit code must not be discarded with `|| true`."""
    path = REPO / relpath
    text = path.read_text()

    # Match the *executable* invocation (start of line), not a mention in a comment.
    m = re.search(rf"^[ \t]*{re.escape(tool)}\b", text, re.M)
    assert m, f"{relpath}: could not find the `{tool}` invocation"

    # The command continues across backslash-continued lines; collect until a line
    # that does not end with a backslash.
    command_lines = []
    for line in text[m.start():].splitlines():
        command_lines.append(line)
        if not line.rstrip().endswith("\\"):
            break
    command = "\n".join(command_lines)

    assert "|| true" not in command, (
        f"{relpath}: `{tool}` exit code is swallowed by `|| true`. "
        "A crashed caller would then be published as an empty result.\n"
        f"Offending command:\n{command}"
    )


def test_melt_does_not_exit_zero_on_missing_install():
    """A missing MELT install is a misconfiguration, not 'zero MEIs found'.

    MELT used to print a WARNING, write a header-only VCF and `exit 0` when
    MELT.jar or the ME references were absent -- making 'MELT was never installed'
    look identical to 'this genome has no mobile element insertions'.
    Skipping MELT on purpose must go through --skip_melt / MELT_STUB instead.
    """
    text = (REPO / "modules" / "melt" / "call.nf").read_text()
    assert "emitting empty VCF" not in text
    assert "--skip_melt" in text, "MELT failure path should point users at --skip_melt"
    assert "exit 1" in text, "MELT must hard-fail on a missing install"


def test_smn_caller_hard_fails_on_missing_output():
    """The originating site of the incident: SMN must fail, not touch an empty TSV."""
    path = REPO / "modules" / "smn_caller" / "call.nf"
    code = "\n".join(_lines(path))
    text = path.read_text()
    assert "touch" not in code, "SMN_CALLER must never touch a placeholder .smn.tsv"
    assert "exit 1" in text, "SMN_CALLER must hard-fail when it produced no result table"
    # It must also reject a header-only table (caller ran but genotyped nothing).
    assert "-lt 2" in text, "SMN_CALLER must reject a header-only (no sample row) TSV"


def test_expansionhunter_hard_fails_on_missing_profile():
    """EH publishes str_profile.json to results/ and OmniGen reads it directly.

    A `{}` profile is what a broken ExpansionHunter looks like, not a genome with no
    repeat expansions, so it must never be published.
    """
    path = REPO / "modules" / "expansionhunter" / "call.nf"
    code = "\n".join(_lines(path))
    assert "echo '{}'" not in code
    assert "exit 1" in path.read_text(), "EXPANSIONHUNTER must hard-fail on a missing profile"


def test_annotsv_distinguishes_empty_input_from_failure():
    """AnnotSV: zero SVs in -> header-only TSV (legit). SVs in, nothing out -> fail."""
    text = (REPO / "modules" / "annotsv" / "annotate.nf").read_text()
    assert "|| touch" not in text
    assert "exit 1" in text, "AnnotSV must hard-fail when it silently produced nothing"
    # The legitimate empty-input path must still emit a real header, not a 0-byte file.
    assert "SV_chrom" in text


def test_svaba_declares_classic_bwa_index_input():
    """SvABA must declare the classic bwa index as a staged input.

    WHY: SvABA calls bwa_idx_load_from_disk internally and needs the CLASSIC bwa index
    (ref_fasta.{amb,ann,bwt,pac,sa}) symlinked next to ${ref_fasta} in the task work dir.
    For the entire history of the repo the SVABA_CALL module declared only `ref_fasta` +
    `ref_fai`, so Nextflow never staged the index and SvABA died with
    "[E::bwa_idx_load_from_disk] fail to locate the index files" -- a failure masked by a
    since-removed `2>&1 || true`. If the index input declaration is dropped again, SvABA
    silently stops staging its index and contributes nothing. This test locks it in.
    """
    text = (REPO / "modules" / "svaba" / "call.nf").read_text()

    # Isolate the SVABA_CALL process (not the STUB) input block.
    call_start = text.index("process SVABA_CALL")
    stub_start = text.index("process SVABA_STUB")
    call_block = text[call_start:stub_start]

    input_start = call_block.index("input:")
    output_start = call_block.index("output:")
    input_block = call_block[input_start:output_start]

    # A `path` input carrying the bwa index must be declared so Nextflow symlinks the
    # .amb/.ann/.bwt/.pac/.sa files into the task work dir alongside ${ref_fasta}.
    assert re.search(r"^\s*path\s+bwa_index\b", input_block, re.M), (
        "SVABA_CALL must declare `path bwa_index` so Nextflow stages the classic bwa "
        "index next to ref_fasta; without it bwa_idx_load_from_disk fails on every run."
    )


def test_main_fails_loud_when_classic_bwa_index_absent():
    """Absent classic bwa index must surface an actionable error, not a silent failure.

    Unless --skip_svaba is set, main.nf must check the index exists and tell the operator
    exactly how to fix it (build with `bwa index` or pass --skip_svaba).
    """
    text = (REPO / "main.nf").read_text()
    assert "skip_svaba" in text, "main.nf must gate the index requirement on --skip_svaba"
    assert "bwa index" in text, "the error must tell operators to run `bwa index <ref>`"
    assert "classic BWA index not found" in text, (
        "main.nf must emit a clear 'classic BWA index not found' error when SvABA runs "
        "without its index."
    )
