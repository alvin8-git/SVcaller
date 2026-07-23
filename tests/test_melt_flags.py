"""Guard the MELT command line against flags MELT does not accept.

WHY. `modules/melt/call.nf` passed `-reads ${params.melt_min_reads}`. MELT v2.2.2
has no such option: it answered "Unrecognized option: -reads", printed its usage
text and exited non-zero, for every ME type, on every sample.

The module's own guard then did the right thing and refused to publish an empty
VCF as "no mobile element insertions" — so the failure was loud. But it was loud
in a way nobody had read: `-reads` landed 2026-06-05 and the run16 F1 baseline is
dated 2026-06-08, so MELT was already failing when that baseline was measured.
The "5-caller ensemble" behind those numbers was in fact four callers, exactly as
SvABA turned out to be absent from a documented six.

A one-token typo in a shell heredoc cost a caller and silently changed what the
published benchmark meant. This test makes the same mistake fail in 0.1 s.

The accepted set is MELT v2.2.2's own `-help` output, captured from a real
container run on 2026-07-22.
"""
import os
import re

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MELT_NF = os.path.join(REPO, "modules", "melt", "call.nf")
MELT_MERGE_NF = os.path.join(REPO, "modules", "melt", "merge.nf")

# Verbatim from MELT v2.2.2 MELT-Single usage.
MELT_SINGLE_FLAGS = {
    "-a", "-ac", "-b", "-bamfile", "-bowtie", "-c", "-cov", "-d", "-e",
    "-exome", "-h", "-j", "-k", "-mcmq", "-n", "-nocleanup", "-q", "-r",
    "-s", "-samtools", "-sr", "-t", "-w", "-z",
}


def _melt_invocation_flags():
    """Flags passed to `java -jar MELT.jar Single` in the module.

    The script block is a Groovy string containing shell, so this reads the
    continuation lines of the java invocation rather than trying to parse either
    language properly.
    """
    src = open(MELT_NF).read()
    # The jar is invoked through a shell variable (`java ... -jar "$melt_jar"
    # Single`), NOT as a literal MELT.jar path — that string only appears in the
    # `find -name "MELT.jar"` lookup above it. Anchoring on the literal made this
    # regex silently match nothing.
    m = re.search(r"java\s+-Xmx.*?-jar\s+\S+\s+Single(.*?)2>&1", src, re.S)
    assert m, "could not locate the MELT Single invocation in modules/melt/call.nf"
    block = m.group(1)
    # a flag is a token starting with a single dash at the head of a line
    return {f for f in re.findall(r"^\s*(-[A-Za-z]\w*)\s", block, re.M)}


def test_melt_flags_are_all_accepted_by_melt():
    used = _melt_invocation_flags()
    assert used, "parsed no flags — the regex has drifted from the module"
    unknown = sorted(used - MELT_SINGLE_FLAGS)
    assert not unknown, (
        f"MELT v2.2.2 does not accept {unknown}. It will print its usage text and "
        f"exit non-zero for every ME type, on every sample. Accepted flags: "
        f"{sorted(MELT_SINGLE_FLAGS)}")


def test_the_specific_regression_is_gone():
    """`-reads` is ExpansionHunter's flag, not MELT's. It is easy to reintroduce
    because the neighbouring module legitimately uses `--reads`."""
    assert "-reads" not in _melt_invocation_flags(), \
        "MELT is being passed -reads again; the correct flag is -sr"


def test_min_reads_maps_to_the_split_read_filter():
    """params.melt_min_reads is documented as split-read support, so it must land
    on -sr. Sending it to -r (read length) or -c (coverage) would be silently
    wrong rather than an error — MELT accepts all three."""
    src = open(MELT_NF).read()
    m = re.search(r"^\s*(-[A-Za-z]\w*)\s+\$\{params\.melt_min_reads\}", src, re.M)
    assert m, "params.melt_min_reads is not passed to any MELT flag"
    assert m.group(1) == "-sr", (
        f"melt_min_reads is passed to {m.group(1)!r}; it is a split-read count "
        "and belongs on -sr")


@pytest.mark.parametrize("param,flag", [("min_depth", "-c")])
def test_other_params_land_on_the_flags_the_docs_claim(param, flag):
    """CLAUDE.md states min_depth doubles as MELT's -c coverage argument."""
    src = open(MELT_NF).read()
    m = re.search(rf"^\s*(-[A-Za-z]\w*)\s+\$\{{params\.{param}\}}", src, re.M)
    assert m, f"params.{param} is not passed to MELT"
    assert m.group(1) == flag, f"params.{param} is on {m.group(1)!r}, expected {flag!r}"


def test_the_output_contract_moved_to_merge_not_lost():
    """MELT parallelisation (2026-07-23) split the serial task into a per-type
    MELT_CALL fan-out plus a MELT_MERGE collector. The per-sample output contract —
    canonical-chr + PASS filter, INFO stripped to SVTYPE/MEITYPE/SVLEN/END, and the
    injected MEITYPE header — used to live at the tail of call.nf and now lives in
    merge.nf. If a future refactor drops it, the merged VCF silently loses MEITYPE
    (breaking bcftools/Truvari, exactly the failure the merge note warns about) or
    stops filtering to canonical/PASS. Guard that the contract still exists and that
    call.nf, now raw per-type, no longer carries it."""
    assert os.path.exists(MELT_MERGE_NF), "modules/melt/merge.nf is missing"
    merge = open(MELT_MERGE_NF).read()
    assert "ID=MEITYPE" in merge, \
        "MELT_MERGE no longer injects the MEITYPE INFO header; the merged VCF will " \
        "lose MEITYPE and bcftools/Truvari will fatal-exit"
    assert "PASS" in merge and "chr([0-9]" in merge, \
        "MELT_MERGE no longer filters to canonical chr + PASS"
    # The final per-sample VCF must be produced by exactly one place. call.nf now emits
    # raw per-type VCFs only; it must NOT also build the .melt.vcf.gz.
    call = open(MELT_NF).read()
    assert "${meta.id}.melt.vcf.gz" not in call, \
        "call.nf still builds the merged .melt.vcf.gz; the merge belongs in merge.nf"


def test_melt_call_is_per_type_fan_out():
    """The parallelisation runs one ME type per task via a `val(metype)` input and a
    tag of the form id:TYPE, mirroring DELLY. If this reverts to a serial loop the
    ~1.2 h → ~30 min win is gone; keep the fan-out shape explicit."""
    call = open(MELT_NF).read()
    assert "val(metype)" in call, \
        "MELT_CALL lost its per-type val(metype) input; it is no longer a fan-out"
    assert re.search(r'tag\s+"\$\{meta\.id\}:\$\{metype\}"', call), \
        "MELT_CALL tag is not id:metype; the per-type fan-out shape has drifted"
