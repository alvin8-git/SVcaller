"""Guard: REPORT normalizes every per-sample meta to [id:] before the join.

WHY. REPORT joins ~18 per-sample channels by their FULL meta map. If an upstream
stage tags some channels with an extra meta key but not others, the exact .join()
silently drops that sample; only the end count-guard notices, and it hard-fails
the whole run instead of producing the report. This happened for real: PICARD
metrics lacked the needs_chr_filter tag its siblings carried, and every
FASTQ-derived sample's report vanished. Keying the joins on a minimal [id:] meta
removes the divergence (nothing in REPORT reads any meta field but id). This test
keeps that normalization in place.
"""
import os
import re

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NF = os.path.join(REPO, "subworkflows", "report.nf")


def _src():
    return open(NF).read()


def _persample_take_channels(src):
    """Per-sample take channels are the ones declared `// [ meta, ...` (path-only
    channels are declared `// path`)."""
    take = src[src.index("take:"):src.index("\n    main:")]
    return re.findall(r"^\s*(ch_\w+)\s*//\s*\[\s*meta", take, re.M)


def test_report_defines_the_id_normalizer():
    assert re.search(r"\bbyId\b\s*=", _src()), \
        "the [id:] meta normalizer (byId) is gone from report.nf"


def test_every_per_sample_channel_is_id_normalized():
    src = _src()
    persample = _persample_take_channels(src)
    assert len(persample) >= 15, \
        f"expected the full set of per-sample channels, found only {persample}"
    missing = [c for c in persample if f"byId({c})" not in src]
    assert not missing, (
        "these per-sample channels are still joined by their full meta map, so an "
        f"added meta key would silently drop the sample: {missing}")


def test_the_landmine_channels_are_covered():
    """The specific channels behind past silent drops: PICARD metrics (missing
    needs_chr_filter), and the newest addition (alpha_globin)."""
    src = _src()
    for ch in ("ch_sv_tsv", "ch_metrics", "ch_alpha_globin"):
        assert f"byId({ch})" in src, f"{ch} must be id-normalized before the join"
