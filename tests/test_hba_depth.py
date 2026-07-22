"""Deterministic unit tests for alpha-globin channel 1 (bin/hba_depth.py).

Pure and fast: no containers, no network, no BAMs, no reference. The calibration
inputs are the committed TSV/BED assets, read at test time, so a future edit to
`validation/giab_alpha_baseline.tsv` or `assets/hba_segments.bed` that breaks the
thresholds fails this suite instead of silently changing every call.
"""
import gzip
import sys
from pathlib import Path

import pytest

REPO = Path("/data/alvin/SVcaller")
sys.path.insert(0, str(REPO / "bin"))

import cnv_traits_common as c  # noqa: E402
import hba_depth  # noqa: E402

SEGMENTS_BED = REPO / "assets" / "hba_segments.bed"
GIAB_BASELINE = REPO / "validation" / "giab_alpha_baseline.tsv"

#: the four segments that carry a baseline; INTER_Z_A is deliberately absent
SCORED_SEGMENTS = ("HBZ", "HBA2", "INTER_A2_A1", "HBA1")

#: Raw depth ratios (segment / chr2:100000000-100020000 control), reproduced
#: 2026-07-22 from /data/alvin/ref/THAL/*_30X.bwa.sortdup.bqsr.bam with
#: `samtools depth -a -r <region>`:
#:   THAL1 ctrl 30.245 · HBZ 26.054 · HBA2 11.086 · INTER 15.273 · HBA1 11.987
#:   THAL2 ctrl 31.363 · HBZ 22.274 · HBA2 25.275 · INTER 31.100 · HBA1 28.224
THAL1_RATIOS = {"HBZ": 0.86, "HBA2": 0.37, "INTER_A2_A1": 0.50, "HBA1": 0.40}
THAL2_RATIOS = {"HBZ": 0.71, "HBA2": 0.81, "INTER_A2_A1": 0.99, "HBA1": 0.90}


# --------------------------------------------------------------------------- #
# Fixtures / helpers
# --------------------------------------------------------------------------- #
def load_segments():
    return {s["name"]: s for s in hba_depth.parse_segments(SEGMENTS_BED)}


def read_giab_baseline():
    """[(sample, {segment: raw_ratio})] from the committed calibration TSV."""
    rows = []
    header = None
    for line in GIAB_BASELINE.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split("\t")
        if header is None:
            header = f
            continue
        rec = dict(zip(header, f))
        rows.append((rec["sample"],
                     {seg: float(rec[seg + "_r"]) for seg in SCORED_SEGMENTS}))
    return rows


def write_depth(tmp_path, segment_means, control_depth=30.0, n_ctrl=8):
    """A mosdepth-style regions.bed.gz carrying alpha segments + CTRL_* windows."""
    coords = {s["name"]: (s["chrom"], s["start"], s["end"])
              for s in hba_depth.parse_segments(SEGMENTS_BED)}
    p = tmp_path / "s.alpha_depth.regions.bed.gz"
    with gzip.open(p, "wt") as fh:
        for label, mean in segment_means.items():
            chrom, start, end = coords.get(label, ("chr16", 1000, 2000))
            fh.write("{}\t{}\t{}\t{}\t{}\n".format(chrom, start, end, label, mean))
        for i in range(1, n_ctrl + 1):
            fh.write("chr{}\t100000000\t100020000\tCTRL_{}\t{}\n".format(
                i + 1, i, control_depth))
    return p


# --------------------------------------------------------------------------- #
# parse_segments
# --------------------------------------------------------------------------- #
def test_parse_segments_reads_the_committed_bed():
    segs = load_segments()
    assert set(segs) == {"HBZ", "INTER_Z_A", "HBA2", "INTER_A2_A1", "HBA1"}
    assert segs["HBA2"]["chrom"] == "chr16"
    assert (segs["HBA2"]["start"], segs["HBA2"]["end"]) == (172875, 173710)
    assert segs["HBA2"]["baseline"] == pytest.approx(0.750)
    assert segs["INTER_A2_A1"]["baseline"] == pytest.approx(1.001)
    assert segs["HBA1"]["baseline"] == pytest.approx(0.964)
    assert segs["HBZ"]["baseline"] == pytest.approx(0.760)
    # trailing '#' notes must not leak into the reliability field
    assert segs["HBA2"]["reliability"] == "needs_own_baseline"
    assert segs["INTER_A2_A1"]["reliability"] == "good"


def test_parse_segments_baseline_na_is_none():
    segs = load_segments()
    assert segs["INTER_Z_A"]["baseline"] is None
    assert segs["INTER_Z_A"]["reliability"] == hba_depth.REL_DO_NOT_AVERAGE


def test_parse_segments_tolerates_na_and_missing_reliability(tmp_path):
    bed = tmp_path / "seg.bed"
    bed.write_text("# header\nchr16\t10\t20\tFOO\tNA\t\t# note\n"
                   "chr16\t20\t30\tBAR\t0.5\tgood\n")
    segs = {s["name"]: s for s in hba_depth.parse_segments(bed)}
    assert segs["FOO"]["baseline"] is None
    assert segs["FOO"]["reliability"] == hba_depth.REL_NO_BASELINE
    assert segs["BAR"]["baseline"] == pytest.approx(0.5)


# --------------------------------------------------------------------------- #
# score / classify
# --------------------------------------------------------------------------- #
def test_score_is_ratio_over_baseline():
    assert hba_depth.score_segment(0.37, 0.750) == pytest.approx(0.4933, abs=1e-4)
    assert hba_depth.score_segment(None, 0.750) is None
    assert hba_depth.score_segment(0.37, None) is None
    assert hba_depth.score_segment(0.37, 0.0) is None


def test_classify_boundaries():
    assert hba_depth.classify(0.0) == "hom_loss"
    assert hba_depth.classify(0.249) == "hom_loss"
    assert hba_depth.classify(0.25) == "het_loss"
    assert hba_depth.classify(0.5) == "het_loss"
    assert hba_depth.classify(0.649) == "het_loss"
    assert hba_depth.classify(0.65) == "intact"
    assert hba_depth.classify(1.0) == "intact"
    assert hba_depth.classify(1.349) == "intact"
    assert hba_depth.classify(1.35) == "gain"
    assert hba_depth.classify(1.5) == "gain"
    assert hba_depth.classify(None) == "uncalibrated"


def test_is_marginal():
    assert hba_depth.is_marginal(0.66) is True     # just above the het/intact cut
    assert hba_depth.is_marginal(0.60) is True
    assert hba_depth.is_marginal(1.00) is False
    assert hba_depth.is_marginal(1.30) is True     # just below the gain cut
    assert hba_depth.is_marginal(None) is False
    # no calibration observation may be marginal, or the thresholds are too tight
    assert hba_depth.is_marginal(0.826) is False   # lowest GIAB normal (HG004 HBZ)
    assert hba_depth.is_marginal(1.189) is False   # highest GIAB normal (HG005 HBZ)
    assert hba_depth.is_marginal(0.500) is False   # highest THAL1 het loss


# --------------------------------------------------------------------------- #
# THE regression that motivated the whole design
# --------------------------------------------------------------------------- #
def test_every_giab_normal_scores_intact_on_all_four_segments():
    """A naive `ratio < 0.8 = loss` calls a het loss in ALL SIX GIAB normals.

    HG001 is excluded: it reads low across the whole locus (0.38/0.67/0.83/0.57
    scored) with a normal chr2 control, in a non-contiguous pattern no deletion
    can produce — technical dropout, per the header of the baseline TSV. It is
    neither a normal nor a carrier, so it must not be asserted either way.
    """
    segs = load_segments()
    rows = [r for r in read_giab_baseline() if r[0] != "HG001"]
    assert len(rows) == 6, "baseline TSV should carry HG002-HG007"
    for sample, ratios in rows:
        for name in SCORED_SEGMENTS:
            score, call = hba_depth.call_segment(segs[name], ratios[name])
            assert call == "intact", (
                "{} {} ratio={} score={:.3f} called {} — a normal must never "
                "be called a carrier".format(sample, name, ratios[name], score, call))
            # and the naive raw-ratio rule this replaces would have been wrong:
            # intact HBA2 sits at 0.665-0.791, i.e. under 0.8 in every normal
            if name == "HBA2":
                assert ratios[name] < 0.8


def test_giab_normals_are_not_marginal():
    """Every GIAB normal must clear the boundary with room, not scrape past it."""
    segs = load_segments()
    for sample, ratios in read_giab_baseline():
        if sample == "HG001":
            continue
        for name in SCORED_SEGMENTS:
            score, _ = hba_depth.call_segment(segs[name], ratios[name])
            assert not hba_depth.is_marginal(score), \
                "{} {} score {:.3f} sits on a boundary".format(sample, name, score)


def test_naive_raw_ratio_threshold_would_false_positive():
    """Guard-rail: documents the bug, so nobody 'simplifies' back into it."""
    naive_losses = [(s, n) for s, r in read_giab_baseline() if s != "HG001"
                    for n in SCORED_SEGMENTS if r[n] < 0.8]
    assert len(naive_losses) >= 6, \
        "the raw-ratio trap must remain demonstrable from the committed data"


# --------------------------------------------------------------------------- #
# THAL1 / THAL2 — the only real positive and negative controls that exist
# --------------------------------------------------------------------------- #
def test_thal1_is_a_het_loss_sparing_hbz():
    segs = load_segments()
    calls = {n: hba_depth.call_segment(segs[n], THAL1_RATIOS[n])
             for n in SCORED_SEGMENTS}
    assert calls["HBA2"][1] == "het_loss"
    assert calls["INTER_A2_A1"][1] == "het_loss"
    assert calls["HBA1"][1] == "het_loss"
    # HBZ spared -> --SEA|--MED, not --FIL|--THAI. This is the discrimination.
    assert calls["HBZ"][1] == "intact"
    assert calls["HBA2"][0] == pytest.approx(0.49, abs=0.01)
    assert calls["INTER_A2_A1"][0] == pytest.approx(0.50, abs=0.01)
    assert calls["HBA1"][0] == pytest.approx(0.41, abs=0.01)
    assert calls["HBZ"][0] == pytest.approx(1.13, abs=0.01)


def test_thal2_has_no_deletion():
    segs = load_segments()
    for name in SCORED_SEGMENTS:
        score, call = hba_depth.call_segment(segs[name], THAL2_RATIOS[name])
        assert call == "intact", \
            "{} called {} (score {:.3f}) in a sample with no deletion".format(
                name, call, score)
    # THAL2's raw HBZ is 0.71 — the exact value a global 0.8 cutoff misreads
    assert THAL2_RATIOS["HBZ"] < 0.8
    assert hba_depth.call_segment(segs["HBZ"], 0.71)[0] == pytest.approx(0.93, abs=0.01)


def test_thal_calls_are_not_marginal():
    segs = load_segments()
    for ratios in (THAL1_RATIOS, THAL2_RATIOS):
        for name in SCORED_SEGMENTS:
            score, _ = hba_depth.call_segment(segs[name], ratios[name])
            assert not hba_depth.is_marginal(score)


# --------------------------------------------------------------------------- #
# INTER_Z_A must never produce a numeric verdict
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("ratio", [0.0, 0.5, 0.99, 1.01, 1.9])
def test_inter_z_a_never_gets_a_numeric_call(ratio):
    """It reads 0.99 ('intact') in THAL1 while a --SEA deletion covers half of it."""
    segs = load_segments()
    score, call = hba_depth.call_segment(segs["INTER_Z_A"], ratio)
    assert score is None
    assert call == "do_not_average"


def test_do_not_average_outranks_a_baseline_if_one_is_ever_added():
    seg = {"name": "INTER_Z_A", "baseline": 1.0,
           "reliability": hba_depth.REL_DO_NOT_AVERAGE}
    assert hba_depth.call_segment(seg, 0.5) == (None, "do_not_average")


def test_no_baseline_is_uncalibrated_not_a_call():
    seg = {"name": "NEW", "baseline": None, "reliability": hba_depth.REL_NO_BASELINE}
    assert hba_depth.call_segment(seg, 0.5) == (None, "uncalibrated")
    seg2 = {"name": "NEW", "baseline": None, "reliability": "good"}
    assert hba_depth.call_segment(seg2, 0.5) == (None, "uncalibrated")


# --------------------------------------------------------------------------- #
# End-to-end through main()
# --------------------------------------------------------------------------- #
def _run(tmp_path, segment_means, control_depth=30.0, sample="S"):
    depth = write_depth(tmp_path, segment_means, control_depth=control_depth)
    out = tmp_path / "{}.alpha_depth.tsv".format(sample)
    hba_depth.main(["--depth", str(depth), "--segments", str(SEGMENTS_BED),
                    "--sample", sample, "--out", str(out)])
    lines = out.read_text().splitlines()
    header = lines[0].lstrip("#").split("\t")
    return lines, [dict(zip(header, ln.split("\t"))) for ln in lines[1:]]


def test_contract_schema_header(tmp_path):
    lines, _ = _run(tmp_path, {"HBA2": 30.0})
    assert lines[0] == ("#sample\tsegment\tchrom\tstart\tend\tmean_depth\t"
                        "control_depth\tratio\tbaseline\tscore\tcall\treliability\t"
                        "marginal")


def test_marginal_reaches_the_file(tmp_path):
    """`marginal` was computed into the row dict but omitted from COLUMNS, so it
    never reached the TSV and bin/alpha_globin.py — which downgrades
    alpha_genes_confidence on it — always saw a confident call."""
    _, rows = _run(tmp_path, {"HBA2": 30.0})
    assert all("marginal" in r for r in rows)
    assert {r["marginal"] for r in rows} <= {"True", "False", "NA"}


def test_end_to_end_thal1_like_sample(tmp_path):
    # ctrl 30 -> the THAL1 raw ratios, as measured on the real BAM
    means = {"HBZ": 0.86 * 30, "INTER_Z_A": 0.99 * 30, "HBA2": 0.37 * 30,
             "INTER_A2_A1": 0.50 * 30, "HBA1": 0.40 * 30}
    lines, rows = _run(tmp_path, means, control_depth=30.0, sample="THAL1")
    assert len(rows) == 5
    by = {r["segment"]: r for r in rows}

    assert by["HBA2"]["sample"] == "THAL1"
    assert by["HBA2"]["chrom"] == "chr16"
    assert by["HBA2"]["start"] == "172875" and by["HBA2"]["end"] == "173710"
    assert float(by["HBA2"]["control_depth"]) == pytest.approx(30.0)
    assert float(by["HBA2"]["ratio"]) == pytest.approx(0.370, abs=1e-3)
    assert float(by["HBA2"]["baseline"]) == pytest.approx(0.750)
    assert float(by["HBA2"]["score"]) == pytest.approx(0.493, abs=1e-3)
    assert by["HBA2"]["call"] == "het_loss"
    assert by["HBA2"]["reliability"] == "needs_own_baseline"

    assert by["INTER_A2_A1"]["call"] == "het_loss"
    assert by["HBA1"]["call"] == "het_loss"
    assert by["HBZ"]["call"] == "intact"

    # the trap: measured ratio present, verdict withheld
    assert by["INTER_Z_A"]["call"] == "do_not_average"
    assert by["INTER_Z_A"]["score"] == "NA"
    assert by["INTER_Z_A"]["baseline"] == "NA"
    assert float(by["INTER_Z_A"]["ratio"]) == pytest.approx(0.990, abs=1e-3)

    # every value rounds to 3 dp
    assert by["HBA2"]["ratio"] == "0.370" and by["HBA2"]["score"] == "0.493"


def test_end_to_end_normal_sample_has_no_losses(tmp_path):
    """A sample sitting exactly on the GIAB baselines must be all-intact."""
    means = {"HBZ": 0.760 * 30, "INTER_Z_A": 1.0 * 30, "HBA2": 0.750 * 30,
             "INTER_A2_A1": 1.001 * 30, "HBA1": 0.964 * 30}
    _, rows = _run(tmp_path, means, control_depth=30.0, sample="NORM")
    calls = {r["segment"]: r["call"] for r in rows}
    assert calls == {"HBZ": "intact", "INTER_Z_A": "do_not_average",
                     "HBA2": "intact", "INTER_A2_A1": "intact", "HBA1": "intact"}
    for r in rows:
        if r["score"] != "NA":
            assert float(r["score"]) == pytest.approx(1.0, abs=0.01)


def test_end_to_end_homozygous_and_gain(tmp_path):
    means = {"HBA2": 0.0, "INTER_A2_A1": 1.5 * 1.001 * 30, "HBA1": 0.964 * 30}
    _, rows = _run(tmp_path, means, control_depth=30.0)
    by = {r["segment"]: r for r in rows}
    assert by["HBA2"]["call"] == "hom_loss"
    assert by["INTER_A2_A1"]["call"] == "gain"     # anti-3.7-like, score ~1.5
    assert by["HBA1"]["call"] == "intact"


def test_missing_segment_and_missing_controls(tmp_path):
    # segment absent from the depth file -> no measurement, no verdict
    _, rows = _run(tmp_path, {"HBA2": 30.0})
    by = {r["segment"]: r for r in rows}
    assert by["HBA1"]["mean_depth"] == "NA"
    assert by["HBA1"]["ratio"] == "NA"
    assert by["HBA1"]["call"] == "uncalibrated"

    # no CTRL_* windows at all -> nothing is callable
    depth = write_depth(tmp_path, {"HBA2": 30.0}, n_ctrl=0)
    out = tmp_path / "noctrl.tsv"
    hba_depth.main(["--depth", str(depth), "--segments", str(SEGMENTS_BED),
                    "--sample", "S", "--out", str(out)])
    body = out.read_text().splitlines()[1:]
    for line in body:
        f = line.split("\t")
        assert f[6] == "NA"                        # control_depth
        assert f[10] in ("uncalibrated", "do_not_average")


def test_module_import_has_no_side_effects(tmp_path):
    """Importing must not read the filesystem or write anything."""
    import importlib
    before = sorted(p.name for p in tmp_path.iterdir())
    importlib.reload(hba_depth)
    assert sorted(p.name for p in tmp_path.iterdir()) == before
    assert hba_depth.HET_LOSS_MAX_SCORE == 0.65
    assert hba_depth.HOM_LOSS_MAX_SCORE == 0.25
    assert hba_depth.GAIN_MIN_SCORE == 1.35


def test_reuses_cnv_traits_common_helpers(tmp_path):
    """Contract decision 2: alpha-globin reuses the depth helpers, not a copy."""
    depth = write_depth(tmp_path, {"HBA2": 22.5}, control_depth=30.0)
    depths = c.read_region_depths(depth)
    assert c.control_baseline(depths) == 30.0
    assert depths["HBA2"] == pytest.approx(22.5)
