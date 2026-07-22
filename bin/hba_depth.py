#!/usr/bin/env python3
"""alpha-globin channel 1 — per-segment copy-number scoring from read depth.

Reads a mosdepth `--by` regions file that contains BOTH the alpha diagnostic
segments of `assets/hba_segments.bed` (HBZ, INTER_Z_A, HBA2, INTER_A2_A1, HBA1)
and the `CTRL_*` windows of `assets/cnv_trait_regions.bed`, and emits one scored
row per alpha segment.

    ratio = segment_depth / control_depth
    score = ratio / baseline        <-- every threshold in this file is on SCORE

THRESHOLDING THE RAW RATIO IS A BUG, NOT A SIMPLIFICATION. Intact depth is not
1.0 across most of this locus: intact HBA2 sits at ratio 0.750 and intact HBZ at
0.760, so `ratio < 0.8 => loss` calls a heterozygous deletion in all six GIAB
normals. The per-segment intact baselines live in col 5 of `hba_segments.bed`
and the raw calibration behind them in `validation/giab_alpha_baseline.tsv`.

Two segments are deliberately never given a numeric verdict:

  * `INTER_Z_A` is `do_not_average`. Mapping inflation over chr16:155000-162000
    (1 kb windows read 1.4-1.9) cancels out the real loss over 164000-172875, so
    the 18 kb mean reads 0.99 -- "intact" -- in THAL1, a sample whose --SEA
    deletion covers half the segment. Its mean is not wrong-by-noise, it is
    meaningless. It is reported with its measured ratio and `call=do_not_average`.
  * anything whose col-6 reliability is `no_baseline_yet` (or whose col-5
    baseline is `NA`) is reported `uncalibrated`. No baseline, no verdict.

This module measures. It does not name alleles and does not interpret; that is
channels 2-4 and OmniGen respectively.

Output: <SAMPLE>.alpha_depth.tsv
    #sample segment chrom start end mean_depth control_depth ratio baseline
    score call reliability
"""
import argparse

import cnv_traits_common as c

# --------------------------------------------------------------------------- #
# Thresholds — derived from committed data, not chosen for roundness.
#
# Calibration inputs (all re-derived, none taken on trust):
#
#   INTACT, n=24 = 6 GIAB samples x 4 segments, scores = raw ratio from
#   validation/giab_alpha_baseline.tsv divided by the col-5 baseline of
#   assets/hba_segments.bed. HG001 is EXCLUDED — it reads low across the whole
#   locus (0.38/0.67/0.83/0.57) with a normal chr2 control, in a NON-CONTIGUOUS
#   pattern no deletion can produce, so it is technical dropout of unknown
#   magnitude, not a normal and not a carrier.
#       per-segment  HBZ  0.826-1.189 (mean 1.019, sd 0.133)
#                    HBA2 0.887-1.055 (mean 0.987, sd 0.069)
#                    INTER_A2_A1 0.973-1.150 (mean 1.028, sd 0.065)
#                    HBA1 0.893-1.173 (mean 1.031, sd 0.115)
#       pooled       0.826-1.189, mean 1.016, sd 0.095
#
#   HET LOSS, n=3 — THAL1 (--SEA het), the ONLY positive control that exists:
#       HBA2 0.493 · INTER_A2_A1 0.505 · HBA1 0.411   (HBZ 1.133 = spared)
#
#   HOM LOSS, n=0.  GAIN, n=0.  Both boundaries below are therefore theoretical
#   on one side; see the notes on each.
#
# HET_LOSS_MAX_SCORE = 0.65
#   The gap to split is [0.500 highest true het loss] .. [0.826 lowest GIAB
#   normal]. 0.65 is the geometric midpoint (sqrt(0.500*0.826) = 0.643) rounded
#   to 2 dp, which balances the two margins instead of favouring one:
#       nearest normal   0.826 / 0.65 = 1.271  -> +27.1% headroom
#       nearest het loss 0.65  / 0.500 = 1.301 -> +30.1% headroom
#   Both margins are real but NOT generous. The binding constraint on the normal
#   side is HG004's HBZ (0.826); HBZ is by far the noisiest segment (sd 0.133,
#   vs 0.065 for INTER_A2_A1), so 0.65 sits only ~2.8 sd below the intact HBZ
#   mean, against ~5.5 sd for INTER_A2_A1. Restated: this boundary is safe for
#   HBA2/INTER_A2_A1/HBA1 and merely adequate for HBZ. And no HBZ het loss has
#   ever been observed here (THAL1's --SEA spares HBZ), so the HBZ half of the
#   --SEA|--MED vs --FIL|--THAI discrimination rests on a threshold calibrated
#   on one side only.
#
# HOM_LOSS_MAX_SCORE = 0.25
#   Expected het = 0.5, expected hom = 0.0; 0.25 is equidistant. No homozygous
#   sample exists to check it against, so the only empirical margin is on the
#   het side: THAL1's lowest het loss 0.411 / 0.25 = 1.64x above the cut. The
#   residual mismapping depth of a true homozygous deletion in this paralogous
#   locus is UNMEASURED; if it exceeds 0.25 a hom would be mis-called het.
#
# GAIN_MIN_SCORE = 1.35
#   An anti-3.7 triplication het gives 3 copies where 2 is normal -> 1.5.
#   Geometric midpoint of [1.189 highest GIAB normal] .. [1.5] is 1.336.
#       highest normal 1.35 / 1.189 = 1.135 -> only +13.5% headroom
#       expected gain  1.5  / 1.35  = 1.111 -> only +11.1% headroom
#   THIS ONE IS UNCOMFORTABLY TIGHT AND IS STATED AS SUCH. The intact upper tail
#   (HG005 HBZ 1.189, HG005 HBA1 1.173) reaches four fifths of the way to a real
#   triplication, there is no triplication sample to calibrate against, and the
#   whole separation is ~26% of signal against a pooled intact sd of ~9.5%.
#   Treat every `gain` from this channel as a lead, not a call.
#
# MARGINAL_WINDOW = 0.08
#   ~0.85x the pooled intact sd (0.095), and strictly below the smallest
#   observed calibration-point-to-boundary distance (0.150 = THAL1's
#   INTER_A2_A1 at 0.500 vs the 0.650 cut), so no GIAB normal and no THAL1
#   observation is flagged marginal, while anything genuinely ambiguous is.
#   The integrator uses it to downgrade `alpha_genes_confidence`.
# --------------------------------------------------------------------------- #
HOM_LOSS_MAX_SCORE = 0.25
HET_LOSS_MAX_SCORE = 0.65
GAIN_MIN_SCORE = 1.35
MARGINAL_WINDOW = 0.08

#: boundaries a score can sit uncomfortably close to
BOUNDARIES = (HOM_LOSS_MAX_SCORE, HET_LOSS_MAX_SCORE, GAIN_MIN_SCORE)

# reliability tokens used in col 6 of assets/hba_segments.bed
REL_DO_NOT_AVERAGE = "do_not_average"
REL_NO_BASELINE = "no_baseline_yet"

# call vocabulary
CALL_INTACT = "intact"
CALL_HET_LOSS = "het_loss"
CALL_HOM_LOSS = "hom_loss"
CALL_GAIN = "gain"
CALL_UNCALIBRATED = "uncalibrated"
CALL_DO_NOT_AVERAGE = "do_not_average"

COLUMNS = ["sample", "segment", "chrom", "start", "end", "mean_depth",
           "control_depth", "ratio", "baseline", "score", "call", "reliability",
           # `marginal` is CONSUMED, not decorative: bin/alpha_globin.py reads it
           # to downgrade alpha_genes_confidence to `low`. It was computed into
           # the row dict but left out of this list, so it never reached the
           # file and every marginal call silently read as confident.
           "marginal"]

NA = "NA"


# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #
def parse_segments(path):
    """Parse `assets/hba_segments.bed` into a list of segment dicts.

    Format is BED6-ish with a trailing `#` note on each data line:
        chrom  start  end  name  baseline  reliability  # free-text note

    `#` header lines are skipped. col 5 is `NA` where no baseline can be given.
    """
    segments = []
    with open(str(path)) as fh:
        for line in fh:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            body, _, note = line.rstrip("\n").partition("#")
            f = body.rstrip().split("\t")
            if len(f) < 4:
                continue
            try:
                start, end = int(f[1]), int(f[2])
            except ValueError:
                continue
            raw_baseline = f[4].strip() if len(f) > 4 else ""
            try:
                baseline = float(raw_baseline)
            except ValueError:
                baseline = None
            segments.append({
                "chrom": f[0],
                "start": start,
                "end": end,
                "name": f[3],
                "baseline": baseline,
                "reliability": (f[5].strip() if len(f) > 5 else "") or REL_NO_BASELINE,
                "note": note.strip(),
            })
    return segments


# --------------------------------------------------------------------------- #
# Scoring
# --------------------------------------------------------------------------- #
def score_segment(ratio, baseline):
    """score = ratio / baseline. None when either input is unusable.

    A None score is not "0" and must never be thresholded — it means the segment
    has no intact-depth calibration, so no verdict is available.
    """
    if ratio is None or baseline is None or baseline <= 0:
        return None
    return ratio / baseline


def classify(score):
    """Map a score onto the call vocabulary. Thresholds are on SCORE, never ratio."""
    if score is None:
        return CALL_UNCALIBRATED
    if score < HOM_LOSS_MAX_SCORE:
        return CALL_HOM_LOSS
    if score < HET_LOSS_MAX_SCORE:
        return CALL_HET_LOSS
    if score < GAIN_MIN_SCORE:
        return CALL_INTACT
    return CALL_GAIN


def is_marginal(score):
    """True when the score sits within MARGINAL_WINDOW of a decision boundary.

    Not a call of its own — the call still stands. It exists so the integrator
    can downgrade `alpha_genes_confidence` rather than present a coin-flip as a
    measurement.
    """
    if score is None:
        return False
    return any(abs(score - b) < MARGINAL_WINDOW for b in BOUNDARIES)


def call_segment(segment, ratio):
    """(score, call) for one segment. Honours the reliability gates first.

    Order matters: `do_not_average` outranks everything, because that segment's
    mean is misleading rather than merely uncertain, and a caller that reaches
    the numeric path for it will confidently emit `intact` over a real deletion.
    """
    if segment.get("reliability") == REL_DO_NOT_AVERAGE:
        return None, CALL_DO_NOT_AVERAGE
    if segment.get("baseline") is None or segment.get("reliability") == REL_NO_BASELINE:
        return None, CALL_UNCALIBRATED
    score = score_segment(ratio, segment["baseline"])
    if score is None:
        return None, CALL_UNCALIBRATED
    return score, classify(score)


def evaluate(sample, segments, depths, control_depth):
    """Score every segment. Returns a list of row dicts keyed by COLUMNS."""
    rows = []
    for seg in segments:
        mean_depth = depths.get(seg["name"])
        ratio = (None if mean_depth is None
                 else c.depth_ratio(mean_depth, control_depth))
        score, call = call_segment(seg, ratio)
        rows.append({
            "sample": sample,
            "segment": seg["name"],
            "chrom": seg["chrom"],
            "start": seg["start"],
            "end": seg["end"],
            "mean_depth": mean_depth,
            "control_depth": control_depth,
            "ratio": ratio,
            "baseline": seg["baseline"],
            "score": score,
            "call": call,
            "reliability": seg["reliability"],
            "marginal": is_marginal(score),
        })
    return rows


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #
def _fmt(value):
    if value is None:
        return NA
    if isinstance(value, float):
        return "{:.3f}".format(value)
    return str(value)


def write_rows(out_path, rows):
    """Write the channel-1 detail TSV: '#'-prefixed header then one row/segment."""
    with open(str(out_path), "w") as fh:
        fh.write("#" + "\t".join(COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(_fmt(r[col]) for col in COLUMNS) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="alpha-globin per-segment depth scoring (channel 1)")
    ap.add_argument("--depth", required=True,
                    help="mosdepth --by regions.bed.gz covering the alpha "
                         "segments AND the CTRL_* control windows")
    ap.add_argument("--segments", required=True, help="assets/hba_segments.bed")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True, help="<SAMPLE>.alpha_depth.tsv")
    args = ap.parse_args(argv)

    depths = c.read_region_depths(args.depth)
    segments = parse_segments(args.segments)
    try:
        control_depth = c.control_baseline(depths)
    except ValueError:
        control_depth = None

    rows = evaluate(args.sample, segments, depths, control_depth)
    write_rows(args.out, rows)
    print("wrote {}: {}".format(
        args.out, " ".join("{}={}".format(r["segment"], r["call"]) for r in rows)))
    return rows


if __name__ == "__main__":
    main()
