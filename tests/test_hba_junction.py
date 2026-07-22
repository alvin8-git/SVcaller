"""Tests for channel 3 (HBA_JUNCTION) — bin/hba_junction.py.

Two halves:

1. Pure-function unit tests (no samtools, no BAM). These run everywhere.
2. End-to-end tests against the committed synthetic fixture
   `tests/fixtures/alpha_junction.bam`, plus a synthetic HOMOZYGOUS-deletion
   BAM built in tmp_path. These need samtools on PATH and are skipped without
   it, exactly as tests/test_junction_fixture.py does.

No containers, no network, no reference genome.

The hom-del control is not decoration. Without it, "calls the fixture het"
is satisfiable by a caller that returns the string "het" unconditionally.
"""
import os
import subprocess
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))

import hba_junction as hj  # noqa: E402

SCRIPT = os.path.join(REPO, "bin", "hba_junction.py")
BAM = os.path.join(REPO, "tests", "fixtures", "alpha_junction.bam")

CHROM = "chr16"
DEL_START, DEL_END = 165001, 185000
LEFT_BP, RIGHT_BP = DEL_START - 1, DEL_END + 1        # 165000 -> 185001
SIZE = DEL_END - DEL_START + 1                        # 20000

HAVE_SAMTOOLS = subprocess.run(["which", "samtools"],
                               capture_output=True).returncode == 0
needs_samtools = pytest.mark.skipif(not HAVE_SAMTOOLS,
                                    reason="samtools not on PATH")


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def run_cli(tmp_path, bam, region, name="S", extra=()):
    """Invoke the CLI the way Nextflow would. Returns (proc, parsed_rows)."""
    out = os.path.join(str(tmp_path), f"{name}.alpha_junction.tsv")
    proc = subprocess.run(
        [sys.executable, SCRIPT, "--bam", bam, "--sample", name,
         "--out", out, "--region", region, *extra],
        capture_output=True, text=True)
    rows = []
    header = None
    if os.path.exists(out):
        with open(out) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                if line.startswith("#"):
                    header = line
                else:
                    rows.append(dict(zip(hj.COLUMNS, line.split("\t"))))
    return proc, header, rows


# ---------------------------------------------------------------------------
# ref_span — the CIGAR trap
# ---------------------------------------------------------------------------

def test_ref_span_plain_match():
    assert hj.ref_span("150M") == 150


def test_ref_span_ignores_soft_and_hard_clips():
    # Clips consume READ bases, never REFERENCE bases.
    assert hj.ref_span("50S100M") == 100
    assert hj.ref_span("100M50S") == 100
    assert hj.ref_span("30H20S100M") == 100


def test_ref_span_insertions_do_not_consume_reference():
    assert hj.ref_span("50M10I90M") == 140


def test_ref_span_deletions_and_skips_do_consume_reference():
    assert hj.ref_span("50M10D90M") == 150
    assert hj.ref_span("50M1000N100M") == 1150


def test_ref_span_beats_the_naive_digit_sum():
    """A naive `sum of every number in the CIGAR` is the bug this guards.

    On these CIGARs the naive answer is wrong by tens of bases, which would
    place the breakend off by the same amount and blow the base-resolution
    requirement.
    """
    import re
    for cigar, truth in [("10S100M20I30D5S", 130),
                         ("5H20M100N30M10S", 150),
                         ("40S110M", 110),
                         ("70S80M", 80)]:
        naive = sum(int(n) for n in re.findall(r"\d+", cigar))
        assert hj.ref_span(cigar) == truth
        assert naive != truth, f"{cigar}: naive sum coincidentally correct"


def test_ref_span_equals_and_x_ops():
    assert hj.ref_span("50=10X90=") == 150


def test_ref_span_rejects_garbage():
    with pytest.raises(ValueError):
        hj.ref_span("100Z")
    assert hj.ref_span("*") == 0
    assert hj.ref_span("") == 0


def test_aln_end_is_last_aligned_base_inclusive():
    # pos 1 with 150M covers 1..150, not 1..151.
    assert hj.aln_end(1, "150M") == 150
    assert hj.aln_end(164891, "110M40S") == LEFT_BP


# ---------------------------------------------------------------------------
# clip_ends
# ---------------------------------------------------------------------------

def test_clip_ends_trailing_clip_gives_the_left_breakend():
    lead, trail = hj.clip_ends(164891, "110M40S", min_clip=15)
    assert lead is None
    assert trail == LEFT_BP


def test_clip_ends_leading_clip_gives_the_right_breakend():
    lead, trail = hj.clip_ends(RIGHT_BP, "40S110M", min_clip=15)
    assert lead == RIGHT_BP
    assert trail is None


def test_clip_ends_unclipped_read_votes_for_nothing():
    assert hj.clip_ends(164880, "150M", min_clip=15) == (None, None)


def test_clip_ends_respects_min_clip():
    # A 5 bp clip is adapter/terminal-mismatch noise, not a junction.
    assert hj.clip_ends(1000, "145M5S", min_clip=15) == (None, None)
    assert hj.clip_ends(1000, "145M5S", min_clip=5)[1] == 1144


def test_clip_ends_both_ends_clipped():
    lead, trail = hj.clip_ends(1000, "20S100M30S", min_clip=15)
    assert (lead, trail) == (1000, 1099)


def test_clip_ends_tolerates_hard_clips_around_the_soft_clip():
    lead, trail = hj.clip_ends(2000, "10H30S100M", min_clip=15)
    assert lead == 2000 and trail is None


def test_clip_ends_uses_ref_span_not_read_length():
    """The breakend must be placed with a reference-aware span."""
    # 50M10D50M40S: reference span 110, so the left breakend is pos+109.
    lead, trail = hj.clip_ends(1000, "50M10D50M40S", min_clip=15)
    assert trail == 1109
    # naive read-length arithmetic (100 aligned bases) would say 1099
    assert trail != 1099


def test_clip_ends_fully_clipped_read_yields_nothing():
    assert hj.clip_ends(1000, "150S") == (None, None)


# ---------------------------------------------------------------------------
# cluster
# ---------------------------------------------------------------------------

def test_cluster_groups_within_window_and_reports_the_mode():
    # 5 reads clip at 165000, jitter of 1-2 bp on two others.
    pos = [165000] * 5 + [165001, 164999]
    assert hj.cluster(pos, 10) == [(165000, 7)]


def test_cluster_representative_is_the_mode_not_the_mean():
    """Mean would give a non-integral, off-base answer."""
    pos = [100, 100, 100, 108]
    (rep, n), = hj.cluster(pos, 10)
    assert rep == 100 and n == 4
    assert rep != sum(pos) / len(pos)


def test_cluster_splits_beyond_the_window():
    assert hj.cluster([100, 101, 500, 501, 502], 10) == [(100, 2), (500, 3)]


def test_cluster_ties_break_leftmost_for_determinism():
    (rep, n), = hj.cluster([100, 100, 105, 105], 10)
    assert (rep, n) == (100, 4)


def test_cluster_empty():
    assert hj.cluster([], 10) == []


def test_cluster_does_not_merge_distinct_alpha_alleles():
    """-a3.7 and --SEA breakends are kilobases apart; window=10 must not
    collapse them into one junction."""
    out = hj.cluster([165000] * 5 + [168700] * 5, hj.DEFAULT_WINDOW)
    assert [p for p, _ in out] == [165000, 168700]


# ---------------------------------------------------------------------------
# zygosity_from_balance — the reason this module exists
# ---------------------------------------------------------------------------

def test_zygosity_balanced_is_het():
    vaf, call = hj.zygosity_from_balance(12, 10)
    assert call == "het"
    assert vaf == 0.545


def test_zygosity_no_spanning_reads_is_hom():
    vaf, call = hj.zygosity_from_balance(12, 0)
    assert (vaf, call) == (1.0, "hom")


def test_zygosity_one_stray_spanning_read_is_still_hom():
    # A single mismapped reference read must not flip a hom-del to het.
    assert hj.zygosity_from_balance(30, 1)[1] == "hom"


def test_zygosity_low_fraction_is_refused_not_guessed():
    """3 split against 40 spanning is ~7% allele fraction. That is not a
    germline heterozygote; say NA rather than invent a genotype."""
    vaf, call = hj.zygosity_from_balance(3, 40)
    assert call == "NA" and vaf < 0.15


def test_zygosity_no_evidence_at_all():
    assert hj.zygosity_from_balance(0, 0) == (None, "NA")


def test_zygosity_rounds_to_three_dp():
    vaf, _ = hj.zygosity_from_balance(1, 2)
    assert vaf == 0.333


# ---------------------------------------------------------------------------
# insert_cutoff / parse_region / confidence
# ---------------------------------------------------------------------------

def test_insert_cutoff_is_derived_from_the_sample_not_hardcoded():
    # median 450 -> 3x = 1350, above the 1000 floor
    assert hj.insert_cutoff([450] * 18) == 1350
    # a long-insert library scales up
    assert hj.insert_cutoff([900] * 10) == 2700


def test_insert_cutoff_floor_protects_against_degenerate_estimates():
    assert hj.insert_cutoff([10, 20, 30]) == hj.MIN_INSERT_CUTOFF
    assert hj.insert_cutoff([]) == 1500       # fallback median 500


def test_parse_region():
    assert hj.parse_region("chr16:1-250000") == ("chr16", 1, 250000)
    assert hj.parse_region("chr16") == ("chr16", 1, None)
    assert hj.parse_region("chr16:150,000-200,000") == ("chr16", 150000, 200000)
    with pytest.raises(ValueError):
        hj.parse_region("chr16:bad-region")


def test_confidence_tiers():
    assert hj.confidence(6, 6, 6, 10, "het") == "high"
    assert hj.confidence(6, 6, 1, 10, "het") == "medium"    # thin discordant
    assert hj.confidence(6, 6, 6, 10, "NA") == "low"
    assert hj.confidence(1, 6, 6, 10, "het") == "low"       # one side too thin


def test_split_only_clip_stack_is_never_high_or_medium():
    """The empirical artifact signature at this locus: plenty of soft clips on
    both sides, zero discordant pairs. THAL2 — which has NO deletion — produces
    exactly this at chr16:~182000 and ~188365. It must not outrank a real
    junction, and it must never be used to collapse a degenerate allele group.
    """
    assert hj.confidence(11, 6, 0, 22, "het") == "low"
    assert hj.confidence(50, 50, 0, 50, "het") == "low"


# ---------------------------------------------------------------------------
# contract / schema
# ---------------------------------------------------------------------------

def test_output_header_is_exactly_the_agreed_schema():
    assert hj.HEADER == (
        "#sample\tchrom\tleft_bp\tright_bp\tsize\tsplit_left\tsplit_right\t"
        "discordant\tspanning\tvaf\tzygosity\tconfidence")


def test_format_row_computes_size_as_right_minus_left_minus_one():
    row = hj.format_row("S", CHROM, LEFT_BP, RIGHT_BP, 6, 6, 6, 10, 0.545,
                        "het", "high").split("\t")
    assert len(row) == len(hj.COLUMNS)
    assert row[hj.COLUMNS.index("size")] == str(SIZE)
    assert row[hj.COLUMNS.index("vaf")] == "0.545"


def test_format_row_renders_missing_vaf_as_NA():
    row = hj.format_row("S", CHROM, 1, 10, 3, 3, 0, 0, None, "NA", "low")
    assert row.split("\t")[hj.COLUMNS.index("vaf")] == "NA"


def test_module_import_has_no_side_effects():
    """Importing must not read a BAM, write a file, or parse argv."""
    import importlib
    importlib.reload(hj)


# ---------------------------------------------------------------------------
# END-TO-END on the committed fixture
# ---------------------------------------------------------------------------

@needs_samtools
def test_fixture_present():
    assert os.path.exists(BAM), (
        "tests/fixtures/alpha_junction.bam is missing. If git dropped it, "
        "check the *.bam rule in .gitignore still has the tests/fixtures/ "
        "exception.")


@needs_samtools
def test_fixture_breakpoint_recovered_to_the_base(tmp_path):
    """THE ACCEPTANCE CRITERION. 165000 -> 185001, exactly, from the data."""
    proc, header, rows = run_cli(tmp_path, BAM, "chr16:1-250000",
                                 name="JUNCTION_FIXTURE")
    assert proc.returncode == 0, proc.stderr
    assert header is not None
    assert len(rows) == 1, f"expected exactly one junction, got {rows}"
    r = rows[0]
    assert r["chrom"] == CHROM
    assert int(r["left_bp"]) == LEFT_BP
    assert int(r["right_bp"]) == RIGHT_BP
    assert int(r["size"]) == SIZE


@needs_samtools
def test_fixture_is_called_heterozygous(tmp_path):
    """The fixture deliberately carries reference-spanning reads. A caller
    that ignores allele balance reports hom-del here and is WRONG."""
    _, _, rows = run_cli(tmp_path, BAM, "chr16:1-250000")
    assert rows[0]["zygosity"] == "het", (
        f"called {rows[0]['zygosity']} — the intact allele's spanning reads "
        "were ignored")
    assert 0.3 < float(rows[0]["vaf"]) < 0.7


@needs_samtools
def test_fixture_support_counts(tmp_path):
    _, _, rows = run_cli(tmp_path, BAM, "chr16:1-250000")
    r = rows[0]
    assert int(r["split_left"]) == 6
    assert int(r["split_right"]) == 6
    # 12 discordant READS = 6 discordant FRAGMENTS; a pair is one observation.
    assert int(r["discordant"]) == 6
    assert int(r["spanning"]) == 10
    assert r["confidence"] == "high"


@needs_samtools
def test_fixture_sample_id_is_propagated(tmp_path):
    _, _, rows = run_cli(tmp_path, BAM, "chr16:1-250000", name="THAL1")
    assert rows[0]["sample"] == "THAL1"


@needs_samtools
def test_min_split_gate_suppresses_the_call(tmp_path):
    """Raising the gate above the available support must yield a clean
    negative, not a crash."""
    proc, header, rows = run_cli(tmp_path, BAM, "chr16:1-250000",
                                 extra=["--min-split", "20"])
    assert proc.returncode == 0
    assert header is not None and rows == []


@needs_samtools
def test_clean_region_produces_a_header_only_file_and_exit_zero(tmp_path):
    """chr16:189000-200000 of the fixture holds only unclipped proper pairs.
    'No junction' is a RESULT, not an error: the file must exist with its
    header and the process must exit 0, so a Nextflow process does not fail
    on every normal sample."""
    proc, header, rows = run_cli(tmp_path, BAM, "chr16:189000-200000")
    assert proc.returncode == 0, proc.stderr
    out = os.path.join(str(tmp_path), "S.alpha_junction.tsv")
    assert os.path.exists(out) and os.path.getsize(out) > 0
    assert header == hj.HEADER
    assert rows == []


# ---------------------------------------------------------------------------
# HOMOZYGOUS-DELETION NEGATIVE CONTROL
# ---------------------------------------------------------------------------

HOM_LEFT, HOM_RIGHT = 100000, 110001
READ_LEN = 150


def _sam_record(qname, flag, pos, cigar, pnext=0, tlen=0, mapq=60):
    seq = "A" * READ_LEN
    return "\t".join([qname, str(flag), CHROM, str(pos), str(mapq), cigar,
                      "=" if pnext else "*", str(pnext), str(tlen),
                      seq, "I" * READ_LEN])


def _build_hom_del_bam(tmp_path):
    """A synthetic HOMOZYGOUS deletion: split reads on both sides, discordant
    pairs across it, and NOT ONE reference-spanning read.

    This is what makes the het logic non-vacuous. If hba_junction.py hardcoded
    "het", or derived zygosity from split count alone, this test fails.
    """
    lines = ["@HD\tVN:1.6\tSO:coordinate",
             f"@SQ\tSN:{CHROM}\tLN:90338345",
             "@CO\tSYNTHETIC hom-del control built by tests/test_hba_junction.py"]
    recs = []
    for i in range(6):
        clip = 40 + i * 5
        m = READ_LEN - clip
        # left splits: aligned block ENDS at HOM_LEFT
        recs.append(_sam_record(f"hsl_{i}", 99, HOM_LEFT - m + 1, f"{m}M{clip}S"))
        # right splits: aligned block STARTS at HOM_RIGHT
        recs.append(_sam_record(f"hsr_{i}", 147, HOM_RIGHT, f"{clip}S{m}M"))
    for i in range(5):
        p1 = HOM_LEFT - 400 - i * 30
        p2 = HOM_RIGHT + 100 + i * 30
        tlen = p2 + READ_LEN - p1
        recs.append(_sam_record(f"hd_{i}", 97, p1, "150M", pnext=p2, tlen=tlen))
        recs.append(_sam_record(f"hd_{i}", 145, p2, "150M", pnext=p1, tlen=-tlen))
    # far-away normal pairs, so an insert-size median can be estimated at all
    for i in range(4):
        p1 = HOM_LEFT - 20000 - i * 500
        p2 = p1 + 300
        recs.append(_sam_record(f"hn_{i}", 99, p1, "150M", pnext=p2, tlen=450))
        recs.append(_sam_record(f"hn_{i}", 147, p2, "150M", pnext=p1, tlen=-450))

    recs.sort(key=lambda r: int(r.split("\t")[3]))
    sam = os.path.join(str(tmp_path), "hom.sam")
    with open(sam, "w") as fh:
        fh.write("\n".join(lines + recs) + "\n")
    bam = os.path.join(str(tmp_path), "hom.bam")
    subprocess.run(["samtools", "sort", "-o", bam, sam], check=True,
                   capture_output=True)
    subprocess.run(["samtools", "index", bam], check=True, capture_output=True)
    return bam


@needs_samtools
def test_hom_del_control_is_called_homozygous(tmp_path):
    bam = _build_hom_del_bam(tmp_path)
    proc, _, rows = run_cli(tmp_path, bam, "chr16:1-250000", name="HOMCTRL")
    assert proc.returncode == 0, proc.stderr
    assert len(rows) == 1, f"expected one junction, got {rows}"
    r = rows[0]
    assert int(r["left_bp"]) == HOM_LEFT
    assert int(r["right_bp"]) == HOM_RIGHT
    assert int(r["spanning"]) == 0, "control must contain no spanning reads"
    assert r["zygosity"] == "hom", (
        f"called {r['zygosity']} on a BAM with zero reference-spanning reads; "
        "the het logic is not actually reading allele balance")
    assert float(r["vaf"]) == 1.0


@needs_samtools
def test_het_and_hom_controls_disagree(tmp_path):
    """Same code path, opposite answers — the discriminating property."""
    _, _, het_rows = run_cli(tmp_path, BAM, "chr16:1-250000", name="HET")
    hom_bam = _build_hom_del_bam(tmp_path)
    _, _, hom_rows = run_cli(tmp_path, hom_bam, "chr16:1-250000", name="HOM")
    assert het_rows[0]["zygosity"] == "het"
    assert hom_rows[0]["zygosity"] == "hom"
    assert float(het_rows[0]["vaf"]) < float(hom_rows[0]["vaf"])


@needs_samtools
def test_supplementary_and_duplicate_reads_are_not_counted(tmp_path):
    """Supplementary alignments are the other half of a split read; counting
    them would double the apparent support from one fragment."""
    bam = _build_hom_del_bam(tmp_path)
    recs = hj.read_sam_records(bam, "chr16:1-250000")
    assert recs
    for r in recs:
        assert not r["flag"] & hj.FLAG_EXCLUDE


@needs_samtools
def test_discordant_pair_must_match_the_junction_size(tmp_path):
    """A pair straddling both breakends is not enough — its insert must also be
    consistent with the deletion's size. One chimeric pair with an absurd TLEN
    otherwise "supports" every candidate junction in the region (THAL1 has such
    a pair, TLEN ~50 Mb)."""
    recs = [
        # straddles 1000..2001 with a sane insert -> supports a 1000 bp deletion
        {"qname": "good", "flag": 0x1, "chrom": CHROM, "pos": 900,
         "cigar": "150M", "rnext": "=", "pnext": 2100, "tlen": 1350,
         "end": 1049},
        # straddles the same breakends but implies a 50 Mb fragment
        {"qname": "chimera", "flag": 0x1, "chrom": CHROM, "pos": 900,
         "cigar": "150M", "rnext": "=", "pnext": 50_000_000,
         "tlen": 50_000_000, "end": 1049},
    ]
    frags = hj.discordant_fragments(recs, CHROM, cutoff=1000)
    assert set(frags) == {"good", "chimera"}, "both are discordant candidates"
    size = 2001 - 1000 - 1
    supporting = [q for q, (s, e) in frags.items()
                  if s <= 1000 and e >= 2001 and (e - s + 1) - size <= 1000]
    assert supporting == ["good"]


@needs_samtools
def test_window_zero_still_recovers_the_exact_breakpoint(tmp_path):
    """The fixture's clips are base-identical, so clustering is not what is
    producing the right answer."""
    _, _, rows = run_cli(tmp_path, BAM, "chr16:1-250000", extra=["--window", "0"])
    assert int(rows[0]["left_bp"]) == LEFT_BP
    assert int(rows[0]["right_bp"]) == RIGHT_BP
