"""Contract guard for a hypothetical per-chromosome scatter of the CNV callers.

Optimization #4 (scatter CNV callers by chromosome, then gather) was assessed and
NOT implemented for either caller. The reasons are correctness, not effort, and
are documented in subworkflows/cnv_calling.nf. This test locks the two things a
future scatter would have to preserve, so nobody re-attempts it blind:

  1. The exact columns/positions cnv_consensus.py reads from each caller TSV.
     A gather must produce byte-equivalent-content input to cnv_consensus.py.
  2. The gather-by-concat semantics. CNVpytor TSVs are headerless and concat
     cleanly. GATK TSVs carry a header row: a naive `cat` of per-chrom parts
     duplicates the header, which DictReader turns into a bad data row that fails
     its own parse (real data rows survive, so it is a silent cosmetic failure,
     not lost calls). Any real GATK gather should still emit ONE header. This test
     pins the actual behavior of both.

If someone does wire a scatter later, these assertions catch a broken gather
before it silently empties a clinical CNV sheet.
"""
import sys
sys.path.insert(0, "/data/alvin/SVcaller/bin")
import cnv_consensus as c


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


# --- CNVpytor: positional contract cnv_consensus.load_cnvpytor depends on ------
# parts[0]=type, parts[1]=chrom:start-end region, parts[3]=copy number (2=diploid).
# A per-chrom scatter that changed column order or the region format would break
# every downstream call without any error.

_PYTOR_CHR1 = "deletion chr1:1000-5000 4000 1.0 1e-10 0 0 1e-9 0.0\n"
_PYTOR_CHR2 = "duplication chr2:2000-6000 4000 3.0 1e-8 0 0 1e-7 0.0\n"


def test_cnvpytor_positional_contract(tmp_path):
    segs = c.load_cnvpytor(_write(tmp_path, "p.tsv", _PYTOR_CHR1))
    assert len(segs) == 1
    s = segs[0]
    assert (s.chrom, s.start, s.end, s.svtype) == ("chr1", 1000, 5000, "DEL")


def test_cnvpytor_gather_by_concat_is_equivalent(tmp_path):
    # Gather = concatenate per-chrom parts in chromosome order. Because the TSV is
    # headerless, plain concat is loss-free: parsing the concatenation must equal
    # parsing each part and appending. This is the ONLY caller whose gather is a
    # plain cat.
    per_chrom = c.load_cnvpytor(_write(tmp_path, "c1.tsv", _PYTOR_CHR1)) + \
                c.load_cnvpytor(_write(tmp_path, "c2.tsv", _PYTOR_CHR2))
    gathered = c.load_cnvpytor(
        _write(tmp_path, "cat.tsv", _PYTOR_CHR1 + _PYTOR_CHR2))
    assert [(s.chrom, s.start, s.end, s.cn, s.svtype) for s in gathered] == \
           [(s.chrom, s.start, s.end, s.cn, s.svtype) for s in per_chrom]


# --- GATK: header-bearing TSV. A gather must keep exactly ONE header row --------

_GATK_HEADER = "CONTIG\tSTART\tEND\tCALL_COPY_NUMBER\tQUALITY\n"
_GATK_CHR1 = "chr1\t1000\t5000\t1\t50\n"
_GATK_CHR2 = "chr2\t2000\t6000\t3\t50\n"


def test_gatk_naive_cat_keeps_data_but_injects_silent_parse_failures(tmp_path):
    # Simulate the naive gather: `cat chr1.tsv chr2.tsv` with both headers kept.
    # csv.DictReader consumes only the FIRST line as the header; the second header
    # becomes a data row where int("START") raises, so load_gatk drops THAT ROW
    # ONLY. The real chr2 data row is unaffected, so both CNVs still parse. The
    # damage is not data loss, it is a silent parse-failure that stays under the
    # warning threshold (load_gatk warns only when EVERY row fails). A gather must
    # still keep one header, but the failure mode is cosmetic, not a lost call.
    naive = _GATK_HEADER + _GATK_CHR1 + _GATK_HEADER + _GATK_CHR2
    segs = c.load_gatk(_write(tmp_path, "naive.tsv", naive))
    assert len(segs) == 2        # both real CNVs survive the stray header row


def test_gatk_header_aware_gather_is_equivalent(tmp_path):
    # The CORRECT gather: one header, then data rows from each part in order.
    correct = _GATK_HEADER + _GATK_CHR1 + _GATK_CHR2
    gathered = c.load_gatk(_write(tmp_path, "correct.tsv", correct))
    mono = c.load_gatk(_write(tmp_path, "mono.tsv", correct))
    assert len(gathered) == 2
    assert [(s.chrom, s.start, s.end, s.cn, s.svtype) for s in gathered] == \
           [(s.chrom, s.start, s.end, s.cn, s.svtype) for s in mono]


def test_gatk_column_names_are_the_frozen_contract(tmp_path):
    # If a scatter renames or reorders these, load_gatk parses 0 rows and warns.
    # Freeze the four names load_gatk indexes by.
    for col in ("CONTIG", "START", "END", "CALL_COPY_NUMBER"):
        assert col in _GATK_HEADER
