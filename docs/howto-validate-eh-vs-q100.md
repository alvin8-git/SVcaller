# ExpansionHunter vs GIAB HG002 assembly truth

Why: ExpansionHunter (EH) on HG002 reports ATXN2 22/62 (CI 52-102) and CNBP 13/51 (CI 46-90). The long
alleles rest only on in-repeat and flanking reads (no spanning read) and exceed the 150 bp read length.
The GIAB HG002 assembly truth shows both are normal: ATXN2 carries a 21 bp insertion on one haplotype
(about 22/29), CNBP a 20 bp deletion; no expansion. Alleles that cannot be sized by short reads are
therefore a standing method limitation (see README, "Sizing rule"), not a per-sample result.

Script: `validation/eh_vs_q100.py <str_profile.json> <truth.vcf.gz> <benchmark.bed> <out.tsv>`.
Outputs: `validation/HG002_eh_vs_T2TQ100.tsv` (T2T-Q100 v1.0 stvar) and `validation/HG002_eh_vs_v5.0q.tsv`
(v5.0q stvar). Input `results/HG002/HG002.str_profile.json` (re-run with `--sex male`, see below), 32 loci.

## Method

Truth allele per haplotype = EH reference repeat count (INFO REF of `HG002.str.vcf.gz`; region length /
motif length if that VCF is absent) + net indel length inside the EH ReferenceRegion on that haplotype
(indels clipped to the region) / motif length, rounded. "no truth" = region not wholly inside the
benchmark bed. Agreement = EH alleles within +-1 repeat of the truth alleles (best pairing): both | one |
none. EH point estimates are used even for unresolved alleles.

## Summary (32 loci)

| Truth file | both alleles agree | one allele | none | no truth | unresolved by sizing rule | of those, truth non-pathogenic |
|---|---|---|---|---|---|---|
| T2T-Q100 v1.0 stvar | 11 | 8 | 12 | 1 | 3 (ATXN2, CNBP, RFC1) | 3/3 |
| v5.0q stvar | 15 | 7 | 7 | 3 | 3 (ATXN2, CNBP, RFC1) | 3/3 |

Unresolved-and-truth-normal loci: ATXN2 (truth 22/29 and 21/29), CNBP (truth 34/39 and 34/36, both far below
the 75-repeat pathogenic minimum), RFC1 (truth ~115-117 repeats of the reference AAAAG motif; the pathogenic
motif is AAGGG and the length-only comparison cannot speak to motif, so read as "long reference-type allele").
None of the three is an expansion. EH reported ATXN2 62 and CNBP 51; the truth maximum is 29-30 and 36-39.

## Per-locus table

EH columns from the HG002 json (sizing = alleleN sizing from `bin/eh_sizing.py`); truth columns are repeat counts.

| Locus | EH GT | EH CI | EH sizing | EH reported as | T2T-Q100 truth | agreement | v5.0q truth | agreement |
|---|---|---|---|---|---|---|---|---|
| FMR1 | 31 | 31-31 | spanned | 31 | 20/31 | both | 31/31 | both |
| AFF2 | 14 | 14-14 | spanned | 14 | 15/15 | both | 14/14 | both |
| AR | 28 | 28-28 | spanned | 28 | 23/23 | none | 27/27 | both |
| ATN1 | 13/15 | 13-13/15-15 | spanned/spanned | 13/15 | 14/14 | both | 12/14 | both |
| ATXN1 | 5/19 | 5-5/19-19 | spanned/spanned | 5/19 | 17/17 | none | 17/17 | none |
| ATXN2 | 22/62 | 20-22/52-102 | flanked/unresolved | 22/>=52 | 22/29 | one | 21/29 | one |
| ATXN3 | 20/23 | 20-20/23-23 | spanned/spanned | 20/23 | 19/23 | both | 19/23 | both |
| ATXN7 | 10/13 | 10-10/13-14 | spanned/spanned | 10/13 | 9/9 | one | 9/9 | one |
| ATXN8OS | 4/16 | 4-4/16-16 | spanned/spanned | 4/16 | 18/18 | none | 18/18 | none |
| ATXN10 | 12/16 | 12-12/16-16 | spanned/spanned | 12/16 | 14/14 | none | 12/16 | both |
| CACNA1A | 1/5 | 1-1/5-5 | spanned/spanned | 1/5 | 6/6 | one | 6/6 | one |
| CNBP | 13/51 | 13-13/46-90 | spanned/unresolved | 13/>=46 | 34/39 | none | 34/36 | none |
| C9ORF72 | 2/2 | 2-2/2-2 | spanned/spanned | 2/2 | 2/2 | both | 2/2 | both |
| DMPK | 11/11 | 11-11/11-11 | spanned/spanned | 11/11 | 10/10 | both | 10/10 | both |
| FXN | 10/11 | 10-10/11-11 | spanned/spanned | 10/11 | 6/6 | none | 6/6 | none |
| HTT | 17/24 | 17-17/24-24 | spanned/spanned | 17/24 | 18/18 | one | 16/23 | both |
| NIPA1 | 4/18 | 3-7/18-18 | spanned/spanned | 4/18 | 23/23 | none | no_truth | no_truth |
| NOTCH2NLC | 20/20 | 20-20/20-20 | spanned/spanned | 20/20 | no_truth | no_truth | no_truth | no_truth |
| NOP56 | 2/5 | 2-2/5-5 | spanned/spanned | 2/5 | 2/2 | one | 1/2 | one |
| PABPN1 | 7/7 | 7-7/7-7 | spanned/spanned | 7/7 | 6/6 | both | 6/6 | both |
| PPP2R2B | 4/12 | 4-4/12-12 | spanned/spanned | 4/12 | 12/12 | one | 12/12 | one |
| RFC1 | 12/33 | 12-12/31-53 | spanned/unresolved | 12/>=31 | 11/117 | one | 12/115 | one |
| TBP | 2/14 | 2-2/14-14 | spanned/spanned | 2/14 | 12/12 | none | 12/12 | none |
| BEAN1 | 9/18 | 9-9/18-18 | spanned/spanned | 9/18 | 11/17 | one | 8/17 | both |
| EIF4A3 | 12/12 | 12-12/12-12 | spanned/spanned | 12/12 | 12/12 | both | 12/12 | both |
| MARCHF6 | 15/15 | 15-15/15-16 | spanned/spanned | 15/15 | 14/14 | both | 14/14 | both |
| SAMD12 | 6/6 | 6-6/6-6 | spanned/spanned | 6/6 | 6/6 | both | 6/6 | both |
| RILPL1 | 1/2 | 0-1/0-2 | spanned/spanned | 1/2 | 6/6 | none | 6/6 | none |
| STARD7 | 4/4 | 4-4/4-4 | spanned/spanned | 4/4 | 7/7 | none | 7/7 | none |
| CSTB | 2/2 | 2-2/2-2 | spanned/spanned | 2/2 | 2/2 | both | 1/1 | both |
| FGF14 | 25/25 | 22-25/25-34 | flanked/flanked | 25/25 | 523/525 | none | no_truth | no_truth |
| DAB1 | 2/14 | 0-2/14-14 | flanked/spanned | 2/14 | 17/17 | none | 14/16 | one |

## Sex fix (2026-11-07)

The first HG002 run had no `--sex`, so EH defaulted to female and called chrX loci diploid for a male
(FMR1 31/31, AR 27/28, AFF2 14/14). Sex is now inferred from depth (`bin/infer_sex.py`: HG002 male,
chrX/auto 0.51, chrY/auto 0.56; HG001 female, 0.99 / 0.01) and passed to EH. After the re-run HG002
reads FMR1 31, AR 28, AFF2 14 (hemizygous); v5.0q truth 31, 27, 14, so FMR1 and AFF2 agree and AR is within 1
repeat. No other locus changed, and HG001 (female control) is identical. Old files are kept as `*.pre-sex`.

## Caveats

- The two truth files disagree at loci where the v1.0 stvar set lacks small indels present in v5.0q (for
  example AR: +12 bp in v5.0q only; FMR1: +33 bp on both haplotypes in v5.0q, on one in v1.0). v5.0q is the
  more complete comparator; v1.0 agreement counts are lower for that reason, not because EH changed.
- Interrupted repeats (ATXN1 CAT, ATXN8OS, ATXN10, FXN, TBP, PPP2R2B, RILPL1, STARD7, CACNA1A, NIPA1) disagree
  with the length-derived truth because reference repeat counts and EH's pure-motif counts differ and a
  length model cannot see interruptions. Treat "none/one" there as a limitation of this comparison; a
  sequence-aware comparison (motif counting on the assembly) is not done.
- FGF14 (1.7 kb region, mixed motif) and NIPA1/NOTCH2NLC: truth is unreliable or outside the benchmark.
- Not tested: HG001, HG003-7 (no assembly truth); the rule only needs the json, so it applies to them unchanged.
