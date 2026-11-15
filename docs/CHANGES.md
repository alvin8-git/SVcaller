# Changes

## 2026-11-15 - code-review fixes on the uncommitted batch (behaviour otherwise unchanged)

HTML report: (1) the STR confidence interval follows the same cited offset as the count (`str_summary` `ci`; the VCF fallback shifts REPCI with `_offset_repci`). (2) `eh_count_offsets.tsv`, `str_retired_loci.tsv`, `str_count_convention.tsv`, `str_not_tested_legacy.tsv` are staged `BUILD_HTML_REPORT` inputs (`--str-offsets/--str-retired/--str-convention/--str-not-tested`), so `-resume` re-renders when one changes. (4) A legacy-catalogue `str_summary.tsv` with no EH VCF now lists its loci (it was an empty section, which reads as "no expansions"). (5) An empty SMN1 count shows SMA as "NOT TYPED" in the known-disease tier (was skipped silently). (9) `str_grading_scope.load_scope()` reads the three tables and the narrow-gap set once per report; `not_graded_reason(locus, scope)` has no module-level cache. `hba_report.py` (6): a depth-only suspected loss makes the Deletion row read "Suspected (read depth only), not confirmed". `eh_sizing.py` (7): an offset with no EH CI keeps the empty CI text instead of writing `None-None`. `nf-run-headless.sh` (8): the default log directory is `$PWD/.nf-headless`, never `/tmp`. Pipeline (3): MELT's `-c` has its own `--melt_coverage` (default 25, the value it ran with before `min_depth` went to 20; `conf/test.config` 5), so lowering `min_depth` does not change MELT. HG002 report re-rendered from stored server outputs before and after: byte-identical (no HG002 result touches these paths).

## 2026-11-15 - mean-depth floor lowered from 25x to 20x (30x stays recommended)

Owner decision. `min_depth` default 25 -> 20 (`nextflow.config`; MELT's `-c` is separate, `melt_coverage` 25). Evidence: truvari 4.3.1 on HG002 raw merged calls vs GIAB HG002
T2T-Q100 v1.0 stvar, downsampled 30/25/22.5/20/17/14/12x (one sample, one replicate, raw calls). Deletion recall 1-10 kb .76 at 30x and 25x, .74 at
20x, .73 at 17x; 10 kb or more .45 -> .41 at 20x; 300 bp-1 kb .62 at 30x, .58 at 20x; no cliff at 25x. Alpha-globin high confidence down to 20x, medium
at 17x; SMN no-call below 18x by design. The abort message is unchanged in form (`ERROR: mean depth N < 20x (required). Aborting.`). The report's QC
card adds, for 20-30x mean depth, "mean depth below the recommended 30x: structural-variant sensitivity is reduced (about 2-4 points fewer deletions
found at 20x than at 30x in GIAB HG002)". README and Documentation.md state 20x (Documentation.md previously said 30 in two places).

## 2026-11-15 - the HTML report follows every owner decision about SVcaller results (never contradicts OmniGen)

**Problem.** `bin/html_report.py` graded PABPN1 from the VCF's raw REPCN (11/11, "INTERMEDIATE") instead of the offset-corrected
`str_summary.tsv` (10/10, normal); it graded retired loci, loci without a cited counting convention and loci of the old catalogue whose
coordinates hold no repeat in GRCh38 (NIPA1 showed EXPANDED); it showed ">=N" alleles as numbers/intermediate; the alpha-globin card
rendered a depth-only suspected loss as "alpha+ deletion, biallelic"; a delly inversion time cap was not mentioned.

**Fix.** (1) STR: `html_report.py` takes `--str-summary` (wired in `subworkflows/report.nf`; without it the VCF count has the cited offset
applied). `bin/str_grading_scope.py` reads three carried curation tables with provenance headers (`assets/str_retired_loci.tsv`,
`assets/str_count_convention.tsv`, `assets/str_not_tested_legacy.tsv`; `tests/test_str_curation_scope.py` compares them with OmniGen's
`panel_sources` when `OMNIGEN_PANEL_SOURCES` is set). A locus set aside shows NOT GRADED with its reason and is never a Finding; an
allele known only as ">=N" shows UNRESOLVED (at least N repeats), ">=0" shows NOT MEASURABLE. (2) `--sv-not-assessed` (published
`<id>.sv_not_assessed.tsv`) shows "Partly screened: delly inversions not assessed"; no record reads "not recorded". (3) `hba_report.py`:
a depth-only suspected -a3.7/-a4.2 reads "suspected, not confirmed; confirm by gap-PCR or MLPA", NA reads "Not determined". (4) The
Clinical Findings card says SV/CNV rows are caller output, not graded findings (the OmniGen Finding gate needs ClinGen dosage files and a
cohort recurrence list that are not carried here); CNV-syndrome rows say the same. (5) An empty SMN1 cell is never a carrier.
Unchanged: SMN no-call below 18x (already honoured by `smn_report.py` and the xlsx), the 25x depth gate (the run aborts, so no report exists).

## 2026-11-15 - tracked files scrubbed of internal paths and private ids

Internal absolute paths (`/data/...`, `/home/...`) replaced: tests derive the repo root from `__file__` and read the reference tree from `REF_ROOT` (default `/path/to/ref`; reference-dependent tests skip when it is absent); validation shell scripts use `${REF_ROOT}`, `${SVCALLER_HOME}`, `${TMPDIR}`; `conf/local.config` reads `SVCALLER_CONDA_ENV`; docs use `/path/to/...` placeholders. A private clinical sample id in this file became "sample A". Git history still contains the old text and is not rewritten.

## 2026-11-15 - PABPN1 repeat count corrected by -1 before grading

**Problem.** STRchive's OPMD_PABPN1 range (normal max 10, intermediate 11, pathogenic min 12) counts GCN codons between ATG and GGG. GRCh38 holds 10 codons there (ATG, (GCG)6(GCA)3(GCG)1, GGG GCT GCG) inside the 31 bp interval chr14:23321472-23321503. ExpansionHunter's degenerate GCN motif absorbs the GGG (one mismatch) as an 11th unit, so the common 10-codon allele is reported as 11 and graded intermediate. In HG002 all 21 spanning reads carry exactly 10 codons while EH calls 11/11; 50 of 64 cohort samples are 11/11.

**Fix (owner decision).** `assets/eh_count_offsets.tsv` (locus, variant_id, offset, evidence) lists PABPN1 as -1; `bin/eh_sizing.py` adds it to genotype, ci, alleleN_est and the displays in `str_summary.tsv`. New trailing columns `eh_raw_genotype` (EH's own value) and `count_offset`. Sizing (spanned/flanked/unresolved) is still decided on EH's raw counts; the spanning/flanking read-count columns stay in EH units. `EH_SIZING` takes the file as an input. No other locus has an offset (NIPA1 checked: reads equal the reference; every other locus has a normal-to-pathogenic gap of 3+ units).

**Not covered.** `bin/html_report.py` grades the STR VCF `REPCN` (EH raw) and does not apply the offset; `str_profile.json` and the VCF stay as EH wrote them. Tests: `tests/test_eh_sizing.py` (reference holds 10 codons, offsets file is PABPN1 -1 only, 11/11 -> 10/10, 13 -> 12 stays pathogenic range). Existing output needs the sizing step re-run (needs `<id>.str_profile.json`, `<id>.sex.json`).

## 2026-11-15 - population-AF filter reads the column of the SV's own type; provenance.json on failed runs

**Problem 1 (audit items 4 and 5).** `GNOMAD_SV_FILTER` used `row.get('B_gain_AFmax', row.get('B_loss_AFmax','0'))`. AnnotSV fills all four `B_gain/B_loss/B_ins/B_inv` column groups for every SV (`AnnotSV-benignsv.tcl`, `benignSVannotation`: benign records with AF >= `-benignAF` 0.01 that completely contain the SV, `bedtools -f 1.0`), so the gain column was read for every type. On HG002 (3.4.6 output, 44,895 rows): 2,725 of 5,466 full-mode DELs were dropped because a common *gain* region contained them (a PASS Delly DEL inside a large gain region vanished), while common *deletions* were kept whenever no gain contained them. INS, INV and TRA had the same flaw.

**Fix.** The filter moved to `bin/sv_af_filter.py` with a type map: DEL->`B_loss_AFmax`, DUP->`B_gain_AFmax`, INS->`B_ins_AFmax`, INV->`B_inv_AFmax`; other types (TRA/BND) have no AF column and are kept. HG002 full rows kept, before -> after: DEL 2741 -> 708 (298 restored, 2331 common deletions now removed), DUP 707 -> 707 (unchanged), INS 943 -> 214, INV 617 -> 234, TRA 1181 -> 2707 (all kept). The other reason a merged-VCF DEL has no row is intended: AnnotSV runs with `-SVminSize 50`, so calls under 50 bp never reach the table (4,759 DELs on HG002). Tests: `tests/test_sv_af_filter.py`.

**Problem 2.** `provenance.json not written: arraycopy: element type mismatch` (commit bdf3625): the `onComplete` command list held GStrings, which `ProcessBuilder(List<String>)` rejects, so no run ever wrote it. After that, on a failed run Nextflow interrupts the `onComplete` thread during shutdown, killing the 8 s writer.

**Fix.** `main.nf`: `write_provenance()` builds a `List<String>`, runs at workflow start and again in `onComplete`.

## 2026-11-15 - repeat-expansion catalog rebuilt from STRchive hg38; compound loci graded on the disease block

**Problem (gap audit 2026-11-15).** `assets/eh_catalog.json` was typed by hand. Four entries used hg19 coordinates (ATXN1, ATXN8OS, PPP2R2B, SAMD12), five matched neither build (CACNA1A, NIPA1, EIF4A3, RILPL1, STARD7), six regions were too narrow (AFF2, AR, ATN1, ATXN3, C9ORF72, TBP), FGF14 too wide, and eight compound repeats were modelled as one pure motif (CNBP, HTT, ATXN7, NOP56, MARCHF6, FXN, RFC1, DAB1), so a count was in the wrong unit (CNBP: whole tract / CCTG; RFC1: AAAAG instead of the pathogenic AAGGG).

**Fix.** The catalog is now generated, not edited: `bin/make_eh_catalog.py` reads the pinned STRchive hg38 catalogs in `assets/sources/` (release v2.28.0, commit `e964fa88df3b2189ee478cbd6eaad15bb28e3a0a`, `assets/sources/PROVENANCE.md`) and `assets/eh_locus_map.tsv` (our LocusId -> STRchive id; the 32 LocusIds are unchanged) and writes `assets/eh_catalog.json` plus `assets/eh_graded_variants.tsv`. Compound loci use STRchive's multi-block `LocusStructure` with its `VariantId`s (ExpansionHunter catalog syntax: Illumina/ExpansionHunter `docs/04_VariantCatalogFiles.md`). 11 loci are multi-block (ATXN7, ATXN8OS, CNBP, FXN, HTT, NOP56, RFC1, MARCHF6, SAMD12, STARD7, DAB1); ExpansionHunter now reports 44 variants for the 32 loci. `eh_graded_variants.tsv` names, per locus, the variant that is the disease block (STRchive `PathologicRegion`), its motif in disease orientation (`DisplayRU`) and STRchive's NormalMax / PathologicMin.

**`bin/eh_sizing.py`.** Reads the graded-variant table (`--graded`, staged by `modules/expansionhunter/sizing.nf`) and writes ONE row per locus, for the graded variant, in repeat units of that block: CNBP is the CAGG (= CCTG) count, RFC1 the AAGGG count. The other blocks stay in `<id>.str_profile.json`. `str_summary.tsv` gains two columns at the END, `variant_id` and `disease_motif`; every existing column keeps its name, position and meaning, so readers that parse by header name are unaffected. A multi-block locus without a graded entry, or a profile that lacks the graded variant (made with an older catalog), is refused with exit 1 instead of guessed.

**Consumers that read `str_profile.json` directly** now see several Variants for the 11 compound loci under `LocusResults.<LocusId>.Variants`; only the graded one (`eh_graded_variants.tsv`) is the disease repeat. Thresholds: see the owner-decision paragraph below.

**Proof.** ExpansionHunter 5.0.0 container (`conf/docker.config`), GIAB HG002 BAM, old and new catalog, then `eh_sizing.py` (also inside `svcaller/utils:1.3`), compared with HG002 v5.0q / T2T-Q100 truth. Of the nine wrong-position loci, six now give measured, spanned calls (ATXN1 31/32, ATXN8OS 16/17, CACNA1A 12/12, PPP2R2B 11/11, NIPA1 9/9, RILPL1 10/11); SAMD12 and STARD7 report 0 copies of the pathogenic block with no reads (shown as `>=0`, unresolved) and EIF4A3's 20-bp motif only bounds the count (`>=7/>=8`); CNBP 13/>=46 became 15/16 (normal); ATXN2 22/>=52 became 23/31 spanned; TBP 2/14 became 37/38.  Tests: `tests/test_eh_catalog.py` (every entry equals its STRchive source, source bytes pinned by sha256, generator output equals the committed files), `tests/test_eh_sizing.py` (compound locus).

**Owner decisions applied (2026-11-15).** (A) AR and AFF2 regions come from the Illumina ExpansionHunter catalog (`variant_catalog/grch38/variant_catalog.json`, commit `7ff9aa5335fdb5b3124cc0f365a64bb7a185803a`, carried in `assets/sources/`), marked as an explicit override in `assets/eh_locus_map.tsv` and checked by `tests/test_eh_catalog.py`: STRchive's 103/149 bp regions made EH count the region length (HG002 AR 40, AFF2 50 flanked). With the override HG002 gives AR 27 (truth 23/27, hemizygous male) and AFF2 19 (low-confidence no longer; truth by region length 19/20, pure CCG run 14/15). (B) `assets/str_disease_ranges.tsv` is now generated from STRchive v2.28.0 (`bin/make_str_ranges.py`, source `assets/sources/STRchive-loci.json`, `tests/test_str_ranges.py`) in repeat units of the graded block: `normal_max` = STRchive benign_max, `path_min` = pathogenic_min, new columns `intermediate_min/max`, `graded_variant`, `strchive_id`, `source` (release, locus page, first references). Existing column names are kept (`premut_max` is blank). BEAN1, MARCHF6, SAMD12, STARD7 have no STRchive normal maximum and the intermediate range is absent for 8 loci; those cells are blank and the `notes` column says so. Disease labels now follow STRchive (RILPL1 = OPDM4, STARD7 = FAME2; FAME4 is YEATS2).

## 2026-11-15 - alpha-globin: HBZ is informational (15 of 64 samples were "not resolved")

**Problem.** HBZ (chr16:152643-154503) read `het_loss` in 12 of 58 independent samples, and in 11 of the 15 `NA` samples it was the only reason
(HBA2, INTER_A2_A1, HBA1 intact, or a -a4.2 het pattern). A het loss of HBZ beside intact genes matches
no allele, so the call was `NA`. This was the other side of 4384ce3, which turned the same inputs from a false
3-gene carrier into `NA`.

**Evidence.** (1) No panel allele removes HBZ alone: --FIL and --THAI take HBZ together with HBA2, INTER_A2_A1 and
HBA1; --SEA and --MED lie inside the cluster and spare HBZ; -a3.7 and -a4.2 are crossovers between the Z and X boxes
of the duplicated alpha genes, downstream of HBZ (Higgs 2013, Cold Spring Harb Perspect Med 3:a011718,
doi 10.1101/cshperspect.a011718, PMC3530043; allele extents in `assets/hba_deletion_alleles.tsv`, basis `literature`).
HBZ therefore matters only to split --SEA|--MED from --FIL|--THAI (South-East Asian alleles). (2) HBZ cannot be
calibrated: across 58 independent samples the ratio is one continuum (min 0.40, p10 0.46, median 0.60, p90 0.85, max
2.5 with a gain; GIAB HG001-HG007 0.47-0.87), against the fixed baseline 0.760 from 6 GIAB samples; a het loss
would read about 0.38, so 12 of 58 normals sit at or below the het-loss cut. HBZ is uncorrelated with the genes
(HG001's HBZ is the lowest GIAB value while its INTER_A2_A1 is intact; the NIST benchmark shows heterozygous
variants across the span).

**Change (`bin/alpha_globin.py`, `_hbz_usable`).** An HBZ `het_loss`, and any marginal HBZ reading, is treated
as not measured; a clear `intact`, `hom_loss` or `gain` still counts (intact excludes --FIL|--THAI; a gained HBZ
vetoes the anti-3.7 triplication, which spares HBZ). A --SEA/--MED carrier whose HBZ reads low is no longer
narrowed to --FIL|--THAI: the reported group is `--SEA|--MED|--FIL|--THAI`. Unchanged: junctions only
corroborate, -a4.2 needs the between-gene drop, marginal and unmeasured segments are not determined, no carrier on
one marginal segment. `assets/hba_segments.bed` is untouched (the depth file still reports HBZ).

**Owner decision (same day): depth-only single-copy -a3.7/-a4.2 is not determined.** Such a call with no
supporting high-confidence junction is emitted as `alpha_genes_called=NA`, confidence `low`, `deletion_alleles=NA`,
`deletion_evidence=depth`, the allele in the new 12th column `suspected_allele` and the segment ratios in the note
("depth pattern suggests -a4.2 (HBA2 x, INTER_A2_A1 y, HBA1 z); no breakpoint junction; confirm by gap-PCR or
MLPA"). Two-copy calls (sample A -a4.2/-a4.2) and --SEA|--MED etc. are unchanged. Replay vs the HBZ-only
version: only three clinical samples changed. OmniGen's `report_alpha_globin.load_contract` reads the new
file unchanged.

**Proof.** `alpha_globin.py` replayed old vs new on the stored depth/junction/sites files of 66 samples (58
independent): `NA` 15 to 5; 10 `NA` samples now 4 genes, 1 now -a4.2/aa (3 genes, medium); 55 unchanged,
including HG002-HG007, THAL1 (--SEA|--MED/aa) and THAL2 (none).

## 2026-11-14 - delly INV step: time cap (`--delly_inv_max_time 2h`), explicit "not assessed" record

**Problem.** On two ~55-60x samples `delly call -t INV` ran 7-8 h in "Split-read assembly" (one thread, no I/O) while DEL/DUP/INS/BND finished in ~80 min; normal samples take ~30 min for INV. Cause of the extra inversion-like split reads not established.

**Fix.** `bin/delly_call_capped.sh` (called by `DELLY_CALL_SVTYPE`) bounds the INV type only by `params.delly_inv_max_time` (default `2h`, `0` = no cap; other types untouched). The cap is a watchdog in the script, not `timeout`: the delly image has BusyBox `timeout`, which returns 143 rather than 124. At the cap the task exits 0 with a WARN in `.command.err`, a header-only `<id>.delly.INV.vcf` (BGZF, sample name from the BAM `SM` tag, contigs from the fai) so `DELLY_MERGE` runs unchanged, and one row in the published `<outdir>/<id>/<id>.sv_not_assessed.tsv` (`sample caller svtype reason limit` = `S1 delly INV time_limit_exceeded 7200s`). The INV task always writes this file: header only = delly INV completed, row = not assessed; file absent = a pipeline version before this change (or a failed task). A delly failure fails the task (no row): the cap counts only when the watchdog fired AND delly died of its TERM/KILL (exit 143/137). The watchdog sends TERM, waits up to 10 s, then KILL. The placeholder's sample column is the SM of the first `@RG` (parsed from the BAM header, `@RG` lines only), else the BAM file name minus extension, which is what delly 1.2.6 itself writes (checked on real delly VCFs, with and without `@RG`), so `DELLY_MERGE`'s concat sees the same sample name. The cap is parsed once by `delly_inv_cap_seconds` (modules/delly/call.nf; bare number = seconds, `0`/`false`/null = none, under 1 s = error at workflow start) and passed to the task as an input value, 0 for every type but INV, so changing the param re-runs only the INV task (and the merge) under `-resume`. Manta (as run here, no `convertInversion`) emits inversions as BND and no INV records; GRIDSS (opt-in) converts BND pairs to INV; so with default callers a capped sample has no INV-typed calls from any caller. Test: `tests/test_delly_inv_cap.py` (fake `delly`); verified by running the real process in the delly 1.2.6 container on a GIAB HG002 slice (normal limit identical to the unmodified process; 1 s cap; forced failure).

## 2026-11-13 — SMN: depth floor, a low-depth result is a no-call (`--smn_min_depth 18`)

**Problem.** A GIAB HG002 downsampling series (truth SMN1=2, SMN2=2) gave a wrong `3/1` flagged `PASS:AgreeWithSum` at 14.7x
(Total_CN_raw 4.06 was right, the split was not); 2/2 returned at 17.9x and above, and `NoCall:InsufficientDepth` only at 9.8x.

**Fix.** New process `SMN_DEPTH_FLOOR` (`bin/smn_depth_floor.py`) after `SMN_CALLER`, fed the mean depth MOSDEPTH already measured. Below
`params.smn_min_depth` (default 18, separate from `min_depth`) `isSMA`/`isCarrier`/`SMN1_CN`/`SMN2_CN` = `None`, `Info` =
`NoCall:BelowValidatedDepth`; the caller's row is kept as `<sample>.smn.raw.tsv`. Unchanged bytes at/above the floor; `NoCall:*` rows pass
through. The 18x value is empirical (one sample per depth), indicative only. No Docker rebuild.

**Follow-up fix (same day).** The process script passed `in/${smn_tsv.name}`, but `stageAs: 'in/*'` already renders the name with its `in/` prefix, so the command read `in/in/<id>.smn.tsv` and every sample that reached `SMN_DEPTH_FLOOR` failed with exit 2. The script now passes `${smn_tsv}`; staging under `in/` is kept. Verified by running the process (above and below the floor) in `svcaller/utils:1.3`; `tests/test_smn_report.py` has a static guard.

## 2026-11-13 — alpha-globin: no more named carriers from one marginal segment, compound genotypes, or noise junctions

**Problem.** About 4 of 45 samples (including GIAB HG001, a European) were reported `--FIL|--THAI/anti-3.7`,
a false carrier call: the only HBZ-loss genotype in the allele table is a two-gene deletion on one
chromosome plus the anti-3.7 triplication on the other, and the whole call rested on HBZ (1.8 kb,
`needs_own_baseline`) scoring 0.62-0.64, just under the het-loss cut and flagged `marginal`. The NIST benchmark shows
balanced heterozygous variants across the supposedly deleted span. Separately, `deletion_evidence` returned
`both`/`junction` whenever ANY junction row existed, and nearly every sample has 1-8 `confidence low` junction rows
in the chr16 subtelomere, so normal samples printed `evidence junction`.

**Fix (`bin/alpha_globin.py`).**
1. *Marginal readings.* Each marginal segment is re-read as its adjacent state (the nearest boundary of `score`,
   or both neighbours when the row has no `score`) and every combination is matched against the allele table.
   (a) a non-wild-type call that turns into wild type under any valid reading is `NA`/`NA` with a note naming the
   segment(s) (not `none`/4: a marginal dip is not evidence of normal); (b) all valid readings non-wild-type but
   with different gene counts: call kept, confidence `medium` (never high; the loss is not in doubt), note gives
   the alternative count; (c) no alternative
   reading gives a different valid genotype: the marginal flag no longer lowers confidence. Confidence is now `low` only for: depth missing, a contradicting
   non-low junction, a called-normal sample whose marginal segment has a valid non-wild-type reading, or `NA` from
   1(a). Junction rows present but none supporting is `medium` (depth only), as when the file is absent
   (previously every marginal segment forced `low`, including clear single-gene carriers, and "channels disagree"
   fired on mere absence of a supporting junction). Note `hba_depth.py`'s `marginal` flag
   only says "near some boundary", hence the use of `score`.
2. *Compound genotype.* A deletion allele on one chromosome plus anti-3.7 on the other is never named from depth
   alone (`NA`/`NA`, note says two-event genotype). **Known blind spot**, alongside the parsimony filter below:
   a real compensated deletion + triplication reads as normal (parsimony) or undetermined (this rule).
3. *Junction channel is corroboration only (`junction_support`).* It can raise confidence (medium to high) and set
   `deletion_evidence=both`; it never lowers confidence, never makes a call `NA`, never contradicts one (only
   two real samples have a high-confidence junction, so a stricter matcher could not be validated). Only
   `confidence == high` rows count (`hba_junction.py`: "use high only"). A row touching no scorable segment is
   ignored. A row SUPPORTS the named genotype if it fully CONTAINS every segment the genotype removes (-1; no
   tolerance) and overlaps no segment it leaves untouched (0 in every named allele; all four for a called-normal sample). Any other high row
   overlapping a scorable segment is only noted ("junction chr16:a-b seen, not used"). `deletion_evidence`: `both`
   (support), `depth`, `junction` (such a high row exists but depth named no deletion: normal or `NA`), `none`.
   The junction file cannot name an allele, so an `NA` stays `NA`. Zygosity and breakpoint tolerances are not used.
4. *Confidence.* `low` needs a positive reason: depth missing or any `NA` genotype. Junction rows present but none
   supporting is `medium` (depth only), as when the file is absent. Marginal readings obey rule 2. `genes` is a
   required argument of `confidence()`. The boundaries and marginal window are copies of `hba_depth.py`'s
   constants, tied by tests. `marginal_alternatives` is judged from the `call` (one state away only): a finite
   score picks the neighbour, a printed score within 0.001 of a boundary that disagrees with the call by one step is
   rounding (trust the call), anything else offers both neighbours. A depth file without a `marginal` column has the
   flag derived from `score` with the same 0.08 window.
5. *Marginal logic, refined.* For a marginal segment called as a loss or gain, an alternative reading toward 2
   copies that gives wild type, or only the banned deletion + anti-3.7 compound, means the call rests on marginal
   evidence: `NA`/`NA` with a note. An alternative that matches no genotype at all caps confidence at `medium`
   with a note (a real two-gene deletion with one borderline segment survives). A marginal segment called intact
   whose alternative matches nothing is ignorable. Valid alternatives differing in gene count or alleles cap at
   `medium` with a note. A called-normal sample with a VALID deletion reading on a marginal segment is `NA`/`NA`,
   never `none`/4. No real sample in the current set is affected by that last rule.
6. *-a4.2 guard (empirical).* HBA2 depth alone dips on some libraries with no deletion (four independent library
   pairs), while every confirmed -a4.2 carrier (three, two confirmed on a second sequencing platform) also shows
   INTER_A2_A1 reduced (ratio 0.47-0.60 het, 0.19-0.21 hom; normals 0.78-1.29). `-a4.2` is named only when
   INTER_A2_A1 is called a loss, or called intact with a marginal flag whose alternative reading is ONE copy
   (toward loss; marginal toward gain does not count); otherwise `NA`/`NA` ("HBA2 depth loss without the expected
   drop between the genes"). The allele table keeps `?`; the guard is a named rule in `bin/alpha_globin.py`.
7. *Unmeasured segments.* If HBA2, INTER_A2_A1 or HBA1 has no usable call the result is `NA`/`NA`, `low`, with a
   note naming it (never `none`/4, never -a4.2 with INTER unmeasured). HBZ unmeasured is tolerated: the matcher
   leaves HBZ unconstrained, so `--SEA|--MED|--FIL|--THAI` is reported as one wider group.
8. *-a4.2 vs -a3.7/-a4.2.* HBA2 and INTER_A2_A1 both het_loss with HBA1 intact keeps `-a4.2/aa` (3 genes) by
   parsimony, with a note that `-a3.7/-a4.2` (2 genes) gives the same depth, and confidence at most `medium`.
9. *Known limit (not fixed).* Homozygous anti-3.7 (6 genes) is not resolved: `gain` maps to 3 copies only, so it is
   emitted as `anti-3.7/aa` (5 genes). Whenever any segment is called `gain` the note says: "gain copy number above
   3 is not resolved; anti-3.7 on both chromosomes (6 genes) reads the same".
10. *-a4.2 with a borderline INTER_A2_A1.* Named at most `medium` (even with a supporting junction) and always with a
   note, whenever INTER_A2_A1 is marginal (marginal loss, or intact marginal toward loss). DECLINED: making a
   marginal het_loss INTER (score ~0.57-0.65) `NA`; confirmed het carriers read INTER 0.47-0.60 (two platforms)
   against 0.78-1.29 in normals, so a marginal loss there is carrier-like, not a coin flip. The guard fires only
   when HBA2 is itself called a loss. A called pure gain (anti-3.7/aa) with one borderline segment whose
   toward-normal reading only fits a compound keeps its call at `medium` with a note; `NA` is reserved for readings
   that reach wild type and for deletion calls resting on marginal segments.

Output columns and the contract are unchanged; `NA` was already a legal value. True positives (--SEA|--MED het,
-a3.7, -a4.2, hom -a4.2-like, anti-3.7) give the same gene counts.

## 2026-11-13 — SMN caller: low-depth input is an explicit no-call, not a crash

On a ~9.8x BAM (mean-depth floor lowered for a depth limit-of-detection series) `SMN_CALLER` died with
`IndexError: list index out of range` in `caller/call_smn12.py` (`get_smn1_call_and_tag`) and ended the pipeline.
Cause: at low depth no per-site call reaches the posterior cutoff, so the `Counter` of site calls is empty and the
upstream code indexes `[0]` (line 137 for the loose cutoff, line 153 for the medium cutoff).

- Fix is a patch to the vendored-by-clone Illumina v1.1 source, `modules/smn_caller/patches/smn_low_depth_nocall.patch`,
  applied in `Dockerfile.smncopynum`; new image `svcaller/smncopynum:1.2` (the `:1.1` tag is untouched). Rebuild with
  `docker build -f Dockerfile.smncopynum -t svcaller/smncopynum:1.2 modules/smn_caller/patches`. Referenced from
  `modules/smn_caller/call.nf` and `conf/docker.config`.
- Only the empty-site-list condition becomes a no-call: `SMN1_CN`/`SMN2_CN`/`isSMA`/`isCarrier` = `None`,
  `Info=NoCall:InsufficientDepth`, exit 0. The TSV gains a trailing `Info` column (existing columns and values are
  unchanged; full-depth output is identical). All other exceptions still propagate, and the `SMN_CALLER` empty/header-only/
  no-JSON guards are unchanged.
- `bin/smn_report.py` renders a no-call (or any row without an SMN1 copy number) as "No call"; before, `None` fell back to
  CN=2 and "Normal".
- Measured on a GRCh38 30x sample (SMN regions only, subsampled): 33x identical to the published call; ~10x `Info=Ambiguous`
  (an existing no-call); ~5x crashed before (line 153) and is now `NoCall:InsufficientDepth`; ~2x crashed before (line 137),
  now `NoCall:InsufficientDepth`.

## 2026-11-11 — Wall-time: scheduling, resource reservations, sharded BAM filter

Goal: cut wall time per 30X sample with the same tools and byte-identical results. Measured on one dev box
(HG002, 115 GB BAM, `--max_cpus 48 --max_memory 100.GB`, Docker, spinning-disk-class shared storage):
**145 min -> 77 min**. All 20 final outputs compared identical (`bin/compare_results.py`, below).

- **`withName` shadowing bug (largest effect).** A `withName` selector in `conf/base.config` whose string is also used
  in `conf/docker.config` (which sets `container` per process) is *replaced*, not merged. The `DELLY_CALL_SVTYPE`
  CPU pin from 2026-07-24 was therefore ignored under `-profile docker` (checked with `nextflow config`): Delly ran with
  `process_medium` (8 CPU / 36 GB), at most two Delly tasks ran at once on a 100 GB budget, and the five types ran
  serially in pairs (~30 min x 3 rounds). Resource overrides now use grouped selectors unique to `base.config`
  (`'DELLY_CALL_SVTYPE|STRLING_CALL'`, `'CNVPYTOR_CALL|GATK_GCNV_CALL'`); all five Delly types now run together.
  `MANTA_CALL`/`SVABA_CALL` blocks in `base.config` are still shadowed under docker (they equal the label values).
- **Right-sized reservations** from measured peak RSS and %CPU (`conf/base.config`): `process_medium` 36 -> 16 GB,
  `process_low` 12 -> 6 GB; Delly 1 CPU / 8 GB; STRling 1 CPU / 8 GB; gCNV and CNVpytor 2 CPU / 16 GB; flagstat and
  mosdepth 4 CPU; Picard insert size 1 CPU. Retries still scale memory (x attempt). No tool argument changed.
- **`SAMTOOLS_FILTER_CHROMS`** takes 3/4 of `--max_cpus` (min 8). `bin/filter_chroms.sh` writes its parts at BGZF
  level 1 (`FILTER_BGZF_LEVEL`, default 1; the compress stage was 67 s vs 29 s for chr21), splits chromosomes
  longer than `FILTER_SHARD_BP` (default 25000000) into position windows (each read is emitted once, original order),
  and indexes with `-@`. Reads and header (minus `@PG`) are identical to the previous script
  (`tests/test_filter_chroms.py::test_sharded_output_equals_unsharded`, plus an md5 of `samtools view` on HG002).
  Level 1 makes the filtered BAM ~8% larger than level 6 (119 GB vs 110 GB for HG002); on slow shared storage
  where five Delly tasks read it together set `FILTER_BGZF_LEVEL=6`.
- Delly BND is submitted first (`subworkflows/sv_calling.nf`); `maxForks 4` removed from `DELLY_CALL_SVTYPE`.
- New `bin/compare_results.py <A> <B> <sample>`: exits 1 if any final output differs (VCF volatile header lines
  such as `##fileDate`/`##reference=` ignored).
- **Evaluated and not changed (results could not be guaranteed or no wall-time gain):**
  Delly sharding by chromosome (BND needs all chromosome pairs and is the longest task, so sharding DEL/DUP/INV/INS
  does not shorten the critical path; sharded IDs would also differ); CNVpytor `-j` (thread output not verified);
  Picard insert size (QC only, feeds the HTML report/MultiQC, no call; runs concurrently from t=0, ~30 min, off the
  critical path); gCNV `CollectReadCounts` sharding (off the critical path).
- Remaining critical path: filter (~5 min of sharded work, then a serial `samtools cat` of ~120 GB and index) ->
  Delly BND (single-threaded, 30-45 min, I/O bound when many readers share the BAM) -> Jasmine/AnnotSV/depth check.

## 2026-11-10 — TRAIT_DEPTH now publishes `<S>.trait_depth.regions.bed.gz`

- `TRAIT_DEPTH` (`modules/traits/depth.nf`) computed the targeted-depth BED but never published it, so no run
  produced the file OmniGen's `report_cnv_traits.py` reads for the `[CNVTRAIT-UNCAL]` AMY1 and LPA KIV-2 lines.
  Added `publishDir "${params.outdir}/${meta.id}", mode: 'copy', pattern: "*.trait_depth.regions.bed.gz"`
  (same pattern as the other per-sample outputs). The mosdepth command and BED contents are unchanged; the `.csi`
  index is not a declared output and is not published.

## 2026-11-07 — ExpansionHunter per-allele sizing (spanned / flanked / unresolved)

- EH reported HG002 ATXN2 22/62 and CNBP 13/51 from in-repeat and flanking reads only (alleles beyond the 150 bp
  read, no spanning read). GIAB HG002 assembly truth shows ATXN2 ~22/29 and CNBP a 20 bp deletion, no expansion.
- New `bin/eh_sizing.py` + process `EH_SIZING` (`modules/expansionhunter/sizing.nf`, wired in
  `subworkflows/sv_calling.nf`, container `svcaller/utils:1.3` in `conf/docker.config`): writes
  `<sample>.str_summary.tsv` with `allele{1,2}_sizing`, `allele{1,2}_display`, `genotype_display`, `conf`.
  Unresolved alleles are `>=<CI lower bound>` and `conf=unresolved`. `str_profile.json` is unchanged.
- ExpansionHunter was run without `--sex` (EH defaults to female): HG002, a male, had chrX loci diploid (FMR1 31/31).
  New `INFER_SEX` process (`bin/infer_sex.py`, depth windows, OmniGen `vcf_qc.py` thresholds) feeds `--sex`; sex, ratios,
  method and any unknown-sex flag go to `<sample>.sex.json` and the `str_summary.tsv` header. Male chrX loci are
  hemizygous (`allele2` = `.`). HG002 re-run: FMR1 31, AR 28, AFF2 14; HG001 (female) unchanged. Old outputs kept as `*.pre-sex`.
- `validation/eh_vs_q100.py` compares EH with Q100 / v5.0q truth; results in `validation/HG002_eh_vs_*.tsv` and
  `docs/howto-validate-eh-vs-q100.md`. Tests: `tests/test_eh_sizing.py`.
- `results/HG001` and `results/HG002` annotated standalone (`str_summary.tsv` is new; no prior summary existed, so no
  `.pre-sizing` copy). Pipeline not re-run.

## 2026-11-07 — Missing local images: svcaller/utils:1.3 and smncopynum:1.1 rebuilt, launch guard added

The HG001_md run died at `ALPHA_GLOBIN:HBA_JUNCTION` (exit 125, "pull access denied for
svcaller/utils") after 42 min. Nothing was wrong with the pin: `svcaller/utils:1.3` is the intended
tag (`nextflow.config`, `conf/docker.config`, ~12 modules, VERSION.md) but the current rootless Docker
daemon (`DOCKER_HOST=unix:///run/user/1001/docker.sock`) had no `svcaller/*` images at all. Built from
the repo's Dockerfiles: `docker build -f Dockerfile.utils -t svcaller/utils:1.3 .` and
`docker build -f Dockerfile.smncopynum -t svcaller/smncopynum:1.1 .` (MELT image not built; opt-in).
No config change. New `tests/test_local_images.py`: one utils tag everywhere, and (when a daemon is
reachable) every pinned `svcaller/*` image must exist locally, failing with the build command.
Public biocontainers (manta, delly, annotsv, ...) are pulled on demand by Nextflow.

## 2026-11-07 — Default callers are Manta + Delly; GRIDSS and MELT opt-in; SUPP_VEC order change

**Downstream parsers (OmniGen) must read the caller order, not assume it.** Owner decision: the
consumer product uses Manta + Delly only. On HG001 read depth refutes 96 % of large DEL/DUP calls,
MELT costs ~8 h per 30x genome with 0 reportable records on the validation samples, GRIDSS has
unbounded run time, and the SMN1/alpha-globin/DMD/HBA tiers already cover the whole-gene deletions
that matter. Expected saving ~8 h (MELT) + 4-6 h (GRIDSS) + ~10 min (Scramble) per genome.

- `nextflow.config`: new `--run_gridss`, `--run_melt` (default false); `skip_gridss`, `skip_melt` now default
  true and still work (`--skip_gridss false` enables). Effective rule: runs iff `run_X || !skip_X`.
  Scramble is also opt-in (`--run_scramble`, `skip_scramble` default true, same rule), so the default SUPP_VEC is always `Manta,Delly` (width 2).
- **SUPP_VEC order** was `Manta,Delly,GRIDSS[,Scramble][,MELT]` (width 3-5). Now one bit per caller that ran:
  `Manta,Delly,[GRIDSS],[Scramble],[MELT]`; default run is `Manta,Delly` (always width 2).
  `bin/caller_set.sh` (new) builds `vcf_list.txt` and writes `<id>.sv_callers.json` (published);
  `JASMINE_MERGE` adds `##SVcaller_callers=` to the VCF header. The GRIDSS-only-TRA filter now uses the real GRIDSS bit.
  `merge.nf` still takes five VCF slots (stubs for callers that did not run); merge tolerances and filters unchanged.
- `filtered.tsv` gains `CALLERS_RAN`, `SUPP_VEC_WIDTH`, `SUPP_CALLERS`, `CALLER_CONSENSUS` (`bin/caller_consensus.py`, in
  `GNOMAD_SV_FILTER`; `DEPTH_*` stay the last seven columns). With only Manta+Delly in the vector, PASS needs both; otherwise
  PASS = >=2 callers as before; single-caller rows are kept as `SINGLE_CALLER`. Old tables lack these columns; a VCF without the header gives `CALLERS_RAN=legacy`, `CALLER_CONSENSUS=NA`.
- `bin/html_report.py` + template name the callers from `CALLERS_RAN` (legacy tables keep the old fallback).
- Tests: `tests/test_caller_set.py` (width-2 vector, provenance list, consensus rule, defaults).

## 2026-11-07 — Read-depth consistency columns on filtered.tsv

HG001 deletions at VHL (chr3, 1.88 Mb), PPARG (chr3, 957 kb) were Manta+Delly-supported but read
1.0x of their flanks. New `DEPTH_CHECK` process (`bin/depth_check.py`, wired in
`subworkflows/annotate.nf` after `GNOMAD_SV_FILTER`) appends `DEPTH_IN`, `DEPTH_FLANK`,
`DEPTH_GLOBAL`, `DEPTH_RATIO`, `DEPTH_RATIO_GLOBAL`, `DEPTH_SUPPORT`, `DEPTH_NOTE` to `<sample>.filtered.tsv` for DEL/DUP >= 1 kb. Annotation only;
no rows dropped, existing columns unchanged. `GNOMAD_SV_FILTER` no longer publishes (DEPTH_CHECK
does). `ANNOTATE` takes the validated BAM as a third input. No new dependency (samtools in the
utils container). Two baselines (local flank and the sample\'s autosomal median): supported needs both, refuted
needs either, so a segdup-inflated flank (HG001 NSF: in 29x, flank 46x, global 31x) is refuted
and noted "flank inflated". Test: `tests/test_depth_check.py`.

## 2026-07-24 — Class A runtime optimisation: right-size DELLY's CPU reservation

Trace-based profiling of the THAL1/THAL2 run showed the wall clock inflated by a
scheduler mismatch, not by any real compute. These are **rigor-preserving** changes:
no caller, threshold, filter, region, or output file changes. VCFs and TSVs are
byte-for-byte identical; only Nextflow's CPU *reservations* move.

**DELLY CPU reservation right-sized (8 → 2).** `DELLY_CALL_SVTYPE` inherited the
shared `process_medium` label (cpus=8) but `delly call` is effectively
single-threaded — it parallelises across the 5 SV *types* (5 separate tasks fanned
out by `subworkflows/sv_calling.nf`), not across threads within a task, and the
script passes no thread flag. On the 64-core THAL host the 5 DELLY tasks reserved
40 cpus which, stacked against Manta (16), Scramble (8) and the CNV callers,
oversubscribed to ~88 cpus and forced DELLY's fan-out to queue serially — its
5-svtype footprint stretched to ~62 min instead of the ~26 min the work actually
takes. Added a `withName: 'DELLY_CALL_SVTYPE'` override in `conf/base.config`
pinning cpus=2 (one worker + headroom), scoped to DELLY *only*. The other
`process_medium` users (samtools sort/subset/filter_chroms, picard markduplicates,
gatk gCNV, cnvpytor, annotsv, strling, scramble, jasmine, melt) are genuinely
multi-threaded and were deliberately left untouched. Because `delly call` receives
no thread argument, its output is unchanged — this frees scheduler slots only.
**Expected saving: ~25–40 min/run.** Validated with `nextflow config` (parses
clean; override resolves to `cpus = { check_max( 2, 'cpus' ) }`).

### Deferred / not changed (with reasons)

**mosdepth pass consolidation — DEFERRED (not safe).** The three full-BAM mosdepth
scans cannot be folded into one pass without changing outputs. Each uses a
*different* `--by` target and mosdepth accepts only one per invocation, emitting one
region-set file per consumer:
- QC `MOSDEPTH` — `--by 50000` (genome-wide 50 kb windows) + `--quantize 0:5:30:500:`
  + the min-depth gate; feeds MultiQC/summary.
- `TRAIT_DEPTH` — `--mapq 0 --by assets/cnv_trait_regions.bed`; feeds `cnv_traits_common.py`.
- `HBA_DEPTH` — `--mapq 0 --by <hba_segments+CTRL bed>`; feeds `hba_depth.py`.

A single pass would require restoring the deliberately-disabled per-base output
(~4.5 GB) and rewriting all three consumer scripts, changing file provenance and the
region-name contracts each script expects — outside the rigor-preserving envelope.
Left as-is.

**Rerun/resume overhead — OPS NOTE, no config change.** ~55 min of the THAL wall was
rerun cost from a transient first-attempt failure, not a config defect.
`errorStrategy` retries only on OOM/kill signals `[143,137,104,134,139]` with
`maxRetries=2` (GRIDSS 3) — a legitimate escalate-memory-on-`task.attempt` policy,
not an over-aggressive blanket retry. `-resume` caching is not defeated by config
(`overwrite=true` is set only on the trace/report/timeline files, which is correct).
No safe config fix exists; operationally, always launch detached and pass an explicit
`-resume <session>` so a transient failure resumes instead of recomputing.

**Heavy-caller guardrail — documented, defaults unchanged.** Current default state
(`conf`/`nextflow.config`, verified 2026-07-24):
- `skip_svaba = true` — SvABA disabled by default (staged but never merged; see the
  2026-07-23 entry). This is the only heavy caller off by default.
- `skip_melt = false` and `skip_gridss = false` — **MELT and GRIDSS RUN by default.**
  The THAL1/THAL2 profiling run must have passed `--skip_melt --skip_gridss` on the
  command line; the pipeline does *not* stub them implicitly. Their skips are already
  explicit params selecting `MELT_STUB`/`GRIDSS_STUB` in `subworkflows/sv_calling.nf`.
  These defaults were **left unchanged** — flipping them to `true` would drop MELT
  MEI calls and GRIDSS SV calls from the standard ensemble, i.e. change outputs, which
  is out of scope for a Class A change. To reproduce the fast THAL profile, pass the
  skip flags at run time; do not bake them into defaults.

## 2026-07-23 — two callers that never ran, the alpha-globin card wired in, and QC made honest

**MELT had never produced a merged call.** `modules/melt/call.nf` passed
`-reads ${params.melt_min_reads}`. MELT v2.2.2 has no such option: it printed its
usage text and exited non-zero for every ME type, on every sample, since the flag
landed on 2026-06-05. The module's own guard correctly refused to publish an empty
VCF as "no mobile element insertions", so the failure was loud, but nobody had read
it. The run16 F1 baseline (2026-06-08) was measured while MELT was failing, so it
was effectively a 4-caller result (Manta + DELLY + GRIDSS + Scramble), not the
5-caller it was labelled. Fixed the flag to `-sr`; MELT now reaches the JASMINE
merge for the first time. `tests/test_melt_flags.py` pins MELT's accepted flag set
so this cannot recur. Also corrected the docs: MELT's own `-sr` default is 0 (no
filtering), so `--melt_min_reads 3` is more stringent than stock MELT, the opposite
of the old "default 5, lowered for recall" claim.

**SvABA correctly labelled and switched off by default.** SvABA has never been
merged (`modules/jasmine/merge.nf` reads `vcfs[0..4]` and never `vcfs[5]`), so it
burned about 4 CPU-hours per sample producing nothing that reached the merge.
`skip_svaba` now defaults to true. The report no longer names SvABA in its caller
list; `tests/test_ensemble_caller_count.py` now scans the report template, not just
the docs, so a stale "+ SvABA" cannot survive there again. Re-enabling needs
`vcfs[5]` wired into the merge first, then a benchmarked HG002 comparison.

**Alpha-globin card wired into the report and redesigned.** `bin/hba_report.py` had
rendered a factual alpha-globin card for a while, but nothing connected it to
`subworkflows/report.nf`, so a `--SEA` carrier's result existed only in
`results/<sample>/alpha_globin/<sample>.alpha_globin.tsv` and was invisible in the
HTML. Wired it using the same pre-render pattern as SMN. Then rebuilt the card
SMN-style: a headline (alpha-gene dosage plus genotype in plain terms), a 3-row body
(deletion / evidence / point-mutation scan), and a muted footer for the panel, its
version, and what was NOT examined. Site hits are named from the panel
(`HBA2:c.377` becomes "Hb Quong Sze") instead of shown as coordinates. Deletion
alleles are reported as a GROUP (`--SEA|--MED`) when depth cannot distinguish them.
The card reports measurements only and makes no clinical call.

**QC published and confirmed honest.** mosdepth summary, flagstat, and insert-size
now publish to `results/<sample>/qc/`, so the report's QC section is reproducible
after a work-dir cleanup (Rule 5). flagstat runs on `ch_final_bam`, which for BAM
inputs is the pre-FILTER_CHROMS input, so mapped% is the true alignment rate
measured against the reads as aligned to full hg38, not a post-filter 100%. A
comment in `subworkflows/preprocess.nf` records that this ordering is deliberate.

**Ops tooling.** `bin/nf-run-headless.sh` launches a run fully detached (`setsid`,
no controlling TTY, records the PID). `bin/nf-status.sh` reports liveness by
`kill -0` on the recorded PID and per-task state from `.exitcode` files, and refuses
a work dir that does not exist rather than reporting a false zero. Both are guarded
by `tests/test_nf_status.py`. FILTER_CHROMS was also measured properly: about 70 min
for a 79 GB BAM at 19 MB/s, awk-bound on one core, so `--max_cpus` does not help it.

Suite grew to 333 tests across 21 modules.

## 2026-07-22 (later) — the contract was wrong three ways; unfrozen and fixed

**Problem.** Three defects in `docs/contracts/alpha_globin_contract.md`, all
raised as escalations by the M8 agent and all genuine. In two of them the
implementation had already diverged from the contract *to stay honest*, which is
the right instinct but leaves spec and code disagreeing.

1. `alpha_genes_called` was declared `int 0-4`, on the unexamined assumption that
   alpha variation only ever REMOVES genes. `anti-3.7` is the reciprocal product
   of the `-a3.7` NAHR and ADDS one, so a carrier has 5 and a homozygote 6. The
   range forced `bin/alpha_globin.py` to emit `NA` for a count it had determined
   perfectly well — reporting a measurement failure that had not occurred.
2. `genotype` was illustrated as `--SEA/aQSa`, writing a site variant INTO a
   haplotype. Short reads at this locus do not establish which chromosome a site
   variant sits on, and on a deletion background that placement is the whole
   clinical question: `--SEA` in trans to Quong Sze is HbH disease, in cis it is
   a carrier who also has a variant on the intact chromosome. Same error class as
   reporting an unphased compound het as "affected".
3. Only ONE fixture existed, showing the *resolved* `--SEA/aa` form. The group
   form is what a depth-only run emits and carries the requirement that OmniGen
   render it verbatim — so the risky path was the untested one. A consumer could
   pass every test while truncating `--SEA|--MED` to its first member, which is
   exactly the bug that later surfaced in the DTC parser.

**Fix.** Range widened to `0-6`; `bin/alpha_globin.py`'s guard now catches only
genuinely impossible counts, and `NA` means only "could not be determined".
`genotype` documented as `<deletion genotype> [ +<site findings> ]`, codifying
what the implementation already did. `validation/examples/` now holds three
fixtures — resolved, group, triplication — and both tracks test against all of
them.

**Consequence.** Suite 299 -> 302. `test_group_form_is_exercised_by_a_fixture`,
`test_triplication_fixture_exceeds_the_old_cap` and
`test_genotype_never_writes_a_site_variant_into_a_haplotype` guard the three.
`test_triplication_is_reported_but_gene_count_is_NA` was renamed and inverted
rather than deleted, so the behaviour change is visible in the diff.

The consumer needed updating in the same breath: OmniGen's `PHENO` map was keyed
0-4, so a count of 5 fell to a default reading "Out-of-range gene count" — in
range, wrongly labelled — and its detail line would have read "5 of 4 alpha
genes". Widening a contract without checking the consumer would have traded one
defect for a worse one. See `OmniGen/docs/CHANGES.md`.

---

## 2026-07-22 — M8 alpha-globin: the four measurement channels

**Problem.** The α-globin locus was a total blind spot. Across HG001–HG007 the pipeline
had never produced a chr16 CNV call below 14.6 Mb and no 3.7 kb deletion in any
`sv_merged.vcf.gz` — including the Han Chinese trio HG005/6/7, where `-α3.7`/`--SEA`
carriers are statistically likely. Meanwhile OmniGen's carrier panel *already claimed*
Alpha Thalassemia as a screened condition and rendered it **negative** for THAL1, a
confirmed `--SEA` carrier, because a gene-based variant lookup structurally cannot see a
deletion. ~80–90% of α-thalassemia is deletional. That is a false negative produced by
scope, not by a crash, so no guard could catch it.

**Fix.** A new M8 subworkflow, `subworkflows/alpha_globin.nf`, with four independent
evidence channels and an integrator:

| | file | what it measures |
|---|---|---|
| ch1 | `bin/hba_depth.py` | per-segment normalized depth over the 5 diagnostic segments |
| ch2 | `bin/alpha_globin.py` | deletion-allele naming from the ch1 signature |
| ch3 | `bin/hba_junction.py` | split-read / discordant-pair breakpoints + zygosity |
| ch4 | `bin/hba_sites.py` | targeted pileup at the pinned pathogenic-site panel |

Output is the frozen contract at `results/<S>/alpha_globin/<S>.alpha_globin.tsv`
(`docs/contracts/alpha_globin_contract.md`), discovered by OmniGen via path convention.
`bin/hba_report.py` renders a thin factual card. Wired into `main.nf` /
`workflows/svcaller.nf` behind `--skip_alpha_globin` (default false); `nextflow run
main.nf --help` compiles all modules and reaches parameter validation.

**Consequence.** SVcaller now measures α-gene dosage, deletion alleles and targeted-site
genotypes. It still does **not** interpret: no HbH/Bart's/trait classification, no
couple-level risk, and `interpretation_complete` is a module constant that no code path
can set true. OmniGen owns every clinical statement, and the rewire of its Alpha
Thalassemia row to consume this contract is still outstanding — until it happens, the
live false negative described above stands.

### Four traps that shaped the implementation

**1. Intact depth is not 1.0, so the threshold goes on `score`, not the ratio.**
`score = (segment_depth / control_depth) / baseline`. Intact HBA2 sits at ratio **0.750**
and HBZ at **0.760** across GIAB HG002–HG007, so a naive `ratio < 0.8 = loss` calls a het
loss in *all six normals*. Baselines are col 5 of `assets/hba_segments.bed`; raw
calibration in `validation/giab_alpha_baseline.tsv`. `bin/hba_depth.py` consumes both and
never thresholds a raw ratio.

**2. `INTER_Z_A` is `do_not_average` and is treated as no observation at all.** It reads
0.99 ("intact") in THAL1, where a `--SEA` deletion covers roughly half of it, because
mapping inflation over chr16:155000-162000 cancels the deletion out. `observed_copies()`
in `bin/alpha_globin.py` drops it rather than letting it vote for "intact" — omission and
a vote are not the same thing, and the difference is a missed 20 kb deletion.

**3. Degenerate groups are emitted as groups.** `--SEA|--MED` and `--FIL|--THAI` have
identical depth signatures. `bin/alpha_globin.py::_collapse_group()` always returns `None`
and says why: the contract permits collapsing only on a junction read or a measured extent
that *excludes* the alternative, and `assets/hba_deletion_alleles.tsv` deliberately carries
no breakpoints (α-cluster NAHR breakpoints sit inside near-identical homology boxes and are
not single-valued per allele), so there is nothing to compare an extent against. Picking
`--SEA` because a sample looks SE Asian is a population inference dressed as a measurement.

**4. Site zygosity is copy-number dependent.** On a `--SEA/αα` background the surviving
HBA2 is hemizygous, so a real variant sits near 100% VAF, not 50%. `bin/hba_sites.py`
takes the α-gene count as an input, emits the raw VAF alongside the call on every row, and
records which rule fired in `zygosity_basis`. With no α-gene count it degrades to a
VAF-only report rather than silently assuming 2 copies. THAL1 exercises the related
`no_call` path for real: its HBA2 is single-copy, so DP at chr16:173208 (Hb Adana) is 4 —
reporting "absent" there would fabricate a negative on a severe allele.

### Defect found in a committed asset: `hba_deletion_alleles.tsv` could not express `-α3.7`

**Problem.** The generated allele table carried copy-change columns for `HBZ`, `HBA2` and
`HBA1` only. But `assets/hba_segments.bed` states plainly that **INTER_A2_A1**, not
HBA1/HBA2, is the diagnostic segment for `-α3.7` — the commonest α-thal deletion
worldwide — and `-α3.7`'s only two entries in the table were both `h` (disrupted/hybrid),
a qualitative marker no threshold can match. The table was therefore unusable by a caller
for that allele: channel 2 would have had to hardcode the `-α3.7` rule in Python, which is
exactly what defining alleles by signature exists to prevent.

**Fix.** Added a `d_INTER_A2_A1` column to `bin/make_hba_deletion_alleles.py` (segments in
genomic order) and regenerated. `-α4.2` gets `h` there rather than `0`: its deletion is
~4.2 kb but HBA2's gene body is only 835 bp, so the X2 crossover lies *inside*
INTER_A2_A1 and part of that segment goes with it — the fraction is not known here, so the
segment must not be thresholded for that allele. `INTER_Z_A` deliberately gets no column
at all, since no allele may be defined against a `do_not_average` segment.

**Consequence.** The two frozen degenerate groups are **unchanged** (`--SEA|--MED`,
`--FIL|--THAI`), so the contract is not affected. `hba_segments.bed` regenerated
byte-identical. A test now rejects an allele table without the column rather than silently
losing the allele.

### `marginal` was computed and thrown away

**Problem.** `bin/hba_depth.py` computed a `marginal` flag for every segment — whether the
score sits near a decision boundary — into its row dict, but `marginal` was missing from
`COLUMNS`, so it never reached the TSV. `bin/alpha_globin.py` reads that flag to downgrade
`alpha_genes_confidence` to `low`, so every marginal call would have been reported as
confident.

**Fix.** Added to `COLUMNS`; test asserts it reaches the file.

### The exit-126 bug recurred, and is now actually guarded

**Problem.** Four of the five new `bin/` scripts were committed mode `100644`. Every module
does `export PATH=${projectDir}/bin:$PATH` and then calls the script by bare name, so the
host checkout's `bin/` shadows the container's `/usr/local/bin` (where `Dockerfile.utils`
runs `chmod +x`). Without the exec bit the process dies with **exit 126**. This is exactly
the incident recorded at `docs/CHANGES.md` 2026-07-13 — "CNV_TRAITS scripts committed
non-executable" — which was fixed at the time but **never guarded**, so it recurred at the
first opportunity.

**Fix.** `chmod +x` on the four, plus
`tests/test_no_empty_placeholders.py::test_every_script_invoked_from_nextflow_is_executable`,
which scans every `.nf` for bare `*.py` invocations and asserts the **git index** mode
(not the working-tree mode — the index is what gets committed) is `100755`. Tool-provided
scripts that live in their own container (`configManta.py`) are excluded. Mutation-tested:
`git update-index --chmod=-x` → fails, restored → passes.

**Consequence.** `bin/cnv_traits_common.py` stays `100644` and correctly so — it is
imported, never invoked.

### Adversarial re-derivation caught a defect the implementation had shipped

Every numeric claim above was re-derived from the raw BAMs by a separate agent
whose brief was to **refute** it, with no sight of the implementers' reasoning.
Most claims survived exactly: all 28 ratios in `validation/giab_alpha_baseline.tsv`
reproduce to 3 dp, the fixture breakpoint and THAL1's real junction were both
re-derived independently to the base, and an exhaustive sweep of all 256
four-segment signatures confirmed a bare `--SEA` is unreachable. Three things did
not survive.

**1. `h` was doing two incompatible jobs — a silent 3-gene carrier was being
called a 2-gene trait.** `_delta('h')` returned `0`, which does not mean "do not
threshold this segment"; it hard-asserts *exactly 2 copies*. That is correct and
load-bearing for `-α3.7`, whose hybrid gene keeps enough of both gene bodies to
read intact. It is wrong for `-α4.2`: that deletion is ~4.2 kb while HBA2's gene
body is only 835 bp, so 2–3.4 kb of the 2969 bp INTER_A2_A1 goes with it and the
segment scores anywhere from ~0.43 to ~0.77 — straddling the 0.65 cutoff. The
note added to `-α4.2` said "do not threshold it for this allele"; the code did
the opposite. Result: `(HBZ 2, HBA2 1, INTER 1, HBA1 2)` came out as
`-a3.7/-a4.2` (2 genes) instead of `-a4.2/aa` (3 genes) — a **silent 3-gene
carrier reported as α-thalassemia trait, flipped by a 0.03 wobble in depth**.

*Fix.* Split the symbol. `h` (hybrid) → `{0}`, an expectation; new `?`
(unconstrained) → `{-1, 0}`, genuinely not compared. `+1` is excluded from `?`
because a deletion allele cannot gain copies. `expected_copies()` now returns a
*set*. Both readings of INTER_A2_A1 now give `-a4.2/aa`, 3 genes.

**2. Group members were joined alphabetically.** The code emitted `--MED|--SEA`;
the asset's own `depth_distinguishable` cells and every occurrence in the
contract write `--SEA|--MED`. OmniGen renders this string, so its spelling is
interface, not cosmetics. Group order now follows the allele table's row order.

**3. `validation/thal_truth_table.tsv` recorded `DP=23`** for THAL2's Hb Quong
Sze site. Not reproducible at any setting — a sweep over `-Q` × `-q` gives 25 /
22 / 20 / 17, never 23 — and internally inconsistent with its own VAF
(12/23 = 0.52, not the 0.55 recorded). Corrected to `DP=22 ref 10 / alt 12
VAF 0.545` with the exact flags, and THAL1's junction extent added as evidence
(method `depth+junction`). The `--SEA` vs `--MED` group is **still not
collapsed**: the breakpoint is now measured, but no GRCh38 coordinates exist for
either allele to compare 19304 bp against, so it excludes neither.

Also hardened: an unknown segment key raised a raw `KeyError`; it now raises
`AlphaGlobinInputError` naming the column and why `INTER_Z_A` has none.

**Two things the verification flagged that were NOT fixed**, and are live risks:

- **The intact/het-loss margin is thin on HBZ.** The lowest GIAB normal score is
  0.826 (HG004 HBZ) against a 0.65 cutoff — only **1.33 sd** of the observed
  intact spread. HBZ alone decides `--SEA|--MED` vs `--FIL|--THAI`, so that
  discrimination rests on the locus's noisiest segment. HBA1 is 2.11 sd; HBA2
  and INTER_A2_A1 are comfortable at 3.45 and 4.95 sd.
- **The gain threshold is not safe.** At 1.35, the parametric false-gain rate is
  ~0.63% per sample on HBZ and ~0.28% on HBA1, and an `anti-3.7` at 1.5 sits
  only ~1.1–1.3 sd above the cut. No triplication control exists. `anti-3.7`
  already yields `alpha_genes_called=NA`, so a gain is a lead, not a call —
  which is the only reason this is tolerable.

### A parsimony rule, and the blind spot it creates

**Problem.** A `--SEA` deletion on one chromosome and an `anti-3.7` triplication on the
other restore every diagnostic segment to 2 copies. Enumerating all genotypes therefore
made **every normal sample** report as ambiguous between "no deletion" and "a deletion
masked by a triplication" — technically true, operationally useless, and a guaranteed
source of false alarms.

**Fix.** `name_alleles()` keeps only the genotypes with the fewest non-wild-type
haplotypes.

**Consequence.** A compensated deletion+triplication carrier reads as normal. That is a
genuine, if remote, blind spot, recorded here and in the module docstring rather than as a
caveat printed on every sample — a caveat on every sample is a caveat nobody reads.

### What was deliberately NOT done

- **No end-to-end pipeline run.** THAL1/THAL2 were exercised only through targeted
  `samtools` region queries. A supervised full run is still required before any use.
- **The card is not wired into `report.nf`.** `REPORT` joins on the full meta map, where
  one divergent key silently drops a sample's entire report with no error; adding a 23rd
  channel without being able to run the pipeline is how that bug gets reintroduced.
- **`svcaller/utils:1.2` was not rebuilt.** The image bakes in `bin/` and `assets/`, so
  until it is rebuilt the pipeline runs the *old* code and none of the above exists at
  runtime.
- **No genome-wide variant calling.** Channel 4 is a pileup at fixed coordinates only.

## 2026-07-16 — Portability: stop hardcoding a host-specific scratch path

**Problem.** Two configs baked a host-specific path into every run, breaking the pipeline
on any other machine:

- `conf/docker.config`: `docker.runOptions = '--rm -v <fixed host scratch path>:/tmp'`
- `nextflow.config`: `env.TMPDIR = '<fixed host scratch path>'`

**Fix.** Introduced a `params.tmp_dir` that defaults to `System.getenv('TMPDIR') ?: '/tmp'`.
Both the Docker bind mount (`--rm -v ${params.tmp_dir}:/tmp`) and `env.TMPDIR = params.tmp_dir`
now derive from it. Behaviour is preserved on hosts that export `TMPDIR`; portable everywhere
else. Override with `--tmp_dir <path>` for a dedicated scratch volume. `nextflow config` and
`nextflow config -profile docker` both parse clean.

## 2026-07-15 — SvABA never actually ran: stage its classic BWA index

**Root cause.** `modules/svaba/call.nf` invoked `svaba run -G ${ref_fasta} ...` while
declaring only `path ref_fasta` + `path ref_fai` as inputs. SvABA internally calls
`bwa_idx_load_from_disk`, which loads the **classic** BWA index
(`hg38.canonical.fa.{amb,ann,bwt,pac,sa}`) from disk next to the reference. Because those
five files were never declared as a staged input, Nextflow never symlinked them into the
isolated task work dir, and SvABA died at startup with:

```
[E::bwa_idx_load_from_disk] fail to locate the index files
Unable to open index file: hg38.canonical.fa
```

**Why nobody noticed.** The `svaba run` line used to end in `2>&1 || true`, so the crash
exited 0 and fell through to the empty-VCF stub — SvABA contributed **zero** SVs on every
run in the repo's history, including the committed HG002 demo (which was really Manta+Delly
+GRIDSS+Scramble+MELT only). Today's fail-loud change (removing `|| true`) unmasked the
latent staging bug. Note the classic BWA index is a **different format** from the bwa-mem2
alignment index (`.0123` / `.bwt.2bit.64`) — the two are not interchangeable, which is why
having the alignment index present did not help SvABA.

**Fix.** `modules/svaba/call.nf` now declares `path bwa_index`, and the index is threaded
`main.nf` → `workflows/svcaller.nf` (`ch_bwa_index`) → `subworkflows/sv_calling.nf` →
`SVABA_CALL`, so Nextflow stages `.amb/.ann/.bwt/.pac/.sa` alongside `${ref_fasta}`. A new
`--bwa_index` param overrides the prefix (defaults to `ref_fasta`). Unless `--skip_svaba` is
set, `main.nf` validates the five files exist and fails loud with an actionable message
(`classic BWA index not found for <ref>; build with 'bwa index <ref>' or pass --skip_svaba`).
`--skip_svaba` (SVABA_STUB) still works and requires no classic index.

**Impact.** HG002/HG003 (and any prior result) gained no SvABA calls; they need a re-run to
pick up SvABA's contribution to the merged SV set.

## 2026-07-15 — Fail-loud guards on empty report inputs + drop dead SURVIVOR_MERGE selector

**Silent-failure class (the OmniGen "clean bill of health" bug), now closed in the reporters.**
A sibling consumer (OmniGen) had a bug where a crashed upstream caller left a 0-byte
placeholder that, being gated only on `os.path.exists()` (True for an empty file), rendered
as a clean *negative* result instead of a failure. Two SVcaller reporters had the same latent
shape:

- `bin/smn_report.py::parse_smn_tsv` opened the SMN TSV with no non-empty check. A 0-byte or
  header-only file (a crashed `SMNCopyNumberCaller`) fell through to `if not lines: return`
  defaults, and `render_html_section` then defaulted `smn1/smn2` to `2` — rendering **"Normal
  (CN=2)"** for a caller that never ran. Now `parse_smn_tsv` raises `SmnInputError` on a
  present-but-empty / header-only-no-data TSV; a `NO_*` sentinel (a legitimate skip) still
  returns the neutral/unknown result.
- `bin/html_report.py::render_report` read its required inputs (`--smn-html`, `--circos-svg`,
  `--sv-tsv`, `--cnv-bed`) directly (`Path(...).read_text()` / `parse_*`), so a 0-byte input
  rendered an empty section / "no findings". Added `UpstreamEmptyError` + `_is_absent` +
  `_require_nonempty`, called at the top of `render_report`. A present-but-empty required
  input now fails loudly.

The distinction mirrors OmniGen `prototype/upstream.py`: **absent != empty != populated.**
ABSENT (a Nextflow `NO_*` sentinel, e.g. `--cnv-bed NO_FILE` for a sample with no CNV data) is
a legitimate skip and passes silently. EMPTY (0 bytes / whitespace — a crash placeholder) is a
hard, visible failure. A **header-only** file ("we looked, found nothing") remains a valid
negative result and is deliberately allowed through — the existing empty-AnnotSV-TSV → merged-VCF
fallback still works. Healthy-path output is unchanged: verified byte-for-byte-shape against the
real `results/HG002/` inputs (report renders at 1.28 MB, matching the committed 1.29 MB), and the
guards do not fire on non-empty inputs.

New regression tests (`tests/test_smn_report.py`, `tests/test_html_report.py`) assert a 0-byte /
header-only required input raises rather than producing blank output, and that a `NO_FILE`
sentinel is still an honest skip. Suite: 67 → **77 passed**.

**Dead config selector removed.** `conf/docker.config` carried
`withName: 'SURVIVOR_MERGE' { container = '...survivor...' }`, but the active pipeline merges SVs
with `JASMINE_MERGE` (`subworkflows/sv_calling.nf`). `SURVIVOR_MERGE` is referenced only by the
separate `workflows/sv_pon_build.nf` panel-of-normals builder, and `modules/survivor/merge.nf`
already declares its container inline — so the selector was dead weight that made Nextflow emit
`WARN: There's no process matching config selector: SURVIVOR_MERGE` on every main-pipeline run.
Removed. `modules/survivor/merge.nf` is **not** orphaned (still used by `sv_pon_build.nf`) and was
left in place.

## 2026-07-13 — CNV_TRAITS scripts committed non-executable (exit 126)

`bin/gst_null.py`, `bin/amy1_cn.py`, `bin/lpa_kiv2.py` and `bin/rh_status.py` are invoked
as bare commands on `PATH` by `subworkflows/cnv_traits.nf`, but were committed with mode
`100644`. On any fresh checkout the CNV_TRAITS processes failed immediately with
`.command.sh: Permission denied` and exit status 126 (observed on a GIAB HG003 run).
Earlier HG001/HG002 runs only worked because the executable bit had been set locally in
that working tree; it was never carried in the git index.

The executable bit is now set in the git index (mode `100755`) for those four scripts.
All four already carried a `#!/usr/bin/env python3` shebang; no script logic changed.
`bin/cnv_traits_common.py` stays `100644` (imported as a Python module, never executed),
as does `bin/nf-cleanup.sh` (invoked as `bash bin/nf-cleanup.sh`).

## 2026-07-13 — Fail loudly instead of publishing empty placeholder outputs

### Why

A production incident traced back to this repo.

SVcaller processes used to end their command blocks with an `|| touch <output>`
fallback (and the JSON equivalent, `echo '{}' > <output>`). When a caller **failed**,
the process still exited 0 and **published a zero-byte output file**. Nextflow saw a
satisfied output declaration, marked the task successful, and the run continued.

Downstream, OmniGen (a consumer genomics report) gated its reads on
`os.path.exists()` — which returns `True` for a 0-byte file. An empty `smn.tsv`
crashed its scan, and the report renderer then produced a **complete-looking consumer
report** whose summary tiles read *"0 Carrier findings, 0 Medication flags, Clear"*.

**A crashed SMN caller was rendered as a clean bill of health for a human being.**

OmniGen has since been hardened to refuse empty artifacts, but SVcaller was the thing
*generating* them. Any consumer — current or future — stayed exposed, and a failed
stage still looked like a completed run on disk. This change fixes the source.

### The distinction that matters

Not every empty output is a bug, and this change is careful not to break working runs:

| Case | Example | Treatment |
|---|---|---|
| **Legitimate empty result** | A VCF with a proper header and zero variant rows — the caller ran, looked, and genuinely found nothing. | **Still allowed.** Preserved as-is. |
| **Explicitly skipped stage** | `--skip_melt`, `--skip_scramble`, `--skip_svaba`, `--skip_gridss` route to a dedicated `*_STUB` process. | **Still allowed.** The skip is a deliberate, recorded choice. |
| **Masked failure** | The caller crashed, or was never installed, and a placeholder file was written so the process could exit 0. | **Now a hard failure.** No placeholder is written. |

The rule: *a crashed caller must never be indistinguishable from a negative result.*

### Sites fixed

| File:line (original) | Original line | Change |
|---|---|---|
| `modules/smn_caller/call.nf:34` | `touch ${meta.id}.smn.tsv` | Removed. Fails with exit 1 if no result table was produced, or if the table has no sample row (header only). The originating site of the incident. |
| `modules/smn_caller/call.nf:38` | `echo '{}' > ${meta.id}.smn_detail.json` | Removed. `smn_caller.py` emits `.tsv` and `.json` together; a missing JSON means an incomplete run. |
| `modules/annotsv/annotate.nf:36` | `[ -n "$f" ] && mv "$f" . \|\| touch ${meta.id}.annotated.tsv` | Now distinguishes the two cases: if the **input VCF has zero SV records**, a header-only TSV is emitted (a legitimate empty annotation). If the input **had** SVs and AnnotSV still produced nothing, the annotation failed silently → exit 1. |
| `modules/expansionhunter/call.nf:38` | `echo '{}' > ${meta.id}.str_profile.json` | Removed → exit 1. This file is published to `results/` and read directly by OmniGen; a `{}` profile reads as "no repeat expansions detected". |
| `modules/svaba/call.nf:23` | `svaba run ... 2>&1 \|\| true` | `\|\| true` removed. A SvABA crash used to fall through to the empty-VCF stub below it and render as "0 SVs found". Also asserts SvABA actually wrote its SV VCF, and that the normalised VCF carries a `#CHROM` header. |
| `modules/scramble/call.nf:26` | `scramble.sh ... --eval-meis \|\| true` | `\|\| true` removed. A SCRAMble crash used to fall through to the header-only-VCF branch and render as "no MEIs detected". The header-only branch is **kept** for the genuine zero-MEI case. |
| `modules/melt/call.nf:25-27` | `echo "WARNING: MELT.jar not found; emitting empty VCF"` + `exit 0` | A missing MELT install is a **misconfiguration**, not a result. Now exit 1, pointing the user at `--skip_melt`. |
| `modules/melt/call.nf:40-42` | `echo "WARNING: no *_MELT.zip found ...; emitting empty VCF"` + `exit 0` | Same: missing ME references → exit 1, not a confident empty VCF. |
| `modules/melt/call.nf:61` | `java ... Single ... 2>&1 \|\| true` | Per-ME-type exit codes are now collected; if any ME type crashed, the process fails rather than reporting zero insertions. |
| `modules/melt/call.nf:71-73` | no `*.final_comp.vcf` → empty VCF + `exit 0` | If every MELT type exited 0 but nothing was written, MELT did not really run → exit 1. |

Also fixed while in `melt/call.nf`: when `--melt_refs` was supplied, `MELT_DIR` was never
set, so the `-n` gene-annotation BED path (`${MELT_DIR}/add_bed_files/...`) silently
resolved to a bogus location. `MELT_DIR` is now set on both branches.

### Deliberately left alone (legitimate, not failure-masking)

- `modules/sv_pon/annotate.nf:21` — `... | cut -f3 > pon_hit_ids.txt || true`. `grep` exits 1
  when there are **no matches**, which is the normal "this sample hits no PON sites" case.
  An empty `pon_hit_ids.txt` here is a real result.
- `modules/truvari/bench.nf:38-44` — four `|| true` on Truvari size-bin benchmarks. Truvari
  legitimately exits non-zero when a size bin contains zero calls. This is the
  **validation/benchmark** path, not the clinical report path.
- `modules/jasmine/merge.nf:93,95` — `grep -cv '^#' ... || echo 0`. Counting variants for a log
  line; `grep -c` exits 1 on zero matches.
- `modules/melt/call.nf:19` (`bowtie2 --version || true`) and `modules/scramble/call.nf:41`
  (version-string fallback) — cosmetic version capture for `versions.yml`.
- `modules/svaba/call.nf` cleanup `rm -f ... 2>/dev/null || true` — intermediate file cleanup.
- `subworkflows/report.nf:130,244` — `.ifEmpty([])` on MultiQC inputs. MultiQC is genuinely
  optional and an absent QC report does not misrepresent a clinical finding.
- The `*_STUB` processes (`melt/stub.nf`, `scramble/stub.nf`, `gridss/stub.nf`,
  `SVABA_STUB`) — these emit header-only VCFs *by design*, and are only reachable via an
  explicit `--skip_*` flag. That is the sanctioned way to express "this stage is optional
  for this sample": a recorded decision, not a fake empty file.

### Blast radius

`conf/base.config:6` is the only `errorStrategy` in the repo:

```groovy
errorStrategy = { task.exitStatus in [143,137,104,134,139] ? 'retry' : 'finish' }
```

There is no `errorStrategy 'ignore'` anywhere. The exit codes above are OOM/kill signals;
a deliberate `exit 1` is **not** in that list, so a newly-failing process will **not** be
retried in a loop. It gets `'finish'`: Nextflow stops launching new tasks, lets in-flight
tasks drain, and then reports the failing process with its work dir, exit status, and the
`.command.err` we now write. That is an actionable error, not a confusing partial run.

**Stages that now hard-fail where they previously passed silently:**

1. `SMN_CALLER` — when `smn_caller.py` produces no result table, or a header-only table.
2. `EXPANSIONHUNTER` — when EH exits 0 without writing its STR profile JSON.
3. `ANNOTSV` — when the input VCF contains SV records but AnnotSV emits no annotation.
   (Zero SVs in → header-only TSV out; still passes.)
4. `SVABA_CALL` — when SvABA itself exits non-zero, or writes no SV VCF.
5. `SCRAMBLE_CALL` — when `scramble.sh` exits non-zero.
6. `MELT_CALL` — when MELT is not installed, ME references are missing, any ME type crashes,
   or no `final_comp.vcf` is produced. Use `--skip_melt` to run without MELT on purpose.

These are all cases that previously produced a **successful-looking run with a fabricated
empty result**. If any of them starts failing on a pipeline that "used to work", that
pipeline was not working — it was silently reporting nothing found.

### Tests

`tests/test_no_empty_placeholders.py` (9 tests) statically asserts the masking patterns stay
dead: no `|| touch` / bare `touch <output>` fallback in any `.nf` file, no `echo '{}' >`
placeholder, and no `|| true` swallowing the exit code of the SvABA or SCRAMble invocations.
Each guard was mutation-tested (bug reintroduced → test fails; reverted → test passes).

Full suite: **67 passed** (58 pre-existing + 9 new).

Not verified by these tests: end-to-end Nextflow behaviour. No pipeline run was performed
for this change (code, docs and unit tests only). The shell blocks of every edited module
were syntax-checked with `bash -n` after rendering the Nextflow interpolation.

### Follow-up fix — Groovy escape bug in `smn_caller/call.nf`

The fail-loud edit introduced a comment in the `script:` block of
`modules/smn_caller/call.nf` (line 31) that contained backslash-escaped backticks. Inside a
Groovy triple-quoted GString those are invalid escapes, so **every** `nextflow run` aborted
before any process started:

```
ERROR ~ Module compilation error
- file : modules/smn_caller/call.nf
- cause: token recognition error at: '`touch\`' @ line 31, column 27.
```

The comment was reworded to use single quotes instead of backticks. No pipeline logic
changed: the `|| touch` fallback stays removed and `errorStrategy` is untouched. Verified
with `nextflow run main.nf --help`, which now compiles all modules and reaches parameter
validation (`ERROR: --input is required`) instead of failing at module compilation.
