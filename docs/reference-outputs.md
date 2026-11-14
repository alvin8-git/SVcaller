# Output file reference

Column-level detail for the per-sample outputs listed in the [README](../README.md#outputs). Moved here from the README on 2026-11-14 without change of content.

## SV caller order (`SUPP_VEC`)

**SUPP_VEC order (downstream parsers: read it, do not hard-code it).** One bit per caller that ran, left to right, in the order `Manta, Delly, [GRIDSS], [Scramble], [MELT]`. GRIDSS, Scramble and MELT bits exist only when opted in (Scramble/MELT only when they have calls). Default run: `SUPP_VEC` = `Manta,Delly`, always width 2. The old default was `Manta,Delly,GRIDSS[,Scramble][,MELT]` (width 3-5), so **old tables have GRIDSS at bit 3; new tables do not**. The actual order is in `<id>.sv_callers.json` (`callers_in_supp_vec_order`, `supp_vec_width`), the VCF header `##SVcaller_callers=Manta,Delly`, and per row in `filtered.tsv` columns `CALLERS_RAN`, `SUPP_VEC_WIDTH`, `SUPP_CALLERS`, `CALLER_CONSENSUS` (`PASS` | `SINGLE_CALLER` | `NA`; `CALLERS_RAN=legacy` on tables written before this change). Print "2 of 2 callers" from `SUPP_VEC_WIDTH`, not "2 of 4". Merge tolerances and filter thresholds are unchanged.

## `provenance.json` (producing-software versions)

At the end of every run (`workflow.onComplete`, `bin/write_provenance.py`) one
`<outdir>/<sample>/provenance.json` is written next to `<sample>.smn.tsv`; OmniGen reads it.
It is non-blocking: any failure only logs a warning and never fails the run.

```json
{"pipeline": "SVcaller <git describe>", "svcaller_commit": "<HEAD, +-dirty>",
 "generated_utc": "...", "sample": "...", "reference": {"path": "<ref_fasta>"},
 "tools": {"manta": "1.6.0", "delly": "1.2.6", "annotsv": "3.4.6", "expansionhunter": "5.0.0",
           "gatk": "...", "samtools": "...", "...": "..."},
 "images": {"manta": "quay.io/biocontainers/manta:1.6.0--py27h9948957_6", "...": "..."}}
```

Each tool version is queried from the tool itself inside its own already-local container
(`docker run --rm <image> <tool> --version`, no pulls); if the tool or image cannot answer, the
value is `"not reported"`. The repo commit is recorded only in `svcaller_commit`, never as a
tool version. Image tags come from `conf/docker.config`. Test: `tests/test_write_provenance.py`.

## STR sizing rule (`str_summary.tsv`)

`<sample>.str_profile.json` is EH's output, untouched. `bin/eh_sizing.py` (process `EH_SIZING`) reads it and
writes `<sample>.str_summary.tsv`, one row per locus, after a `#` header note. Columns:
`locus, repeat_unit, genotype, ci, read_length` (existing EH fields), then `allele1_est, allele1_ci,
allele1_sizing, allele1_display, allele2_est, allele2_ci, allele2_sizing, allele2_display`
(allele2 is `.` for haploid calls), `genotype_display, conf, spanning_reads, flanking_reads, inrepeat_reads`.

| `alleleN_sizing` | Meaning |
|---|---|
| `spanned` | at least 2 spanning reads at exactly that repeat count |
| `flanked` | flanking/in-repeat reads only; the allele (repeats x motif length) fits within the read length |
| `unresolved` | estimate longer than the read length with no spanning-read support, or CI wider than 30 repeats |

An `unresolved` allele is shown as `>=<lower CI bound>` in `alleleN_display` and `genotype_display` (for example
ATXN2 `22/>=52`), and the row gets `conf=unresolved`. Otherwise `conf` is `high` (CI width of the larger
allele <= max(6, 15% of its size)) or `low`. Read length is `ReadLength` in the json (150 for these runs);
the thresholds (2 spanning reads, CI width 30) are constants at the top of `bin/eh_sizing.py`.

Sex and chrX loci: `INFER_SEX` (`bin/infer_sex.py`, `<sample>.sex.json`) infers sex from `samtools coverage` depth
(chrY X-degenerate 6.5-10.5 Mb and chrX 40-44 Mb versus chr1/2/3 20-24 Mb; female if chrY/auto < 0.10 and chrX/auto
>= 0.75, male if chrY/auto >= 0.25 and chrX/auto < 0.75, else unknown, the same rule as OmniGen `vcf_qc.py`) and passes
`--sex male|female` to ExpansionHunter (EH defaults to female, which genotypes chrX loci as diploid in a male).
Unknown sex runs as female and sets `sex_flag` in the json and a `# FLAG` line in the summary. The inferred sex,
ratios and method are in `<sample>.sex.json` and in the `str_summary.tsv` header comment. For a male, chrX loci
(FMR1, AR, AFF2) are hemizygous: `allele2_*` is `.` and `genotype` has one number.

Limitation statement for the report: **alleles longer than the read length cannot be sized from short reads;
such alleles are reported as a lower bound.** It is a standing method limitation, not a per-sample result.
Validation against GIAB HG002 assembly truth: `docs/howto-validate-eh-vs-q100.md`.

## Depth consistency check (`filtered.tsv`)

Every DEL and DUP of at least 1 kb in `<sample>.filtered.tsv` is checked against read depth
(`bin/depth_check.py`, the `DEPTH_CHECK` process, the last step before the table is published).
Seven columns are appended after the existing ones; no row is dropped.

| Column | Meaning |
|--------|---------|
| `DEPTH_IN` | Mean depth inside the SV, excluding 500 bp at each breakpoint |
| `DEPTH_FLANK` | Mean depth of two flanks, each min(SV length, 100 kb), either side |
| `DEPTH_GLOBAL` | The sample's autosomal median depth (same on every row; see below) |
| `DEPTH_RATIO` | `DEPTH_IN / DEPTH_FLANK`, 2 decimals; `.` when the SV is < 1 kb or flank depth < 5x |
| `DEPTH_RATIO_GLOBAL` | `DEPTH_IN / DEPTH_GLOBAL`, 2 decimals |
| `DEPTH_SUPPORT` | `supported` / `refuted` / `inconclusive` / `na` (other SV types, or no ratio) |
| `DEPTH_NOTE` | `flank inflated` if `DEPTH_FLANK / DEPTH_GLOBAL` >= 1.3, `flank depressed` if <= 0.7 (the segmental-duplication signature), else `.` |

Per-ratio bands (one config block at the top of `bin/depth_check.py`): DEL `supported` at
ratio <= 0.65, `refuted` at >= 0.85; DUP `supported` at >= 1.35, `refuted` at <= 1.15; between
is inconclusive. A heterozygous deletion reads near 0.5x, a homozygous one near 0, a
heterozygous duplication near 1.5x. Two baselines are used because a segmental duplication
inflates the local flank and makes a normal-depth region look deleted: `DEPTH_SUPPORT` is
`supported` only if **both** ratios are in the supported band, `refuted` if **either** is in the
refuted band, otherwise `inconclusive` (if the flank is too shallow for a ratio, the global
ratio alone decides). **`refuted` means the read depth contradicts the call** (for example a
"deletion" that reads 1.0x of its flanks or of the genome-wide depth); it is evidence against
the call, not proof it is false, and a caller-supported call is never removed.

`DEPTH_GLOBAL` is the median of `samtools coverage` mean depth over 300 fixed 20 kb windows
spaced evenly along chr1-22 (deterministic; zero-coverage windows skipped). It is printed on
stderr and cached in a sidecar `<out.tsv>.global_depth` so re-annotation does not recompute.
Depth uses `samtools coverage` default read filters (no duplicates, secondary, QC-fail). To
annotate an existing result without a rerun:

```bash
bin/depth_check.py results/S/S.filtered.tsv S.bam results/S/S.filtered.depth.tsv
```

## Copy-number / blood-group trait contracts

Small, stable per-sample TSVs emitted by the `CNV_TRAITS` subworkflow for
downstream interpretation (e.g. OmniGen). Each is a targeted, normalized
read-depth estimate from the aligned BAM (mosdepth over `assets/cnv_trait_regions.bed`,
`--mapq 0`, normalized to copy-number-stable `CTRL_*` control regions); the CNV
consensus BED corroborates called deletions but does not drive the call.

| File | Columns |
|------|---------|
| `<S>.trait_depth.regions.bed.gz` | headerless 5-column BED: `chrom`, `start`, `end`, `region` (RHD, AMY1_CLUSTER, GSTM1, GSTT1, LPA_KIV2, CTRL_1..8), `mean_depth` (mosdepth `--mapq 0`) |
| `bloodgroup/<S>.rh_status.tsv` | `sample`, `RHD_copies`, `Rh_status` (pos/neg), `confidence` |
| `cnv_traits/<S>.amy1.tsv` | `sample`, `AMY1_copies`, `method` |
| `cnv_traits/<S>.gst_null.tsv` | `sample`, `GSTM1` (null/present), `GSTT1` (null/present) |
| `cnv_traits/<S>.lpa_kiv2.tsv` | `sample`, `KIV2_copies`, `method` |

See [docs/omnigen-additions-plan.md](omnigen-additions-plan.md) for design,
coordinates, thresholds, and calibration notes.
