# Changes

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

## 2026-07-16 — Portability: stop hardcoding `/data/alvin/tmp`

**Problem.** Two configs baked a host-specific path into every run, breaking the pipeline
on any other machine:

- `conf/docker.config`: `docker.runOptions = '--rm -v /data/alvin/tmp:/tmp'`
- `nextflow.config`: `env.TMPDIR = '/data/alvin/tmp'`

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
