# Agent brief — SVcaller α-globin module

**Repo:** `/data/alvin/SVcaller` · **Branch from:** `main` (currently `586f3a9`)
**Companion track:** OmniGen (see `omnigen-thalassemia.md`). You share a frozen
interface with it. **You may not change that interface unilaterally.**

Read first, in order: `docs/contracts/alpha_globin_contract.md`, then
`docs/superpowers/plans/2026-07-21-alpha-globin-thalassemia.md`.

---

## 1. Mission

Build the α-globin measurement module: four channels emitting one contract file
per sample.

```
subworkflows/alpha_globin.nf     new — its own subworkflow, NOT a 5th cnv_trait
bin/alpha_globin.py              integrates channels -> the contract TSV
bin/hba_report.py                thin HTML card (mirror bin/smn_report.py)
```

| Ch | What | Inputs already committed |
|---|---|---|
| 1 | α gene count from depth | `assets/hba_segments.bed`, CTRL windows in `assets/cnv_trait_regions.bed` |
| 2 | deletion allele naming | `assets/hba_deletion_alleles.tsv` |
| 3 | junction / breakpoint | `tests/fixtures/alpha_junction.bam` |
| 4 | targeted pathogenic sites | `assets/hba_pathogenic_sites.tsv` |

Output: `results/<S>/alpha_globin/<S>.alpha_globin.tsv`, exactly the 11 columns
in the contract, validated against `validation/examples/SAMPLE.alpha_globin.tsv`.

## 2. Non-goals — do not do these

- **Do not touch β-globin.** It needs a VCF, not reads, and OmniGen already
  handles it. Adding an HBB channel here duplicates working code.
- **Do not run a genome-wide variant caller.** Channel 4 is a pileup at fixed
  coordinates. Genome-wide SNV calling is an explicit spec non-goal, and the
  pileup is what keeps channel 4 inside it.
- **Do not interpret.** No HbH / Bart's / trait classification, no couple risk,
  no "carrier"/"clear" language. You measure; OmniGen interprets.
  `interpretation_complete` is always `false`.
- **Do not edit generated assets by hand.** `assets/hba_*.tsv|bed` are emitted by
  `bin/make_globin_panels.py` and `bin/make_hba_deletion_alleles.py`. Change the
  curation in the generator and re-run.
- **Do not widen the contract** to make your implementation easier. Escalate (§7).

## 3. Traps that already cost time — each was a real bug

**Threshold on `score`, never on the raw ratio.**
```
ratio = segment_depth / control_depth
score = ratio / baseline          <-- threshold on THIS
```
Intact depth is **not 1.0** across most of this locus. HBA2's intact baseline is
0.750 — a naive `ratio < 0.8 = loss` calls a het loss in **all six GIAB
normals**. Baselines live in col 5 of `hba_segments.bed`.

**`INTER_Z_A` is `do_not_average`.** It reads 0.99 ("intact") in a sample where a
`--SEA` deletion covers half of it, because mapping inflation over
chr16:155000-162000 cancels the deletion out. Use mappable sub-windows or skip
it. Do not compute its mean and believe the answer.

**Never name one allele from a degenerate group.** `--SEA|--MED` and
`--FIL|--THAI` have *identical* depth signatures. Emit the group. Collapse only
on a channel-3 junction or an extent that excludes the alternative, and say
which in `deletion_evidence`. Picking `--SEA` because the sample looks SE Asian
is a population inference dressed as a measurement.

**Zygosity in channel 4 depends on channel 1.** On a `--SEA` background the
surviving HBA2 is hemizygous, so a real variant sits near 100% VAF, not 50%.
Emit the raw VAF *and* the copy-number-aware call, never the latter alone.

**Rebuild the container after editing `bin/`.** `svcaller/utils:1.2` bakes in
`bin/` and `assets/`. Edit a script, skip the rebuild, and the pipeline silently
runs the old code. Also: `html_report.py` resolves its template via
`Path(__file__).parent.parent / "assets"` — never change `bin/`'s depth relative
to `assets/`.

**Meta-map joins drop samples silently.** `REPORT` joins on the *full* meta map;
a channel whose meta differs by one key produces no report and no error. See
CLAUDE.md "Meta-map consistency". Tag every sibling channel or normalise to
`meta.id`.

**Shared reference files must be `Channel.value()`.** A queue channel is consumed
by the first subworkflow and every later one receives nothing.

**Never `docker rmi` / `prune -a`.** `svcaller/{melt:2.2.2,utils,smncopynum}` are
built locally and not re-pullable.

## 4. Verification — what "done" means

1. `pytest tests/` green. Currently **129**; you only add.
2. Tests stay pure and fast (~2 s, no containers, no network). Follow
   `tests/test_cnv_traits.py` for depth channels and
   `tests/test_junction_fixture.py` for channel 3.
3. Channel 3 recovers the fixture breakpoint **to the base**: `165000 -> 185001`,
   and calls it **heterozygous** (the fixture deliberately contains
   reference-spanning reads; a hom-del call is a failure).
4. Re-scored against the committed baselines, THAL1 and THAL2 must come out as
   `validation/thal_truth_table.tsv` says. That file is machine-read by
   `tests/test_thal_truth_table.py` — keep it that way.
5. Your emitted TSV validates against the committed example fixture.
6. Run end-to-end on THAL1 and THAL2 before declaring done. Path variables and
   the run command are in CLAUDE.md. Use `NXF_ANSI_LOG=false` for background
   runs and a per-sample `-work-dir`.

## 5. Documentation — required, not optional

Append to **`docs/CHANGES.md`** following its existing shape: a dated `##`
section, **Problem → Fix → Consequence**, naming files and line numbers. It is a
narrative record of *why*, not a commit list.

Also update, in the same commit as the change they describe:
- `docs/superpowers/plans/2026-07-21-alpha-globin-thalassemia.md` — tick off
  what you completed and **correct anything you find to be wrong**. Several
  claims in it were already wrong and were fixed in place; that is the expected
  behaviour, not an exception.
- `CLAUDE.md` — if you add params, modules, or `bin/` scripts. Its tables are
  routinely stale; keep yours accurate.
- `validation/thal_truth_table.tsv` — if you establish a genotype. Never set
  `orthogonal=yes` without naming the assay.

**Write down what you tried that did not work.** The measurement dead-ends in
this project were as valuable as the successes.

## 6. Subagents — use them, with these rules

The channels are genuinely independent, so parallelism is real. Robustness comes
from *separating implementation from verification*, not from more implementers.

**Fan out (parallel, independent):**
- one agent per channel: 1 (depth), 3 (junction), 4 (sites)
- channel 2 depends on 1 and 3 — start it after both land

**Then verify adversarially — this is the part that matters.** For anything
numeric — a threshold, a coordinate, a baseline, a breakpoint — spawn a
*separate* agent whose brief is to **refute** the result, with no sight of the
implementer's reasoning. Give it the raw data and ask it to derive the number
independently. Accept only on agreement.

Rules that keep it robust:
- **Subagents do not commit.** You commit, after verification.
- **Subagents do not edit the contract or the generated assets.** If one wants
  to, that is an escalation, not a change.
- **Every numeric claim carries its derivation** in the agent's report, so you
  can check it rather than trust it.
- **One channel per worktree** if agents run concurrently, so parallel edits
  cannot collide.
- Give each subagent this brief plus its channel's section — not a summary.
  Summaries drop the traps in §3, and the traps are the point.

## 7. Escalate instead of deciding — stop and ask if

- the contract seems wrong or insufficient (it is frozen; the other track depends on it)
- a baseline looks wrong for a sample you did not calibrate against
- you want to collapse a degenerate group on anything other than a junction or extent
- a real sample disagrees with `validation/thal_truth_table.tsv`
- you are about to write anything a reader could take as clinical advice

**Known-open, do not treat as settled:** THAL1/THAL2 supplier labels are still
reversed and unconfirmed; the GIAB baseline assumes six untyped samples are
α-normal; HG001 is an unresolved outlier (probably technical); `--SEA` vs
`--MED` is unresolved by design until channel 3 lands.
