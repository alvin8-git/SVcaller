# Agent brief — OmniGen thalassemia screen

**Repo:** `/data/alvin/OmniGen` · **Branch from:** `main` (currently `9be43d3`)
**Companion track:** SVcaller (see `svcaller-alpha-globin.md`). You share a frozen
interface with it. **You may not change that interface unilaterally.**

Read first: `SVcaller/docs/contracts/alpha_globin_contract.md`. Then
`prototype/entrypoint.py:16-40` — the hard boundary that shapes this whole track.

> ⚠️ **The working tree already has uncommitted changes that are not yours** —
> six modified `examples/*.html`, plus untracked `.stitch/`, `examples/CLAUDE.md`
> and HG004–007 report HTML. Stage only files you touched. Do not `git add -A`.

---

## 1. Mission

Own the *screen* and the *interpretation*. Three pieces, smallest first.

**(a) β-globin — mostly already works. Verify, don't rebuild.**
The evidence DB holds **431 HBB P/LP variants**, and 11 of the 12 curated alleles
in `SVcaller/assets/hbb_pathogenic_sites.tsv` are already there with canonical
rsIDs (HbE `rs33950507`, IVS-II-654 `rs34451549`, CD39 `rs11549407`, CD41-42
`rs36029927`, −28 `rs33931746`, …). Your job:

- confirm the panel row renders correctly for a real carrier
- add the one gap: **CD71-72** (`c.216_217insA`, chr11:5226676)
- validate against the 48 known genotypes in
  `SVcaller/validation/hbb_1000g_carriers.tsv` (23 carriers) and
  `hbb_1000g_controls.tsv` (25 confirmed non-carriers) — sensitivity *and*
  specificity, since a caller that flags everyone scores perfectly on carriers

**(b) α-globin — consume the new contract.**
```
report_alpha_globin.py            new, mirroring report_cnv_traits.py
prototype/omnigen_report.py:543   register it in the hardcoded script list
prototype/entrypoint.py           add the tier, manifest key None
```
Read `${SV}/results/<S>/alpha_globin/<S>.alpha_globin.tsv` **by path
convention**, matching `report_cnv_traits.py`. *Not* a manifest key — that
decision is settled in the contract with its reasoning.

**(c) Replace the `partial` stopgap with real results.**
`report_carrier_panel.py` currently marks Alpha Thalassemia `partial` because a
variant scan cannot see deletions (commit `9be43d3`). Once (b) lands, those rows
should report an actual result — but `partial` must remain the fallback whenever
the contract file is absent. Never let a missing file render as negative.

## 2. Non-goals

- **Do not make OmniGen read BAMs.** `entrypoint.py:16` is explicit: *"a VCF is a
  list of variants; it is not the reads."* That boundary is load-bearing — it is
  why the α false negative was diagnosable at all. Depth, copy number and
  breakpoints come from SVcaller.
- **Do not re-call variants.** README: *"grade, don't re-call."*
- **Do not rebuild the β panel.** 11/12 alleles already exist in ClinVar.
- **Do not duplicate SVcaller's measurement** — no depth logic, no CNV inference.

## 3. Traps that already cost time

**Producer and parsers must change together.** `[PANEL …]` statuses are
hardcoded in three places: the producer (`report_carrier_panel.py`) and two
regexes (`omnigen_dtc.py:249`, `omnigen_tech.py:116`). A status added to one and
not the others is **silently dropped** — the condition vanishes from the report,
which is worse than a wrong label. `test_panel_consumers_accept_all_statuses`
guards this; do not weaken it.

**Styling has semantics.** In `omnigen_dtc.py`, `red` means "has a card above" —
a positive finding. `partial` is neither positive nor a clean negative, which is
why the membership test is explicit rather than `level != "neg"`. Keep new
statuses out of the positive set unless they really are positive.

**Render degenerate groups verbatim.** The contract may deliver
`--SEA|--MED/aa`. Display the group. Silently showing the first member invents
precision the measurement does not have, and `--SEA` vs `--MED` is a population
inference, not a result.

**`required_cols` must match the contract exactly.** OmniGen's `_load()` idiom
(see the SMN tier, `omnigen_report.py:144`) fails the run on a column
mismatch. That is the desired behaviour — do not soften it to a warning.

**`not_screened` must be rendered, not just present.** The whole failure mode
this project exists to fix is a consumer gating on `os.path.exists()` and
printing "Clear". `run_tests.py:840` already proves the fail-loud pattern for
SMN — a 0-byte file makes the report exit. Give α the same.

**RUO relaxes evidence, not truthfulness.** This is an educational DTC report, so
n=1 validation graded C/D is acceptable. Reporting "negative" for something the
method cannot see is *not* — that is a false statement, and in DTC no clinician
mediates it.

## 4. Verification — what "done" means

1. `python3 prototype/run_tests.py` green. Currently **56**; you only add.
2. β validated against all 48 known 1000G genotypes, reported as sensitivity
   **and** specificity.
3. A `--SEA` carrier no longer renders as a bare negative anywhere in the report.
4. A missing / empty / malformed contract file makes the run **fail loudly** —
   never render as clear. Test all three cases.
5. **Regenerate `examples/omnigen_HG00*_bam_report.html`** — they predate commit
   `9be43d3` and still show Alpha Thalassemia as a plain negative. They are what
   a reviewer will look at.

## 5. Documentation — required, not optional

OmniGen has no changelog. **Create `docs/CHANGES.md`** following SVcaller's
convention (`SVcaller/docs/CHANGES.md`): dated `##` sections, **Problem → Fix →
Consequence**, naming files and line numbers.

Also update in the same commit as the change:
- `docs/FINDINGS.md` — the validation record. Add α and β coverage with measured
  numbers, and state plainly what is unmeasured.
- `docs/input-coverage.md` — the VCF-vs-BAM table. α-globin is a new BAM-only
  tier; the count changes.
- `docs/ConditionsCovered.md` — if the thalassemia rows change meaning.

Record dead ends too, and anything you find to be **wrong** in these documents —
correcting them in place is expected.

## 6. Subagents — use them, with these rules

This track is smaller than SVcaller's, so parallelism matters less than
independent checking.

**Fan out:** β validation (a) and α consumption (b) are independent — run
concurrently. (c) depends on (b).

**Then verify adversarially.** For the β concordance run, spawn a *separate*
agent that re-derives sensitivity and specificity from the raw TSVs without
sight of the first agent's numbers. Concordance claims are exactly the kind that
look right and are off by a denominator.

**A dedicated false-negative hunter is worth one whole agent.** Brief: *"find any
input for which this report says negative, clear, or nothing-found about a
condition it cannot actually assess."* That is the bug class that started this
work; it deserves an adversary rather than a checklist.

Rules:
- **Subagents do not commit** — you commit, after verification.
- **Subagents never `git add -A`** (see the pre-existing dirty tree above).
- **No subagent changes the contract or SVcaller's assets.** Cross-repo edits are
  an escalation.
- Give each subagent this brief plus the contract — not a summary. The traps in
  §3 are the content.

## 7. Escalate instead of deciding — stop and ask if

- the contract is insufficient for what the report needs to say
- SVcaller emits a group you want to collapse for display
- β concordance is below ~95% on the 1000G set (suspect the harness, not the
  callers — HBB is not in a segdup and these calls should be near-perfect)
- you are about to write anything a reader could take as clinical advice

**Known-open, do not treat as settled:** THAL1/THAL2 supplier labels are
reversed and unconfirmed; THAL1 is a double heterozygote (α `--SEA` + β
IVS-II-654) resting on one VCF record at depth 19; HG02379's compound-het phasing
rests on 4 fragments and its **phenotype is unknown** — never render it as
affected.
