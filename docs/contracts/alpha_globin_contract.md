# Contract: `alpha_globin.tsv` — SVcaller → OmniGen

**Status:** frozen 2026-07-22 (phase 0). Change only by editing this file *and*
the example, together, in one commit.

This exists because the α-globin work is being built by two independent tracks —
measurement in SVcaller, screen ownership in OmniGen. They must agree on the
interface *before* either starts, or integration fails at the end when it is
most expensive. Both sides test against `validation/examples/SAMPLE.alpha_globin.tsv`.

---

## Decision 1 — discovery is by path convention, not manifest key

OmniGen ingests SVcaller output two different ways today, so "follow the
existing pattern" is ambiguous. It is resolved here:

| Existing tier | Mechanism |
|---|---|
| SMN | manifest key — `config.resolve(SAMPLE, "smn")`, read from `samples/<S>.json` |
| CNV traits | path convention — `${SV}/results/<S>/cnv_traits/<S>.gst_null.tsv` |

**α-globin uses the path convention.** `config.py` states that manifest keys
exist for *"per-sample, irregularly-named paths"*. SVcaller's own output paths
are regular and predictable, so a manifest key would add a per-sample JSON edit
that buys nothing. The tier's entry in `entrypoint.py:READ_DEPENDENT` therefore
carries `None` as its manifest key, exactly like CNV traits and Rh(D).

```
${SV}/results/<SAMPLE>/alpha_globin/<SAMPLE>.alpha_globin.tsv
```

## Decision 2 — module shape

An earlier note called α-globin "a fifth CNV trait caller". That was only half
right and the two tracks would have read it differently, so:

- **Its own subworkflow**, `subworkflows/alpha_globin.nf` — the four channels
  include junction detection and a targeted pileup, which do not fit
  `TRAIT_DEPTH`'s one-depth-pass-then-interpret shape.
- **But it reuses `bin/cnv_traits_common.py`** for normalized-depth helpers
  rather than growing a second copy of that logic.
- It consumes the same `CTRL_*` windows from `assets/cnv_trait_regions.bed`.

---

## Format

A **single-row TSV** with a header — matching how OmniGen reads `smn.tsv`
(`_load(..., required_cols=[...], expect_rows=True)` then takes row `[0]`).
Tab-separated, no quoting, `NA` for genuinely unknown.

### Required columns

OmniGen declares these in `required_cols`; a rename on either side fails the
run rather than silently mis-rendering.

| Column | Type | Meaning |
|---|---|---|
| `sample` | string | sample id, matches `meta.id` |
| `alpha_genes_called` | int 0–6, or `NA` | functional α-gene count. **Not capped at 4** — see below |
| `alpha_genes_confidence` | `high`\|`medium`\|`low` | `high` = depth and a supporting junction agree · `medium` = depth only (no junction file, or none supporting), a marginal segment changing only the genotype/gene count, or the -a4.2 vs -a3.7/-a4.2 ambiguity · `low` = depth missing or any `NA` genotype. A junction never lowers it |
| `deletion_alleles` | string | e.g. `--SEA\|--MED/aa`, `none`, `NA`. `/`-separated haplotypes, `a` = intact α. **`\|` inside a haplotype means an unresolved degenerate group — see below** |
| `deletion_evidence` | `depth`\|`junction`\|`both`\|`none` | `both` only when a `confidence == high` junction fully contains every segment the named genotype removes (no tolerance) and overlaps none it leaves untouched; `depth` = depth-named deletion alone; `junction` = a high junction overlapping a scorable segment exists but depth named no deletion (normal or `NA`); `none` |
| `site_genotypes` | string | `;`-joined `GENE:HGVS:zygosity`, or `none` |
| `site_panel_version` | string | `hba_pathogenic_sites.tsv@<sha1>` — the panel actually used |
| `genotype` | string | `<deletion genotype>` optionally ` +<site findings>`, e.g. `--SEA\|--MED/aa +HBA2:c.377T>C:het`. **Never writes a site variant into a haplotype** — see below |
| `screened` | string | `,`-joined tier ids actually run |
| `not_screened` | string | `,`-joined tier ids explicitly NOT run |
| `interpretation_complete` | `false` | always false from SVcaller; it measures, it does not interpret |

### Depth cannot name every allele — report the group, never pick one

`assets/hba_deletion_alleles.tsv` carries a `depth_distinguishable` column, and
two groups are flagged `no:`

```
--SEA | --MED     both remove HBA1+HBA2 and spare HBZ
--FIL | --THAI    both remove HBA1+HBA2 and HBZ
```

Within a group the copy-number signature is **identical**, so depth alone cannot
choose between them. A caller that emits `--SEA` from depth is inventing
precision it does not have — and `--SEA` vs `--MED` is a population inference,
not a measurement, so guessing it will look right in SE Asian samples and be
wrong elsewhere.

Emit the group: `--SEA|--MED/aa`. Collapse to a single allele **only** when
channel 3 supplies a junction read or the measured extent excludes the
alternative, and record which in `deletion_evidence`.

OmniGen must render the group as-is. Do not silently display the first member.

### When the genotype is `NA` (added 2026-11-13)

`deletion_alleles` and `alpha_genes_called` are both `NA` (never `none`/4) when: (1) a non-wild-type call rests on
marginal segments (read toward 2 copies it is wild type or only a deletion + anti-3.7 compound); (2) the only
fit is a deletion allele plus the anti-3.7 triplication (known blind spot); (3) a called-normal sample has a valid
deletion reading on a marginal segment; (4) **-a4.2 guard**, empirical: HBA2 loss with INTER_A2_A1 not called a
loss and not marginal toward loss (HBA2 alone dips on some libraries without a deletion; every confirmed -a4.2
carrier also showed INTER_A2_A1 reduced, ratio 0.47-0.60 het / 0.19-0.21 hom against 0.78-1.29 in normals);
(5) HBA2, INTER_A2_A1 or HBA1 has no usable call (HBZ unmeasured instead widens the degenerate group). The run
note names the segment(s). The junction channel only corroborates and cannot change any of this.

The guard fires only when HBA2 is itself called a loss. It is deliberately NOT extended to a marginal het_loss INTER_A2_A1 (score ~0.57-0.65): `-a4.2` is then still named, at `medium` at most (even with a supporting junction) and with a note that INTER_A2_A1 was borderline, because confirmed het carriers read INTER 0.47-0.60 (two platforms) against 0.78-1.29 in normals, so a marginal loss there is carrier-like, not a coin flip.

**Known limit:** homozygous anti-3.7 (6 genes) is not resolved: `gain` maps to 3 copies only, so it is emitted as `anti-3.7/aa` (5 genes), and any `gain` call adds the note "gain copy number above 3 is not resolved; anti-3.7 on both chromosomes (6 genes) reads the same".

### The α-gene count is not capped at 4 (widened 2026-07-22)

`alpha_genes_called` was originally declared `int 0–4`, on the unexamined
assumption that α-globin variation only ever removes genes. It does not:
`anti-3.7` is the **reciprocal product** of the same NAHR that creates `-α3.7`,
and it *adds* a gene. An `anti-3.7` carrier genuinely has **5**; homozygous
`anti-3.7` would have 6 but is emitted as 5 today (see the known limit below). The range is therefore **0–6**.

The old range forced the implementation to emit `NA` for a real, countable
result — the allele survived in `deletion_alleles` so nothing was lost, but the
count field lied by omission. A consumer that renders "α genes: NA" for a
sample whose gene count is perfectly well determined is reporting a measurement
failure that did not occur.

`NA` remains valid, and now means what it should: the count could not be
determined. It must not be used for counts the field merely could not hold.

### `genotype` does not assert phase (corrected 2026-07-22)

This field previously illustrated `--SEA/aQSa` — a site variant written **into**
a haplotype. That was wrong, and the implementation was right to diverge from it
rather than obey it.

Short reads at this locus do not establish which chromosome a site variant sits
on. Writing `aQSa` into a haplotype claims exactly that, and on a deletion
background the placement *is* the clinical question: `--SEA` in trans to a
Quong Sze allele is HbH disease, whereas in cis the person is a `--SEA` carrier
who also happens to carry a site variant on the intact chromosome. Asserting the
first from data that cannot distinguish them is the same class of error as
reporting an unphased compound het as "affected".

The format is therefore two separable measurements:

```
<deletion genotype> [ +<site findings> ]

--SEA|--MED/aa                        deletion only
--SEA|--MED/aa +HBA2:c.377T>C:het     deletion AND a site finding, phase unstated
aa/aa +HBA2:c.427T>C:het              no deletion, site finding only
```

Neither half implies the other. A consumer must not render the ` +` form as a
compound genotype, and must not reorder it into a haplotype string.

### Zygosity is copy-number dependent — do not drop the raw VAF

`site_genotypes` zygosity **must** be computed against `alpha_genes_called`. On a
`--SEA/αα` background the surviving HBA2 is hemizygous, so a real variant sits
near 100% VAF, not 50%. A naive `0/1` call at that site is wrong, and it is
wrong precisely at the compound heterozygotes that matter clinically.

Emit the raw VAF alongside the zygosity call in the per-site detail file
(`<SAMPLE>.alpha_sites.tsv`, free-form); never let the contract's zygosity
string be the only surviving evidence.

### `not_screened` is load-bearing

A consumer must never be able to gate on `os.path.exists()` and render "Clear".
At minimum `not_screened` carries `beta_globin` and
`alpha_nondeletional_outside_panel`. OmniGen must render its contents, not just
check the file is present.

`run_tests.py:794` already proves the fail-loud path for SMN — a 0-byte file
makes `omnigen_report` exit rather than render. α-globin gets the same test.

---

## Examples — three fixtures, not one

`validation/examples/` holds **three** fixtures. Both tracks must test against
all of them.

One fixture was a mistake. The original showed only the *resolved* form
(`--SEA/aa`, collapsed on junction evidence), which is the rarer and easier case.
The **group** form is what a depth-only run actually emits, and it carries the
subtle requirement that OmniGen render it verbatim rather than showing the first
member — so the risky path was the one nothing tested. A consumer could have
passed every test while silently truncating `--SEA|--MED` to `--SEA`, which is
exactly the bug that later turned up in the DTC parser.

| Fixture | Case | Exercises |
|---|---|---|
| `SAMPLE.alpha_globin.tsv` | `--SEA/aa`, evidence `both` | group collapsed by a junction read |
| `SAMPLE_group.alpha_globin.tsv` | `--SEA\|--MED/aa`, evidence `depth` | **unresolved group — must render verbatim** |
| `SAMPLE_triplication.alpha_globin.tsv` | `anti-3.7/aa`, 5 α genes | count above 4; the widened range |

```
sample	alpha_genes_called	alpha_genes_confidence	deletion_alleles	deletion_evidence	site_genotypes	site_panel_version	genotype	screened	not_screened	interpretation_complete
SAMPLE	2	high	--SEA/aa	both	none	hba_pathogenic_sites.tsv@<sha>	--SEA/aa	alpha_deletional,alpha_targeted_sites	beta_globin,alpha_nondeletional_outside_panel	false
SAMPLE_GROUP	2	medium	--SEA|--MED/aa	depth	none	hba_pathogenic_sites.tsv@<sha>	--SEA|--MED/aa	alpha_deletional,alpha_targeted_sites	beta_globin,alpha_nondeletional_outside_panel	false
SAMPLE_TRIP	5	high	anti-3.7/aa	both	HBA2:c.377T>C:het	hba_pathogenic_sites.tsv@<sha>	anti-3.7/aa +HBA2:c.377T>C:het	alpha_deletional,alpha_targeted_sites	beta_globin,alpha_nondeletional_outside_panel	false
```

Note the group fixture's `alpha_genes_confidence` is `medium`, not `high`:
depth alone resolved the gene *count* but not the allele *identity*. And the
triplication fixture shows the ` +` site form, so the phase-neutral `genotype`
format is covered too.

---

## Open items this contract does NOT settle

- **Research-use vs clinical** is undecided, and it changes report wording and
  whether the module may gate reproductive decisions. Both tracks should build
  the fields; the wording decision can land later.
- **Site panel is not frozen for content.** Coordinates are derived and
  cross-checked (`bin/make_globin_panels.py`), but the *allele list* is curated
  and incomplete — Hb Constant Spring has no validation sample, and α coverage
  beyond `--SEA` is unmeasured.
- **α validation is n=1** (THAL1) and there is no compound heterozygote.
