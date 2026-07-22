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
| `alpha_genes_called` | int 0–4, or `NA` | functional α-gene count |
| `alpha_genes_confidence` | `high`\|`medium`\|`low` | depth+junction agree / depth only / marginal |
| `deletion_alleles` | string | e.g. `--SEA\|--MED/aa`, `none`, `NA`. `/`-separated haplotypes, `a` = intact α. **`\|` inside a haplotype means an unresolved degenerate group — see below** |
| `deletion_evidence` | `depth`\|`junction`\|`both`\|`none` | which channels supported the call |
| `site_genotypes` | string | `;`-joined `GENE:HGVS:zygosity`, or `none` |
| `site_panel_version` | string | `hba_pathogenic_sites.tsv@<sha1>` — the panel actually used |
| `genotype` | string | integrated, e.g. `--SEA/aQSa` |
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

## Example

`validation/examples/SAMPLE.alpha_globin.tsv` is the canonical fixture. It
encodes THAL1's expected output (`--SEA` heterozygote, no site hits):

```
sample	alpha_genes_called	alpha_genes_confidence	deletion_alleles	deletion_evidence	site_genotypes	site_panel_version	genotype	screened	not_screened	interpretation_complete
SAMPLE	2	high	--SEA/aa	both	none	hba_pathogenic_sites.tsv@0000000	--SEA/aa	alpha_deletional,alpha_targeted_sites	beta_globin,alpha_nondeletional_outside_panel	false
```

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
