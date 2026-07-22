# Plan — α-globin copy-number genotyping for thalassemia carrier screening

**Status:** proposal, not started. Written 2026-07-21.
**Scope (approved):** α-deletional **plus** a targeted pathogenic-site check at a
fixed list of HBA1/HBA2 positions. β-thalassemia measurement remains out of scope
for SVcaller.
**Intended use:** couples carrier screening.
**Validation:** THAL1 (`--SEA`) and THAL2 (Hb Quong Sze, HBA2 c.377T>C). Two
samples only; both purchased and sequenced as standard samples.
**⚠️ Sample identity corrected 2026-07-22** — the supplied labels were the
reverse of what the reads show. See § Sample identity below before using either
sample as truth.
**DRAGEN:** unavailable — the caller must be written.

---

## P0 — CONFIRMED 2026-07-22. This is a live false negative.

**Ran it. Alpha Thalassemia reports negative for a confirmed `--SEA` carrier.**

`report_carrier_panel.py:111-115` renders a row negative whenever the
genome-wide ClinVar P/LP **variant** scan returns nothing for that gene:

```python
recs = by_gene.get(gene)
if not recs:
    lines.append(f"  [PANEL neg] {cond} | {gene} |  | {d}")
```

THAL1 carries a heterozygous ~20 kb `--SEA` deletion of HBA1+HBA2 (verified from
its reads, § Sample identity). Its VCF holds 6 variants across HBA1/HBA2, all
common polymorphisms; the deletion appears nowhere, because deleting a gene
*removes* variants rather than adding one. The panel therefore prints:

```
[PANEL neg] Alpha Thalassemia | HBA1 |  | 138
[PANEL neg] Alpha Thalassemia | HBA2 |  | 167
```

**The trailing counts make it worse.** 138 and 167 are the P/LP variants
screened per gene — the report tells the reader 305 α-globin variants were
examined and none was found. That reads as thoroughness. What it cannot say is
that ~80–90% of α-thalassemia is deletional and none of those 305 records
describes a deletion. The generic footnote (*"negative = no known P/LP variant
found … not a guarantee"*) is true for the other 163 panel rows, where SNV
scanning is the right method; for α-globin the dominant disease mechanism is
structurally invisible, and the footnote does not distinguish those cases.

### The same sample proves the beta path already works

THAL1 also carries IVS-II-654 (§ THAL1 also carries a beta-thalassemia allele).
That variant **is** in the evidence catalog — `chr11:5225923 G>A, rs34451549,
classification P` — so the existing scan would flag
`Beta Chain-Related Hemoglobinopathy` as a carrier hit.

One sample, one report, both outcomes:

| | Mechanism | Catalog | Result |
|---|---|---|---|
| β — IVS-II-654 | SNV | in ClinVar as P | **caught** |
| α — `--SEA` | 20 kb deletion | not representable | **missed, rendered negative** |

That is the architecture split of this plan, demonstrated rather than argued:
β is a VCF lookup OmniGen already does correctly; α needs reads and needs the
module below.

(Caveat on the β half: THAL1 is not wired into OmniGen's sample config, so the
scan was not executed end-to-end. The catalog entry and THAL1's genotype match
on chromosome, position and both alleles, so a position-matching scan hits — but
that is inference, not a run.)

### Mitigate now, independently of building anything

The module is weeks away; the false negative is live today. Either suppress the
two Alpha Thalassemia rows or mark them not-covered, stating that deletional
α-thal is not assessed from a variant scan. That is a small change to
`report_carrier_panel.py` and does not wait on any of the work below.

---

## Original P0 note — check this before any build work

`OmniGen/prototype/carrier_panel.tsv:13-14` already lists **Alpha Thalassemia →
HBA1, HBA2** as a screened condition, alongside Beta Chain-Related
Hemoglobinopathy → HBB (`:31`). The carrier section renders *every* screened
condition, negatives included, with a running "N negative" count
(`omnigen_dtc.py:698,706-723`), fed by `[PANEL neg|carrier|atrisk|affected]`
lines parsed at `:247`.

Whatever produces those PANEL lines is doing **gene-based variant lookup**, which
structurally cannot see a deletion. So OmniGen may today be listing *Alpha
Thalassemia — negative* for a genuine `--SEA` carrier: a false negative on a
consumer carrier screen, from a claim the panel already makes.

**This is testable in about an hour and you already have the sample.** Run THAL1
— the `--SEA` deletion carrier — through OmniGen and read the Alpha Thalassemia
row. If it reports negative, that
is a live shipping defect, and it can be mitigated immediately — mark α-thal
not-covered until the module below exists — entirely independently of building
anything. Do not let this wait on the project.

I could not confirm it from `prototype/samples/*.json`; those files carry no
PANEL lines, so the producer is elsewhere and must be traced.

---

## Context

SVcaller today has no hemoglobin content of any kind — a repo-wide grep for
`HBA`, `HBA1`, `HBA2`, `HBB`, `thalass`, `globin`, `hemoglobin` returns nothing
in code, assets, or docs. The question is whether thalassemia detection and
classification belong here.

Thalassemia splits into two diseases with different variant classes, and only
one of them is this pipeline's kind of problem:

| | Locus | Variant class | In scope for an SV/CNV pipeline? |
|---|---|---|---|
| **α-thalassemia** | HBA1/HBA2, chr16p13.3 | ~80–90% **deletions** (`-α3.7`, `-α4.2`, `--SEA`, `--MED`, `--FIL`, `--THAI`) | **Yes** |
| **β-thalassemia** | HBB, chr11p15.4 | ~230 point mutations vs ~18 large deletions of ~280 known alleles | **No** |

The design spec (`docs/superpowers/specs/2026-05-09-svcnv-caller-design.md` §3)
lists "SNP/indel calling (handled by separate pipeline)" as an explicit
non-goal, which β-thalassemia squarely is.

---

## Recommendation

**Add α-globin copy-number genotyping. Do not have SVcaller classify
"thalassemia carrier status."**

Two separate answers to the two halves of the question:

- **Detection — yes**, as a dedicated paralog-aware caller (new M8 subworkflow),
  *not* through the existing generic SV/CNV path and *not* through the
  `cnv_traits` depth framework. Both are demonstrably incapable of it; evidence
  below.
- **Classification — no.** SVcaller should emit α-gene dosage, deletion-allele
  identity and targeted-site genotypes as a *measurement*. Turning that into a
  couples risk statement still requires β-globin input this pipeline will not
  have, plus non-deletional α alleles outside the targeted panel.

### Why classification must live elsewhere

For couples carrier screening specifically, α-deletional status alone is
necessary but not sufficient, and reporting it as a thalassemia result is
actively unsafe:

- **β-thal is invisible here** (~82% of known β alleles are point mutations).
- **Non-deletional α is invisible here.** Hb Constant Spring (HBA2 c.427T>C) is
  common in SE Asia and, combined with `--SEA`, produces HbH-CS — *more* severe
  than deletional HbH.
- Concrete failure: partner A `--SEA/αα` (detected); partner B `αα/αα` by
  deletion analysis but a silent HbCS or β-thal carrier (not detected). Reported
  as low risk; true risk is 25% HbH-CS, or 25% Bart's hydrops if B also carries
  `--SEA`. Bart's is lethal in utero and carries maternal morbidity.

This is the same failure class as the documented production incident in
`TODO.md`: an empty `smn.tsv` rendered by OmniGen as "0 Carrier findings …
Clear" — a crashed caller shown as a clean bill of health. Here it would arise
from **scope** rather than a crash, which is worse, because nothing is broken
for a guard to catch.

### Strategic caveat worth deciding before building

Thalassemia carrier screening already has a mature, cheap, accurate pathway:
CBC/MCV → HPLC → gap-PCR/MLPA. That pathway is the clinical reference standard
and detects both α and β. A WGS-derived α-CN call is only worth building if WGS
is *already* being run for other reasons and you want incremental value from it
— it should not displace first-line screening, and it would need to demonstrate
non-inferiority before informing reproductive decisions.

---

## Why the existing pipeline cannot do this today

Not a matter of tuning. Four independent mechanisms each remove α-thal variants,
and empirically the pipeline finds nothing at the locus.

**Empirical baseline.** Across HG001–HG007 there is **no chr16 CNV call below
14.6 Mb in any sample**, and no 3.7 kb deletion in any `sv_merged.vcf.gz` —
including the Han Chinese trio HG005/6/7, where `-α3.7`/`--SEA` carriers are
statistically likely. The locus is currently a blind spot.

| # | Mechanism | Effect on α-thal |
|---|---|---|
| 1 | `bin/cnv_consensus.py:92` `cnvpytor_only_min_bp = 1_000_000` | CNVpytor-only calls <1 Mb dropped. A 3.7 kb DEL survives only if GATK gCNV independently calls it. |
| 2 | `modules/jasmine/merge.nf:116-141` — post-merge awk drops single-caller non-INS SVs with `\|SVLEN\| > 10000` | **Deletes `--SEA` (20 kb), `--MED` (17 kb), `--FIL` (31 kb), `--THAI` (34 kb)** when only one caller finds them. These are the most clinically critical SE Asian alleles. |
| 3 | `modules/annotsv/annotate.nf:70-103` `GNOMAD_SV_FILTER`, AF 0.01, **hard-drops rows** | `-α3.7` exceeds 1% AF globally. Currently spared **only by a bug**: `row.get('B_gain_AFmax', row.get('B_loss_AFmax','0'))` — `B_gain_AFmax` is always present in the header, so `B_loss_AFmax` is never read. Fix that `.get()` and every common α-thal deletion is silently deleted. **This is a landmine.** |
| 4 | `bin/html_report.py:516,544,612-641` — CNV reaches Clinical Findings only on ≥50% overlap with a region in `assets/cnv_syndromes.tsv`, whose chr16 entries are 16p13.11/16p12.1/16p11.2 only | A correctly called `-α3.7` **would never be displayed.** |

Also: `pon/sv_pon/giab_sv_pon.bed` contains `chr16:182288-183574`, inside the
α-globin cluster ~5 kb 3′ of HBA1/HBQ1 — the locus is already partially
blacklisted (flagged, not removed; `modules/sv_pon/annotate.nf:20-35`).

**Resolution is marginal even before filtering.** GATK gCNV uses
`--bin-length 1000` (`modules/gatk/gcnv_pon.nf:23`) and CNVpytor bins at 1 kb —
3.7 kb is 3–4 bins, in a segmental duplication, with a PON built without GC or
mappability correction.

**Depth alone cannot separate HBA1 from HBA2.** This is the same problem SMN
has, and the reason `SMNCopyNumberCaller` emits a 16-value per-PSV vector
(`SMN1_CN_raw`) and an explicit paralogous-sequence-variant column
(`g.27134T>G_CN`). The `cnv_traits` framework has no equivalent machinery: it
runs `mosdepth --mapq 0` and normalizes against 8 control windows
(`bin/cnv_traits_common.py:64-74`), handling paralogy only by threshold fudging
(`rh_status.py:12` allows RHCE cross-mapping to lift a true 0-copy to ~0.3).
That is adequate for present/absent (GSTM1 null) but not for resolving 4 vs 3
α genes — a ratio of 1.0 vs 0.75.

---

## Tooling

No open-source short-read α-globin caller exists. Illumina's **DRAGEN HBA
caller** solves this and validates well against long reads, but is proprietary
and requires the DRAGEN platform. Long-read (PacBio HiFi / SMRT) is the accurate
route for rare and complex alleles. Clinical reference remains gap-PCR + MLPA.

**Conclusion: the caller must be written.** It is the largest single addition
the pipeline has taken on — budget accordingly, and confirm the DRAGEN route is
genuinely unavailable before committing.

---

## Design

New M8 subworkflow, following the pattern the repo itself names the "contract
emitter model" (`docs/omnigen-additions-plan.md` §0): `subworkflows/smn_calling.nf`
→ a per-locus tool → a published TSV contract.

Four independent evidence channels, because no single one is sufficient:

1. **`HBA_DEPTH` — total α dosage.** Normalized depth over HBA1, HBA2, HBQ1 and
   the HS-40/MCS-R2 enhancer, against control windows. Gives α gene count
   (0–4). Reuse the `modules/traits/depth.nf` mosdepth `--by` pattern, but with
   its own scaling — `estimate_copies()` assumes `2 * ratio` and α is a 4-copy
   baseline.
2. **`HBA_PSV` — paralog discrimination.** Count reads at paralogous sequence
   variants distinguishing HBA1 from HBA2. This is the part depth cannot do and
   the part with no existing code in this repo.
3. **`HBA_JUNCTION` — allele identity.** The common deletions have fixed,
   published breakpoints (that is why gap-PCR works). Targeted split-read and
   discordant-pair search at the documented coordinates for `-α3.7`, `-α4.2`,
   `--SEA`, `--MED`, `--FIL`, `--THAI`, `-(α)20.5`, plus the `anti-3.7`
   triplication. Junction evidence gives the *allele*, not just a count, and
   sidesteps the segdup depth problem entirely.

4. **`HBA_SITES` — targeted non-deletional alleles.** Pileup at a fixed,
   version-controlled list of ~10–15 known pathogenic HBA1/HBA2 positions in a
   new `assets/hba_pathogenic_sites.tsv` (`chrom, pos, ref, alt, gene, hgvs,
   allele_name, clinical_note`). Minimum panel:

   | Allele | HGVS | Note |
   |---|---|---|
   | Hb Constant Spring | HBA2 c.427T>C | Most common non-deletional allele in SE Asia |
   | Hb Quong Sze | HBA2 c.377T>C | THAL2; major cause of non-deletional HbH in southern China |
   | Hb Adana | HBA2 c.179G>A | Severe; also occurs on HBA1 |
   | Hb Paksé, poly-A signal variants | — | Confirm final list against HbVar/IthaGenes before freezing |

   This is a pileup at known coordinates, **not variant discovery** — the same
   kind of operation as `HBA_JUNCTION`, which is why it does not breach the
   spec's SNP/indel-calling non-goal. Implement with `pysam` (already in the
   `svcaller/smncopynum` image) or `samtools mpileup`, reporting ref/alt depth
   and VAF per site.

   **This channel must consume channel 1's copy number.** VAF interpretation is
   copy-number dependent: on a `--SEA/αα` background the surviving HBA2 is
   hemizygous, so a real variant sits near 100% VAF, not 50%. Calling zygosity
   without the α-gene count produces wrong genotypes at exactly the compound
   heterozygotes that matter (`--SEA`/Quong Sze → HbH disease). Emit the raw VAF
   **and** the copy-number-aware zygosity call, never the latter alone.

`bin/alpha_globin.py` integrates the four into a contract TSV: α gene count,
called deletion alleles, targeted-site genotypes, an integrated genotype string
(e.g. `--SEA/αQSα`), per-channel confidence, and the scope manifest below. A
separate `bin/hba_report.py` builds the HTML card, mirroring `smn_report.py`
(which shares zero code with the trait callers).

### The safety contract — the non-negotiable part

The output must carry a **machine-readable** declaration of what was and was not
screened, so no downstream consumer can gate on `os.path.exists()` and render
"Clear". Minimum fields:

```
alpha_genes_called      0-4
deletion_alleles        e.g. --SEA/-a3.7
site_genotypes          e.g. HBA2:c.377T>C:het
genotype                integrated, e.g. --SEA/aQSa
screened                alpha_deletional,alpha_targeted_sites
site_panel_version      hba_pathogenic_sites.tsv @ <sha>
not_screened            beta_globin,alpha_nondeletional_outside_panel
interpretation_complete false
```

`site_panel_version` matters: "no pathogenic site found" means nothing without
knowing which sites were interrogated. Pin the panel file's hash in the contract
so a historical report can be re-read against the panel that actually ran.

`interpretation_complete=false` must be structurally impossible for SVcaller to
set true. The HTML card must render the limitation **inline**, not as a
footnote, and must never use the unqualified word "thalassemia" in a result
line. OmniGen must be updated to refuse to emit a thalassemia carrier statement
from this file alone — that is a coordinated change, not a SVcaller-only one.

### Guard wiring

`_require_nonempty()` in `bin/html_report.py:40-55` is generic but called on a
hardcoded four-item list (`:1408-1411`); the four CNV-trait TSVs get no guard at
all. A new caller inherits **nothing** automatically — the call site must be
added explicitly. Reuse the three-state model already established: *absent ≠
empty ≠ populated*, where a header-only file is a valid negative
(`bin/smn_report.py:60-73`).

Join with `remainder: true` + `?: file("NO_FILE")`, **not** an exact `.join()`.
`report.nf:163` uses an exact join for SMN, and per CLAUDE.md any meta-map key
divergence silently drops the whole sample's report.

---

## Filter changes required

These must land regardless, or a working caller's output dies downstream:

1. **Exempt the α-globin locus from the size and PON filters** — a region
   allowlist consulted by `modules/jasmine/merge.nf` before the `>10000` drop,
   so `--SEA`/`--MED`/`--FIL`/`--THAI` survive single-caller support.
2. **Fix `GNOMAD_SV_FILTER`'s `.get()` chain** (`annotsv/annotate.nf:94`) *and*
   exempt known-pathogenic common variants from AF-based dropping. Today
   deletions escape by accident; someone will fix that bug and silently break
   this. Add a regression test pinning the behaviour.
3. **Add α-globin regions to `assets/cnv_syndromes.tsv`** or the caller's output
   will never reach Clinical Findings.
4. Consider `bin/cnv_consensus.py:23-28` — `reciprocal_overlap()` divides by
   `min(len)`, so it is not true reciprocal overlap. Out of scope here, but it
   affects small-CNV behaviour generally.

---

## Validation

Two WGS thalassemia samples are on hand, and with the targeted-site channel
approved **both are now in scope** — one per half of the module. The fact that
one of them was *out* of scope under the original plan is what drove the scope
change, and remains the clearest argument in this document.

| Sample | Variant | Class | Channel that must detect it |
|---|---|---|---|
| **THAL1** | `--SEA`, HBA1+HBA2 deletion (~20 kb), het | Deletional | 1/3 — depth + junction |
| **THAL2** | Hb Quong Sze, HBA2 c.377T>C (p.Leu126Pro), het | Non-deletional point mutation | **4 — targeted sites** (in scope as of 2026-07-21) |

With the targeted-site channel approved, **both samples are now true positives**,
and each exercises a different half of the module. That is a much better
validation set than one positive plus one documented blind spot — but note the
two are *independent* tests: THAL1 says nothing about the site channel and THAL2
says nothing about deletion calling.

### Sample identity — labels were reversed (corrected 2026-07-22)

**The sample names as supplied are the opposite of what the reads contain.** This
was caught by inspecting the BAMs directly, before any code was written. Earlier
revisions of this plan had the truth set backwards; a correct caller validated
against those labels would have scored 0/2.

Evidence, from `/data/alvin/ref/THAL/*_30X.bwa.sortdup.bqsr.bam`
(~76 GB each, `samtools quickcheck` OK, full hg38, DNBSEQ/MGI reads —
`@RG PL:COMPLETE`, aligned with BGI MegaBOLT bwa):

**Depth**, 1 kb windows, normalized to a chr2:100,000,000-100,020,000 control:

| Region | THAL1 | THAL2 |
|---|---|---|
| HBZ, chr16:142–155 kb | 1.00 | 1.00 |
| chr16:164–186 kb (HBA1+HBA2) | **0.45–0.60** | ~1.00 |
| chr16:187 kb onward | 1.05 | 1.05 |

THAL1 carries a clean, sustained half-depth block from ~164 kb to ~186 kb that
spans HBA1+HBA2 and spares HBZ — the `--SEA` signature, heterozygous. THAL2 is
flat across the whole cluster: **no deletion**. (The near-zero window at 186 kb
and the dip at 172 kb appear in *both* samples — unmappable regions, not
variants. Do not let a caller report them.)

**Point mutation**, `samtools mpileup` over chr16:172,800-177,600:

`chr16:173,548 T>C, VAF 0.55, DP=22` — **present in THAL2 only**, absent from
THAL1. HBA2 spans chr16:172,876-173,710; walking its exon structure places c.377
at ≈chr16:173,547. THAL2 is the Hb Quong Sze heterozygote.

Two caveats on that coordinate, both to be closed before the site panel is
frozen (§ Open questions 5):

- ≈173,547 is arithmetic from exon boundaries, **not** an authoritative
  transcript alignment. It lands within ~1 bp and the substitution is T>C as
  expected, but confirm against the HBA2 annotation.
- The deletion is "consistent with `--SEA`" by size and position only. Exact
  breakpoints need split-read analysis — which channel 3 (`HBA_JUNCTION`) will
  do anyway, so this doubles as its first real test.

One encouraging signal: the Quong Sze VAF is 0.55, not the ~0.25 expected if
HBA1/HBA2 reads collapsed onto one paralog. Reads are segregating correctly at
that site, which is mild evidence the targeted-site channel is feasible at 30×.

### THAL1 also carries a beta-thalassemia allele (found 2026-07-22)

Running the newly derived HBB panel against THAL1's germline VCF hit
immediately:

```
IVS-II-654  c.316-197  chr11:5225923  G>A  GT=0/1  AD=9,10  (VAF 0.53)
```

IVS-II-654 is a common Chinese β⁺-thalassemia allele. **THAL1 is therefore a
double heterozygote — `--SEA` α-thal trait *and* β-thal trait.** Nothing in the
purchase description said so, and nothing in this plan assumed it.

Two consequences:

- As a validation sample THAL1 is *more* valuable than assumed: it exercises the
  α deletional channel and the β panel simultaneously, and it is the only sample
  on hand that tests whether the two screens stay independent and neither
  suppresses the other's finding.
- It sharpens the false-reassurance risk. A report that measures only α would
  describe THAL1 by its `--SEA` allele alone, which is an incomplete and
  clinically misleading picture of that individual.

Caveats: this rests on **one VCF record at depth 19**, not orthogonally
confirmed, and the allele is intronic (pathogenic via cryptic splicing) so it
cannot be read off the protein consequence. Confirm from the reads and, ideally,
by the supplier's own assay before entering it in the truth table.

**Action:** confirm with the supplier whether the mix-up is in the filenames or
in our notes, and record the resolution in `validation/thal_truth_table.tsv`.
Until then treat the *sequence* as authoritative and the *filename* as suspect.
Re-verify identity for MITO too — a third BAM,
`MITO_30X.bwa.sortdup.bqsr.bam` (78 GB), sits in the same directory and its
provenance has not been checked.

> **Nomenclature — resolved.** The allele is **HBA2:c.377T>C** (ClinVar,
> p.Leu126Pro); the earlier `c.377C>T` was a typo. Use `T>C` when the site panel
> is written — a targeted assay keyed on the wrong base calls nothing.

Hb Quong Sze is not an edge case: it is one of the major alleles causing
**non-deletional HbH disease in southern Chinese populations**, and the
hyper-unstable α chain it produces is degraded post-translationally. A
`--SEA` × Quong Sze couple — precisely THAL1 × THAL2 — carries a 25% risk of
HbH disease. Under deletional-only scope the module would have detected one
partner and reported the other as having no α-globin deletions; **that pairing
is why the targeted-site channel was added.**

### The validation set does not support the stated intended use

This needs saying plainly. **Two purchased samples are the whole positive set,
and no more are available.** They cover one deletion allele (`--SEA`) and one
site allele (Quong Sze). That means at launch:

- `-α3.7` — the most common α-thal deletion worldwide — **unvalidated**
- `-α4.2`, `--MED`, `--FIL`, `--THAI`, `anti-3.7` triplication — **unvalidated**
- Hb Constant Spring — the most common SE Asian non-deletional allele — **unvalidated**
- Copy-number-aware zygosity on a deletion background — **unvalidated** (neither
  sample is a compound heterozygote)

Sensitivity measured on n=1 per channel cannot be generalised to a panel. For
**couples carrier screening**, where a false negative can lead to an affected
pregnancy, that is not an adequate evidence base. Two honest options:

1. **Label it research-use**, exclude it from reproductive decision-making, and
   state the unvalidated alleles explicitly in the contract and the report; or
2. **Strengthen the evidence base** before clinical use — orthogonally type the
   GIAB samples already on hand (gap-PCR on HG005/6/7, the Han Chinese trio,
   where `-α3.7`/`--SEA` carriers are plausible) to add real controls cheaply,
   and/or source a `-α3.7` carrier and a compound heterozygote.

In-silico validation — simulated reads carrying known breakpoints — should be
built regardless. It validates the *detection logic* across all target alleles,
which is genuinely useful, but it does not measure real-world sensitivity and
must not be reported as if it did.

Consequences for the validation plan:

1. **THAL1 is the only positive control for deletional calling.** One positive,
   covering one of seven-plus target alleles, is thin. Sensitivity for `-α3.7`,
   `-α4.2`, `--MED`, `--FIL`, `--THAI` and the `anti-3.7` triplication is
   **unmeasured**, and must be declared as such rather than assumed. Sourcing at
   least a `-α3.7` carrier is a prerequisite for clinical use — it is the most
   common allele worldwide.
2. **THAL2 is the only positive control for the site channel**, and covers one
   of ~10–15 panel positions. Hb Constant Spring in particular is the highest-
   frequency SE Asian non-deletional allele and has **no validation sample** —
   its sensitivity will be unmeasured at launch. Say so explicitly in the
   contract rather than implying panel-wide sensitivity from one allele.
3. **THAL2 doubles as the module's hardest safety test.** It must (a) yield a
   correct 4-α-gene *deletional* result while (b) still flagging a pathogenic
   site, and (c) never produce a string a clinician could read as "not a
   thalassemia carrier" — β-thal remains untested for it. Automate this: a
   carrier whose only finding lies outside the deletional channel must never
   yield a reassuring report.
4. **Test the copy-number-aware zygosity path directly.** Neither sample is a
   compound heterozygote, so the `--SEA`-background hemizygous-VAF logic — the
   most error-prone part of channel 4 — has **no real sample to validate it**.
   Cover it with a synthetic/downsampled BAM fixture, and treat a genuine
   compound-heterozygote sample as the highest-value future acquisition.
3. Build `validation/thal_truth_table.tsv` mirroring `validation/smn_truth_table.tsv`
   (`sample, alpha_genes, alleles, class, method, source`), recording variant
   **class** so non-deletional entries are explicitly out-of-scope-by-design
   rather than looking like misses.
4. **Write a script that actually reads it.** `smn_truth_table.tsv` is currently
   decorative — nothing in `bin/`, `tests/`, or `validation/` consumes it, and
   SMN accuracy rests on a human eyeballing TSVs. Do not repeat that.
   `validation/thal_benchmark.sh` producing a concordance table is part of the
   deliverable, not a follow-up.
5. Report per-allele sensitivity/specificity and the α-gene-count confusion
   matrix, with unmeasured alleles shown as unmeasured. Gate clinical use on it.

---

## Scope decision — targeted sites included (approved 2026-07-21)

Scope was initially set to α-deletional only, on the reasoning that
non-deletional α requires SNV calling, an explicit spec non-goal. THAL2 (Hb
Quong Sze) forced a narrow revisit, and the targeted-site channel was
**approved**.

The non-goal is respected because this is a pileup at a fixed, version-controlled
coordinate list — the same kind of operation as the `HBA_JUNCTION` breakpoint
check, not variant discovery. No caller is run genome-wide; nothing is
discovered that was not already named in `assets/hba_pathogenic_sites.tsv`.

Effect: α-thal coverage rises from ~80–90% (deletional) to ~95%+, THAL2 becomes
a true positive rather than a documented blind spot, and the false-reassurance
risk is reduced structurally rather than managed purely through report wording.

**It does not solve β-thalassemia**, which remains ~82% point mutations across
~280 known alleles and out of scope *for SVcaller*. The classification split
therefore stands unchanged: SVcaller measures α, something else owns the
couple-level verdict. β is not blocked on this pipeline at all — see
§ β-globin below.

### Where this lives: measurement in SVcaller, screen ownership in OmniGen

Reconsidered directly, because the objection is a fair one: none of this work
reuses existing SVcaller output, so why not build it all in OmniGen, or
standalone?

**Three findings settle it.**

**1. OmniGen already owns the thalassemia claim.** This is not a new screen to
site — `carrier_panel.tsv:13-14,31` already lists Alpha Thalassemia and Beta
Chain-Related Hemoglobinopathy as screened conditions. The goal "OmniGen owns
the total thalassemia screen" is already true on paper. What is missing is the
capability behind the α half. Nothing needs relocating for OmniGen to own the
screen; it does.

**2. β-thalassemia is probably already largely working, and needs no new
measurement.** β-thal is ~82% point mutations in HBB, and **HBB is not in a
segmental duplication** — so a standard germline VCF plus gene-panel lookup,
exactly what OmniGen already does, should capture most β alleles. The gap in the
"total screen" is therefore **specifically α, not β**. Verify allele-level
coverage against a β-thal panel (IVS-I-110, cd39, IVS-I-1, IVS-II-654, cd41/42,
−28, cd17), but the mechanism exists and costs nothing new.

**3. α-thalassemia cannot be done in OmniGen's shape at all.** 80–90% is
deletional and invisible to any VCF lookup; the rest sits in a segdup where
standard SNV calls are unreliable. Both need the BAM, plus copy-number-aware
interpretation. And OmniGen is architecturally a *detector of what ran*:
`entrypoint._resolves()` (`:83`) checks which upstream artifacts exist and
`detect()` (`:88`) infers BAM- vs VCF-derived from that. If OmniGen starts
consuming BAMs itself it becomes both the detector and a producer of coverage,
and that model stops meaning anything. It would also need containers, resource
tiers, retry and per-sample compute orchestration — all of which SVcaller has
and OmniGen deliberately does not.

**Against standalone:** a third pipeline to version, deploy and operate for ~15
positions and a handful of breakpoints, which would still need everything
SVcaller already provides. The organization rules just adopted favour fewer
moving parts.

**This is not duplicated work — it is the SMA split, unchanged.**
SMNCopyNumberCaller runs inside SVcaller and emits `smn.tsv`; OmniGen presents
SMA carrier status. Nobody regards that as split ownership. α-globin is
structurally identical. The rule "SVcaller does SV/CNV, OmniGen parses results"
is *preserved*, not broken: a `--SEA` allele is a 20 kb deletion, which is
squarely SV/CNV work.

**Where the objection is right, and what changes because of it:** thalassemia
*interpretation* — HbH/Bart's/trait classification, couple-level risk — must not
be duplicated in SVcaller's HTML. SVcaller's card should be thin and factual
(gene count, alleles, site genotypes, scope), and OmniGen owns the clinical
narrative. The plan is amended accordingly.

### Why the site check specifically cannot be split off

Its **interpretation depends on the α-globin copy number the same module
produces**:

> On a `--SEA/αα` background, both α genes on one chromosome are gone. An HBA2
> variant on the surviving chromosome is **hemizygous** and appears at ~100% VAF,
> not ~50%. A caller without the copy-number context reads that as homozygous —
> a materially wrong genotype, and precisely the `--SEA` × Quong Sze compound
> heterozygote (HbH disease) that matters most clinically.

Copy number and point mutation are one analysis at this locus. Splitting them
across repos means shipping the VAF without the denominator.

The architectural facts confirm it:

- **OmniGen does not analyse BAMs per sample.** `prototype/entrypoint.py:83`
  `_resolves(sample, key)` detects which upstream *artifacts exist*; `detect()`
  (`:88`) infers whether a report is BAM- or VCF-derived from that. Its BAM
  references elsewhere are reference-DB builders, not per-sample callers. Adding
  a pileup engine would turn the reporting layer into a caller.
- **OmniGen aggregates several pipelines** (SVcaller, pgx-suite, poly-suite) and
  computes coverage "from what ACTUALLY resolves" (`entrypoint.py:41-43`, with
  HG003 as the worked example). It is a consumer of contracts.
- pgx-suite is the other BAM-consuming pipeline, but HBA is not pharmacogenomics.

**So: the fourth evidence channel goes in SVcaller's α-globin module, sharing
state with the depth, PSV and junction channels.**

### What OmniGen should own instead

OmniGen already has exactly the machinery this plan's safety contract needs —
better than anything hand-rolled in SVcaller's HTML:

- `entrypoint.not_covered(sample, entry)` (`prototype/entrypoint.py:106`) returns
  the tiers an input could not test, rendered as an explicit amber banner
  (`omnigen_dtc.py:1226-1244`): *"their absence is **not a negative result**.
  Nothing was tested, so nothing can be ruled out."*
- A separate red banner names screens that **failed** to run (`:1215-1226`),
  written in response to the same incident this plan cites.

The β-thalassemia and beyond-panel non-deletional α limitations should be
registered as first-class **never-covered tiers** in that model, so the existing
banner names them automatically. That is strictly better than prose in the
SVcaller report, because it survives the report being regenerated, and it is
enforced in the layer the couple actually reads.

Division of responsibility:

| Layer | Owns | Status |
|---|---|---|
| **SVcaller** | **Measures α-globin only**: deletion alleles, gene count, targeted HBA1/HBA2 site genotypes, copy-number-aware zygosity. One contract file with scope metadata. Thin factual HTML card. **Never says "thalassemia"; no HbH/Bart's/trait classification; no couple-level risk.** | To build |
| **OmniGen** | **Owns the total thalassemia screen.** Rewires its existing Alpha Thalassemia panel row to consume SVcaller's α contract instead of VCF gene lookup; keeps β-thal on the existing HBB/VCF path; registers beyond-panel α as never-covered; owns all interpretation and any couple-level statement. | Claim already exists (`carrier_panel.tsv:13-14`); α capability missing |
| **Existing SNV pipeline** | Supplies the germline VCF that OmniGen already reads for HBB. | Verify β allele coverage |

The single most important integration point: **OmniGen's Alpha Thalassemia row
must stop being fed by HBA1/HBA2 gene lookup and start being fed by the α
contract.** Until that rewire happens, the row is answering a question it cannot
answer — see P0 above.

---

## Verification

- `pytest tests/` — new `tests/test_alpha_globin.py` covering the pure functions:
  dosage → gene count, junction evidence → allele, site pileup → genotype,
  **copy-number-aware zygosity** (the `--SEA`-background hemizygous case), and
  integration conflicts (e.g. a site variant on a gene the depth channel says is
  deleted). Plus a contract-schema test like
  `test_cnv_traits.py::test_contract_schemas`.
- Add an entry to `tests/test_no_empty_placeholders.py` mirroring
  `test_smn_caller_hard_fails_on_missing_output` (`:142`) — asserts no `touch`
  fallback, presence of `exit 1`, and a header-only rejection. A header-only
  site table is a **valid negative** (no pathogenic site found) and must not be
  treated as failure; an absent one must.
- Synthetic BAM fixture for the hemizygous-VAF path, since neither THAL1 nor
  THAL2 is a compound heterozygote.
- Run THAL1 + THAL2 end-to-end; compare against the truth table via the new
  benchmark script. THAL1 must give 2 α genes and `--SEA`; THAL2 must give 4 α
  genes **and** a heterozygous HBA2 c.377T>C. Require both before any
  carrier-screening use.
- Re-run HG001–HG007 and confirm the α locus now yields calls where expected,
  that the Han Chinese trio is consistent with population frequency, and that
  all seven are negative across the site panel.
- Confirm the HTML card renders the not-screened declaration and
  `site_panel_version` inline; that a deliberately emptied contract fails the run
  rather than rendering "Clear"; and that THAL2's report contains no string a
  clinician could read as "not a thalassemia carrier".

---

## β-globin — the VCF path (added 2026-07-22)

**β-thalassemia never needed SVcaller.** Earlier revisions treated β as "out of
scope, needs its own plan", which implied a missing capability. It isn't one —
β is a *different input*, not a harder problem, and the input already exists.

HBB (chr11p15.4) is **not in a segmental duplication**. Ordinary short-read SNV
calling is reliable there, so the ~82% of β-thal that is point mutations is
readable from a plain germline VCF: no depth, no PSV discrimination, no reads.
That places it on OmniGen's side of the hard boundary in `entrypoint.py:16`, not
SVcaller's. OmniGen already claims Beta Chain-Related Hemoglobinopathy → HBB in
`carrier_panel.tsv:31` and already parses VCFs.

### The VCFs exist and are on disk

Confirmed 2026-07-22 in `/data/alvin/ref/THAL/`:

```
THAL1_30X.bwa.sortdup.bqsr.hc4.vcf.gz   287 MB  + .tbi
THAL2_30X.bwa.sortdup.bqsr.hc4.vcf.gz   288 MB  + .tbi
```

`hc4` = GATK HaplotypeCaller 4.x, same MegaBOLT delivery as the BAMs. Equivalent
`.hc.vcf.gz` files already exist for HG001–HG007 in `/data/alvin/ref/GIAB/`, so
the GIAB set doubles as a negative control cohort for the β panel at no cost.

**HBB is well covered.** `tabix chr11:5225000-5232000` returns 17 records for
THAL1 and 10 for THAL2, with sound QUAL (191–1321), sensible genotypes and
balanced AD. Callability at the locus is not in question.

### What is actually missing

Not a caller — an **annotated β-thal allele list**, and the means to apply it:

- **No ClinVar VCF on disk.**
- **VEP is installed** (`/data/alvin/ensembl-vep`) but **has no cache**, so it
  cannot annotate offline as-is.

Until one of those is fixed, the 17/10 HBB records above are uninterpreted. They
are *presumed* common polymorphisms — neither sample was purchased as a β
carrier — but nothing has been checked, and "presumed benign" must not be
reported as "β screened, negative".

The deliverable is therefore a curated `hbb_pathogenic_sites.tsv` (CD41-42
`-TTCT`, IVS-II-654 C>T, −28 A>G, CD17 A>T, CD71-72 `+A`, CD26/HbE, HbS, plus
population-appropriate additions), version-pinned exactly like the α site panel.
β alleles cluster strongly by population, so ~20–30 entries cover >90% in most
groups — the same fixed-coordinate lookup pattern, applied to a VCF instead of a
pileup.

### Neither THAL sample validates the β path

Both are α samples. There is **no β-thal positive control on hand**, so β
sensitivity will be entirely unmeasured at launch — the same n=0 problem the α
channels have at n=1. This must appear in the `not_screened` contract, not be
quietly implied as covered.

### Bonus finding — the VCF also calls the α site

`chr16:173,548 T>C` appears in **THAL2's** VCF at `GT=0/1, QUAL=330, AD=12,13`
(VAF 0.52) and is correctly absent from THAL1 — matching the pileup result
(0.55) closely.

This does **not** retire channel 4. Three reasons the pileup stays:

1. The VCF cannot supply α copy number, which channel 4's zygosity call requires.
2. `GT=0/1` is *correct here only because THAL2 has no deletion*. On a
   `--SEA/αQSα` background the surviving HBA2 is hemizygous and the VCF genotype
   would be wrong — precisely the compound heterozygote that matters clinically,
   and precisely the case neither sample tests.
3. One concordant site in a segdup is not evidence of general reliability there.

It is genuinely useful as a **cross-check**: agreement between pileup and VCF at
panel sites is a cheap correctness signal, and disagreement is a strong hint the
copy-number-aware logic has a bug. Wire it into `validation/thal_benchmark.sh`.

### Revised ownership

| Layer | Input | Home |
|---|---|---|
| α deletional (gene count, alleles) | BAM — depth + junction | **SVcaller** (to build) |
| α non-deletional, named sites | BAM — targeted pileup | **SVcaller** (to build) |
| β point mutations (~82% of β) | **germline VCF** | **OmniGen** (panel content) |
| β large deletions (~6% of β alleles) | BAM — SV ensemble | SVcaller, incidental, unvalidated |
| Classification, couple risk | all of the above | **OmniGen** |

β work is therefore unblocked and independent of this plan — it needs no
Nextflow change, and can proceed in parallel.

---

## Phase 0 — shared interface, completed 2026-07-22

Built before fanning out to independent SVcaller and OmniGen tracks, because
both consume the same artifacts and would otherwise each invent their own.

| Artifact | What it settles |
|---|---|
| `docs/contracts/alpha_globin_contract.md` | file path, format, exact column list, both open decisions |
| `validation/examples/SAMPLE.alpha_globin.tsv` | canonical fixture both sides test against |
| `assets/hba_pathogenic_sites.tsv` | 4 α sites |
| `assets/hbb_pathogenic_sites.tsv` | 12 β sites |
| `assets/hba_segments.bed` | 5 diagnostic segments, derived from gene models |
| `assets/hba_deletion_alleles.tsv` | 7 deletion/triplication alleles by CN signature |
| `bin/hgvs_map.py` | HGVS c. → GRCh38, with self-test |
| `bin/make_globin_panels.py` | regenerates both site panels; refuses to emit a bad row |
| `bin/make_hba_deletion_alleles.py` | regenerates the segments + allele signatures |
| `tests/test_globin_panels.py` | 28 tests; suite now 107 |

**The β site panel turned out to be largely redundant.** OmniGen's evidence DB
already holds 431 HBB P/LP variants, and **11 of the 12 curated sites are
already there with canonical rsIDs** (HbE `rs33950507`, IVS-I-5 `rs33915217`,
IVS-II-654 `rs34451549`, CD39 `rs11549407`, CD41-42 `rs36029927`, −28
`rs33931746`, …); only CD71-72 is absent. So OmniGen's β path works as-is and
that track is much smaller than first scoped. The panel's residual value is that
it records *which* alleles matter by population — which a flat 431-row list does
not — and its coordinates are what found IVS-II-654 in THAL1. Those rsIDs also
confirm the derived coordinates a third time, independently of population
structure and the FASTA check.

### Deletion alleles are defined by signature, not by breakpoint

The site panel covers only point mutations. The deletional alleles — 80–90% of
α-thal, and the module's primary job — had **no definition anywhere**, which
blocked channels 2 and 3.

They are not derivable the way sites are: α-cluster NAHR breakpoints sit inside
near-identical homology boxes, are published against varying builds, and are not
even single-valued for a given allele. Hand-typing GRCh38 coordinates would be
precisely the failure the site generator exists to prevent.

So `hba_deletion_alleles.tsv` defines each allele by **which diagnostic segments
it removes** — what a depth caller measures anyway — over the five segments in
`hba_segments.bed`, whose boundaries *are* derived from the RefSeq models.
`approx_size` is documentary and a test asserts it can never be mistaken for a
coordinate.

**Two degenerate groups fall out, and they change the contract:**

```
--SEA | --MED     remove HBA1+HBA2, spare HBZ
--FIL | --THAI    remove HBA1+HBA2 and HBZ
```

Depth cannot separate members of a group. A caller emitting `--SEA` from depth
alone is inventing precision — and choosing `--SEA` over `--MED` is a population
inference, which would look correct in SE Asian samples and be wrong elsewhere.
The contract now requires the group form (`--SEA|--MED/aa`), collapsing only on
a junction read or an extent that excludes the alternative.

`anti-3.7` (triplication) is included deliberately: a caller written to look
only for losses would never report it, and it modifies β-thal severity.

**Two decisions that were genuinely ambiguous** (details in the contract):

1. **Discovery is by path convention**, not manifest key. OmniGen does it both
   ways today — SMN uses `config.resolve(SAMPLE,"smn")`, CNV traits use a
   hardcoded `${SV}/results/...` path — so "follow the existing pattern" would
   have produced two incompatible integrations. `config.py` reserves manifest
   keys for *irregularly-named* paths; SVcaller's are regular.
2. **α-globin is its own subworkflow**, not a fifth CNV trait caller, because
   the junction and pileup channels do not fit `TRAIT_DEPTH`'s shape — but it
   reuses `bin/cnv_traits_common.py` for normalized depth rather than
   duplicating it. An earlier note in this plan said "fifth trait caller"; that
   was half right and is superseded here.

**Coordinates are derived, never hand-typed.** `hgvs_map.py` reproduces six
coordinates established independently of it (HbS, HbE, CD41-42, CD17 from 1000G
population structure; Hb Quong Sze from THAL2's reads; IVS-II-654 cross-checked
two ways), and `make_globin_panels.py` verifies every row against the reference
FASTA — HBB is minus-strand, so a panel that forgets to complement calls
nothing while looking healthy.

The guards earned their keep immediately: generation **failed** on the `-28`
allele, because that name is numbered from the CAP site, not the ATG. HBB has a
50 nt 5′ UTR, so `-28` is `c.-78` (chr11:5227099, inside the `CATAAAA` TATA box).
Hand-curation would have shipped a site that silently never matched.

---

## Open questions

**Resolved 2026-07-21:** HGVS is `c.377T>C` (typo corrected). DRAGEN unavailable
— caller must be written. Placement: measurement in SVcaller, screen ownership
in OmniGen (§ Where this lives). Site check in scope.

**Resolved 2026-07-22:** Sample identity — the supplied THAL1/THAL2 labels are
reversed. THAL1 is the `--SEA` deletion; THAL2 is Hb Quong Sze. Verified from
the reads (§ Sample identity). Both BAMs are on disk, complete, and aligned to
full hg38. **No FASTQ exists for either sample** — BAM is the only input, so
both always run `FILTER_CHROMS`. Germline VCFs (`*.hc4.vcf.gz`) for both
samples were copied the same day and HBB is confirmed callable, which resolves
the old "is a germline VCF routinely produced?" question and unblocks β
independently of this plan (§ β-globin).

Still open:

1. ~~**P0 — does OmniGen currently report Alpha Thalassemia as negative?**~~
   **RESOLVED 2026-07-22: yes, it does.** Confirmed by running
   `report_carrier_panel.py`. Producer traced to `report_carrier_panel.py:111-115`.
   See top — the remaining action is the mitigation, not the investigation.
2. **Does OmniGen's HBB panel row detect the common β alleles?** The *input* is
   now confirmed present and callable (§ β-globin); what remains is whether
   OmniGen's panel content actually recognises the common β-thal alleles, and
   what it renders when it finds none. Blocked on annotation — no ClinVar VCF on
   disk and VEP has no cache.
3. ~~**Research-use or clinical?**~~ **RESOLVED 2026-07-22: research-use.** The
   consumer is OmniGen, a **RUO / educational DTC** report — its README states
   *"Nothing it produces is medical advice or a diagnosis"* and everything is
   graded A–D with caveats. This substantially relaxes § The validation set does
   not support the stated intended use: that section was written against
   "couples carrier screening" driving reproductive decisions, and n=1 per
   channel *is* acceptable for a graded educational report that states its own
   confidence.

   **It does not relax the false-negative obligation, and DTC sharpens it.**
   Two different risks were being conflated:

   - *insufficient evidence for a clinical claim* — relaxed; grade it C/D and say so.
   - *reporting "Alpha Thalassemia — negative" when the caller structurally
     cannot see the deletion* — *not* relaxed. That is a false statement, not an
     ungraded one, and in DTC there is no clinician mediating it. P0 stands.

   See `docs/what-a-30x-bam-yields.md` for the full scoping note.
4. **Freeze the site panels' CONTENT.** Coordinates are settled (phase 0 —
   derived and FASTA-checked), but the *allele list* is curated and incomplete:
   4 α sites, 12 β sites. Review against HbVar/IthaGenes and decide what a
   population-appropriate panel must contain. Reviewers should check the HGVS
   and allele names in `bin/make_globin_panels.py`, never the coordinates.
5. **Get annotation working.** No ClinVar VCF on disk; VEP is installed with no
   cache. No longer blocks the c. → genomic mapping (phase 0 solved that from
   the AnnotSV RefSeq models), but still blocks classifying anything *outside*
   the curated panels — i.e. any variant we did not already name.
6. **Are THAL1/THAL2's orthogonal results gap-PCR, MLPA, or Sanger?** Each
   validates a different channel; Quong Sze implies Sanger or a targeted assay,
   which does not validate deletion calling at all. **Now also the route to
   settling the label swap** — ask the supplier for the assay reports keyed to
   the tube IDs, not the filenames.
7. Does OmniGen's `not_covered` model, keyed on entry point (VCF vs BAM), need
   extending to express a tier that is **never** coverable by any input?
8. **What is `MITO_30X`?** A third 78 GB BAM in `/data/alvin/ref/THAL/` whose
   identity and purpose are unverified. Given the THAL1/THAL2 swap, do not
   assume its label is right either.
9. **Source a β-thal positive control.** β sensitivity is currently n=0. Same
   acquisition conversation as the `-α3.7` carrier and the compound
   heterozygote.
