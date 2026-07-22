# What a 30× WGS BAM actually yields — scoping note for OmniGen (RUO/DTC)

Written 2026-07-22, prompted by the thalassemia work. Belongs beside
`OmniGen/docs/blueprint-expansion.md`; kept here because that is where the
caller-side evidence lives. Move it if it reads better there.

## The motivation, as the project itself states it

- `OmniGen/README.md` — *"A graded, literature-backed interpretation layer for
  whole-genome sequencing."* Design principle: **grade, don't re-call**. OmniGen
  consumes specialist-caller output and applies uniform A–D grading.
- Scope badge and prose: **educational, non-clinical, RUO.** *"Nothing it
  produces is medical advice or a diagnosis."*
- `docs/omnigen-design-doc.md` — BAM tier is the higher-fee product because it
  is genuinely more compute and more value; consumer ARPU is capped, so
  **virality is the growth engine**; the social layer is deliberately confined
  to ancestry + benign traits, health stays private.
- `docs/input-coverage.md` — the containment hierarchy, already measured:

  ```
  VCF  8/13 reports  ~5–15 min   <$1
  BAM 13/13 reports  ~4–5 h      ~$14 on-demand / ~$6 spot
  FASTQ 13/13        ~8–10 h     ~$25 / ~$11    (adds nothing over BAM but alignment)
  ```

  The five BAM-only tiers today: **PGx star alleles, HLA, SMA/SMN, STR repeat
  expansions, SV/CNV.**

**So the BAM tier's entire commercial justification is the read-only tiers.**
Anything that widens that gap widens the reason to buy the more expensive
product. That is the lens for everything below.

## What this changes about the thalassemia plan

The α-globin plan was written against *"intended use: couples carrier
screening"*, and its harshest section — *"The validation set does not support
the stated intended use"* — argues n=2 cannot support reproductive decisions.

**RUO/DTC answers that.** Open question 3 (research-use vs clinical) resolves to
research-use, and the validation-sufficiency objection largely dissolves: n=1
per channel is acceptable for a graded, educational report that states its
confidence.

**It does not dissolve the false-negative obligation, and arguably sharpens it.**
Those are different risks:

| Risk | Relaxed by RUO? |
|---|---|
| Insufficient evidence to support a clinical claim | **Yes** — grade it C/D and say so |
| Reporting "Alpha Thalassemia — negative" when the caller structurally cannot see the deletion | **No** — that is a false statement, not an ungraded one |

In DTC it is worse than in clinic: there is no clinician mediating, and a
consumer reads "negative" as reassurance regardless of an RUO badge. The P0 item
stands unchanged and is the highest-priority thing in the plan.

## Extractable from a 30× BAM but not currently produced

Checked against `blueprint-expansion.md` Part 1 and the docs set — items with
zero mentions anywhere are marked **new**.

### Tier 1 — reuses `cnv_traits.nf` wholesale (targeted depth + a ~50-line interpreter)

The existing four trait callers (AMY1, GSTM1/T1, LPA KIV-2, RHD) already prove
this shape. Marginal architecture ≈ zero; each is a BED region plus an
interpreter.

| Signal | Why it fits DTC | Status |
|---|---|---|
| **α-globin / thalassemia** | carrier interest, high SE-Asian prevalence | in progress |
| **C4A/C4B copy number** | complement CNV; strong schizophrenia association (Sekar 2016), also SLE | **new** |
| **UGT2B17 deletion** | very common in East Asians; testosterone metabolism + drug glucuronidation — the "sports/doping" angle DTC loves | **new** |
| **Haptoglobin Hp1/Hp2** | CNV with cardiovascular/diabetes interaction; popular nutrition-DTC item | **new** |
| **Sex-chromosome karyotype** (XXY/XYY/XO) | listed "Medium" in Part 1, but is nearly free from existing MOSDEPTH output | listed, not built |

### Tier 2 — new method, high consumer appeal

- **mtDNA heteroplasmy** — at 30× nuclear, chrM sits at roughly 1000–5000×, so
  heteroplasmy down to ~1–2% is measurable. `entrypoint.py:37` already flags
  heteroplasmy as reads-only and explicitly distinct from the haplogroup call,
  but nothing produces it. `MITO_30X.bwa.sortdup.bqsr.bam` is already on disk.
  Strongest single candidate: cheap, genuinely BAM-only, unambiguously novel
  versus a VCF product.
- **Telomere length** (TelSeq-style, telomeric-repeat read fraction) — a
  headline DTC item ("biological age"). Cheap. Honest grade is D: noisy,
  batch-sensitive, and not clinically meaningful at the individual level. The
  A–D framework can carry it truthfully where a competitor's cannot.
- **Viral reads / integration** (HPV, HBV, EBV) from the unmapped fraction —
  reads that are currently discarded. Novel for consumer reports and cheap.

### Tier 3 — real work, flagged for completeness

- **HLA class II** (DQ2/DQ8 → celiac, DQB1*06:02 → narcolepsy) — already
  identified as a gap; needs a class-II typer in pgx-suite.
- **KIR haplotypes** — structural/copy-number immune locus.
- **Mosaic chromosomal alterations / clonal haematopoiesis** from BAF+LRR —
  30× is marginal, and it is a health finding that may not belong in a DTC
  product at all.

## The honest ceiling — what 30× short reads cannot do

Worth stating plainly, because competitors advertise some of it:

- **Methylation / "epigenetic age"** — impossible from standard WGS. Needs
  bisulfite or nanopore. This is the single biggest DTC claim the product
  structurally cannot make.
- **Phasing beyond fragment length** — cis vs trans is often unresolvable
  without parents. Directly seen this session: calling HG02379's HbE/CD41-42
  compound heterozygote in *trans* rested on 4 informative fragments, and the
  two sites are only 181 bp apart. A pair 5 kb apart would have been
  unresolvable.
- **Large repeat expansions** — ExpansionHunter estimates beyond read length; it
  does not measure a full expansion.
- **Somatic mosaicism below ~5%** — depth-limited.
- **De novo assembly of hyperpolymorphic/structurally complex loci.**

## Suggested order

1. **P0** — confirm whether α-thal currently reports negative. Independent of everything, ~1 h.
2. **mtDNA heteroplasmy** — best novelty-per-effort, sample already on disk.
3. **C4A/C4B + UGT2B17 + haptoglobin** — one batch through the existing trait machinery.
4. **Sex-chromosome karyotype** — nearly free from data already produced.
5. Telomere length and viral reads as grade-D/novelty items once the above land.
