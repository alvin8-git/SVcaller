# WES (whole-exome capture) support — limitations

**Verdict: do NOT run SVcaller on a whole-exome (WES / hybrid-capture) BAM.**
The pipeline is a **human-WGS** pipeline (`nextflow.config` manifest:
`description = 'SV/CNV/SMN1-SMN2 pipeline for human WGS'`). On a WES BAM it
either **aborts at the coverage QC gate**, or — if that gate is bypassed —
emits **silently-wrong** depth/CNV results and **misses most SV breakpoints**.
There is no capture-BED code path that makes any of the depth/SV/QC logic
exome-aware.

Why WES breaks the assumptions: exome capture covers only exons (~1–2 % of the
genome). Off-target regions (introns, intergenic) sit at ~0 depth. Every
read-depth step here assumes a genome-wide / fixed-window baseline that has no
valid analogue on WES, and SV breakpoints almost always fall in introns that
capture does not cover.

---

## Per-module verdict

| Module / output | File:line (evidence) | Runs on WES? | Result validity |
|---|---|---|---|
| **Coverage QC gate** (mosdepth `--by 50000`, mean-depth threshold) | `modules/mosdepth/coverage.nf:24,30,33-37`; threshold `min_depth = 25` `nextflow.config:18`; called in `subworkflows/preprocess.nf:50` | **ABORTS** (`exit 1`) | Halts the whole run before any caller. Genome-wide mean depth on WES ≪ 25×. |
| **α-globin / HBA depth** (thalassemia) | `bin/hba_depth.py:286,289` → `control_baseline()` on `CTRL_*` | Runs only if gate bypassed | **SILENTLY WRONG / no-call** (CTRL baseline ≈ 0) |
| **AMY1 copy number** | `bin/amy1_cn.py:21` → `control_baseline()` | if gate bypassed | **SILENTLY WRONG / no-call** |
| **LPA / KIV-2** | `bin/lpa_kiv2.py` → `control_baseline()` | if gate bypassed | **SILENTLY WRONG / no-call** |
| **RHD status** | `bin/rh_status.py` → `control_baseline()` | if gate bypassed | **SILENTLY WRONG / no-call** |
| **GSTM1 / GSTT1 null** | `bin/gst_null.py` → `control_baseline()` | if gate bypassed | **SILENTLY WRONG / no-call** |
| **GATK gCNV** | `modules/gatk/gcnv_call.nf:28` (`-L intervals`); PON `workflows/pon_build.nf` | if gate bypassed | **SILENTLY WRONG** — WGS-binned PON; WES on/off-target bias not modeled |
| **CNVpytor** | `modules/cnvpytor/call.nf` | if gate bypassed | **SILENTLY WRONG** — genome-wide RD histogram baseline collapses off-target |
| **SMN1/SMN2 caller** | `subworkflows/smn_calling.nf`, `modules/smn_caller` | if gate bypassed | **UNRELIABLE** — SMNCopyNumberCaller normalization is WGS-trained |
| **Manta** (SV) | `modules/manta/call.nf` | Runs | **DEGRADED** — whole-genome, no capture BED; intronic breakpoints uncaptured → missed |
| **DELLY** (SV) | `modules/delly/call.nf` | Runs | **DEGRADED** — same; misses intronic breakpoints |
| **GRIDSS** (SV) | `modules/gridss/call.nf` | Runs | **DEGRADED** — same |
| **Scramble** (MEI) | `modules/scramble/call.nf` | Runs | **DEGRADED** — same |
| **MELT / STRling / ExpansionHunter** | `modules/melt`, `modules/strling`, `modules/expansionhunter` | Runs | **DEGRADED** — depth/repeat models assume WGS coverage |

Legend: **ABORTS** = process exits non-zero and stops the pipeline. **SILENTLY
WRONG** = produces a number/call that looks normal but rests on an invalid
baseline. **no-call** = the `baseline <= 0 → None` guard suppresses output
(`bin/cnv_traits_common.py:104-108`). **DEGRADED** = runs and emits a VCF, but
sensitivity is silently reduced; nothing flags the loss.

---

## Detailed findings

### 1. It ABORTS at the coverage gate (default behaviour)

`modules/mosdepth/coverage.nf` runs `mosdepth --no-per-base --quantize ...
--by 50000` and then:

```
MEAN_DEPTH=$(awk '/^total/{print $4; exit}' <sample>.mosdepth.summary.txt)
if awk -v d="$MEAN_DEPTH" -v m="${min_depth}" 'BEGIN{exit (d >= m) ? 0 : 1}'; then
    echo "PASS: mean depth ..."
else
    echo "ERROR: mean depth $MEAN_DEPTH < ${min_depth}x (required). Aborting." >&2
    exit 1
fi
```
(`coverage.nf:30-37`)

- **`min_depth` default = 25** (`nextflow.config:18`; `conf/test.config:5`
  lowers it to 5 for CI only).
- The `MEAN_DEPTH` is the **`total` row** of the mosdepth summary — the
  genome-wide mean across all 50 kb windows. On WES the overwhelming majority
  of 50 kb windows are off-target (~0×), so this genome-wide mean is a few ×
  at most, far below 25 → the `else` branch fires.
- The abort is **unconditional** given a failing depth: `exit 1` with no
  override flag or bypass.
- **No partial output escapes first.** MOSDEPTH is invoked inside
  `PREPROCESS` (`subworkflows/preprocess.nf:50`), which runs **before** any SV
  / CNV / trait / SMN caller in `workflows/svcaller.nf`. A failed gate kills the
  run before a single variant is produced.

So on a WES BAM with default parameters, SVcaller **produces nothing** — it
stops at QC.

### 2. If the gate is tuned/bypassed, the depth traits are silently wrong

All the copy-number trait interpreters normalize to the **`CTRL_*` control
windows** of `assets/cnv_trait_regions.bed`:

```
chr2  100000000 100020000 CTRL_1
chr3  100000000 100020000 CTRL_2
...  (chr4,5,7,8,11,12)   CTRL_8
```

These are **fixed intergenic coordinates at ~100 Mb** — pure off-target on any
exome capture, i.e. ~0 reads. `control_baseline()` takes the **median** of the
eight CTRL depths (`bin/cnv_traits_common.py:64-73`), which collapses to ≈ 0.
`depth_ratio()` / `estimate_copies()` then either:

- return **`None`** via the `baseline <= 0` guard
  (`cnv_traits_common.py:104-108`) → **no-call**, or
- (a few stray off-target reads make the baseline a tiny non-zero number) →
  **wildly inflated ratios / copy numbers** presented as real.

This invalidates every consumer of `control_baseline()`:
`bin/hba_depth.py:289` (α-globin/thalassemia), `bin/amy1_cn.py:21` (AMY1),
`bin/lpa_kiv2.py` (LPA KIV-2), `bin/rh_status.py` (RHD), `bin/gst_null.py`
(GSTM1/GSTT1). The α-globin caller additionally carries per-segment WGS
calibration baselines (`assets/hba_segments.bed` col 5, from
`validation/giab_alpha_baseline.tsv`) that are not valid for capture data.

The genome-wide CNV callers fail for the same reason at a different layer:
- **GATK gCNV** denoises against a **panel of normals** built from WGS BAMs over
  the same intervals (`workflows/pon_build.nf`); the WES on-/off-target
  coverage profile is nothing like the WGS PON, so denoised copy ratios are
  meaningless. `modules/gatk/gcnv_call.nf:28`.
- **CNVpytor** builds a **genome-wide read-depth histogram baseline**
  (`modules/cnvpytor/call.nf`); off-target near-zero depth destroys that
  baseline.
- **SMN1/SMN2** (`subworkflows/smn_calling.nf`, `modules/smn_caller`) uses
  Illumina's SMNCopyNumberCaller, whose normalization is trained/calibrated for
  WGS depth profiles → unreliable on WES.

### 3. SV callers miss intronic breakpoints; nothing flags it

Manta (`modules/manta/call.nf`), DELLY (`modules/delly/call.nf`), GRIDSS
(`modules/gridss/call.nf`) and Scramble (`modules/scramble/call.nf`) all run
**whole-genome** with **no target/capture BED** restricting or informing them.
They will run, but SV breakpoints overwhelmingly land in introns and
intergenic sequence that exome capture does not cover, so most events have no
supporting reads and are **missed**. The only BED produced in the SV path is
`modules/manta/residual_regions.nf` (Manta's own non-PASS regions for tiered
GRIDSS) — not a capture BED. **Nothing in the pipeline detects or reports the
WES sensitivity loss**; the VCF simply looks sparse.

### 4. `params.intervals` does NOT make the pipeline exome-aware

`params.intervals` defaults to `null` (`nextflow.config`). It is wired into
**exactly one place**: GATK read-counting.

- `main.nf:62-63` → `workflows/svcaller.nf:22,57` →
  `subworkflows/cnv_calling.nf:60,66` → `GATK_PREPROCESS_INTERVALS` and
  `GATK_GCNV_CALL -L ${intervals}` (`modules/gatk/gcnv_call.nf:28`).
- Also `workflows/pon_build.nf:36,55` when **building** the gCNV PON.

It is **not** consumed by:
- the **mosdepth coverage QC** gate (`--by 50000`, whole-genome — the abort in §1),
- **any SV caller** (Manta/DELLY/GRIDSS/Scramble),
- **CNVpytor**,
- the **depth traits** (which use the hard-coded `assets/cnv_trait_regions.bed`
  with its off-target CTRL windows, not `params.intervals`),
- **SMN**.

So even supplying a correct exome capture BED as `--intervals`:
- does **not** stop the QC gate from aborting,
- does **not** fix the trait baselines,
- does **not** restore SV sensitivity,
- and only feeds gCNV read counting — which is still meaningless because its PON
  is WGS-derived, and matching a WES capture BED to a WGS PON is invalid anyway.

---

## Bottom line

**SVcaller produces no reliable output on a WES BAM.** With defaults it aborts
at the coverage QC gate (`coverage.nf:30`, `min_depth = 25`) before any caller
runs. If that gate is lowered or removed, the depth-based results
(α-globin/thalassemia, AMY1, LPA, RHD, GST, gCNV, CNVpytor, SMN) become
silently wrong because their baselines (fixed intergenic CTRL windows and a
WGS panel-of-normals) have no valid signal on exome data, and the SV callers
silently miss the intronic breakpoints that make up most structural variation.
`params.intervals` only reaches GATK gCNV and cannot rescue the QC, SV, or
trait-depth paths.

## Recommendation

- **Do not run this pipeline on WES / hybrid-capture BAMs.** It is designed and
  validated for ~30× human WGS (see `docs/what-a-30x-bam-yields.md`).
- If WES support is ever required, it is a **feature project**, not a parameter
  change: a capture-BED-aware coverage QC (on-target mean depth, not
  genome-wide), capture-matched trait baselines / control windows, a
  WES-specific gCNV PON, and explicit intronic-breakpoint sensitivity caveats
  on the SV output.
- **Minimum near-term guard** (separate from this docs-only change): have the
  QC/entry point detect a likely-WES BAM (e.g. large fraction of near-zero
  50 kb windows / very low genome-wide-to-on-target depth ratio) and **refuse
  with a clear "WES not supported" error** rather than aborting with a generic
  low-depth message or, worse, running past a lowered `min_depth` into
  silently-wrong output.
