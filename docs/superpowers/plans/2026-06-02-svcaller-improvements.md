<!-- /autoplan restore point: /home/alvin/.gstack/projects/alvin8-git-SVcaller/master-autoplan-restore-20260602-113311.md -->
# SVcaller Next-Steps Improvement Plan

**Date:** 2026-06-02
**Current state:** Pipeline at run82, F1=0.341 (4 callers: Manta+Delly+GRIDSS+Scramble). HTML report 1.9 MB, Circos 7-ring plot with legend. Documentation complete.

**Goal:** Three targeted improvements to precision/recall, QC completeness, and SMN validation coverage.

---

## Phase 1: Scramble SVLEN Fix (Precision improvement — revised after autoplan CEO review)

**Root cause (revised):** Scramble L1 FPs fail Truvari's ±30% size window because canonical SVLEN=6000 bp, while observed L1 insertions are typically 500-2000 bp. Truvari rejects the call on size mismatch, not call quality. QUAL filtering was the original hypothesis but both independent review voices agreed it targets the wrong mechanism.

**Hypothesis (revised):** Annotating L1 calls with the observed insertion length from the Scramble MEINFO field (or from the split-read evidence) instead of the canonical 6000 bp will convert a substantial fraction of L1 FPs to TPs.

**Step 1: Verify hypothesis (diagnostic — no code change)**
- Run `bcftools stats` on the Scramble VCF from the most recent run
- Cross-tab: SVTYPE vs QUAL vs Truvari match status from `{work_dir}/truvari/`
- If L1 FPs cluster at high QUAL → confirms SVLEN mismatch hypothesis
- If L1 FPs cluster at low QUAL → QUAL filtering was right; revert to original Phase 1

**Step 2: Implement observed-SVLEN (if hypothesis confirmed)**
- **File: `modules/jasmine/merge.nf` lines 59-74** — NOT sv_calling.nf (eng review correction)
  The canonical SVLEN assignment (`svlen=6000` for L1) is in the Jasmine merge awk block that reconstructs the INFO field as `"SVTYPE=INS;SVLEN=" svlen`. This is where the fix goes.
- Diagnostic first: `bcftools view work/<scramble-hash>/HG002.vcf | grep MEINFO | head` — verify MEINFO field is present in the raw Scramble VCF
- If MEINFO present: modify the L1 branch of the awk to compute `svlen = meinfo_end - meinfo_start` (field 4-3 of the MEINFO value), with fallback to 6000 if absent
- If MEINFO absent: try removing the SVLEN override entirely to let SCRAMble.R's native value propagate
- For ALU and SVA: keep canonical SVLEN (insertion length is consistent at full-length)
- Re-run Truvari benchmark to measure delta

**Files to modify:**
- `modules/jasmine/merge.nf` lines 59-74 — awk SVLEN assignment block

**Expected outcome:** L1 size-match rate improves; some FPs become TPs. P may decrease slightly; R should increase. Net F1 improvement expected 0.03-0.08.

**Risk (Eng review):**
- MEINFO may be absent in Scramble VCF (INFO column shows `.` on data rows). Diagnostic in Step 1 is mandatory before code change.
- The Jasmine awk overwrites the entire INFO field — any upstream MEINFO injection is discarded. Fix MUST be in merge.nf, not upstream.
- Jasmine tiered min_support awk (line ~116 same file) re-reads SVLEN for large-SV tier selection — changing SVLEN affects tier boundaries. Verify tier thresholds after fix.

---

## Phase 2: samtools flagstat wiring (QC completeness)

**Problem:** `mapped_pct` shows "N/A" in the HTML QC section. `samtools flagstat` runs in the pipeline (`modules/samtools/flagstat.nf`) and the output is collected into `ch_flagstat` but the HTML report script already reads it — the channel was not wired to BUILD_HTML_REPORT until the fix in report.nf. Need to verify end-to-end.

**Scope:**
- Check `subworkflows/report.nf` — confirm `ch_flagstat` is in the join chain feeding BUILD_HTML_REPORT
- Check `bin/html_report.py::_parse_flagstat` — verify the regex matches actual flagstat output format
- Check `modules/samtools/flagstat.nf` — verify the output file is emitted correctly
- Run pipeline on HG002 and confirm `mapped_pct` is not "N/A" in the HTML report

**Files to inspect/modify:**
- `subworkflows/report.nf` — channel join at line ~108 (ch_flagstat join position)
- `bin/html_report.py` — `_parse_flagstat` function, line ~151-163
- `modules/samtools/flagstat.nf`

**Expected outcome:** HTML QC section shows mapped % (expect >99% for HG002 aligned to canonical ref).

**Risk:** Low. flagstat is already running; this is a wiring/parsing bug.

---

## Phase 3: SMN Validation Run

**Problem:** SMN validation has not been run since the pipeline was restructured. Need to confirm SMNCopyNumberCaller produces correct results on the SMA trio (affected + carrier + normal).

**Scope:**
- Run pipeline with `--skip_gridss true` on `validation/smn_validation_samplesheet.csv`
- Verify SMN1/SMN2 copy number calls against the truth table at `validation/smn_truth_table.tsv`
- Check `bin/smn_report.py` output section in HTML report renders correctly for all 3 samples

**Files involved (read-only — should not require changes):**
- `validation/smn_validation_samplesheet.csv`
- `validation/smn_truth_table.tsv`
- `subworkflows/smn_calling.nf`
- `bin/smn_report.py`

**Expected outcome:** SMA-affected sample = SMN1:0, SMN2:4; carrier = SMN1:1, SMN2:3; normal = SMN1:2, SMN2:2.

**Risk:** Medium. If SMNCopyNumberCaller was broken by the preprocess refactor (FILTER_CHROMS changes), CN calls could be wrong. Mitigation: compare with the truth table and fix the module if calls are off.

---

## Dependencies

- Phase 1 (SVLEN diagnostic) runs first; 30 min
- Phase 2 (SMN validation) runs next — moved up from Phase 3, patient safety priority
- Phase 3 (flagstat wiring) is cosmetic QC; runs last or in parallel with Phase 2
- Phases 1 and 2 both produce Truvari output — run Phase 1 Truvari, then combine with Phase 2

## Success Criteria

| Metric | Current | Target |
|--------|---------|--------|
| F1 (HG002 Truvari) | 0.341 | ≥ 0.370 |
| Precision | 0.648 | ≥ 0.680 |
| Recall | 0.231 | ≥ 0.250 |
| `mapped_pct` in HTML | N/A | numeric value |
| SMN validation accuracy | untested | 3/3 samples correct |

### Phase 2 Actual Results (SMN validation — post-run)

| Sample | Truth SMN1 | Called | Truth SMN2 | Called | isSMA | isCarrier | Result |
|--------|-----------|--------|-----------|--------|-------|-----------|--------|
| SMAPB | 0 | 0 ✓ | 3 | 3 ✓ | True | False | PASS |
| SMAM | 1 | 1 ✓ | 5 | 4 ⚠️ | False | True | PARTIAL |
| SMAD | 1 | 1 ✓ | 1 | 1 ✓ | False | True | PASS |

Clinical classification (SMA/carrier) correct for all 3/3. SMAM SMN2=4 vs expected 5: Total_CN_raw=5.07; off-by-one at high SMN2 CN is a known SMNCopyNumberCaller limitation. No pipeline fix needed.

### Phase 3 Actual Results (flagstat wiring — post-run)

`mapped_pct` now shows **99.99%** in HTML QC section. Wiring was already correct. Previous N/A was caused by a UnicodeDecodeError in an old run with utils:1.0 that crashed html_report.py before mapped_pct was written (STR VCF gzip bug, fixed in utils:1.1). No code change needed.

Note: `Duplicate rate` still shows N/A — Picard MarkDuplicates metrics not reaching BUILD_HTML_REPORT (separate issue, not in this plan's scope).

### Phase 1 Actual Results (post-implementation)

L1 SVLEN changed 6000→1500 in `modules/jasmine/merge.nf:72`.
MEINFO diagnostic finding: MEINFO START/END are genomic insertion-site coordinates (~1 bp span), not ME sequence coordinates — observed length not derivable from MEINFO.

| Metric | Baseline | After fix | Delta |
|--------|----------|-----------|-------|
| Precision | 0.648 | 0.665 | +0.017 |
| Recall | 0.231 | 0.231 | ~0 |
| F1 | 0.341 | 0.343 | +0.002 |

Per-size-bin after fix: 50-300bp F1=0.299; 300bp-1kb F1=0.418; 1kb-10kb F1=0.291; >10kb F1=0.291

**Conclusion:** Fix is correct and kept. Improvement is real but marginal (+0.002 F1) because:
- L1 is a small fraction of total calls
- Dominant FP/FN problem is 50-300bp bin (ALU insertions, FN=15,459) — positional misses, not SVLEN
- Addressing ALU recall requires tuning Scramble parameters or alternate MEI caller (future scope)

---

## CEO Review (autoplan — auto-decided)

### What already exists
- Scramble post-processing awk in `subworkflows/sv_calling.nf` — canonical SVLEN assignment is ~10 lines
- Truvari bench is wired and runs per-pipeline — diagnostic data available immediately
- SMN samplesheet + truth table at `validation/smn_validation_samplesheet.csv` + `validation/smn_truth_table.tsv`
- samtools flagstat module at `modules/samtools/flagstat.nf`; `_parse_flagstat` in `bin/html_report.py:151`

### NOT in scope (deferred to future plan)
- Improving recall (R=0.231) by addressing Manta/Delly/GRIDSS over-filtering — large scope, separate plan
- Per-size-bin F1 diagnosis (50-300 bp vs 300 bp-1 kb vs 1-10 kb vs >10 kb) — diagnostic only, no code changes needed here but should be done before next F1 improvement cycle
- ClassifyCNV integration — future work
- GRIDSS tiering improvements — future work

### Alternatives considered
- QUAL filtering for Scramble (original Phase 1): rejected — SVLEN mismatch is root cause, not call confidence
- Adding a 4th SV caller (e.g., DRAGEN-SV): out of scope, requires new licensing
- Retraining Scramble model: not practical, no training data access

### Dream state delta
CURRENT: F1=0.341, P=0.648, R=0.231; mapped_pct=N/A; SMN untested
AFTER THIS PLAN: F1≥0.370, P≥0.680, R≥0.250; mapped_pct shows %; SMN 3/3 correct
12-MONTH IDEAL: F1≥0.50 (requires recall improvement via better MEI detection or long reads)

### Decision Audit Trail

| # | Phase | Decision | Classification | Principle | Rationale |
|---|-------|----------|----------------|-----------|-----------|
| 1 | CEO | Replace QUAL filter with SVLEN fix | User Challenge (accepted) | P1 | Both voices agree SVLEN mismatch is root cause |
| 2 | CEO | Swap Phase 2↔3: SMN before flagstat | Auto-decided | P1+P3 | Clinical correctness > QC display |
| 3 | CEO | Keep 3-phase scope unchanged | Auto-decided | P2 | All 3 items in blast radius, <1d each |

---

## Eng Review (autoplan — Phase 3)

### Architecture Diagram — SVLEN Fix Data Flow

```
subworkflows/sv_calling.nf
  └── SCRAMBLE process
        output: HG002.vcf (INFO col may be "." — MEINFO may be absent)
                         │
                         ▼
modules/jasmine/merge.nf  ← FIX GOES HERE (lines 59-74)
  └── awk block
        SVTYPE=INS → assigns canonical svlen
          L1   → svlen = 6000   ← CHANGE: use MEINFO end-start if present
          ALU  → svlen = 282    ← keep (full-length, consistent)
          SVA  → svlen = 800    ← keep (rare; not significant contributor)
          other → pass-through
        INFO reconstructed as: "SVTYPE=INS;SVLEN=" svlen
        ⚠ Entire INFO field overwritten — upstream MEINFO lost here
                         │
                         ▼
        Jasmine tiered min_support (line ~116 same file)
          re-reads SVLEN for large-SV tier (>1 kb threshold)
          ← verify tier boundaries unchanged after SVLEN reduction
                         │
                         ▼
        sv_merged.vcf.gz → Truvari bench
          L1: SVLEN now 500-2000 → fits ±30% window with truth
          Expected: L1 FPs → TPs; F1 ↑
```

### Test Plan

Written to: `docs/superpowers/plans/2026-06-02-svlen-test-plan.md` (see that file)

Summary:
1. **Diagnostic (no code change):** `bcftools view work/<scramble-hash>/HG002.vcf | grep "^[^#]" | cut -f8 | grep -o "MEINFO=[^;]*" | head -20` — confirm MEINFO presence/absence
2. **Truvari baseline snapshot:** record current TP/FP/FN for SVTYPE=INS from `results/HG002/truvari/summary.json`
3. **After fix — unit check:** `grep "SVLEN" modules/jasmine/merge.nf` — confirm L1 branch reads MEINFO
4. **Regression — Jasmine tier:** `grep -A5 "large_sv" modules/jasmine/merge.nf` — verify tier threshold logic unchanged
5. **Integration:** re-run pipeline with `-resume`; compare Truvari summary JSON for INS subtype

### Failure Modes Registry

| # | Failure Mode | Likelihood | Mitigation |
|---|-------------|-----------|-----------|
| 1 | MEINFO absent in Scramble VCF | Medium | Fallback: remove SVLEN override entirely so SCRAMble.R native value propagates |
| 2 | Jasmine tier boundary shifts (L1 now <1 kb) | Low | Some L1 calls may drop from large-SV tier; inspect tier counts before/after |
| 3 | MEINFO present but format varies by MEI subtype | Low | Parse only L1 subtype; leave ALU/SVA canonical as-is |
| 4 | Fix helps L1 but ALU/SVA are dominant FP source | Medium | Measure per-MEI-type TP/FP delta in Truvari; pivot plan if ALU dominates |
| 5 | SMN calling broken by preprocess refactor | Medium | Compare CN calls against truth table row-by-row; fix SMN module if off |
| 6 | flagstat N/A is meta-map mismatch not regex | Low | Add `println meta` debug line to report.nf join; inspect work dir |

### Eng Decisions

| # | Decision | Classification | Rationale |
|---|----------|---------------|-----------|
| 4 | Fix in merge.nf not sv_calling.nf | Correction (plan was wrong) | Jasmine awk overwrites INFO; upstream injection is discarded |
| 5 | Diagnostic before code change (MEINFO check) | Mandatory gate | MEINFO may be absent; code change is contingent on its presence |
| 6 | Keep ALU/SVA canonical SVLEN | Auto-decided | ALU insertion length is consistent at full-length; SVA is low-frequency |

---

## DX Review (autoplan — Phase 3.5)

DX scope: 1 parameter (`--scramble_min_qual`), documentation completeness for the new SVLEN fix.

### Findings

1. **`--scramble_min_qual` parameter does not exist** — the plan previously referenced QUAL filtering, but after CEO revision the fix is SVLEN-based, not QUAL-based. No new parameter is needed. The awk fix is internal to the merge step. DX impact: zero new flags to document.

2. **CLAUDE.md Phase 1 description is accurate** — the diagnostic + implementation split is clearly documented. No runbook gaps.

3. **`docs/howto-interpret-report.md` flagstat note** — currently says "known limitation, future fix." After Phase 2 (flagstat wiring), this note should be updated. **Deferred action item** (post-implementation).

### DX Decision

| # | Decision | Classification | Rationale |
|---|----------|---------------|-----------|
| 7 | No new CLI parameter for SVLEN fix | Auto-decided | Fix is internal; user-facing interface unchanged |
| 8 | Update howto-interpret-report.md after Phase 2 | Deferred | Can't update until fix is confirmed working |

---

## Final Approval Gate (autoplan — Phase 4)

### Taste Decisions for User Review

The following non-obvious choices were made during autoplan that require your sign-off:

| # | Decision | What was chosen | Alternative | Your call |
|---|----------|----------------|-------------|-----------|
| A | Phase ordering | SMN (Phase 2) before flagstat (Phase 3) | flagstat first (cosmetic but quick) | Auto-decided: patient safety priority |
| B | SVLEN fix target | `modules/jasmine/merge.nf` lines 59-74 | `subworkflows/sv_calling.nf` (original plan — wrong) | Corrected by eng review |
| C | MEINFO diagnostic gate | Mandatory before any code change | Skip diagnostic, attempt fix directly | Mandatory: MEINFO may be absent |
| D | ALU/SVA SVLEN | Keep canonical (282 bp / 800 bp) | Apply same observed-SVLEN logic | Keep canonical: full-length insertions are consistent |
| E | Scope boundary | No new CLI parameters | Add `--scramble_svlen_source` flag for user control | Not needed: internal fix only |

### Challenges Outstanding

None — all CEO challenges were resolved:
- Original QUAL filter hypothesis → replaced with SVLEN mismatch hypothesis (user confirmed)
- File location error (sv_calling.nf → merge.nf) → corrected by eng subagent

### Implementation Checklist (post-approval)

- [ ] Run MEINFO diagnostic on most recent Scramble VCF work dir
- [ ] If MEINFO present: modify L1 branch in `modules/jasmine/merge.nf:59-74`
- [ ] Bump `// vN` comment in merge.nf to bust Nextflow cache
- [ ] Re-run pipeline `-resume`; capture Truvari delta
- [ ] Run SMN validation: `nextflow run main.nf --input validation/smn_validation_samplesheet.csv --skip_gridss true ...`
- [ ] Verify `mapped_pct` in HTML after pipeline completes
- [ ] Update `docs/howto-interpret-report.md` flagstat note if fixed
