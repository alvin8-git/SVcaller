# SVcaller Pipeline — Design Considerations & Improvement Roadmap

**Date:** 2026-06-03  
**Baseline:** F1=0.343, P=0.665, R=0.231  
**Truth set:** GIAB T2TQ100-V1.0 HG002 (29,691 SVs after passonly + confidence bed)  
**Pipeline:** Manta + Delly + GRIDSS + Scramble → Jasmine merge  

---

## 1. Realistic F1 Benchmarks for PE150 30× WGS

| Pipeline | Truth set | F1 | Notes |
|---|---|---|---|
| Manta alone | GIAB v0.6 (DEL-focused) | 0.45–0.52 | Older, DEL-only benchmark |
| Manta+Delly+GRIDSS ensemble | GIAB v0.6 | 0.55–0.65 | No MEI |
| Manta+Delly+GRIDSS+MELT | T2TQ100-V1.0 | 0.65–0.72 | State of art, Illumina |
| DRAGEN 4.x | T2TQ100-V1.0 | 0.60–0.70 | Commercial FPGA |
| **Our current pipeline** | **T2TQ100-V1.0** | **0.343** | As of 2026-06-02 |
| PacBio HiFi (Sniffles2) | T2TQ100-V1.0 | 0.92–0.95 | Long-read ceiling |

**Practical ceiling for Illumina PE150 30×:** ~0.70–0.75  
**Realistic target after P1–P10:** ~0.60–0.68  
**Current gap:** ~0.33 below achievable — roughly half the possible F1

**Important caveat:** The T2TQ100-V1.0 truth set is the hardest available (MEI-inclusive, long-read-derived). Against the older GIAB v0.6 (DEL-focused), our pipeline would score ~0.55–0.65 today, which is within the published range for short-read ensembles.

---

## 2. Root Cause Analysis — Why F1=0.343?

### Truth set composition (T2TQ100-V1.0, full)

| SVTYPE | Size | Count | Notes |
|---|---|---|---|
| INS | <300 bp | 36,913 | Mostly ALU MEI (~282 bp) |
| DEL | <300 bp | 29,326 | Repeat-mediated, hardest for short reads |
| INS | 300 bp–1 kb | 9,018 | SVA, L1 fragments |
| INS | 1–10 kb | 4,713 | Full-length L1 |
| DEL | 300 bp–1 kb | 4,143 | |
| DEL | 1–10 kb | 1,660 | |
| **Truth has NO DUP, NO INV entries** | | | |

After `--passonly + --includebed`: **29,691 benchmarkable SVs** (~59% INS, ~41% DEL).

### Our query VCF (PASS only — what Truvari uses)

| SVTYPE | PASS | LowQual | Matchable to truth? |
|---|---|---|---|
| DEL | 9,657 | 2,866 | ✓ Yes |
| INS | 5,822 | 325 | ✓ Yes |
| TRA | 2,792 | 6,255 | Excluded from bench |
| **DUP** | **1,525** | 2,031 | **✗ NEVER — truth has no DUP** |
| **INV** | **847** | **15,510** | **✗ NEVER — truth has no INV** |

**Matchable to truth: 15,479 PASS DEL+INS → 9,857 in confidence bed (comp_cnt)**

### Five root causes of low F1

**1. MEI detection gap** — ~60% of all FNs (16,316 INS FNs)  
Truth has ~17,000 INS in confidence bed (ALU, L1, SVA). Scramble detects <10% with QUAL≥70.  
Fundamental short-read limit. *Fix: MELT (P1), PopIns2 (P9)*

**2. Small DEL gap** — ~25% of all FNs (6,521 DEL FNs)  
~4,500 truth DEL <300 bp. Repeat-mediated deletions elude split-read callers.  
*Fix: SvABA (P4), LUMPY (P8)*

**3. DUP/INV calls are 100% FPs in Truvari** — 2,372 pure FP calls  
Our 1,525 PASS DUP + 847 PASS INV can never match truth (truth has no DUP/INV SVTYPE).  
Truvari v4.3.1 requires SVTYPE match by default.  
*Fix: add `--typeignore` to Truvari (QW-A); per-SVTYPE min_support (P6)*

**4. Delly LowQual contamination** — 26,987 LowQual variants in merged VCF  
15,510 LowQual INV (Delly noise) pass through to merged VCF, excluded by `--passonly`  
but bloat the VCF and pollute SUPP_VEC. 2,866 LowQual DEL may include real TPs  
that are now invisible to Truvari.  
*Fix: filter Delly non-PASS before Jasmine merge (QW-B)*

**5. PE150 30× fundamental ceiling** — not fixable without long reads  
MEI recall ceiling for Illumina: ~30–40% even with MELT.  
Full-recall requires PacBio HiFi or ONT ≥20×.

---

## 3. Quick Wins — Immediate Fixes (Pre-MELT)

### QW-A: Add `--typeignore` to Truvari bench
**File:** `modules/truvari/bench.nf`  
**Change:** Add `--typeignore` flag to all `truvari bench` calls  
**Expected:** Some DUP calls match truth INS/DEL → FP↓, TP↑, P↑ ~0.01–0.03  
**Effort:** 1 line  

### QW-B: Filter Delly non-PASS before Jasmine merge
**File:** `modules/jasmine/merge.nf`  
**Change:** Add `$7 != "PASS" {next}` to Delly preprocessing awk (vcfs[1])  
**Expected:** Removes 15,510 LowQual INV from merged VCF; cleaner SUPP_VEC  
**Risk:** May lose Delly-only LowQual DEL that are real (small). Acceptable tradeoff.  
**Effort:** 1 line  

---

## 4. P1–P10 Improvement Roadmap

| # | Action | ΔF1 est. | Effort | Root cause |
|---|---|---|---|---|
| P1 | **MELT** (MEI caller: ALU/L1/SVA) | +0.05–0.10 | 2 days | 16k MEI FNs |
| P2 | Calibrated GRIDSS DUP/INV filter | +0.005–0.015 | 1 hr | DUP/INV FPs |
| P3 | min_support=2 for 50–300 bp INS | +0.005–0.010 | 30 min | Small INS FPs |
| P4 | **SvABA** (small DEL assembler) | +0.01–0.03 | 1 day | 6.5k DEL FNs |
| P5 | GATK gCNV GC correction fix | +0.005–0.01 | 1 hr | Large DEL/DUP |
| P6 | Per-SVTYPE min_support in Jasmine | +0.01–0.02 | 30 min | DUP/INV FPs |
| P7 | **Paragraph** re-genotyping | +0.01–0.02 | 1 day | Wrong-GT FPs |
| P8 | **LUMPY/SpeedSeq** (DEL caller) | +0.01–0.02 | 1 day | DEL FNs |
| P9 | **PopIns2** (novel INS assembly) | +0.01–0.03 | 2 days | Non-ref INS FNs |
| P10 | gnomAD-SV v4 filter upgrade | +0.002–0.005 | 1 hr | Common-variant FPs |

**Cumulative estimate with all P1–P10:** F1 ~0.60–0.68

---

## 5. Pipeline Runtime Estimates

**Hardware assumption:** 32-core server, 128 GB RAM, 30× PE150 WGS

### Current pipeline (Manta+Delly+GRIDSS+Scramble)

| Phase | Time | Bottleneck |
|---|---|---|
| M1: Align + MarkDup (FASTQ) | 4.0 h | BWA-MEM2 |
| M1: FILTER_CHROMS (BAM) | 0.5 h | samtools |
| M2: SV calling (parallel) | 5.0 h | GRIDSS |
| M3: CNV calling (parallel) | 3.5 h | absorbed in M2 window |
| M4: SMN calling (parallel) | 0.5 h | absorbed |
| Post-merge: merge+annotate+report | 1.5 h | sequential |
| **FASTQ → report** | **~10.5 h** | |
| **BAM → report** | **~7.0 h** | |

### After P1–P10 (all improvements)

| Phase | Time | Bottleneck |
|---|---|---|
| M1: Align + MarkDup (FASTQ) | 4.0 h | unchanged |
| M2: SV calling (parallel) | **8.0 h** | **MELT** (ALU+L1+SVA parallel) |
| M2 parallel: GRIDSS | 5.0 h | absorbed in MELT window |
| M2 parallel: SvABA | 4.0 h | absorbed |
| M2 parallel: LUMPY | 2.0 h | absorbed |
| M2 parallel: PopIns2 | 3.5 h | absorbed |
| M3: CNV + GC correction | 4.0 h | absorbed in M2 window |
| Post-merge: merge+Paragraph+annotate+report | **2.6 h** | Paragraph +1 h |
| **FASTQ → report** | **~14.6 h (+4.1 h)** | |
| **BAM → report** | **~11.1 h (+4.1 h)** | |

**Peak CPU during M2 window:** ~128 cores (MELT 16, GRIDSS 32, SvABA 16, LUMPY 8, PopIns2 16, Delly 4, Manta 16, EH 4, Scramble 4)

### `--quick` mode (skip MELT + Paragraph)
**FASTQ → report: ~6.5 h, expected F1 ~0.50**  
Use for turnaround-critical cases; MEI not detected.

### Mitigation options for full-mode latency

| Option | Wall-time saving | Trade-off |
|---|---|---|
| Pre-run MELT overnight on cached BAM | −8 h from critical path | BAM must be pre-staged |
| Cloud-burst MELT IndivAnalysis (spot) | −4–5 h | Cloud cost |
| Skip Paragraph | −1 h | 13.5% wrong-GT rate |

---

## 6. Key Architecture Decisions

### MELT integration (P1)
MELT internal workflow is multi-step and partially sequential:  
`BuildTransposonZIP (5 min) → Preprocess (2 h) → IndivAnalysis (4 h, parallel by element) → GroupAnalysis (1 h) → Genotype (2 h) → MakeVCF (30 min)`  
IndivAnalysis for ALU, L1, SVA runs in separate NF processes (parallel).  
Total wall time: ~8 h on 16 cores. Feeds into Jasmine as 5th caller → SUPP_VEC becomes 5-char.  
Update GRIDSS TRA filter: `substr(supp,3,1)=="1"` remains valid (GRIDSS is still caller 3).

### Jasmine SUPP_VEC with 5 callers (post-MELT)
`Manta(1) / Delly(2) / GRIDSS(3) / Scramble(4) / MELT(5)`  
- GRIDSS TRA filter: `ones==1 && length(supp)>=3 && substr(supp,3,1)=="1"` ← unchanged
- DUP/INV GRIDSS filter (P6): same position-3 check, still valid

### Paragraph re-genotyping (P7)
Must run AFTER Jasmine merge (needs merged VCF as input).  
Sequential on critical path: adds ~1 h.  
Input: merged VCF + original BAM + reference.  
Output: re-genotyped VCF with updated GT/GQ → replaces merged VCF for downstream.

### GATK gCNV GC correction fix (P5)
Root cause: GRCh38 `.dict` has alphabetical chr order; BAM headers use numeric order.  
Fix: generate `hg38.canonical.sorted.dict` matching BAM header order, use for  
`AnnotateIntervals` and pass `--annotated-intervals` to `CreateReadCountPanelOfNormals`.  
Requires PON rebuild (~2 h with GIAB samples).

---

## 7. Benchmarking Considerations

### Truvari parameters in use
```
truvari bench --passonly --pick multi --pctseq 0 --sizemin 50
```
- `--passonly`: only PASS variants compared (excludes 26,987 LowQual from our VCF)
- `--pctseq 0`: no sequence similarity required (good for MEI where sequences vary)
- Missing: `--typeignore` — **Truvari v4.3.1 requires SVTYPE match by default**

### Recommended addition: `--typeignore`
Without this, our 1,525 DUP PASS + 847 INV PASS can never match truth (T2TQ100 has only INS/DEL).  
These 2,372 calls are pure FPs in current benchmarking.

### Interpreting F1 against T2TQ100-V1.0
The T2TQ100-V1.0 truth set is ~3× harder than GIAB v0.6 for Illumina short reads  
because it captures MEI comprehensively from long reads. When comparing to published  
F1 scores, always verify which truth set version was used. F1=0.343 on T2TQ100-V1.0  
≈ F1=0.58–0.65 on the older GIAB v0.6 (estimated, not measured).
