# Per-Size-Bin F1 Diagnosis

**Date:** 2026-06-02  
**Run:** post-SVLEN-fix (L1 1500 bp, commit e2ce614)  
**Overall:** F1=0.343, P=0.665, R=0.231

## Results by Size Bin

| Bin | P | R | F1 | TP | FP | FN |
|-----|---|---|----|----|----|-----|
| 50-300bp | 0.674 | 0.192 | 0.299 | 3,682 | 1,635 | 15,459 |
| 300bp-1kb | 0.551 | 0.337 | 0.418 | 2,479 | 1,892 | 4,888 |
| 1kb-10kb | 0.452 | 0.215 | 0.291 | 652 | 784 | 2,382 |
| >10kb | 0.618 | 0.190 | 0.291 | 34 | 21 | 145 |

## FN/FP SVTYPE Breakdown

```
50-300bp  FN: 10,012 INS  + 5,447 DEL
          FP:  1,116 INS  +   318 DUP  +  163 DEL  +  38 INV

300bp-1kb FN:  4,023 INS  +   865 DEL
          FP:    627 INS  +   598 INV  +  455 DUP  + 212 DEL

1kb-10kb  FN:  2,194 INS  +   188 DEL
          FP:    388 INS  +   139 DUP  +  138 DEL  + 119 INV

>10kb     FN:    110 INS  +    35 DEL
          FP:     13 DUP  +     8 DEL

Overall   FN: 16,316 INS  + 6,521 DEL
          FP:  1,027 INS  +   925 DUP  + 751 INV   + 509 DEL
```

## Root Cause Analysis

### 1. INS FN = 16,316 (71% of all FNs) — MEI detection gap
- ALU (~282 bp, 50-300bp bin): 10,012 FNs. Scramble misses ~73% of GIAB ALU truth insertions.
- L1 (1500 bp, 1kb-10kb bin): 2,194 FNs. Scramble misses most L1 truth insertions.
- SVA/other (300bp-1kb): 4,023 FNs.
- **Root cause**: Scramble is limited by Illumina short-read sensitivity for MEI. QUAL≥70 filter may be over-filtering (currently 70; lowering to 50 might recover ALUs).

### 2. DEL FN = 6,521 (29% of all FNs) — small deletion gap
- 50-300bp bin: 5,447 DEL FNs. Manta, Delly, and GRIDSS all miss small deletions.
- These are likely repeat-mediated deletions <300 bp where split-read evidence is ambiguous.
- **Root cause**: Multi-caller min_support=1 (Jasmine), so even single-caller calls pass. Manta is the best small DEL caller here — its recall at <300 bp is known to drop off.

### 3. DUP FP = 925 + INV FP = 751 — GRIDSS artifact calls
- 1,676 combined DUP/INV FPs, mostly from GRIDSS.
- Current filter removes single-caller SVs only if SVLEN >10 kb. Many GRIDSS DUP/INV calls are <10 kb and pass this gate.
- **Root cause**: GRIDSS DUP/INV calls without corroboration from Manta/Delly are unreliable.

## Improvement Roadmap (ranked by expected F1 delta)

| Priority | Action | Expected F1 delta | Effort |
|----------|--------|------------------|--------|
| P1 | Lower Scramble QUAL threshold 70→50 | +0.02-0.05 | Low (1-line awk change) |
| P2 | Filter GRIDSS-only DUP/INV (not just >10kb, all sizes) | +0.01-0.02 | Low (awk filter change) |
| P3 | Add MELT as 5th MEI caller | +0.05-0.10 | High (new module) |
| P4 | Add Lumpy/SURVIVOR for small DEL recall | +0.01-0.03 | High (new module) |

## Quick Wins (can be done now)

### Quick Win 1: Lower Scramble QUAL threshold to 50
File: `modules/jasmine/merge.nf` line 70 — `if(\$6+0 < 70) next` → `if(\$6+0 < 50) next`  
Risk: More ALU calls → precision may drop in 50-300bp bin.  
Measure: Re-run Truvari and check 50-300bp bin F1 delta.

### Quick Win 2: Extend GRIDSS-only DUP/INV filter to all sizes
File: `modules/jasmine/merge.nf` lines 118-122 — post-merge awk filter  
Currently: `if(ones==1 && svlen>10000) next` removes single-caller calls only if >10kb  
Change: Also remove single-caller DUP and INV regardless of size (Manta/Delly corroboration required)  
Risk: May remove a few legitimate single-caller DUPs. Low impact on recall (TP DUP count is small).  
Measure: Check TP-base delta for DUP before/after.
