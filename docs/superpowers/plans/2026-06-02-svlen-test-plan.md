# SVLEN Fix — Test Plan

**Plan:** 2026-06-02-svcaller-improvements.md Phase 1
**Target file:** `modules/jasmine/merge.nf` lines 59-74

## Pre-condition: Diagnostic Gate (mandatory — no code change before this)

```bash
# 1. Find most recent Scramble work dir
ls -lt /data/alvin/SVcaller/work/*/HG002.vcf 2>/dev/null | head -5

# 2. Check for MEINFO field in Scramble VCF output
bcftools view /data/alvin/SVcaller/work/<scramble-hash>/HG002.vcf \
  | grep "^[^#]" | cut -f8 | grep -c "MEINFO="
# If count > 0: MEINFO present → implement fix
# If count = 0: MEINFO absent → use fallback (remove override entirely)

# 3. Sample MEINFO values for L1 entries
bcftools view /data/alvin/SVcaller/work/<scramble-hash>/HG002.vcf \
  | grep "MEINFO=" | grep "L1" | cut -f8 | head -10
# Expected format: MEINFO=L1HS,start,end,strand
# Observed SVLEN = end - start (field 3 minus field 2 of MEINFO value)
```

## Truvari Baseline Snapshot

Before any code change, record current INS metrics:

```bash
python3 -c "
import json
with open('/data/alvin/SVcaller/results/HG002/truvari/summary.json') as f:
    d = json.load(f)
print('TP:', d.get('TP-base', '?'), 'FP:', d.get('FP', '?'),
      'FN:', d.get('FN', '?'), 'F1:', d.get('f1', '?'))
"
```

## Implementation Verification (after code change)

```bash
# Unit: confirm awk L1 branch reads MEINFO
grep -n "MEINFO\|L1\|svlen" modules/jasmine/merge.nf | head -20

# Regression: Jasmine tier thresholds unchanged
grep -n "large_sv\|tier\|1000\|svlen" modules/jasmine/merge.nf | head -20
```

## Integration Test

```bash
# Re-run pipeline with -resume (only JASMINE_MERGE and downstream re-run)
nextflow run main.nf -profile docker \
  --input validation/validation_samplesheet.csv \
  --ref_fasta /data/alvin/ref/GRCh38/hg38.canonical.fa \
  --intervals /data/alvin/ref/GRCh38/wgs_autosomal.bed \
  --pon /data/alvin/SVcaller/pon/pon/giab_cnv_pon.hdf5 \
  --giab_truth /data/alvin/ref/GIAB/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz \
  --eh_catalog assets/eh_catalog.json \
  --annotsv_db /data/alvin/ref/annotsv/Annotations_Human \
  --outdir /data/alvin/SVcaller/results \
  -work-dir /data/alvin/SVcaller/work \
  -resume > /data/alvin/tmp/svlen_fix_run.log 2>&1 &

# Compare Truvari before/after
python3 -c "
import json
after = json.load(open('/data/alvin/SVcaller/results/HG002/truvari/summary.json'))
print('After — TP:', after.get('TP-base'), 'FP:', after.get('FP'),
      'F1:', after.get('f1'))
"
```

## Pass Criteria

| Check | Pass condition |
|-------|---------------|
| MEINFO diagnostic | Field present in ≥1 L1 row OR fallback path chosen |
| Truvari F1 | ≥ 0.370 (delta ≥ +0.029 from 0.341) |
| Truvari Precision | ≥ 0.680 |
| Truvari Recall | ≥ 0.250 |
| Jasmine tier regression | `grep large_sv` shows same threshold as before fix |
| Report renders | HTML opens, Circos visible, no N/A in coverage |
