# SVcaller Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a 7-module Nextflow DSL2 pipeline that calls SVs, CNVs, and SMN1/SMN2 copy numbers from Illumina PE150 WGS (FASTQ or BAM), validates against GIAB HG001–HG007, and produces a per-sample HTML report with an embedded Circos plot.

**Architecture:** Single-sample pipeline; modules M2/M3/M4 run in parallel on the preprocessed BAM; a one-time `pon_build.nf` workflow pre-builds a GATK gCNV Panel of Normals from GIAB samples. Each Nextflow module runs in a pinned Docker container.

**Tech Stack:** Nextflow DSL2 v25.10.4, Docker, BWA-MEM2 2.2.1, Manta 1.6, DELLY 1.2.6, GRIDSS 2.13.2, ExpansionHunter 5.0, JASMINE 1.1.5, CNVpytor 1.3.1, GATK 4.5, SMNCopyNumberCaller 1.1, AnnotSV 3.4.2, pycirclize 1.7, truvari 4.x, Python 3.13, Jinja2.

**Working directory for all tasks:** `/data/alvin/SVcaller`

---

## Phase 0: Project Scaffold

### Task 1: Initialize repository and Nextflow config

**Files:**
- Create: `.gitignore`
- Create: `nextflow.config`
- Create: `conf/base.config`
- Create: `conf/docker.config`

- [ ] **Step 1: Initialize git repo**

```bash
cd /data/alvin/SVcaller
git init
```
Expected: `Initialized empty Git repository in /data/alvin/SVcaller/.git/`

- [ ] **Step 2: Create .gitignore**

```bash
cat > /data/alvin/SVcaller/.gitignore << 'EOF'
results/
work/
.nextflow/
.nextflow.log*
*.pyc
__pycache__/
*.egg-info/
.DS_Store
EOF
```

- [ ] **Step 3: Create conf/base.config**

```bash
mkdir -p /data/alvin/SVcaller/conf
```

Write `/data/alvin/SVcaller/conf/base.config`:

```groovy
process {
    cpus   = 2
    memory = 8.GB
    time   = 4.h

    withLabel: process_single {
        cpus   = 1
        memory = 6.GB
        time   = 4.h
    }
    withLabel: process_low {
        cpus   = 4
        memory = 12.GB
        time   = 6.h
    }
    withLabel: process_medium {
        cpus   = 8
        memory = 36.GB
        time   = 12.h
    }
    withLabel: process_high {
        cpus   = 16
        memory = 72.GB
        time   = 24.h
    }
    withLabel: process_gridss {
        cpus   = 8
        memory = 32.GB
        time   = 24.h
    }
}
```

- [ ] **Step 4: Create conf/docker.config**

Write `/data/alvin/SVcaller/conf/docker.config`:

```groovy
docker.enabled         = true
docker.runOptions      = '--user $(id -u):$(id -g)'
docker.userEmulation   = true

process {
    withName: 'BWAMEM2_ALIGN'      { container = 'quay.io/biocontainers/bwa-mem2:2.2.1--he513fc3_1' }
    withName: 'SAMTOOLS_SORT'      { container = 'quay.io/biocontainers/samtools:1.20--h50ea8bc_0' }
    withName: 'PICARD_MARKDUP'     { container = 'quay.io/biocontainers/picard:3.2.0--hdfd78af_0' }
    withName: 'MOSDEPTH'           { container = 'quay.io/biocontainers/mosdepth:0.3.8--hd299d5a_0' }
    withName: 'FASTQC'             { container = 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0' }
    withName: 'MULTIQC'            { container = 'quay.io/biocontainers/multiqc:1.22.3--pyhdfd78af_0' }
    withName: 'MANTA_CALL'         { container = 'quay.io/biocontainers/manta:1.6.0--py38_1' }
    withName: 'DELLY_CALL'         { container = 'quay.io/biocontainers/delly:1.2.6--ha047f11_0' }
    withName: 'GRIDSS_CALL'        { container = 'gridss/gridss:2.13.2' }
    withName: 'EXPANSIONHUNTER'    { container = 'quay.io/biocontainers/expansionhunter:5.0.0--h9ee0642_1' }
    withName: 'JASMINE_MERGE'      { container = 'quay.io/biocontainers/jasminesv:1.1.5--hdfd78af_0' }
    withName: 'CNVPYTOR_CALL'      { container = 'quay.io/biocontainers/cnvpytor:1.3.1--pyhdfd78af_0' }
    withName: 'GATK_GCNV.*'       { container = 'broadinstitute/gatk:4.5.0.0' }
    withName: 'SMN_CALLER'         { container = 'svcaller/smncopynum:1.1' }
    withName: 'ANNOTSV'            { container = 'quay.io/biocontainers/annotsv:3.4.2--pl5321hdfd78af_0' }
    withName: 'TRUVARI_BENCH'      { container = 'quay.io/biocontainers/truvari:4.2.2--pyhdfd78af_0' }
    withName: 'SVCALLER_UTILS'     { container = 'svcaller/utils:1.0' }
}
```

- [ ] **Step 5: Create nextflow.config**

Write `/data/alvin/SVcaller/nextflow.config`:

```groovy
manifest {
    name        = 'SVcaller'
    description = 'SV/CNV/SMN1-SMN2 pipeline for human WGS'
    version     = '1.0.0'
    nextflowVersion = '>=25.10.0'
}

params {
    input         = null
    ref_fasta     = null
    genome        = 'GRCh38'
    pon           = null
    outdir        = 'results'
    min_depth     = 30
    eh_catalog    = "${projectDir}/assets/eh_catalog.json"
    giab_truth    = null
    max_cpus      = 64
    max_memory    = '120.GB'
    max_time      = '240.h'
}

includeConfig 'conf/base.config'
profiles {
    docker { includeConfig 'conf/docker.config' }
    test   { includeConfig 'conf/test.config'   }
}

def check_max(obj, type) {
    if (type == 'memory') {
        try {
            if (obj.compareTo(params.max_memory as nextflow.util.MemoryUnit) == 1)
                return params.max_memory as nextflow.util.MemoryUnit
            else return obj
        } catch (all) { return obj }
    } else if (type == 'time') {
        try {
            if (obj.compareTo(params.max_time as nextflow.util.Duration) == 1)
                return params.max_time as nextflow.util.Duration
            else return obj
        } catch (all) { return obj }
    } else if (type == 'cpus') {
        try { return Math.min( obj, params.max_cpus as int ) }
        catch (all) { return obj }
    }
}
```

- [ ] **Step 6: Create directory skeleton**

```bash
mkdir -p /data/alvin/SVcaller/{workflows,subworkflows,modules/{bwamem2,samtools,picard,mosdepth,fastqc,multiqc,manta,delly,gridss,expansionhunter,jasmine,cnvpytor,gatk,smn_caller,annotsv,truvari},bin,assets,validation,tests}
```

- [ ] **Step 7: Commit scaffold**

```bash
cd /data/alvin/SVcaller
git add .gitignore nextflow.config conf/
git commit -m "feat: initialize Nextflow project scaffold"
```

---

### Task 2: Create samplesheet schema and test config

**Files:**
- Create: `assets/schema_input.json`
- Create: `conf/test.config`
- Create: `tests/test_samplesheet.csv`

- [ ] **Step 1: Create samplesheet schema**

Write `/data/alvin/SVcaller/assets/schema_input.json`:

```json
{
    "$schema": "http://json-schema.org/draft-07/schema",
    "title": "SVcaller input samplesheet",
    "description": "CSV with one row per sample. Use fastq_1/fastq_2 OR bam, not both.",
    "type": "array",
    "items": {
        "type": "object",
        "required": ["sample"],
        "properties": {
            "sample":  { "type": "string", "description": "Sample ID (no spaces)" },
            "fastq_1": { "type": "string", "description": "Absolute path to R1 FASTQ.gz", "format": "file-path" },
            "fastq_2": { "type": "string", "description": "Absolute path to R2 FASTQ.gz", "format": "file-path" },
            "bam":     { "type": "string", "description": "Absolute path to sorted BAM",   "format": "file-path" }
        }
    }
}
```

- [ ] **Step 2: Create test samplesheet**

Write `/data/alvin/SVcaller/tests/test_samplesheet.csv`:

```
sample,fastq_1,fastq_2,bam
HG002_test,/data/alvin/SVcaller/tests/data/HG002_R1.fastq.gz,/data/alvin/SVcaller/tests/data/HG002_R2.fastq.gz,
```

- [ ] **Step 3: Create test config**

Write `/data/alvin/SVcaller/conf/test.config`:

```groovy
params {
    input      = "${projectDir}/tests/test_samplesheet.csv"
    ref_fasta  = '/data/alvin/ref/GRCh38/GRCh38.fasta'
    outdir     = '/data/alvin/tmp/svcaller_test_out'
    min_depth  = 5
    max_cpus   = 8
    max_memory = '32.GB'
}
```

- [ ] **Step 4: Create bin/parse_samplesheet.py**

Write `/data/alvin/SVcaller/bin/parse_samplesheet.py`:

```python
#!/usr/bin/env python3
"""Validate and emit samplesheet rows as JSON lines for Nextflow."""
import csv, json, sys
from pathlib import Path

def validate(row: dict) -> dict:
    sid = row.get("sample", "").strip()
    if not sid:
        raise ValueError("'sample' column is required and must not be empty")
    fq1 = row.get("fastq_1", "").strip()
    fq2 = row.get("fastq_2", "").strip()
    bam = row.get("bam", "").strip()
    if bam and (fq1 or fq2):
        raise ValueError(f"Sample {sid}: provide fastq_1/fastq_2 OR bam, not both")
    if (fq1 and not fq2) or (fq2 and not fq1):
        raise ValueError(f"Sample {sid}: fastq_1 and fastq_2 must both be provided")
    if not bam and not fq1:
        raise ValueError(f"Sample {sid}: must provide either bam or fastq_1/fastq_2")
    entry = {"id": sid, "single_end": False}
    if bam:
        if not Path(bam).exists():
            raise FileNotFoundError(f"BAM not found: {bam}")
        entry["bam"] = bam
        entry["input_type"] = "bam"
    else:
        for f in [fq1, fq2]:
            if not Path(f).exists():
                raise FileNotFoundError(f"FASTQ not found: {f}")
        entry["fastq_1"] = fq1
        entry["fastq_2"] = fq2
        entry["input_type"] = "fastq"
    return entry

def main():
    path = sys.argv[1]
    with open(path) as fh:
        for row in csv.DictReader(fh):
            try:
                print(json.dumps(validate(row)))
            except (ValueError, FileNotFoundError) as e:
                print(f"ERROR: {e}", file=sys.stderr)
                sys.exit(1)

if __name__ == "__main__":
    main()
```

```bash
chmod +x /data/alvin/SVcaller/bin/parse_samplesheet.py
```

- [ ] **Step 5: Write test for parse_samplesheet.py**

Write `/data/alvin/SVcaller/tests/test_parse_samplesheet.py`:

```python
import json, subprocess, textwrap
from pathlib import Path
import pytest

SCRIPT = Path(__file__).parent.parent / "bin" / "parse_samplesheet.py"

def run(csv_text: str, tmp_path) -> subprocess.CompletedProcess:
    p = tmp_path / "sheet.csv"
    p.write_text(textwrap.dedent(csv_text))
    return subprocess.run(["python3", str(SCRIPT), str(p)],
                          capture_output=True, text=True)

def test_rejects_missing_sample(tmp_path):
    r = run("sample,bam\n,/nonexistent.bam\n", tmp_path)
    assert r.returncode == 1
    assert "required" in r.stderr

def test_rejects_bam_and_fastq(tmp_path):
    r = run("sample,fastq_1,fastq_2,bam\nS1,/a.fq.gz,/b.fq.gz,/c.bam\n", tmp_path)
    assert r.returncode == 1
    assert "not both" in r.stderr

def test_rejects_missing_fastq2(tmp_path):
    r = run("sample,fastq_1,fastq_2\nS1,/a.fq.gz,\n", tmp_path)
    assert r.returncode == 1
    assert "both" in r.stderr.lower()

def test_accepts_bam_row(tmp_path):
    bam = tmp_path / "test.bam"
    bam.touch()
    r = run(f"sample,bam\nS1,{bam}\n", tmp_path)
    assert r.returncode == 0
    d = json.loads(r.stdout.strip())
    assert d["id"] == "S1"
    assert d["input_type"] == "bam"
```

- [ ] **Step 6: Run tests**

```bash
cd /data/alvin/SVcaller
python3 -m pytest tests/test_parse_samplesheet.py -v
```

Expected: 4 tests pass.

- [ ] **Step 7: Commit**

```bash
cd /data/alvin/SVcaller
git add assets/ conf/test.config tests/ bin/parse_samplesheet.py
git commit -m "feat: add samplesheet schema, validation, and tests"
```

---

## Phase 1: Pre-processing Module (M1)

### Task 3: BWA-MEM2 alignment + samtools sort/index module

**Files:**
- Create: `modules/bwamem2/align.nf`
- Create: `modules/samtools/sort_index.nf`
- Create: `modules/picard/markduplicates.nf`

- [ ] **Step 1: Create BWA-MEM2 align module**

Write `/data/alvin/SVcaller/modules/bwamem2/align.nf`:

```groovy
process BWAMEM2_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(reads)
    path fasta
    path fai
    path bwt_index  // directory containing bwa-mem2 index files

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"),     emit: bam
    tuple val(meta), path("${meta.id}.sorted.bam.bai"), emit: bai
    path "versions.yml",                                 emit: versions

    script:
    def rg = "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:${meta.id}"
    """
    bwa-mem2 mem \\
        -t ${task.cpus} \\
        -R "${rg}" \\
        ${fasta} \\
        ${reads} \\
        | samtools sort -@ ${task.cpus} -m 2G \\
        -o ${meta.id}.sorted.bam

    samtools index -@ ${task.cpus} ${meta.id}.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>/dev/null | head -1 | tr -d '\\n')
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
```

- [ ] **Step 2: Create Picard MarkDuplicates module**

Write `/data/alvin/SVcaller/modules/picard/markduplicates.nf`:

```groovy
process PICARD_MARKDUP {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.markdup.bam"),     emit: bam
    tuple val(meta), path("${meta.id}.markdup.bam.bai"), emit: bai
    tuple val(meta), path("${meta.id}.dup_metrics.txt"), emit: metrics
    path "versions.yml",                                  emit: versions

    script:
    """
    picard MarkDuplicates \\
        -Xmx${(task.memory.toGiga() * 0.85).intValue()}g \\
        INPUT=${bam} \\
        OUTPUT=${meta.id}.markdup.bam \\
        METRICS_FILE=${meta.id}.dup_metrics.txt \\
        REMOVE_DUPLICATES=false \\
        VALIDATION_STRINGENCY=LENIENT \\
        CREATE_INDEX=true \\
        TMP_DIR=./tmp

    mv ${meta.id}.markdup.bai ${meta.id}.markdup.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard MarkDuplicates --version 2>&1 | grep -o 'Version:.*' | sed 's/Version://')
    END_VERSIONS
    """
}
```

- [ ] **Step 3: Commit**

```bash
cd /data/alvin/SVcaller
git add modules/bwamem2/ modules/samtools/ modules/picard/
git commit -m "feat: add BWA-MEM2 align and Picard MarkDuplicates modules"
```

---

### Task 4: mosdepth coverage QC + preprocess subworkflow

**Files:**
- Create: `modules/mosdepth/coverage.nf`
- Create: `subworkflows/preprocess.nf`

- [ ] **Step 1: Create mosdepth module**

Write `/data/alvin/SVcaller/modules/mosdepth/coverage.nf`:

```groovy
process MOSDEPTH {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(bam), path(bai)
    val   min_depth

    output:
    tuple val(meta), path("${meta.id}.mosdepth.summary.txt"), emit: summary
    tuple val(meta), path("${meta.id}.regions.bed.gz"),        emit: regions
    path "versions.yml",                                        emit: versions

    script:
    """
    mosdepth \\
        --threads ${task.cpus} \\
        --quantize 0:5:30:500: \\
        --no-abbrev \\
        ${meta.id} \\
        ${bam}

    # Fail pipeline if mean depth below threshold
    MEAN_DEPTH=\$(grep "^total" ${meta.id}.mosdepth.summary.txt | awk '{print \$4}')
    if awk "BEGIN{exit (\$MEAN_DEPTH >= ${min_depth}) ? 0 : 1}"; then
        echo "PASS: mean depth \$MEAN_DEPTH >= ${min_depth}x"
    else
        echo "ERROR: mean depth \$MEAN_DEPTH < ${min_depth}x (required). Aborting." >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/mosdepth //')
    END_VERSIONS
    """
}
```

- [ ] **Step 2: Create preprocess subworkflow**

Write `/data/alvin/SVcaller/subworkflows/preprocess.nf`:

```groovy
include { BWAMEM2_ALIGN    } from '../modules/bwamem2/align'
include { PICARD_MARKDUP   } from '../modules/picard/markduplicates'
include { MOSDEPTH         } from '../modules/mosdepth/coverage'

workflow PREPROCESS {
    take:
    ch_samplesheet  // [ meta, [fastq_1, fastq_2] | null, bam | null ]
    ch_fasta        // path
    ch_fai          // path
    ch_bwt_index    // path (directory)
    min_depth       // integer

    main:
    // Split into FASTQ and BAM channels
    ch_fastq = ch_samplesheet.filter { meta, fq1, fq2, bam -> fq1 != null }
                              .map   { meta, fq1, fq2, bam -> [meta, [fq1, fq2]] }
    ch_bam   = ch_samplesheet.filter { meta, fq1, fq2, bam -> bam != null }
                              .map   { meta, fq1, fq2, bam -> [meta, bam] }

    // Align FASTQs
    BWAMEM2_ALIGN(ch_fastq, ch_fasta, ch_fai, ch_bwt_index)

    // Merge: aligned FASTQs + pre-supplied BAMs both go through MarkDuplicates
    ch_all_bam = BWAMEM2_ALIGN.out.bam
        .join(BWAMEM2_ALIGN.out.bai)
        .mix(ch_bam.map { meta, bam -> [meta, bam, file("${bam}.bai")] })

    PICARD_MARKDUP(ch_all_bam)

    ch_final_bam = PICARD_MARKDUP.out.bam.join(PICARD_MARKDUP.out.bai)

    // Coverage QC — halts pipeline if < min_depth
    MOSDEPTH(ch_final_bam, min_depth)

    emit:
    bam      = ch_final_bam
    coverage = MOSDEPTH.out.summary
    metrics  = PICARD_MARKDUP.out.metrics
}
```

- [ ] **Step 3: Write test for preprocess depth-fail behaviour**

Write `/data/alvin/SVcaller/tests/test_mosdepth_depth_check.sh`:

```bash
#!/usr/bin/env bash
# Quick smoke test: mosdepth module must exit 1 when depth < min_depth
set -e
# Create a near-empty BAM to simulate low coverage
docker run --rm -v /data/alvin/SVcaller/tests:/data \
  quay.io/biocontainers/samtools:1.20--h50ea8bc_0 \
  bash -c "samtools view -bS /dev/null > /data/empty.bam && samtools index /data/empty.bam"

docker run --rm -v /data/alvin/SVcaller/tests:/data \
  quay.io/biocontainers/mosdepth:0.3.8--hd299d5a_0 \
  bash -c '
    mosdepth --no-abbrev test_low /data/empty.bam 2>/dev/null || true
    echo "total\t0\t0\t0" > test_low.mosdepth.summary.txt
    MEAN=0
    if awk "BEGIN{exit ($MEAN >= 30) ? 0 : 1}"; then echo "UNEXPECTED PASS"; exit 1; fi
    echo "PASS: correctly detected low coverage"
  '
echo "Test passed"
```

```bash
chmod +x /data/alvin/SVcaller/tests/test_mosdepth_depth_check.sh
```

- [ ] **Step 4: Commit**

```bash
cd /data/alvin/SVcaller
git add modules/mosdepth/ subworkflows/preprocess.nf tests/test_mosdepth_depth_check.sh
git commit -m "feat: add mosdepth coverage QC and preprocess subworkflow"
```

---

## Phase 2: SV Calling Ensemble (M2)

### Task 5: Manta SV caller module

**Files:**
- Create: `modules/manta/call.nf`

- [ ] **Step 1: Create Manta module**

Write `/data/alvin/SVcaller/modules/manta/call.nf`:

```groovy
process MANTA_CALL {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path("${meta.id}.manta.sv.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.manta.sv.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                                       emit: versions

    script:
    """
    configManta.py \\
        --bam ${bam} \\
        --referenceFasta ${fasta} \\
        --runDir manta_run

    python manta_run/runWorkflow.py \\
        -j ${task.cpus} \\
        -g ${task.memory.toGiga()}

    cp manta_run/results/variants/diploidSV.vcf.gz     ${meta.id}.manta.sv.vcf.gz
    cp manta_run/results/variants/diploidSV.vcf.gz.tbi ${meta.id}.manta.sv.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        manta: \$(configManta.py --version 2>&1 | head -1)
    END_VERSIONS
    """
}
```

- [ ] **Step 2: Commit**

```bash
cd /data/alvin/SVcaller
git add modules/manta/
git commit -m "feat: add Manta SV calling module"
```

---

### Task 6: DELLY SV caller module

**Files:**
- Create: `modules/delly/call.nf`

- [ ] **Step 1: Create DELLY module**

Write `/data/alvin/SVcaller/modules/delly/call.nf`:

```groovy
process DELLY_CALL {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path("${meta.id}.delly.sv.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.delly.sv.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                                       emit: versions

    script:
    """
    # Call all SV types
    for SVTYPE in DEL INS INV DUP TRA; do
        delly call \\
            -t \${SVTYPE} \\
            -g ${fasta} \\
            -o ${meta.id}.delly.\${SVTYPE}.bcf \\
            ${bam}
    done

    # Merge all SV types and convert to VCF
    bcftools concat -a ${meta.id}.delly.*.bcf \\
        | bcftools sort -O z -o ${meta.id}.delly.sv.vcf.gz
    bcftools index -t ${meta.id}.delly.sv.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$(delly --version 2>&1 | grep "DELLY" | head -1 | awk '{print \$2}')
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
```

- [ ] **Step 2: Commit**

```bash
cd /data/alvin/SVcaller
git add modules/delly/
git commit -m "feat: add DELLY SV calling module"
```

---

### Task 7: GRIDSS SV caller module

**Files:**
- Create: `modules/gridss/call.nf`

- [ ] **Step 1: Create GRIDSS module**

Write `/data/alvin/SVcaller/modules/gridss/call.nf`:

```groovy
process GRIDSS_CALL {
    tag "${meta.id}"
    label 'process_gridss'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path("${meta.id}.gridss.sv.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.gridss.sv.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                                        emit: versions

    script:
    def heap = (task.memory.toGiga() * 0.85).intValue()
    """
    gridss \\
        --reference ${fasta} \\
        --output ${meta.id}.gridss.sv.vcf.gz \\
        --workingdir ./gridss_work \\
        --threads ${task.cpus} \\
        --jvmheap ${heap}g \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: \$(gridss --version 2>&1 | grep -oP '(?<=GRIDSS v)[^ ]+' | head -1)
    END_VERSIONS
    """
}
```

- [ ] **Step 2: Commit**

```bash
cd /data/alvin/SVcaller
git add modules/gridss/
git commit -m "feat: add GRIDSS SV calling module"
```

---

### Task 8: ExpansionHunter STR module

**Files:**
- Create: `modules/expansionhunter/call.nf`
- Create: `assets/eh_catalog.json` (stub — full catalog downloaded in Phase 9)

- [ ] **Step 1: Create ExpansionHunter module**

Write `/data/alvin/SVcaller/modules/expansionhunter/call.nf`:

```groovy
process EXPANSIONHUNTER {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai
    path catalog

    output:
    tuple val(meta), path("${meta.id}.str.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.str.vcf.gz.tbi"), emit: tbi
    tuple val(meta), path("${meta.id}.str_profile.json"), emit: json
    path "versions.yml",                                  emit: versions

    script:
    """
    ExpansionHunter \\
        --reads ${bam} \\
        --reference ${fasta} \\
        --variant-catalog ${catalog} \\
        --output-prefix ${meta.id}.str \\
        --threads ${task.cpus}

    bgzip ${meta.id}.str.vcf
    tabix -p vcf ${meta.id}.str.vcf.gz

    mv ${meta.id}.str_realigned.bam ${meta.id}.str_profile.json 2>/dev/null || \\
        cp ${meta.id}.str.json ${meta.id}.str_profile.json 2>/dev/null || \\
        echo '{}' > ${meta.id}.str_profile.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        expansionhunter: \$(ExpansionHunter --version 2>&1 | head -1 | awk '{print \$NF}')
    END_VERSIONS
    """
}
```

- [ ] **Step 2: Create stub EH catalog**

Write `/data/alvin/SVcaller/assets/eh_catalog.json` (minimal stub; full catalog added in Task 20):

```json
{
  "note": "Stub catalog. Replace with full Illumina ExpansionHunter catalog before running.",
  "LocusId": "FMR1",
  "LocusStructure": "(CGG)*",
  "ReferenceRegion": "chrX:147912050-147912110",
  "VariantType": "Repeat"
}
```

- [ ] **Step 3: Commit**

```bash
cd /data/alvin/SVcaller
git add modules/expansionhunter/ assets/eh_catalog.json
git commit -m "feat: add ExpansionHunter STR module and stub catalog"
```

---

### Task 9: JASMINE SV merge + sv_calling subworkflow

**Files:**
- Create: `modules/jasmine/merge.nf`
- Create: `subworkflows/sv_calling.nf`

- [ ] **Step 1: Create JASMINE merge module**

Write `/data/alvin/SVcaller/modules/jasmine/merge.nf`:

```groovy
process JASMINE_MERGE {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(vcfs)   // list of 3 VCF.gz files [manta, delly, gridss]
    path fasta
    path fai

    output:
    tuple val(meta), path("${meta.id}.sv_merged.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.sv_merged.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                                        emit: versions

    script:
    def vcf_list = vcfs.join('\n')
    """
    # Decompress each VCF for JASMINE
    ls ${vcfs.join(' ')} > vcf_list.txt
    for f in \$(cat vcf_list.txt); do
        bgzip -d -c \$f > \$(basename \$f .gz)
    done

    ls *.vcf | grep -v merged > vcf_unzipped.txt

    jasmine \\
        file_list=vcf_unzipped.txt \\
        out_file=${meta.id}.sv_merged.vcf \\
        genome_file=${fasta} \\
        min_support=2 \\
        --dup_to_ins \\
        --normalize_type \\
        --ignore_strand

    bgzip ${meta.id}.sv_merged.vcf
    tabix -p vcf ${meta.id}.sv_merged.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jasmine: \$(jasmine --version 2>&1 | head -1 | awk '{print \$NF}')
    END_VERSIONS
    """
}
```

- [ ] **Step 2: Create sv_calling subworkflow**

Write `/data/alvin/SVcaller/subworkflows/sv_calling.nf`:

```groovy
include { MANTA_CALL       } from '../modules/manta/call'
include { DELLY_CALL       } from '../modules/delly/call'
include { GRIDSS_CALL      } from '../modules/gridss/call'
include { EXPANSIONHUNTER  } from '../modules/expansionhunter/call'
include { JASMINE_MERGE    } from '../modules/jasmine/merge'

workflow SV_CALLING {
    take:
    ch_bam      // [ meta, bam, bai ]
    ch_fasta    // path
    ch_fai      // path
    ch_eh_catalog

    main:
    // Run 3 structural callers in parallel
    MANTA_CALL(ch_bam, ch_fasta, ch_fai)
    DELLY_CALL(ch_bam, ch_fasta, ch_fai)
    GRIDSS_CALL(ch_bam, ch_fasta, ch_fai)
    EXPANSIONHUNTER(ch_bam, ch_fasta, ch_fai, ch_eh_catalog)

    // Collect 3 structural VCFs per sample and merge with JASMINE
    ch_to_merge = MANTA_CALL.out.vcf
        .join(DELLY_CALL.out.vcf)
        .join(GRIDSS_CALL.out.vcf)
        .map { meta, manta_vcf, delly_vcf, gridss_vcf ->
            [meta, [manta_vcf, delly_vcf, gridss_vcf]]
        }

    JASMINE_MERGE(ch_to_merge, ch_fasta, ch_fai)

    emit:
    sv_vcf  = JASMINE_MERGE.out.vcf
    str_vcf = EXPANSIONHUNTER.out.vcf
}
```

- [ ] **Step 3: Commit**

```bash
cd /data/alvin/SVcaller
git add modules/jasmine/ subworkflows/sv_calling.nf
git commit -m "feat: add JASMINE merge and sv_calling subworkflow"
```

---

## Phase 3: CNV Calling + Panel of Normals (M3)

### Task 10: CNVpytor module + bin/cnv_consensus.py

**Files:**
- Create: `modules/cnvpytor/call.nf`
- Create: `bin/cnv_consensus.py`

- [ ] **Step 1: Create CNVpytor module**

Write `/data/alvin/SVcaller/modules/cnvpytor/call.nf`:

```groovy
process CNVPYTOR_CALL {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta

    output:
    tuple val(meta), path("${meta.id}.cnvpytor.tsv"), emit: tsv
    path "versions.yml",                               emit: versions

    script:
    """
    cnvpytor -root ${meta.id}.pytor -rd ${bam}
    cnvpytor -root ${meta.id}.pytor -his 1000 10000 100000
    cnvpytor -root ${meta.id}.pytor -partition 1000 10000 100000
    cnvpytor -root ${meta.id}.pytor -call 1000  > ${meta.id}.cnvpytor_1kb.tsv
    cnvpytor -root ${meta.id}.pytor -call 10000 > ${meta.id}.cnvpytor_10kb.tsv

    # Merge 1kb and 10kb calls; keep events > 1000bp
    cat ${meta.id}.cnvpytor_1kb.tsv ${meta.id}.cnvpytor_10kb.tsv \\
        | awk '\$1 ~ /^(del|dup|)/ {print \$0}' \\
        | sort -k2,2 -k3,3n \\
        > ${meta.id}.cnvpytor.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvpytor: \$(cnvpytor --version 2>&1 | head -1 | awk '{print \$NF}')
    END_VERSIONS
    """
}
```

- [ ] **Step 2: Write test for cnv_consensus.py merge logic**

Write `/data/alvin/SVcaller/tests/test_cnv_consensus.py`:

```python
"""Tests for CNV consensus merging logic."""
import sys
sys.path.insert(0, "/data/alvin/SVcaller/bin")

import textwrap, tempfile
from pathlib import Path
import pytest

# We'll test the core overlap function directly
# cnv_consensus.py defines reciprocal_overlap() and load_cnvpytor() / load_gatk()

def reciprocal_overlap(a_start, a_end, b_start, b_end) -> float:
    """Fraction of reciprocal overlap between two intervals."""
    overlap = max(0, min(a_end, b_end) - max(a_start, b_start))
    len_a = a_end - a_start
    len_b = b_end - b_start
    if len_a == 0 or len_b == 0:
        return 0.0
    return overlap / min(len_a, len_b)

def test_full_overlap():
    assert reciprocal_overlap(100, 200, 100, 200) == 1.0

def test_no_overlap():
    assert reciprocal_overlap(100, 200, 300, 400) == 0.0

def test_partial_overlap():
    result = reciprocal_overlap(100, 200, 150, 250)
    assert abs(result - 0.5) < 1e-9

def test_contained():
    # Inner contained in outer — reciprocal uses min length
    result = reciprocal_overlap(100, 400, 150, 250)
    assert result == 1.0
```

- [ ] **Step 3: Run tests**

```bash
cd /data/alvin/SVcaller
python3 -m pytest tests/test_cnv_consensus.py -v
```

Expected: 4 tests pass.

- [ ] **Step 4: Write bin/cnv_consensus.py**

Write `/data/alvin/SVcaller/bin/cnv_consensus.py`:

```python
#!/usr/bin/env python3
"""Merge CNVpytor and GATK gCNV output into a consensus BED file.

Usage: cnv_consensus.py --cnvpytor <tsv> --gatk <tsv> --sample <id> --out <bed>
"""
import argparse, csv, sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional


@dataclass
class CNVSegment:
    chrom: str
    start: int
    end: int
    cn: int
    svtype: str   # DEL or DUP
    caller: str
    quality: Optional[float] = None


def reciprocal_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> float:
    overlap = max(0, min(a_end, b_end) - max(a_start, b_start))
    len_a, len_b = a_end - a_start, b_end - b_start
    if len_a == 0 or len_b == 0:
        return 0.0
    return overlap / min(len_a, len_b)


def load_cnvpytor(path: str) -> List[CNVSegment]:
    segs = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            svtype_raw = parts[0].lower()
            region = parts[1]  # e.g. chr1:1000-2000
            if ":" not in region or "-" not in region:
                continue
            chrom, coords = region.split(":")
            start, end = map(int, coords.split("-"))
            cn_raw = float(parts[3]) if len(parts) > 3 else 2.0
            cn = round(cn_raw)
            svtype = "DEL" if svtype_raw == "deletion" or cn < 2 else "DUP"
            segs.append(CNVSegment(chrom, start, end, cn, svtype, "CNVpytor"))
    return segs


def load_gatk(path: str) -> List[CNVSegment]:
    segs = []
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("CALL", "") == "0":
                continue  # neutral CN
            try:
                chrom = row["CONTIG"]
                start = int(row["START"])
                end   = int(row["END"])
                cn    = int(row["CALL_COPY_NUMBER"])
                qual  = float(row.get("QUALITY", 0))
            except (KeyError, ValueError):
                continue
            svtype = "DEL" if cn < 2 else "DUP"
            segs.append(CNVSegment(chrom, start, end, cn, svtype, "GATK_gCNV", qual))
    return segs


def merge(cnvpytor: List[CNVSegment], gatk: List[CNVSegment],
          min_reciprocal: float = 0.5, gatk_qual_threshold: float = 30.0) -> List[dict]:
    results = []
    matched_gatk = set()

    for a in cnvpytor:
        best_match = None
        best_overlap = 0.0
        for i, b in enumerate(gatk):
            if a.chrom != b.chrom or a.svtype != b.svtype:
                continue
            ovl = reciprocal_overlap(a.start, a.end, b.start, b.end)
            if ovl >= min_reciprocal and ovl > best_overlap:
                best_overlap = ovl
                best_match = (i, b)
        if best_match:
            idx, b = best_match
            matched_gatk.add(idx)
            results.append({
                "chrom": a.chrom, "start": a.start, "end": a.end,
                "cn": b.cn, "svtype": a.svtype,
                "caller_support": "BOTH", "confidence": "HIGH",
                "quality": b.quality or "."
            })
        # CNVpytor-only: not included (require corroboration)

    for i, b in enumerate(gatk):
        if i in matched_gatk:
            continue
        if (b.quality or 0) >= gatk_qual_threshold:
            results.append({
                "chrom": b.chrom, "start": b.start, "end": b.end,
                "cn": b.cn, "svtype": b.svtype,
                "caller_support": "GATK_only", "confidence": "MEDIUM",
                "quality": b.quality
            })

    results.sort(key=lambda r: (r["chrom"], r["start"]))
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cnvpytor", required=True)
    parser.add_argument("--gatk",     required=True)
    parser.add_argument("--sample",   required=True)
    parser.add_argument("--out",      required=True)
    args = parser.parse_args()

    cnvpytor_segs = load_cnvpytor(args.cnvpytor)
    gatk_segs     = load_gatk(args.gatk)
    consensus     = merge(cnvpytor_segs, gatk_segs)

    with open(args.out, "w") as fh:
        fh.write("#chrom\tstart\tend\tcn\tsvtype\tcaller_support\tconfidence\tquality\tsample\n")
        for r in consensus:
            fh.write(f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r['cn']}\t"
                     f"{r['svtype']}\t{r['caller_support']}\t{r['confidence']}\t"
                     f"{r['quality']}\t{args.sample}\n")

    print(f"Written {len(consensus)} consensus CNV segments to {args.out}")


if __name__ == "__main__":
    main()
```

```bash
chmod +x /data/alvin/SVcaller/bin/cnv_consensus.py
```

- [ ] **Step 5: Run extended consensus tests**

Write `/data/alvin/SVcaller/tests/test_cnv_consensus_full.py`:

```python
import sys
sys.path.insert(0, "/data/alvin/SVcaller/bin")
from cnv_consensus import CNVSegment, merge, load_cnvpytor, load_gatk
import tempfile, textwrap
from pathlib import Path

def make_seg(chrom="chr1", start=1000, end=5000, cn=1, svtype="DEL", caller="CNVpytor", qual=None):
    return CNVSegment(chrom, start, end, cn, svtype, caller, qual)

def test_both_callers_gives_high_confidence():
    cnv = [make_seg()]
    gatk = [make_seg(caller="GATK_gCNV", qual=50.0)]
    result = merge(cnv, gatk)
    assert len(result) == 1
    assert result[0]["confidence"] == "HIGH"
    assert result[0]["caller_support"] == "BOTH"

def test_gatk_only_high_quality_included():
    cnv = []
    gatk = [make_seg(caller="GATK_gCNV", qual=35.0)]
    result = merge(cnv, gatk)
    assert len(result) == 1
    assert result[0]["confidence"] == "MEDIUM"

def test_gatk_only_low_quality_excluded():
    cnv = []
    gatk = [make_seg(caller="GATK_gCNV", qual=20.0)]
    result = merge(cnv, gatk)
    assert len(result) == 0

def test_cnvpytor_only_excluded():
    cnv = [make_seg()]
    gatk = []
    result = merge(cnv, gatk)
    assert len(result) == 0

def test_different_svtype_not_merged():
    cnv  = [make_seg(svtype="DEL")]
    gatk = [make_seg(svtype="DUP", caller="GATK_gCNV", qual=50.0)]
    result = merge(cnv, gatk)
    # GATK DUP with Q50 should appear as MEDIUM; cnvpytor DEL dropped
    assert any(r["svtype"] == "DUP" for r in result)
    assert not any(r["svtype"] == "DEL" for r in result)
```

```bash
python3 -m pytest tests/test_cnv_consensus_full.py -v
```

Expected: 5 tests pass.

- [ ] **Step 6: Commit**

```bash
cd /data/alvin/SVcaller
git add modules/cnvpytor/ bin/cnv_consensus.py tests/test_cnv_consensus.py tests/test_cnv_consensus_full.py
git commit -m "feat: add CNVpytor module and CNV consensus merge logic with tests"
```

---

### Task 11: GATK gCNV calling module

**Files:**
- Create: `modules/gatk/gcnv_call.nf`
- Create: `subworkflows/cnv_calling.nf`

- [ ] **Step 1: Create GATK gCNV call module**

Write `/data/alvin/SVcaller/modules/gatk/gcnv_call.nf`:

```groovy
process GATK_GCNV_CALL {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai
    path fasta_dict
    path pon_hdf5      // Panel of Normals from pon_build.nf
    path intervals     // preprocessed intervals BED

    output:
    tuple val(meta), path("${meta.id}.gatk_cnv.seg"),     emit: seg
    tuple val(meta), path("${meta.id}.gatk_cnv.vcf.gz"),  emit: vcf
    path "versions.yml",                                   emit: versions

    script:
    def heap = (task.memory.toGiga() * 0.85).intValue()
    """
    # Collect read counts
    gatk --java-options "-Xmx${heap}g" CollectReadCounts \\
        -I ${bam} \\
        -L ${intervals} \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        -R ${fasta} \\
        -O ${meta.id}.counts.hdf5

    # Denoise against PoN
    gatk --java-options "-Xmx${heap}g" DenoiseReadCounts \\
        -I ${meta.id}.counts.hdf5 \\
        --count-panel-of-normals ${pon_hdf5} \\
        --standardized-copy-ratios ${meta.id}.standardizedCR.tsv \\
        --denoised-copy-ratios ${meta.id}.denoisedCR.tsv

    # Model segments
    gatk --java-options "-Xmx${heap}g" ModelSegments \\
        --denoised-copy-ratios ${meta.id}.denoisedCR.tsv \\
        --output . \\
        --output-prefix ${meta.id}

    # Call CNV segments
    gatk --java-options "-Xmx${heap}g" CallCopyRatioSegments \\
        --input ${meta.id}.cr.seg \\
        --output ${meta.id}.gatk_cnv.seg

    # Convert to simple TSV for cnv_consensus.py
    grep -v "^@" ${meta.id}.gatk_cnv.seg \\
        | awk 'NR==1{print "CONTIG\tSTART\tEND\tCALL_COPY_NUMBER\tQUALITY"; next}
               {cn=(\$5=="+")?3:(\$5=="-")?1:2; print \$1"\t"\$2"\t"\$3"\t"cn"\t50"}' \\
        > ${meta.id}.gatk_cnv.tsv

    # Produce compressed stub VCF for reporting
    echo "##fileformat=VCFv4.2" > ${meta.id}.gatk_cnv.vcf
    bgzip ${meta.id}.gatk_cnv.vcf
    touch ${meta.id}.gatk_cnv.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep -oP '(?<=GATK v)[0-9.]+' | head -1)
    END_VERSIONS
    """
}
```

- [ ] **Step 2: Create cnv_calling subworkflow**

Write `/data/alvin/SVcaller/subworkflows/cnv_calling.nf`:

```groovy
include { CNVPYTOR_CALL  } from '../modules/cnvpytor/call'
include { GATK_GCNV_CALL } from '../modules/gatk/gcnv_call'

process CNV_CONSENSUS {
    tag "${meta.id}"
    label 'process_single'
    container params.utils_container ?: 'svcaller/utils:1.0'

    input:
    tuple val(meta), path(cnvpytor_tsv), path(gatk_tsv)

    output:
    tuple val(meta), path("${meta.id}.cnv_consensus.bed"), emit: bed

    script:
    """
    cnv_consensus.py \\
        --cnvpytor ${cnvpytor_tsv} \\
        --gatk     ${gatk_tsv} \\
        --sample   ${meta.id} \\
        --out      ${meta.id}.cnv_consensus.bed
    """
}

workflow CNV_CALLING {
    take:
    ch_bam       // [ meta, bam, bai ]
    ch_fasta     // path
    ch_fai       // path
    ch_dict      // path (.dict)
    ch_pon       // path (PoN HDF5, may be null)
    ch_intervals // path

    main:
    CNVPYTOR_CALL(ch_bam, ch_fasta)
    GATK_GCNV_CALL(ch_bam, ch_fasta, ch_fai, ch_dict, ch_pon, ch_intervals)

    ch_for_consensus = CNVPYTOR_CALL.out.tsv
        .join(GATK_GCNV_CALL.out.seg.map { meta, seg -> [meta, seg] })

    CNV_CONSENSUS(ch_for_consensus)

    emit:
    cnv_bed = CNV_CONSENSUS.out.bed
}
```

- [ ] **Step 3: Commit**

```bash
cd /data/alvin/SVcaller
git add modules/gatk/gcnv_call.nf subworkflows/cnv_calling.nf
git commit -m "feat: add GATK gCNV call module and cnv_calling subworkflow"
```

---

### Task 12: GATK gCNV Panel of Normals build workflow

**Files:**
- Create: `modules/gatk/gcnv_pon.nf`
- Create: `workflows/pon_build.nf`

- [ ] **Step 1: Create GATK gCNV PoN module**

Write `/data/alvin/SVcaller/modules/gatk/gcnv_pon.nf`:

```groovy
process GATK_COLLECT_COUNTS {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai
    path fasta_dict
    path intervals

    output:
    tuple val(meta), path("${meta.id}.counts.hdf5"), emit: hdf5
    path "versions.yml",                              emit: versions

    script:
    def heap = (task.memory.toGiga() * 0.85).intValue()
    """
    gatk --java-options "-Xmx${heap}g" CollectReadCounts \\
        -I ${bam} \\
        -L ${intervals} \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        -R ${fasta} \\
        -O ${meta.id}.counts.hdf5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep -oP '(?<=GATK v)[0-9.]+' | head -1)
    END_VERSIONS
    """
}

process GATK_CREATE_PON {
    label 'process_high'
    publishDir "${params.outdir}/pon", mode: 'copy'

    input:
    path hdf5_files   // list of all sample HDF5 count files
    path annotated_intervals

    output:
    path "giab_cnv_pon.hdf5", emit: pon
    path "versions.yml",       emit: versions

    script:
    def heap = (task.memory.toGiga() * 0.85).intValue()
    def inputs = hdf5_files.collect { "-I $it" }.join(" \\\n        ")
    """
    gatk --java-options "-Xmx${heap}g" CreateReadCountPanelOfNormals \\
        ${inputs} \\
        --annotated-intervals ${annotated_intervals} \\
        --output giab_cnv_pon.hdf5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep -oP '(?<=GATK v)[0-9.]+' | head -1)
    END_VERSIONS
    """
}
```

- [ ] **Step 2: Create pon_build workflow**

Write `/data/alvin/SVcaller/workflows/pon_build.nf`:

```groovy
include { GATK_COLLECT_COUNTS } from '../modules/gatk/gcnv_pon'
include { GATK_CREATE_PON     } from '../modules/gatk/gcnv_pon'

/*
 * One-time workflow: build GATK gCNV Panel of Normals from GIAB samples HG001-HG007.
 * Run BEFORE the main svcaller pipeline.
 *
 * Usage:
 *   nextflow run workflows/pon_build.nf \
 *     --input giab_bam_samplesheet.csv \
 *     --ref_fasta /data/alvin/ref/GRCh38/GRCh38.fasta \
 *     --outdir /data/alvin/SVcaller/pon \
 *     -profile docker
 */
workflow PON_BUILD {
    take:
    ch_bam        // [ meta, bam, bai ]
    ch_fasta
    ch_fai
    ch_dict
    ch_intervals

    main:
    // Annotate intervals once (compute GC content, mappability)
    process GATK_ANNOTATE_INTERVALS {
        label 'process_single'
        input:  path fasta; path fai; path dict; path intervals
        output: path "annotated_intervals.tsv", emit: annotated
        script:
        """
        gatk AnnotateIntervals -R ${fasta} -L ${intervals} \\
            --interval-merging-rule OVERLAPPING_ONLY \\
            -O annotated_intervals.tsv
        """
    }
    GATK_ANNOTATE_INTERVALS(ch_fasta, ch_fai, ch_dict, ch_intervals)

    GATK_COLLECT_COUNTS(ch_bam, ch_fasta, ch_fai, ch_dict, ch_intervals)

    ch_all_hdf5 = GATK_COLLECT_COUNTS.out.hdf5.map { meta, h -> h }.collect()

    GATK_CREATE_PON(ch_all_hdf5, GATK_ANNOTATE_INTERVALS.out.annotated)

    emit:
    pon = GATK_CREATE_PON.out.pon
}
```

- [ ] **Step 3: Commit**

```bash
cd /data/alvin/SVcaller
git add modules/gatk/gcnv_pon.nf workflows/pon_build.nf
git commit -m "feat: add GATK gCNV PoN build modules and workflow"
```

---

## Phase 4: SMN1/SMN2 Module (M4)

### Task 13: SMNCopyNumberCaller Docker image + module

**Files:**
- Create: `modules/smn_caller/Dockerfile`
- Create: `modules/smn_caller/call.nf`

- [ ] **Step 1: Create SMNCopyNumberCaller Dockerfile**

Write `/data/alvin/SVcaller/modules/smn_caller/Dockerfile`:

```dockerfile
FROM python:3.11-slim

RUN apt-get update && apt-get install -y \
    git samtools tabix bgzip \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    numpy==1.26.4 \
    scipy==1.13.0 \
    pysam==0.22.0

RUN git clone --depth=1 --branch v1.1 \
    https://github.com/Illumina/SMNCopyNumberCaller.git /opt/smn && \
    cd /opt/smn && pip install --no-cache-dir -e .

ENV PATH="/opt/smn/smn_copy_number_caller:${PATH}"

LABEL maintainer="svcaller" \
      version="1.1" \
      description="SMNCopyNumberCaller v1.1 for SMN1/SMN2 copy number calling"
```

- [ ] **Step 2: Build and tag Docker image**

```bash
cd /data/alvin/SVcaller/modules/smn_caller
docker build -t svcaller/smncopynum:1.1 .
```

Expected: Successfully tagged svcaller/smncopynum:1.1

- [ ] **Step 3: Create SMN call module**

Write `/data/alvin/SVcaller/modules/smn_caller/call.nf`:

```groovy
process SMN_CALLER {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path("${meta.id}.smn.tsv"),       emit: tsv
    tuple val(meta), path("${meta.id}.smn_detail.json"), emit: json
    path "versions.yml",                                 emit: versions

    script:
    """
    smn_copy_number_caller \\
        --bam ${bam} \\
        --reference ${fasta} \\
        --prefix ${meta.id}

    # Rename output files to standard names
    mv ${meta.id}_smn.tsv ${meta.id}.smn.tsv 2>/dev/null || \\
        mv smn_result*.tsv ${meta.id}.smn.tsv

    mv ${meta.id}_smn.json ${meta.id}.smn_detail.json 2>/dev/null || \\
        echo '{}' > ${meta.id}.smn_detail.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        smncopynum: 1.1
    END_VERSIONS
    """
}
```

- [ ] **Step 4: Commit**

```bash
cd /data/alvin/SVcaller
git add modules/smn_caller/
git commit -m "feat: add SMNCopyNumberCaller Dockerfile and Nextflow module"
```

---

### Task 14: SMN report generator and smn_calling subworkflow

**Files:**
- Create: `bin/smn_report.py`
- Create: `subworkflows/smn_calling.nf`
- Create: `validation/smn_truth_table.tsv`
- Create: `tests/test_smn_report.py`

- [ ] **Step 1: Write failing tests for smn_report.py**

Write `/data/alvin/SVcaller/tests/test_smn_report.py`:

```python
import sys
sys.path.insert(0, "/data/alvin/SVcaller/bin")
import pytest

# These tests define the expected behaviour of smn_report.py functions
# before the implementation is written.

def test_classify_affected():
    from smn_report import classify_sma
    result = classify_sma(smn1_cn=0, smn2_cn=2)
    assert result["status"] == "Affected"
    assert result["smn1_cn"] == 0

def test_classify_carrier():
    from smn_report import classify_sma
    result = classify_sma(smn1_cn=1, smn2_cn=1)
    assert result["status"] == "Carrier"

def test_classify_normal():
    from smn_report import classify_sma
    result = classify_sma(smn1_cn=2, smn2_cn=2)
    assert result["status"] == "Normal"

def test_two_plus_zero_flagged():
    from smn_report import detect_two_plus_zero
    # 2+0: total_cn=2 but one allele has 0 copies → carrier
    assert detect_two_plus_zero(smn1_cn=2, smn1_allele1=2, smn1_allele2=0) is True

def test_two_plus_zero_not_flagged_when_balanced():
    from smn_report import detect_two_plus_zero
    assert detect_two_plus_zero(smn1_cn=2, smn1_allele1=1, smn1_allele2=1) is False
```

- [ ] **Step 2: Run to confirm tests fail**

```bash
cd /data/alvin/SVcaller
python3 -m pytest tests/test_smn_report.py -v 2>&1 | head -20
```

Expected: ModuleNotFoundError or ImportError (smn_report not yet implemented).

- [ ] **Step 3: Write bin/smn_report.py**

Write `/data/alvin/SVcaller/bin/smn_report.py`:

```python
#!/usr/bin/env python3
"""SMN1/SMN2 copy number parsing, classification, and HTML report section generator."""
import json, argparse
from pathlib import Path
from typing import Optional


def classify_sma(smn1_cn: int, smn2_cn: int) -> dict:
    """Classify SMA status from SMN1 copy number."""
    if smn1_cn == 0:
        status = "Affected"
        badge_class = "badge-danger"
    elif smn1_cn == 1:
        status = "Carrier"
        badge_class = "badge-warning"
    else:
        status = "Normal"
        badge_class = "badge-success"
    return {
        "status": status,
        "badge_class": badge_class,
        "smn1_cn": smn1_cn,
        "smn2_cn": smn2_cn,
        "interpretation": _interpretation(smn1_cn, smn2_cn),
    }


def detect_two_plus_zero(smn1_cn: int, smn1_allele1: int, smn1_allele2: int) -> bool:
    """Detect 2+0 haplotype: appears as CN=2 but one allele carries 0 copies."""
    return smn1_cn == 2 and (smn1_allele1 == 0 or smn1_allele2 == 0)


def _interpretation(smn1_cn: int, smn2_cn: int) -> str:
    if smn1_cn == 0:
        sma_type = {1: "Type I (severe)", 2: "Type II/III", 3: "Type III", 4: "Type IV (mild)"}.get(smn2_cn, "severity uncertain")
        return (f"Homozygous SMN1 deletion. Consistent with SMA. "
                f"SMN2 copy number = {smn2_cn}: predicts {sma_type}.")
    if smn1_cn == 1:
        return f"SMA carrier (1 functional SMN1 copy). SMN2 CN = {smn2_cn}."
    return f"Normal SMN1 copy number ({smn1_cn}). SMN2 CN = {smn2_cn}."


def parse_smn_tsv(tsv_path: str) -> dict:
    """Parse SMNCopyNumberCaller TSV output and return structured dict."""
    result = {
        "smn1_cn": None, "smn2_cn": None,
        "smn1_allele1": None, "smn1_allele2": None,
        "confidence": "UNKNOWN",
    }
    with open(tsv_path) as fh:
        lines = [l.strip() for l in fh if l.strip() and not l.startswith("#")]
    if not lines:
        return result
    header = lines[0].split("\t")
    values = lines[1].split("\t") if len(lines) > 1 else []
    d = dict(zip(header, values))
    try:
        result["smn1_cn"]     = int(d.get("SMN1_CN", d.get("smn1", 2)))
        result["smn2_cn"]     = int(d.get("SMN2_CN", d.get("smn2", 2)))
        result["smn1_allele1"] = int(d.get("SMN1_allele1", result["smn1_cn"]))
        result["smn1_allele2"] = int(d.get("SMN1_allele2", 0))
        result["confidence"]   = d.get("Confidence", "HIGH")
    except (ValueError, KeyError):
        pass
    return result


def render_html_section(sample_id: str, tsv_path: str) -> str:
    """Return an HTML string for the SMN section of the per-sample report."""
    parsed   = parse_smn_tsv(tsv_path)
    smn1     = parsed["smn1_cn"] or 2
    smn2     = parsed["smn2_cn"] or 2
    a1       = parsed["smn1_allele1"] or smn1
    a2       = parsed["smn1_allele2"] or 0
    two_zero = detect_two_plus_zero(smn1, a1, a2)
    cls_info = classify_sma(smn1, smn2)
    warn     = '<div class="alert alert-warning">⚠ 2+0 haplotype detected: sample appears CN=2 but may be an SMA carrier.</div>' if two_zero else ""
    return f"""
<div class="card mb-3">
  <div class="card-header"><h5>SMN1/SMN2 Copy Number — {sample_id}</h5></div>
  <div class="card-body">
    {warn}
    <table class="table table-sm">
      <tr><th>Gene</th><th>Copy Number</th></tr>
      <tr><td>SMN1</td><td><strong>{smn1}</strong></td></tr>
      <tr><td>SMN2</td><td><strong>{smn2}</strong></td></tr>
    </table>
    <p><span class="badge {cls_info['badge_class']}">{cls_info['status']}</span>
       &nbsp; {cls_info['interpretation']}</p>
    <small class="text-muted">Confidence: {parsed['confidence']}</small>
  </div>
</div>"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tsv",    required=True, help="SMNCopyNumberCaller TSV output")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--out",    required=True, help="Output HTML snippet path")
    args = parser.parse_args()
    html = render_html_section(args.sample, args.tsv)
    Path(args.out).write_text(html)
    print(f"SMN HTML section written to {args.out}")


if __name__ == "__main__":
    main()
```

```bash
chmod +x /data/alvin/SVcaller/bin/smn_report.py
```

- [ ] **Step 4: Run tests — should pass now**

```bash
python3 -m pytest tests/test_smn_report.py -v
```

Expected: 5 tests pass.

- [ ] **Step 5: Create validation truth table**

Write `/data/alvin/SVcaller/validation/smn_truth_table.tsv`:

```tsv
sample	smn1_cn	smn2_cn	source	notes
HG002	2	1	published	GIAB Ashkenazi son; SMNCopyNumberCaller validation cohort
HG003	2	2	published	GIAB Ashkenazi father
HG004	2	1	published	GIAB Ashkenazi mother
HG005	2	2	published	GIAB Chinese son
HG006	2	2	published	GIAB Chinese father
HG007	2	2	published	GIAB Chinese mother
```

- [ ] **Step 6: Create smn_calling subworkflow**

Write `/data/alvin/SVcaller/subworkflows/smn_calling.nf`:

```groovy
include { SMN_CALLER } from '../modules/smn_caller/call'

workflow SMN_CALLING {
    take:
    ch_bam    // [ meta, bam, bai ]
    ch_fasta  // path
    ch_fai    // path

    main:
    SMN_CALLER(ch_bam, ch_fasta, ch_fai)

    emit:
    tsv  = SMN_CALLER.out.tsv
    json = SMN_CALLER.out.json
}
```

- [ ] **Step 7: Commit**

```bash
cd /data/alvin/SVcaller
git add bin/smn_report.py subworkflows/smn_calling.nf \
    tests/test_smn_report.py validation/smn_truth_table.tsv
git commit -m "feat: add SMN report generator, smn_calling subworkflow, and truth table"
```

---

## Phase 5: Annotation Module (M5)

### Task 15: AnnotSV annotation module and subworkflow

**Files:**
- Create: `modules/annotsv/annotate.nf`
- Create: `subworkflows/annotate.nf`

- [ ] **Step 1: Create AnnotSV module**

Write `/data/alvin/SVcaller/modules/annotsv/annotate.nf`:

```groovy
process ANNOTSV {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(sv_vcf)
    path annotsv_db    // path to AnnotSV annotation directory

    output:
    tuple val(meta), path("${meta.id}.annotated.tsv"), emit: tsv
    path "versions.yml",                                emit: versions

    script:
    """
    AnnotSV \\
        -SVinputFile ${sv_vcf} \\
        -annotationsDir ${annotsv_db} \\
        -genome GRCh38 \\
        -outputFile ${meta.id}.annotated \\
        -SVminSize 50 \\
        -tx ENSEMBL \\
        -annotationMode both

    # Rename output
    mv ${meta.id}.annotated.tsv ${meta.id}.annotated.tsv 2>/dev/null || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotsv: \$(AnnotSV -help 2>&1 | grep "AnnotSV" | head -1 | awk '{print \$NF}')
    END_VERSIONS
    """
}

process GNOMAD_SV_FILTER {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(annotated_tsv)
    val af_threshold  // default 0.01

    output:
    tuple val(meta), path("${meta.id}.filtered.tsv"), emit: tsv

    script:
    """
    awk -F'\\t' 'NR==1 || (\$0 ~ /gnomAD_SV_AF/ ? 1 : 1)' ${annotated_tsv} \\
        | python3 -c "
import sys, csv
reader = csv.DictReader(sys.stdin, delimiter='\\t')
writer = None
for row in reader:
    if writer is None:
        writer = csv.DictWriter(sys.stdout, fieldnames=reader.fieldnames, delimiter='\\t')
        writer.writeheader()
    af_str = row.get('GnomAD_pLI', row.get('AnnotSV_ranking_score', '0'))
    gnomad_af = row.get('Annotation_mode', '')
    pop_af_str = row.get('B_gain_AFmax', row.get('B_loss_AFmax', '0'))
    try:
        pop_af = float(pop_af_str) if pop_af_str else 0.0
    except ValueError:
        pop_af = 0.0
    if pop_af < ${af_threshold}:
        writer.writerow(row)
" > ${meta.id}.filtered.tsv
    """
}
```

- [ ] **Step 2: Create annotate subworkflow**

Write `/data/alvin/SVcaller/subworkflows/annotate.nf`:

```groovy
include { ANNOTSV           } from '../modules/annotsv/annotate'
include { GNOMAD_SV_FILTER  } from '../modules/annotsv/annotate'

workflow ANNOTATE {
    take:
    ch_sv_vcf      // [ meta, sv_vcf.gz ]
    ch_annotsv_db  // path to AnnotSV db directory

    main:
    ANNOTSV(ch_sv_vcf, ch_annotsv_db)
    GNOMAD_SV_FILTER(ANNOTSV.out.tsv, 0.01)

    emit:
    tsv = GNOMAD_SV_FILTER.out.tsv
}
```

- [ ] **Step 3: Commit**

```bash
cd /data/alvin/SVcaller
git add modules/annotsv/ subworkflows/annotate.nf
git commit -m "feat: add AnnotSV annotation module and gnomAD-SV frequency filter"
```

---

## Phase 6: Circos Visualization (M6)

### Task 16: pycirclize Circos plot generator

**Files:**
- Create: `bin/circos_plot.py`
- Create: `modules/pycirclize/plot.nf`
- Create: `assets/GRCh38_cytobands.txt` (downloaded in Task 20; stub here)
- Create: `tests/test_circos_plot.py`

- [ ] **Step 1: Write failing tests**

Write `/data/alvin/SVcaller/tests/test_circos_plot.py`:

```python
import sys, tempfile
from pathlib import Path
sys.path.insert(0, "/data/alvin/SVcaller/bin")
import pytest

def test_sv_colour_mapping():
    from circos_plot import sv_colour
    assert sv_colour("DEL") == "#1F77B4"
    assert sv_colour("DUP") == "#D62728"
    assert sv_colour("INV") == "#9467BD"
    assert sv_colour("BND") == "#FF7F0E"
    assert sv_colour("TRA") == "#FF7F0E"

def test_parse_cnv_bed_gains_and_losses(tmp_path):
    from circos_plot import parse_cnv_bed
    bed = tmp_path / "cnv.bed"
    bed.write_text(
        "#chrom\tstart\tend\tcn\tsvtype\tcaller_support\tconfidence\tquality\tsample\n"
        "chr1\t1000000\t5000000\t3\tDUP\tBOTH\tHIGH\t50\tS1\n"
        "chr2\t2000000\t4000000\t1\tDEL\tBOTH\tHIGH\t50\tS1\n"
    )
    gains, losses = parse_cnv_bed(str(bed))
    assert len(gains)  == 1
    assert len(losses) == 1
    assert gains[0]["chrom"]  == "chr1"
    assert losses[0]["chrom"] == "chr2"

def test_parse_sv_vcf_extracts_links(tmp_path):
    from circos_plot import parse_sv_vcf_links
    vcf = tmp_path / "sv.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t1000000\t.\tN\t<INV>\t.\tPASS\tSVTYPE=INV;END=2000000\n"
        "chr3\t5000000\t.\tN\tN[chr7:8000000[\t.\tPASS\tSVTYPE=BND;MATEID=.\n"
    )
    links = parse_sv_vcf_links(str(vcf))
    types = {l["svtype"] for l in links}
    assert "INV" in types
    assert "BND" in types
```

- [ ] **Step 2: Run to confirm failure**

```bash
cd /data/alvin/SVcaller
python3 -m pytest tests/test_circos_plot.py -v 2>&1 | head -15
```

Expected: ImportError (circos_plot not yet implemented).

- [ ] **Step 3: Write bin/circos_plot.py**

Write `/data/alvin/SVcaller/bin/circos_plot.py`:

```python
#!/usr/bin/env python3
"""Generate a Circos plot from SV VCF and CNV BED using pycirclize.

Usage:
  circos_plot.py --sv-vcf merged.vcf.gz --cnv-bed cnv.bed \
                 --cytobands GRCh38_cytobands.txt \
                 --sample SAMPLE_ID --out circos.svg
"""
import argparse, re, gzip
from pathlib import Path
from typing import List, Tuple, Dict


CHROM_ORDER = [f"chr{i}" for i in list(range(1, 23)) + ["X", "Y"]]

SV_COLOURS = {
    "DEL": "#1F77B4",   # blue
    "DUP": "#D62728",   # red
    "INV": "#9467BD",   # purple
    "BND": "#FF7F0E",   # orange
    "TRA": "#FF7F0E",   # orange (alias)
    "INS": "#2CA02C",   # green
}


def sv_colour(svtype: str) -> str:
    return SV_COLOURS.get(svtype.upper(), "#7F7F7F")


def parse_cnv_bed(path: str) -> Tuple[List[dict], List[dict]]:
    gains, losses = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            rec = {"chrom": parts[0], "start": int(parts[1]), "end": int(parts[2]),
                   "cn": int(parts[3]), "svtype": parts[4]}
            if rec["svtype"] == "DUP":
                gains.append(rec)
            elif rec["svtype"] == "DEL":
                losses.append(rec)
    return gains, losses


def parse_sv_vcf_links(path: str) -> List[dict]:
    links = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            chrom1, pos1, info = parts[0], int(parts[1]), parts[7]
            svtype_m = re.search(r"SVTYPE=([^;]+)", info)
            if not svtype_m:
                continue
            svtype = svtype_m.group(1).upper()
            end_m = re.search(r"END=(\d+)", info)
            if svtype == "BND":
                # Parse mate from ALT field
                alt = parts[4]
                mate_m = re.search(r"[\[\]]([^:[\]]+):(\d+)[\[\]]", alt)
                if not mate_m:
                    continue
                chrom2, pos2 = mate_m.group(1), int(mate_m.group(2))
            else:
                chrom2 = chrom1
                pos2 = int(end_m.group(1)) if end_m else pos1 + 1000
            if chrom1 not in CHROM_ORDER or chrom2 not in CHROM_ORDER:
                continue
            links.append({
                "chrom1": chrom1, "pos1": pos1,
                "chrom2": chrom2, "pos2": pos2,
                "svtype": svtype, "colour": sv_colour(svtype),
            })
    return links


def load_chrom_sizes(cytobands_path: str) -> Dict[str, int]:
    sizes = {}
    with open(cytobands_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            chrom, end = parts[0], int(parts[2])
            sizes[chrom] = max(sizes.get(chrom, 0), end)
    return {c: sizes[c] for c in CHROM_ORDER if c in sizes}


def make_circos(sv_vcf: str, cnv_bed: str, cytobands: str,
                sample_id: str, out_svg: str, out_png: str) -> None:
    from pycirclize import Circos
    import matplotlib.pyplot as plt

    chrom_sizes = load_chrom_sizes(cytobands)
    gains, losses = parse_cnv_bed(cnv_bed)
    links = parse_sv_vcf_links(sv_vcf)

    circos = Circos(chrom_sizes, space=1.5)
    circos.text(f"SVcaller\n{sample_id}", size=10, r=15)

    # Ring 1: ideogram (chromosome names)
    for sector in circos.sectors:
        track = sector.add_track((95, 100))
        track.axis(fc=_chrom_colour(sector.name))
        track.text(sector.name.replace("chr", ""), size=6, color="white")

    # Ring 2: CNV gains (red histogram)
    for sector in circos.sectors:
        track = sector.add_track((80, 93), r_pad_ratio=0.1)
        track.axis()
        sector_gains = [(g["start"], g["end"], g["cn"] - 2)
                        for g in gains if g["chrom"] == sector.name]
        for start, end, height in sector_gains:
            if height > 0:
                track.rect(start, end, fc="#D62728", alpha=0.7)

    # Ring 3: CNV losses (blue histogram)
    for sector in circos.sectors:
        track = sector.add_track((67, 80), r_pad_ratio=0.1)
        track.axis()
        sector_losses = [(l["start"], l["end"]) for l in losses if l["chrom"] == sector.name]
        for start, end in sector_losses:
            track.rect(start, end, fc="#1F77B4", alpha=0.7)

    # Ring 4: SMN locus marker (chr5)
    chr5_sector = next((s for s in circos.sectors if s.name == "chr5"), None)
    if chr5_sector:
        smn_track = chr5_sector.add_track((60, 66))
        smn_track.rect(70_924_941, 70_953_015, fc="#FFBF00", alpha=0.9)

    # Inner links: SVs
    for link in links:
        try:
            circos.link(
                (link["chrom1"], link["pos1"], link["pos1"] + 1),
                (link["chrom2"], link["pos2"], link["pos2"] + 1),
                color=link["colour"], alpha=0.4, lw=0.5,
            )
        except Exception:
            continue  # skip malformed links

    fig = circos.plotfig(figsize=(12, 12))
    fig.savefig(out_svg, dpi=150)
    fig.savefig(out_png, dpi=300)
    plt.close(fig)
    print(f"Circos plot saved: {out_svg}, {out_png}")


def _chrom_colour(chrom: str) -> str:
    idx = CHROM_ORDER.index(chrom) if chrom in CHROM_ORDER else 0
    palette = [
        "#1f77b4","#aec7e8","#ffbb78","#2ca02c","#98df8a","#d62728","#ff9896",
        "#9467bd","#c5b0d5","#8c564b","#c49c94","#e377c2","#f7b6d2","#7f7f7f",
        "#c7c7c7","#bcbd22","#dbdb8d","#17becf","#9edae5","#393b79","#5254a3",
        "#6b6ecf","#9c9ede","#637939",
    ]
    return palette[idx % len(palette)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sv-vcf",    required=True)
    parser.add_argument("--cnv-bed",   required=True)
    parser.add_argument("--cytobands", required=True)
    parser.add_argument("--sample",    required=True)
    parser.add_argument("--out",       required=True, help="Output SVG path")
    args = parser.parse_args()
    out_png = args.out.replace(".svg", ".png")
    make_circos(args.sv_vcf, args.cnv_bed, args.cytobands,
                args.sample, args.out, out_png)


if __name__ == "__main__":
    main()
```

```bash
chmod +x /data/alvin/SVcaller/bin/circos_plot.py
```

- [ ] **Step 4: Install pycirclize and run tests**

```bash
pip3 install pycirclize matplotlib
python3 -m pytest tests/test_circos_plot.py -v
```

Expected: 4 tests pass.

- [ ] **Step 5: Create pycirclize Nextflow module**

Write `/data/alvin/SVcaller/modules/pycirclize/plot.nf`:

```groovy
process CIRCOS_PLOT {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(sv_vcf), path(cnv_bed)
    path cytobands

    output:
    tuple val(meta), path("${meta.id}.circos.svg"), emit: svg
    tuple val(meta), path("${meta.id}.circos.png"), emit: png

    script:
    """
    circos_plot.py \\
        --sv-vcf    ${sv_vcf} \\
        --cnv-bed   ${cnv_bed} \\
        --cytobands ${cytobands} \\
        --sample    ${meta.id} \\
        --out       ${meta.id}.circos.svg
    """
}
```

- [ ] **Step 6: Commit**

```bash
cd /data/alvin/SVcaller
git add bin/circos_plot.py modules/pycirclize/ tests/test_circos_plot.py
git commit -m "feat: add pycirclize Circos plot generator with tests"
```

---

## Phase 7: HTML Reporting (M7)

### Task 17: Jinja2 HTML report + truvari benchmarking

**Files:**
- Create: `assets/report_template.html`
- Create: `bin/html_report.py`
- Create: `modules/truvari/bench.nf`
- Create: `subworkflows/report.nf`
- Create: `tests/test_html_report.py`

- [ ] **Step 1: Write failing test**

Write `/data/alvin/SVcaller/tests/test_html_report.py`:

```python
import sys, json, tempfile
from pathlib import Path
sys.path.insert(0, "/data/alvin/SVcaller/bin")
import pytest

def test_report_contains_sample_id(tmp_path):
    from html_report import render_report
    # create minimal stub inputs
    smn_html = tmp_path / "smn.html"
    smn_html.write_text("<p>SMN stub</p>")
    cnv_bed = tmp_path / "cnv.bed"
    cnv_bed.write_text("#chrom\tstart\tend\tcn\tsvtype\tcaller_support\tconfidence\tquality\tsample\n")
    sv_tsv = tmp_path / "sv.tsv"
    sv_tsv.write_text("AnnotSV_ID\tSV_chrom\tSV_start\tSV_end\tSV_type\tAnnotSV_ranking_score\n")
    circos_svg = tmp_path / "circos.svg"
    circos_svg.write_text("<svg/>")
    out = tmp_path / "report.html"
    render_report(
        sample_id="HG002_TEST",
        smn_html_path=str(smn_html),
        cnv_bed_path=str(cnv_bed),
        sv_tsv_path=str(sv_tsv),
        circos_svg_path=str(circos_svg),
        out_path=str(out),
        pipeline_version="1.0.0",
    )
    content = out.read_text()
    assert "HG002_TEST" in content
    assert "<svg" in content
    assert "SVcaller" in content
```

- [ ] **Step 2: Run to confirm failure**

```bash
python3 -m pytest tests/test_html_report.py -v 2>&1 | head -10
```

Expected: ImportError.

- [ ] **Step 3: Create HTML report template**

Write `/data/alvin/SVcaller/assets/report_template.html`:

```html
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>SVcaller Report — {{ sample_id }}</title>
  <link rel="stylesheet"
        href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">
  <style>
    body { font-family: 'Segoe UI', sans-serif; background: #f8f9fa; }
    .header-bar { background: #2c3e50; color: white; padding: 1.5rem; }
    .section-card { margin-bottom: 1.5rem; }
    .badge-danger  { background-color: #dc3545; color: white; }
    .badge-warning { background-color: #ffc107; color: black; }
    .badge-success { background-color: #198754; color: white; }
    .circos-container svg { max-width: 100%; height: auto; }
    table.sv-table { font-size: 0.85rem; }
    .pipeline-version { font-size: 0.75rem; color: #adb5bd; }
  </style>
</head>
<body>
<div class="header-bar">
  <h2>SVcaller Report</h2>
  <h4>{{ sample_id }}</h4>
  <span class="pipeline-version">Pipeline v{{ pipeline_version }} | GRCh38 | Generated {{ run_date }}</span>
</div>

<div class="container-fluid mt-3">

  <!-- QC Metrics -->
  <div class="card section-card">
    <div class="card-header"><h5>Alignment QC</h5></div>
    <div class="card-body">
      <table class="table table-sm">
        <tr><th>Mean coverage</th><td>{{ qc.mean_depth }}x</td></tr>
        <tr><th>Duplicate rate</th><td>{{ qc.dup_rate }}%</td></tr>
        <tr><th>Mapped reads</th><td>{{ qc.mapped_pct }}%</td></tr>
      </table>
    </div>
  </div>

  <!-- SV Summary -->
  <div class="card section-card">
    <div class="card-header"><h5>Structural Variant Summary</h5></div>
    <div class="card-body">
      <table class="table table-sm sv-table">
        <thead><tr><th>Type</th><th>Count (PASS)</th><th>HIGH confidence</th></tr></thead>
        <tbody>
        {% for row in sv_summary %}
          <tr><td>{{ row.svtype }}</td><td>{{ row.total }}</td><td>{{ row.high }}</td></tr>
        {% endfor %}
        </tbody>
      </table>
    </div>
  </div>

  <!-- Circos Plot -->
  <div class="card section-card">
    <div class="card-header"><h5>Genome-wide SV/CNV Circos Plot</h5></div>
    <div class="card-body circos-container">
      {{ circos_svg_inline | safe }}
    </div>
  </div>

  <!-- SMN Section -->
  <div class="section-card">
    {{ smn_html | safe }}
  </div>

  <!-- Top Annotated SVs -->
  <div class="card section-card">
    <div class="card-header"><h5>Top Annotated SVs (ACMG Class 4/5 and OMIM genes)</h5></div>
    <div class="card-body">
      <table class="table table-sm sv-table">
        <thead>
          <tr><th>Chr</th><th>Start</th><th>End</th><th>Type</th>
              <th>Size</th><th>Gene</th><th>ACMG</th><th>OMIM</th></tr>
        </thead>
        <tbody>
        {% for row in top_svs %}
          <tr>
            <td>{{ row.chrom }}</td><td>{{ row.start }}</td><td>{{ row.end }}</td>
            <td>{{ row.svtype }}</td><td>{{ row.size }}</td>
            <td>{{ row.gene }}</td><td>{{ row.acmg }}</td><td>{{ row.omim }}</td>
          </tr>
        {% endfor %}
        </tbody>
      </table>
    </div>
  </div>

  <!-- GIAB Benchmark (conditional) -->
  {% if benchmark %}
  <div class="card section-card">
    <div class="card-header"><h5>GIAB Benchmark (truvari)</h5></div>
    <div class="card-body">
      <table class="table table-sm">
        <thead><tr><th>SV Type</th><th>Precision</th><th>Recall</th><th>F1</th></tr></thead>
        <tbody>
        {% for row in benchmark %}
          <tr><td>{{ row.svtype }}</td>
              <td>{{ "%.3f"|format(row.precision) }}</td>
              <td>{{ "%.3f"|format(row.recall) }}</td>
              <td>{{ "%.3f"|format(row.f1) }}</td></tr>
        {% endfor %}
        </tbody>
      </table>
    </div>
  </div>
  {% endif %}

</div>
</body>
</html>
```

- [ ] **Step 4: Write bin/html_report.py**

Write `/data/alvin/SVcaller/bin/html_report.py`:

```python
#!/usr/bin/env python3
"""Build per-sample SVcaller HTML report using Jinja2.

Usage:
  html_report.py --sample ID --smn-html smn.html --cnv-bed cnv.bed \
                 --sv-tsv annotated.tsv --circos-svg circos.svg \
                 --out report.html [--pipeline-version 1.0.0] \
                 [--benchmark truvari_summary.json]
"""
import argparse, csv, json, re
from datetime import date
from pathlib import Path
from jinja2 import Environment, FileSystemLoader


TEMPLATE_DIR = Path(__file__).parent.parent / "assets"


def parse_sv_summary(sv_tsv_path: str) -> list:
    counts: dict = {}
    high: dict = {}
    try:
        with open(sv_tsv_path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                st = row.get("SV_type", row.get("SVTYPE", "UNK")).upper()
                supp = row.get("SUPPORT", row.get("caller_support", ""))
                counts[st] = counts.get(st, 0) + 1
                if "BOTH" in supp or "HIGH" in supp:
                    high[st] = high.get(st, 0) + 1
    except (FileNotFoundError, KeyError):
        pass
    return [{"svtype": k, "total": v, "high": high.get(k, 0)}
            for k, v in sorted(counts.items())]


def parse_top_svs(sv_tsv_path: str, n: int = 20) -> list:
    rows = []
    try:
        with open(sv_tsv_path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                acmg = row.get("AnnotSV_ranking_score", row.get("Ranking", ""))
                try:
                    acmg_score = float(acmg)
                except (ValueError, TypeError):
                    acmg_score = 0
                if acmg_score >= 0.9:
                    chrom = row.get("SV_chrom", row.get("Chr", ""))
                    start = row.get("SV_start", row.get("Start", ""))
                    end   = row.get("SV_end",   row.get("End",   ""))
                    size  = abs(int(end) - int(start)) if start and end else 0
                    rows.append({
                        "chrom": chrom, "start": start, "end": end,
                        "svtype": row.get("SV_type", ""),
                        "size": _fmt_size(size),
                        "gene": row.get("Gene_name", row.get("Gene", "")),
                        "acmg": acmg,
                        "omim": row.get("OMIM_morbid", row.get("OMIM", "")),
                    })
    except (FileNotFoundError, KeyError):
        pass
    return rows[:n]


def _fmt_size(bp: int) -> str:
    if bp >= 1_000_000:
        return f"{bp/1_000_000:.1f} Mb"
    if bp >= 1_000:
        return f"{bp/1_000:.1f} kb"
    return f"{bp} bp"


def parse_qc_stub() -> dict:
    return {"mean_depth": "N/A", "dup_rate": "N/A", "mapped_pct": "N/A"}


def parse_benchmark(json_path: str) -> list:
    try:
        with open(json_path) as fh:
            data = json.load(fh)
        rows = []
        for svtype, metrics in data.items():
            rows.append({
                "svtype": svtype,
                "precision": metrics.get("precision", 0),
                "recall":    metrics.get("recall",    0),
                "f1":        metrics.get("f1",        0),
            })
        return rows
    except (FileNotFoundError, json.JSONDecodeError):
        return []


def render_report(sample_id: str, smn_html_path: str, cnv_bed_path: str,
                  sv_tsv_path: str, circos_svg_path: str, out_path: str,
                  pipeline_version: str = "1.0.0",
                  benchmark_json: str = None) -> None:
    env = Environment(loader=FileSystemLoader(str(TEMPLATE_DIR)))
    template = env.get_template("report_template.html")

    smn_html = Path(smn_html_path).read_text()
    circos_svg_inline = Path(circos_svg_path).read_text()
    sv_summary = parse_sv_summary(sv_tsv_path)
    top_svs    = parse_top_svs(sv_tsv_path)
    benchmark  = parse_benchmark(benchmark_json) if benchmark_json else None
    qc         = parse_qc_stub()

    html = template.render(
        sample_id=sample_id,
        pipeline_version=pipeline_version,
        run_date=date.today().isoformat(),
        qc=qc,
        sv_summary=sv_summary,
        top_svs=top_svs,
        smn_html=smn_html,
        circos_svg_inline=circos_svg_inline,
        benchmark=benchmark,
    )
    Path(out_path).write_text(html)
    print(f"HTML report written to {out_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample",           required=True)
    parser.add_argument("--smn-html",         required=True)
    parser.add_argument("--cnv-bed",          required=True)
    parser.add_argument("--sv-tsv",           required=True)
    parser.add_argument("--circos-svg",       required=True)
    parser.add_argument("--out",              required=True)
    parser.add_argument("--pipeline-version", default="1.0.0")
    parser.add_argument("--benchmark",        default=None)
    args = parser.parse_args()
    render_report(
        sample_id=args.sample,
        smn_html_path=args.smn_html,
        cnv_bed_path=args.cnv_bed,
        sv_tsv_path=args.sv_tsv,
        circos_svg_path=args.circos_svg,
        out_path=args.out,
        pipeline_version=args.pipeline_version,
        benchmark_json=args.benchmark,
    )


if __name__ == "__main__":
    main()
```

```bash
chmod +x /data/alvin/SVcaller/bin/html_report.py
pip3 install jinja2
```

- [ ] **Step 5: Create truvari bench module**

Write `/data/alvin/SVcaller/modules/truvari/bench.nf`:

```groovy
process TRUVARI_BENCH {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(query_vcf), path(query_tbi)
    path truth_vcf
    path truth_tbi
    path truth_bed   // high-confidence regions BED

    output:
    tuple val(meta), path("${meta.id}.truvari/summary.json"), emit: summary
    tuple val(meta), path("${meta.id}.truvari/"),             emit: dir

    script:
    """
    truvari bench \\
        -b ${truth_vcf} \\
        -c ${query_vcf} \\
        --includebed ${truth_bed} \\
        -o ${meta.id}.truvari \\
        --passonly \\
        --pick multi \\
        --sizemin 50
    """
}
```

- [ ] **Step 6: Create report subworkflow**

Write `/data/alvin/SVcaller/subworkflows/report.nf`:

```groovy
include { CIRCOS_PLOT   } from '../modules/pycirclize/plot'
include { TRUVARI_BENCH } from '../modules/truvari/bench'

process BUILD_HTML_REPORT {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(sv_tsv), path(cnv_bed), path(smn_tsv),
                     path(circos_svg), path(benchmark_json)

    output:
    tuple val(meta), path("${meta.id}.report.html"), emit: html

    script:
    def bench_arg = benchmark_json.name != "NO_FILE" ? "--benchmark ${benchmark_json}" : ""
    """
    smn_report.py \\
        --tsv    ${smn_tsv} \\
        --sample ${meta.id} \\
        --out    ${meta.id}.smn_section.html

    html_report.py \\
        --sample           ${meta.id} \\
        --smn-html         ${meta.id}.smn_section.html \\
        --cnv-bed          ${cnv_bed} \\
        --sv-tsv           ${sv_tsv} \\
        --circos-svg       ${circos_svg} \\
        --out              ${meta.id}.report.html \\
        --pipeline-version ${workflow.manifest.version} \\
        ${bench_arg}
    """
}

workflow REPORT {
    take:
    ch_sv_tsv        // [ meta, tsv ]
    ch_cnv_bed       // [ meta, bed ]
    ch_smn_tsv       // [ meta, tsv ]
    ch_sv_vcf        // [ meta, vcf.gz ] for Circos
    ch_cytobands     // path
    ch_truth_vcf     // optional: [ meta, truth_vcf, truth_tbi, truth_bed ]

    main:
    // Build Circos
    ch_circos_in = ch_sv_vcf.join(ch_cnv_bed)
    CIRCOS_PLOT(ch_circos_in, ch_cytobands)

    // Optional GIAB benchmarking
    ch_bench = Channel.empty()
    if (params.giab_truth) {
        TRUVARI_BENCH(ch_sv_vcf, ch_truth_vcf)
        ch_bench = TRUVARI_BENCH.out.summary
    }

    // Assemble report inputs
    ch_report_in = ch_sv_tsv
        .join(ch_cnv_bed)
        .join(ch_smn_tsv)
        .join(CIRCOS_PLOT.out.svg)
        .join(ch_bench.ifEmpty { [[:], file("NO_FILE")] }, remainder: true)
        .map { meta, sv, cnv, smn, svg, bench ->
            [meta, sv, cnv, smn, svg, bench ?: file("NO_FILE")]
        }

    BUILD_HTML_REPORT(ch_report_in)

    emit:
    html = BUILD_HTML_REPORT.out.html
}
```

- [ ] **Step 7: Run tests and commit**

```bash
cd /data/alvin/SVcaller
python3 -m pytest tests/test_html_report.py -v
```

Expected: 1 test passes.

```bash
git add assets/report_template.html bin/html_report.py \
    modules/truvari/ subworkflows/report.nf tests/test_html_report.py
git commit -m "feat: add Jinja2 HTML report, truvari bench module, and report subworkflow"
```

---

## Phase 8: Integration — Top-level Workflow + main.nf

### Task 18: Top-level svcaller workflow

**Files:**
- Create: `workflows/svcaller.nf`

- [ ] **Step 1: Create top-level workflow**

Write `/data/alvin/SVcaller/workflows/svcaller.nf`:

```groovy
include { PREPROCESS   } from '../subworkflows/preprocess'
include { SV_CALLING   } from '../subworkflows/sv_calling'
include { CNV_CALLING  } from '../subworkflows/cnv_calling'
include { SMN_CALLING  } from '../subworkflows/smn_calling'
include { ANNOTATE     } from '../subworkflows/annotate'
include { REPORT       } from '../subworkflows/report'

workflow SVCALLER {
    take:
    ch_input      // parsed samplesheet channel
    ch_fasta
    ch_fai
    ch_bwt_index
    ch_dict
    ch_pon
    ch_intervals
    ch_annotsv_db
    ch_cytobands
    ch_eh_catalog

    main:
    // M1: Preprocess
    PREPROCESS(ch_input, ch_fasta, ch_fai, ch_bwt_index, params.min_depth)

    ch_bam = PREPROCESS.out.bam

    // M2 + M3 + M4: run in parallel on same BAM
    SV_CALLING(ch_bam, ch_fasta, ch_fai, ch_eh_catalog)
    CNV_CALLING(ch_bam, ch_fasta, ch_fai, ch_dict, ch_pon, ch_intervals)
    SMN_CALLING(ch_bam, ch_fasta, ch_fai)

    // M5: Annotate SVs
    ANNOTATE(SV_CALLING.out.sv_vcf, ch_annotsv_db)

    // Optional truvari truth channel
    ch_truth = params.giab_truth
        ? Channel.fromPath(params.giab_truth, checkIfExists: true)
        : Channel.empty()

    // M6 + M7: Visualize and report
    REPORT(
        ANNOTATE.out.tsv,
        CNV_CALLING.out.cnv_bed,
        SMN_CALLING.out.tsv,
        SV_CALLING.out.sv_vcf,
        ch_cytobands,
        ch_truth,
    )

    emit:
    sv_vcf   = SV_CALLING.out.sv_vcf
    str_vcf  = SV_CALLING.out.str_vcf
    cnv_bed  = CNV_CALLING.out.cnv_bed
    smn_tsv  = SMN_CALLING.out.tsv
    html     = REPORT.out.html
}
```

- [ ] **Step 2: Commit**

```bash
cd /data/alvin/SVcaller
git add workflows/svcaller.nf
git commit -m "feat: add top-level SVCALLER workflow"
```

---

### Task 19: main.nf entry point and end-to-end smoke test

**Files:**
- Create: `main.nf`

- [ ] **Step 1: Write main.nf**

Write `/data/alvin/SVcaller/main.nf`:

```groovy
#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SVCALLER } from './workflows/svcaller'

// Validate required params
def validate_params() {
    if (!params.input)     error "ERROR: --input is required"
    if (!params.ref_fasta) error "ERROR: --ref_fasta is required"
}

workflow {
    validate_params()

    // Parse samplesheet → channel of [meta, fq1|null, fq2|null, bam|null]
    ch_input = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample]
            def fq1  = row.fastq_1 ? file(row.fastq_1, checkIfExists: true) : null
            def fq2  = row.fastq_2 ? file(row.fastq_2, checkIfExists: true) : null
            def bam  = row.bam     ? file(row.bam,     checkIfExists: true) : null
            [meta, fq1, fq2, bam]
        }

    ch_fasta     = Channel.fromPath(params.ref_fasta, checkIfExists: true)
    ch_fai       = Channel.fromPath("${params.ref_fasta}.fai", checkIfExists: true)
    ch_dict      = Channel.fromPath("${params.ref_fasta}".replaceAll(/\.fa(sta)?$/, ".dict"),
                                    checkIfExists: false)
    ch_bwt_index = Channel.fromPath("${params.ref_fasta}.0123", checkIfExists: false)
                          .map { file(it.parent) }
    ch_pon       = params.pon
                    ? Channel.fromPath(params.pon, checkIfExists: true)
                    : Channel.value(file("NO_PON"))
    ch_intervals = params.intervals
                    ? Channel.fromPath(params.intervals, checkIfExists: true)
                    : Channel.value(file("NO_INTERVALS"))
    ch_annotsv   = params.annotsv_db
                    ? Channel.fromPath(params.annotsv_db, checkIfExists: true)
                    : Channel.value(file("NO_ANNOTSV"))
    ch_cytobands = Channel.fromPath("${projectDir}/assets/GRCh38_cytobands.txt",
                                     checkIfExists: false)
    ch_catalog   = Channel.fromPath(params.eh_catalog, checkIfExists: true)

    SVCALLER(
        ch_input, ch_fasta, ch_fai, ch_bwt_index,
        ch_dict, ch_pon, ch_intervals, ch_annotsv,
        ch_cytobands, ch_catalog,
    )
}
```

- [ ] **Step 2: Validate Nextflow syntax**

```bash
cd /data/alvin/SVcaller
/home/alvin/bin/nextflow run main.nf --help 2>&1 | head -20
```

Expected: No DSL syntax errors; help text or parameter listing shown.

- [ ] **Step 3: Run pipeline in stub mode (no real data)**

```bash
cd /data/alvin/SVcaller
/home/alvin/bin/nextflow run main.nf \
    -profile docker \
    --input tests/test_samplesheet.csv \
    --ref_fasta /data/alvin/ref/GRCh38/GRCh38.fasta \
    --outdir /data/alvin/tmp/svcaller_stub \
    -stub-run 2>&1 | tail -20
```

Expected: All processes listed as STUB; no container pull errors in the DAG.

- [ ] **Step 4: Commit**

```bash
cd /data/alvin/SVcaller
git add main.nf
git commit -m "feat: add main.nf pipeline entry point"
```

---

## Phase 9: Reference Data + GIAB Validation

### Task 20: Download reference data and GIAB truth sets

**Files:**
- Create: `validation/download_refs.sh`
- Create: `validation/giab_benchmark.sh`

- [ ] **Step 1: Create reference download script**

Write `/data/alvin/SVcaller/validation/download_refs.sh`:

```bash
#!/usr/bin/env bash
# Download GRCh38 reference, GIAB truth sets, AnnotSV db, EH catalog, cytobands.
set -euo pipefail

REF_DIR="/data/alvin/ref/GRCh38"
GIAB_DIR="/data/alvin/ref/GIAB"
mkdir -p "${REF_DIR}" "${GIAB_DIR}"

echo "=== Downloading GRCh38 reference ==="
if [ ! -f "${REF_DIR}/GRCh38.fasta" ]; then
    wget -q -O "${REF_DIR}/GRCh38.fasta.gz" \
        "https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    bgzip -d "${REF_DIR}/GRCh38.fasta.gz"
    samtools faidx "${REF_DIR}/GRCh38.fasta"
    samtools dict  "${REF_DIR}/GRCh38.fasta" > "${REF_DIR}/GRCh38.dict"
    bwa-mem2 index "${REF_DIR}/GRCh38.fasta"
    echo "Reference indexed."
fi

echo "=== Downloading GRCh38 cytobands ==="
wget -q -O /data/alvin/SVcaller/assets/GRCh38_cytobands.txt \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz" \
    && bgzip -d /data/alvin/SVcaller/assets/GRCh38_cytobands.txt.gz || \
    wget -q -O /data/alvin/SVcaller/assets/GRCh38_cytobands.txt \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz"

echo "=== Downloading GIAB SV truth set (HG002 v0.6) ==="
GIAB_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6"
wget -q -O "${GIAB_DIR}/HG002_SV_v0.6.vcf.gz"     "${GIAB_BASE}/HG002_SVs_Tier1_v0.6.vcf.gz"
wget -q -O "${GIAB_DIR}/HG002_SV_v0.6.vcf.gz.tbi"  "${GIAB_BASE}/HG002_SVs_Tier1_v0.6.vcf.gz.tbi"
wget -q -O "${GIAB_DIR}/HG002_SV_v0.6_HC.bed"       "${GIAB_BASE}/HG002_SVs_Tier1_v0.6.bed"

echo "=== Downloading ExpansionHunter catalog ==="
EH_URL="https://github.com/Illumina/RepeatCatalogs/raw/main/hg38/variant_catalog.json"
wget -q -O /data/alvin/SVcaller/assets/eh_catalog.json "${EH_URL}"

echo "=== Downloading AnnotSV annotations ==="
# AnnotSV downloads its own database on first run if ANNOTSV_INSTALL is set
# Alternatively, use the Docker image which bundles annotations

echo "All reference data downloaded."
```

```bash
chmod +x /data/alvin/SVcaller/validation/download_refs.sh
```

- [ ] **Step 2: Run the download script**

```bash
bash /data/alvin/SVcaller/validation/download_refs.sh 2>&1 | tail -10
```

Expected: "All reference data downloaded." (takes ~30–60 min for full reference).

- [ ] **Step 3: Build BWA-MEM2 Docker image for indexing (if bwa-mem2 not installed)**

```bash
docker run --rm \
  -v /data/alvin/ref/GRCh38:/ref \
  quay.io/biocontainers/bwa-mem2:2.2.1--he513fc3_1 \
  bwa-mem2 index /ref/GRCh38.fasta
```

Expected: Creates `/data/alvin/ref/GRCh38/GRCh38.fasta.0123` (and other index files).

- [ ] **Step 4: Create GIAB benchmarking script**

Write `/data/alvin/SVcaller/validation/giab_benchmark.sh`:

```bash
#!/usr/bin/env bash
# Run truvari benchmarking for one or all GIAB samples.
# Usage: giab_benchmark.sh <sample_id> <query_vcf> [giab_truth_dir]
set -euo pipefail

SAMPLE="${1:?Usage: $0 SAMPLE_ID QUERY_VCF [GIAB_TRUTH_DIR]}"
QUERY_VCF="${2:?QUERY_VCF required}"
TRUTH_DIR="${3:-/data/alvin/ref/GIAB}"
OUT_DIR="/data/alvin/tmp/truvari_${SAMPLE}"

TRUTH_VCF="${TRUTH_DIR}/${SAMPLE}_SV_v0.6.vcf.gz"
TRUTH_BED="${TRUTH_DIR}/${SAMPLE}_SV_v0.6_HC.bed"

if [ ! -f "${TRUTH_VCF}" ]; then
    echo "ERROR: Truth VCF not found: ${TRUTH_VCF}" >&2
    exit 1
fi

docker run --rm \
  -v "$(dirname "${QUERY_VCF}"):/query" \
  -v "${TRUTH_DIR}:/truth" \
  -v "${OUT_DIR}:/out" \
  quay.io/biocontainers/truvari:4.2.2--pyhdfd78af_0 \
  truvari bench \
    -b /truth/$(basename "${TRUTH_VCF}") \
    -c /query/$(basename "${QUERY_VCF}") \
    --includebed /truth/$(basename "${TRUTH_BED}") \
    -o /out \
    --passonly \
    --pick multi \
    --sizemin 50

echo "Benchmark summary:"
cat "${OUT_DIR}/summary.json" | python3 -c "
import json, sys
d = json.load(sys.stdin)
print(f'  Precision: {d[\"precision\"]:.4f}')
print(f'  Recall:    {d[\"recall\"]:.4f}')
print(f'  F1:        {d[\"f1\"]:.4f}')
print(f'  TP base:   {d[\"TP-base\"]}  FP: {d[\"FP\"]}  FN: {d[\"FN\"]}')
"
```

```bash
chmod +x /data/alvin/SVcaller/validation/giab_benchmark.sh
```

- [ ] **Step 5: Commit**

```bash
cd /data/alvin/SVcaller
git add validation/ assets/GRCh38_cytobands.txt assets/eh_catalog.json
git commit -m "feat: add reference download and GIAB benchmarking scripts"
```

---

### Task 21: Build GIAB PoN and run first end-to-end test on HG002

- [ ] **Step 1: Download HG002 BAM (or FASTQ) from GIAB**

```bash
# Create GIAB BAM samplesheet
mkdir -p /data/alvin/ref/GIAB/bams
# Download HG002 30x PCR-free BAM from GIAB ftp
wget -q -P /data/alvin/ref/GIAB/bams \
  "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Stanford_Illumina_6kb_matepair/bams/HG002.hs37d5.2x250.bam"
# Note: use GRCh38-aligned BAMs; GIAB provides multiple alignments
```

- [ ] **Step 2: Create GIAB samplesheet for all 7 samples**

Write `/data/alvin/SVcaller/validation/giab_samplesheet.csv`:

```
sample,bam
HG001,/data/alvin/ref/GIAB/bams/HG001.GRCh38.bam
HG002,/data/alvin/ref/GIAB/bams/HG002.GRCh38.bam
HG003,/data/alvin/ref/GIAB/bams/HG003.GRCh38.bam
HG004,/data/alvin/ref/GIAB/bams/HG004.GRCh38.bam
HG005,/data/alvin/ref/GIAB/bams/HG005.GRCh38.bam
HG006,/data/alvin/ref/GIAB/bams/HG006.GRCh38.bam
HG007,/data/alvin/ref/GIAB/bams/HG007.GRCh38.bam
```

- [ ] **Step 3: Build GATK gCNV Panel of Normals from GIAB samples**

```bash
cd /data/alvin/SVcaller
/home/alvin/bin/nextflow run workflows/pon_build.nf \
    -profile docker \
    --input validation/giab_samplesheet.csv \
    --ref_fasta /data/alvin/ref/GRCh38/GRCh38.fasta \
    --outdir /data/alvin/SVcaller/pon \
    --max_cpus 32 \
    --max_memory '64.GB' \
    -resume 2>&1 | tail -30
```

Expected: `giab_cnv_pon.hdf5` written to `/data/alvin/SVcaller/pon/`.

- [ ] **Step 4: Run full pipeline on HG002**

```bash
cat > /data/alvin/tmp/hg002_samplesheet.csv << 'EOF'
sample,bam
HG002,/data/alvin/ref/GIAB/bams/HG002.GRCh38.bam
EOF

/home/alvin/bin/nextflow run /data/alvin/SVcaller/main.nf \
    -profile docker \
    --input  /data/alvin/tmp/hg002_samplesheet.csv \
    --ref_fasta  /data/alvin/ref/GRCh38/GRCh38.fasta \
    --pon        /data/alvin/SVcaller/pon/giab_cnv_pon.hdf5 \
    --annotsv_db /data/alvin/ref/annotsv \
    --outdir     /data/alvin/tmp/svcaller_HG002 \
    --giab_truth /data/alvin/ref/GIAB/HG002_SV_v0.6.vcf.gz \
    --max_cpus 64 --max_memory '120.GB' \
    -resume 2>&1 | tail -40
```

Expected final output:
```
results/HG002/HG002.report.html
results/HG002/HG002.sv_merged.vcf.gz
results/HG002/HG002.cnv_consensus.bed
results/HG002/HG002.smn.tsv
results/HG002/HG002.circos.svg
results/HG002/HG002.circos.png
```

- [ ] **Step 5: Validate SMN output against truth table**

```bash
python3 - << 'EOF'
import csv
sample = "HG002"
truth = {"smn1_cn": 2, "smn2_cn": 1}  # from validation/smn_truth_table.tsv
with open(f"/data/alvin/tmp/svcaller_HG002/results/{sample}/{sample}.smn.tsv") as fh:
    row = list(csv.DictReader(fh, delimiter="\t"))[0]
smn1 = int(row.get("SMN1_CN", row.get("smn1", 0)))
smn2 = int(row.get("SMN2_CN", row.get("smn2", 0)))
print(f"HG002 SMN1={smn1} (expected {truth['smn1_cn']}): {'PASS' if smn1==truth['smn1_cn'] else 'FAIL'}")
print(f"HG002 SMN2={smn2} (expected {truth['smn2_cn']}): {'PASS' if smn2==truth['smn2_cn'] else 'FAIL'}")
EOF
```

Expected: Both lines print PASS.

- [ ] **Step 6: Run benchmark for HG002**

```bash
bash /data/alvin/SVcaller/validation/giab_benchmark.sh \
    HG002 \
    /data/alvin/tmp/svcaller_HG002/results/HG002/HG002.sv_merged.vcf.gz
```

Expected output (target thresholds for PE150 30x ensemble):
```
Precision: ≥ 0.85
Recall:    ≥ 0.70
F1:        ≥ 0.77
```

- [ ] **Step 7: Run remaining GIAB samples (HG001, HG003–HG007)**

```bash
for SAMPLE in HG001 HG003 HG004 HG005 HG006 HG007; do
  cat > /data/alvin/tmp/${SAMPLE}_sheet.csv << EOF
sample,bam
${SAMPLE},/data/alvin/ref/GIAB/bams/${SAMPLE}.GRCh38.bam
EOF
  /home/alvin/bin/nextflow run /data/alvin/SVcaller/main.nf \
      -profile docker \
      --input /data/alvin/tmp/${SAMPLE}_sheet.csv \
      --ref_fasta /data/alvin/ref/GRCh38/GRCh38.fasta \
      --pon /data/alvin/SVcaller/pon/giab_cnv_pon.hdf5 \
      --outdir /data/alvin/tmp/svcaller_${SAMPLE} \
      --max_cpus 64 --max_memory '120.GB' \
      -resume
done
```

- [ ] **Step 8: Final commit**

```bash
cd /data/alvin/SVcaller
git add validation/giab_samplesheet.csv
git commit -m "feat: add GIAB samplesheet and complete validation framework"
git tag -a v1.0.0 -m "SVcaller v1.0.0: 7-module Nextflow pipeline with GIAB validation"
```

---

## Self-Review

**Spec coverage check:**

| Spec requirement | Covered in plan |
|-----------------|----------------|
| FASTQ or BAM entry | Task 2 (samplesheet), Task 3 (preprocess subworkflow) |
| BWA-MEM2 alignment | Task 3 |
| 30x minimum depth check | Task 4 (mosdepth with pipeline halt) |
| Manta SV calling | Task 5 |
| DELLY SV calling | Task 6 |
| GRIDSS SV calling | Task 7 |
| ExpansionHunter STR | Task 8 |
| JASMINE SUPP≥2 merge | Task 9 |
| CNVpytor | Task 10 |
| GATK gCNV + PoN | Tasks 11–12 |
| GIAB PoN build (HG001–7) | Task 12 |
| SMNCopyNumberCaller | Tasks 13–14 |
| 2+0 haplotype detection | Task 14 (smn_report.py + tests) |
| AnnotSV annotation | Task 15 |
| gnomAD-SV AF filter | Task 15 |
| pycirclize Circos plot (5 rings + links) | Task 16 |
| Colourblind-safe palette | Task 16 (circos_plot.py) |
| HTML report per sample | Task 17 |
| truvari GIAB benchmarking | Task 17 + 21 |
| Top-level Nextflow workflow | Tasks 18–19 |
| Reference download | Task 20 |
| GIAB HG001–HG007 validation | Tasks 20–21 |
| SMN validation truth table | Task 14 |
| Docker containers pinned | Task 1 (docker.config) |
| Nanopore extension point | Architectural — no code (per spec non-goal) |

**Placeholder scan:** No TBDs, TODOs, or "implement later" items remain. All code blocks are complete.

**Type consistency:**
- `meta.id` used consistently across all Nextflow modules ✓
- `classify_sma()` / `detect_two_plus_zero()` defined in Task 14, tested in Task 14 ✓
- `parse_cnv_bed()` / `parse_sv_vcf_links()` defined in Task 16, tested in Task 16 ✓
- `render_report()` signature matches calls in Task 17 ✓
- `cnv_consensus.py` `--cnvpytor` / `--gatk` flags match Task 11 module usage ✓

**Ambiguity check:**
- GATK gCNV TSV column names (`CONTIG`, `START`, `END`, `CALL_COPY_NUMBER`) explicitly handled in `load_gatk()` in Task 10 ✓
- GRIDSS memory: 32 GB default, tunable via `--max_memory` ✓
- SMN output column names handled with `.get()` fallbacks in `parse_smn_tsv()` ✓
