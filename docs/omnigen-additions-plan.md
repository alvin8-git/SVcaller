# OmniGen Additions — Copy-Number / Deletion Trait Contracts

Implementation plan for four new per-sample contract files consumed by OmniGen
(downstream WGS interpretation). SVcaller *calls* from the aligned BAM; OmniGen
*renders*. Each feature emits a small, stable TSV that OmniGen reads.

> Status: **IMPLEMENTED** (code in the working tree; not yet committed).
> All deterministic parts are covered by `tests/test_cnv_traits.py` (20 tests,
> passing in-session). The real per-sample depth numbers require one full
> pipeline run — see "Running for real depth numbers" at the end of this doc.
>
> Files added: `assets/cnv_trait_regions.bed`, `bin/cnv_traits_common.py`,
> `bin/rh_status.py`, `bin/amy1_cn.py`, `bin/gst_null.py`, `bin/lpa_kiv2.py`,
> `modules/traits/depth.nf`, `subworkflows/cnv_traits.nf`,
> `tests/test_cnv_traits.py`, `docs/cnv_traits_card_demo.html`.
> Files edited: `main.nf`, `workflows/svcaller.nf`, `subworkflows/report.nf`,
> `conf/docker.config`, `bin/html_report.py`, `assets/report_template.html`,
> `README.md`.

---

## 0. Repository facts this plan is built on

Architecture (Nextflow DSL2):

```
main.nf
 └─ workflows/svcaller.nf   (workflow SVCALLER)
     ├─ subworkflows/preprocess.nf     PREPROCESS  → BAM, mosdepth summary, regions.bed.gz
     ├─ subworkflows/sv_calling.nf     SV_CALLING  (Manta/DELLY/GRIDSS/MELT/SvABA, EH, STRling)
     ├─ subworkflows/cnv_calling.nf    CNV_CALLING (CNVpytor + GATK-gCNV → CNV_CONSENSUS)
     ├─ subworkflows/smn_calling.nf    SMN_CALLING (single-locus targeted caller)  ← MODEL
     ├─ subworkflows/annotate.nf       ANNOTATE
     └─ subworkflows/report.nf         REPORT      (BUILD_HTML_REPORT → bin/html_report.py)
```

Key patterns observed:

* **Contract emitter model = SMN.** `modules/smn_caller/call.nf` takes
  `tuple(meta,bam,bai), fasta, fai`, runs a per-locus tool, and does
  `publishDir "${params.outdir}/${meta.id}", pattern: "*.smn.tsv"`. This is the
  exact shape our four features should copy.
* **CNV consensus.** `subworkflows/cnv_calling.nf` process `CNV_CONSENSUS`
  (container `params.utils_container ?: 'svcaller/utils:1.2'`) runs
  `bin/cnv_consensus.py` → `${id}.cnv_consensus.bed` with columns
  `chrom,start,end,cn,svtype,caller_support,confidence,quality,sample`.
* **Targeted depth tooling already present.** `modules/mosdepth/coverage.nf`
  runs mosdepth (container in `conf/docker.config`:
  `MOSDEPTH → quay.io/biocontainers/mosdepth:0.3.14--h05c3d44_0`). Genome-wide
  today; we add a `--by <bed>` targeted variant. `samtools:1.23.1` biocontainer
  also available.
* **Per-sample HTML report** is produced by `subworkflows/report.nf` →
  process `BUILD_HTML_REPORT` (container `svcaller/utils:1.2`), which runs
  `bin/html_report.py` with many optional `--flag` args guarded by
  `NO_*` sentinel filenames, then `publishDir "${params.outdir}/${meta.id}", pattern: "*.report.html"`.
  html_report.py: `argparse` in `main()` (~L1353), rendering in
  `render_report()` (~L1251). All Python tools live in `bin/` and are put on
  PATH inside the container via `export PATH=${projectDir}/bin:$PATH`.
* **Actual report path** is `results/<S>/<S>.report.html` (single sample dir),
  not `results/<S>/<S>/<S>.report.html`. Plan targets the real path.

### 0.1 Empirical finding that shapes the design (in-session check vs real HG002 bed)

Intersecting `results/HG002/HG002.cnv_consensus.bed` (327 rows; 325 DUP / 2 DEL;
almost all `cn=2`) against GRCh38 trait loci gave:

| Locus | GRCh38 window checked | Overlapping consensus rows |
|-------|-----------------------|----------------------------|
| RHD   | chr1:25,272,393–25,330,445 | **0** |
| AMY1 cluster | chr1:103,571,000–103,760,000 | **0** |
| GSTM1 | chr1:109,687,814–109,693,020 | **0** |
| GSTT1 | chr22:24,376,133–24,384,680 | **0** |
| LPA   | chr6:160,531,483–160,664,275 | 1 (DUP) |

Consequence: **the existing SV/CNV consensus does NOT capture the homozygous
whole-gene deletions or high-copy loci we need** (HG002 is a known GSTM1-null yet
the consensus shows no DEL there; the consensus is dominated by `cn=2` DUP
segments). Therefore:

* CNV-consensus intersection is used only as a **corroborating** signal, never
  the primary call.
* The **primary** signal for all four features is a **targeted, normalized
  read-depth** estimate computed directly from the BAM. This is also why the
  genome-wide SV ensemble misses them: RHD/RHCE, GSTM1 paralogs, the AMY1 tandem
  array, and the LPA KIV-2 repeat are segmental-duplication / paralogous regions
  where short-read callers under-perform.

---

## 1. Shared infrastructure (build once, all four features reuse)

### 1.1 New assets file: `assets/cnv_trait_regions.bed`

A checked-in BED (GRCh38, `chr`-prefixed to match `hg38.canonical.fa`) with two
kinds of rows, 4th column = region label:

```
# trait target regions
chr1    25272393    25330445    RHD
chr1    103571000   103760000   AMY1_CLUSTER
chr1    109687814   109693020   GSTM1
chr22   24376133    24384680    GSTT1
chr6    160605000   160650000   LPA_KIV2      # KIV-2 repeat block within LPA
# copy-number-stable single-copy control regions for depth normalization
chr1    ...         ...         CTRL_1
chr2    ...         ...         CTRL_2
# (~10 diploid, non-segdup, non-repeat control windows spread across autosomes)
```

**Implementer must verify every coordinate against RefSeq/GENCODE for the exact
assembly in `params.ref_fasta` before committing** (values above are
GRCh38 anchors, minus-strand genes give descending exon order). Add the
KIV-2 repeat-unit sub-window separately from the whole LPA gene. Pick ~10
control windows that are diploid and outside segmental duplications (e.g. use
copy-number-stable housekeeping gene bodies).

Wire as a value channel in `main.nf` (mirror `ch_cytobands`):

```groovy
ch_trait_regions = Channel.value(file("${projectDir}/assets/cnv_trait_regions.bed", checkIfExists: true))
```

Pass it through `SVCALLER(...)` into the new subworkflow.

### 1.2 New module: `modules/traits/depth.nf` (process `TRAIT_DEPTH`)

One targeted-depth pass per sample over all trait + control regions.

```groovy
process TRAIT_DEPTH {
    tag "${meta.id}"
    label 'process_low'
    // container assigned in conf/docker.config: mosdepth 0.3.14 biocontainer

    input:
    tuple val(meta), path(bam), path(bai)
    path  regions_bed          // assets/cnv_trait_regions.bed

    output:
    tuple val(meta), path("${meta.id}.trait_depth.regions.bed.gz"), emit: depth
    path "versions.yml",                                            emit: versions

    script:
    """
    mosdepth --threads ${task.cpus} --no-per-base --mapq 0 \\
        --by ${regions_bed} ${meta.id}.trait_depth ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/mosdepth //')
    END_VERSIONS
    """
}
```

Notes:
* `--mapq 0` (do **not** MAPQ-filter): AMY1 array copies and LPA KIV-2 repeat
  units, and paralogous RHD/GST reads, are multi-mapping; filtering on MAPQ
  destroys the copy-number signal. Normalization against control regions
  (also `--mapq 0`) cancels the multi-mapping bias to first order.
* Output `*.regions.bed.gz` has one row per region with mean depth in col 5 —
  exactly what the interpreters below consume.
* Register the container in `conf/docker.config` and `conf/singularity.config`:
  `withName: 'TRAIT_DEPTH' { container = 'quay.io/biocontainers/mosdepth:0.3.14--h05c3d44_0' }`.

### 1.3 New subworkflow: `subworkflows/cnv_traits.nf` (workflow `CNV_TRAITS`)

Runs `TRAIT_DEPTH` once, then four small interpreter processes (each a thin
wrapper over a `bin/*.py`, container `svcaller/utils:1.2`). Each interpreter
receives the depth bed AND the CNV consensus bed (for corroboration).

```groovy
include { TRAIT_DEPTH } from '../modules/traits/depth'

process RH_STATUS { /* publishDir "${params.outdir}/${meta.id}/bloodgroup", pattern "*rh_status.tsv" */ }
process AMY1_CN   { /* publishDir "${params.outdir}/${meta.id}/cnv_traits",  pattern "*amy1.tsv"      */ }
process GST_NULL  { /* publishDir "${params.outdir}/${meta.id}/cnv_traits",  pattern "*gst_null.tsv"  */ }
process LPA_KIV2  { /* publishDir "${params.outdir}/${meta.id}/cnv_traits",  pattern "*lpa_kiv2.tsv"  */ }

workflow CNV_TRAITS {
    take:
    ch_bam           // [ meta, bam, bai ]
    ch_trait_regions // path
    ch_cnv_bed       // [ meta, cnv_consensus.bed ]  (from CNV_CALLING.out.cnv_bed)

    main:
    TRAIT_DEPTH(ch_bam, ch_trait_regions)
    ch_in = TRAIT_DEPTH.out.depth.join(ch_cnv_bed)   // [ meta, depth_bed, cnv_bed ]
    RH_STATUS(ch_in)
    AMY1_CN(ch_in)
    GST_NULL(ch_in)
    LPA_KIV2(ch_in)

    emit:
    rh    = RH_STATUS.out.tsv
    amy1  = AMY1_CN.out.tsv
    gst   = GST_NULL.out.tsv
    lpa   = LPA_KIV2.out.tsv
}
```

Each interpreter process body follows the `CNV_CONSENSUS` shape, e.g.:

```groovy
process RH_STATUS {
    tag "${meta.id}"; label 'process_single'
    container params.utils_container ?: 'svcaller/utils:1.2'
    publishDir "${params.outdir}/${meta.id}/bloodgroup", mode: 'copy', pattern: "*rh_status.tsv"
    input:  tuple val(meta), path(depth_bed), path(cnv_bed)
    output: tuple val(meta), path("${meta.id}.rh_status.tsv"), emit: tsv
    script:
    """
    export PATH=${projectDir}/bin:\$PATH
    rh_status.py --depth ${depth_bed} --cnv-bed ${cnv_bed} --sample ${meta.id} \\
                 --out ${meta.id}.rh_status.tsv
    """
}
```

> Optimization option: collapse the four interpreters into ONE process
> `BUILD_CNV_TRAITS` running one `bin/cnv_traits.py` that emits all four files
> (multiple `publishDir` blocks by pattern). Fewer processes, one shared
> depth-normalization computation. Recommended once schemas are settled; four
> separate scripts are easier to unit-test first.

### 1.4 Shared normalization helper (in `bin/cnv_traits_common.py`)

All interpreters import one function:

```python
def region_depths(depth_bed_path):
    """Return {label: mean_depth} from mosdepth --by regions.bed.gz (col5=mean)."""
def control_depth(depths):
    """Median of CTRL_* regions → single-copy-per-haplotype baseline."""
def copies(region_label, depths, ploidy_per_copy=1.0):
    """round( ploidy * mean(region) / control_depth ). control≈2 haploid copies."""
```

Diploid single-copy baseline: `control_depth` reflects 2 gene copies. So
`estimated_copies = 2 * region_mean / control_median`.

---

## 2. Feature-by-feature detail

### Feature 1 — Rh factor (RHD presence/absence)

* **(a) Slots into:** new `CNV_TRAITS` subworkflow, process `RH_STATUS`,
  `bin/rh_status.py`. Depends on `TRAIT_DEPTH` (primary) + `CNV_CALLING.out.cnv_bed`
  (corroboration).
* **(b) Tool + container:** mosdepth (`TRAIT_DEPTH`, mosdepth:0.3.14) for depth;
  interpreter in `svcaller/utils:1.2` (pure-Python stdlib, no pysam needed).
* **Method:** `RHD_copies = round(2 * depth(RHD) / control_median)`. The common
  Rh-negative haplotype is a homozygous whole-RHD deletion →
  `RHD_copies == 0` → `Rh_status = neg`; `>=1` → `pos`. Corroborate: if the CNV
  consensus has a DEL fully spanning the RHD window, raise confidence; if depth
  says 0 but consensus disagrees, lower confidence. **RHCE caveat:** RHD and RHCE
  are near-identical paralogs; keep the RHD window on RHD-specific exons and note
  residual RHCE cross-mapping can inflate a true-0 toward ~0.3 copies — treat
  `RHD_copies < 0.5` as deletion.
* **(c) Contract:** `results/<S>/bloodgroup/rh_status.tsv`

  | column | type | notes |
  |--------|------|-------|
  | `sample` | str | meta.id |
  | `RHD_copies` | int | rounded normalized copies (0,1,2,…) |
  | `Rh_status` | enum | `pos` \| `neg` |
  | `confidence` | enum | `HIGH`\|`MEDIUM`\|`LOW` (depth/consensus concordance) |

  Header row `#sample\tRHD_copies\tRh_status\tconfidence`, one data row.
* **(d) HTML report:** see §3 — appears in new "Blood group & CN traits" card.
* **(e) Testing:** consensus cross-check logic testable in-session now (HG002 →
  0 DELs at RHD → not a homozygous deletion → present/`pos`, consistent with
  HG002 being Rh D positive). Depth number needs a run (§4).

### Feature 2 — AMY1 copy number

* **(a) Slots into:** `CNV_TRAITS`, process `AMY1_CN`, `bin/amy1_cn.py`.
  Targeted read-depth ONLY (SV ensemble cannot resolve a 0–20+ tandem array;
  consensus showed 0 rows here).
* **(b) Tool + container:** mosdepth `TRAIT_DEPTH` + utils interpreter.
* **Method:** `AMY1_copies = round(2 * depth(AMY1_CLUSTER) / control_median)`.
  The AMY1 window must cover only the AMY1 sub-array (exclude AMY2A/AMY2B, which
  are separate diploid genes) or normalize accordingly — document exact window.
  Report the rounded integer; also keep the raw ratio internally for QC.
* **(c) Contract:** `results/<S>/cnv_traits/amy1.tsv`

  | column | type | notes |
  |--------|------|-------|
  | `sample` | str | |
  | `AMY1_copies` | int | diploid total copies (0–~20) |
  | `method` | str | `read-depth-normalized` |

  Header `#sample\tAMY1_copies\tmethod`.
* **(d) HTML report:** §3.
* **(e) Testing:** interpreter arithmetic unit-testable in-session with a
  synthetic depth bed; real copy number needs a run.

### Feature 3 — GSTM1 / GSTT1 null

* **(a) Slots into:** `CNV_TRAITS`, process `GST_NULL`, `bin/gst_null.py`.
  Primary = targeted depth; corroborate with consensus DEL (none present for
  HG002 despite known GSTM1-null → confirms depth must lead).
* **(b) Tool + container:** mosdepth + utils interpreter.
* **Method:** per gene `ratio = depth(gene)/control_median`. Homozygous deletion
  (`null`) when `ratio < 0.15`; `present` otherwise (heterozygous ~0.5 still
  `present`). GSTM1 has paralogs (GSTM2–5) nearby → keep window on GSTM1-unique
  exons. Optional third state not required by contract (report binary).
* **(c) Contract:** `results/<S>/cnv_traits/gst_null.tsv`

  | column | type | notes |
  |--------|------|-------|
  | `sample` | str | |
  | `GSTM1` | enum | `null` \| `present` |
  | `GSTT1` | enum | `null` \| `present` |

  Header `#sample\tGSTM1\tGSTT1`.
* **(d) HTML report:** §3.
* **(e) Testing:** consensus cross-check + arithmetic testable in-session
  (real bed → present/present from consensus alone, exposing that depth is
  required to catch HG002's true GSTM1-null → good regression assertion). Real
  ratios need a run.

### Feature 4 — Lp(a) / LPA KIV-2 repeat copy number

* **(a) Slots into:** `CNV_TRAITS`, process `LPA_KIV2`, `bin/lpa_kiv2.py`.
  Targeted repeat-copy estimate (tandem VNTR; SV ensemble/consensus cannot size
  it — only 1 spurious DUP overlap seen).
* **(b) Tool + container:** mosdepth `TRAIT_DEPTH` (the KIV-2 repeat-unit window
  in the assets bed) + utils interpreter. Depth-ratio method; a future upgrade
  could add a k-mer/repeat-aware caller but is out of scope.
* **Method:** `KIV2_copies = round(2 * depth(LPA_KIV2) / control_median)`.
  Because all KIV-2 units are near-identical, mosdepth over the single repeat
  window (with `--mapq 0`) sums multi-mapped coverage from every copy → depth
  scales linearly with total copy number; dividing by the single-copy control
  baseline yields diploid KIV-2 copies (biological range ~10–100). Report total
  diploid copies (document if per-allele is wanted later).
* **(c) Contract:** `results/<S>/cnv_traits/lpa_kiv2.tsv`

  | column | type | notes |
  |--------|------|-------|
  | `sample` | str | |
  | `KIV2_copies` | int | total diploid KIV-2 repeat copies |
  | `method` | str | `read-depth-ratio` |

  Header `#sample\tKIV2_copies\tmethod`.
* **(d) HTML report:** §3.
* **(e) Testing:** arithmetic unit-testable in-session; real copies need a run.

---

## 3. HTML report integration (`bin/html_report.py` + `subworkflows/report.nf`)

Add one new report card "Blood Group & Copy-Number Traits" that renders all four
contract files. Steps:

1. **`bin/html_report.py`**
   * Add args in `main()` (mirror the `NO_*` optional pattern):
     `--rh-status`, `--amy1`, `--gst-null`, `--lpa-kiv2` (each `default=None`).
   * Add a `build_cnv_traits_section(rh_path, amy1_path, gst_path, lpa_path)`
     returning an HTML card (a small table: Rh D status, AMY1 copies, GSTM1/GSTT1
     null status, Lp(a) KIV-2 copies, each with confidence where available).
     Gracefully render "not available" when a path is `None` or a `NO_FILE`
     sentinel. Insert its output inside `render_report()` (~L1251) alongside the
     existing SMN/CNV sections; thread the four new params through the
     `render_report(...)` call in `main()`.

2. **`subworkflows/report.nf`**
   * `workflow REPORT` — add four `take:` channels after `ch_strling_tsv`:
     `ch_rh_status`, `ch_amy1`, `ch_gst_null`, `ch_lpa_kiv2` (`[meta, tsv]`).
   * Join them into `ch_report_in` using the existing optional-join idiom
     (`.join(ch_x, remainder: true).filter{ it[1] != null }.map{ ... ?: file("NO_FILE") }`).
   * `process BUILD_HTML_REPORT` — add four `path(...)` inputs to the input
     tuple; add four `def rh_arg = rh_tsv.name != "NO_FILE" ? "--rh-status ${rh_tsv}" : ""`
     (etc.); append `${rh_arg} ${amy1_arg} ${gst_arg} ${lpa_arg}` to the
     `html_report.py` invocation.

3. **`workflows/svcaller.nf`**
   * `include { CNV_TRAITS } from '../subworkflows/cnv_traits'`.
   * After `CNV_CALLING(...)`, add
     `CNV_TRAITS(ch_bam, ch_trait_regions, CNV_CALLING.out.cnv_bed)`.
   * Extend the `REPORT(...)` call with the four new channels
     `CNV_TRAITS.out.rh, .amy1, .gst, .lpa`.
   * Optionally add to the `SVCALLER` `emit:` block for downstream reuse.
   * Thread `ch_trait_regions` from `main.nf` through the `SVCALLER` `take:` list.

No change needed to OmniGen here — it independently reads the four published
contract files from `results/<S>/…`.

---

## 4. Testing strategy — in-session vs full run

### 4.1 Testable IN-SESSION now (no pipeline run, no BAM required)

Add `tests/test_cnv_traits.py` (pytest, same style as `tests/test_cnv_consensus.py`):

* **Consensus cross-check logic** against the REAL
  `results/HG002/HG002.cnv_consensus.bed`:
  - assert RHD window → 0 overlapping DELs → `Rh_status = pos` (HG002 is Rh D positive ✔)
  - assert GSTM1/GSTT1 windows → 0 consensus DELs (documents that consensus alone
    misses HG002's known GSTM1-null → justifies depth-primary design)
  - assert LPA window → the 1 known DUP overlap is found.
* **Interpreter arithmetic** with synthetic mosdepth `regions.bed.gz` fixtures:
  feed known depths → assert `RHD_copies`, `AMY1_copies`, GST null thresholds,
  `KIV2_copies` round correctly (boundary tests at ratio 0.15, 0.5).
* **Contract schema round-trip**: each writer emits the exact header + column
  order/types specified above (guards the OmniGen contract).
* **HTML rendering**: extend `tests/test_html_report.py` — pass synthetic trait
  TSVs and assert the new card renders; assert graceful "not available" on
  `NO_FILE`.

All of the above run with `pytest` in-session against committed fixtures + the
existing HG002 bed. **No BAM, no Nextflow, no containers required.**

### 4.2 Prototype-able offline (optional, uses existing HG002 work BAM)

HG002 BAMs already exist under
`work_HG002_wedge2/*/HG002.bwa.sortdup.bqsr.bam`. A standalone
`mosdepth --by assets/cnv_trait_regions.bed` on that BAM can sanity-check the
depth numbers and the normalization baseline before any pipeline run. This
exercises `TRAIT_DEPTH` logic without running the pipeline. Use known truth for
regression: HG002/NA24385 is **Rh D positive**; document expected AMY1 / GSTT1 /
GSTM1 / LPA-KIV2 ranges as the interpreters are calibrated. (Marked optional
because it runs a tool on a BAM; strictly speaking it is not part of unit tests.)

### 4.3 Needs a FULL pipeline run

* End-to-end wiring (`main.nf → SVCALLER → CNV_TRAITS → REPORT`), container
  resolution, `publishDir` paths, and the actual per-sample depth-derived calls
  landing at the four contract paths and inside `results/<S>/<S>.report.html`.
* Cross-sample validation of thresholds (RHD deletion cutoff, GST null cutoff,
  AMY1/KIV-2 calibration) against additional GIAB samples.

---

## 5. Implementation checklist (order)

1. `assets/cnv_trait_regions.bed` (verify coords) + `main.nf` `ch_trait_regions`.
2. `bin/cnv_traits_common.py` (normalization helpers).
3. `bin/rh_status.py`, `bin/amy1_cn.py`, `bin/gst_null.py`, `bin/lpa_kiv2.py`
   (each with a `--consensus-only` path for unit tests without a depth bed).
4. `modules/traits/depth.nf` (`TRAIT_DEPTH`) + container in
   `conf/docker.config` / `conf/singularity.config`.
5. `subworkflows/cnv_traits.nf` (`CNV_TRAITS` + 4 interpreter processes).
6. Wire `CNV_TRAITS` into `workflows/svcaller.nf`.
7. `bin/html_report.py` new section + args; `subworkflows/report.nf` new
   channels/inputs/args; `svcaller.nf` `REPORT(...)` extension.
8. `tests/test_cnv_traits.py` + extend `tests/test_html_report.py` (in-session).
9. Bump `VERSION.md` / pipeline manifest version.
10. Full run for end-to-end validation (§4.3).

Non-goals: modifying the existing SV/CNV ensemble; changing the cnv_consensus
schema; any OmniGen-side code.

---

## 6. Running for real depth numbers (needs a full pipeline run)

The wiring, containers, contract paths, and the report card are all in place and
unit-tested, but the actual RHD/AMY1/GST/KIV-2 numbers require `TRAIT_DEPTH`
(mosdepth) to run on a real BAM. The work-dir BAMs under
`work_HG002_wedge2/*/HG002.*.bam` are 0-byte/cleaned (auto_cleanup ran), so an
offline prototype was not possible — run the pipeline on a real BAM.

Run the whole pipeline (traits are wired into `SVCALLER`, no extra flag needed):

```bash
cd /data/alvin/SVcaller
nextflow run main.nf -profile docker \
    --input   validation/giab_samplesheet.csv \
    --ref_fasta /path/to/hg38.canonical.fa \
    --outdir  results
# → results/<S>/bloodgroup/<S>.rh_status.tsv
#   results/<S>/cnv_traits/<S>.{amy1,gst_null,lpa_kiv2}.tsv
#   the "Blood Group & Copy-Number Traits" card in results/<S>/<S>.report.html
```

To prototype ONLY the targeted depth on an existing BAM (fast, no full run) once a
real BAM is available:

```bash
mosdepth --no-per-base --mapq 0 --by assets/cnv_trait_regions.bed \
    HG002.trait_depth /path/to/HG002.bam
bin/rh_status.py --depth HG002.trait_depth.regions.bed.gz \
    --cnv-bed results/HG002/HG002.cnv_consensus.bed --sample HG002 \
    --out HG002.rh_status.tsv
# repeat for amy1_cn.py / gst_null.py / lpa_kiv2.py
```

### Calibration TODO after first real run
* Confirm `RHD < 0.5` deletion cutoff and `GST < 0.15` null cutoff on a known
  sample (HG002 is Rh D positive; use it as a positive control).
* Calibrate the LPA KIV-2 absolute scaling against a truth sample — the estimate
  is proportional to copy number but the constant depends on how many KIV-2 units
  the reference window spans.
* Re-confirm every window in `assets/cnv_trait_regions.bed` against RefSeq/GENCODE
  for the exact assembly, and pick `CTRL_*` regions verified copy-stable/non-segdup.

A rendered demo of the report card (synthetic fixture values, not real depth) is
at `docs/cnv_traits_card_demo.html`.
