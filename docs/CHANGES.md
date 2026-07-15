# Changes

## 2026-07-15 — Fail-loud guards on empty report inputs + drop dead SURVIVOR_MERGE selector

**Silent-failure class (the OmniGen "clean bill of health" bug), now closed in the reporters.**
A sibling consumer (OmniGen) had a bug where a crashed upstream caller left a 0-byte
placeholder that, being gated only on `os.path.exists()` (True for an empty file), rendered
as a clean *negative* result instead of a failure. Two SVcaller reporters had the same latent
shape:

- `bin/smn_report.py::parse_smn_tsv` opened the SMN TSV with no non-empty check. A 0-byte or
  header-only file (a crashed `SMNCopyNumberCaller`) fell through to `if not lines: return`
  defaults, and `render_html_section` then defaulted `smn1/smn2` to `2` — rendering **"Normal
  (CN=2)"** for a caller that never ran. Now `parse_smn_tsv` raises `SmnInputError` on a
  present-but-empty / header-only-no-data TSV; a `NO_*` sentinel (a legitimate skip) still
  returns the neutral/unknown result.
- `bin/html_report.py::render_report` read its required inputs (`--smn-html`, `--circos-svg`,
  `--sv-tsv`, `--cnv-bed`) directly (`Path(...).read_text()` / `parse_*`), so a 0-byte input
  rendered an empty section / "no findings". Added `UpstreamEmptyError` + `_is_absent` +
  `_require_nonempty`, called at the top of `render_report`. A present-but-empty required
  input now fails loudly.

The distinction mirrors OmniGen `prototype/upstream.py`: **absent != empty != populated.**
ABSENT (a Nextflow `NO_*` sentinel, e.g. `--cnv-bed NO_FILE` for a sample with no CNV data) is
a legitimate skip and passes silently. EMPTY (0 bytes / whitespace — a crash placeholder) is a
hard, visible failure. A **header-only** file ("we looked, found nothing") remains a valid
negative result and is deliberately allowed through — the existing empty-AnnotSV-TSV → merged-VCF
fallback still works. Healthy-path output is unchanged: verified byte-for-byte-shape against the
real `results/HG002/` inputs (report renders at 1.28 MB, matching the committed 1.29 MB), and the
guards do not fire on non-empty inputs.

New regression tests (`tests/test_smn_report.py`, `tests/test_html_report.py`) assert a 0-byte /
header-only required input raises rather than producing blank output, and that a `NO_FILE`
sentinel is still an honest skip. Suite: 67 → **77 passed**.

**Dead config selector removed.** `conf/docker.config` carried
`withName: 'SURVIVOR_MERGE' { container = '...survivor...' }`, but the active pipeline merges SVs
with `JASMINE_MERGE` (`subworkflows/sv_calling.nf`). `SURVIVOR_MERGE` is referenced only by the
separate `workflows/sv_pon_build.nf` panel-of-normals builder, and `modules/survivor/merge.nf`
already declares its container inline — so the selector was dead weight that made Nextflow emit
`WARN: There's no process matching config selector: SURVIVOR_MERGE` on every main-pipeline run.
Removed. `modules/survivor/merge.nf` is **not** orphaned (still used by `sv_pon_build.nf`) and was
left in place.

## 2026-07-13 — CNV_TRAITS scripts committed non-executable (exit 126)

`bin/gst_null.py`, `bin/amy1_cn.py`, `bin/lpa_kiv2.py` and `bin/rh_status.py` are invoked
as bare commands on `PATH` by `subworkflows/cnv_traits.nf`, but were committed with mode
`100644`. On any fresh checkout the CNV_TRAITS processes failed immediately with
`.command.sh: Permission denied` and exit status 126 (observed on a GIAB HG003 run).
Earlier HG001/HG002 runs only worked because the executable bit had been set locally in
that working tree; it was never carried in the git index.

The executable bit is now set in the git index (mode `100755`) for those four scripts.
All four already carried a `#!/usr/bin/env python3` shebang; no script logic changed.
`bin/cnv_traits_common.py` stays `100644` (imported as a Python module, never executed),
as does `bin/nf-cleanup.sh` (invoked as `bash bin/nf-cleanup.sh`).

## 2026-07-13 — Fail loudly instead of publishing empty placeholder outputs

### Why

A production incident traced back to this repo.

SVcaller processes used to end their command blocks with an `|| touch <output>`
fallback (and the JSON equivalent, `echo '{}' > <output>`). When a caller **failed**,
the process still exited 0 and **published a zero-byte output file**. Nextflow saw a
satisfied output declaration, marked the task successful, and the run continued.

Downstream, OmniGen (a consumer genomics report) gated its reads on
`os.path.exists()` — which returns `True` for a 0-byte file. An empty `smn.tsv`
crashed its scan, and the report renderer then produced a **complete-looking consumer
report** whose summary tiles read *"0 Carrier findings, 0 Medication flags, Clear"*.

**A crashed SMN caller was rendered as a clean bill of health for a human being.**

OmniGen has since been hardened to refuse empty artifacts, but SVcaller was the thing
*generating* them. Any consumer — current or future — stayed exposed, and a failed
stage still looked like a completed run on disk. This change fixes the source.

### The distinction that matters

Not every empty output is a bug, and this change is careful not to break working runs:

| Case | Example | Treatment |
|---|---|---|
| **Legitimate empty result** | A VCF with a proper header and zero variant rows — the caller ran, looked, and genuinely found nothing. | **Still allowed.** Preserved as-is. |
| **Explicitly skipped stage** | `--skip_melt`, `--skip_scramble`, `--skip_svaba`, `--skip_gridss` route to a dedicated `*_STUB` process. | **Still allowed.** The skip is a deliberate, recorded choice. |
| **Masked failure** | The caller crashed, or was never installed, and a placeholder file was written so the process could exit 0. | **Now a hard failure.** No placeholder is written. |

The rule: *a crashed caller must never be indistinguishable from a negative result.*

### Sites fixed

| File:line (original) | Original line | Change |
|---|---|---|
| `modules/smn_caller/call.nf:34` | `touch ${meta.id}.smn.tsv` | Removed. Fails with exit 1 if no result table was produced, or if the table has no sample row (header only). The originating site of the incident. |
| `modules/smn_caller/call.nf:38` | `echo '{}' > ${meta.id}.smn_detail.json` | Removed. `smn_caller.py` emits `.tsv` and `.json` together; a missing JSON means an incomplete run. |
| `modules/annotsv/annotate.nf:36` | `[ -n "$f" ] && mv "$f" . \|\| touch ${meta.id}.annotated.tsv` | Now distinguishes the two cases: if the **input VCF has zero SV records**, a header-only TSV is emitted (a legitimate empty annotation). If the input **had** SVs and AnnotSV still produced nothing, the annotation failed silently → exit 1. |
| `modules/expansionhunter/call.nf:38` | `echo '{}' > ${meta.id}.str_profile.json` | Removed → exit 1. This file is published to `results/` and read directly by OmniGen; a `{}` profile reads as "no repeat expansions detected". |
| `modules/svaba/call.nf:23` | `svaba run ... 2>&1 \|\| true` | `\|\| true` removed. A SvABA crash used to fall through to the empty-VCF stub below it and render as "0 SVs found". Also asserts SvABA actually wrote its SV VCF, and that the normalised VCF carries a `#CHROM` header. |
| `modules/scramble/call.nf:26` | `scramble.sh ... --eval-meis \|\| true` | `\|\| true` removed. A SCRAMble crash used to fall through to the header-only-VCF branch and render as "no MEIs detected". The header-only branch is **kept** for the genuine zero-MEI case. |
| `modules/melt/call.nf:25-27` | `echo "WARNING: MELT.jar not found; emitting empty VCF"` + `exit 0` | A missing MELT install is a **misconfiguration**, not a result. Now exit 1, pointing the user at `--skip_melt`. |
| `modules/melt/call.nf:40-42` | `echo "WARNING: no *_MELT.zip found ...; emitting empty VCF"` + `exit 0` | Same: missing ME references → exit 1, not a confident empty VCF. |
| `modules/melt/call.nf:61` | `java ... Single ... 2>&1 \|\| true` | Per-ME-type exit codes are now collected; if any ME type crashed, the process fails rather than reporting zero insertions. |
| `modules/melt/call.nf:71-73` | no `*.final_comp.vcf` → empty VCF + `exit 0` | If every MELT type exited 0 but nothing was written, MELT did not really run → exit 1. |

Also fixed while in `melt/call.nf`: when `--melt_refs` was supplied, `MELT_DIR` was never
set, so the `-n` gene-annotation BED path (`${MELT_DIR}/add_bed_files/...`) silently
resolved to a bogus location. `MELT_DIR` is now set on both branches.

### Deliberately left alone (legitimate, not failure-masking)

- `modules/sv_pon/annotate.nf:21` — `... | cut -f3 > pon_hit_ids.txt || true`. `grep` exits 1
  when there are **no matches**, which is the normal "this sample hits no PON sites" case.
  An empty `pon_hit_ids.txt` here is a real result.
- `modules/truvari/bench.nf:38-44` — four `|| true` on Truvari size-bin benchmarks. Truvari
  legitimately exits non-zero when a size bin contains zero calls. This is the
  **validation/benchmark** path, not the clinical report path.
- `modules/jasmine/merge.nf:93,95` — `grep -cv '^#' ... || echo 0`. Counting variants for a log
  line; `grep -c` exits 1 on zero matches.
- `modules/melt/call.nf:19` (`bowtie2 --version || true`) and `modules/scramble/call.nf:41`
  (version-string fallback) — cosmetic version capture for `versions.yml`.
- `modules/svaba/call.nf` cleanup `rm -f ... 2>/dev/null || true` — intermediate file cleanup.
- `subworkflows/report.nf:130,244` — `.ifEmpty([])` on MultiQC inputs. MultiQC is genuinely
  optional and an absent QC report does not misrepresent a clinical finding.
- The `*_STUB` processes (`melt/stub.nf`, `scramble/stub.nf`, `gridss/stub.nf`,
  `SVABA_STUB`) — these emit header-only VCFs *by design*, and are only reachable via an
  explicit `--skip_*` flag. That is the sanctioned way to express "this stage is optional
  for this sample": a recorded decision, not a fake empty file.

### Blast radius

`conf/base.config:6` is the only `errorStrategy` in the repo:

```groovy
errorStrategy = { task.exitStatus in [143,137,104,134,139] ? 'retry' : 'finish' }
```

There is no `errorStrategy 'ignore'` anywhere. The exit codes above are OOM/kill signals;
a deliberate `exit 1` is **not** in that list, so a newly-failing process will **not** be
retried in a loop. It gets `'finish'`: Nextflow stops launching new tasks, lets in-flight
tasks drain, and then reports the failing process with its work dir, exit status, and the
`.command.err` we now write. That is an actionable error, not a confusing partial run.

**Stages that now hard-fail where they previously passed silently:**

1. `SMN_CALLER` — when `smn_caller.py` produces no result table, or a header-only table.
2. `EXPANSIONHUNTER` — when EH exits 0 without writing its STR profile JSON.
3. `ANNOTSV` — when the input VCF contains SV records but AnnotSV emits no annotation.
   (Zero SVs in → header-only TSV out; still passes.)
4. `SVABA_CALL` — when SvABA itself exits non-zero, or writes no SV VCF.
5. `SCRAMBLE_CALL` — when `scramble.sh` exits non-zero.
6. `MELT_CALL` — when MELT is not installed, ME references are missing, any ME type crashes,
   or no `final_comp.vcf` is produced. Use `--skip_melt` to run without MELT on purpose.

These are all cases that previously produced a **successful-looking run with a fabricated
empty result**. If any of them starts failing on a pipeline that "used to work", that
pipeline was not working — it was silently reporting nothing found.

### Tests

`tests/test_no_empty_placeholders.py` (9 tests) statically asserts the masking patterns stay
dead: no `|| touch` / bare `touch <output>` fallback in any `.nf` file, no `echo '{}' >`
placeholder, and no `|| true` swallowing the exit code of the SvABA or SCRAMble invocations.
Each guard was mutation-tested (bug reintroduced → test fails; reverted → test passes).

Full suite: **67 passed** (58 pre-existing + 9 new).

Not verified by these tests: end-to-end Nextflow behaviour. No pipeline run was performed
for this change (code, docs and unit tests only). The shell blocks of every edited module
were syntax-checked with `bash -n` after rendering the Nextflow interpolation.

### Follow-up fix — Groovy escape bug in `smn_caller/call.nf`

The fail-loud edit introduced a comment in the `script:` block of
`modules/smn_caller/call.nf` (line 31) that contained backslash-escaped backticks. Inside a
Groovy triple-quoted GString those are invalid escapes, so **every** `nextflow run` aborted
before any process started:

```
ERROR ~ Module compilation error
- file : modules/smn_caller/call.nf
- cause: token recognition error at: '`touch\`' @ line 31, column 27.
```

The comment was reworded to use single quotes instead of backticks. No pipeline logic
changed: the `|| touch` fallback stays removed and `errorStrategy` is untouched. Verified
with `nextflow run main.nf --help`, which now compiles all modules and reaches parameter
validation (`ERROR: --input is required`) instead of failing at module compilation.
