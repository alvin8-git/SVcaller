# SVcaller directory organization + rules

Supersedes the root `ORGANIZATION.md` (which covered only the 2026-07-12
`results/` consolidation). Written 2026-07-21 against the tree as it stands after
the git cleanup and SMA BAM relocation.

---

## Part 1 — What is actually in here

Every directory in `/data/alvin/SVcaller`, grouped by lifecycle. Lifecycle is the
useful axis for a Nextflow project, because it determines the three things you
need to know about any path: is it versioned, is it safe to delete, and who
regenerates it.

### A. Pipeline source — tracked, reviewed, small

The pipeline itself. Nothing here is machine-generated.

| Path | Size | Contents |
|---|---|---|
| `main.nf` | 5 KB | Entry point: samplesheet parsing, channel setup, `SVCALLER` call |
| `nextflow.config` | 3 KB | `params` block, profiles, `check_max()` |
| `conf/` | 20 KB | `base` (resource tiers, OOM retry), `docker`, `local`, `test` |
| `workflows/` | 16 KB | 3 entry workflows: `svcaller`, `pon_build`, `sv_pon_build` |
| `subworkflows/` | 48 KB | 7 module-level stages (M1–M6) |
| `modules/` | 272 KB | 24 tool wrappers, one dir per tool, `<tool>/<verb>.nf` |
| `bin/` | 412 KB | 12 Python CLIs + `nf-cleanup.sh`; auto-staged onto `PATH` |
| `assets/` | 344 KB | Static inputs baked into the image: report template, EH catalog, cytobands, trait BEDs, syndrome/STR tables |

### B. Tests — tracked

| Path | Size | Contents |
|---|---|---|
| `tests/` | 228 KB | 9 pytest modules, 79 tests, ~0.5 s, no containers needed |

### C. Run definitions — tracked, small, but they encode absolute paths

| Path | Size | Contents |
|---|---|---|
| `validation/` | 84 KB | 13 samplesheets, `smn_truth_table.tsv`, 4 driver/download shell scripts |

### D. Documentation — tracked

| Path | Size | Contents |
|---|---|---|
| `docs/` | 5.5 MB | Diátaxis-ish set (tutorial / howto ×5 / reference ×2 / explanation), `CHANGES.md`, `superpowers/{specs,plans}`, `demo/` (3 committed HTML reports), `img/` |
| `README.md`, `CLAUDE.md` | 53 KB | Entry doc; agent instructions |
| `Documentation.md`, `ORGANIZATION.md`, `TODO.md`, `VERSION.md` | 74 KB | Root-level sprawl — see Rule 9 |

### E. Container definitions — tracked (except the blob)

| Path | Size | Contents |
|---|---|---|
| `Dockerfile.{utils,melt,smncopynum}` | 4 KB | The three locally-built, non-pullable images |
| `environment.yml` | 645 B | Conda env (secondary to Docker) |
| `MELTv2.2.2.tar.gz` | 90 MB | Registration-walled MELT installer; gitignored, required for the image build |

### F. External input data — **does not belong in the repo**

| Path | Size | Status |
|---|---|---|
| ~~`benchmarks_data/`~~ | 260 GB | **Moved 2026-07-21** → `/data/alvin/ref/benchmarks/<sample>/fastq/` |
| `Benchmarks/` | 168 MB | Still here — GIAB truth BEDs/VCFs (HG002 v5.0q, HG003/HG004 v4.2.1) |

Downloaded external artifacts, not products of this pipeline. The COLO829 and
HCC1395 reads now sit beside the truth sets they are benchmarked against, under
one `/data/alvin/ref/benchmarks/` tree. `Benchmarks/` is the same category and
should follow, into `/data/alvin/ref/GIAB/truth/` — its near-identical name to
the old `benchmarks_data/` was its own hazard. See Rule 2.

### G. Built panels — regenerable, expensive

| Path | Size | Contents |
|---|---|---|
| `pon/` | 525 MB | `pon/pon/giab_cnv_pon.hdf5` + interval list; `pon/sv_pon/giab_sv_pon.bed` |

Products of `workflows/pon_build.nf` and `sv_pon_build.nf`. Rebuilt only if the
GIAB BAMs change. The `pon/pon/` double-nesting is an artifact of
`--outdir pon` (see Rule 6).

### H. Pipeline outputs — regenerable, the deliverable

| Path | Size | Contents |
|---|---|---|
| `results/` | 356 GB | 9 sample dirs (~40–103 MB each: VCF, BED, circos, `report.html`, xlsx), `SMN/` trio, `_experiments/`, plus the two entries below that are *not* outputs |
| `results/cache`, `results/.cache` | 11 GB | `storeDir` caches: `gridss_ref`, `gatk_preprocess`, `filter_chroms` |
| `results/pon` | 79 MB | A *copy* of the panel that also lives in `pon/` |

Actual per-sample outputs total under 1 GB. The 356 GB is dominated by
`_experiments/` and the caches.

### I. Scratch — always disposable

| Path | Size | Contents |
|---|---|---|
| `work/` | 1.1 GB | Nextflow task dirs — current shared scratch |
| ~~`work_HG001`, `work_HG004…HG007`~~ | 20 GB | **Deleted 2026-07-21** once each report was published |
| `.nextflow/`, `.nextflow.log{,.1–.9}` | ~600 KB | Session state + 10 rotated logs |
| `pytest-of-alvin/` | 2.6 MB | pytest `tmp_path` factory output — see Rule 4 |

---

## Part 2 — Organization rules

### Rule 1 — Four lifecycles, and every path belongs to exactly one

`source` (tracked, reviewed) · `input` (external, read-only, central) ·
`output` (regenerable, gitignored) · `scratch` (deletable at any moment).

If you cannot say which one a new directory is, it does not have a home yet.
Never let one directory hold two lifecycles — that is precisely why `results/`
is 356 GB: it mixes outputs, an 11 GB regenerable cache, and a duplicate panel.

### Rule 2 — Input data lives in `/data/alvin/ref`, never in the pipeline repo

The repo holds *code that reads data*, not the data. Reference genomes, truth
sets, PON source BAMs, and benchmark reads all belong to the shared store, where
they are reused by other pipelines and backed up once:

```
/data/alvin/ref/GRCh38/      reference FASTA + indices
/data/alvin/ref/GIAB/        GIAB BAMs + truth sets
/data/alvin/ref/SMA/         SMA trio BAMs        (moved 2026-07-21)
/data/alvin/ref/benchmarks/  COLO829, HCC1395 — truth sets AND reads
/data/alvin/ref/annotsv/     AnnotSV annotation DB
```

**Outstanding:** `benchmarks_data/` (260 GB) and `Benchmarks/` (168 MB) still
violate this. Both are same-filesystem renames, so the move is instant:

```bash
mv benchmarks_data/COLO829/* /data/alvin/ref/benchmarks/COLO829/
mv benchmarks_data/HCC1395/* /data/alvin/ref/benchmarks/HCC1395/
mkdir -p /data/alvin/ref/GIAB/truth && mv Benchmarks/* /data/alvin/ref/GIAB/truth/
```

Test for this rule: *if I deleted this directory, would I re-download it or
re-compute it?* Re-download ⇒ it is input ⇒ it belongs in `ref`.

### Rule 3 — Absolute paths appear only in samplesheets and `--` arguments

Never hardcode a data path in `.nf`, `.config`, or `bin/`. Samplesheets under
`validation/` are the one place absolute paths are allowed, and even there prefer
a `${REF}`-style variable so a sheet stays runnable on another machine. Moving
the SMA BAMs required editing 4 samplesheets and 1 doc; had they been
`${REF}`-relative it would have required editing nothing.

### Rule 4 — the environment must supply `TMPDIR`; scratch never lands on `/` or in the repo

```
nextflow.config:10    tmp_dir = System.getenv('TMPDIR') ?: '/tmp'
conf/docker.config:2  docker.runOptions = "--rm -v ${params.tmp_dir}:/tmp"
```

`TMPDIR` is **unset** in a login or interactive shell here, and is set nowhere in
`~/.bashrc`, `~/.profile`, `/etc/environment`, or `/etc/profile.d/`. So a real run
takes the `?: '/tmp'` fallback and bind-mounts `/tmp` — which lives on
`/dev/mapper/vgubuntu-root`, the 915 GB **OS partition**, not on `/data`. That
contradicts this machine's standing convention:

> Always use `/data/alvin/tmp` instead of `/tmp` … `/tmp` is on a small root partition.

Commit `f4a5ac0` did the right thing by removing the hardcoded `/data/alvin/tmp`
from `conf/docker.config` — the path is host-specific and does not belong in a
portable pipeline. It moved the responsibility to the environment. The
environment then never took it up. **Fix the environment, not the config:**

```bash
echo 'export TMPDIR=/data/alvin/tmp' >> ~/.bashrc     # persistent
nextflow run … --tmp_dir /data/alvin/tmp              # per-run override
```

Two guards worth adding: refuse to start if `params.tmp_dir` resolves inside
`projectDir` (a container would then bind-mount the repo over `/tmp`), and warn
if it resolves onto the root filesystem.

*Caveat on `pytest-of-alvin/`:* it is real and is now gitignored, but it was
produced by an agent sandbox that exports `TMPDIR=$PWD`, not by normal use. It is
evidence that an unset `TMPDIR` lets **whatever the caller happens to set**
decide where scratch lands — which is the argument for pinning it explicitly.

### Rule 5 — One work directory per sample or batch, deleted once results land

`work_<sampleId>` for single samples, `work_<batchName>` for batches. This is
deliberate (see `docs/explanation-design.md`) — a shared `work/` causes session
lock conflicts and makes targeted cleanup impossible. The corollary is the half
people skip: **delete it as soon as `results/<sample>/<sample>.report.html`
exists.** Do not consolidate stale work dirs into `work/`; remove them.

Outstanding: `work_HG001` + `work_HG004…HG007` ≈ 20 GB, all five published.

### Rule 6 — `--outdir` is a publish target, and holds only published outputs

`results/` should contain sample directories and nothing else. Today it also
holds an 11 GB `storeDir` cache and a duplicate panel. Give each its own root:

```
--outdir      results/          per-sample deliverables only
storeDir      cache/            regenerable, reference-derived, shared across samples
pon/                            built panels, one canonical copy
```

Never let `storeDir` live under `--outdir`: the cache is keyed to the *reference*,
not the sample, so it survives cleanup that outputs should not — and burying it in
`results/` makes "how big is my output?" unanswerable. Also avoid
`--outdir pon` when the workflow publishes into `pon/`; that is what produced
`pon/pon/`.

### Rule 7 — Module layout is `modules/<tool>/<verb>.nf`, one process per file

Already followed by all 24 tool dirs. Each file: `input` tuple → `script` →
`output` with named `emit` + a `versions.yml`. A process gets a resource *label*
(`process_low`…`process_gridss`), never inline `cpus`/`memory` — tiers live in
`conf/base.config` so they can be retuned in one place.

### Rule 8 — Gitignore by lifecycle, and state the reason

Ignore whole categories (outputs, scratch, input data, binaries), not individual
files, so a new sample or run never needs a gitignore edit. Keep the comment
explaining *why* each block exists and where the data actually lives — the
current `.gitignore` does this well and it is why the 884 GB of packed BAM blobs
never became a tracked-file problem.

The single intentional exception is `docs/demo/**`, negated back in on purpose.

### Rule 9 — Documentation lives in `docs/`; the root holds `README` + `CLAUDE` only

`docs/` already follows a clean tutorial / how-to / reference / explanation split.
The root does not: `Documentation.md`, `ORGANIZATION.md`, `TODO.md`, `VERSION.md`
should move to `docs/` (`VERSION.md` → `docs/CHANGELOG.md`), and this file
replaces `ORGANIZATION.md`.

### Rule 10 — Ship the files a public repo is expected to have

Missing: `LICENSE` (while `README.md` renders an MIT badge linking to it — the
link 404s and the project is therefore unlicensed), `CONTRIBUTING.md`,
`CITATION.cff`, and `.github/workflows/` (the 0.5 s test suite is a three-line CI
job). `LICENSE` is the urgent one: it is a claim currently contradicted by the
tree.

---

## Outstanding violations, by payoff

| # | Rule | Action | Reclaims |
|---|---|---|---|
| 1 | 4 | `export TMPDIR=/data/alvin/tmp` in `~/.bashrc` — currently unset, so runs bind-mount the OS partition's `/tmp` | correctness |
| 2 | 2 | `Benchmarks/` → `/data/alvin/ref/GIAB/truth/` | 168 MB |
| 3 | 6 | Move `storeDir` cache + `results/pon` out of `results/` | clarity (11 GB) |
| 4 | 10 | Add `LICENSE` | correctness |
| 5 | 9 | Move 4 root `.md` files into `docs/` | tidiness |

## Completed 2026-07-21

- **Rule 2** — SMA trio BAMs (397 GB) → `/data/alvin/ref/SMA/`; 4 samplesheets +
  `Documentation.md` repointed.
- **Rule 2** — `benchmarks_data/` (260 GB) → `/data/alvin/ref/benchmarks/<sample>/fastq/`,
  reuniting the COLO829/HCC1395 reads with the truth sets already stored there.
  `urls.txt` → `download_urls.txt` beside the data, its 8 `dir=` lines repointed;
  `validation/colo829_samplesheet.csv` updated.
- **Rule 5** — `work_HG001`, `work_HG004…HG007` (20 GB) deleted after confirming
  each sample had a published `report.html` and no session lock was held.
- **Rule 8** — 884 GB of abandoned git pack debris removed (see
  `reorg-plan-2026-07-21.md`); `pytest-of-alvin/` gitignored.

Repo went from 637 GB to ~357 GB, of which 356 GB is `results/`. The tracked
source tree is about 6 MB.
