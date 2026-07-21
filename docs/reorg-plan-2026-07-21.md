# Repo reorganization + SvABA scope decision — 2026-07-21

**Status: §1 (git cleanup) and §2 (SMA BAM move) executed 2026-07-21. §3–§6 are
still open recommendations.**

Result: `.git` 884 GB → 2.7 MB, `/data` free space 1.3 TB → 2.1 TB, SMA BAMs
relocated to the central store. All 174 commits, both rescue branches, and
`origin/main` parity intact; `git fsck` clean; 79/79 tests pass. Pre-change
backup bundle (all refs, 2.5 MB) at
`/data/alvin/tmp/SVcaller-prehousekeep-2026-07-21.bundle`.

---

## 1. `.git` is 884 GB of abandoned repack debris — reclaim ~880 GB

`git count-objects -v -H` on 2026-07-21:

```
count: 1427          size: 500.59 MiB     # the real history, loose
in-pack: 3           size-pack: 396.25 GiB
garbage: 3           size-garbage: 486.93 GiB
```

`.git/objects/pack/` holds:

| File | Size | What it is |
|---|---|---|
| `pack-e5ef4758….pack` | 425 GB | `.idx` is **1172 bytes** and indexes **3 objects** |
| `tmp_pack_AzWxH1` | 1.5 GB | aborted repack, Jun 5 10:52 |
| `tmp_pack_dHcW8g` | 107 GB | aborted repack, Jun 5 12:59 |
| `tmp_pack_kbMcEt` | 414 GB | aborted repack, Jun 5 17:14 |

Three successive failed repacks on Jun 5, almost certainly a `git add` of
`results/` / `ValidationBAM/` before `.gitignore` covered them. Git itself
classifies the `tmp_pack_*` files as *garbage*.

**The genuine repository is a few MB.** 170 commits, 1345 blobs, and the largest
reachable blob in all of history is `docs/demo/COLO829_report.html` at 2.2 MB.
`main` is level with `origin/main`; nothing is unpushed.

Reachability was verified independently. The 3 objects in the orphan pack are
blobs of 132 GB, 134 GB and 130 GB — the three SMA BAMs, from the Jun 5 `git add`.
All three are dangling; `git fsck --connectivity-only` reports zero missing or
broken objects, so nothing reachable depends on any of it.

### ⚠️ Two orphaned stashes were found — already rescued

`git stash list` reported nothing and `refs/stash` did not exist, but two real
stash commits survived as dangling objects, absent from all 184 reflog entries.
`git gc --prune=now` would have destroyed both permanently:

| Rescue branch | Commit | WIP on |
|---|---|---|
| `rescue/stash-2026-06-07` | `e8e7d21` | `9ab6757` fix(report): caller list + gnomAD text for P5-P7 |
| `rescue/stash-2026-06-04` | `f97d0ee` | `ffe24a6` fix(P2): filter single-caller DUP/INV in Jasmine post-merge awk |

Both touch `bin/html_report.py`, `assets/report_template.html` and
`subworkflows/report.nf`, and differ substantially from HEAD — this is not
already-merged work. **Branches were created on 2026-07-21, so both are now
reachable and gc-safe.** Review them and either merge or delete the branches
before pruning; do not drop them silently.

**Action (safe to run now — no dangling commits remain):**

```bash
cd /data/alvin/SVcaller
rm -f .git/objects/pack/tmp_pack_*          # 487 GB, git-declared garbage
git gc --prune=now                          # drops the 3-object orphan pack
```

Skip `--aggressive`: it buys nothing on a few-MB history and forces a full
re-delta. Confirm no Nextflow run is active first.

Belt-and-braces alternative, given the remote holds everything: re-clone into a
fresh directory and move the untracked working data across. Do this only while no
Nextflow run is active (verified idle 2026-07-21).

---

## 2. Move the SMA BAMs to the central `/data/alvin/ref` — 397 GB out of the repo

The repo is inconsistent about where input BAMs live. HG002 already reads from
the central store, but the three SMA BAMs sit *inside the working tree*:

```
validation/validation_samplesheet.csv  →  /data/alvin/ref/GIAB/HG002…bam        ✅
validation/smn_*_samplesheet.csv       →  /data/alvin/SVcaller/ValidationBAM/…  ❌ 397 GB
```

`ValidationBAM/SMA_BAM/` is 3 × ~140 GB BAM + `.bai`, owned by `root`.
`/data/alvin/ref` and `/data/alvin/SVcaller` are the same filesystem (`/dev/sdb`),
so the move is a rename — instant, no copy, no extra space needed.

```bash
sudo mkdir -p /data/alvin/ref/SMA
sudo mv /data/alvin/SVcaller/ValidationBAM/SMA_BAM/* /data/alvin/ref/SMA/
rmdir /data/alvin/SVcaller/ValidationBAM/SMA_BAM /data/alvin/SVcaller/ValidationBAM
```

Then update the four samplesheets that reference the old path
(`smn_SMAD`, `smn_SMAM`, `smn_SMAPB`, `smn_validation`) and drop `ValidationBAM/`
from `.gitignore`. Consider `${REF}`-style indirection so samplesheets stop
hardcoding absolute paths — this is the same portability issue commit `f4a5ac0`
addressed for the tmp dir.

---

## 3. `work_HG00*` — delete rather than nest

Five stale scratch dirs, ~20 GB total. All five samples have published reports in
`results/<sample>/<sample>.report.html`, which is exactly the condition the
project's own convention requires before deletion.

```bash
cd /data/alvin/SVcaller && rm -rf work_HG001 work_HG004 work_HG005 work_HG006 work_HG007
```

**Do not** consolidate them into `work/`. Per-sample work directories are a
deliberate design decision documented in `docs/explanation-design.md`
("Per-sample work directories") — a shared `work/` reintroduces the session-lock
conflicts and all-or-nothing cleanup that layout was chosen to avoid. They are
also already gitignored (`work_*/`), so they are invisible on GitHub. The
untidiness is local disk hygiene, not a repository-structure problem, and the fix
is removal, not relocation.

---

## 4. Missing GitHub-standard files

`README.md` renders a `License: MIT` badge linking to `LICENSE`, and §License
says "MIT — see LICENSE". **No LICENSE file exists**, so both links 404 and the
project is technically unlicensed (default copyright, not MIT).

Missing: `LICENSE`, `CONTRIBUTING.md`, `CITATION.cff`, `.github/` (no CI —
`pytest tests/` runs in 0.5 s and would be a trivial workflow).

Adding `LICENSE` is the only urgent one; it is a stated claim that is currently false.

---

## 5. Root-level Markdown sprawl (cosmetic)

Six `.md` files at root: `README.md`, `CLAUDE.md`, `Documentation.md`,
`ORGANIZATION.md`, `TODO.md`, `VERSION.md`. Convention keeps `README` + `CLAUDE`
at root and moves the rest under `docs/` (`VERSION.md` → `CHANGELOG.md`).
`ORGANIZATION.md` largely duplicates what this file covers and could be folded in.

Also: `pytest-of-alvin/` is untracked scratch at root — delete it and add
`pytest-of-alvin/` to `.gitignore`.

---

## 6. SvABA is *not* part of the minimum wedge — do not fix the merge bug

The design spec (`docs/superpowers/specs/2026-05-09-svcnv-caller-design.md` §5,
M2) defines the SV ensemble as exactly **three** structural callers plus a STR
caller:

> JASMINE merges all three structural callers (Manta + DELLY + GRIDSS) into a
> consensus VCF.

SvABA appears nowhere in the spec's goals, non-goals, or tool tables. It arrived
later via `docs/superpowers/plans/2026-06-03-pipeline-design-considerations.md`
as priority **P4**, justified by a *projection* — "+0.01–0.03 F1, 1 day,
6.5k DEL FNs" — that has never been measured, because:

1. Until 2026-07-15 SvABA crashed at startup on every run (its classic BWA index
   was never staged), and a `2>&1 || true` swallowed the failure. See
   `docs/CHANGES.md` 2026-07-15.
2. It *still* contributes nothing: `subworkflows/sv_calling.nf:135` hands six
   VCFs to `JASMINE_MERGE`, but `modules/jasmine/merge.nf` only decompresses
   `vcfs[0..4]` and builds `vcf_list.txt` from those. `vcfs[5]` (SvABA) is staged
   into the task dir and ignored.

This was already suspected — `results/COLO829/BND_validation.md` line 17 flagged
it as "(verify)" on 2026-06-18 and it sat unchecked for a month. Now confirmed.

Meanwhile `SVABA_CALL` is `process_high` (16 CPUs pinned) at roughly 4 h/sample,
and `--skip_svaba` defaults to `false` — so every run has paid full freight for a
caller that has never produced a single variant.

**Recommendation: flip the default to skip, don't wire it in.**

```groovy
skip_svaba = true   // never validated; see docs/reorg-plan-2026-07-21.md §6
```

One line, reversible, and it stops burning ~4 CPU-hours per sample immediately.
Wire SvABA into `vcf_list.txt` only as a deliberate, benchmarked experiment: run
HG002 with and without, and keep it only if the measured ΔF1 justifies the cost.
Until then it is a 4-hour no-op inside the critical path.

**Correct the false claims either way** — `README.md` (lines 21, 362) and the
CLAUDE.md F1 table advertise a "6-caller" ensemble. Every published benchmark
number, including the run16 baseline, is a 5-caller result.
