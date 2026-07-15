# Results directory organization

_Last reorganized: 2026-07-12._

The pipeline previously wrote each run to its own top-level `results_<tag>/` and
`work_<tag>/` directory. Those have been **consolidated into a single `results/`
tree** (one subdirectory per sample) and a single disposable `work/` scratch dir.
Everything under `results*/` and `work*/` is gitignored — only source and the
committed demo reports under `docs/demo/` are tracked.

## Current layout

```
results/
├── HG002/                         # GIAB reference — SV/CNV/STR/SMN + blood-group/CNV traits
├── HG001/                         # GIAB NA12878 (30X GRCh38) + blood-group/CNV traits
├── COLO829/                       # COLO829_tumor — somatic SV/CNV (depth traits omitted, no BAM)
├── SMN/
│   ├── SMAD/                      # SMA trio
│   ├── SMAM/
│   └── SMAPB/
├── _experiments/
│   ├── HG002_wedge/               # retained experimental HG002 run
│   └── svpon/                     # retained experimental HG002 SV-PON run
├── pon/                           # shared Panel of Normals (single copy; all runs used the same PON)
├── cache/  .cache/                # shared reference-derived pipeline cache (regenerable, gitignored)
```

## What changed on 2026-07-12

- **Consolidated** `results_HG001/HG001`, `results_colo829/COLO829_tumor`,
  `results_smn/{SMAD,SMAM,SMAPB}`, `results_HG002_wedge/HG002`, and
  `results_svpon/HG002` into the single `results/` tree above (moves on the same
  filesystem — instant, no copy). Real outputs were verified present at the
  destination before any source directory was removed.
- **De-duplicated** the Panel of Normals: all six per-run `pon/` copies were
  byte-identical (`giab_cnv_pon.hdf5`, 81,904,419 bytes), so a single shared
  `results/pon/` is kept.
- **Reclaimed ~772 GB** by deleting the disposable `work_*` scratch dirs
  (`work_HG001`, `work_HG002`, `work_HG002_wedge`, `work_HG002_wedge2`, `work_smn`)
  and the now-empty `results_*` shells (which held only regenerable cache + the
  duplicate PON). Going forward the pipeline uses a single `work/` scratch dir.
- **Added** HG001 blood-group/CNV depth traits + `HG001.bam_stats.json` and
  regenerated the HG001 and COLO829 HTML reports (published to `docs/demo/`).

No tracked pipeline logic changed as part of the move; the only code edit was a
backward-compatible `--traits-note` option added to `bin/html_report.py` so a
BAM-less sample (COLO829) can state explicitly that depth traits were not computed.

## What changed on 2026-07-13

Reclaiming `work_HG002` also removed the only copies of the HG002 Truvari
`benchmark.json` and ExpansionHunter STR VCF, so those two report sections had
dropped out of the demo. Both stages were **re-run from the surviving inputs**
(`results/HG002/HG002.sv_merged.vcf.gz` + the 30X GRCh38 HG002 BAM) with the
same commands as `modules/truvari/bench.nf` and `modules/expansionhunter/call.nf`,
and `docs/demo/HG002_report.html` was regenerated with `bin/html_report.py`.
The regenerated demo is a strict superset of the 2026-07-12 one (circos with the
CNV-trait ring, trait card and QC all retained). Their outputs are now kept in
`results/` so they survive the next scratch reclaim:

```
results/HG002/HG002.str.vcf.gz{,.tbi}        # ExpansionHunter 5.0.0, 32-locus catalog
results/HG002/HG002.str_profile.json
results/HG002/benchmark/HG002.{T2T,v5q}.truvari_{summary,sizebin}.json
```

Numbers are recorded in [TODO.md](TODO.md) ("HG002 demo report: Truvari + STR
sections restored — 2026-07-13").
