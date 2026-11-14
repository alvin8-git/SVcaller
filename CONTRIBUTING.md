# Contributing to SVcaller

Thank you for helping. SVcaller is a research and educational pipeline; please keep changes and
wording non-clinical.

## Before you open a pull request

1. Run the unit tests: `pytest tests/` (no containers or network needed; a few tests skip without
   `samtools` on `PATH`).
2. Run the smoke test: `bash tests/smoke/run.sh` (Docker and the locally built
   `svcaller/utils:1.3` image are needed; see the README).
3. If you change `bin/` or `assets/`, rebuild `svcaller/utils` from `Dockerfile.utils`; the pipeline
   runs the scripts inside that image.
4. Record any change in behaviour in `docs/CHANGES.md`: what was wrong, what changed, how it was
   checked.

## Data rules

- Only public GIAB samples (HG001 to HG007) may appear in code, tests, fixtures, screenshots or
  example output.
- Do not commit BAM, CRAM, VCF or reference files outside `tests/fixtures/`, and keep fixtures small.
- No hostnames, IP addresses, personal paths or credentials in any committed file.

## Style

- One process per module under `modules/<tool>/`, with a named `emit`.
- A failed step must exit non-zero; never publish an empty placeholder output
  (`tests/test_no_empty_placeholders.py` checks this).
- Documentation figures are Mermaid sources in `docs/img/*.mmd`, rendered with
  `docs/img/render.sh`.

## Reporting problems

Use the issue templates. For security problems, do not open a public issue; contact the maintainer
through GitHub instead.
