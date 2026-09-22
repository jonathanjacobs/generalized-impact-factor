# Changelog

Material changes to this project, and their verification limits. Use Git history for detailed implementation changes; this file is for changes that matter to people developing, reviewing, running, or maintaining the project.

## Unreleased

- Restructured `impact_lookup.py` into a reusable library (`src/impact_factor/`) plus a thin CLI, so the pipeline functions can be imported and reused without running the CLI (`docs/REQUIREMENTS.md` R4). No change to the CLI's argument names or defaults.
- Fixed a bug where the reported "Elapsed time" was always near zero (`start_clock` was never assigned; `start_time` was used instead).
- Fixed `dump_term_citations_tsv`'s `request_delay` parameter, which was previously accepted but unused (the function used a module-level variable instead); it now actually controls the pacing between citation-fetch batches.
- Added a pytest suite (`tests/`) covering `compute_h_g_i10`, `generate_fuzzy_queries`, `smart_tc`, `chunks`, and `read_terms_from_file`.
- Added `environments/requirements.txt`, `.gitignore`, `.editorconfig`, `docs/DESIGN.md`, `docs/REQUIREMENTS.md`, `docs/TESTING.md`, `docs/DATA_GOVERNANCE.md`, `docs/THIRD_PARTY.md`, `docs/CODE_STYLE.md`, `docs/history/VALIDATION_HISTORY.md`, `AGENTS.md`, `CONTRIBUTING.md`, adapted from [atcc-dev-template](https://github.com/ATCC-Bioinformatics/atcc-dev-template).
- Added `ruff.toml` and an advisory GitHub Actions style-check workflow (not required for merging yet).

Verification limits as of this entry: the restructuring was checked with the unit test suite and a partial live run against NCBI/iCite (see `docs/history/VALIDATION_HISTORY.md`); no full end-to-end run producing both output TSVs has completed yet, and the fuzzy-matching accuracy requirement (R5) is still an open gap.
