# Validation history

- **Audience:** Anyone deciding whether to rely on a past result from this project.
- **Use this when:** Checking what has already been validated, and under what conditions.
- **Update this when:** A new dated conclusion is reached. Add a new entry; do not edit or remove an earlier one — add a correction or superseding entry instead.
- **Do not update for:** A routine passing test run; that belongs in CI output, not here.

## 2026-09-22 — Restructuring into src/impact_factor did not change pipeline behavior

- **Code or workflow version:** commit at the time of the `src/impact_factor` library split (see `git log` around this date); pre-split baseline was the single-file `impact_lookup.py` as of commit `c9e1731`.
- **Settings, input, and reference identifiers:** terms `ATCC CRL-1585` and `SARS-CoV-2 spike protein`; `--window_length_years 2` or `3`, small `--retmax` (5–15); live NCBI E-utilities and NIH iCite APIs, `NCBI_API_KEY` from the environment.
- **Environment and checks performed:** ran on the `an internal ATCC HPC host` host, Python 3.9.19. Confirmed: (1) `python -m pytest tests/` — 17/17 pass, covering `compute_h_g_i10`, `generate_fuzzy_queries`, `smart_tc`, `chunks`, `read_terms_from_file`; (2) the refactored `search_pmc` returns real PMC IDs for a known term (`ATCC CRL-1585` → 5 PMC IDs on a direct call); (3) `impact_lookup.py --help` and the missing-`--email` error path behave identically to the original CLI's documented interface; (4) direct package import (`from impact_factor import ...`) works without going through the CLI, confirming R4.
- **Observed result:** the library split, `main()` extraction, and the two bug fixes (dead `elapsed_sec`/`start_clock`, and `dump_term_citations_tsv`'s previously-unused `request_delay` parameter now actually controls its pacing) did not change the CLI's argument surface or the PMC query construction. A full end-to-end run producing both output TSVs was not completed in this session — live PMC search calls began returning `HTTP 400: Bad Request` partway through testing, traced to NCBI-side throttling on the shared `NCBI_API_KEY` (confirmed by reproducing the same 400 on an unmodified, minimal `Entrez.esearch` call with no code from this project involved).
- **Limitations and skipped or inconclusive checks:** no check yet of fuzzy-matching false-positive/false-negative rate (`R5`, still open — see `docs/REQUIREMENTS.md` and `docs/TESTING.md`); no completed full pipeline run (PMC search through both TSV outputs) in this session due to the NCBI throttling above; no check of `dump_term_citations_tsv` or `get_citation_counts_icite` against live data in this session.
- **Evidence location:** this entry; the pytest run is reproducible via `python -m pytest tests/` in this repository.
- **Reviewer and decision:** not yet reviewed by a person. The checks above were performed during the template-adoption restructuring but have not been confirmed by Jonathan Jacobs (Project Lead / Technical Lead). A full live end-to-end run and the R5 fuzzy-matching check are still needed before treating this tool's output as validated for a scientific or business conclusion.
