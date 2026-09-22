# Testing and validation

- **Audience:** Anyone deciding whether to trust this tool's output, or contributing a change.
- **Use this when:** Deciding what evidence is needed, performing checks, or evaluating readiness.
- **Update this when:** The planned checks, environments, acceptance criteria, or responsible people change.
- **Do not update for:** One observed result; retain that evidence in `docs/history/VALIDATION_HISTORY.md`.

## What is covered today

`tests/` has pytest unit tests for the pure, deterministic parts of `src/impact_factor/`:

| Check | What it covers |
| --- | --- |
| `tests/test_metrics.py` | `compute_h_g_i10` against hand-verified citation-count cases, including empty input, all-zero citations, one highly-cited paper, and input-order independence (R1). |
| `tests/test_terms.py` | `generate_fuzzy_queries` variation generation for terms with spaces, hyphens, both, or neither; `read_terms_from_file` reading and blank-line handling. |
| `tests/test_text.py` | `smart_tc` acronym/alphanumeric-code preservation; `chunks` batching. |

Run with `python3 -m pytest tests/`.

## What is not covered

- **The network-calling functions** (`search_pmc`, `pmc_to_pmid_idconv`, `dump_term_citations_tsv`, `get_citation_counts_icite`) have no automated test — they were checked manually against live NCBI/iCite services during the `src/` restructuring (see `docs/history/VALIDATION_HISTORY.md` once that entry is added), but there is no regression test that would catch, for example, an NCBI response-format change.
- **Fuzzy-matching accuracy** (R5 in `docs/REQUIREMENTS.md`): nobody has checked what fraction of matched PMC hits for a real term are false positives, or estimated the false-negative rate from phrasing variations the fuzzy matcher doesn't generate. This is the most important open gap before trusting these numbers in a report.
- **iCite coverage limits**: no check exists for how citation counts behave for pre-1980 publications (iCite's stated coverage start).

## Connect requirements to evidence

| Requirement | Evidence needed | Responsible person | Status |
| --- | --- | --- | --- |
| R1 | `tests/test_metrics.py` | Jonathan Jacobs | done |
| R2 | Manual CLI run against a zero-hit term | Jonathan Jacobs | done (see `docs/history/VALIDATION_HISTORY.md`) |
| R3 | Covered incidentally by `pmc_to_pmid_idconv`'s use of a `set`; no dedicated regression test | TBD | gap |
| R4 | `tests/` (every test imports the package directly, not through the CLI) | Jonathan Jacobs | done |
| R5 | Hand-reviewed sample of real search results | TBD | gap — not started |

## Small test inputs and environments

All current unit tests use small literal values defined directly in the test files — no test data files are needed, and none are versioned in `tests/`.

## Who decides evidence is sufficient

TBD — until named otherwise, treat this as Jonathan Jacobs (Project Lead / Technical Lead, per `README.md`).
