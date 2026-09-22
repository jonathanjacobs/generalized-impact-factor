# Tests

Pytest unit tests for the pure functions in `src/impact_factor/`: `compute_h_g_i10` (H-index/G-index/i10-index math, with hand-verified reference values), `generate_fuzzy_queries` (term-variation generation), `smart_tc` (acronym/alphanumeric-preserving titlecase), `chunks`, and `read_terms_from_file`.

Run with:

```bash
python3 -m pytest tests/
```

No test data files are used; all inputs are small literals defined in the test files. This suite does not exercise the network-calling functions (`search_pmc`, `pmc_to_pmid_idconv`, `dump_term_citations_tsv`, `get_citation_counts_icite`) — see `docs/TESTING.md` for how those are checked.
