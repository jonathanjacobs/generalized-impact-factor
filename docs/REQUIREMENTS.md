# Requirements

- **Audience:** Anyone deciding whether a change is acceptable, or reviewing this project's output before it supports a conclusion.
- **Use this when:** Deciding whether the tool's behavior or a proposed change to it is acceptable.
- **Update this when:** What counts as acceptable behavior changes.
- **Do not update for:** Current tasks or status — those belong in the work tracker.

## Requirements

### R1 — Deterministic index math
Given a fixed list of `{PMID, CitedByCount}` pairs, `compute_h_g_i10` must always return the same H-index, G-index, and i10-index, independent of input order.
**We will know this is met when:** the unit tests in `tests/test_metrics.py` pass, including the input-order test.

### R2 — No missing rows
A term with zero PMC hits or zero convertible PMIDs must still produce a row in `terms_metrics.tsv` with `Total Mentions`, `H-index`, `G-index`, and `i10-index` all `0` — it must not be silently dropped from the output.
**We will know this is met when:** running the CLI on a term known to have no matches produces a row for that term, not a missing one (see the "no PMIDs found" branch in `impact_lookup.py`).

### R3 — No double-counted PMIDs
A PMID matched via more than one fuzzy query variation for the same term must be counted once, not once per variation, in `Total Mentions` and in the citation list used for H/G/i10.
**We will know this is met when:** `pmc_to_pmid_idconv` returns a de-duplicated PMID list (it already does, via a `set`); a regression test should be added if this behavior is ever changed.

### R4 — Standalone, reusable library
The functions behind the CLI must be importable and independently callable without invoking `argparse` or running the CLI's `main()`.
**We will know this is met when:** `from impact_factor import ...` works from a script or notebook that only adds `src/` to its path, as demonstrated in `src/README.md` and exercised by every test in `tests/`.

### R5 — Fuzzy-matching evidence, not just plausibility
Before this tool's output is used to support a scientific or business conclusion, the fuzzy-matching false-positive and false-negative rate must be checked against at least one hand-reviewed sample of real search results, not just assumed to be acceptable.
**We will know this is met when:** a dated entry exists in `docs/history/VALIDATION_HISTORY.md` describing that check and its result.

## Out of scope

- A web interface or scheduled/automated runs.
- Explaining *why* a term is or isn't cited — the tool counts and scores, it does not interpret.
- Full-text mining of PMC articles; only bibliographic metadata (title, journal, year, authors, affiliations) is retrieved.

## Maintaining the agreement

Keep this list short enough that everyone involved can hold it in their head. Add a requirement when a real disagreement or surprise reveals that one was missing; do not add one speculatively. Update the acceptance criterion, not just the requirement text, when the evidence location changes (e.g. once `docs/TESTING.md` and `docs/history/VALIDATION_HISTORY.md` are more built out).
