# Design

- **Audience:** Anyone deciding whether this project's approach fits their need, or reviewing a proposed change to it.
- **Use this when:** Understanding why this exists and what it does and does not try to do.
- **Update this when:** The problem, scope, or proposed approach changes.
- **Do not update for:** Implementation details that belong in `ARCHITECTURE.md` or the code itself.

## Problem and importance

- **Problem, need, or opportunity:** ATCC catalog items — cell lines, strains, instruments, reagents — appear in the biomedical literature at very different rates and with very different downstream influence, but there is no existing measure of that influence comparable across items. Author-level bibliometrics (H-index, G-index, i10-index) already solve an analogous problem for people; this project applies the same math to product terms instead of author names.
- **Why it matters:** these numbers are intended to support scientific or business conclusions (e.g. which products are most cited, or most influential once cited) — see the audience and stakes recorded in `README.md` and the project's adoption-interview answers.
- **Who is affected:** whoever uses the output to justify a claim or a decision, and downstream readers of that report or conclusion.

**Project description:** a command-line tool that, given a list of terms, searches PubMed Central for mentions of each term (with fuzzy matching for spacing/hyphenation), converts matches to PMIDs, retrieves citation counts from NIH iCite, and reports H-index/G-index/i10-index and total mention count per term.

## Intended outcome and scope

**Outcome and scope:** produce a reusable, importable library (`src/impact_factor/`) plus a standalone CLI (`impact_lookup.py`) that other bioinformaticians and data scientists can run themselves or import functions from directly, without requiring a package install. Out of scope for now: a web interface, scheduled/automated runs, or any claim about *why* a term is or isn't cited (the tool counts and scores; it does not explain).

## Proposed approach

A term-by-term pipeline: (1) build fuzzy query variations of the term and search PMC within a configurable publication-date window; (2) convert the resulting PMC IDs to PMIDs via NCBI's ID Converter API; (3) fetch per-PMID metadata (title, journal, year, authors, affiliations) from PubMed; (4) fetch per-PMID citation counts from NIH iCite; (5) compute H-index, G-index, and i10-index from those counts. Steps 1–4 are implemented as small, independently testable functions in `src/impact_factor/`; `impact_lookup.py` only handles argument parsing and wiring the steps together, so the same functions can be reused from another script or a notebook.

## Important considerations

- **Scientific assumptions, references, and interpretation limits:**
  - Fuzzy matching (space/hyphen variations of a term) can both over-match — an unrelated paper happens to contain the literal string — and under-match — a paper phrases the term in a way no generated variation covers. Nobody has yet measured the project's false-positive/false-negative rate against a hand-checked sample; see `docs/REQUIREMENTS.md` and `docs/TESTING.md`.
  - NIH iCite's citation-count coverage begins in 1980; earlier publications will be undercounted or show `0` citations, which affects H/G/i10 for terms with an older literature history.
  - PMC (full text index) and PubMed (metadata index) do not perfectly overlap; a term matched in PMC is assumed to also resolve to a PMID via the ID Converter, which is not guaranteed for every record.
- **Data-handling and security:** covered in `docs/DATA_GOVERNANCE.md` — in particular, whether a real input term list (actual ATCC SKUs/product names) is safe to commit to this repository's current, non-ATCC-org-hosted remote.
- **Computing-environment, portability, and reproducibility:** tested only on Python 3.9.19 so far (see `environments/README.md`); network calls to NCBI/iCite mean a run is not exactly reproducible byte-for-byte over time (citation counts and PMC's index both change), so a run record (`docs/RUN_PROVENANCE.md`) matters more here than for a purely deterministic tool.
- **Operational, compatibility, and support considerations:** not yet an operational service — see `README.md` status. Revisit this document if that changes.
- **Known limitations, alternatives, and open questions:** no alternative approach (e.g. a different citation database, or exact-match-only search) has been evaluated yet; this is a reasonable next investigation if the fuzzy-matching false-positive rate turns out to be a problem.

Considerations required for acceptance belong in `docs/REQUIREMENTS.md`. How the parts connect is only in `docs/ARCHITECTURE.md` if that document is ever added — for now, `src/README.md` covers the two-part CLI/library split, which is simple enough not to need a separate architecture document.
