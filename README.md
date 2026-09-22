# generalized-impact-factor

Measures the "scientific impact" of technologies and products used in the life sciences — cell lines, bacterial strains, instruments, reagents, and similar terms — using the same bibliometric indices normally used to score author impact (H-index, G-index, i10-index), applied instead to how often and how influentially a term appears in the biomedical literature.

- **Status:** in development. The core pipeline works and has unit-test coverage (see `tests/`); it has not yet had a scientific review of its output against known cases, and end-to-end network behavior has not been checked into `docs/history/VALIDATION_HISTORY.md`. Do not treat its numbers as validated until that review has happened.
- **Started from:** [atcc-dev-template](https://github.com/ATCC-Bioinformatics/atcc-dev-template) at commit `b36531f` (the template has no version tags yet).
- **Project Lead / Technical Lead:** Jonathan L. Jacobs (jjacobs@atcc.org).
- **Work tracker:** TBD — link the approved tracker here once one exists for this project.

See [`docs/DESIGN.md`](docs/DESIGN.md) for why this exists and the approach, [`docs/REQUIREMENTS.md`](docs/REQUIREMENTS.md) for what counts as an acceptable result, and [`docs/TESTING.md`](docs/TESTING.md) for what evidence exists so far.

## Features
- **Term-based search:** searches PubMed Central for publications mentioning specific terms.
- **Fuzzy matching:** automatically generates variations of input terms (with/without spaces, hyphens) to broaden the search.
- **Date windowing:** filters search results to a specified publication date range.
- **PMID conversion:** converts PMC IDs to PubMed IDs (PMIDs) for consistent citation tracking.
- **Citation metrics:** fetches citation counts for PMIDs using the NIH iCite API.
- **Bibliometric indices:** calculates H-index, G-index, and i10-index for each term based on citation data.
- **Detailed output:** generates TSV files with per-term metrics and detailed citation metadata.

## Repository layout
- `impact_lookup.py` — the standalone CLI. Run it directly; it adds `src/` to its own import path, so no install step is required.
- `src/impact_factor/` — the reusable library behind the CLI (term matching, PMC/PubMed/iCite calls, index math). See [`src/README.md`](src/README.md). Other scripts can import from it the same way.
- `tests/` — pytest unit tests for the library's pure functions.
- `docs/` — project definition, data-handling and third-party terms, and testing evidence.
- `environments/` — pinned Python dependencies.

## Install
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r environments/requirements.txt
```

An NCBI API key is highly recommended for higher request limits and improved reliability, though not strictly required for basic use. You can obtain one from your [NCBI account](https://ncbi.nlm.nih.gov/account/).

## Usage

### 1. Create an input file
Create a plain text file (e.g., `terms.txt`) where each line contains a single term you want to analyze. If this file is not found, the script uses two default terms for demonstration.

Example `terms.txt`:
```
ATCC CRL-1585
VR-1490
Human Papillomavirus
SARS-CoV-2
```

Before committing a real `terms.txt`, see the note on the input term list in [`docs/DATA_GOVERNANCE.md`](docs/DATA_GOVERNANCE.md).

### 2. Run the script

```bash
python impact_lookup.py \
    --input_file terms.txt \
    --output_metrics_file my_metrics.tsv \
    --output_citations_file my_citations.tsv \
    --email your.email@example.com \
    --api_key YOUR_NCBI_API_KEY \
    --window_length_years 10 \
    --end_offset_years 0 \
    --additional_query_term "(\"cell line\" OR virus)"
```

**Arguments:**
- `--input_file` (str, default: `terms.txt`): path to the file containing terms, one per line.
- `--output_metrics_file` (str, default: `terms_metrics.tsv`): path for the output TSV file with per-term citation metrics.
- `--output_citations_file` (str, default: `terms_citations.tsv`): path for the output TSV file with detailed per-citation metadata.
- `--email` (str, **required**): your email address for NCBI Entrez identification. Defaults to the `NCBI_EMAIL` environment variable if set.
- `--tool_name` (str, default: `ATCC-Term-CitationMetrics`): tool name for NCBI Entrez identification.
- `--api_key` (str, optional): your NCBI API key. If not provided, checks the `NCBI_API_KEY` environment variable.
- `--window_length_years` (int, default: 5): number of years to look back for publications.
- `--end_offset_years` (int, default: 0): offset from the current year for the end of the search window (0 for current year, 1 for last year).
- `--retmax` (int, default: 200): maximum number of PMC hits per term search.
- `--idconv_batch_size` (int, default: 100): batch size for NCBI ID Converter API calls.
- `--chunksize` (int, default: 20): chunk size for Entrez queries (NCBI recommends <= 20).
- `--timedelay` (float, default: 0.34): delay in seconds between API requests to be polite to NCBI services.
- `--maxdelay` (int, default: 20): maximum delay for failed retries of NCBI requests.
- `--icite_batch_size` (int, default: 200): safe batch size for iCite API calls.
- `--additional_query_term` (str, optional): an additional term ANDed with every search query (e.g., `"(\"cell line\" OR virus)"` for more focused results).

### Google Colab execution
If running in Google Colab, add `src/` to the path the same way `impact_lookup.py` does, then call its `main()`, or copy the argument-parsing and pipeline code directly into a cell. `parser.parse_args([])` is used automatically whenever `ipykernel` is detected, so command-line syntax is not required in a notebook — set values on the `args` object instead:

```python
# args.input_file = 'my_custom_terms.txt'
# args.email = 'your.email@example.com'
# args.window_length_years = 10
# args.api_key = "YOUR_NCBI_API_KEY"  # or set as an environment variable in Colab secrets
```

## Output files
The script generates two TSV (tab-separated value) files. **These are run results, not repository content — see [`docs/DATA_GOVERNANCE.md`](docs/DATA_GOVERNANCE.md); they are `.gitignore`d and should not be committed.**

1. **`terms_metrics.tsv`** (or your specified `--output_metrics_file`)
   Aggregated bibliometric indices for each term:
   - `Term`: the analyzed term.
   - `Total Mentions`: total number of unique PubMed IDs found.
   - `H-index`: Hirsch index.
   - `G-index`: G-index.
   - `i10-index`: number of publications with at least 10 citations.

2. **`terms_citations.tsv`** (or your specified `--output_citations_file`)
   Detailed metadata for each found citation:
   - `Term`: the term associated with this citation.
   - `PubMedID`: the PubMed ID of the article.
   - `Year`: publication year.
   - `Title`: article title.
   - `Journal`: journal name.
   - `AuthorList`: semicolon-separated list of authors.
   - `Institution`: pipe-separated list of affiliated institutions.

## How it works
1. **Read terms:** loads terms from the specified input file.
2. **PMC search:** for each term, searches PubMed Central using `Entrez.esearch`. The query includes fuzzy variations of the term and a date filter, plus an optional `additional_query_term`.
3. **ID conversion:** converts the retrieved PMC IDs to PMIDs using the NCBI PMC ID Converter API, for compatibility with PubMed and iCite.
4. **Fetch citation details:** for each PMID, uses `Entrez.efetch` to retrieve article metadata (title, authors, journal, year, institutions) from PubMed.
5. **iCite metrics:** queries the NIH iCite API for citation counts for each PMID.
6. **Calculate indices:** computes the H-index, G-index, and i10-index from the citation counts.
7. **Output:** writes the aggregated term metrics and detailed citation data to separate TSV files.

Known limitations: fuzzy term matching can both over-match (an unrelated paper happens to contain the exact variant string) and under-match (a paper phrases the term differently than any generated variation); iCite's citation-count coverage starts in 1980; and PMC's index of full text lags and does not perfectly overlap with PubMed's index of metadata-only records. See `docs/DESIGN.md` for how these affect the intended use.

## License
This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
