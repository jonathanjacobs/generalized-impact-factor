# generalized-impact-factor
Measures the "scientific impact" of technologies and products used in the life sciences, such as cell lines, bacterial strains, instruments, reagents, etc. using metrics similar determining author impact factors (i.e. H-index, G-index, i10-index).

## Table of Contents
- [Features](#features)
- [Requirements](#requirements)
- [Usage](#usage)
  - [Input File Format](#input-file-format)
  - [Command-line Arguments](#command-line-arguments)
- [Output Files](#output-files)
- [How it Works](#how-it-works)
- [License](#license)

## Features
- **Term-based Search**: Searches PubMed Central for publications mentioning specific terms.
- **Fuzzy Matching**: Automatically generates variations of input terms (e.g., with/without spaces, hyphens) to ensure comprehensive search results.
- **Date Windowing**: Filters search results to a specified publication date range.
- **PMID Conversion**: Converts PMC IDs to PubMed IDs (PMIDs) for consistent citation tracking.
- **Citation Metrics**: Fetches citation counts for PMIDs using the iCite API.
- **Bibliometric Indices**: Calculates H-index, G-index, and i10-index for each term based on citation data.
- **Detailed Output**: Generates TSV files with per-term metrics and detailed citation metadata.

## Requirements
- Python 3.x
- Biopython library (`pip install biopython`)
- requests library (`pip install requests`)
- titlecase library (`pip install titlecase`)

An NCBI API key is highly recommended for higher request limits and improved reliability, though not strictly required for basic use. You can obtain one from your [NCBI account](https://ncbi.nlm.nih.gov/account/).

## Usage

### 1. Create an Input File
Create a plain text file (e.g., `terms.txt`) where each line contains a single term you want to analyze. If this file is not found, the script will use default terms for demonstration.

Example `terms.txt`:
```
ATCC CRL-1585
VR-1490
Human Papillomavirus
SARS-CoV-2
```

### 2. Run the Script

You can run the script from the command line or within an environment like Google Colab. The script uses `argparse` to handle configuration.

**Command-line execution:**

```bash
python your_script_name.py \
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
- `--input_file` (str, default: `terms.txt`): Path to the file containing terms, one per line.
- `--output_metrics_file` (str, default: `terms_metrics.tsv`): Path for the output TSV file with per-term citation metrics.
- `--output_citations_file` (str, default: `terms_citations.tsv`): Path for the output TSV file with detailed per-citation metadata.
- `--email` (str, **required**): Your email address for NCBI Entrez identification.
- `--tool_name` (str, default: `ATCC-Term-CitationMetrics`): Tool name for NCBI Entrez identification.
- `--api_key` (str, optional): Your NCBI API key. If not provided, it will check the `NCBI_API_KEY` environment variable.
- `--window_length_years` (int, default: 5): Number of years to look back for publications.
- `--end_offset_years` (int, default: 0): Offset from the current year for the end of the search window (0 for current year, 1 for last year).
- `--retmax` (int, default: 200): Maximum number of PMC hits per term search.
- `--idconv_batch_size` (int, default: 100): Batch size for NCBI ID Converter API calls.
- `--chunksize` (int, default: 20): Chunk size for Entrez queries (NCBI recommends <= 20).
- `--timedelay` (float, default: 0.34): Delay in seconds between API requests to be polite to NCBI services.
- `--maxdelay` (int, default: 20): Maximum delay for failed retries of NCBI requests.
- `--icite_batch_size` (int, default: 200): Safe batch size for iCite API calls.
- `--additional_query_term` (str, optional): An additional term to be ANDed with every search query (e.g., `"(\"cell line\" OR virus)"` for more focused results).

### Google Colab execution:
If running in Google Colab, the `parser.parse_args([])` line is set to handle this by passing an empty list. You can modify the `args` object directly after its creation if you need to override defaults programmatically without using command-line syntax.

```python
# ... (existing code) ...

args = parser.parse_args([]) if "ipykernel" in sys.modules else parser.parse_args()

# Example of overriding arguments in Colab:
# args.input_file = 'my_custom_terms.txt'
# args.email = 'your.email@example.com'
# args.window_length_years = 10
# args.api_key = "YOUR_NCBI_API_KEY" # Or set as an environment variable in Colab secrets

# ... (rest of the script) ...
```

## Output Files
The script generates two TSV (Tab Separated Value) files:

1.  **`terms_metrics.tsv`** (or your specified `--output_metrics_file`)
    Contains aggregated bibliometric indices for each term:
    -   `Term`: The analyzed term.
    -   `Total Mentions`: Total number of unique PubMed IDs found.
    -   `H-index`: Hirsch index.
    -   `G-index`: G-index.
    -   `i10-index`: Number of publications with at least 10 citations.

2.  **`terms_citations.tsv`** (or your specified `--output_citations_file`)
    Contains detailed metadata for each found citation:
    -   `Term`: The term associated with this citation.
    -   `PubMedID`: The PubMed ID of the article.
    -   `Year`: Publication year.
    -   `Title`: Article title.
    -   `Journal`: Journal name.
    -   `AuthorList`: Semicolon-separated list of authors.
    -   `Institution`: Pipe-separated list of affiliated institutions.

## How it Works
1.  **Read Terms**: Loads terms from the specified input file.
2.  **PMC Search**: For each term, performs a search on PubMed Central using `Entrez.esearch`. The search query includes fuzzy variations of the term and a date filter. An optional `additional_query_term` can be appended.
3.  **ID Conversion**: Converts the PMC IDs retrieved from PMC into PMIDs using the NCBI ID Converter API. This step ensures compatibility with PubMed and iCite services.
4.  **Fetch Citation Details**: For each PMID, uses `Entrez.efetch` to retrieve detailed article metadata (title, authors, journal, year, institutions) from PubMed.
5.  **iCite Metrics**: Queries the NIH iCite API to get citation counts for each PMID.
6.  **Calculate Indices**: Computes the H-index, G-index, and i10-index based on the collected citation counts.
7.  **Output**: Writes the aggregated term metrics and detailed citation data to separate TSV files.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

```
# MIT License
#
# Copyright (c) 2024 Jonathan L. Jacobs
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
```
