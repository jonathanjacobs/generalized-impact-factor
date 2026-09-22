#!/usr/bin/env python3
"""ATCC Term citation analysis pipeline: fuzzy-search PMC, fetch citation counts from iCite,
and compute per-term H-index, G-index, and i10-index.

Run directly (no installation required): python impact_lookup.py --email you@example.com
"""

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

import argparse
import logging
import os
import sys
import time
from datetime import datetime
from pathlib import Path

# Make the reusable library importable without installing the package, so this
# stays a standalone script per docs/DESIGN.md.
sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from Bio import Entrez  # noqa: E402

from impact_factor import (  # noqa: E402
    compute_h_g_i10,
    dump_term_citations_tsv,
    get_citation_counts_icite,
    pmc_to_pmid_idconv,
    read_terms_from_file,
    search_pmc,
    write_term_metrics_tsv,
)


def build_arg_parser() -> argparse.ArgumentParser:
    """Build the CLI argument parser."""
    parser = argparse.ArgumentParser(
        description="ATCC Term citation analysis pipeline."
    )

    parser.add_argument(
        "--input_file",
        type=str,
        default="terms.txt",
        help="Path to the file containing terms (e.g., SKUs), one per line.",
    )
    parser.add_argument(
        "--output_metrics_file",
        type=str,
        default="terms_metrics.tsv",
        help="Path to the TSV file for outputting per-term citation metrics.",
    )
    parser.add_argument(
        "--output_citations_file",
        type=str,
        default="terms_citations.tsv",
        help="Path to the TSV file for outputting per-citation metadata.",
    )
    parser.add_argument(
        "--email",
        type=str,
        default=os.getenv("NCBI_EMAIL"),
        help="Email address required by NCBI Entrez for identification. Default is to pull NCBI_EMAIL from the ENV",
    )
    parser.add_argument(
        "--tool_name",
        type=str,
        default="ATCC-Term-CitationMetrics",
        help="Tool name required by NCBI Entrez for identification.",
    )
    parser.add_argument(
        "--api_key",
        type=str,
        default=None,
        help="Optional: NCBI API key. If not provided, will check environment variable NCBI_API_KEY.",
    )
    parser.add_argument(
        "--window_length_years",
        type=int,
        default=5,
        help="Number of years to look back for publications.",
    )
    parser.add_argument(
        "--end_offset_years",
        type=int,
        default=0,
        help=(
            "Offset from the current year for the end of the search window (e.g., 0 for current "
            "year, 1 for last year). Default is 0 for current year."
        ),
    )
    parser.add_argument(
        "--retmax",
        type=int,
        default=200,
        help="Maximum number of PMC hits per term search.",
    )
    parser.add_argument(
        "--idconv_batch_size",
        type=int,
        default=100,
        help="Maximum number of IDs to convert in a single NCBI ID Converter API call.",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=20,
        help="Chunk size for Entrez queries (NCBI likes <= 20).",
    )
    parser.add_argument(
        "--timedelay",
        type=float,
        default=0.34,
        help="Delay in seconds between API requests (conservative: ~3 requests/sec without API key).",
    )
    parser.add_argument(
        "--maxdelay",
        type=int,
        default=20,
        help="Maximum delay in seconds for failed retries of NCBI requests.",
    )
    parser.add_argument(
        "--icite_batch_size",
        type=int,
        default=200,
        help="Safe batch size for iCite API calls.",
    )
    parser.add_argument(
        "--additional_query_term",
        type=str,
        default=None,
        help="Optional: An additional term to be ANDed with every search query.",
    )
    return parser


def configure_logging() -> None:
    """Configure root logging to stdout, replacing any existing handlers."""
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[logging.StreamHandler(sys.stdout)],
    )


def main() -> None:
    configure_logging()

    parser = build_arg_parser()
    # parser.parse_args([]) is Colab-friendly: pass an empty list when running under ipykernel.
    args = (
        parser.parse_args([])
        if "ipykernel" in sys.modules
        else parser.parse_args()
    )

    Entrez.email = args.email
    Entrez.tool = args.tool_name

    if not args.email:
        parser.error("--email is required (or set NCBI_EMAIL)")

    if args.api_key:
        Entrez.api_key = args.api_key
    elif os.getenv("NCBI_API_KEY"):
        Entrez.api_key = os.getenv("NCBI_API_KEY")

    # Patient with transient NCBI failures (Biopython will retry automatically)
    Entrez.max_tries = 5
    Entrez.sleep_between_tries = 30

    start_time = datetime.now()
    start_clock = time.time()
    logging.info(f"Script start: {start_time}")

    terms = read_terms_from_file(args.input_file)

    if not terms:
        logging.info(
            f"No terms found in '{args.input_file}'. Using default terms!"
        )
        terms = ["ATCC CRL-1585", "VR-1490"]

    logging.info(f"Terms to process: {terms}")

    all_terms_metrics = []

    for term in terms:
        logging.info(f"\n--- Processing Term: {term} ---")

        pmc_ids = search_pmc(
            term,
            args.window_length_years,
            args.end_offset_years,
            args.retmax,
            args.additional_query_term,
        )
        logging.info(
            f"Found {len(pmc_ids)} PMC IDs for {term} from a {args.window_length_years}-year "
            f"window ending {args.end_offset_years} year(s) ago."
        )

        pubmed_ids = pmc_to_pmid_idconv(
            pmc_ids,
            Entrez.tool,
            Entrez.email,
            args.idconv_batch_size,
            Entrez.max_tries,
            args.timedelay,
            args.maxdelay,
        )

        logging.info(
            f"Converted {len(pmc_ids)} PMCIDs to {len(pubmed_ids)} unique PMIDs for {term}."
        )
        total_mentions = len(pubmed_ids)

        if not pubmed_ids:
            logging.info(f"No PMIDs found for Term {term}.")
            all_terms_metrics.append(
                {
                    "Term": term,
                    "Total Mentions": 0,
                    "H-index": 0,
                    "G-index": 0,
                    "i10-index": 0,
                }
            )
            continue

        dump_term_citations_tsv(
            term=term,
            pmids=pubmed_ids,
            out_tsv_path=args.output_citations_file,
            batch_size=args.retmax,
            request_delay=args.timedelay,
        )

        citation_data = get_citation_counts_icite(
            pubmed_ids,
            args.icite_batch_size,
            Entrez.max_tries,
            args.timedelay,
            args.maxdelay,
        )

        h, g, i10 = compute_h_g_i10(citation_data)

        logging.info(f"\nMetrics for Term '{term}':")
        logging.info(f"  H-index: {h}")
        logging.info(f"  G-index: {g}")
        logging.info(f"  i10-index: {i10}")

        all_terms_metrics.append(
            {
                "Term": term,
                "Total Mentions": total_mentions,
                "H-index": h,
                "G-index": g,
                "i10-index": i10,
            }
        )

    logging.info("\n--- Summary of All Terms ---")
    for entry in all_terms_metrics:
        logging.info(entry)

    write_term_metrics_tsv(
        all_terms_metrics, out_tsv_path=args.output_metrics_file
    )

    end_time = datetime.now()
    elapsed_sec = time.time() - start_clock

    logging.info(f"Script end: {end_time}")
    logging.info(
        f"Elapsed time: {elapsed_sec:.1f} seconds ({elapsed_sec / 60:.2f} minutes)"
    )


if __name__ == "__main__":
    main()
