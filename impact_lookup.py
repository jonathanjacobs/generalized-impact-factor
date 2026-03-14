import sys
import argparse
import random
import time
from datetime import datetime
from Bio import Entrez
import requests
import csv
import os
import re
from titlecase import titlecase
import logging

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

# --- Configure Logging ---
# Remove any existing handlers from the root logger to prevent duplicate output
for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)

logging.basicConfig(
    level=logging.INFO, # Set to logging.DEBUG for more verbose output
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.StreamHandler(sys.stdout) # Direct output to stdout for Colab visibility
    ]
)

# --- Argument Parsing ---
parser = argparse.ArgumentParser(description="ATCC Term citation analysis pipeline.")

parser.add_argument(
    "--input_file",
    type=str,
    default='terms.txt',
    help="Path to the file containing terms (e.g., SKUs), one per line."
)
parser.add_argument(
    "--output_metrics_file",
    type=str,
    default="terms_metrics.tsv",
    help="Path to the TSV file for outputting per-term citation metrics."
)
parser.add_argument(
    "--output_citations_file",
    type=str,
    default="terms_citations.tsv",
    help="Path to the TSV file for outputting per-citation metadata."
)

parser.add_argument(
    "--email",
    type=str,
    #default=os.getenv("NCBI_EMAIL"),
    default="jjacobs@atcc.org",
    help="Email address required by NCBI Entrez for identification. Default is is to pull NCBI_EMAIL from the ENV"
)

parser.add_argument(
    "--tool_name",
    type=str,
    default="ATCC-Term-CitationMetrics", # Default tool name
    help="Tool name required by NCBI Entrez for identification."
)
parser.add_argument(
    "--api_key",
    type=str,
    default=None, # It's better to keep API key from env/secrets rather than default in code
    help="Optional: NCBI API key. If not provided, will check environment variable NCBI_API_KEY."
)
parser.add_argument(
    "--window_length_years",
    type=int,
    default=5,
    help="Number of years to look back for publications."
)
parser.add_argument(
    "--end_offset_years",
    type=int,
    default=0,
    help="""Offset from the current year for the end of the search window (e.g., 0 for current year, 1 for last year).
    Default is 0 for current year."""
)
parser.add_argument(
    "--retmax",
    type=int,
    default=200,
    help="Maximum number of PMC hits per term search."
)
parser.add_argument(
    "--idconv_batch_size",
    type=int,
    default=100,
    help="Maximum number of IDs to convert in a single NCBI ID Converter API call."
)
parser.add_argument(
    "--chunksize",
    type=int,
    default=20,
    help="Chunk size for Entrez queries (NCBI likes <= 20)."
)
parser.add_argument(
    "--timedelay",
    type=float,
    default=0.34,
    help="Delay in seconds between API requests (conservative: ~3 requests/sec without API key)."
)
parser.add_argument(
    "--maxdelay",
    type=int,
    default=20,
    help="Maximum delay in seconds for failed retries of NCBI requests."
)
parser.add_argument(
    "--icite_batch_size",
    type=int,
    default=200,
    help="Safe batch size for iCite API calls."
)
parser.add_argument(
    "--additional_query_term",
    type=str,
    default=None,
    help="Optional: An additional term to be ANDed with every search query."
)

args = parser.parse_args([]) if "ipykernel" in sys.modules else parser.parse_args()
# Above is Google colab friendly. For Colab execution, pass an empty list if not running from CLI

# Globals (now derived from args or defaults)
input_terms_file = args.input_file
out_metrics_tsv = args.output_metrics_file
out_citations_tsv = args.output_citations_file
additional_query_term = args.additional_query_term

# REQUIRED by NCBI (identify yourself)
Entrez.email = args.email
Entrez.tool = args.tool_name

# EMAIL IS REQUIRED
if not args.email:
    parser.error("--email is required (or set NCBI_EMAIL)")

# OPTIONAL: if you have one, set in env var NCBI_API_KEY and uncomment:
if args.api_key:
    Entrez.api_key = args.api_key
elif os.getenv("NCBI_API_KEY"):
    Entrez.api_key = os.getenv("NCBI_API_KEY")

# Patient with transient NCBI failures (Biopython will retry automatically)
Entrez.max_tries = 5
Entrez.sleep_between_tries = 30

# search parameters
window_length_years = args.window_length_years
end_offset_years = args.end_offset_years
retmax = args.retmax
IDCONV_BATCH_SIZE = args.idconv_batch_size

# pacing
chunksize = args.chunksize
timedelay = args.timedelay
maxdelay = args.maxdelay

# ICITE GLOBALS
ICITE_BASE = "https://icite.od.nih.gov/api/pubs"
ICITE_BATCH_SIZE = args.icite_batch_size # Re-added here

# Utilities

# titlecase tool
def smart_tc(s: str) -> str:
  if not s:
    return ""
  #preserve biomedical/sciency stuff and acrynyms...
  def keep_acronyms(word,**kwargs):
    # preserve ALLCAPS, etc and common patterns with digits
    if word.isupper():
      return word
    if any(ch.isdigit() for ch in word) and any(ch.isalpha() for ch in word):
      return word
    return None

  return titlecase(s,callback=keep_acronyms)

# chunking tool
def chunks(lst, n):
        for i in range(0, len(lst), n):
            yield lst[i:i+n]

#deal with NCBI being pissy
def backoff_sleep(attempt, base=timedelay, cap=maxdelay):
    # Exponential backoff with jazzy jittering.
    # this should now only be called if first attempt fails / throws exception and we want a second try
    delay = min(cap, (base * (2 ** attempt)) + random.uniform(0, 1.0))
    time.sleep(delay)

# Helper for fuzzy matching
def _generate_fuzzy_queries(term: str) -> list[str]:
    """Generates a list of query variations for fuzzy matching."""
    variations = {term} # Start with the original term

    # Remove all spaces
    no_spaces = term.replace(' ', '')
    if no_spaces != term:
        variations.add(no_spaces)

    # Replace spaces with hyphens
    if ' ' in term:
        hyphenated = term.replace(' ', '-')
        if hyphenated != term:
            variations.add(hyphenated)

    # Remove all hyphens
    no_hyphens = term.replace('-', '')
    if no_hyphens != term:
        variations.add(no_hyphens)

    return sorted(list(variations))

# -----------------------------
# Term input / output
# -----------------------------

def read_terms_from_file(filename=input_terms_file):
    """Reads terms from a text file, one term per line."""
    terms = []
    try:
        with open(filename, 'r') as f:
            for line in f:
                term = line.strip()
                if term:
                    terms.append(term)
    except FileNotFoundError:
        logging.error(f"Input file not found: '{filename}'.")
        # No re-raising, as the main block handles default terms if the list is empty
    return terms

def write_term_metrics_tsv(metrics, out_tsv_path):
    """
    Write per-term citation metrics to a TSV.
    Expected keys per entry:
      Term, Total Mentions, H-index, G-index, i10-index
    """
    if not metrics:
        logging.info("No term metrics to write.")
        return

    fieldnames = ["Term", "Total Mentions", "H-index", "G-index", "i10-index"]

    try:
        with open(out_tsv_path, "w", newline="\n", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            w.writeheader()
            for row in metrics:
                w.writerow(row)
        logging.info(f"Wrote term metrics table to {out_tsv_path}")
    except IOError as e:
        logging.error(f"Error writing term metrics to {out_tsv_path}: {e}")

# -----------------------------
# PubMed Central search (unchanged)
# -----------------------------

def search_pmc(query, window_length_years, end_offset_years, retmax, additional_query_term):
  # Note:Searches PubMed Central for articles matching the query within a specified date window.

  current_year = datetime.now().year
  end_year_query = current_year - end_offset_years
  start_year_query = end_year_query - window_length_years + 1

  date_since = datetime(start_year_query, 1, 1).strftime("%Y/%m/%d")
  date_until = datetime(end_year_query, 12, 31).strftime("%Y/%m/%d")

  # Generate fuzzy query variations
  query_variations = _generate_fuzzy_queries(query)
  # Wrap each variation in quotes for phrase searching and join with OR
  fuzzy_query_part = "(" + " OR ".join([f'"{v}"' for v in query_variations]) + ")"

  full_query = f'{fuzzy_query_part} AND ({date_since}[pdat] : {date_until}[pdat])'

  if additional_query_term:
      full_query += f' AND {additional_query_term}'

  logging.info(f"PMC search query: {full_query}") # Log the actual query being used

  try:
    handle = Entrez.esearch(db="pmc", term=full_query, retmax=retmax)
    record = Entrez.read(handle)
    handle.close()
    return record.get("IdList", [])
  except Exception as e:
    logging.error(f"Error searching PMC for query '{full_query}': {e}")
    return []

def pmc_to_pmid_idconv(pmc_ids, tool, email, batch_size, retries):
    """
    Convert PMC IDs (either 'PMCxxxx' or numeric strings) -> PMIDs using NCBI PMC ID Converter API.
    Returns a de-duplicated list of PMID strings.
    """
    if not pmc_ids:
        return []

    # normalize: ensure numeric strings, since your PMC esearch returns numeric UIDs
    pmc_numeric = []
    for x in pmc_ids:
        s = str(x).strip()
        if s.upper().startswith("PMC"):
            s = s[3:]
        if s.isdigit():
            pmc_numeric.append(s)

    pmids = set()

    for i in range(0, len(pmc_numeric), batch_size):
        batch = pmc_numeric[i:i+batch_size]
        ids_param = ",".join(batch)

        # idtype=pmcid lets us pass numeric PMCIDs without the 'PMC' prefix
        url = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
        params = {
            "ids": ids_param,
            "idtype": "pmcid",
            "format": "json",
            "tool": tool,
            "email": email,
            "versions": "no",
        }

        for attempt in range(retries):
          try:
            r = requests.get(url, params=params, timeout=maxdelay)
            r.raise_for_status()
            payload = r.json()
            for rec in payload.get("records", []):
              if rec.get("pmid"):pmids.add(str(rec["pmid"]))
            break # success

          except requests.exceptions.RequestException as e:
            logging.warning(f"NCBI ID conversion failed (attempt {attempt+1}): {e}; backing off...")
            backoff_sleep(attempt, timedelay, maxdelay)
            if attempt == retries - 1:
              logging.error(f"Final attempt failed for NCBI ID conversion: {e}")
              raise # Re-raise if all retries fail
          except Exception as e:
            logging.error(f"Unexpected error during NCBI ID conversion for batch {ids_param}: {e}")
            if attempt == retries - 1:
                raise
          time.sleep(timedelay)

    return sorted(pmids)

def dump_term_citations_tsv(
    term: str,
    pmids: list[str],
    out_tsv_path: str,
    email: str,
    tool: str,
    batch_size: int,
    request_delay: float
):
    """
    Write a TSV with columns to a single file:
      Term, PubMedID, Title, AuthorList, Institution, etc.

    Notes:
    - Institution is taken from PubMed Author->AffiliationInfo entries (free text).
      Multiple affiliations are joined with ' | ' and de-duplicated.
    - AuthorList is "LastName Initials" joined with '; '.
    - Requires Entrez.email to be set
    """

    # Defensive: if no PMIDs, still write header-only TSV.
    fieldnames = ["Term", "PubMedID", "Year", "Title", "Journal", "AuthorList", "Institution"]

    file_exists = os.path.exists(out_tsv_path)

    try:
        with open(out_tsv_path, "a", encoding="utf-8", newline="\n") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")

            # write header only once
            if not file_exists:
                w.writeheader()

            if not pmids:
                return

            # Fetch in batches to reduce calls
            for batch in chunks([str(p) for p in pmids], batch_size):
                handle = Entrez.efetch(db="pubmed", id=",".join(batch), retmode="xml")
                records = Entrez.read(handle)
                handle.close()

                for pubmed_article in records.get("PubmedArticle", []):
                    medline = pubmed_article.get("MedlineCitation", {})
                    article = medline.get("Article", {})

                    pmid = str(medline.get("PMID", ""))

                    # Title (may be a StringElement)
                    title = article.get("ArticleTitle", "")
                    title = str(title) if title is not None else ""

                    # get journal name
                    journal = ""
                    j = article.get("Journal", {})
                    if j and "Title" in j:
                        journal = str(j["Title"])

                    # get year
                    year = ""
                    try:
                      pubdate = article.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {})
                      if "Year" in pubdate and pubdate["Year"]:
                        year = str(pubdate["Year"])
                      elif "MedlineDate" in pubdate and pubdate["MedlineDate"]:
                        m = re.search(r"\b(19|20)\d{2}\b", str(pubdate["MedlineDate"]))
                        year = int(m.group(0)) if m else ""
                      else:
                        ad = article.get("ArticleDate", [])
                        if ad and ad[0].get("Year"):
                          year = str(ad[0]["Year"])
                    except Exception:
                      logging.warning(f"Could not parse publication year for PMID {pmid}.")

                    # Authors and affiliations
                    authors_out = []
                    affs = []

                    for author in article.get("AuthorList", []):
                        last = author.get("LastName")
                        initials = author.get("Initials")
                        collective = author.get("CollectiveName")

                        if collective:
                            authors_out.append(str(collective))
                        elif last:
                            if initials:
                                authors_out.append(f"{last} {initials}")
                            else:
                                authors_out.append(str(last))

                        for aff in author.get("AffiliationInfo", []):
                            aff_text = aff.get("Affiliation")
                            if aff_text:
                                affs.append(str(aff_text))

                    # De-duplicate affiliations but preserve order
                    seen = set()
                    affs_dedup = []
                    for a in affs:
                        if a not in seen:
                            affs_dedup.append(a)
                            seen.add(a)

                    w.writerow({
                        "Term": term, # Changed from SKU
                        "PubMedID": pmid,
                        "Year": year,
                        "Title": smart_tc(title.replace("\t", " ").replace("\n", " ").strip()),
                        "Journal": smart_tc(journal.replace("\t", " ").replace("\n", " ").strip()),
                        "AuthorList": smart_tc("; ".join(authors_out).replace("\t", " ").replace("\n", " ").strip()),
                        "Institution": smart_tc(" | ".join(affs_dedup).replace("\t", " ").replace("\n", " ").strip())
                    })

                # polite pacing between requests (esp. for large jobs)
                time.sleep(timedelay)
    except IOError as e:
        logging.error(f"Error writing citation data for term {term} to {out_tsv_path}: {e}")


# -----------------------------
# iCite citation counts (NEW)
# -----------------------------

def get_citation_counts_icite(pubmed_id_list):
     # Retrieves citation counts for each PMID using NIH iCite.
     # iCite supports PMID lookups and returns citation_counts

    if not pubmed_id_list:
        return []

    # iCite coverage is 1980-present, and RCR availability is 1980-2026; older articles may be missing.
    citation_counts = {pmid: 0 for pmid in pubmed_id_list}

    # Request only what we need (citation_count). Field naming differs between legacy vs non-legacy,
    # but citation_count is present in the support doc example.
    params_template = {
        "format": "json",
        "fl": "pmid,citation_count",
        # optionally: "legacy": "false"
    }

    for batch in chunks(pubmed_id_list, ICITE_BATCH_SIZE):
        pmids_csv = ",".join(batch)
        params = dict(params_template)
        params["pmids"] = pmids_csv

        for attempt in range(Entrez.max_tries):
            try:
                r = requests.get(ICITE_BASE, params=params, timeout=maxdelay)
                r.raise_for_status()
                payload = r.json()

                # payload["data"] is a list of records
                for rec in payload.get("data", []):
                  pmid = str(rec.get("pmid"))
                  # iCite uses 'citation_count' in examples.
                  c = rec.get("citation_count")
                  if c is None:
                    citation_counts[pmid] = 0
                  else:
                    citation_counts[pmid] = int(c)
                break  # success

            except requests.exceptions.RequestException as e:
                logging.warning(f"iCite batch failure (attempt {attempt+1}): {e}; backing off...")
                backoff_sleep(attempt, timedelay, maxdelay)
                if attempt == Entrez.max_tries - 1:
                    logging.error(f"Final attempt failed for iCite batch: {e}")
                    raise # Re-raise if all retries fail
            except Exception as e:
                logging.error(f"Unexpected error during iCite API call for batch {pmids_csv}: {e}")
                if attempt == Entrez.max_tries - 1:
                    raise

        # be polite to the service; iCite doesn't publish a strict rps spec in these docs,
        # but batching reduces overall load.
        time.sleep(timedelay)

    return [{"PMID": pmid, "CitedByCount": citation_counts.get(pmid, 0)} for pmid in pubmed_id_list]

# -----------------------------
# Index calculations
# -----------------------------

def compute_h_g_i10(citation_data):
    # Given list of {'PMID','CitedByCount'}, compute H, G, i10.
    if not citation_data:
        return 0, 0, 0

    citations_sorted = sorted(citation_data, key=lambda x: x["CitedByCount"], reverse=True)

    # H-index
    h_index = 0
    for i, pub_info in enumerate(citations_sorted):
        if pub_info["CitedByCount"] >= (i + 1):
            h_index = i + 1
        else:
            break

    # G-index
    g_index = 0
    cumulative = 0
    for i, pub_info in enumerate(citations_sorted):
        cumulative += pub_info["CitedByCount"]
        if cumulative >= (i + 1) ** 2:
            g_index = i + 1
        else:
            break

    # i10-index
    i10_index = sum(1 for pub_info in citations_sorted if pub_info["CitedByCount"] >= 10)

    return h_index, g_index, i10_index

# -----------------------------
# Main
# -----------------------------

if __name__ == "__main__":
    start_time = datetime.now()
    logging.info(f"Script start: {start_time}")

    # To run this in Colab, you might need to mock or remove sys.argv parsing
    # or explicitly pass args=["..."] to parser.parse_args()

    terms = read_terms_from_file(args.input_file)

    # Add default terms if the input file was not found or was empty
    if not terms:
        logging.info(f"No terms found in '{args.input_file}'. Using default terms!")
        terms = ['ATCC CRL-1585', 'VR-1490']

    logging.info(f"Terms to process: {terms}")

    all_terms_metrics = []

    for term in terms:
        logging.info(f"\n--- Processing Term: {term} ---")

        # 1) Search PMC for the term
        # The query is now constructed with fuzzy variations inside search_pmc
        pmc_ids = search_pmc(term, window_length_years, end_offset_years, retmax, additional_query_term)
        logging.info(f"Found {len(pmc_ids)} PMC IDs for {term} from a {window_length_years}-year window ending {end_offset_years} year(s) ago.")

        # 2) Convert PMCIDs -> PMIDs
        pubmed_ids = pmc_to_pmid_idconv(pmc_ids, Entrez.tool, Entrez.email, IDCONV_BATCH_SIZE, Entrez.max_tries)

        logging.info(f"Converted {len(pmc_ids)} PMCIDs to {len(pubmed_ids)} unique PMIDs for {term}.")
        total_mentions = len(pubmed_ids)

        if not pubmed_ids:
            logging.info(f"No PMIDs found for Term {term}.")
            all_terms_metrics.append({"Term": term, "Total Mentions": 0, "H-index": 0, "G-index": 0, "i10-index": 0})
            continue

        # dump citations for this term to TSV
        dump_term_citations_tsv(
            term=term,
            pmids=pubmed_ids,
            out_tsv_path=out_citations_tsv,
            email=Entrez.email,
            tool=Entrez.tool,
            batch_size=retmax,
            request_delay=timedelay
        )

        # 3) Get citation counts via iCite
        citation_data = get_citation_counts_icite(pubmed_ids)

        # 4) Compute metrics
        h, g, i10 = compute_h_g_i10(citation_data)

        logging.info(f"\nMetrics for Term '{term}':")
        logging.info(f"  H-index: {h}")
        logging.info(f"  G-index: {g}")
        logging.info(f"  i10-index: {i10}")

        all_terms_metrics.append({"Term": term, "Total Mentions": total_mentions, "H-index": h, "G-index": g, "i10-index": i10})

    logging.info("\n--- Summary of All Terms ---")
    for entry in all_terms_metrics:
      logging.info(entry)

    write_term_metrics_tsv(all_terms_metrics, out_tsv_path=out_metrics_tsv)

    end_time = datetime.now()
    # Ensure start_clock is defined, in case the script is run directly outside __main__
    if 'start_clock' not in locals():
        start_clock = time.time() # Define if not already defined
    elapsed_sec = time.time() - start_clock

    logging.info(f"Script end: {end_time}")
    logging.info(f"Elapsed time: {elapsed_sec:.1f} seconds ({elapsed_sec/60:.2f} minutes)")
    
