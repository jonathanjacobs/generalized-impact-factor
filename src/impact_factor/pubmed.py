"""PubMed Central search, PMC-to-PMID conversion, and citation metadata retrieval."""

from __future__ import annotations

import csv
import logging
import os
import re
import time
from datetime import datetime

import requests
from Bio import Entrez

from .retry import backoff_sleep
from .terms import generate_fuzzy_queries
from .text import chunks, smart_tc


def search_pmc(
    query: str,
    window_length_years: int,
    end_offset_years: int,
    retmax: int,
    additional_query_term: str | None,
) -> list[str]:
    """Search PubMed Central for ``query`` (with fuzzy variations) within a date window, returning PMC IDs.

    Requires ``Entrez.email`` (and optionally ``Entrez.api_key``) to already be set by the caller.
    """
    current_year = datetime.now().year
    end_year_query = current_year - end_offset_years
    start_year_query = end_year_query - window_length_years + 1

    date_since = datetime(start_year_query, 1, 1).strftime("%Y/%m/%d")
    date_until = datetime(end_year_query, 12, 31).strftime("%Y/%m/%d")

    query_variations = generate_fuzzy_queries(query)
    fuzzy_query_part = (
        "(" + " OR ".join([f'"{v}"' for v in query_variations]) + ")"
    )

    full_query = (
        f"{fuzzy_query_part} AND ({date_since}[pdat] : {date_until}[pdat])"
    )

    if additional_query_term:
        full_query += f" AND {additional_query_term}"

    logging.info(f"PMC search query: {full_query}")

    try:
        handle = Entrez.esearch(db="pmc", term=full_query, retmax=retmax)
        record = Entrez.read(handle)
        handle.close()
        return record.get("IdList", [])
    except Exception as e:
        logging.error(f"Error searching PMC for query '{full_query}': {e}")
        return []


def pmc_to_pmid_idconv(
    pmc_ids: list[str],
    tool: str,
    email: str,
    batch_size: int,
    retries: int,
    timedelay: float,
    maxdelay: float,
) -> list[str]:
    """Convert PMC IDs to a de-duplicated, sorted list of PMID strings using the NCBI PMC ID Converter API."""
    if not pmc_ids:
        return []

    pmc_numeric = []
    for x in pmc_ids:
        s = str(x).strip()
        if s.upper().startswith("PMC"):
            s = s[3:]
        if s.isdigit():
            pmc_numeric.append(s)

    pmids = set()

    for batch in chunks(pmc_numeric, batch_size):
        ids_param = ",".join(batch)

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
                    if rec.get("pmid"):
                        pmids.add(str(rec["pmid"]))
                break

            except requests.exceptions.RequestException as e:
                logging.warning(
                    f"NCBI ID conversion failed (attempt {attempt + 1}): {e}; backing off..."
                )
                backoff_sleep(attempt, timedelay, maxdelay)
                if attempt == retries - 1:
                    logging.error(
                        f"Final attempt failed for NCBI ID conversion: {e}"
                    )
                    raise
            except Exception as e:
                logging.error(
                    f"Unexpected error during NCBI ID conversion for batch {ids_param}: {e}"
                )
                if attempt == retries - 1:
                    raise
            time.sleep(timedelay)

    return sorted(pmids)


def dump_term_citations_tsv(
    term: str,
    pmids: list[str],
    out_tsv_path: str,
    batch_size: int,
    request_delay: float,
) -> None:
    """Append citation metadata (title, journal, authors, institutions) for ``pmids`` to a per-citation TSV.

    Requires ``Entrez.email`` to already be set by the caller. Writes a header the first time the file
    is created; subsequent calls for other terms append to the same file.
    """
    fieldnames = [
        "Term",
        "PubMedID",
        "Year",
        "Title",
        "Journal",
        "AuthorList",
        "Institution",
    ]

    file_exists = os.path.exists(out_tsv_path)

    try:
        with open(out_tsv_path, "a", encoding="utf-8", newline="\n") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")

            if not file_exists:
                w.writeheader()

            if not pmids:
                return

            for batch in chunks([str(p) for p in pmids], batch_size):
                handle = Entrez.efetch(
                    db="pubmed", id=",".join(batch), retmode="xml"
                )
                records = Entrez.read(handle)
                handle.close()

                for pubmed_article in records.get("PubmedArticle", []):
                    medline = pubmed_article.get("MedlineCitation", {})
                    article = medline.get("Article", {})

                    pmid = str(medline.get("PMID", ""))

                    title = article.get("ArticleTitle", "")
                    title = str(title) if title is not None else ""

                    journal = ""
                    j = article.get("Journal", {})
                    if j and "Title" in j:
                        journal = str(j["Title"])

                    year = ""
                    try:
                        pubdate = (
                            article.get("Journal", {})
                            .get("JournalIssue", {})
                            .get("PubDate", {})
                        )
                        if "Year" in pubdate and pubdate["Year"]:
                            year = str(pubdate["Year"])
                        elif (
                            "MedlineDate" in pubdate and pubdate["MedlineDate"]
                        ):
                            m = re.search(
                                r"\b(19|20)\d{2}\b", str(pubdate["MedlineDate"])
                            )
                            year = int(m.group(0)) if m else ""
                        else:
                            ad = article.get("ArticleDate", [])
                            if ad and ad[0].get("Year"):
                                year = str(ad[0]["Year"])
                    except Exception:
                        logging.warning(
                            f"Could not parse publication year for PMID {pmid}."
                        )

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

                    seen = set()
                    affs_dedup = []
                    for a in affs:
                        if a not in seen:
                            affs_dedup.append(a)
                            seen.add(a)

                    w.writerow(
                        {
                            "Term": term,
                            "PubMedID": pmid,
                            "Year": year,
                            "Title": smart_tc(
                                title.replace("\t", " ")
                                .replace("\n", " ")
                                .strip()
                            ),
                            "Journal": smart_tc(
                                journal.replace("\t", " ")
                                .replace("\n", " ")
                                .strip()
                            ),
                            "AuthorList": smart_tc(
                                "; ".join(authors_out)
                                .replace("\t", " ")
                                .replace("\n", " ")
                                .strip()
                            ),
                            "Institution": smart_tc(
                                " | ".join(affs_dedup)
                                .replace("\t", " ")
                                .replace("\n", " ")
                                .strip()
                            ),
                        }
                    )

                time.sleep(request_delay)
    except IOError as e:
        logging.error(
            f"Error writing citation data for term {term} to {out_tsv_path}: {e}"
        )
