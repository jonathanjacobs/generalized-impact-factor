"""Citation counts from the NIH iCite API."""

from __future__ import annotations

import logging
import time

import requests

from .retry import backoff_sleep
from .text import chunks

ICITE_BASE = "https://icite.od.nih.gov/api/pubs"


def get_citation_counts_icite(
    pubmed_id_list: list[str],
    icite_batch_size: int,
    max_tries: int,
    timedelay: float,
    maxdelay: float,
    icite_base: str = ICITE_BASE,
) -> list[dict]:
    """Fetch citation counts for ``pubmed_id_list`` from NIH iCite, returning ``{"PMID", "CitedByCount"}`` dicts.

    iCite coverage is 1980-present; older PMIDs are returned with a citation count of 0 rather than
    being omitted, so every input PMID always has a corresponding output entry.
    """
    if not pubmed_id_list:
        return []

    citation_counts = {pmid: 0 for pmid in pubmed_id_list}

    params_template = {
        "format": "json",
        "fl": "pmid,citation_count",
    }

    for batch in chunks(pubmed_id_list, icite_batch_size):
        pmids_csv = ",".join(batch)
        params = dict(params_template)
        params["pmids"] = pmids_csv

        for attempt in range(max_tries):
            try:
                r = requests.get(icite_base, params=params, timeout=maxdelay)
                r.raise_for_status()
                payload = r.json()

                for rec in payload.get("data", []):
                    pmid = str(rec.get("pmid"))
                    c = rec.get("citation_count")
                    citation_counts[pmid] = 0 if c is None else int(c)
                break

            except requests.exceptions.RequestException as e:
                logging.warning(
                    f"iCite batch failure (attempt {attempt + 1}): {e}; backing off..."
                )
                backoff_sleep(attempt, timedelay, maxdelay)
                if attempt == max_tries - 1:
                    logging.error(f"Final attempt failed for iCite batch: {e}")
                    raise
            except Exception as e:
                logging.error(
                    f"Unexpected error during iCite API call for batch {pmids_csv}: {e}"
                )
                if attempt == max_tries - 1:
                    raise

        time.sleep(timedelay)

    return [
        {"PMID": pmid, "CitedByCount": citation_counts.get(pmid, 0)}
        for pmid in pubmed_id_list
    ]
