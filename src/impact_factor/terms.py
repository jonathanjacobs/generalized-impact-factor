"""Reading search terms and generating fuzzy-matching query variations for them."""

from __future__ import annotations

import logging


def generate_fuzzy_queries(term: str) -> list[str]:
    """Return sorted query variations of ``term`` with spaces/hyphens added or removed, for fuzzy PMC matching."""
    variations = {term}

    no_spaces = term.replace(" ", "")
    if no_spaces != term:
        variations.add(no_spaces)

    if " " in term:
        hyphenated = term.replace(" ", "-")
        if hyphenated != term:
            variations.add(hyphenated)

    no_hyphens = term.replace("-", "")
    if no_hyphens != term:
        variations.add(no_hyphens)

    return sorted(variations)


def read_terms_from_file(filename: str) -> list[str]:
    """Read one term per line from ``filename``, skipping blank lines. Returns an empty list if the file is missing."""
    terms = []
    try:
        with open(filename, "r") as f:
            for line in f:
                term = line.strip()
                if term:
                    terms.append(term)
    except FileNotFoundError:
        logging.error(f"Input file not found: '{filename}'.")
    return terms
