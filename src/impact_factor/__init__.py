"""Reusable functions behind impact_lookup.py: term matching, PubMed/iCite lookups, and bibliometric indices."""

from .icite import get_citation_counts_icite
from .metrics import compute_h_g_i10, write_term_metrics_tsv
from .pubmed import dump_term_citations_tsv, pmc_to_pmid_idconv, search_pmc
from .retry import backoff_sleep
from .terms import generate_fuzzy_queries, read_terms_from_file
from .text import chunks, smart_tc

__all__ = [
    "backoff_sleep",
    "chunks",
    "compute_h_g_i10",
    "dump_term_citations_tsv",
    "generate_fuzzy_queries",
    "get_citation_counts_icite",
    "pmc_to_pmid_idconv",
    "read_terms_from_file",
    "search_pmc",
    "smart_tc",
    "write_term_metrics_tsv",
]
