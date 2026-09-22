"""Bibliometric indices (H-index, G-index, i10-index) and writing per-term metrics output."""

from __future__ import annotations

import csv
import logging


def compute_h_g_i10(citation_data: list[dict]) -> tuple[int, int, int]:
    """Compute the H-index, G-index, and i10-index from a list of ``{"PMID", "CitedByCount"}`` dicts."""
    if not citation_data:
        return 0, 0, 0

    citations_sorted = sorted(
        citation_data, key=lambda x: x["CitedByCount"], reverse=True
    )

    h_index = 0
    for i, pub_info in enumerate(citations_sorted):
        if pub_info["CitedByCount"] >= (i + 1):
            h_index = i + 1
        else:
            break

    g_index = 0
    cumulative = 0
    for i, pub_info in enumerate(citations_sorted):
        cumulative += pub_info["CitedByCount"]
        if cumulative >= (i + 1) ** 2:
            g_index = i + 1
        else:
            break

    i10_index = sum(
        1 for pub_info in citations_sorted if pub_info["CitedByCount"] >= 10
    )

    return h_index, g_index, i10_index


def write_term_metrics_tsv(metrics: list[dict], out_tsv_path: str) -> None:
    """Write per-term bibliometric metrics to a TSV with columns Term, Total Mentions, H-index, G-index, i10-index."""
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
