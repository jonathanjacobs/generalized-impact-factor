"""Text helpers for titlecasing scientific text and batching sequences."""

from __future__ import annotations

from titlecase import titlecase


def smart_tc(s: str) -> str:
    """Titlecase ``s`` while preserving ALLCAPS acronyms and alphanumeric codes (e.g. SARS-CoV-2, ATCC CRL-1585)."""
    if not s:
        return ""

    def keep_acronyms(word: str, **kwargs) -> str | None:
        if word.isupper():
            return word
        if any(ch.isdigit() for ch in word) and any(
            ch.isalpha() for ch in word
        ):
            return word
        return None

    return titlecase(s, callback=keep_acronyms)


def chunks(lst: list, n: int):
    """Yield successive ``n``-sized slices of ``lst``."""
    for i in range(0, len(lst), n):
        yield lst[i : i + n]
