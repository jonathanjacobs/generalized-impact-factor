"""Backoff helper for retrying NCBI/iCite requests after a failed attempt."""

from __future__ import annotations

import random
import time


def backoff_sleep(attempt: int, base: float, cap: float) -> None:
    """Sleep with exponential backoff and jitter before retrying a failed request.

    Call this only after a first attempt has failed; ``attempt`` is the zero-based retry count.
    """
    delay = min(cap, (base * (2**attempt)) + random.uniform(0, 1.0))
    time.sleep(delay)
