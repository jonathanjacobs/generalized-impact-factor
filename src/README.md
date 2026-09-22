# src/impact_factor

The reusable library behind `impact_lookup.py`: fuzzy term matching, PubMed Central search, PMC-to-PMID conversion, citation metadata retrieval, NIH iCite lookups, and H-index/G-index/i10-index calculation.

`impact_lookup.py` at the repository root is a thin CLI that imports this package (it adds `src/` to `sys.path` at startup, so no install step is required) and is the only file that touches `argparse` or global `Entrez` configuration. Other scripts or notebooks can import from `impact_factor` directly the same way, after adding `src/` to their path, without going through the CLI.

Needs the environment in `environments/requirements.txt`. Tested with Python 3.9.19; see `tests/` for what is covered by automated checks and `docs/TESTING.md` for what is not.
