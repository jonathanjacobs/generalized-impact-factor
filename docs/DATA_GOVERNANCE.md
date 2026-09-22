# Working with data safely

- **Audience:** Anyone working with data, sample information, settings, logs, or results.
- **Use this when:** Deciding what may be added to Git and where project data or results should be stored.
- **Update this when:** The project's approved data-handling or storage rules change.
- **Do not update for:** One workflow run; record that run in the project's approved location.

This document is a starting point. Follow all applicable ATCC, project, legal, and information-security rules. When two rules differ, follow the stricter one.

## The basic rule

Use Git for code, documentation, and small approved examples. Do not use it as the main storage location for scientific data or workflow results.

## Keep these out of Git

This applies to any public-facing surface tied to this repository, not only committed files: GitHub issues, pull-request descriptions, comments, and CI logs are just as visible as the repository itself once this remote is public, so treat them with the same care as a commit.

- Raw sequencing data, production sample data, large reference files or indexes, working directories, and production results.
- Confidential, controlled, proprietary, regulated, or otherwise unapproved data and sample information.
- Passwords, tokens, private keys, connection details, and other credentials.
- Logs, reports, or screenshots that reveal restricted information or internal locations that should not be shared — this includes internal hostnames, server or node names, IP addresses, and internal domain names (for example, a specific ATCC HPC cluster node name). When a document needs to describe where something ran, describe the computing environment generically (e.g. "an internal ATCC HPC host," "a shared compute node") instead of naming the specific machine, even in a place that feels as routine as a validation-history entry, a GitHub issue, or a log excerpt pasted for debugging.

If restricted material is added accidentally, stop sharing it and follow the appropriate ATCC response process. Deleting the latest copy does not remove it from Git history.

## What is usually appropriate

- Code, documentation, software-environment files, and examples with sensitive values removed.
- Small synthetic, public, or explicitly approved de-identified test data.
- Approved sample or dataset identifiers, reference versions, checksums, and links to results stored in the correct system.

## Settings and run records

Use obvious placeholders in example settings and keep local or secret settings out of Git. For an important shared or production run, record the versions, inputs, settings, computing environment, checks, result location, and review described in `RUN_PROVENANCE.md`.

## Specific to this project

- `impact_lookup.py`'s two generated outputs, `terms_metrics.tsv` and `terms_citations.tsv`, are run results, not code — they are listed in `.gitignore` and must not be committed. `terms_citations.tsv` contains public PubMed metadata (titles, journals, author names, institutional affiliations); that data is public, but it is still a run result, not a repository artifact.
- **Open question about the input term list:** this repository's remote (`github.com/jonathanjacobs/generalized-impact-factor`) is not under the `ATCC-Bioinformatics` GitHub org and appears to be a personal, likely-public account. `terms.txt` (the input file of product/SKU terms) is not currently tracked in Git, and it should stay that way unless someone confirms that the specific terms in it are fine to publish — a list of real ATCC catalog numbers or unreleased product names is internal information even though each individual term, once searched, only touches public PubMed/PMC data. Decide this before adding any real `terms.txt` to the repository, and before deciding whether this repository should move under the ATCC-Bioinformatics org.

## Decisions for each project

- Who will review data handling? **TBD** — likely the project's Technical Lead (see `README.md`); confirm.
- Which data types and storage systems are approved? **TBD** for where `terms_metrics.tsv`/`terms_citations.tsv` outputs should be kept once generated (not Git — see above). Confirm an approved location (e.g. a shared drive, an internal system) before this tool is used to support a report.
- How long should records be kept, and how should the team link to them? **TBD**.
