# Using other people's software, data, and code

- **Audience:** Anyone adding a tool, container, reference, database, dataset, or borrowed code to a project.
- **Use this when:** Choosing a dependency, adapting published code, or preparing work that will leave the team.
- **Update this when:** The project's approved tools or references change, or their terms change.
- **Do not update for:** One workflow run, or questions about data sensitivity and storage.

This document is a starting point and is not legal advice. Follow all applicable ATCC, legal, procurement, and information-security rules, and ask when the answer is not clear. When two rules differ, follow the stricter one.

[Working with data safely](DATA_GOVERNANCE.md) covers whether material is sensitive and where it may be stored. This document covers a different question: whether the project is permitted to use something, and what it owes in return.

## The basic rule

Check the terms while you are still choosing the tool or reference, when it is easy to pick something else. A tool that has to be removed after an analysis has been built on it costs far more than one ruled out in an afternoon.

Many widely used bioinformatics tools and reference databases are free to download but are licensed for academic or non-commercial research only, and commercial use requires a separate license. ATCC is a commercial organization, so the fact that something is free to obtain does not answer whether the project may use it. Popularity, publication, availability from a public package channel, and use by other groups in the field do not answer it either.

## Four questions to ask

- **Use:** may we run it for our purpose, in our setting?
- **Modify:** may we change it, and does changing it create obligations?
- **Share:** may we pass it on — in a container image, a workflow, a deliverable, a publication, or to a customer?
- **Credit:** what attribution, citation, or notice is required?

If any answer is unclear, record it as unresolved and ask before building on it.

## What to check

| What is being added | Watch for |
| --- | --- |
| Software tool or package | Academic or non-commercial-only licenses; per-seat, per-site, or per-server terms |
| Container or environment image | Everything inside it is redistributed with it, under its own terms |
| Reference genome, annotation, or database | Use restrictions, redistribution limits, required citation, version-specific terms |
| Public dataset | Data use agreements, controlled access, publication and attribution conditions |
| Code from a repository, paper, or forum | License compatibility and required attribution; a public repository with no license grants no rights |
| Another team's or a vendor's workflow | Internal approval, support expectations, and limits on sharing it onward |

Copyleft licenses such as GPL and AGPL are usually straightforward for internal use but create obligations as soon as the work is shared outside the team, including inside a container image. Identify that situation early.

## Studied compared with adapted

Reading someone's method or code to understand their approach is a different situation from copying it into this project, and the two are recorded differently.

- **Studied:** note the reference in the project's design, architecture, or investigation record so the influence stays traceable, and reimplement the approach independently.
- **Adapted or copied:** record it in the register below, with its source, license, and required attribution, before the work is shared.

Public availability does not grant redistribution rights.

## Register of third-party material

| Item and version | Source | Terms, permitted use, and required attribution | Checked by and when |
| --- | --- | --- | --- |
| Biopython 1.78 (`Bio.Entrez`) | https://biopython.org/ | Biopython License Agreement (permissive: free use, copy, modify, distribute, including commercial use, no fee; must keep the copyright/permission notice; no redistribution obligations beyond that). No attribution required beyond the notice. | Jonathan Jacobs, 2026-09-22, from the package's bundled `LICENSE.rst` |
| NCBI E-utilities / PMC ID Converter API (esearch, efetch, idconv) | https://www.ncbi.nlm.nih.gov/home/develop/api/ | Free to use with required identification (`email`, `tool` params — already set in `impact_lookup.py`). Rate limit is 3 requests/second without an API key, up to 10/second with one (`NCBI_API_KEY`, already supported). No redistribution restriction found for retrieved bibliographic metadata (titles, authors, journal, affiliations); this tool does not fetch full text. Commercial use is not restricted by the E-utilities usage policy itself. | Jonathan Jacobs, 2026-09-22, from NCBI's published E-utilities usage guidelines |
| NIH iCite API (citation counts) | https://icite.od.nih.gov/api | NIH Open Citation Collection data is released under the CC0 public domain dedication — no restriction on use, modification, or redistribution, and no attribution legally required (citing iCite as the source is good practice but not mandatory). | Jonathan Jacobs, 2026-09-22, from iCite's published documentation |
| `titlecase` (Python package) | https://github.com/ppannuto/python-titlecase | MIT License — permissive, commercial use permitted, requires keeping the license/copyright notice if redistributed. | Jonathan Jacobs, 2026-09-22, from PyPI package metadata |
| `requests` (Python package) | https://requests.readthedocs.io | Apache License 2.0 — permissive, commercial use permitted, requires keeping the license notice and stating changes if the code itself is modified and redistributed (this project does not modify it). | Jonathan Jacobs, 2026-09-22, from PyPI package metadata |

Unresolved: whether an internal ATCC report built from PMC/PubMed-derived metadata and iCite citation counts needs any additional disclosure or citation convention beyond what's noted above — ask whoever reviews the report before it goes out.

## Decisions for each project

- Who reviews terms before a dependency is adopted? **TBD**
- Which tools, references, and databases are already approved for this work? Biopython, NCBI E-utilities/PMC ID Converter, NIH iCite, `titlecase`, and `requests`, per the register above.
- Is the result intended for internal use, a customer deliverable, a product, or publication? Internal use supporting a scientific or business conclusion (see `docs/DESIGN.md`). Revisit this register if that changes — the last question changes the answers to most of the others.
