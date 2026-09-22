# Repository guidance for automated contributors

## Scope

These instructions apply throughout the repository unless a more specific `AGENTS.md` is added in a subdirectory.

## Working rules

- Preserve unrelated work already present in the checkout.
- Consult the "Where information belongs" list below before creating documentation. (This project has not adopted the template's full `docs/DOCUMENTATION_OWNERSHIP.md` matrix — see `docs/GETTING_STARTED.md`'s proportionality guidance in `atcc-dev-template` for why.)
- Follow `docs/DATA_GOVERNANCE.md` before writing data, metadata, settings, logs, reports, or generated files — this applies to GitHub issues, pull-request text, and comments too, not only files being committed.
- Follow `docs/THIRD_PARTY.md` before adding a tool, container, reference, database, dataset, or code taken from outside the project. Do not copy outside code into the repository without recording its source and terms.
- Do not change scientific thresholds, reference choices, identifier mappings, or workflow behavior unless explicitly requested.
- Keep examples sanitized; never introduce credentials, internal paths, internal hostnames or server/node names, restricted data, or production outputs — including when pasting a log excerpt or describing a computing environment for a validation-history entry, issue, or commit message. Before publishing any of those, check the text for anything that identifies specific internal infrastructure rather than describing it generically.
- Do not claim portability, validation, reproducibility, or release readiness without recorded evidence.
- State skipped checks as plainly as checks that passed.
- Treat an actual run, log, or tool output as stronger evidence than remembered tool behavior, documentation, or an earlier assertion in conversation.

## Writing for scientific teams

- Write for bioinformaticians, biologists, scientists, and workflow users first.
- Prefer familiar scientific and operational language, concrete examples, and direct questions over specialized software-engineering terminology.
- When a technical term is necessary, explain it in plain language the first time it appears. Do not sacrifice accuracy to avoid a useful term.
- Lead with the document's practical purpose: what question it answers, when someone should use it, and what decision or work it supports.
- Keep guidance proportionate. Do not add fields, process, or separate documents unless they make the work easier to understand, review, reproduce, operate, or maintain.
- Avoid or explain terms such as artifact, contract, orchestration, boundary, idempotency, side effect, fixture, and provenance when plainer wording would communicate the same meaning.

## Sentence-level writing style

Write documentation, commit messages, pull request text, and issue comments as direct statements of fact. Report what is true and what was observed. Do not argue a case, build toward a conclusion, or promote a recommendation.

Do not use these constructions:

- Antithesis, such as "X, not Y", "it isn't A, it's B", or "this is A rather than B". Write what is true and continue.
- Short fragments used for emphasis, and parallel clauses stacked for rhythm.
- A bold one-line pronouncement to open a section.
- Intensifiers that carry no information, such as actually, genuinely, properly, real, and truly.
- Business idiom, such as "pays for itself", "worth the investment", or "high value per effort".
- Developer and startup slang, such as dogfooding, foot-gun, load-bearing, blast radius, surface area, happy path, source of truth, and north star.

Those are examples, and the list will never be complete. The general rule is that a phrase drawn from software-industry or startup culture does not belong in writing aimed at scientists. Write the plain meaning instead.

Name the concrete thing instead of describing the shape of an argument. A sentence such as "the failure modes are asymmetric" tells a reader nothing about what happens; state the input, the result, and the consequence. Before keeping a summary sentence, check whether someone who does not already know the answer could work out what it refers to.

## Markdown line breaks

Do not hard-wrap markdown at any column width. Write each paragraph and each list item as one long line and let the renderer wrap it. This applies to repository files, GitHub issue bodies, issue comments, and pull request descriptions.

Do not insert manual line breaks with `<br>` or with two trailing spaces. When several short labeled lines should render separately, such as the audience and purpose block at the top of each document, write them as a bullet list.

Table rows stay on one line, and text inside fenced code blocks keeps whatever line structure it needs.

## Where information belongs

- `README.md`: orientation and project status.
- `CHANGELOG.md`: material project changes and verification limits.
- `CONTRIBUTING.md`: contribution and review expectations.
- `docs/DESIGN.md`: project rationale, goals, and proposed solution.
- `docs/REQUIREMENTS.md`: acceptable behavior, scope, and completion criteria.
- `docs/CODE_STYLE.md`: readable Python conventions and the role of automated style checks.
- `docs/DATA_GOVERNANCE.md`: rules for data, settings, logs, and results.
- `docs/THIRD_PARTY.md`: terms for outside tools, references, data, and code.
- `docs/TESTING.md`: repeatable quality procedures.
- `docs/history/VALIDATION_HISTORY.md`: dated observed outcomes.
