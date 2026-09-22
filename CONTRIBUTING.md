# Contributing

Make each change understandable to another bioinformatician and proportionate to its scientific, data, and operational risk.

## Before changing the repository

1. Read the nearest README and the relevant requirement or work item.
2. Check [data-handling rules](docs/DATA_GOVERNANCE.md) before adding inputs, settings, logs, or results, and [third-party terms](docs/THIRD_PARTY.md) before adding an outside tool, reference, dataset, or borrowed code.
3. Work on a branch and give the change one clear purpose.
4. Use [the code style guide](docs/CODE_STYLE.md) for Python.

## Make the change reviewable

- Explain what changed, why it was needed, and who or what is affected.
- Make changes to scientific assumptions, references, thresholds, formats, or supported behavior explicit.
- Keep local settings and credentials out of code and Git.
- Update the document that owns information changed by the work.
- Report what was checked, what was skipped, what was observed, and important limitations. A successful style check is not scientific validation.
- Request the technical, data, operational, user, or scientific review that matches the change. Routine work does not need every kind of review.

A change is ready when its intended outcome and applicable requirements are met, the evidence and limitations are clear, affected guidance is current, and the appropriate review has occurred.

## Marking a version

Run records and validation history both refer to a version, so the project needs a way to assign one. Use a git tag of the form `vX.Y.Z` and record what changed under that heading in [CHANGELOG.md](CHANGELOG.md). The scheme is usually called semantic versioning.

Standard semantic versioning treats a major increment as a breaking change to a software interface. For a workflow or analysis, choose the increment by what happens to the results instead:

| Increment | When to use it | Effect on earlier results |
| --- | --- | --- |
| `X`, major | A different reference, changed algorithm, or moved threshold | No longer directly comparable; continuing work must re-run or explain the difference |
| `Y`, minor | New capability, output, or option, with existing behavior unchanged | Still valid |
| `Z`, patch | A fix or speed change that does not alter a valid result | Still valid, unless the bug affected them |

Start at `0.1.0`. Below `1.0.0` the project is still changing and its behavior may change without a major increment. Move to `1.0.0` when other people can depend on the work.

A major increment means earlier validation conclusions no longer apply automatically; record what was re-checked in [validation history](docs/history/VALIDATION_HISTORY.md). If a workflow reports its own version at runtime, keep the number in one file that both the code and the release process read.
