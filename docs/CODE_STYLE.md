# Writing readable scientific code

- **Audience:** Anyone writing or reviewing Python code in this project.
- **Use this when:** Adding code, reviewing a change, or deciding how a script should be organized.
- **Update this when:** The team agrees to change a coding convention or an automated style check.
- **Do not update for:** A one-time exception; explain that exception near the code or in the pull request.

The goal of a shared style is to make scientific code easier to read, check, reuse, and maintain. A teammate should be able to understand the work without having to learn the original author's personal conventions.

Automated checks can identify many formatting problems and common coding mistakes. They cannot decide whether a biological assumption, reference, threshold, statistical method, or interpretation is scientifically sound.

This project is Python-only today; the R and Bash sections of the source guide are kept out of this copy for that reason. If R or Bash code is added later, restore them from `atcc-dev-template`'s `docs/CODE_STYLE.md`.

## Practices for every language

- Use names that describe the biological, analytical, or operational meaning of a value. Familiar scientific abbreviations are appropriate when their meaning is clear to the intended readers.
- Keep file locations, credentials, system-specific values, and adjustable settings out of the code. Put safe shared examples in `config/` and keep local values in ignored files such as `settings.local.yaml`.
- Make important thresholds, reference versions, identifier mappings, units, coordinate systems, and assumptions visible and documented.
- Prefer small, focused functions or steps. A reader should be able to tell what each one receives, what it produces, and what can go wrong.
- Explain *why* a non-obvious scientific or technical choice was made. Do not add comments that repeat what a line of code already says.
- Cite a paper, method, or other source when variable names or calculations rely on specialized notation that may not be familiar to the team.
- Provide useful error messages that identify the affected input and expected condition without revealing restricted data or credentials.
- Record or control random seeds when randomness affects a result that must be repeatable.
- Add checks for important expected results, meaningful limits, missing data, and malformed inputs when applicable. Follow `docs/TESTING.md` for the broader testing and validation approach.
- Follow the surrounding project's established style when it does not conflict with an agreed rule in this document.

## Python

- Use `snake_case` for files, functions, methods, and variables; use `CapWords` for classes; and use `UPPER_CASE` for constants.
- Use four spaces for indentation and aim for lines no longer than 80 characters. Readability is more important than forcing an awkward line break for a URL, path, or similar indivisible value.
- Keep imports at the top of the file, organize them consistently, and avoid imports whose origin is unclear.
- Do not use mutable values such as lists or dictionaries as default function arguments.
- Use type hints for public functions and for code where the expected type is not obvious. Do not add complicated annotations that make simple scientific code harder to understand. This project targets Python 3.9, so use `from __future__ import annotations` in any module using `X | Y` union syntax, rather than requiring Python 3.10+.
- Use a context manager, such as `with open(...)`, for files and other resources that must be closed.
- Put executable behavior in a `main()` function and protect it with `if __name__ == "__main__":` when a file can be run directly.
- Use docstrings for public functions, classes, modules, and non-obvious behavior. Describe information a caller needs rather than restating the implementation.

## Automated checks

The pull-request workflow (`.github/workflows/code-style-checks.yml`) runs Ruff on Python files when a pull request touches them. It also flags common local-settings filenames that should not be committed; that check recognizes only common filenames, it does not inspect file contents and is not a secret scanner.

The workflow reports a failed check but is not a required condition for merging — repository owners make that decision through GitHub branch protection. During trial use, correct a finding when appropriate or explain why an exception is reasonable.

Run the same checks locally when the tools are available:

```text
ruff check .
ruff format --check .
```

Formatting and code-quality checks are supporting evidence only. They do not replace the project checks, representative-user testing, scientific review, or validation described in `docs/TESTING.md`.

## Sources

This guide is adapted from `atcc-dev-template`'s `docs/CODE_STYLE.md`, which was itself informed by the [Google Python Style Guide](https://google.github.io/styleguide/pyguide.html) (Creative Commons Attribution 3.0 license) and condensed for ATCC bioinformatics work.
