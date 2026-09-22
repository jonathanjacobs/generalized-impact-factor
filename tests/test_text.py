"""Tests for smart_tc (acronym/alphanumeric-preserving titlecase) and chunks."""

from impact_factor import chunks, smart_tc


def test_empty_string_returns_empty_string():
    assert smart_tc("") == ""


def test_preserves_uppercase_acronyms():
    assert smart_tc("DNA extraction protocol") == "DNA Extraction Protocol"


def test_preserves_multiple_acronyms_and_alphanumeric_codes():
    assert (
        smart_tc("expression of GAPDH in HEK293 cells")
        == "Expression of GAPDH in HEK293 Cells"
    )


def test_preserves_alphanumeric_catalog_number():
    assert (
        smart_tc("ATCC CRL-1585 cell line characterization")
        == "ATCC CRL-1585 Cell Line Characterization"
    )


def test_chunks_splits_into_expected_sizes():
    assert list(chunks([1, 2, 3, 4, 5], 2)) == [[1, 2], [3, 4], [5]]


def test_chunks_empty_list_yields_nothing():
    assert list(chunks([], 3)) == []
