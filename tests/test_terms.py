"""Tests for generate_fuzzy_queries and read_terms_from_file."""

from impact_factor import generate_fuzzy_queries, read_terms_from_file


def test_term_with_space_and_hyphen_generates_all_variations():
    assert generate_fuzzy_queries("ATCC CRL-1585") == [
        "ATCC CRL-1585",
        "ATCC CRL1585",
        "ATCC-CRL-1585",
        "ATCCCRL-1585",
    ]


def test_term_with_only_hyphen():
    assert generate_fuzzy_queries("VR-1490") == ["VR-1490", "VR1490"]


def test_term_with_only_space():
    assert generate_fuzzy_queries("Human Papillomavirus") == [
        "Human Papillomavirus",
        "Human-Papillomavirus",
        "HumanPapillomavirus",
    ]


def test_term_with_neither_space_nor_hyphen_has_one_variation():
    assert generate_fuzzy_queries("SARSCoV2") == ["SARSCoV2"]


def test_read_terms_from_file_skips_blank_lines(tmp_path):
    terms_file = tmp_path / "terms.txt"
    terms_file.write_text("ATCC CRL-1585\n\nVR-1490\n   \nSARS-CoV-2\n")

    assert read_terms_from_file(str(terms_file)) == [
        "ATCC CRL-1585",
        "VR-1490",
        "SARS-CoV-2",
    ]


def test_read_terms_from_file_missing_file_returns_empty_list(tmp_path):
    missing = tmp_path / "does_not_exist.txt"

    assert read_terms_from_file(str(missing)) == []
