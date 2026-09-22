"""Tests for compute_h_g_i10 (H-index, G-index, i10-index)."""

from impact_factor import compute_h_g_i10


def make_citations(counts):
    return [{"PMID": str(i), "CitedByCount": c} for i, c in enumerate(counts)]


def test_empty_citation_list_returns_zeros():
    assert compute_h_g_i10([]) == (0, 0, 0)


def test_all_zero_citations_returns_zeros():
    assert compute_h_g_i10(make_citations([0, 0, 0])) == (0, 0, 0)


def test_known_citation_counts():
    # Sorted desc: 10, 8, 5, 4, 3
    # H: 10>=1, 8>=2, 5>=3, 4>=4, 3>=5(no) -> h=4
    # G: cum 10>=1, 18>=4, 23>=9, 27>=16, 30>=25 -> g=5
    # i10: only the 10-citation paper -> i10=1
    assert compute_h_g_i10(make_citations([10, 8, 5, 4, 3])) == (4, 5, 1)


def test_single_highly_cited_paper():
    assert compute_h_g_i10(make_citations([1000])) == (1, 1, 1)


def test_input_order_does_not_matter():
    # Sorted desc: 10, 3, 1
    # H: 10>=1, 3>=2, 1>=3(no) -> h=2
    # G: cum 10>=1, 13>=4, 14>=9 -> g=3
    # i10: only the 10-citation paper -> i10=1
    unsorted = make_citations([3, 10, 1])
    assert compute_h_g_i10(unsorted) == (2, 3, 1)
