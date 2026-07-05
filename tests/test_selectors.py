import pytest

from histo_com.selectors import (
    DomainPart,
    DomainRef,
    ResidueRef,
    SelectorError,
    format_domain,
    format_residue,
    parse_domains,
    parse_residues,
)


def test_parse_domains_single_chain():
    assert parse_domains("P") == [DomainRef(parts=(DomainPart(chain="P"),))]


def test_parse_domains_chain_list_is_one_combined_domain():
    # A comma-separated string is ONE domain: its parts are combined into
    # a single centre of mass, not treated as separate domains.
    assert parse_domains("A,B") == [
        DomainRef(parts=(DomainPart(chain="A"), DomainPart(chain="B")))
    ]


def test_parse_domains_chain_range():
    assert parse_domains("L:1-180") == [DomainRef(parts=(DomainPart(chain="L", start=1, end=180),))]


def test_parse_domains_bare_range():
    assert parse_domains("1-180") == [DomainRef(parts=(DomainPart(chain=None, start=1, end=180),))]


def test_parse_domains_from_iterable_gives_separate_domains():
    # An iterable of separate elements describes multiple, separate
    # domains — unlike a single comma-joined string.
    assert parse_domains(["A", "B"]) == [
        DomainRef(parts=(DomainPart(chain="A"),)),
        DomainRef(parts=(DomainPart(chain="B"),)),
    ]


def test_parse_domains_iterable_element_can_be_compound():
    assert parse_domains(["A,B", "C"]) == [
        DomainRef(parts=(DomainPart(chain="A"), DomainPart(chain="B"))),
        DomainRef(parts=(DomainPart(chain="C"),)),
    ]


def test_parse_domains_invalid_range():
    with pytest.raises(SelectorError):
        parse_domains("A:180-1")


def test_parse_domains_empty():
    with pytest.raises(SelectorError):
        parse_domains("")


def test_parse_residues_single():
    assert parse_residues("5") == [ResidueRef(chain=None, resseq=5)]


def test_parse_residues_range_expands():
    refs = parse_residues("1-9")
    assert refs == [ResidueRef(chain=None, resseq=i) for i in range(1, 10)]


def test_parse_residues_mixed_list():
    refs = parse_residues("A:1-3,B:12")
    assert refs == [
        ResidueRef(chain="A", resseq=1),
        ResidueRef(chain="A", resseq=2),
        ResidueRef(chain="A", resseq=3),
        ResidueRef(chain="B", resseq=12),
    ]


def test_parse_residues_from_int_iterable():
    refs = parse_residues(range(1, 4))
    assert refs == [
        ResidueRef(chain=None, resseq=1),
        ResidueRef(chain=None, resseq=2),
        ResidueRef(chain=None, resseq=3),
    ]


def test_parse_residues_from_tuple():
    assert parse_residues([("A", 5)]) == [ResidueRef(chain="A", resseq=5)]


@pytest.mark.parametrize(
    "ref,expected",
    [
        (DomainRef(parts=(DomainPart(chain="P"),)), "P"),
        (DomainRef(parts=(DomainPart(chain="A", start=1, end=180),)), "A:1-180"),
        (DomainRef(parts=(DomainPart(chain=None, start=1, end=180),)), "1-180"),
        (DomainRef(parts=(DomainPart(chain="A", start=12, end=12),)), "A:12"),
        (
            DomainRef(parts=(DomainPart(chain="A"), DomainPart(chain="B"))),
            "A,B",
        ),
    ],
)
def test_format_domain(ref, expected):
    assert format_domain(ref) == expected


@pytest.mark.parametrize(
    "ref,expected",
    [
        (ResidueRef(chain=None, resseq=5), "5"),
        (ResidueRef(chain="A", resseq=12), "A:12"),
    ],
)
def test_format_residue(ref, expected):
    assert format_residue(ref) == expected


def test_format_domain_round_trips_through_parse():
    for spec in ["P", "A,B", "L:1-180", "1-180", "A:12"]:
        (ref,) = parse_domains(spec)
        assert format_domain(ref) == spec
