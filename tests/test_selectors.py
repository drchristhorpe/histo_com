import pytest

from histo_com.selectors import (
    DomainRef,
    ResidueRef,
    SelectorError,
    parse_domains,
    parse_residues,
)


def test_parse_domains_single_chain():
    assert parse_domains("P") == [DomainRef(chain="P")]


def test_parse_domains_chain_list():
    assert parse_domains("A,B") == [DomainRef(chain="A"), DomainRef(chain="B")]


def test_parse_domains_chain_range():
    assert parse_domains("L:1-180") == [DomainRef(chain="L", start=1, end=180)]


def test_parse_domains_bare_range():
    assert parse_domains("1-180") == [DomainRef(chain=None, start=1, end=180)]


def test_parse_domains_from_iterable():
    assert parse_domains(["A", "B"]) == [DomainRef(chain="A"), DomainRef(chain="B")]


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
