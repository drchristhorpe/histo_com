from pathlib import Path

import numpy as np
import pytest

from histo_com import HistoCom, StructureError, com_all
from histo_com.core import load_structure

FIXTURES = Path(__file__).parent / "fixtures"
ABD = FIXTURES / "1hhk_1_abd.cif"
PEPTIDE = FIXTURES / "1hhk_1_peptide.cif"
COMPLEX = FIXTURES / "8gvi_1_aligned.cif"
ALIGNED_PDB = FIXTURES / "7r7y_1_aligned.pdb"


def _is_finite_vec3(v):
    arr = np.asarray(v)
    return arr.shape == (3,) and np.all(np.isfinite(arr))


def test_load_structure_cif():
    structure = load_structure(ABD)
    assert len(structure) >= 1
    assert len(list(structure[0])) == 1


def test_load_structure_missing_file():
    with pytest.raises(StructureError):
        load_structure(FIXTURES / "does_not_exist.cif")


def test_load_structure_unknown_extension(tmp_path):
    bogus = tmp_path / "structure.xyz"
    bogus.write_text("not a structure")
    with pytest.raises(StructureError):
        load_structure(bogus)


# -- 1hhk_1_abd.cif : mode all -----------------------------------------------


def test_abd_com_all():
    com = com_all(ABD)
    assert _is_finite_vec3(com)


def test_abd_com_all_matches_manual_average():
    h = HistoCom(ABD)
    com = h.com()
    atoms = list(h.structure[0].get_atoms())
    coords = np.array([a.coord for a in atoms], dtype=np.float64)
    masses = np.array([a.mass for a in atoms], dtype=np.float64)
    expected = np.average(coords, axis=0, weights=masses)
    assert np.allclose(com, expected)


# -- 1hhk_1_peptide.cif : mode all, and residues 1-9 -------------------------


def test_peptide_com_all():
    com = com_all(PEPTIDE)
    assert _is_finite_vec3(com)


def test_peptide_com_by_residues_range_object():
    h = HistoCom(PEPTIDE)
    coms = h.com_by_residues(range(1, 10))
    assert len(coms) == 9
    for com in coms:
        assert _is_finite_vec3(com)
    # Residues are distinct positions along the peptide backbone.
    assert len({tuple(np.round(c, 3)) for c in coms}) == 9


def test_peptide_com_by_residues_string_selector():
    h = HistoCom(PEPTIDE)
    coms_str = h.com_by_residues("1-9")
    coms_iter = h.com_by_residues(range(1, 10))
    assert len(coms_str) == len(coms_iter) == 9
    for a, b in zip(coms_str, coms_iter):
        assert np.allclose(a, b)


def test_peptide_single_chain_loaded_once():
    h = HistoCom(PEPTIDE)
    structure_before = h.structure
    h.com()
    h.com_by_residues(range(1, 10))
    h.com_by_domains(["C"])
    # Same in-memory structure object reused across all calls.
    assert h.structure is structure_before


# -- 8gvi_1_aligned.cif : domains --------------------------------------------


def test_complex_domain_single_chain_P():
    h = HistoCom(COMPLEX)
    coms = h.com_by_domains(["P"])
    assert len(coms) == 1
    assert _is_finite_vec3(coms[0])


def test_complex_domain_chain_L_range_1_180():
    h = HistoCom(COMPLEX)
    coms = h.com_by_domains(["L:1-180"])
    assert len(coms) == 1
    assert _is_finite_vec3(coms[0])
    # Chain L only has residues 1-99; the range should include all of them
    # and equal the whole-chain centre of mass.
    whole_chain = h.com_by_domains(["L"])[0]
    assert np.allclose(coms[0], whole_chain)


def test_complex_domain_chains_A_B():
    h = HistoCom(COMPLEX)
    coms = h.com_by_domains(["A", "B"])
    assert len(coms) == 2
    for com in coms:
        assert _is_finite_vec3(com)
    assert not np.allclose(coms[0], coms[1])


def test_complex_domain_string_selector_matches_iterable():
    h = HistoCom(COMPLEX)
    from_string = h.com_by_domains("A,B")
    from_iter = h.com_by_domains(["A", "B"])
    for a, b in zip(from_string, from_iter):
        assert np.allclose(a, b)


def test_complex_domain_unknown_chain_raises():
    h = HistoCom(COMPLEX)
    with pytest.raises(StructureError):
        h.com_by_domains(["Z"])


def test_complex_bare_range_ambiguous_chain_raises():
    h = HistoCom(COMPLEX)
    with pytest.raises(StructureError):
        h.com_by_domains(["1-180"])


# -- PDB format support -------------------------------------------------------


def test_pdb_format_matches_cif(tmp_path):
    from Bio.PDB import PDBIO

    pdb_path = tmp_path / "peptide.pdb"
    io = PDBIO()
    io.set_structure(load_structure(PEPTIDE))
    io.save(str(pdb_path))

    cif_com = com_all(PEPTIDE)
    pdb_com = com_all(pdb_path)
    assert np.allclose(cif_com, pdb_com)


# -- 7r7y_1_aligned.pdb : native PDB file, genuine partial-range domain -----


def test_native_pdb_domain_chain_A_partial_range():
    h = HistoCom(ALIGNED_PDB)
    coms = h.com_by_domains(["A:1-180"])
    assert len(coms) == 1
    assert _is_finite_vec3(coms[0])
    # Chain A has 276 residues, so 1-180 is a genuine truncation and
    # should differ from the whole-chain centre of mass.
    whole_chain = h.com_by_domains(["A"])[0]
    assert not np.allclose(coms[0], whole_chain)


def test_native_pdb_domain_matches_manual_average():
    # This structure has alternate-location (altloc) atoms in the 1-180
    # range, so the "expected" atom list must resolve disorder the same
    # way Biopython's own center_of_mass()/get_unpacked_list() does,
    # rather than naively iterating residues (which only sees the
    # currently-selected altloc per atom).
    h = HistoCom(ALIGNED_PDB)
    com = h.com_by_domains(["A:1-180"])[0]
    chain_a = h.model["A"]
    residues = [r for r in chain_a if 1 <= r.id[1] <= 180]
    atoms = [a for r in residues for a in r.get_unpacked_list()]
    coords = np.array([a.coord for a in atoms], dtype=np.float64)
    masses = np.array([a.mass for a in atoms], dtype=np.float64)
    expected = np.average(coords, axis=0, weights=masses)
    assert np.allclose(com, expected)


# -- write_com_pdb() ----------------------------------------------------------


def _read_marker_atoms(path):
    structure = load_structure(path)
    return list(structure[0].get_atoms())


def test_write_com_pdb_all_mode(tmp_path):
    h = HistoCom(ABD)
    out = h.write_com_pdb(tmp_path / "markers.pdb", mode="all")
    assert out.is_file()

    atoms = _read_marker_atoms(out)
    assert len(atoms) == 1
    assert atoms[0].get_parent().resname == "COM"
    assert np.allclose(atoms[0].coord, h.com(), atol=1e-3)


def test_write_com_pdb_domains_mode(tmp_path):
    h = HistoCom(COMPLEX)
    domains = "P,L:1-180,A,B"
    expected = h.com_by_domains(domains)
    out = h.write_com_pdb(tmp_path / "markers.pdb", mode="domains", domains=domains)

    atoms = _read_marker_atoms(out)
    assert len(atoms) == len(expected) == 4
    for atom, exp in zip(atoms, expected):
        assert np.allclose(atom.coord, exp, atol=1e-3)
    # Each domain marker is placed on its own chain for traceability.
    assert [a.get_parent().get_parent().id for a in atoms] == ["P", "L", "A", "B"]


def test_write_com_pdb_residues_mode(tmp_path):
    h = HistoCom(PEPTIDE)
    expected = h.com_by_residues(range(1, 10))
    out = h.write_com_pdb(tmp_path / "markers.pdb", mode="residues", residues="1-9")

    atoms = _read_marker_atoms(out)
    assert len(atoms) == len(expected) == 9
    for atom, exp in zip(atoms, expected):
        assert np.allclose(atom.coord, exp, atol=1e-3)
    # Residue markers keep the original residue numbering.
    assert [a.get_parent().id[1] for a in atoms] == list(range(1, 10))


def test_write_com_pdb_handles_duplicate_residue_request(tmp_path):
    # Requesting the same residue twice would otherwise collide on
    # (chain, resSeq) when building the marker structure.
    h = HistoCom(PEPTIDE)
    out = h.write_com_pdb(tmp_path / "markers.pdb", mode="residues", residues="5,5")

    atoms = _read_marker_atoms(out)
    assert len(atoms) == 2
    assert np.allclose(atoms[0].coord, atoms[1].coord, atol=1e-3)


def test_write_com_pdb_returns_path(tmp_path):
    h = HistoCom(ABD)
    out_path = tmp_path / "subdir_does_not_need_to_exist_check.pdb"
    result = h.write_com_pdb(out_path, mode="all")
    assert result == out_path


def test_write_com_pdb_unknown_mode_raises(tmp_path):
    h = HistoCom(ABD)
    with pytest.raises(ValueError):
        h.write_com_pdb(tmp_path / "markers.pdb", mode="bogus")
