"""Structure loading and centre-of-mass calculations."""

from __future__ import annotations

from collections import deque
from pathlib import Path
from typing import Iterable

import numpy as np
from Bio.PDB.Chain import Chain
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Structure import Structure

from histo_com.selectors import DomainRef, ResidueRef, parse_domains, parse_residues

_CIF_SUFFIXES = {".cif", ".mmcif"}
_PDB_SUFFIXES = {".pdb", ".ent"}

# Bio.PDB entity "levels" that may contain disorder duplicates and thus
# need get_unpacked_list() rather than raw child_list traversal.
_DISORDERABLE_LEVELS = {"R", "C"}


class StructureError(ValueError):
    """Raised for problems loading or querying a structure."""


def load_structure(path: str | Path, structure_id: str | None = None) -> Structure:
    """Parse a PDB or mmCIF file into a Bio.PDB Structure.

    Format is chosen from the file extension (case-insensitive):
    ``.cif``/``.mmcif`` -> mmCIF, ``.pdb``/``.ent`` -> legacy PDB.
    """
    path = Path(path)
    if not path.is_file():
        raise StructureError(f"No such file: {path}")

    suffix = path.suffix.lower()
    sid = structure_id or path.stem

    if suffix in _CIF_SUFFIXES:
        parser = MMCIFParser(QUIET=True)
    elif suffix in _PDB_SUFFIXES:
        parser = PDBParser(QUIET=True)
    else:
        raise StructureError(
            f"Unrecognised structure file extension {suffix!r} for {path}; "
            "expected one of .cif, .mmcif, .pdb, .ent"
        )

    structure = parser.get_structure(sid, str(path))
    if len(structure) == 0:
        raise StructureError(f"No models found in {path}")
    return structure


def centre_of_mass(entities: Iterable) -> np.ndarray:
    """Mass-weighted centre of mass over a homogeneous iterable of Bio.PDB
    entities (all Residues, all Chains, ...) or Atoms.

    Correctly resolves disordered residues/atoms (altlocs) by using the
    same unpacking approach as ``Bio.PDB.Entity.center_of_mass``,
    generalised to an arbitrary starting collection rather than a single
    entity's children.
    """
    queue = deque(entities)
    if not queue:
        raise StructureError("No atoms/residues to compute a centre of mass over")

    while {e.level for e in queue} != {"A"}:
        e = queue.popleft()
        if e.level in _DISORDERABLE_LEVELS:
            queue.extend(e.get_unpacked_list())
        else:
            queue.extend(e.child_list)

    atoms = list(queue)
    coords = np.asarray([a.coord for a in atoms], dtype=np.float64)
    weights = np.asarray([a.mass for a in atoms], dtype=np.float64)
    return np.average(coords, axis=0, weights=weights)


def _get_chain(model, chain_id: str) -> Chain:
    try:
        return model[chain_id]
    except KeyError:
        available = ", ".join(c.id for c in model)
        raise StructureError(
            f"No such chain {chain_id!r}; available chains: {available}"
        ) from None


def _default_chain(model) -> Chain:
    chains = list(model)
    if len(chains) != 1:
        ids = ", ".join(c.id for c in chains)
        raise StructureError(
            "A chain-less selector is only valid for single-chain structures; "
            f"this structure has {len(chains)} chains ({ids}). "
            "Prefix the selector with a chain id, e.g. 'A:1-180'."
        )
    return chains[0]


def _residues_in_range(chain: Chain, start: int, end: int) -> list:
    residues = [r for r in chain if start <= r.id[1] <= end]
    if not residues:
        raise StructureError(
            f"No residues in range {start}-{end} found in chain {chain.id!r}"
        )
    return residues


def _residue_by_number(chain: Chain, resseq: int):
    for r in chain:
        if r.id[1] == resseq:
            return r
    raise StructureError(f"No residue {resseq} found in chain {chain.id!r}")


class HistoCom:
    """Loads a structure once and computes centres of mass from it.

    >>> h = HistoCom("structure.cif")
    >>> h.com()
    >>> h.com_by_domains(["A", "B"])
    >>> h.com_by_residues(range(1, 10))
    """

    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.structure = load_structure(self.path)
        self.model = self.structure[0]

    def com(self) -> np.ndarray:
        """Centre of mass over every atom in the first model."""
        return centre_of_mass([self.model])

    def _domain_com(self, ref: DomainRef) -> np.ndarray:
        chain = _get_chain(self.model, ref.chain) if ref.chain else _default_chain(self.model)
        if ref.start is None:
            return centre_of_mass([chain])
        return centre_of_mass(_residues_in_range(chain, ref.start, ref.end))

    def com_by_domains(self, domains) -> list[np.ndarray]:
        """Centre of mass per domain (chain, or chain + residue range).

        ``domains`` may be a comma-separated selector string or an
        iterable of chain letters / selector tokens.
        """
        refs = parse_domains(domains)
        return [self._domain_com(ref) for ref in refs]

    def _residue_com(self, ref: ResidueRef) -> np.ndarray:
        chain = _get_chain(self.model, ref.chain) if ref.chain else _default_chain(self.model)
        residue = _residue_by_number(chain, ref.resseq)
        return centre_of_mass([residue])

    def com_by_residues(self, residues) -> list[np.ndarray]:
        """Centre of mass per residue (one per residue number, ranges expand).

        ``residues`` may be a comma-separated selector string or an
        iterable of residue numbers / selector tokens.
        """
        refs = parse_residues(residues)
        return [self._residue_com(ref) for ref in refs]


def com_all(path: str | Path) -> np.ndarray:
    """Convenience wrapper: centre of mass of an entire structure file."""
    return HistoCom(path).com()


def com_domains(path: str | Path, domains) -> list[np.ndarray]:
    """Convenience wrapper: centre of mass per domain for a structure file."""
    return HistoCom(path).com_by_domains(domains)


def com_residues(path: str | Path, residues) -> list[np.ndarray]:
    """Convenience wrapper: centre of mass per residue for a structure file."""
    return HistoCom(path).com_by_residues(residues)
