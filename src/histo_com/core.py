"""Structure loading and centre-of-mass calculations."""

from __future__ import annotations

from collections import deque
from pathlib import Path
from typing import Iterable

import numpy as np
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

from histo_com.selectors import (
    DomainRef,
    ResidueRef,
    format_domain,
    format_residue,
    parse_domains,
    parse_residues,
)

# Residue/atom naming used for the pseudo-atoms written by write_com_pdb().
_MARKER_RESNAME = "COM"
_MARKER_CHAIN_ID = "Z"

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

    def _resolve_domain(self, ref: DomainRef) -> tuple[str, np.ndarray]:
        """Resolve a (possibly multi-part) domain to its combined centre
        of mass, by gathering every part's residues into one homogeneous
        list before averaging."""
        residues: list = []
        chain_ids: list[str] = []
        for part in ref.parts:
            chain = _get_chain(self.model, part.chain) if part.chain else _default_chain(self.model)
            chain_ids.append(chain.id)
            if part.start is None:
                residues.extend(chain.get_unpacked_list())
            else:
                residues.extend(_residues_in_range(chain, part.start, part.end))
        coord = centre_of_mass(residues)
        unique_chain_ids = list(dict.fromkeys(chain_ids))
        chain_id = unique_chain_ids[0] if len(unique_chain_ids) == 1 else _MARKER_CHAIN_ID
        return chain_id, coord

    def com_by_domains(self, domains) -> list[np.ndarray]:
        """Centre of mass per domain.

        ``domains`` may be a single selector string — a chain, a
        chain+range, or several comma-joined chains/ranges combined into
        one domain (e.g. ``"A,B"``) — or an iterable of such selectors,
        each describing its own, separate domain.
        """
        refs = parse_domains(domains)
        return [self._resolve_domain(ref)[1] for ref in refs]

    def _resolve_residue(self, ref: ResidueRef) -> tuple[str, np.ndarray]:
        chain = _get_chain(self.model, ref.chain) if ref.chain else _default_chain(self.model)
        residue = _residue_by_number(chain, ref.resseq)
        return chain.id, centre_of_mass([residue])

    def com_by_residues(self, residues) -> list[np.ndarray]:
        """Centre of mass per residue (one per residue number, ranges expand).

        ``residues`` may be a comma-separated selector string or an
        iterable of residue numbers / selector tokens.
        """
        refs = parse_residues(residues)
        return [self._resolve_residue(ref)[1] for ref in refs]

    def _com_markers(self, mode: str, domains=None, residues=None) -> list[tuple[str, int, str, np.ndarray]]:
        """Resolve (chain_id, resseq, label, coord) rows for write_com_pdb()."""
        if mode == "all":
            return [(_MARKER_CHAIN_ID, 1, "all", self.com())]

        if mode == "domains":
            refs = parse_domains(domains)
            rows = []
            for i, ref in enumerate(refs, start=1):
                chain_id, coord = self._resolve_domain(ref)
                rows.append((chain_id, i, format_domain(ref), coord))
            return rows

        if mode == "residues":
            refs = parse_residues(residues)
            rows = []
            for ref in refs:
                chain_id, coord = self._resolve_residue(ref)
                rows.append((chain_id, ref.resseq, format_residue(ref), coord))
            return rows

        raise ValueError(f"Unknown mode {mode!r}; expected 'all', 'domains', or 'residues'")

    def write_com_pdb(self, output_path: str | Path, mode: str = "all", domains=None, residues=None) -> Path:
        """Write a PDB file containing one pseudo-atom per computed centre
        of mass, so it can be viewed alongside the structure in a
        molecular viewer.

        Each marker is a ``HETATM`` residue named ``COM`` with a single
        atom also named ``COM``. ``mode``/``domains``/``residues`` follow
        the same semantics as :meth:`com`, :meth:`com_by_domains`, and
        :meth:`com_by_residues`.
        """
        rows = self._com_markers(mode, domains=domains, residues=residues)
        output_path = Path(output_path)
        _write_marker_pdb(rows, output_path)
        return output_path


def _write_marker_pdb(rows: list[tuple[str, int, str, np.ndarray]], output_path: Path) -> None:
    """Write ``(chain_id, resseq, label, coord)`` rows as HETATM pseudo-atoms."""
    if not rows:
        raise StructureError("No centre-of-mass markers to write")

    structure = Structure("com_markers")
    model = Model(0)
    structure.add(model)

    chains: dict[str, Chain] = {}
    seen_ids: set[tuple[str, int, str]] = set()
    for serial, (chain_id, resseq, _label, coord) in enumerate(rows, start=1):
        chain = chains.get(chain_id)
        if chain is None:
            chain = Chain(chain_id)
            model.add(chain)
            chains[chain_id] = chain

        # Two selector tokens can resolve to the same (chain, residue
        # number), e.g. overlapping domain ranges or a repeated residue
        # token; disambiguate with an insertion code rather than crash.
        icode = " "
        for suffix in " ABCDEFGHIJKLMNOPQRSTUVWXYZ":
            icode = suffix
            if (chain_id, resseq, icode) not in seen_ids:
                break
        seen_ids.add((chain_id, resseq, icode))

        residue = Residue((f"H_{_MARKER_RESNAME}", resseq, icode), _MARKER_RESNAME, "")
        atom = Atom(
            name=_MARKER_RESNAME,
            coord=np.asarray(coord, dtype=np.float64),
            bfactor=0.0,
            occupancy=1.0,
            altloc=" ",
            fullname=f" {_MARKER_RESNAME}",
            serial_number=serial,
            element="C",
        )
        residue.add(atom)
        chain.add(residue)

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_path))


def com_all(path: str | Path) -> np.ndarray:
    """Convenience wrapper: centre of mass of an entire structure file."""
    return HistoCom(path).com()


def com_domains(path: str | Path, domains) -> list[np.ndarray]:
    """Convenience wrapper: centre of mass per domain for a structure file."""
    return HistoCom(path).com_by_domains(domains)


def com_residues(path: str | Path, residues) -> list[np.ndarray]:
    """Convenience wrapper: centre of mass per residue for a structure file."""
    return HistoCom(path).com_by_residues(residues)
