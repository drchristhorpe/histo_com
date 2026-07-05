"""histo_com: centre of mass calculations for 3D biological structures."""

from histo_com.core import (
    HistoCom,
    StructureError,
    centre_of_mass,
    com_all,
    com_domains,
    com_residues,
    load_structure,
)
from histo_com.selectors import DomainRef, ResidueRef, SelectorError

__version__ = "0.1.0"

__all__ = [
    "HistoCom",
    "StructureError",
    "SelectorError",
    "DomainRef",
    "ResidueRef",
    "centre_of_mass",
    "com_all",
    "com_domains",
    "com_residues",
    "load_structure",
    "__version__",
]
