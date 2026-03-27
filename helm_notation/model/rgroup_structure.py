"""Data model for chemical structures with R-groups."""

from __future__ import annotations

from dataclasses import dataclass, field

from rdkit import Chem


@dataclass
class RgroupStructure:
    """A chemical structure with R-group attachment points.

    molecule: RDKit Mol object
    rgroup_map: Maps R-group label (e.g., "R1") to atom index in the molecule
    """

    molecule: Chem.Mol | None = None
    rgroup_map: dict[str, int] = field(default_factory=dict)
