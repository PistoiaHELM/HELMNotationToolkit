"""Data model for a peptide fragment during recursive bond-breaking."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from rdkit import Chem


@dataclass
class PeptideFragment:
    """Represents a molecular fragment during peptide structure parsing.

    Attributes:
        id: Unique identifier assigned during bond-breaking.
        level: Depth level in the fragmentation tree.
        molecule: RDKit Mol object for this fragment.
        monomer_ids_list: List of linear peptide chains, each a list of monomer IDs.
        connection_maps_list: Parallel structure to monomer_ids_list; for each monomer,
            a dict mapping a connection key (e.g. ``"R3"``) to the R-group label.
    """

    id: int = 0
    level: int = 0
    molecule: Optional[Chem.Mol] = None
    monomer_ids_list: Optional[list[list[str]]] = None
    connection_maps_list: Optional[list[list[dict[str, str]]]] = None
