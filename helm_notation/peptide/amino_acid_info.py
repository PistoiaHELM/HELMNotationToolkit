"""Data model for amino acid identification results."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class AminoAcidInfo:
    """Holds the monomer ID and R-group mapping for an identified amino acid fragment."""

    monomer_id: str = ""
    rgroup_map: dict[str, str] = field(default_factory=dict)
