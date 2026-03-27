"""Caching mechanism for monomer data."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from helm_notation.model.attachment import Attachment
    from helm_notation.model.monomer import Monomer


@dataclass
class MonomerCache:
    """Cache for monomer database, attachment database, and SMILES-monomer mapping."""

    monomer_db: dict[str, dict[str, Monomer]] = field(default_factory=dict)
    attachment_db: dict[str, Attachment] = field(default_factory=dict)
    smiles_monomer_db: dict[str, Monomer] = field(default_factory=dict)
