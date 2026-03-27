"""Data model for polymer nodes in complex polymer notation."""

from __future__ import annotations

import re
from dataclasses import dataclass, field


@dataclass
class PolymerNode:
    """A node in the polymer graph (e.g., RNA1, PEPTIDE2, CHEM1)."""

    id: str = ""
    label: str = ""
    anotation: str = ""

    @property
    def type(self) -> str:
        return get_polymer_type(self.id)


def get_polymer_type(node_id: str) -> str:
    """Extract polymer type from node ID by stripping trailing digits."""
    result: list[str] = []
    for ch in node_id:
        if ch.isdigit():
            break
        result.append(ch)
    return "".join(result)
