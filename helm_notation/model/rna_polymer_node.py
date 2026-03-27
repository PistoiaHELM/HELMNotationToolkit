"""Data model for RNA polymer nodes."""

from __future__ import annotations

from dataclasses import dataclass

from helm_notation.model.polymer_node import PolymerNode


@dataclass
class RNAPolymerNode(PolymerNode):
    """Specialized polymer node for RNA structures."""

    sequence: str = ""
    modified_sequence: str = ""

    @property
    def type(self) -> str:
        return "RNA"

    @classmethod
    def from_polymer_node(cls, node: PolymerNode) -> RNAPolymerNode:
        return cls(
            id=node.id,
            label=node.label,
            anotation=node.anotation,
        )

    def __str__(self) -> str:
        return f"{self.id}\t{self.label}\t{self.anotation}\t{self.sequence}\t{self.modified_sequence}"
