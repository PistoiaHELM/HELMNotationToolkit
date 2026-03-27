"""Data model for complex polymers."""

from __future__ import annotations

from dataclasses import dataclass, field

from helm_notation.model.polymer_node import PolymerNode
from helm_notation.model.polymer_edge import PolymerEdge


@dataclass
class ComplexPolymer:
    """Root model containing lists of nodes, edges, and base pairs."""

    polymer_node_list: list[PolymerNode] = field(default_factory=list)
    polymer_edge_list: list[PolymerEdge] = field(default_factory=list)
    base_pair_list: list[PolymerEdge] = field(default_factory=list)
    polymer_node_annotation_map: dict[str, str] = field(default_factory=dict)
