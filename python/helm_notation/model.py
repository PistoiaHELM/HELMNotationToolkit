"""Data model classes ported from org.helm.notation.model.*"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


# ---------------------------------------------------------------------------
# MoleculeInfo  (already partially in types.py, but we keep a richer version)
# ---------------------------------------------------------------------------


@dataclass
class MoleculeInfo:
    """Port of org.helm.notation.model.MoleculeInfo."""

    molecular_weight: float = 0.0
    exact_mass: float = 0.0
    molecular_formula: str = ""


# ---------------------------------------------------------------------------
# Attachment
# ---------------------------------------------------------------------------

# Well-known constants
BACKBONE_MONOMER_LEFT_ATTACHMENT = "R1"
BACKBONE_MONOMER_RIGHT_ATTACHMENT = "R2"
BACKBONE_MONOMER_BRANCH_ATTACHMENT = "R3"
BRANCH_MONOMER_ATTACHMENT = "R1"
STARTING_MONOMER_ATTACHMENT = "Starting"
CAP_GROUP_OH = "OH"
CAP_GROUP_H = "H"
PAIR_ATTACHMENT = "pair"


@dataclass
class Attachment:
    """Port of org.helm.notation.model.Attachment."""

    id: int = 0
    alternate_id: str = ""
    label: str = ""
    cap_group_name: str = ""
    cap_group_smiles: str = ""

    def __str__(self) -> str:
        return f"{self.label or ' '}: {self.cap_group_name}"


# ---------------------------------------------------------------------------
# Monomer
# ---------------------------------------------------------------------------

NUCLIEC_ACID_POLYMER_TYPE = "RNA"
PEPTIDE_POLYMER_TYPE = "PEPTIDE"
CHEMICAL_POLYMER_TYPE = "CHEM"
SUPPORTED_POLYMER_TYPES = [NUCLIEC_ACID_POLYMER_TYPE, PEPTIDE_POLYMER_TYPE, CHEMICAL_POLYMER_TYPE]

BACKBONE_MONOMER_TYPE = "Backbone"
BRANCH_MONOMER_TYPE = "Branch"
UNDEFINED_MONOMER_TYPE = "Undefined"

ID_P = "P"

MODIFICATION_START_SYMBOL = "["
MODIFICATION_END_SYMBOL = "]"

ATTACHMENT_LIST_DELIMITER = "$"


@dataclass
class Monomer:
    """Port of org.helm.notation.model.Monomer."""

    id: int = 0
    alternate_id: str = ""
    natural_analog: str = ""
    name: str = ""
    can_smiles: str = ""
    molfile: str = ""
    monomer_type: str = ""
    polymer_type: str = ""
    attachment_list: List[Attachment] = field(default_factory=list)
    new_monomer: bool = False
    ad_hoc_monomer: bool = False

    def get_attachment(self, label: str) -> Optional[Attachment]:
        for att in self.attachment_list:
            if att.label.lower() == label.lower():
                return att
        return None

    def get_cap_molecule_info(self, label: str) -> Optional[MoleculeInfo]:
        """Return MoleculeInfo for the cap group of the requested R-group."""
        from helm_notation.chem import ChemEngine

        att = self.get_attachment(label)
        if att and att.cap_group_smiles:
            return ChemEngine.get_molecule_info(att.cap_group_smiles)
        return None

    def is_modified(self) -> bool:
        return len(self.alternate_id) > 1

    def is_same_type(self, other: "Monomer") -> bool:
        a1 = self.natural_analog or self.alternate_id
        a2 = other.natural_analog or other.alternate_id
        return (
            self.monomer_type == other.monomer_type
            and self.polymer_type == other.polymer_type
            and a1 == a2
        )

    def get_attachment_list_string(self) -> str:
        ids = sorted(a.alternate_id for a in self.attachment_list)
        return ATTACHMENT_LIST_DELIMITER.join(ids)

    def attachment_equals(self, other: "Monomer") -> bool:
        return self.get_attachment_list_string() == other.get_attachment_list_string()

    def attachment_contains(self, other: "Monomer") -> bool:
        return other.get_attachment_list_string() in self.get_attachment_list_string()

    def contains_any_atom(self) -> bool:
        """Check if the monomer SMILES contains a generic 'A' (any) atom."""
        from helm_notation.chem import ChemEngine  # noqa: E501 circular-safe

        if not self.can_smiles:
            return False
        mol = ChemEngine.get_molecule(self.can_smiles)
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "A":
                return True
        return False


# ---------------------------------------------------------------------------
# MonomerCache
# ---------------------------------------------------------------------------


@dataclass
class MonomerCache:
    """Port of org.helm.notation.model.MonomerCache."""

    monomer_db: Dict[str, Dict[str, Monomer]] = field(default_factory=dict)
    smiles_monomer_db: Dict[str, Monomer] = field(default_factory=dict)
    attachment_db: Dict[str, Attachment] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Nucleotide
# ---------------------------------------------------------------------------

STARTING_POSITION_TYPE = 1
MIDDLE_POSITION_TYPE = 2
ENDING_POSITION_TYPE = 3


@dataclass
class Nucleotide:
    """Port of org.helm.notation.model.Nucleotide."""

    symbol: str = ""
    notation: str = ""
    natural_analog: str = ""
    position_type: int = MIDDLE_POSITION_TYPE

    def is_modified(self) -> bool:
        if self.notation:
            if "[" in self.notation:
                return True
            if "(X)" in self.notation:
                return True
            if self.notation.endswith(")"):
                return True
        return False

    def unmodified_without_phosphate(self) -> bool:
        if not self.notation:
            return True
        return "[" not in self.notation and "X" not in self.notation


# ---------------------------------------------------------------------------
# Fragment
# ---------------------------------------------------------------------------


@dataclass
class Fragment:
    """Port of org.helm.notation.model.Fragment."""

    source: str = ""
    type: str = ""
    symbol: str = ""
    notation: str = ""


# ---------------------------------------------------------------------------
# PolymerNode
# ---------------------------------------------------------------------------

_NODE_ID_RE = re.compile(r"^([A-Za-z]+)(\d+)$")


@dataclass
class PolymerNode:
    """Port of org.helm.notation.model.PolymerNode."""

    id: str = ""
    label: str = ""
    annotation: str = ""

    @property
    def type(self) -> str:
        return get_polymer_type(self.id)


def get_polymer_type(node_id: str) -> str:
    """Extract the polymer type prefix from a node id like 'RNA1'."""
    result: list[str] = []
    for ch in node_id:
        if ch.isdigit():
            break
        result.append(ch)
    return "".join(result)


# ---------------------------------------------------------------------------
# RNAPolymerNode
# ---------------------------------------------------------------------------


@dataclass
class RNAPolymerNode(PolymerNode):
    """Port of org.helm.notation.model.RNAPolymerNode."""

    sequence: str = ""
    modified_sequence: str = ""

    @property
    def type(self) -> str:
        return NUCLIEC_ACID_POLYMER_TYPE


# ---------------------------------------------------------------------------
# PolymerEdge
# ---------------------------------------------------------------------------

EDGE_COMPONENT_DELIMITER = ","
CONNECTION_SEPARATOR = "-"
MONOMER_ATTACHMENT_SEPARATOR = ":"
NODE_CONCATENATOR = "+"

STANDARD_EDGE = 1
PAIR_EDGE = 2
GENERIC_EDGE = 3

STANDARD_EDGE_ATTACHMENT = 1
PAIR_EDGE_ATTACHMENT = 2
GENERIC_EDGE_ATTACHMENT = 3

PAIR_EDGE_KEY = "pair"
GENERIC_EDGE_KEY = "generic"
UNKNOWN_ATTACHMENT = "?"


@dataclass
class PolymerEdge:
    """Port of org.helm.notation.model.PolymerEdge."""

    edge_notation: str = ""
    source_node: str = ""
    target_node: str = ""
    connection: str = ""
    edge_type: int = 0
    source_attachment_type: int = 0
    target_attachment_type: int = 0

    @property
    def source_nodes(self) -> List[str]:
        return self.source_node.split("+") if self.source_node else []

    @property
    def target_nodes(self) -> List[str]:
        return self.target_node.split("+") if self.target_node else []

    @property
    def source_connection(self) -> Optional[str]:
        if not self.connection:
            return None
        return self.connection.split(CONNECTION_SEPARATOR)[0]

    @property
    def target_connection(self) -> Optional[str]:
        if not self.connection:
            return None
        parts = self.connection.split(CONNECTION_SEPARATOR)
        return parts[1] if len(parts) == 2 else None

    @property
    def reverse_connection(self) -> str:
        return f"{self.target_connection}{CONNECTION_SEPARATOR}{self.source_connection}"

    @property
    def source_uid(self) -> str:
        return f"{self.source_node}{MONOMER_ATTACHMENT_SEPARATOR}{self.source_connection}"

    @property
    def target_uid(self) -> str:
        return f"{self.target_node}{MONOMER_ATTACHMENT_SEPARATOR}{self.target_connection}"

    def to_reverse_string(self) -> str:
        return (
            f"{self.target_node}{EDGE_COMPONENT_DELIMITER}"
            f"{self.source_node}{EDGE_COMPONENT_DELIMITER}"
            f"{self.target_connection}{CONNECTION_SEPARATOR}{self.source_connection}"
        )

    @property
    def source_r(self) -> Optional[str]:
        sc = self.source_connection
        if sc and self.source_attachment_type == STANDARD_EDGE_ATTACHMENT:
            parts = sc.split(MONOMER_ATTACHMENT_SEPARATOR)
            return parts[1] if len(parts) == 2 else None
        return None

    @property
    def target_r(self) -> Optional[str]:
        tc = self.target_connection
        if tc and self.target_attachment_type == STANDARD_EDGE_ATTACHMENT:
            parts = tc.split(MONOMER_ATTACHMENT_SEPARATOR)
            return parts[1] if len(parts) == 2 else None
        return None

    def __str__(self) -> str:
        return (
            f"{self.source_node}{EDGE_COMPONENT_DELIMITER}"
            f"{self.target_node}{EDGE_COMPONENT_DELIMITER}"
            f"{self.source_connection}{CONNECTION_SEPARATOR}{self.target_connection}"
        )


# ---------------------------------------------------------------------------
# ComplexPolymer
# ---------------------------------------------------------------------------


@dataclass
class ComplexPolymer:
    """Port of org.helm.notation.model.ComplexPolymer."""

    polymer_node_list: List[PolymerNode] = field(default_factory=list)
    polymer_edge_list: List[PolymerEdge] = field(default_factory=list)
    base_pair_list: List[PolymerEdge] = field(default_factory=list)
    polymer_node_annotation_map: Dict[str, str] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# RgroupStructure  (chemistry-dependent: uses Chem.Mol)
# ---------------------------------------------------------------------------


@dataclass
class RgroupStructure:
    """Port of org.helm.notation.model.RgroupStructure.

    molecule: an rdkit Chem.Mol (or None)
    rgroup_map: maps R-group label (e.g. "R1") → atom index in molecule
    """

    molecule: Any = None
    rgroup_map: Dict[str, int] = field(default_factory=dict)
