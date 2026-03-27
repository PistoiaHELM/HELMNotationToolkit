"""Parser for complex polymer notation in HELM.

Complex Polymer Notation is based on graph theory, and uses simple polymer
notation and chemical structure as building blocks.

Complex notation follows the format:
    NodeList$EdgeList$BasePairList$NodeLabelList$OtherTBD

where EdgeList, BasePairList, NodeLabelList and OtherTBD could be empty,
but still keep the positional $.

Copyright (C) 2012, The Pistoia Alliance

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from __future__ import annotations

import logging
import re
from collections import defaultdict
from typing import TYPE_CHECKING

from rdkit import Chem

from helm_notation.constants import MONOMER_COUNT_THRESHOLD
from helm_notation.exceptions import (
    MonomerException,
    NotationException,
    StructureException,
)
from helm_notation.model.complex_polymer import ComplexPolymer
from helm_notation.model.molecule_info import MoleculeInfo
from helm_notation.model.monomer import (
    CHEMICAL_POLYMER_TYPE,
    NUCLEIC_ACID_POLYMER_TYPE,
    PEPTIDE_POLYMER_TYPE,
    Monomer,
)
from helm_notation.model.polymer_edge import (
    CONNECTION_SEPARATOR,
    EDGE_COMPONENT_DELIMITER,
    GENERIC_EDGE,
    GENERIC_EDGE_ATTACHMENT,
    MONOMER_ATTACHMENT_SEPARATOR,
    NODE_CONCATENATOR,
    NODE_CONCATENATOR_REGEX,
    PAIR_EDGE,
    STANDARD_EDGE,
    STANDARD_EDGE_ATTACHMENT,
    PolymerEdge,
)
from helm_notation.model.polymer_node import PolymerNode, get_polymer_type
from helm_notation.model.rgroup_structure import RgroupStructure
from helm_notation.model.rna_polymer_node import RNAPolymerNode
from helm_notation.model.nucleotide import Nucleotide
from helm_notation.tools import edge_parser
from helm_notation.tools import structure_parser
from helm_notation.tools import monomer_parser as monomer_parser_module

if TYPE_CHECKING:
    from helm_notation.monomer_store import MonomerStore

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
TOP_LEVEL_STOPPER = "}$"
TOP_LEVEL_DELIMITER = "$"
LIST_LEVEL_DELIMITER = "|"
LIST_LEVEL_DELIMITER_REGEX = r"\|"
NODE_LABEL_START_SYMBOL = "{"
NODE_LABEL_END_SYMBOL = "}"
INVALID_POLYMER_NODE = (
    "Polymer node notation must be in the format of PolymerTypeNumber{Notation}"
)
INVALID_NODE_ID = "Polymer node ID must be in the format of PolymerTypeNumber"
DEFAULT_PADDING_CHAR = " "
DEFAULT_BASE_PAIR_CHAR = "|"


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _check_for_monomer_store(monomer_store: MonomerStore | None) -> MonomerStore | None:
    """Return *monomer_store* unchanged when it is not ``None``/empty,
    otherwise fetch a fresh store from the global :class:`MonomerFactory`.
    """
    if monomer_store is not None:
        return monomer_store
    try:
        from helm_notation.monomer_factory import get_instance as get_monomer_factory
        factory = get_monomer_factory()
        return factory.monomer_store
    except Exception:
        logger.exception("Unable to get default monomer store")
        return None


# ---------------------------------------------------------------------------
# Public API -- SMILES helpers
# ---------------------------------------------------------------------------

def get_complex_polymer_smiles(
    extended_notation: str,
    monomer_store: MonomerStore | None = None,
) -> str:
    """Return the unique SMILES string for a complex polymer notation.

    All monomers must have specific structures.
    """
    monomer_store = _check_for_monomer_store(monomer_store)
    mol_list = get_complex_polymer_structure(extended_notation, monomer_store)
    parts: list[str] = []
    for mol in mol_list:
        smi = Chem.MolToSmiles(mol)
        parts.append(smi)
    mixture_smiles = ".".join(parts)
    mol = structure_parser.get_molecule(mixture_smiles)
    return Chem.MolToSmiles(mol)


def get_complex_polymer_canonical_smiles(
    extended_notation: str,
    monomer_store: MonomerStore | None = None,
) -> str:
    """Return the canonical unique SMILES for a complex polymer notation."""
    monomer_store = _check_for_monomer_store(monomer_store)
    mol_list = get_complex_polymer_structure(extended_notation, monomer_store)
    parts: list[str] = []
    for mol in mol_list:
        smi = structure_parser.get_unique_smiles(Chem.MolToSmiles(mol))
        parts.append(smi)
    mixture_smiles = ".".join(parts)
    return structure_parser.get_unique_smiles(mixture_smiles)


# ---------------------------------------------------------------------------
# Public API -- generic structure check
# ---------------------------------------------------------------------------

def contains_generic_structure(
    extended_notation: str,
    monomer_store: MonomerStore | None = None,
) -> bool:
    """Return ``True`` if the notation contains a generic (non-SMILES) structure."""
    monomer_store = _check_for_monomer_store(monomer_store)

    # Lazy import to avoid circular dependency
    from helm_notation.tools import simple_notation_parser

    all_node_string = get_all_node_string(extended_notation)
    node_list = get_polymer_node_list(all_node_string, monomer_store)
    for node in node_list:
        if node.id.startswith(CHEMICAL_POLYMER_TYPE):
            monomer_id = node.label
            monomer = simple_notation_parser.get_monomer(
                monomer_id, CHEMICAL_POLYMER_TYPE, monomer_store
            )
            if monomer.can_smiles is None or len(monomer.can_smiles) == 0:
                return True
            if monomer.contain_any_atom():
                return True
    return False


# ---------------------------------------------------------------------------
# Public API -- structure generation
# ---------------------------------------------------------------------------

def get_complex_polymer_structure(
    extended_notation: str,
    monomer_store: MonomerStore | None = None,
) -> list[Chem.Mol]:
    """Return a list of RDKit ``Mol`` objects with all R-groups filled."""
    from helm_notation.tools import simple_notation_parser

    monomer_store = _check_for_monomer_store(monomer_store)

    total_monomer_count = get_total_monomer_count(extended_notation, monomer_store)
    if total_monomer_count > MONOMER_COUNT_THRESHOLD:
        raise NotationException(
            f"Total monomer count [{total_monomer_count}] is above support "
            f"threshold [{MONOMER_COUNT_THRESHOLD}]"
        )

    complex_polymer = parse(extended_notation, monomer_store)
    validate_complex_polymer(complex_polymer, monomer_store)

    node_list = complex_polymer.polymer_node_list
    edge_list = complex_polymer.polymer_edge_list

    if edge_list is None:
        edge_list = []
    else:
        for edge in edge_list:
            if edge.edge_type != STANDARD_EDGE:
                raise NotationException(
                    "Polymer notation contains non-standard connection"
                )

    # Convert each polymer node into RgroupStructure, prefix rgroup keys with nodeID
    node_struc_map = _get_polymer_node_structure_map(node_list, monomer_store)
    for node_id, struc in node_struc_map.items():
        rg_map = struc.rgroup_map
        new_rg_map: dict[str, int] = {}
        for key, val in rg_map.items():
            new_rg_map[f"{node_id}:{key}"] = val
        struc.rgroup_map = new_rg_map

    # Group connected nodes
    group_node_list_map = _get_group_node_list_map(node_list, edge_list)

    # Convert each connected group into a single RgroupStructure
    group_structure_map: dict[int, RgroupStructure] = {}
    for group, n_list in group_node_list_map.items():
        if len(n_list) == 1:
            struc = node_struc_map[n_list[0]]
            group_structure_map[group] = struc
        else:
            conn_edge_list = _get_connected_edge_list(n_list, edge_list)
            merged_node_struc_map: dict[tuple[str, ...], RgroupStructure] = {}

            while len(conn_edge_list) > 0:
                edge = conn_edge_list[0]
                source_id = edge.source_node
                target_id = edge.target_node

                source_id_list: list[str] = []
                source_struc: RgroupStructure
                source_entry = _get_structure_entry(source_id, merged_node_struc_map)
                if source_entry is None:
                    source_struc = node_struc_map[source_id]
                else:
                    source_id_list = list(source_entry[0])
                    source_struc = source_entry[1]

                target_id_list: list[str] = []
                target_struc: RgroupStructure
                target_entry = _get_structure_entry(target_id, merged_node_struc_map)
                if target_entry is None:
                    target_struc = node_struc_map[target_id]
                else:
                    target_id_list = list(target_entry[0])
                    target_struc = target_entry[1]

                source_r = edge.source_uid  # e.g. RNA1:5:R1
                target_r = edge.target_uid
                source_atom_map = source_struc.rgroup_map
                target_atom_map = target_struc.rgroup_map

                temp_mol = source_struc.molecule
                structure_parser.merge(
                    temp_mol,
                    source_atom_map[source_r],
                    target_struc.molecule,
                    target_atom_map[target_r],
                )

                # Build new structure after connection
                dirty_struc = RgroupStructure()
                dirty_struc.molecule = temp_mol
                rgroup_map: dict[str, int] = {}
                for k, v in source_atom_map.items():
                    rgroup_map[k] = v
                for k, v in target_atom_map.items():
                    rgroup_map[k] = v
                rgroup_map.pop(source_r, None)
                rgroup_map.pop(target_r, None)
                dirty_struc.rgroup_map = rgroup_map

                # Build new ID list
                new_id_list: list[str] = []
                if not source_id_list and not target_id_list:
                    new_id_list = [source_id, target_id]
                elif not source_id_list and target_id_list:
                    new_id_list = list(target_id_list)
                    if source_id not in new_id_list:
                        new_id_list.append(source_id)
                    merged_node_struc_map.pop(tuple(target_id_list), None)
                elif source_id_list and not target_id_list:
                    new_id_list = list(source_id_list)
                    if target_id not in new_id_list:
                        new_id_list.append(target_id)
                    merged_node_struc_map.pop(tuple(source_id_list), None)
                else:
                    new_id_list = list(source_id_list)
                    for tid in target_id_list:
                        if tid not in new_id_list:
                            new_id_list.append(tid)
                    merged_node_struc_map.pop(tuple(source_id_list), None)
                    merged_node_struc_map.pop(tuple(target_id_list), None)

                merged_node_struc_map[tuple(new_id_list)] = dirty_struc
                conn_edge_list.remove(edge)

            # Only one entry finally
            merged_struc = next(iter(merged_node_struc_map.values()))
            group_structure_map[group] = merged_struc

    # Remove all remaining R groups not used in inter-polymer connection
    return _get_molecule_list(node_list, group_structure_map, monomer_store)


def _get_molecule_list(
    node_list: list[PolymerNode],
    group_structure_map: dict[int, RgroupStructure],
    monomer_store: MonomerStore | None,
) -> list[Chem.Mol]:
    """Cap unused R-groups and return the list of final molecules."""
    from helm_notation.tools import simple_notation_parser

    result: list[Chem.Mol] = []
    for struc in group_structure_map.values():
        mol = struc.molecule
        rgroup_map = struc.rgroup_map

        for key, r_atom_idx in list(rgroup_map.items()):
            # key format: NodeID:MonomerNumber:R#
            key_items = key.split(MONOMER_ATTACHMENT_SEPARATOR)
            node_id = key_items[0]
            polymer_type = get_polymer_type(node_id)
            node_notation = _get_polymer_notation(node_id, node_list)
            monomer_id_list = simple_notation_parser.get_monomer_id_list(
                node_notation, polymer_type, monomer_store
            )

            monomer_number = int(key_items[1])
            monomer_id = monomer_id_list[monomer_number - 1]
            monomer = simple_notation_parser.get_monomer(
                monomer_id, polymer_type, monomer_store
            )
            attachments = monomer.attachment_list

            rlabel = key_items[2]
            for att in attachments:
                if att.cap_group_smiles is None:
                    monomer_parser_module.fill_attachment_info(att)
                if att.get_label() == rlabel:
                    att_mol = structure_parser.get_molecule(att.cap_group_smiles)
                    rgroup_id = int(rlabel[1:])
                    att_atom = structure_parser.find_rgroup_atom(att_mol, rgroup_id)
                    structure_parser.merge(
                        mol, rgroup_map[key], att_mol, att_atom
                    )

        result.append(mol)
    return result


def _get_structure_entry(
    node_id: str,
    merged_node_struc_map: dict[tuple[str, ...], RgroupStructure],
) -> tuple[tuple[str, ...], RgroupStructure] | None:
    """Find an entry in *merged_node_struc_map* whose key-tuple contains *node_id*."""
    for ids, struc in merged_node_struc_map.items():
        if node_id in ids:
            return (ids, struc)
    return None


# ---------------------------------------------------------------------------
# Group / edge helpers
# ---------------------------------------------------------------------------

def _get_group_node_list_map(
    all_node_list: list[PolymerNode],
    all_edge_list: list[PolymerEdge],
) -> dict[int, list[str]]:
    """Break ``PolymerNode`` lists into connected groups based on edge info.

    The list values may contain duplicated entries for self-connection or
    multiple connections; a single-element list indicates a standalone node.
    """
    group_node_list_map: dict[int, list[str]] = {}
    group_id = 0

    # Put connected nodes into groups
    for edge in all_edge_list:
        source_id = edge.source_node
        target_id = edge.target_node
        overlap_groups: list[int] = []
        for gid, node_ids in group_node_list_map.items():
            if source_id in node_ids or target_id in node_ids:
                overlap_groups.append(gid)

        if not overlap_groups:
            group_id += 1
            group_node_list_map[group_id] = [source_id, target_id]
        else:
            if len(overlap_groups) > 1:
                # Merge all overlapping groups into the first one
                first_group = overlap_groups[0]
                for j in range(1, len(overlap_groups)):
                    other_group = overlap_groups[j]
                    group_node_list_map[first_group].extend(
                        group_node_list_map[other_group]
                    )
                    del group_node_list_map[other_group]

            first_group = overlap_groups[0]
            group_node_list_map[first_group].append(source_id)
            group_node_list_map[first_group].append(target_id)

    # Put standalone nodes into separate groups
    for node in all_node_list:
        node_id = node.id
        found_parent = False
        for node_ids in group_node_list_map.values():
            if node_id in node_ids:
                found_parent = True
                break
        if not found_parent:
            group_id += 1
            group_node_list_map[group_id] = [node_id]

    return group_node_list_map


def _get_connected_edge_list(
    connected_node_list: list[str],
    all_edge_list: list[PolymerEdge],
) -> list[PolymerEdge]:
    """Return edges whose source and target are both in *connected_node_list*."""
    result: list[PolymerEdge] = []
    for edge in all_edge_list:
        if (
            edge.source_node in connected_node_list
            and edge.target_node in connected_node_list
        ):
            result.append(edge)
    return result


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate_complex_notation(
    extended_notation: str,
    monomer_store: MonomerStore | None = None,
) -> bool:
    """Validate the complex polymer notation. Returns ``True`` on success."""
    monomer_store = _check_for_monomer_store(monomer_store)
    validate_notation_format(extended_notation)
    cp = parse(extended_notation, monomer_store)
    validate_complex_polymer(cp, monomer_store)
    return True


def validate_complex_polymer(
    complex_polymer: ComplexPolymer,
    monomer_store: MonomerStore | None = None,
) -> bool:
    """Check that nodes and edges match, and validate each node label."""
    monomer_store = _check_for_monomer_store(monomer_store)

    node_list = complex_polymer.polymer_node_list
    edge_list = complex_polymer.polymer_edge_list
    bp_list = complex_polymer.base_pair_list
    annotation_map = complex_polymer.polymer_node_annotation_map

    # Must have at least one node
    if not node_list:
        raise NotationException(
            "Complex notation must have at least one polymer node"
        )
    node_map: dict[str, str] = {}
    for node in node_list:
        node_map[node.id] = node.label
    if len(node_list) > len(node_map):
        raise NotationException("Polymer node IDs are not unique")

    # This call also validates simple polymers
    polymer_node_structure_map = _get_polymer_node_structure_map(
        node_list, monomer_store
    )

    if edge_list:
        for edge in edge_list:
            if edge.edge_type == STANDARD_EDGE:
                _validate_standard_edge(edge, node_map, polymer_node_structure_map)
            elif edge.edge_type == GENERIC_EDGE:
                _validate_generic_edge(edge, node_map, polymer_node_structure_map)
            elif edge.edge_type == PAIR_EDGE:
                _validate_pair_edge(edge, node_map)
            else:
                raise NotationException(f"Invalid edge type: {edge.edge_type}")

    if bp_list:
        for edge in bp_list:
            _validate_pair_edge(edge, node_map)

    if annotation_map:
        for aid in annotation_map:
            if aid not in node_map:
                raise NotationException(
                    "Polymer annotation contains unknown polymer node ID"
                )

    return True


def _validate_standard_edge(
    edge: PolymerEdge,
    node_map: dict[str, str],
    polymer_node_structure_map: dict[str, RgroupStructure],
) -> None:
    node = edge.source_node
    connection = edge.source_connection
    _validate_standard_attachment(node, connection, node_map, polymer_node_structure_map)
    node = edge.target_node
    connection = edge.target_connection
    _validate_standard_attachment(node, connection, node_map, polymer_node_structure_map)


def _validate_pair_edge(
    edge: PolymerEdge,
    node_map: dict[str, str],
) -> None:
    if edge.edge_type != PAIR_EDGE:
        raise NotationException(
            f"Invalid base pair edge: {edge.edge_notation} => {edge}"
        )
    if edge.source_node not in node_map:
        raise NotationException(
            "Polymer edge contains unknown polymer node ID"
        )
    if edge.target_node not in node_map:
        raise NotationException(
            "Polymer edge contains unknown polymer node ID"
        )


def _validate_generic_edge(
    edge: PolymerEdge,
    node_map: dict[str, str],
    polymer_node_structure_map: dict[str, RgroupStructure],
) -> None:
    # Source node check
    att_type = edge.source_attachment_type
    if att_type == STANDARD_EDGE_ATTACHMENT:
        _validate_standard_attachment(
            edge.source_node,
            edge.source_connection,
            node_map,
            polymer_node_structure_map,
        )
    elif att_type == GENERIC_EDGE_ATTACHMENT:
        _validate_generic_attachment(edge.source_nodes, node_map)

    # Target node check
    att_type = edge.target_attachment_type
    if att_type == STANDARD_EDGE_ATTACHMENT:
        _validate_standard_attachment(
            edge.target_node,
            edge.target_connection,
            node_map,
            polymer_node_structure_map,
        )
    elif att_type == GENERIC_EDGE_ATTACHMENT:
        _validate_generic_attachment(edge.target_nodes, node_map)


def _validate_standard_attachment(
    node: str,
    connection: str | None,
    node_map: dict[str, str],
    polymer_node_structure_map: dict[str, RgroupStructure],
) -> None:
    rstructure = polymer_node_structure_map.get(node)
    if node not in node_map:
        raise NotationException(
            "Polymer edge contains unknown polymer node ID"
        )
    if rstructure is None:
        raise NotationException(
            "Polymer edge contains polymer node without connection point"
        )
    r_map = rstructure.rgroup_map
    if r_map is None or connection not in r_map:
        raise NotationException(
            "Polymer edge contains polymer node without connection point"
        )


def _validate_generic_attachment(
    nodes: list[str],
    node_map: dict[str, str],
) -> None:
    for node in nodes:
        if node not in node_map:
            raise NotationException(
                "Polymer edge contains unknown polymer node ID"
            )


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse(
    extended_notation: str,
    monomer_store: MonomerStore | None = None,
) -> ComplexPolymer:
    """Parse *extended_notation* into a :class:`ComplexPolymer`."""
    if monomer_store is None:
        from helm_notation.monomer_factory import get_instance as get_monomer_factory
        try:
            factory = get_monomer_factory()
        except Exception as ex:
            raise NotationException(
                "Unable to initialize monomer factory"
            ) from ex
        monomer_store = factory.monomer_store

    cp = ComplexPolymer()

    all_node_string = get_all_node_string(extended_notation)
    nodes = get_polymer_node_list(all_node_string, monomer_store)
    cp.polymer_node_list = nodes

    all_edge_string = get_all_edge_string(extended_notation)
    edges = _get_polymer_edge_list(all_edge_string)
    cp.polymer_edge_list = edges

    all_base_pair_string = get_all_base_pair_string(extended_notation)
    base_pair_edges = _get_polymer_edge_list(all_base_pair_string)
    cp.base_pair_list = base_pair_edges

    all_node_label_string = get_all_node_label_string(extended_notation)
    annotation_map = _get_polymer_node_id_annotation_map(all_node_label_string)
    cp.polymer_node_annotation_map = annotation_map

    return cp


# ---------------------------------------------------------------------------
# Component-string extraction
# ---------------------------------------------------------------------------

def get_all_node_string(extended_notation: str) -> str:
    return _get_component_string(extended_notation, 0, False)


def get_all_edge_string(extended_notation: str) -> str:
    return _get_component_string(extended_notation, 1, False)


def get_all_base_pair_string(extended_notation: str) -> str:
    return _get_component_string(extended_notation, 2, False)


def get_all_node_label_string(extended_notation: str) -> str:
    return _get_component_string(extended_notation, 3, False)


def get_other_string(extended_notation: str) -> str:
    return _get_component_string(extended_notation, 4, True)


def _get_component_string(
    extended_notation: str,
    component_index: int,
    is_last: bool,
) -> str:
    """Return the string for the specified positional component."""
    validate_notation_format(extended_notation)
    count = 0
    start_pos = 0
    end_pos = 0

    first_pos = extended_notation.index(TOP_LEVEL_STOPPER) + 1

    if component_index == 0:
        end_pos = first_pos
    else:
        end_pos = first_pos
        while count < component_index:
            start_pos = end_pos + 1
            end_pos = extended_notation.index(TOP_LEVEL_DELIMITER, start_pos)
            count += 1

    if is_last:
        return extended_notation[start_pos:]
    else:
        return extended_notation[start_pos:end_pos]


def validate_notation_format(extended_notation: str) -> bool:
    """Validate that the notation has the required ``$$$$`` delimiters."""
    pos = extended_notation.find(TOP_LEVEL_STOPPER)
    if pos < 0:
        raise NotationException(
            "Invalid complex notation format,missing first component stopper }$"
        )
    start_pos = pos + len(TOP_LEVEL_STOPPER)

    count = 0
    while count < 3:
        pos = extended_notation.find(TOP_LEVEL_DELIMITER, start_pos)
        if pos < 0:
            raise NotationException(
                "Invalid complex notation format, must have four positional "
                "delimieters $$$$"
            )
        start_pos = pos + 1
        count += 1

    return True


# ---------------------------------------------------------------------------
# Node-list parsing
# ---------------------------------------------------------------------------

def get_polymer_node_list(
    all_node_string: str,
    monomer_store: MonomerStore | None = None,
) -> list[PolymerNode]:
    """Convert the *all_node_string* to a list of :class:`PolymerNode`.

    CHEM nodes may contain adhoc monomers with SMILES that include ``|``.
    """
    monomer_store = _check_for_monomer_store(monomer_store)

    result: list[PolymerNode] = []
    if not all_node_string:
        return result

    nodes = all_node_string
    if TOP_LEVEL_STOPPER in all_node_string:
        nodes = get_all_node_string(all_node_string)

    node_start_pos = 0
    delimiter_start_pos = nodes.find(NODE_LABEL_START_SYMBOL, node_start_pos)
    end_pos = nodes.find(NODE_LABEL_END_SYMBOL, node_start_pos)
    if end_pos >= 0:
        end_pos += 1

    while delimiter_start_pos >= 0 and end_pos > 0:
        node_string = nodes[node_start_pos:end_pos]
        node = _get_polymer_node(node_string, monomer_store)
        result.append(node)

        node_start_pos = end_pos + 1
        delimiter_start_pos = nodes.find(NODE_LABEL_START_SYMBOL, node_start_pos)
        end_pos = nodes.find(NODE_LABEL_END_SYMBOL, node_start_pos)
        if end_pos >= 0:
            end_pos += 1

    return result


def _get_polymer_node(
    node_string: str,
    monomer_store: MonomerStore | None = None,
) -> PolymerNode:
    """Convert a single node string into a :class:`PolymerNode`."""
    from helm_notation.tools import simple_notation_parser

    monomer_store = _check_for_monomer_store(monomer_store)

    start_symbol_pos = node_string.find(NODE_LABEL_START_SYMBOL)
    end_symbol_pos = node_string.find(NODE_LABEL_END_SYMBOL)

    if start_symbol_pos <= 0 or end_symbol_pos <= 0:
        raise NotationException(f"{INVALID_POLYMER_NODE}: {node_string}")
    if end_symbol_pos != len(node_string) - 1:
        raise NotationException(f"{INVALID_POLYMER_NODE}: {node_string}")

    node_id = node_string[:start_symbol_pos]
    _validate_node_id(node_id)
    label = node_string[start_symbol_pos + 1 : end_symbol_pos]

    # Add adhoc chem monomer if applicable
    if node_id.startswith(CHEMICAL_POLYMER_TYPE):
        label = simple_notation_parser.process_node(
            label, CHEMICAL_POLYMER_TYPE, "X", monomer_store
        )

    node = PolymerNode()
    node.id = node_id
    node.label = label
    return node


def _get_polymer_notation(
    node_id: str,
    node_list: list[PolymerNode],
) -> str | None:
    """Return the notation label for the given *node_id*."""
    for node in node_list:
        if node.id == node_id:
            return node.label
    return None


def _validate_node_id(node_id: str) -> bool:
    """Validate that *node_id* matches ``[A-Za-z]+[0-9]+``."""
    chars = list(node_id)
    if not chars[0].isalpha():
        raise NotationException(f"{INVALID_NODE_ID}: {node_id}")
    if not chars[-1].isdigit():
        raise NotationException(f"{INVALID_NODE_ID}: {node_id}")
    found_num = False
    for ch in chars:
        if ch.isdigit():
            found_num = True
        if found_num and ch.isalpha():
            raise NotationException(f"{INVALID_NODE_ID}: {node_id}")
    return True


# ---------------------------------------------------------------------------
# Edge-list parsing
# ---------------------------------------------------------------------------

def _get_polymer_edge_list(all_edge_string: str) -> list[PolymerEdge]:
    """Convert *all_edge_string* to a list of :class:`PolymerEdge`."""
    result: list[PolymerEdge] = []
    if all_edge_string:
        all_edges = re.split(LIST_LEVEL_DELIMITER_REGEX, all_edge_string)
        for edge_str in all_edges:
            edge = edge_parser.parse(edge_str)
            result.append(edge)
    return result


# ---------------------------------------------------------------------------
# RNA helpers
# ---------------------------------------------------------------------------

def get_rna_polymer_node_list(
    complex_notation: str,
    monomer_store: MonomerStore | None = None,
) -> list[RNAPolymerNode]:
    """Return the :class:`RNAPolymerNode` entries found in the notation."""
    from helm_notation.tools import simple_notation_parser

    monomer_store = _check_for_monomer_store(monomer_store)
    node_list = get_polymer_node_list(get_all_node_string(complex_notation))
    all_node_annotation_string = get_all_node_label_string(complex_notation)
    annotation_map = _get_polymer_node_id_annotation_map(all_node_annotation_string)
    result: list[RNAPolymerNode] = []
    for node in node_list:
        if node.type == NUCLEIC_ACID_POLYMER_TYPE:
            rna_node = RNAPolymerNode.from_polymer_node(node)
            nuc_list = simple_notation_parser.get_nucleotide_list(
                node.label, monomer_store=monomer_store
            )
            seq = simple_notation_parser.get_nucleotide_sequence(nuc_list)
            rna_node.sequence = seq
            mod_seq = simple_notation_parser.get_modified_nucleotide_sequence(nuc_list)
            rna_node.modified_sequence = mod_seq
            if rna_node.id in annotation_map:
                rna_node.anotation = annotation_map[rna_node.id]
            result.append(rna_node)
    return result


# ---------------------------------------------------------------------------
# siRNA formatting
# ---------------------------------------------------------------------------

def get_formated_sirna_sequences(
    complex_notation: str,
    padding_char: str = DEFAULT_PADDING_CHAR,
    base_pair_char: str = DEFAULT_BASE_PAIR_CHAR,
) -> list[str]:
    """Generate formatted siRNA sequences with padding and base-pair chars.

    Returns a list of strings (0, 1, or 3 elements).
    """
    from helm_notation.tools import nucleotide_sequence_parser

    if padding_char is None or len(padding_char) != 1:
        raise NotationException("Padding string must be single character")
    if base_pair_char is None or len(base_pair_char) != 1:
        raise NotationException("Base pair string must be single character")

    rna_list = get_rna_polymer_node_list(complex_notation)
    count = len(rna_list)
    if count == 0:
        return []
    elif count == 1:
        return [rna_list[0].sequence]
    elif count == 2:
        rna1_seq: str | None = None
        rna2_seq: str | None = None
        rna1_annotation: str | None = None
        rna2_annotation: str | None = None

        for node in rna_list:
            if node.id == "RNA1":
                rna1_seq = node.sequence
                rna1_annotation = node.anotation
            elif node.id == "RNA2":
                rna2_seq = node.sequence
                rna2_annotation = node.anotation

        reverse_rna2_seq = nucleotide_sequence_parser.get_reverse_sequence(rna2_seq)

        bp_string = get_all_base_pair_string(complex_notation)
        if not bp_string:
            return [rna1_seq, rna2_seq]
        else:
            monomer_position_map = _get_sirna_monomer_position_map(bp_string)
            seq_position_map: dict[int, int] = {}
            for key, value in monomer_position_map.items():
                seq_key = key // 3 + 1
                seq_value = value // 3 + 1
                seq_position_map[seq_key] = seq_value

            seq_list = sorted(seq_position_map.keys())
            rna1_first = seq_list[0]
            rna2_last = seq_position_map[seq_list[0]]
            rna1_last = seq_list[-1]
            rna2_first = seq_position_map[seq_list[-1]]

            if (rna1_last - rna1_first) != (rna2_last - rna2_first):
                raise NotationException("siRNA matching lengths are different")

            rna1_left_overhang = rna1_first - 1
            rna1_right_overhang = len(rna1_seq) - rna1_last
            rna2_left_overhang = len(rna2_seq) - rna2_last
            rna2_right_overhang = rna2_first - 1
            sbs: list[list[str]] = [[], [], []]

            if rna1_left_overhang >= rna2_left_overhang:
                sbs[0].append(rna1_seq)

                for _ in range(rna1_left_overhang):
                    sbs[1].append(padding_char)
                for i in range(rna1_first, rna1_last + 1):
                    if i in seq_position_map:
                        sbs[1].append(base_pair_char)
                    else:
                        sbs[1].append(padding_char)

                for _ in range(rna1_left_overhang - rna2_left_overhang):
                    sbs[2].append(padding_char)
                sbs[2].append(reverse_rna2_seq)
            else:
                for _ in range(rna2_left_overhang - rna1_left_overhang):
                    sbs[0].append(padding_char)
                sbs[0].append(rna1_seq)

                for _ in range(rna2_left_overhang):
                    sbs[1].append(padding_char)
                for i in range(rna1_first, rna1_last + 1):
                    if i in seq_position_map:
                        sbs[1].append(base_pair_char)
                    else:
                        sbs[1].append(padding_char)

                sbs[2].append(reverse_rna2_seq)

            if rna1_right_overhang >= rna2_right_overhang:
                for _ in range(rna1_right_overhang):
                    sbs[1].append(padding_char)
                for _ in range(rna1_right_overhang - rna2_right_overhang):
                    sbs[2].append(padding_char)
            else:
                for _ in range(rna2_right_overhang - rna1_right_overhang):
                    sbs[0].append(padding_char)
                for _ in range(rna2_right_overhang - rna1_right_overhang):
                    sbs[1].append(padding_char)

            row0 = "".join(sbs[0])
            row1 = "".join(sbs[1])
            row2 = "".join(sbs[2])

            if (
                (rna1_annotation and rna1_annotation.upper() == "AS")
                or (rna2_annotation and rna2_annotation.upper() == "SS")
            ):
                return [row2[::-1], row1[::-1], row0[::-1]]
            else:
                return [row0, row1, row2]
    else:
        raise NotationException(
            "Structure contains more than two RNA sequences"
        )


def _get_sirna_monomer_position_map(base_pair_string: str) -> dict[int, int]:
    pairs = re.split(LIST_LEVEL_DELIMITER_REGEX, base_pair_string)
    result: dict[int, int] = {}
    for pair in pairs:
        edge_components = pair.split(EDGE_COMPONENT_DELIMITER)
        pair_info = edge_components[2]
        edge_properties = pair_info.split(CONNECTION_SEPARATOR)

        str1 = edge_properties[0].split(MONOMER_ATTACHMENT_SEPARATOR)[0]
        int1 = int(str1)

        str2 = edge_properties[1].split(MONOMER_ATTACHMENT_SEPARATOR)[0]
        int2 = int(str2)

        if pair.startswith("RNA1"):
            result[int1] = int2
        elif pair.startswith("RNA2"):
            result[int2] = int1
        else:
            raise NotationException(
                "Structure contains more than two RNA sequences"
            )
    return result


def _get_polymer_node_id_annotation_map(
    all_node_annotation_string: str,
) -> dict[str, str]:
    result: dict[str, str] = {}
    if all_node_annotation_string and all_node_annotation_string.strip():
        items = re.split(LIST_LEVEL_DELIMITER_REGEX, all_node_annotation_string)
        for item in items:
            start = item.index(NODE_LABEL_START_SYMBOL)
            end = item.index(NODE_LABEL_END_SYMBOL)
            node_id = item[:start]
            annotation = item[start + 1 : end]
            result[node_id] = annotation
    return result


# ---------------------------------------------------------------------------
# Canonical notation
# ---------------------------------------------------------------------------

def get_canonical_notation(
    complex_notation: str,
    include_validation: bool = False,
    monomer_store: MonomerStore | None = None,
) -> str:
    """Generate canonical notation for *complex_notation*.

    When *include_validation* is ``True`` the notation is also validated.
    """
    from helm_notation.tools import simple_notation_parser
    from helm_notation.tools import permutation_expansion

    monomer_store = _check_for_monomer_store(monomer_store)

    validate_notation_format(complex_notation)
    cp = parse(complex_notation, monomer_store)

    if include_validation:
        validate_complex_polymer(cp, monomer_store)

    node_list = cp.polymer_node_list
    edge_list = cp.polymer_edge_list

    # Deal with ad-hoc CHEM monomers -- use SMILES instead of temp ID
    for node in node_list:
        if (
            node.type == CHEMICAL_POLYMER_TYPE
            and node.label.startswith(
                simple_notation_parser.get_ad_hoc_monomer_id_prefix(
                    CHEMICAL_POLYMER_TYPE
                )
            )
        ):
            m = monomer_store.get_monomer(CHEMICAL_POLYMER_TYPE, node.label)
            smiles = m.can_smiles
            uni_smi = structure_parser.get_unique_extended_smiles(smiles)

            rgroups = structure_parser.get_rgroups_from_extended_smiles(uni_smi)
            old_new_rgroup_map: dict[str, str] = {}
            for i, old_id in enumerate(rgroups):
                new_id = f"Q{i + 1}"
                old_new_rgroup_map[old_id] = new_id
                uni_smi = uni_smi.replace(old_id, new_id)
            uni_smi = uni_smi.replace("Q", "R")
            node.label = uni_smi

            node_id = node.id

            # Deal with edges
            for edge in edge_list:
                if edge.source_node == node_id:
                    old_r = edge.source_r
                    new_r = old_new_rgroup_map.get(old_r)
                    if new_r is not None:
                        new_r = new_r.replace("Q", "R")
                        edge.source_r = new_r

                if edge.target_node == node_id:
                    old_r = edge.target_r
                    new_r = old_new_rgroup_map.get(old_r)
                    if new_r is not None:
                        new_r = new_r.replace("Q", "R")
                        edge.target_r = new_r

    # Identify all backbone-cycle and branch-cycle nodes
    backbone_cycle_nodes: list[str] = []
    branch_cycle_nodes: list[str] = []
    for edge in edge_list:
        if edge.is_backbone_self_cycle:
            backbone_cycle_nodes.append(edge.source_node)
        elif edge.is_self_cycle:
            branch_cycle_nodes.append(edge.source_node)

    # Deal with backbone cycles
    for node_id in backbone_cycle_nodes:
        offset = 0
        monomer_count = 0
        for node in node_list:
            if node.id == node_id:
                can_entry = simple_notation_parser.get_simple_canonical_notation_map_entry(
                    node.label, node.type, monomer_store
                )
                node.label = can_entry[1]
                offset = can_entry[0]
                monomer_count = simple_notation_parser.get_monomer_count(
                    node.label, node.type, monomer_store
                )
                break

        for edge in edge_list:
            if (
                edge.source_node == node_id
                or edge.target_node == node_id
            ):
                if edge.is_backbone_self_cycle:
                    connections = sorted(
                        [edge.connection, edge.reverse_connection]
                    )
                    edge.connection = connections[0]
                else:
                    if edge.source_node == node_id:
                        old_src = edge.source_monomer_number
                        if old_src > offset:
                            new_src = old_src - offset
                        else:
                            new_src = monomer_count + old_src - offset
                        edge.source_monomer_number = new_src

                    if edge.target_node == node_id:
                        old_tar = edge.target_monomer_number
                        if old_tar > offset:
                            new_tar = old_tar - offset
                        else:
                            new_tar = monomer_count + old_tar - offset
                        edge.target_monomer_number = new_tar

    # Deal with branch self-cycles
    for node_id in branch_cycle_nodes:
        for edge in edge_list:
            if (
                edge.is_self_cycle
                and not edge.is_backbone_self_cycle
                and edge.source_node == node_id
            ):
                connections = sorted(
                    [edge.connection, edge.reverse_connection]
                )
                edge.connection = connections[0]

    id_label_map: dict[str, str] = {}
    label_id_map: dict[str, list[str]] = {}  # sorted iteration via TreeMap equivalent

    for node in node_list:
        id_label_map[node.id] = node.label
        if node.label in label_id_map:
            label_id_map[node.label].append(node.id)
        else:
            label_id_map[node.label] = [node.id]

    # Sort by label
    sorted_labels = sorted(label_id_map.keys())
    lol: list[list[list[str]]] = []
    for key in sorted_labels:
        value = label_id_map[key]
        al: list[list[str]] = [value]
        permutation_expansion.expand(lol, al)

    node_id_permutations = permutation_expansion.linearize(lol)
    notation_list: list[str] = []
    for sorted_id_list in node_id_permutations:
        notation = _generate_notation_based_node_order(
            sorted_id_list, node_list, edge_list
        )
        notation_list.append(notation)
    notation_list.sort()
    return notation_list[0]


def _generate_notation_based_node_order(
    sorted_node_id_list: list[str],
    node_list: list[PolymerNode],
    edge_list: list[PolymerEdge],
) -> str:
    id_label_map: dict[str, str] = {}
    for node in node_list:
        id_label_map[node.id] = node.label

    type_count_map: dict[str, int] = {}
    old_new_id_map: dict[str, str] = {}

    sb_parts: list[str] = []
    for old_id in sorted_node_id_list:
        label = id_label_map[old_id]
        ptype = re.split(r"\d", old_id)[0]
        if ptype in type_count_map:
            type_count_map[ptype] += 1
        else:
            type_count_map[ptype] = 1
        counter = type_count_map[ptype]
        new_id = f"{ptype}{counter}"
        old_new_id_map[old_id] = new_id
        sb_parts.append(f"{new_id}{{{label}}}")
    sb = LIST_LEVEL_DELIMITER.join(sb_parts) + TOP_LEVEL_DELIMITER

    if edge_list:
        for edge in edge_list:
            edge_parser.canonicalize(edge, sorted_node_id_list)

        sorted_edge_list = edge_parser.sort(edge_list, sorted_node_id_list)

        for edge in sorted_edge_list:
            old_source = edge.source_node
            new_source = _replace_node_id(old_source, old_new_id_map)
            edge.source_node = new_source

            old_target = edge.target_node
            new_target = _replace_node_id(old_target, old_new_id_map)
            edge.target_node = new_target

        edge_parts: list[str] = []
        for edge in sorted_edge_list:
            edge_parts.append(str(edge))
        sb += LIST_LEVEL_DELIMITER.join(edge_parts)

    sb += "$$$"
    return sb


def _replace_node_id(
    old_node_id: str,
    old_new_id_map: dict[str, str],
) -> str:
    old_nodes = re.split(NODE_CONCATENATOR_REGEX, old_node_id)
    new_nodes = [old_new_id_map[n] for n in old_nodes]
    return NODE_CONCATENATOR.join(new_nodes)


# ---------------------------------------------------------------------------
# Combining notations
# ---------------------------------------------------------------------------

def get_combined_complex_notation(
    complex_notation1: str,
    complex_notation2: str,
    monomer_store: MonomerStore | None = None,
) -> str:
    """Combine two complex notations with no chemical connection or base pairing."""
    monomer_store = _check_for_monomer_store(monomer_store)
    if not complex_notation1 or not complex_notation1.strip():
        return complex_notation2
    if not complex_notation2 or not complex_notation2.strip():
        return complex_notation1

    rna_count = 0
    peptide_count = 0
    chem_count = 0

    node_string1 = get_all_node_string(complex_notation1)
    edge_string1 = get_all_edge_string(complex_notation1)
    base_pair_string1 = get_all_base_pair_string(complex_notation1)
    node_label_string1 = get_all_node_label_string(complex_notation1)
    node_list1 = get_polymer_node_list(node_string1, monomer_store)
    for node in node_list1:
        ptype = node.type
        if ptype == NUCLEIC_ACID_POLYMER_TYPE:
            rna_count += 1
        elif ptype == PEPTIDE_POLYMER_TYPE:
            peptide_count += 1
        elif ptype == CHEMICAL_POLYMER_TYPE:
            chem_count += 1
        else:
            raise NotationException(f"Unknown polymer type found [{ptype}]")

    node_id_map: dict[str, str] = {}
    node_string2 = get_all_node_string(complex_notation2)
    edge_string2 = get_all_edge_string(complex_notation2)
    base_pair_string2 = get_all_base_pair_string(complex_notation2)
    node_label_string2 = get_all_node_label_string(complex_notation2)
    node_list2 = get_polymer_node_list(node_string2, monomer_store)
    for node in node_list2:
        ptype = node.type
        old_id = node.id
        new_id: str | None = None
        if ptype == NUCLEIC_ACID_POLYMER_TYPE:
            rna_count += 1
            new_id = f"{ptype}{rna_count}"
        elif ptype == PEPTIDE_POLYMER_TYPE:
            peptide_count += 1
            new_id = f"{ptype}{peptide_count}"
        elif ptype == CHEMICAL_POLYMER_TYPE:
            chem_count += 1
            new_id = f"{ptype}{chem_count}"
        else:
            raise NotationException(f"Unknown polymer type found [{ptype}]")
        node_id_map[old_id] = new_id

    for old_id, new_id in node_id_map.items():
        node_string2 = node_string2.replace(old_id, new_id)
        edge_string2 = edge_string2.replace(old_id, new_id)
        base_pair_string2 = base_pair_string2.replace(old_id, new_id)
        node_label_string2 = node_label_string2.replace(old_id, new_id)

    sb: list[str] = []
    # Build all nodes
    sb.append(node_string1)
    sb.append("|")
    sb.append(node_string2)
    sb.append("$")

    # Build all edges
    sb.append(edge_string1)
    if edge_string1:
        sb.append("|")
    sb.append(edge_string2)
    sb.append("$")

    # Build all base pairs
    sb.append(base_pair_string1)
    if base_pair_string1:
        sb.append("|")
    sb.append(base_pair_string2)
    sb.append("$")

    # Build all node labels
    sb.append(node_label_string1)
    if node_label_string1:
        sb.append("|")
    sb.append(node_label_string2)
    sb.append("$")

    return "".join(sb)


# ---------------------------------------------------------------------------
# Hybridize
# ---------------------------------------------------------------------------

def hybridize(complex_notation: str) -> str:
    """Auto-add base pair info when the notation has exactly two RNA nodes
    and no existing base-pairing.
    """
    rna_list = get_rna_polymer_node_list(complex_notation)
    base_pair_string = get_all_base_pair_string(complex_notation)
    if len(rna_list) == 2 and not base_pair_string:
        node_string = get_all_node_string(complex_notation)
        edge_string = get_all_edge_string(complex_notation)
        label_string = get_all_node_label_string(complex_notation)
        base_pair_string = _get_base_pair_string(rna_list[0], rna_list[1])
        return f"{node_string}${edge_string}${base_pair_string}${label_string}$"
    else:
        return complex_notation


def _get_base_pair_string(
    node1: RNAPolymerNode,
    node2: RNAPolymerNode,
    monomer_store: MonomerStore | None = None,
) -> str:
    from helm_notation.tools import simple_notation_parser
    from helm_notation.tools import nucleotide_sequence_parser

    monomer_store = _check_for_monomer_store(monomer_store)

    seq1 = node1.sequence.replace("T", "U")
    seq2 = node2.sequence.replace("T", "U")
    # Reverse complement of seq2
    chars = list(seq2)
    comp_parts: list[str] = []
    for i in range(len(chars) - 1, -1, -1):
        symbol = chars[i]
        comp_symbol = nucleotide_sequence_parser.complement_map.get(symbol, symbol)
        comp_parts.append(comp_symbol)
    comp_seq2 = "".join(comp_parts)

    max_seq_match = nucleotide_sequence_parser.get_max_match_fragment(seq1, comp_seq2)
    seq_match_length = len(max_seq_match)
    seq1_nuc_start = -1
    seq1_monomer_start = 0
    seq2_nuc_start = -1
    seq2_monomer_start = 0

    seq1_nuc_list = simple_notation_parser.get_nucleotide_list(
        node1.label, include_last=False, monomer_store=monomer_store
    )
    seq2_nuc_list = simple_notation_parser.get_nucleotide_list(
        node2.label, include_last=False, monomer_store=monomer_store
    )

    base_pair = ""
    if seq_match_length > 0:
        # Get starting monomer position for sequence 1
        seq1_nuc_start = seq1.index(max_seq_match)
        for i in range(seq1_nuc_start):
            nuc = seq1_nuc_list[i]
            mc = simple_notation_parser.get_monomer_count_for_rna(
                nuc.notation, monomer_store
            )
            seq1_monomer_start += mc

        # Get starting monomer position for sequence 2
        comp_seq2_nuc_start = comp_seq2.index(max_seq_match)
        seq2_nuc_start = len(seq2) - seq_match_length - comp_seq2_nuc_start
        for i in range(seq2_nuc_start):
            nuc = seq2_nuc_list[i]
            mc = simple_notation_parser.get_monomer_count_for_rna(
                nuc.notation, monomer_store
            )
            seq2_monomer_start += mc

        # Build matching monomer positions
        for i in range(seq_match_length):
            nuc1 = seq1_nuc_list[i + seq1_nuc_start]
            if nuc1.get_base_monomer(monomer_store) is None:
                raise NotationException(
                    "Nucleotide without base cannot be hybridized with others"
                )
            if i == 0:
                seq1_monomer_start += 2
            else:
                seq1_monomer_start += 3

            nuc2 = seq2_nuc_list[i + seq2_nuc_start]
            if nuc2.get_base_monomer(monomer_store) is None:
                raise NotationException(
                    "Nucleotide without base cannot be hybridized with others"
                )
            if i == 0:
                seq2_monomer_start += 2
            else:
                seq2_monomer_start += 3

        # Build base pair string, offset by 3
        for i in range(seq_match_length, 0, -1):
            seq1_monomer_pos = seq1_monomer_start - (i - 1) * 3
            seq2_monomer_pos = seq2_monomer_start - (seq_match_length - i) * 3

            if base_pair:
                base_pair += "|"
            base_pair += (
                f"RNA1,RNA2,{seq1_monomer_pos}:pair-{seq2_monomer_pos}:pair"
            )

    return base_pair


# ---------------------------------------------------------------------------
# Decompose
# ---------------------------------------------------------------------------

def decompose(
    complex_notation: str,
    monomer_store: MonomerStore | None = None,
) -> list[str]:
    """Decompose *complex_notation* into notations of covalently connected
    monomers. Base-pairing info is lost.
    """
    monomer_store = _check_for_monomer_store(monomer_store)
    cp = parse(complex_notation, monomer_store)
    annotation_map = cp.polymer_node_annotation_map

    grouped_nodes = _get_group_node_list_map(
        cp.polymer_node_list, cp.polymer_edge_list
    )
    if len(grouped_nodes) == 1:
        return [complex_notation]

    group_notation_map: dict[int, str] = {}
    for group_id, all_nodes in grouped_nodes.items():
        # Remove duplicated nodes
        connected_nodes: list[str] = []
        for n in all_nodes:
            if n not in connected_nodes:
                connected_nodes.append(n)

        connected_edges = _get_connected_edge_list(
            connected_nodes, cp.polymer_edge_list
        )

        node_sb_parts: list[str] = []
        for nid in connected_nodes:
            poly_notation = _get_polymer_notation(nid, cp.polymer_node_list)
            node_sb_parts.append(f"{nid}{{{poly_notation}}}")

        edge_sb_parts: list[str] = []
        for edge in connected_edges:
            edge_sb_parts.append(str(edge))

        annotation_sb_parts: list[str] = []
        for nid in connected_nodes:
            if nid in annotation_map:
                annotation_sb_parts.append(f"{nid}{{{annotation_map[nid]}}}")

        result_notation = (
            LIST_LEVEL_DELIMITER.join(node_sb_parts)
            + "$"
            + LIST_LEVEL_DELIMITER.join(edge_sb_parts)
            + "$$"
            + LIST_LEVEL_DELIMITER.join(annotation_sb_parts)
            + "$"
        )
        group_notation_map[group_id] = result_notation

    # Sort based on polymerNode list order
    results: list[str] = []
    used_nodes: list[str] = []
    for polymer_node in cp.polymer_node_list:
        nid = polymer_node.id
        if nid not in used_nodes:
            for group_id, nodes in grouped_nodes.items():
                if nid in nodes:
                    notation = group_notation_map[group_id]
                    results.append(notation)
                    used_nodes.extend(nodes)
                    break

    return results


# ---------------------------------------------------------------------------
# Replace monomer
# ---------------------------------------------------------------------------

def replace_monomer(
    complex_notation: str,
    polymer_type: str,
    existing_monomer_id: str,
    new_monomer_id: str,
    monomer_store: MonomerStore | None = None,
    validate: bool = True,
) -> str:
    """Replace *existing_monomer_id* with *new_monomer_id* in the notation."""
    from helm_notation.tools import simple_notation_parser

    monomer_store = _check_for_monomer_store(monomer_store)

    if validate:
        simple_notation_parser.validate_monomer_replacement(
            polymer_type, existing_monomer_id, new_monomer_id, monomer_store
        )
    if not complex_notation:
        return complex_notation
    if existing_monomer_id == new_monomer_id:
        return complex_notation

    all_node_string = get_all_node_string(complex_notation)
    rest_of_notation = complex_notation[len(all_node_string):]
    polymers = get_polymer_node_list(all_node_string, monomer_store)
    sb_parts: list[str] = []
    for polymer in polymers:
        parts: list[str] = [polymer.id, NODE_LABEL_START_SYMBOL]
        if polymer.type == polymer_type:
            notation = polymer.label
            result = simple_notation_parser.replace_monomer(
                notation,
                polymer_type,
                existing_monomer_id,
                new_monomer_id,
                monomer_store,
                validate,
            )
            parts.append(result)
        else:
            parts.append(polymer.label)
        parts.append(NODE_LABEL_END_SYMBOL)
        sb_parts.append("".join(parts))

    return LIST_LEVEL_DELIMITER.join(sb_parts) + rest_of_notation


# ---------------------------------------------------------------------------
# Standardize
# ---------------------------------------------------------------------------

def standardize(
    notation: str,
    monomer_store: MonomerStore | None = None,
) -> str:
    """Standardize all forms of HELM notation into complex notation."""
    from helm_notation.tools import simple_notation_parser
    from helm_notation.monomer_factory import get_instance as get_monomer_factory

    if not notation:
        return notation

    if monomer_store is None:
        monomer_store = get_monomer_factory().monomer_store

    # If it has $ delimiters, it's already complex notation
    if TOP_LEVEL_DELIMITER in notation:
        validate_complex_notation(notation, monomer_store)
        return notation
    else:
        factory = get_monomer_factory()
        polymer_types = sorted(factory.monomer_db.keys())
        for polymer_type in polymer_types:
            try:
                simple_notation_parser.validate_simple_notation(
                    notation, polymer_type, monomer_store
                )
                return simple_notation_parser.get_complex_notation(
                    notation, polymer_type, monomer_store
                )
            except Exception:
                pass
        raise NotationException(
            f"Invalid HELM notation for simple polymer: {notation}"
        )


# ---------------------------------------------------------------------------
# Nucleotide modification check
# ---------------------------------------------------------------------------

def has_nucleotide_modification(
    polymer_nodes: list[PolymerNode],
    monomer_store: MonomerStore | None = None,
) -> bool:
    """Return ``True`` if any RNA node contains a modified nucleotide."""
    from helm_notation.tools import simple_notation_parser

    monomer_store = _check_for_monomer_store(monomer_store)

    for node in polymer_nodes:
        one_node_id = node.id
        notation = node.label

        if one_node_id.startswith(CHEMICAL_POLYMER_TYPE):
            pass
        elif one_node_id.startswith(PEPTIDE_POLYMER_TYPE):
            pass
        else:
            nucleotide_list = simple_notation_parser.get_strict_nucleotide_list(
                notation, include_last=False, monomer_store=monomer_store
            )
            for nuc in nucleotide_list:
                if nuc.is_modified:
                    return True
    return False


# ---------------------------------------------------------------------------
# Total monomer count
# ---------------------------------------------------------------------------

def get_total_monomer_count(
    notation: str,
    monomer_store: MonomerStore | None = None,
) -> int:
    """Return the total monomer count across all polymer nodes."""
    from helm_notation.tools import simple_notation_parser

    monomer_store = _check_for_monomer_store(monomer_store)

    total = 0
    nodes = get_polymer_node_list(notation, monomer_store)
    for node in nodes:
        polymer_type = node.type
        label = node.label
        mc = simple_notation_parser.get_monomer_count(
            label, polymer_type, monomer_store
        )
        total += mc
    return total


# ---------------------------------------------------------------------------
# MoleculeInfo calculation
# ---------------------------------------------------------------------------

def get_molecule_info(
    extended_notation: str,
    include_validation: bool = False,
    monomer_store: MonomerStore | None = None,
) -> MoleculeInfo:
    """Return :class:`MoleculeInfo` using the divide-and-conquer approach.

    Ignores fuzzy / generic edges in the calculation (good enough for large
    structures).
    """
    from helm_notation.tools import simple_notation_parser

    monomer_store = _check_for_monomer_store(monomer_store)

    cp = parse(extended_notation, monomer_store)

    if include_validation:
        validate_complex_polymer(cp, monomer_store)

    node_list = cp.polymer_node_list
    edge_list = cp.polymer_edge_list

    # Deal with all connections between polymers
    caps: list[MoleculeInfo] = []
    for edge in edge_list:
        if edge.edge_type == GENERIC_EDGE:
            continue

        source_id = edge.source_node
        source_node = _get_polymer_node_from_id(node_list, source_id)
        source_monomer_number = edge.source_monomer_number
        source_r = edge.source_r
        source_cap_mi = _get_cap_molecule_info(
            source_node, source_monomer_number, source_r, monomer_store
        )
        caps.append(source_cap_mi)

        target_id = edge.target_node
        target_node = _get_polymer_node_from_id(node_list, target_id)
        target_monomer_number = edge.target_monomer_number
        target_r = edge.target_r
        target_cap_mi = _get_cap_molecule_info(
            target_node, target_monomer_number, target_r, monomer_store
        )
        caps.append(target_cap_mi)

    # Deal with all simple polymers
    chunks: list[MoleculeInfo] = []
    for node in node_list:
        polymer_type = node.type
        label = node.label
        tmp_mi = simple_notation_parser.get_molecule_info_for_polymer(
            label, polymer_type, monomer_store
        )
        chunks.append(tmp_mi)

    return structure_parser.process_molecule_info(chunks, caps)


def _get_cap_molecule_info(
    node: PolymerNode,
    monomer_number: int,
    rgroup: str,
    monomer_store: MonomerStore | None = None,
) -> MoleculeInfo:
    from helm_notation.tools import simple_notation_parser

    monomer_store = _check_for_monomer_store(monomer_store)
    polymer_type = node.type
    notation = node.label

    monomer_id_list = simple_notation_parser.get_monomer_id_list(
        notation, polymer_type, monomer_store
    )
    monomer = simple_notation_parser.get_monomer(
        monomer_id_list[monomer_number - 1], polymer_type, monomer_store
    )
    return monomer.get_cap_molecule_info(rgroup)


def _get_polymer_node_from_id(
    nodes: list[PolymerNode],
    node_id: str,
) -> PolymerNode | None:
    for node in nodes:
        if node.id == node_id:
            return node
    return None


def _get_polymer_node_structure_map(
    node_list: list[PolymerNode],
    monomer_store: MonomerStore | None,
) -> dict[str, RgroupStructure]:
    from helm_notation.tools import simple_notation_parser

    monomer_store = _check_for_monomer_store(monomer_store)

    node_struc_map: dict[str, RgroupStructure] = {}
    for node in node_list:
        node_id = node.id
        node_label = node.label
        polymer_type = get_polymer_type(node_id)
        struc = simple_notation_parser.get_simple_polymer_structure(
            node_label, polymer_type, monomer_store
        )
        if struc.molecule is None:
            raise NotationException(
                "Polymer notation contains non-specific monomer structure"
            )
        node_struc_map[node_id] = struc
    return node_struc_map


# ---------------------------------------------------------------------------
# SMILES replacement helper
# ---------------------------------------------------------------------------

def get_notation_by_replacing_smiles(
    helm_string: str,
    monomer_store: MonomerStore,
) -> str:
    """Replace SMILES in complex notation with temporary IDs."""
    from helm_notation.tools import simple_notation_parser

    all_node_string = get_all_node_string(helm_string)
    rest_of_notation = helm_string[len(all_node_string):]
    polymers = get_polymer_node_list(all_node_string, monomer_store)
    sb_parts: list[str] = []
    for polymer in polymers:
        parts: list[str] = [polymer.id, NODE_LABEL_START_SYMBOL]

        notation = polymer.label
        result = simple_notation_parser.get_notation_by_replacing_smiles(
            notation, polymer.type, monomer_store
        )
        parts.append(result)

        parts.append(NODE_LABEL_END_SYMBOL)
        sb_parts.append("".join(parts))

    return LIST_LEVEL_DELIMITER.join(sb_parts) + rest_of_notation
