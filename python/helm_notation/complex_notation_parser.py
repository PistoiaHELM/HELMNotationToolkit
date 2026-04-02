"""Phase 3 port of core org.helm.notation.tools.ComplexNotationParser APIs.

Current scope:
- standard complex notation with node list + edge list
- canonicalization for common cyclic peptide case
- complex polymer smiles and molecule info

This is an incremental implementation aimed at baseline parity checkpoints.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from rdkit import Chem

from helm_notation.chem import ChemEngine
from helm_notation.exceptions import NotationException, StructureException
from helm_notation.model import MoleculeInfo
from helm_notation.monomer_factory import MonomerFactory
from helm_notation.monomer_store import MonomerStore
from helm_notation.simple_notation_parser import (
    SimpleNotationParser,
    _find_r_atom,
    _merge_mols_on_r_atoms,
    cyclize_on_r_atoms,
)


@dataclass
class _Node:
    node_id: str
    polymer_type: str
    label: str


@dataclass
class _Edge:
    source: str
    target: str
    connection: str


class ComplexNotationParser:
    TOP_LEVEL_DELIMITER = "$"
    LIST_LEVEL_DELIMITER = "|"

    @staticmethod
    def get_complex_polymer_smiles(
        extended_notation: str,
        monomer_store: Optional[MonomerStore] = None,
    ) -> str:
        store = monomer_store or MonomerFactory.get_instance().get_monomer_store()
        components = _build_connected_components(extended_notation, store)

        capped_components: List[str] = []
        for mol in components:
            capped = _cap_remaining_rgroups(mol, store)
            capped_components.append(ChemEngine.get_unique_smiles(capped))

        if not capped_components:
            raise NotationException("No structure generated from complex notation")

        mixture = ".".join(capped_components)
        return ChemEngine.get_unique_smiles(mixture)

    @staticmethod
    def get_molecule_info(
        extended_notation: str,
        monomer_store: Optional[MonomerStore] = None,
    ) -> MoleculeInfo:
        smiles = ComplexNotationParser.get_complex_polymer_smiles(
            extended_notation, monomer_store
        )
        return ChemEngine.get_molecule_info(smiles)

    @staticmethod
    def get_canonical_notation(
        extended_notation: str,
        monomer_store: Optional[MonomerStore] = None,
    ) -> str:
        nodes, edges, bp, annotations, other = _parse_extended_notation(extended_notation)

        # Incremental canonicalization: rotate cyclic peptide labels to lexicographically
        # minimal rotation when a self-cycle edge 1:R1-n:R2 exists.
        edge_set = {(e.source, e.target, e.connection) for e in edges}

        canon_nodes: List[_Node] = []
        for node in nodes:
            if node.polymer_type == "PEPTIDE":
                tokens = [t for t in node.label.split(".") if t]
                n = len(tokens)
                if n > 1:
                    has_cycle = (
                        (node.node_id, node.node_id, f"1:R1-{n}:R2") in edge_set
                        or (node.node_id, node.node_id, f"{n}:R2-1:R1") in edge_set
                    )
                    if has_cycle:
                        rotations = [tokens[i:] + tokens[:i] for i in range(n)]
                        best = min(rotations)
                        node = _Node(node.node_id, node.polymer_type, ".".join(best))
            canon_nodes.append(node)

        return _render_extended_notation(canon_nodes, edges, bp, annotations, other)


# ---------------------------------------------------------------------------
# Parsing / rendering
# ---------------------------------------------------------------------------


def _parse_extended_notation(
    extended_notation: str,
) -> Tuple[List[_Node], List[_Edge], str, str, str]:
    parts = extended_notation.split("$")
    while len(parts) < 5:
        parts.append("")

    nodes_part, edges_part, basepair_part, annotations_part, other_part = parts[:5]

    nodes: List[_Node] = []
    for token in [t for t in nodes_part.split("|") if t]:
        if "{" not in token or not token.endswith("}"):
            raise NotationException(f"Invalid polymer node notation: {token}")
        node_id, label = token.split("{", 1)
        label = label[:-1]
        polymer_type = _polymer_type(node_id)
        nodes.append(_Node(node_id=node_id, polymer_type=polymer_type, label=label))

    edges: List[_Edge] = []
    if edges_part:
        for token in [t for t in edges_part.split("|") if t]:
            comps = token.split(",")
            if len(comps) != 3:
                raise NotationException(f"Invalid edge notation: {token}")
            edges.append(_Edge(comps[0], comps[1], comps[2]))

    return nodes, edges, basepair_part, annotations_part, other_part


def _render_extended_notation(
    nodes: List[_Node],
    edges: List[_Edge],
    basepair_part: str,
    annotations_part: str,
    other_part: str,
) -> str:
    nodes_part = "|".join([f"{n.node_id}{{{n.label}}}" for n in nodes])
    edges_part = "|".join([f"{e.source},{e.target},{e.connection}" for e in edges])
    return "$".join([nodes_part, edges_part, basepair_part, annotations_part, other_part])


def _polymer_type(node_id: str) -> str:
    out = []
    for ch in node_id:
        if ch.isdigit():
            break
        out.append(ch)
    return "".join(out)


# ---------------------------------------------------------------------------
# Structure generation
# ---------------------------------------------------------------------------


def _build_connected_components(
    extended_notation: str,
    monomer_store: MonomerStore,
) -> List[Chem.Mol]:
    nodes, edges, _, _, _ = _parse_extended_notation(extended_notation)

    # One molecule per node initially
    components: Dict[str, Chem.Mol] = {}
    node_component: Dict[str, str] = {}

    for node in nodes:
        structure = SimpleNotationParser.get_simple_polymer_structure(
            node.label,
            node.polymer_type,
            monomer_store,
            node_id=node.node_id,
        )
        if structure.molecule is None:
            raise NotationException(f"Unable to build structure for node {node.node_id}")
        components[node.node_id] = structure.molecule
        node_component[node.node_id] = node.node_id

    for edge in edges:
        if "generic" in edge.connection or "pair" in edge.connection:
            raise NotationException(
                "Phase 3 parser currently supports standard edges only"
            )

        left, right = edge.connection.split("-")
        src_pos_str, src_r = left.split(":")
        tgt_pos_str, tgt_r = right.split(":")
        src_pos = int(src_pos_str)
        tgt_pos = int(tgt_pos_str)

        src_comp = node_component[edge.source]
        tgt_comp = node_component[edge.target]

        if src_comp == tgt_comp:
            mol = components[src_comp]
            a = _find_r_atom(mol, monomer_position=src_pos, r_label=src_r, node_id=edge.source)
            b = _find_r_atom(mol, monomer_position=tgt_pos, r_label=tgt_r, node_id=edge.target)
            components[src_comp] = cyclize_on_r_atoms(mol, a, b)
            continue

        left_mol = components[src_comp]
        right_mol = components[tgt_comp]

        left_idx = _find_r_atom(
            left_mol,
            monomer_position=src_pos,
            r_label=src_r,
            node_id=edge.source,
        )
        right_idx = _find_r_atom(
            right_mol,
            monomer_position=tgt_pos,
            r_label=tgt_r,
            node_id=edge.target,
        )

        merged = _merge_mols_on_r_atoms(left_mol, left_idx, right_mol, right_idx)

        # Keep src component id, absorb tgt component id
        components[src_comp] = merged
        del components[tgt_comp]
        for node_id, comp_id in list(node_component.items()):
            if comp_id == tgt_comp:
                node_component[node_id] = src_comp

    return list(components.values())


# ---------------------------------------------------------------------------
# Capping remaining R-groups
# ---------------------------------------------------------------------------


def _cap_remaining_rgroups(mol: Chem.Mol, monomer_store: MonomerStore) -> Chem.Mol:
    current = Chem.Mol(mol)

    while True:
        r_idx = _first_open_r_atom(current)
        if r_idx is None:
            break

        atom = current.GetAtomWithIdx(r_idx)
        monomer_id = atom.GetProp("helm_monomer_id") if atom.HasProp("helm_monomer_id") else ""
        polymer_type = atom.GetProp("helm_polymer_type") if atom.HasProp("helm_polymer_type") else ""
        r_label = atom.GetProp("helm_r_label") if atom.HasProp("helm_r_label") else ""

        monomer = monomer_store.get_monomer_db().get(polymer_type, {}).get(monomer_id)
        if monomer is None:
            # fallback: just remove atom and satisfy valence implicitly
            current = _strip_r_atom(current, r_idx)
            continue

        attachment = monomer.get_attachment(r_label)
        if attachment is None or not attachment.cap_group_smiles:
            current = _strip_r_atom(current, r_idx)
            continue

        cap = ChemEngine.get_molecule(attachment.cap_group_smiles)
        _annotate_cap(cap, r_label)

        cap_r_idx = _find_cap_r_atom(cap, r_label)
        current = _merge_mols_on_r_atoms(current, r_idx, cap, cap_r_idx)

    return current


def _first_open_r_atom(mol: Chem.Mol) -> Optional[int]:
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            return atom.GetIdx()
    return None


def _strip_r_atom(mol: Chem.Mol, r_idx: int) -> Chem.Mol:
    rw = Chem.RWMol(mol)
    rw.RemoveAtom(r_idx)
    out = rw.GetMol()
    Chem.SanitizeMol(out)
    return out


def _annotate_cap(cap_mol: Chem.Mol, r_label: str) -> None:
    for atom in cap_mol.GetAtoms():
        if atom.GetAtomicNum() != 0:
            continue
        label = ChemEngine._get_r_label(atom)
        if not label:
            isotope = atom.GetIsotope()
            amap = atom.GetAtomMapNum()
            if isotope > 0:
                label = f"R{isotope}"
            elif amap > 0:
                label = f"R{amap}"
        if label:
            atom.SetProp("helm_r_label", label)


def _find_cap_r_atom(cap_mol: Chem.Mol, label: str) -> int:
    # preferred exact label
    for atom in cap_mol.GetAtoms():
        if atom.GetAtomicNum() == 0 and atom.HasProp("helm_r_label"):
            if atom.GetProp("helm_r_label") == label:
                return atom.GetIdx()

    # fallback first R atom
    for atom in cap_mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            return atom.GetIdx()

    raise StructureException("Cap group does not contain an R-group")
