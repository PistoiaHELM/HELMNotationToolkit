"""Phase 3 port of core org.helm.notation.tools.SimpleNotationParser APIs.

Current scope:
- linear PEPTIDE simple notation (e.g. A.G.C)
- linear RNA simple notation (e.g. R(A)P.R(U)P)
- single-monomer CHEM notation (e.g. MCC)

This is intentionally incremental and geared to baseline parity checkpoints.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from rdkit import Chem

from helm_notation.exceptions import MonomerException, NotationException, StructureException
from helm_notation.model import Monomer, RgroupStructure
from helm_notation.monomer_store import MonomerStore
from helm_notation.monomer_factory import MonomerFactory
from helm_notation.chem import ChemEngine

_RNA_TOKEN_RE = re.compile(r"^\[?([A-Za-z0-9]+)\]?\(([A-Za-z0-9]+)\)(P?)$")


@dataclass
class _MonomerToken:
    monomer_id: str
    polymer_type: str
    position: int


class SimpleNotationParser:
    GROUP_LEVEL_DELIMITER = "."

    @staticmethod
    def get_monomer(
        monomer_id: str,
        polymer_type: str,
        monomer_store: Optional[MonomerStore] = None,
    ) -> Monomer:
        store = monomer_store or MonomerFactory.get_instance().get_monomer_store()
        poly_db = store.get_monomer_db().get(polymer_type, {})
        monomer = poly_db.get(monomer_id)
        if monomer is None:
            raise MonomerException(
                f"Monomer '{monomer_id}' not found in polymer type '{polymer_type}'"
            )
        return monomer

    @staticmethod
    def get_simple_polymer_structure(
        polymer_notation: str,
        polymer_type: str,
        monomer_store: Optional[MonomerStore] = None,
        node_id: Optional[str] = None,
    ) -> RgroupStructure:
        store = monomer_store or MonomerFactory.get_instance().get_monomer_store()

        polymer_type = polymer_type.upper()
        if polymer_type == "PEPTIDE":
            mol = _build_linear_peptide(polymer_notation, store, node_id=node_id)
        elif polymer_type == "RNA":
            mol = _build_linear_rna(polymer_notation, store, node_id=node_id)
        elif polymer_type == "CHEM":
            mol = _build_chem(polymer_notation, store, node_id=node_id)
        else:
            raise NotationException(f"Unsupported polymer type: {polymer_type}")

        rgroup_map: Dict[str, int] = {}
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 0:
                continue
            label = atom.GetProp("helm_r_label") if atom.HasProp("helm_r_label") else None
            if label:
                # keep the first one, matching Java behavior at map granularity
                rgroup_map.setdefault(label, atom.GetIdx())

        return RgroupStructure(molecule=mol, rgroup_map=rgroup_map)

    @staticmethod
    def get_simple_polymer_smiles(
        polymer_notation: str,
        polymer_type: str,
        monomer_store: Optional[MonomerStore] = None,
    ) -> str:
        structure = SimpleNotationParser.get_simple_polymer_structure(
            polymer_notation, polymer_type, monomer_store
        )
        if structure.molecule is None:
            raise NotationException("Unable to build polymer structure")
        clean = _strip_internal_atom_props(structure.molecule)
        return ChemEngine.get_unique_extended_smiles(clean)


def _build_chem(polymer_notation: str, store: MonomerStore, node_id: Optional[str]) -> Chem.Mol:
    # For Phase 3, support a single monomer ID in CHEM notation.
    monomer_id = polymer_notation.strip()
    monomer = SimpleNotationParser.get_monomer(monomer_id, "CHEM", store)
    mol = ChemEngine.get_molecule(monomer.can_smiles)
    _annotate_rgroups(mol, position=1, monomer=monomer, node_id=node_id)
    return mol


def _build_linear_peptide(polymer_notation: str, store: MonomerStore, node_id: Optional[str]) -> Chem.Mol:
    tokens = [token.strip() for token in polymer_notation.split(".") if token.strip()]
    if not tokens:
        raise NotationException("Empty peptide notation")

    acc: Optional[Chem.Mol] = None
    prev_position: Optional[int] = None
    position = 0

    for monomer_id in tokens:
        position += 1
        monomer = SimpleNotationParser.get_monomer(monomer_id, "PEPTIDE", store)
        mol = ChemEngine.get_molecule(monomer.can_smiles)
        _annotate_rgroups(mol, position=position, monomer=monomer, node_id=node_id)

        if acc is None:
            acc = mol
            prev_position = position
            continue

        left_idx = _find_r_atom(acc, monomer_position=prev_position, r_label="R2")
        right_idx = _find_r_atom(mol, monomer_position=position, r_label="R1")
        acc = _merge_mols_on_r_atoms(acc, left_idx, mol, right_idx)
        prev_position = position

    assert acc is not None
    return acc


def _build_linear_rna(polymer_notation: str, store: MonomerStore, node_id: Optional[str]) -> Chem.Mol:
    tokens = [token.strip() for token in polymer_notation.split(".") if token.strip()]
    if not tokens:
        raise NotationException("Empty RNA notation")

    acc: Optional[Chem.Mol] = None
    prev_phosphate_position: Optional[int] = None
    monomer_position = 0

    for nt in tokens:
        match = _RNA_TOKEN_RE.match(nt)
        if not match:
            raise NotationException(f"Unsupported RNA token format: {nt}")

        sugar_id, base_id, has_phosphate = match.groups()

        monomer_position += 1
        sugar_pos = monomer_position
        sugar = SimpleNotationParser.get_monomer(sugar_id, "RNA", store)
        sugar_mol = ChemEngine.get_molecule(sugar.can_smiles)
        _annotate_rgroups(sugar_mol, position=sugar_pos, monomer=sugar, node_id=node_id)

        if acc is None:
            acc = sugar_mol
        else:
            assert prev_phosphate_position is not None
            left_idx = _find_r_atom(acc, monomer_position=prev_phosphate_position, r_label="R2")
            right_idx = _find_r_atom(sugar_mol, monomer_position=sugar_pos, r_label="R1")
            acc = _merge_mols_on_r_atoms(acc, left_idx, sugar_mol, right_idx)

        monomer_position += 1
        base_pos = monomer_position
        base = SimpleNotationParser.get_monomer(base_id, "RNA", store)
        base_mol = ChemEngine.get_molecule(base.can_smiles)
        _annotate_rgroups(base_mol, position=base_pos, monomer=base, node_id=node_id)

        sugar_r3_idx = _find_r_atom(acc, monomer_position=sugar_pos, r_label="R3")
        base_r1_idx = _find_r_atom(base_mol, monomer_position=base_pos, r_label="R1")
        acc = _merge_mols_on_r_atoms(acc, sugar_r3_idx, base_mol, base_r1_idx)

        if has_phosphate:
            monomer_position += 1
            phosphate_pos = monomer_position
            phosphate = SimpleNotationParser.get_monomer("P", "RNA", store)
            phosphate_mol = ChemEngine.get_molecule(phosphate.can_smiles)
            _annotate_rgroups(
                phosphate_mol,
                position=phosphate_pos,
                monomer=phosphate,
                node_id=node_id,
            )

            sugar_r2_idx = _find_r_atom(acc, monomer_position=sugar_pos, r_label="R2")
            p_r1_idx = _find_r_atom(phosphate_mol, monomer_position=phosphate_pos, r_label="R1")
            acc = _merge_mols_on_r_atoms(acc, sugar_r2_idx, phosphate_mol, p_r1_idx)
            prev_phosphate_position = phosphate_pos
        else:
            prev_phosphate_position = None

    if acc is None:
        raise NotationException("Unable to construct RNA structure")
    return acc


def _annotate_rgroups(mol: Chem.Mol, position: int, monomer: Monomer, node_id: Optional[str]) -> None:
    for atom in mol.GetAtoms():
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
            else:
                label = "R?"

        atom.SetProp("helm_r_label", label)
        atom.SetProp("helm_monomer_pos", str(position))
        atom.SetProp("helm_monomer_id", monomer.alternate_id)
        atom.SetProp("helm_polymer_type", monomer.polymer_type)
        if node_id:
            atom.SetProp("helm_node_id", node_id)


def _find_r_atom(
    mol: Chem.Mol,
    monomer_position: Optional[int],
    r_label: str,
    node_id: Optional[str] = None,
) -> int:
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 0:
            continue
        if not atom.HasProp("helm_r_label"):
            continue
        if atom.GetProp("helm_r_label") != r_label:
            continue
        if monomer_position is not None:
            if not atom.HasProp("helm_monomer_pos"):
                continue
            if atom.GetProp("helm_monomer_pos") != str(monomer_position):
                continue
        if node_id is not None:
            if not atom.HasProp("helm_node_id"):
                continue
            if atom.GetProp("helm_node_id") != node_id:
                continue
        return atom.GetIdx()
    raise StructureException(
        f"Unable to locate R-group atom for pos={monomer_position} label={r_label} node={node_id}"
    )


def _merge_mols_on_r_atoms(
    left_mol: Chem.Mol,
    left_r_idx: int,
    right_mol: Chem.Mol,
    right_r_idx: int,
) -> Chem.Mol:
    combined = Chem.CombineMols(left_mol, right_mol)
    rw = Chem.RWMol(combined)

    right_r_idx_offset = left_mol.GetNumAtoms() + right_r_idx

    left_neighbors = [n.GetIdx() for n in rw.GetAtomWithIdx(left_r_idx).GetNeighbors()]
    right_neighbors = [n.GetIdx() for n in rw.GetAtomWithIdx(right_r_idx_offset).GetNeighbors()]

    if len(left_neighbors) != 1 or len(right_neighbors) != 1:
        raise StructureException("Only one bond is allowed from Rgroup")

    left_anchor = left_neighbors[0]
    right_anchor = right_neighbors[0]

    for atom_idx in sorted([left_r_idx, right_r_idx_offset], reverse=True):
        rw.RemoveAtom(atom_idx)
        if left_anchor > atom_idx:
            left_anchor -= 1
        if right_anchor > atom_idx:
            right_anchor -= 1

    rw.AddBond(left_anchor, right_anchor, Chem.BondType.SINGLE)
    merged = rw.GetMol()
    Chem.SanitizeMol(merged)
    return merged


def cyclize_on_r_atoms(mol: Chem.Mol, r_idx_a: int, r_idx_b: int) -> Chem.Mol:
    rw = Chem.RWMol(mol)

    neighbors_a = [n.GetIdx() for n in rw.GetAtomWithIdx(r_idx_a).GetNeighbors()]
    neighbors_b = [n.GetIdx() for n in rw.GetAtomWithIdx(r_idx_b).GetNeighbors()]
    if len(neighbors_a) != 1 or len(neighbors_b) != 1:
        raise StructureException("Only one bond is allowed from Rgroup")

    a_anchor = neighbors_a[0]
    b_anchor = neighbors_b[0]

    for atom_idx in sorted([r_idx_a, r_idx_b], reverse=True):
        rw.RemoveAtom(atom_idx)
        if a_anchor > atom_idx:
            a_anchor -= 1
        if b_anchor > atom_idx:
            b_anchor -= 1

    rw.AddBond(a_anchor, b_anchor, Chem.BondType.SINGLE)
    cycled = rw.GetMol()
    Chem.SanitizeMol(cycled)
    return cycled


def _strip_internal_atom_props(mol: Chem.Mol) -> Chem.Mol:
    clone = Chem.Mol(mol)
    keys = (
        "helm_r_label",
        "helm_monomer_pos",
        "helm_monomer_id",
        "helm_polymer_type",
        "helm_node_id",
    )
    for atom in clone.GetAtoms():
        for key in keys:
            if atom.HasProp(key):
                atom.ClearProp(key)
    return clone
