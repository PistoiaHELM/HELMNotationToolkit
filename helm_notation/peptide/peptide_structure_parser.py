"""Convert peptide chemical structures into HELM polymer notation.

Handles amide bonds (C(=O)-N) and disulfide bonds (S-S).

Steps:
1. Find all amide and disulfide bonds in the structure.
2. Exclude preserved amide bonds (urea, urethane, imide).
3. Fragment the molecule by breaking these bonds, tracking history in a tree.
4. Close opened lactam rings in any fragment.
5. Match each fragment against the monomer core structure database.
6. Generate extended SMILES for unmatched fragments (inline HELM notation).
7. Combine broken side-chain amine from Asn/Gln with Asp/Glu.
8. Generate HELM notation from monomer IDs and tree model.

Limitations:
- Only supports R1, R2, and R3 attachments.
- R1/R2/R3 assignment may be imprecise (based on alpha-amino-acid configuration).
"""

from __future__ import annotations

import collections
import logging
import random
from typing import Optional

from rdkit import Chem
from rdkit.Chem import RWMol

from helm_notation.exceptions import MonomerException, StructureException
from helm_notation.model.monomer import PEPTIDE_POLYMER_TYPE, Monomer
from helm_notation.peptide.amino_acid_info import AminoAcidInfo
from helm_notation.peptide.connection_info import ConnectionInfo
from helm_notation.peptide.peptide_fragment import PeptideFragment
from helm_notation.peptide.tree import Tree
from helm_notation.tools import structure_parser
from helm_notation.tools.simple_notation_parser import (
    GROUP_LEVEL_DELIMITER,
    MODIFICATION_END_SYMBOL,
    MODIFICATION_START_SYMBOL,
)

logger = logging.getLogger(__name__)

# Attachment-type labels stored on dummy atoms via ``SetProp``
AMINE_ATTACHMENT_LABEL = "R1-R3"
CARBONYL_ATTACHMENT_LABEL = "R2-R3"
SULFIDE_ATTACHMENT_LABEL = "R3"

_ATTACHMENT_PROP = "_peptide_attachment_label"
_RGROUP_PROP = "_peptide_rgroup_id"


# ---------------------------------------------------------------------------
# Module-level singleton
# ---------------------------------------------------------------------------

class PeptideStructureParser:
    """Singleton that converts peptide structures into HELM notation."""

    _instance: PeptideStructureParser | None = None

    def __init__(self) -> None:
        self._id_list: list[str] = []
        self._non_capped: list[str] = []
        self._r1_capped: list[str] = []
        self._r2_capped: list[str] = []
        self._r3_capped: list[str] = []
        self._r1r2_capped: list[str] = []
        self._r1r3_capped: list[str] = []
        self._r2r3_capped: list[str] = []
        self._all_capped: list[str] = []
        self.seed_id: int = 1

    @classmethod
    def get_instance(cls) -> PeptideStructureParser:
        if cls._instance is None:
            cls._instance = PeptideStructureParser()
        return cls._instance

    # ------------------------------------------------------------------
    # Initialisation
    # ------------------------------------------------------------------

    def init_amino_acid_lists(self) -> None:
        """Pre-compute capped canonical SMILES for all known peptide monomers."""
        from helm_notation.monomer_factory import get_instance as get_factory

        id_monomer_map = get_factory().get_monomer_db().get(PEPTIDE_POLYMER_TYPE, {})

        for mid, monomer in id_monomer_map.items():
            self._id_list.append(mid)
            self._non_capped.append(_get_capped_unique_smiles(monomer, []))
            self._r1_capped.append(_get_capped_unique_smiles(monomer, [1]))
            self._r2_capped.append(_get_capped_unique_smiles(monomer, [2]))
            self._r3_capped.append(_get_capped_unique_smiles(monomer, [3]))
            self._r1r2_capped.append(_get_capped_unique_smiles(monomer, [1, 2]))
            self._r1r3_capped.append(_get_capped_unique_smiles(monomer, [1, 3]))
            self._r2r3_capped.append(_get_capped_unique_smiles(monomer, [2, 3]))
            self._all_capped.append(_get_capped_unique_smiles(monomer, [1, 2, 3]))

    # ------------------------------------------------------------------
    # Public entry points
    # ------------------------------------------------------------------

    def molfile2notation(self, molfile: str) -> str:
        """Convert a molfile string to HELM notation."""
        mol = structure_parser.get_molecule(molfile)
        if mol is None:
            raise StructureException("Cannot parse molfile input")
        return self.molecule2notation(mol)

    def smiles2notation(self, smiles: str) -> str:
        """Convert a SMILES string to HELM notation."""
        mol = structure_parser.get_molecule(smiles)
        if mol is None:
            raise StructureException(f"Cannot parse SMILES: {smiles}")
        return self.molecule2notation(mol)

    def molecule2notation(self, molecule: Chem.Mol) -> str:
        """Convert an RDKit Mol to HELM notation."""
        rwmol = RWMol(molecule)
        frags = Chem.GetMolFrags(rwmol, asMols=True)
        if len(frags) > 1:
            raise StructureException("Input structure can't be a mixture")

        self.seed_id = 1

        pf = PeptideFragment(molecule=rwmol, level=1)
        frag_tree: Tree[PeptideFragment] = Tree(pf)

        self.break_peptide_bonds(frag_tree)
        _close_opened_lactam_tree(frag_tree)
        _roll_up(self, frag_tree)

        monomer_ids_list = frag_tree.head.monomer_ids_list
        connection_maps_list = frag_tree.head.connection_maps_list

        if monomer_ids_list is None or connection_maps_list is None:
            raise StructureException("Failed to parse peptide structure")

        _combine_opened_side_chain_amide(monomer_ids_list, connection_maps_list)
        notation = _generate_notation(monomer_ids_list, connection_maps_list)
        return notation

    # ------------------------------------------------------------------
    # Bond-breaking methods
    # ------------------------------------------------------------------

    def break_peptide_bonds(self, tree: Tree[PeptideFragment]) -> None:
        """Recursively break all amide and disulfide bonds."""
        frag = tree.head
        if frag.id <= 0:
            frag.id = self.seed_id
            self.seed_id += 1
        mol = frag.molecule
        if mol is None:
            return

        current_level = frag.level
        current_id = frag.id

        disulfide = _find_disulfide_bond(mol)
        amide = _find_amide_bond(mol)

        if disulfide is not None or amide is not None:
            if disulfide is not None:
                mols = _break_disulfide_bond(mol, disulfide, current_id)
            else:
                mols = _break_amide_bond(mol, amide, current_id)  # type: ignore[arg-type]

            next_level = current_level + 1

            if len(mols) == 1:
                child_pf = PeptideFragment(molecule=mols[0], level=next_level)
                child_tree = tree.add_leaf(child_pf)
                self.break_peptide_bonds(child_tree)
            else:
                pf1 = PeptideFragment(molecule=mols[0], level=next_level)
                t1 = tree.add_leaf(pf1)
                self.break_peptide_bonds(t1)

                pf2 = PeptideFragment(molecule=mols[1], level=next_level)
                t2 = tree.add_leaf(pf2)
                self.break_peptide_bonds(t2)

    def break_disulfide_bonds(self, tree: Tree[PeptideFragment]) -> None:
        """Recursively break only disulfide bonds."""
        frag = tree.head
        if frag.id <= 0:
            frag.id = self.seed_id
            self.seed_id += 1
        mol = frag.molecule
        if mol is None:
            return

        current_level = frag.level
        current_id = frag.id
        bond_idx = _find_disulfide_bond(mol)
        if bond_idx is not None:
            mols = _break_disulfide_bond(mol, bond_idx, current_id)
            next_level = current_level + 1

            if len(mols) == 1:
                child_pf = PeptideFragment(molecule=mols[0], level=next_level)
                child_tree = tree.add_leaf(child_pf)
                self.break_disulfide_bonds(child_tree)
            else:
                pf1 = PeptideFragment(molecule=mols[0], level=next_level)
                t1 = tree.add_leaf(pf1)
                self.break_disulfide_bonds(t1)
                pf2 = PeptideFragment(molecule=mols[1], level=next_level)
                t2 = tree.add_leaf(pf2)
                self.break_disulfide_bonds(t2)

    def break_amide_bonds(self, tree: Tree[PeptideFragment]) -> None:
        """Recursively break only amide bonds."""
        frag = tree.head
        if frag.id <= 0:
            frag.id = self.seed_id
            self.seed_id += 1
        mol = frag.molecule
        if mol is None:
            return

        current_level = frag.level
        current_id = frag.id
        bond_idx = _find_amide_bond(mol)
        if bond_idx is not None:
            mols = _break_amide_bond(mol, bond_idx, current_id)
            next_level = current_level + 1

            if len(mols) == 1:
                child_pf = PeptideFragment(molecule=mols[0], level=next_level)
                child_tree = tree.add_leaf(child_pf)
                self.break_amide_bonds(child_tree)
            else:
                pf1 = PeptideFragment(molecule=mols[0], level=next_level)
                t1 = tree.add_leaf(pf1)
                self.break_amide_bonds(t1)
                pf2 = PeptideFragment(molecule=mols[1], level=next_level)
                t2 = tree.add_leaf(pf2)
                self.break_amide_bonds(t2)


# ---------------------------------------------------------------------------
# Module-level helper functions (private)
# ---------------------------------------------------------------------------


def _get_capped_unique_smiles(monomer: Monomer, rgroup_ids: list[int]) -> str:
    """Compute the canonical SMILES for a monomer with certain R-groups capped."""
    mol_smi = monomer.can_smiles
    mol = structure_parser.get_molecule(mol_smi)
    if mol is None:
        return ""

    rwmol = RWMol(mol)
    for rid in rgroup_ids:
        try:
            r_idx = structure_parser.find_rgroup_atom(rwmol, rid)
        except StructureException:
            continue
        att = monomer.get_attachment(f"R{rid}")
        if att is None:
            continue
        att_smi = att.cap_group_smiles
        att_mol = structure_parser.get_molecule(att_smi)
        if att_mol is None:
            continue
        att_rwmol = RWMol(att_mol)
        try:
            att_r_idx = structure_parser.find_rgroup_atom(att_rwmol, rid)
        except StructureException:
            continue
        rwmol = structure_parser.merge(rwmol, r_idx, att_rwmol, att_r_idx)

    return structure_parser.get_unique_smiles_from_mol(rwmol)


# ---------------------------------------------------------------------------
# Atom / bond classification helpers (RDKit equivalents of ChemAxon checks)
# ---------------------------------------------------------------------------


def _is_carbonyl_carbon_atom(mol: Chem.Mol, atom_idx: int) -> bool:
    """Check if atom at *atom_idx* is a carbonyl carbon: C(=O) with 3 total connections."""
    atom = mol.GetAtomWithIdx(atom_idx)
    if atom.GetAtomicNum() != 6:
        return False
    total_degree = atom.GetDegree() + atom.GetTotalNumHs()
    if total_degree != 3:
        return False
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        if bond.GetBondTypeAsDouble() == 2.0 and other.GetAtomicNum() == 8:
            return True
    return False


def _is_amino_nitrogen_atom(mol: Chem.Mol, atom_idx: int) -> bool:
    """Check if atom is an amino nitrogen (N with total degree 3)."""
    atom = mol.GetAtomWithIdx(atom_idx)
    if atom.GetAtomicNum() != 7:
        return False
    total_degree = atom.GetDegree() + atom.GetTotalNumHs()
    return total_degree == 3


def _is_disulfide_sulfur_atom(mol: Chem.Mol, atom_idx: int) -> bool:
    """Check if atom is a sulfur participating in a disulfide (R-S-S-R, degree 2)."""
    atom = mol.GetAtomWithIdx(atom_idx)
    return atom.GetAtomicNum() == 16 and atom.GetDegree() == 2


def _is_disulfide_bond(mol: Chem.Mol, bond_idx: int) -> bool:
    """Check if a bond is a disulfide bond (S-S, single)."""
    bond = mol.GetBondWithIdx(bond_idx)
    if bond.GetBondTypeAsDouble() != 1.0:
        return False
    a1 = bond.GetBeginAtomIdx()
    a2 = bond.GetEndAtomIdx()
    return _is_disulfide_sulfur_atom(mol, a1) and _is_disulfide_sulfur_atom(mol, a2)


def _is_amide_bond(mol: Chem.Mol, bond_idx: int) -> bool:
    """Check if a bond is a cleavable amide bond (C(=O)-N, not preserved)."""
    bond = mol.GetBondWithIdx(bond_idx)
    if bond.GetBondTypeAsDouble() != 1.0:
        return False
    a1 = bond.GetBeginAtomIdx()
    a2 = bond.GetEndAtomIdx()

    found = False
    if _is_carbonyl_carbon_atom(mol, a1) and _is_amino_nitrogen_atom(mol, a2):
        found = True
    elif _is_carbonyl_carbon_atom(mol, a2) and _is_amino_nitrogen_atom(mol, a1):
        found = True

    if found and not _is_preserved_amide_bond(mol, bond_idx):
        return True
    return False


def _is_preserved_amide_bond(mol: Chem.Mol, bond_idx: int) -> bool:
    """Check if an amide bond should be preserved (urea, urethane, or imide)."""
    bond = mol.GetBondWithIdx(bond_idx)
    a1_idx = bond.GetBeginAtomIdx()
    a2_idx = bond.GetEndAtomIdx()

    a1 = mol.GetAtomWithIdx(a1_idx)
    a2 = mol.GetAtomWithIdx(a2_idx)

    if a1.GetAtomicNum() == 6:
        c_idx, n_idx = a1_idx, a2_idx
    else:
        c_idx, n_idx = a2_idx, a1_idx

    c_atom = mol.GetAtomWithIdx(c_idx)
    n_atom = mol.GetAtomWithIdx(n_idx)

    # Urea: N-C(=O)-N
    match = 0
    for b in c_atom.GetBonds():
        other = b.GetOtherAtom(c_atom)
        if b.GetBondTypeAsDouble() == 2.0 and other.GetAtomicNum() == 8:
            match += 1
        elif b.GetBondTypeAsDouble() == 1.0 and other.GetAtomicNum() == 7:
            match += 1
    if match == 3:
        return True

    # Urethane: O-C(=O)-N
    match = 0
    for b in c_atom.GetBonds():
        other = b.GetOtherAtom(c_atom)
        if b.GetBondTypeAsDouble() == 2.0 and other.GetAtomicNum() == 8:
            match += 1
        elif b.GetBondTypeAsDouble() == 1.0 and other.GetAtomicNum() == 7:
            match += 1
        elif b.GetBondTypeAsDouble() == 1.0 and other.GetAtomicNum() == 8:
            match += 1
    if match == 3:
        return True

    # Imide: C(=O)-N-C(=O)
    match = 0
    for b in n_atom.GetBonds():
        other = b.GetOtherAtom(n_atom)
        if b.GetBondTypeAsDouble() == 1.0 and other.GetAtomicNum() == 6 and _is_carbonyl_carbon_atom(mol, other.GetIdx()):
            match += 1
    if match == 2:
        return True

    return False


def _find_disulfide_bond(mol: Chem.Mol) -> int | None:
    """Return the index of the first disulfide bond, or ``None``."""
    for bond in mol.GetBonds():
        if _is_disulfide_bond(mol, bond.GetIdx()):
            return bond.GetIdx()
    return None


def _find_amide_bond(mol: Chem.Mol) -> int | None:
    """Return the index of the first cleavable amide bond, or ``None``."""
    for bond in mol.GetBonds():
        if _is_amide_bond(mol, bond.GetIdx()):
            return bond.GetIdx()
    return None


# ---------------------------------------------------------------------------
# Bond-breaking helpers
# ---------------------------------------------------------------------------


def _add_rgroup_dummy(rwmol: RWMol, neighbor_idx: int, rgroup_id: int, label: str) -> int:
    """Add a dummy atom (R-group) bonded to *neighbor_idx* and return its index."""
    dummy_idx = rwmol.AddAtom(Chem.Atom(0))
    dummy_atom = rwmol.GetAtomWithIdx(dummy_idx)
    dummy_atom.SetIsotope(rgroup_id)
    dummy_atom.SetProp(_ATTACHMENT_PROP, label)
    dummy_atom.SetIntProp(_RGROUP_PROP, rgroup_id)
    rwmol.AddBond(neighbor_idx, dummy_idx, Chem.BondType.SINGLE)
    return dummy_idx


def _get_attachment_label(mol: Chem.Mol, atom_idx: int) -> str:
    """Retrieve the attachment label from a dummy atom."""
    atom = mol.GetAtomWithIdx(atom_idx)
    if atom.HasProp(_ATTACHMENT_PROP):
        return atom.GetProp(_ATTACHMENT_PROP)
    return ""


def _get_rgroup_id(mol: Chem.Mol, atom_idx: int) -> int:
    """Retrieve the R-group id stored on a dummy atom."""
    atom = mol.GetAtomWithIdx(atom_idx)
    if atom.HasProp(_RGROUP_PROP):
        return atom.GetIntProp(_RGROUP_PROP)
    return atom.GetIsotope()


def _break_disulfide_bond(mol: Chem.Mol, bond_idx: int, rgroup_id: int) -> list[RWMol]:
    """Break a disulfide bond, replacing it with two R-group dummy atoms."""
    rwmol = RWMol(mol)
    bond = rwmol.GetBondWithIdx(bond_idx)
    a1_idx = bond.GetBeginAtomIdx()
    a2_idx = bond.GetEndAtomIdx()

    rwmol.RemoveBond(a1_idx, a2_idx)
    _add_rgroup_dummy(rwmol, a1_idx, rgroup_id, SULFIDE_ATTACHMENT_LABEL)
    _add_rgroup_dummy(rwmol, a2_idx, rgroup_id, SULFIDE_ATTACHMENT_LABEL)

    frags = Chem.GetMolFrags(rwmol.GetMol(), asMols=True)
    return [RWMol(f) for f in frags]


def _break_amide_bond(mol: Chem.Mol, bond_idx: int, rgroup_id: int) -> list[RWMol]:
    """Break an amide bond, labelling the carbonyl side as R2-type and amino as R1-type."""
    rwmol = RWMol(mol)
    bond = rwmol.GetBondWithIdx(bond_idx)
    a1_idx = bond.GetBeginAtomIdx()
    a2_idx = bond.GetEndAtomIdx()

    if _is_carbonyl_carbon_atom(rwmol, a1_idx):
        carbonyl_idx, amino_idx = a1_idx, a2_idx
    else:
        carbonyl_idx, amino_idx = a2_idx, a1_idx

    rwmol.RemoveBond(a1_idx, a2_idx)
    carbonyl_r_idx = _add_rgroup_dummy(rwmol, carbonyl_idx, rgroup_id, CARBONYL_ATTACHMENT_LABEL)
    amino_r_idx = _add_rgroup_dummy(rwmol, amino_idx, rgroup_id, AMINE_ATTACHMENT_LABEL)

    frags = Chem.GetMolFrags(rwmol.GetMol(), asMols=True)
    results = [RWMol(f) for f in frags]

    if len(results) == 1:
        return results

    # Ensure the carbonyl-containing fragment comes first
    # Check which fragment contains the carbonyl R-group dummy
    for i, frag in enumerate(results):
        for atom in frag.GetAtoms():
            if atom.GetAtomicNum() == 0 and atom.HasProp(_ATTACHMENT_PROP):
                if atom.GetProp(_ATTACHMENT_PROP) == CARBONYL_ATTACHMENT_LABEL:
                    if i == 0:
                        return results
                    else:
                        return [results[1], results[0]]

    return results


# ---------------------------------------------------------------------------
# Lactam closure
# ---------------------------------------------------------------------------


def _close_opened_lactam_mol(mol: RWMol) -> RWMol:
    """Close lactam rings in a single molecule fragment.

    If two dummy atoms share the same R-group id (isotope), replace them
    with a direct bond between their neighbours.
    """
    rwmol = RWMol(mol)
    level_map: dict[int, list[int]] = collections.defaultdict(list)

    for atom in rwmol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            rg = atom.GetIsotope()
            if rg > 0:
                level_map[rg].append(atom.GetIdx())

    # Collect pairs that need closing (sorted descending so removal is safe)
    pairs_to_close: list[tuple[int, int]] = []
    for rg, idxs in level_map.items():
        if len(idxs) == 2:
            pairs_to_close.append((max(idxs), min(idxs)))

    if not pairs_to_close:
        return rwmol

    # Process each pair: remove dummy atoms and add bond between neighbours
    # Sort by descending first index so we remove from the end
    pairs_to_close.sort(reverse=True)
    for idx_high, idx_low in pairs_to_close:
        atom_high = rwmol.GetAtomWithIdx(idx_high)
        atom_low = rwmol.GetAtomWithIdx(idx_low)

        # Get the neighbor of each dummy
        neighbor_high_idx = atom_high.GetNeighbors()[0].GetIdx()
        neighbor_low_idx = atom_low.GetNeighbors()[0].GetIdx()

        # Remove higher index first
        rwmol.RemoveAtom(idx_high)
        # Adjust lower index if needed
        adjusted_low = idx_low if idx_low < idx_high else idx_low - 1
        # Adjust neighbor indices
        adj_neighbor_high = neighbor_high_idx
        if neighbor_high_idx > idx_high:
            adj_neighbor_high -= 1

        rwmol.RemoveAtom(adjusted_low)
        # Adjust both neighbor indices for the second removal
        if adj_neighbor_high > adjusted_low:
            adj_neighbor_high -= 1
        adj_neighbor_low = neighbor_low_idx
        if adj_neighbor_low > idx_high:
            adj_neighbor_low -= 1
        if adj_neighbor_low > adjusted_low:
            adj_neighbor_low -= 1

        rwmol.AddBond(adj_neighbor_high, adj_neighbor_low, Chem.BondType.SINGLE)

    return rwmol


def _close_opened_lactam_tree(tree: Tree[PeptideFragment]) -> None:
    """Recursively close lactam rings in all leaf fragments."""
    subtrees = tree.get_sub_trees()
    if not subtrees:
        mol = tree.head.molecule
        if mol is not None:
            tree.head.molecule = _close_opened_lactam_mol(RWMol(mol))
    else:
        for subtree in subtrees:
            _close_opened_lactam_tree(subtree)


# ---------------------------------------------------------------------------
# Alpha amino acid detection helpers
# ---------------------------------------------------------------------------


def _is_alpha_amino_acid_nitrogen_atom(mol: Chem.Mol, atom_idx: int) -> bool:
    """Check if the nitrogen is an alpha-amino-acid N (N-C-C(=O) pattern)."""
    atom = mol.GetAtomWithIdx(atom_idx)
    if atom.GetAtomicNum() != 7:
        return False
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        if bond.GetBondTypeAsDouble() == 1.0 and other.GetAtomicNum() == 6:
            # Check if this carbon is alpha-carbon (bonded to a carbonyl carbon)
            for cbond in other.GetBonds():
                neighbor = cbond.GetOtherAtom(other)
                if (
                    cbond.GetBondTypeAsDouble() == 1.0
                    and neighbor.GetAtomicNum() == 6
                    and _is_carbonyl_carbon_atom(mol, neighbor.GetIdx())
                ):
                    return True
    return False


def _is_alpha_amino_acid_carbonyl_carbon_atom(mol: Chem.Mol, atom_idx: int) -> bool:
    """Check if the carbon is the carbonyl C of an alpha amino acid."""
    atom = mol.GetAtomWithIdx(atom_idx)
    if atom.GetAtomicNum() != 6:
        return False
    found_oxygen = False
    found_alpha_nitrogen = False
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        if bond.GetBondTypeAsDouble() == 2.0 and other.GetAtomicNum() == 8:
            found_oxygen = True
        elif bond.GetBondTypeAsDouble() == 1.0 and other.GetAtomicNum() == 6:
            for cbond in other.GetBonds():
                neighbor = cbond.GetOtherAtom(other)
                if cbond.GetBondTypeAsDouble() == 1.0 and neighbor.GetAtomicNum() == 7:
                    found_alpha_nitrogen = True
                    break
    return found_oxygen and found_alpha_nitrogen


# ---------------------------------------------------------------------------
# Fragment -> AminoAcidInfo identification
# ---------------------------------------------------------------------------


def _molecule2amino_acid_info(
    parser: PeptideStructureParser, mol: Chem.Mol
) -> AminoAcidInfo:
    """Identify a molecular fragment as a known amino acid monomer."""
    can_smi = structure_parser.get_unique_smiles_from_mol(mol)
    monomer_id: str | None = None
    rgroup_map: dict[str, str] = {}

    # Try matching against all capped SMILES lists
    for smiles_list in [
        parser._non_capped,
        parser._r1_capped,
        parser._r2_capped,
        parser._r3_capped,
        parser._r1r2_capped,
        parser._r1r3_capped,
        parser._r2r3_capped,
        parser._all_capped,
    ]:
        if can_smi in smiles_list:
            idx = smiles_list.index(can_smi)
            monomer_id = parser._id_list[idx]
            break

    # Map R-group dummy atoms to R1/R2/R3 labels
    contains_sulfide = False
    r_atoms: list[int] = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0 and atom.GetIsotope() > 0:
            r_atoms.append(atom.GetIdx())
            label = _get_attachment_label(mol, atom.GetIdx())
            if label == SULFIDE_ATTACHMENT_LABEL:
                contains_sulfide = True

    for r_idx in r_atoms:
        atom = mol.GetAtomWithIdx(r_idx)
        rg_id = _get_rgroup_id(mol, r_idx)
        label = _get_attachment_label(mol, r_idx)
        key = f"R{rg_id}"

        if label.startswith(SULFIDE_ATTACHMENT_LABEL):
            rgroup_map[key] = "R3"
        elif label.startswith(AMINE_ATTACHMENT_LABEL):
            # Could be R1 or R3
            neighbor = atom.GetNeighbors()[0]
            n_idx = neighbor.GetIdx()
            if (
                _is_alpha_amino_acid_nitrogen_atom(mol, n_idx)
                or len(r_atoms) == 1
                or contains_sulfide
            ):
                rgroup_map[key] = "R1"
            else:
                rgroup_map[key] = "R3"
        elif label.startswith(CARBONYL_ATTACHMENT_LABEL):
            # Could be R2 or R3
            neighbor = atom.GetNeighbors()[0]
            n_idx = neighbor.GetIdx()
            if (
                _is_alpha_amino_acid_carbonyl_carbon_atom(mol, n_idx)
                or len(r_atoms) == 1
                or contains_sulfide
            ):
                rgroup_map[key] = "R2"
            else:
                rgroup_map[key] = "R3"

    if monomer_id is None:
        # Generate inline extended SMILES
        monomer_id = structure_parser.get_unique_extended_smiles_from_mol(mol)
        r_numbers = [int(k[1:]) for k in rgroup_map]
        random.shuffle(r_numbers)

        for rnum in reversed(r_numbers):
            old_r = f"R{rnum}"
            new_r = rgroup_map[old_r]
            new_x = new_r.replace("R", "X")
            monomer_id = monomer_id.replace(old_r, new_x)
        monomer_id = monomer_id.replace("X", "R")
    else:
        if monomer_id == "Ggu":
            # gamma glutamic acid => same structure as E with flipped R2/R3
            monomer_id = "E"

    return AminoAcidInfo(monomer_id=monomer_id, rgroup_map=rgroup_map)


# ---------------------------------------------------------------------------
# Roll-up: build monomer ID lists from the fragmentation tree
# ---------------------------------------------------------------------------


def _roll_up(parser: PeptideStructureParser, tree: Tree[PeptideFragment]) -> None:
    """Repeatedly perform single roll-up until the tree is fully resolved."""
    depth = tree.get_depth()
    for _ in range(depth):
        _single_roll_up(parser, tree)


def _single_roll_up(parser: PeptideStructureParser, tree: Tree[PeptideFragment]) -> None:
    """Attempt to resolve one level of the fragmentation tree."""
    frag = tree.head
    mol = frag.molecule
    monomer_ids_list = frag.monomer_ids_list

    if monomer_ids_list is None or not monomer_ids_list:
        if mol is not None and mol.GetNumAtoms() > 0:
            # Leaf node: identify the monomer
            mi = _molecule2amino_acid_info(parser, mol)
            frag.monomer_ids_list = [[mi.monomer_id]]
            frag.connection_maps_list = [[dict(mi.rgroup_map)]]
        else:
            subtrees = tree.get_sub_trees()
            if len(subtrees) == 1:
                subtree1 = subtrees[0]
                frag1 = subtree1.head
                if frag1.monomer_ids_list and frag1.monomer_ids_list:
                    frag.monomer_ids_list = frag1.monomer_ids_list
                    frag.connection_maps_list = frag1.connection_maps_list
                else:
                    _single_roll_up(parser, subtree1)
            elif len(subtrees) == 2:
                parent_id = frag.id
                subtree1 = subtrees[0]
                frag1 = subtree1.head
                ids_list1 = frag1.monomer_ids_list
                conn_list1 = frag1.connection_maps_list

                subtree2 = subtrees[1]
                frag2 = subtree2.head
                ids_list2 = frag2.monomer_ids_list
                conn_list2 = frag2.connection_maps_list

                if ids_list1 and ids_list1:
                    if ids_list2 and ids_list2:
                        assert conn_list1 is not None
                        assert conn_list2 is not None
                        ci1 = _get_connection_info(parent_id, ids_list1, conn_list1)
                        ci2 = _get_connection_info(parent_id, ids_list2, conn_list2)

                        parent_ids_list: list[list[str]] = []
                        parent_conn_list: list[list[dict[str, str]]] = []

                        if (
                            (ci1.rgroup == "R2" and ci2.rgroup == "R1")
                            or (ci2.rgroup == "R2" and ci1.rgroup == "R1")
                        ):
                            keep_order = not (ci2.rgroup == "R2" and ci1.rgroup == "R1")

                            # Extract the matching lists
                            monomer_ids_1: list[str] | None = None
                            conn_maps_1: list[dict[str, str]] | None = None
                            for i in range(len(ids_list1)):
                                if i == ci1.list_index:
                                    monomer_ids_1 = ids_list1[i]
                                    conn_maps_1 = conn_list1[i]
                                else:
                                    parent_ids_list.append(ids_list1[i])
                                    parent_conn_list.append(conn_list1[i])

                            assert monomer_ids_1 is not None and conn_maps_1 is not None
                            conn_maps_1[ci1.monomer_index].pop(f"R{parent_id}", None)

                            monomer_ids_2: list[str] | None = None
                            conn_maps_2: list[dict[str, str]] | None = None
                            for i in range(len(ids_list2)):
                                if i == ci2.list_index:
                                    monomer_ids_2 = ids_list2[i]
                                    conn_maps_2 = conn_list2[i]
                                else:
                                    parent_ids_list.append(ids_list2[i])
                                    parent_conn_list.append(conn_list2[i])

                            assert monomer_ids_2 is not None and conn_maps_2 is not None
                            conn_maps_2[ci2.monomer_index].pop(f"R{parent_id}", None)

                            if keep_order:
                                new_ids = monomer_ids_1 + monomer_ids_2
                                new_conns = conn_maps_1 + conn_maps_2
                            else:
                                new_ids = monomer_ids_2 + monomer_ids_1
                                new_conns = conn_maps_2 + conn_maps_1

                            parent_ids_list.append(new_ids)
                            parent_conn_list.append(new_conns)
                        else:
                            parent_ids_list.extend(ids_list1)
                            parent_ids_list.extend(ids_list2)
                            parent_conn_list.extend(conn_list1)
                            parent_conn_list.extend(conn_list2)

                        frag.monomer_ids_list = parent_ids_list
                        frag.connection_maps_list = parent_conn_list
                    else:
                        _single_roll_up(parser, subtree2)
                else:
                    _single_roll_up(parser, subtree1)
                    if ids_list2 is None or not ids_list2:
                        _single_roll_up(parser, subtree2)


def _get_connection_info(
    parent_id: int,
    monomer_ids_list: list[list[str]],
    connection_map_list: list[list[dict[str, str]]],
) -> ConnectionInfo:
    """Find the connection entry that references the given parent fragment ID."""
    connection_key = f"R{parent_id}"
    list_index = -1
    monomer_index = -1
    rgroup = ""
    monomer_id = ""

    for i, conn_maps in enumerate(connection_map_list):
        for j, cmap in enumerate(conn_maps):
            if connection_key in cmap:
                list_index = i
                monomer_index = j
                rgroup = cmap[connection_key]
                break

    if list_index < 0 or monomer_index < 0 or not rgroup:
        raise StructureException("Unable to find connection info")

    monomer_id = monomer_ids_list[list_index][monomer_index]
    return ConnectionInfo(
        list_index=list_index,
        monomer_index=monomer_index,
        monomer_id=monomer_id,
        rgroup=rgroup,
    )


# ---------------------------------------------------------------------------
# Side-chain amide combination (Asp/Glu + am -> Asn/Gln)
# ---------------------------------------------------------------------------


def _combine_opened_side_chain_amide(
    monomer_ids_list: list[list[str]],
    connection_maps_list: list[list[dict[str, str]]],
) -> None:
    """Merge free amide fragments back into Asn (N) or Gln (Q)."""
    amine_indices: list[int] = []
    other_indices: list[int] = []

    for idx, ids in enumerate(monomer_ids_list):
        if len(ids) == 1 and ids[0] == "am":
            amine_indices.append(idx)
        else:
            other_indices.append(idx)

    removal_indices: list[int] = []

    for am_idx in amine_indices:
        conn_maps = connection_maps_list[am_idx]
        am_connections = conn_maps[0]
        level = next(iter(am_connections.keys()))

        found_match = False
        for other_idx in other_indices:
            ids = monomer_ids_list[other_idx]
            for mon_idx, mon_id in enumerate(ids):
                if mon_id in ("D", "E"):
                    other_conns = connection_maps_list[other_idx][mon_idx]
                    connection_r = other_conns.get(level)
                    if connection_r == "R3":
                        removal_indices.append(am_idx)
                        new_id = "N" if mon_id == "D" else "Q"
                        ids[mon_idx] = new_id
                        other_conns.pop(level, None)
                        found_match = True
                        break
            if found_match:
                break

    # Remove last-first so indices stay valid
    for i in sorted(removal_indices, reverse=True):
        monomer_ids_list.pop(i)
        connection_maps_list.pop(i)


# ---------------------------------------------------------------------------
# Notation generation
# ---------------------------------------------------------------------------


def _get_simple_notation(ids: list[str]) -> str:
    """Build a simple notation string from a list of monomer IDs."""
    parts: list[str] = []
    for mid in ids:
        if len(mid) > 1:
            parts.append(f"{MODIFICATION_START_SYMBOL}{mid}{MODIFICATION_END_SYMBOL}")
        else:
            parts.append(mid)
    return GROUP_LEVEL_DELIMITER.join(parts)


def _generate_notation(
    monomer_ids_list: list[list[str]],
    connection_maps_list: list[list[dict[str, str]]],
) -> str:
    """Generate a complete HELM notation string from parsed data."""
    simple_notations: list[str] = []
    for ids in monomer_ids_list:
        simple_notations.append(_get_simple_notation(ids))

    # Build connection map: key = level id, value = list of "polymerIndex:monomerIndex:rgroup"
    connection_map: dict[str, list[str]] = {}
    for i, conn_maps in enumerate(connection_maps_list):
        for j, cmap in enumerate(conn_maps):
            for level_id, rgroup in cmap.items():
                con_string = f"{i + 1}:{j + 1}:{rgroup}"
                if level_id not in connection_map:
                    connection_map[level_id] = []
                connection_map[level_id].append(con_string)

    # Build connection strings
    connections: list[str] = []
    for level_id in sorted(connection_map.keys()):
        entries = connection_map[level_id]
        if len(entries) != 2:
            raise StructureException("Connection must be paired")

        s1, s2 = entries
        t1 = s1.split(":")
        t2 = s2.split(":")
        rnum1 = int(t1[2].replace("R", ""))
        rnum2 = int(t2[2].replace("R", ""))

        if rnum1 > rnum2:
            connection = (
                f"PEPTIDE{t1[0]},PEPTIDE{t2[0]},{t1[1]}:{t1[2]}-{t2[1]}:{t2[2]}"
            )
        else:
            connection = (
                f"PEPTIDE{t2[0]},PEPTIDE{t1[0]},{t2[1]}:{t2[2]}-{t1[1]}:{t1[2]}"
            )
        connections.append(connection)

    # Build final notation
    sb: list[str] = []
    for i, notation in enumerate(simple_notations):
        if i > 0:
            sb.append("|")
        sb.append(f"PEPTIDE{i + 1}{{{notation}}}")

    sb.append("$")
    if connections:
        sb.append("|".join(connections))
    sb.append("$$$")

    return "".join(sb)
