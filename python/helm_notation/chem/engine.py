from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple, Union

from helm_notation.model import MoleculeInfo

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
except Exception as exc:  # pragma: no cover
    raise ImportError(
        "RDKit is required for helm_notation. Install with `pip install rdkit-pypi`."
    ) from exc


class EngineError(ValueError):
    pass


MolLike = Union[str, "Chem.Mol"]


@dataclass
class RemovedRGroup:
    molecule: "Chem.Mol"
    attachment_atom_idx: int


class ChemEngine:
    SMILES_EXTENSION_SEPARATOR_REGEX = "|"
    EXTENSION_COMPONENT_SEPARATOR_REGEX = "$"
    SMILES_MIXTURE_DELIMITER = "."

    @staticmethod
    def get_molecule_info(smiles: str) -> MoleculeInfo:
        mol = ChemEngine.get_molecule(smiles)
        return MoleculeInfo(
            molecular_weight=Descriptors.MolWt(mol),
            exact_mass=Descriptors.ExactMolWt(mol),
            molecular_formula=rdMolDescriptors.CalcMolFormula(mol),
        )

    @staticmethod
    def validate_smiles(smiles: str) -> bool:
        try:
            ChemEngine.get_molecule(smiles)
            return True
        except EngineError:
            return False

    @staticmethod
    def get_simple_extended_smiles(extended_smiles: str) -> str:
        if "|" not in extended_smiles:
            raise EngineError("Invalid extended SMILES format")
        parts = extended_smiles.split("|")
        if len(parts) != 3:
            raise EngineError("Invalid extended SMILES format")

        smiles = parts[0]
        extension = parts[1]
        ext_parts = extension.split("$", -1)
        if len(ext_parts) != 3:
            raise EngineError("Invalid extended SMILES format")

        atom_positions = ext_parts[1]
        return f"{smiles}|${atom_positions}$|"

    @staticmethod
    def get_unique_extended_smiles(mol_or_smiles: MolLike) -> str:
        mol = ChemEngine.get_molecule(mol_or_smiles)
        return Chem.MolToCXSmiles(mol, canonical=True)

    @staticmethod
    def get_unique_smiles(mol_or_smiles: MolLike) -> str:
        mol = ChemEngine.get_molecule(mol_or_smiles)
        return Chem.MolToSmiles(mol, canonical=True)

    @staticmethod
    def smiles_to_pdb(smiles: str) -> str:
        mol = Chem.AddHs(ChemEngine.get_molecule(smiles))
        return Chem.MolToPDBBlock(mol)

    @staticmethod
    def get_molecule(structure_input: MolLike) -> "Chem.Mol":
        if isinstance(structure_input, Chem.Mol):
            return Chem.Mol(structure_input)
        if not isinstance(structure_input, str) or not structure_input.strip():
            raise EngineError("Structure input is empty")

        text = structure_input.strip()
        mol = None

        if "M  END" in text or "\n" in text:
            mol = Chem.MolFromMolBlock(text, sanitize=True, removeHs=False)
        if mol is None:
            mol = Chem.MolFromSmiles(text, sanitize=True)

        if mol is None:
            raise EngineError("Unable to parse structure input as molfile or SMILES")
        return mol

    @staticmethod
    def decompose_mixture(mixture_extended_smiles: str) -> List[str]:
        mixture = ChemEngine.get_simple_extended_smiles(mixture_extended_smiles)
        smiles_part, extension_part, _ = mixture.split("|")
        smiles_tokens = [token for token in smiles_part.replace(" ", "").split(".") if token]

        ext_tokens = [token for token in extension_part.split("$") if token]
        if len(ext_tokens) != 1:
            raise EngineError("Invalid extended SMILES format")

        atom_tags = [tag for tag in ext_tokens[0].split(";")]

        pos = 0
        components: List[str] = []
        for token in smiles_tokens:
            mol = ChemEngine.get_molecule(token)
            atom_count = mol.GetNumAtoms()
            chunk = atom_tags[pos : pos + atom_count]
            pos += atom_count
            components.append(f"{token}|${';'.join(chunk)}$|")

        return components

    @staticmethod
    def get_rgroup_atom_idx(mol_or_smiles: MolLike, rgroup_id: int) -> int:
        mol = ChemEngine.get_molecule(mol_or_smiles)
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 0:
                continue
            label = ChemEngine._get_r_label(atom)
            if label == f"R{rgroup_id}":
                return atom.GetIdx()
            if atom.GetIsotope() == rgroup_id:
                return atom.GetIdx()
            if atom.GetAtomMapNum() == rgroup_id:
                return atom.GetIdx()
        raise EngineError("Rgroup does not exist in the structure")

    @staticmethod
    def remove_rgroup(mol_or_smiles: MolLike, rgroup_id: int) -> RemovedRGroup:
        mol = ChemEngine.get_molecule(mol_or_smiles)
        idx = ChemEngine.get_rgroup_atom_idx(mol, rgroup_id)
        return ChemEngine.remove_rgroup_by_atom_idx(mol, idx)

    @staticmethod
    def remove_rgroup_by_atom_idx(mol_or_smiles: MolLike, r_atom_idx: int) -> RemovedRGroup:
        mol = ChemEngine.get_molecule(mol_or_smiles)
        rw = Chem.RWMol(mol)
        r_atom = rw.GetAtomWithIdx(r_atom_idx)
        neighbors = [n.GetIdx() for n in r_atom.GetNeighbors()]
        if len(neighbors) != 1:
            raise EngineError("Only one bond is allowed from Rgroup")
        attachment_atom_idx = neighbors[0]

        rw.RemoveAtom(r_atom_idx)
        if r_atom_idx < attachment_atom_idx:
            attachment_atom_idx -= 1

        new_mol = rw.GetMol()
        Chem.SanitizeMol(new_mol)
        return RemovedRGroup(new_mol, attachment_atom_idx)

    @staticmethod
    def merge(
        left_mol_or_smiles: MolLike,
        left_rgroup_id: int,
        right_mol_or_smiles: MolLike,
        right_rgroup_id: int,
    ) -> "Chem.Mol":
        left_removed = ChemEngine.remove_rgroup(left_mol_or_smiles, left_rgroup_id)
        right_removed = ChemEngine.remove_rgroup(right_mol_or_smiles, right_rgroup_id)

        left = left_removed.molecule
        right = right_removed.molecule

        combined = Chem.CombineMols(left, right)
        rw = Chem.RWMol(combined)

        left_idx = left_removed.attachment_atom_idx
        right_idx = left.GetNumAtoms() + right_removed.attachment_atom_idx

        rw.AddBond(left_idx, right_idx, Chem.BondType.SINGLE)
        merged = rw.GetMol()
        Chem.SanitizeMol(merged)
        return merged

    @staticmethod
    def _get_r_label(atom: "Chem.Atom") -> Optional[str]:
        for key in ("atomLabel", "_MolFileRLabel", "_TriposAtomName"):
            if atom.HasProp(key):
                value = atom.GetProp(key)
                if not value:
                    continue
                normalized = value.lstrip("_")
                if normalized.startswith("R"):
                    return normalized
        return None

    @staticmethod
    def get_rgroups_from_extended_smiles(extended_smiles: str) -> List[str]:
        pattern = re.compile(r"R(\d+)")
        matches = pattern.findall(extended_smiles)
        unique_rgroups = sorted(set(f"R{m}" for m in matches), key=lambda x: int(x[1]))
        return unique_rgroups

    @staticmethod
    def process_molecule_info(
        chunks: List[MoleculeInfo],
        caps: List[MoleculeInfo],
    ) -> MoleculeInfo:
        total_mw = 0.0
        total_exact = 0.0
        atom_number_map: Dict[str, int] = {}

        for mi in chunks:
            total_mw += mi.molecular_weight
            total_exact += mi.exact_mass
            ChemEngine._process_molecular_formula(
                atom_number_map, mi.molecular_formula, 1
            )

        for mi in caps:
            total_mw -= mi.molecular_weight
            total_exact -= mi.exact_mass
            ChemEngine._process_molecular_formula(
                atom_number_map, mi.molecular_formula, -1
            )

        formula_parts = []
        for atom in sorted(atom_number_map.keys()):
            count = atom_number_map[atom]
            formula_parts.append(atom)
            if count > 1:
                formula_parts.append(str(count))

        return MoleculeInfo(
            molecular_weight=total_mw,
            exact_mass=total_exact,
            molecular_formula="".join(formula_parts),
        )

    @staticmethod
    def _process_molecular_formula(
        atom_number_map: Dict[str, int],
        formula: str,
        operation: int,
    ) -> None:
        import re
        pattern = re.compile(r"([A-Z][a-z]?)(\d*)")
        for match in pattern.finditer(formula):
            atom = match.group(1)
            count_str = match.group(2)
            count = int(count_str) if count_str else 1
            count = count * operation

            if atom in atom_number_map:
                atom_number_map[atom] += count
            else:
                atom_number_map[atom] = count
