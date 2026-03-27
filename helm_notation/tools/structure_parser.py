"""Chemical structure handling using RDKit (replaces ChemAxon MarvinBeans)."""

from __future__ import annotations

import re
from typing import Optional

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, RWMol, rdMolDescriptors

# Suppress "not removing hydrogen atom with dummy atom neighbors" warnings
# that occur during CXSMILES parsing — these are expected and harmless.
RDLogger.logger().setLevel(RDLogger.ERROR)

from helm_notation.exceptions import NotationException, StructureException
from helm_notation.model.molecule_info import MoleculeInfo

# Extended SMILES format constants
SMILES_EXTENSION_SEPARATOR_REGEX = r"\|"
SMILES_EXTENSION_SEPARATOR_SYMBOL = "|"
EXTENSION_COMPONENT_SEPARATOR_REGEX = r"[$]"
EXTENSION_COMPONENT_SEPARATOR_SYMBOL = "$"
SMILES_MIXTURE_DELIMITER_REGEX = r"\."
ATOM_POSITION_DELIMITER_REGEX = r"\$|;"
ATOM_POSITION_DELIMITER_SYMBOL = ";"
WHITE_SPACE_REGEX = r"\s"

CANONICALIZATION_ROUND = 3

_ADDITION = 1
_SUBTRACTION = 2


def get_molecule_info(smiles: str) -> MoleculeInfo:
    """Calculate molecular weight, formula, and exact mass from SMILES."""
    mol = get_molecule(smiles)
    if mol is None:
        raise StructureException(f"Cannot parse SMILES: {smiles}")
    mi = MoleculeInfo()
    mi.molecular_formula = rdMolDescriptors.CalcMolFormula(mol)
    mi.molecular_weight = round(Descriptors.MolWt(mol), 2)
    mi.exact_mass = round(Descriptors.ExactMolWt(mol), 2)
    return mi


def validate_smiles(smiles: str) -> bool:
    """Validate if the provided SMILES string is valid."""
    mol = get_molecule(smiles)
    if mol is None:
        return False
    problems = Chem.DetectChemistryProblems(mol)
    return len(problems) == 0


def get_simple_extended_smiles(extended_smiles: str) -> str:
    """Remove everything in ChemAxon Extended SMILES extension except atom positions.

    Input:  [*][H].O[*].OC1=CC=C(...) |r,$_R1;;;_R2;...$,c:14,t:3,5|
    Output: [*][H].O[*].OC1=CC=C(...) |$_R1;;;_R2;...$|
    """
    components = re.split(SMILES_EXTENSION_SEPARATOR_REGEX, extended_smiles)

    # Split on | gives: [smiles_part, extension_content, trailing_empty]
    # or [smiles_part, extension_content] depending on trailing |
    if len(components) < 2:
        raise StructureException("Invalid extended SMILES format")

    smi = components[0]
    # Rejoin in case there were extra | chars, strip trailing empty
    extension = "|".join(components[1:]).rstrip("|").strip()

    ext_components = re.split(EXTENSION_COMPONENT_SEPARATOR_REGEX, extension)
    if len(ext_components) < 3:
        raise StructureException("Invalid extended SMILES format")

    ext = ext_components[1]
    return (
        f"{smi}{SMILES_EXTENSION_SEPARATOR_SYMBOL}"
        f"{EXTENSION_COMPONENT_SEPARATOR_SYMBOL}{ext}"
        f"{EXTENSION_COMPONENT_SEPARATOR_SYMBOL}"
        f"{SMILES_EXTENSION_SEPARATOR_SYMBOL}"
    )


def get_unique_extended_smiles(extended_smiles: str) -> str:
    """Return unique canonical extended SMILES (with R-group atom labels)."""
    mol = get_molecule(extended_smiles)
    if mol is None:
        raise StructureException(f"Cannot parse SMILES: {extended_smiles}")
    return get_unique_extended_smiles_from_mol(mol)


def get_unique_extended_smiles_from_mol(mol: Chem.Mol) -> str:
    """Return unique canonical extended SMILES from an RDKit Mol."""
    return _get_unique_smiles_impl(mol, use_extended=True)


def get_unique_smiles(smiles: str) -> str:
    """Return unique canonical SMILES."""
    mol = get_molecule(smiles)
    if mol is None:
        raise StructureException(f"Cannot parse SMILES: {smiles}")
    return get_unique_smiles_from_mol(mol)


def get_unique_smiles_from_mol(mol: Chem.Mol) -> str:
    """Return unique canonical SMILES from an RDKit Mol."""
    return _get_unique_smiles_impl(mol, use_extended=False)


def _get_unique_smiles_impl(mol: Chem.Mol, use_extended: bool = False) -> str:
    """Internal: canonicalize SMILES, optionally preserving CXSMILES extensions."""
    mol = Chem.RemoveHs(mol)
    if use_extended:
        new_smi = Chem.MolToCXSmiles(mol, Chem.SmilesWriteParams())
    else:
        new_smi = Chem.MolToSmiles(mol, canonical=True)

    old_smi = ""
    count = 0
    while old_smi != new_smi and count < CANONICALIZATION_ROUND:
        count += 1
        old_smi = new_smi
        reparsed = get_molecule(old_smi)
        if reparsed is None:
            break
        reparsed = Chem.RemoveHs(reparsed)
        if use_extended:
            new_smi = Chem.MolToCXSmiles(reparsed, Chem.SmilesWriteParams())
        else:
            new_smi = Chem.MolToSmiles(reparsed, canonical=True)

    return new_smi


def get_molecule(structure_input: str) -> Chem.Mol | None:
    """Convert SMILES or molfile string to an RDKit Mol object.

    Handles both plain SMILES, ChemAxon extended SMILES (cxsmiles with
    |$...$| atom labels), and V2000/V3000 molfile blocks.
    """
    if structure_input is None:
        return None

    stripped = structure_input.strip()

    # Detect molfile format
    if _is_molfile(stripped):
        mol = Chem.MolFromMolBlock(stripped, sanitize=True)
        if mol is not None:
            AllChem.Compute2DCoords(mol)
        return mol

    # Try parsing as CXSMILES first (handles |$...$| atom labels)
    mol = _parse_cxsmiles(stripped)
    if mol is not None:
        AllChem.Compute2DCoords(mol)
        return mol

    # Fall back to plain SMILES
    mol = Chem.MolFromSmiles(stripped, sanitize=True)
    if mol is not None:
        AllChem.Compute2DCoords(mol)
    return mol


def _is_molfile(text: str) -> bool:
    """Check if the input looks like a molfile/SDF block."""
    lines = text.split("\n")
    # V2000/V3000 molfiles have "V2000" or "V3000" in the counts line (line 4)
    if len(lines) >= 4:
        if "V2000" in lines[3] or "V3000" in lines[3]:
            return True
    # Also check for "M  END" which is present in all molfiles
    return "M  END" in text


def _parse_cxsmiles(smiles_str: str) -> Chem.Mol | None:
    """Parse ChemAxon extended SMILES with R-group labels.

    ChemAxon format: SMILES |$_R1;;;_R2;...$|
    The $...$ section maps atom labels to positions.

    Uses RDKit's native CXSMILES support to correctly resolve label positions,
    then converts _R1/_R2/_R3 atomLabel properties to isotope-labeled dummy
    atoms [1*], [2*], [3*].
    """
    # Check if this contains a CXSMILES extension
    if "|" not in smiles_str:
        return None

    # Use RDKit's native CXSMILES parser which handles position mapping
    mol = Chem.MolFromSmiles(smiles_str, sanitize=True)
    if mol is None:
        # Try without sanitization
        mol = Chem.MolFromSmiles(smiles_str, sanitize=False)
        if mol is None:
            return None
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            return None

    # Find atoms with atomLabel properties matching _R<n> pattern
    rw_mol = Chem.RWMol(mol)
    found_rgroup = False
    for atom in rw_mol.GetAtoms():
        if atom.HasProp("atomLabel"):
            label = atom.GetProp("atomLabel")
            if label.startswith("_R") and label[2:].isdigit():
                rgroup_id = int(label[2:])
                found_rgroup = True
                # If this atom is a dummy, set its isotope directly
                if atom.GetAtomicNum() == 0:
                    atom.SetIsotope(rgroup_id)
                else:
                    # The label is on a real atom adjacent to a dummy;
                    # find the neighboring dummy and label it instead.
                    # This handles cases where CXSMILES position mapping
                    # lands on the neighbor rather than the dummy.
                    labeled_dummy = False
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 0 and neighbor.GetIsotope() == 0:
                            neighbor.SetIsotope(rgroup_id)
                            labeled_dummy = True
                            break
                    if not labeled_dummy:
                        # As fallback, set isotope on the labeled atom itself
                        atom.SetIsotope(rgroup_id)

    if not found_rgroup:
        # No R-group labels found; this may not be a CXSMILES or labels
        # are in a format we don't handle. Return None to fall through.
        return None

    return rw_mol.GetMol()


def find_rgroup_atom(mol: Chem.Mol, rgroup_id: int) -> int:
    """Find the atom index of the dummy atom representing the given R-group.

    R-groups are represented as dummy atoms (atomic num 0) with isotope = rgroup_id.
    """
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0 and atom.GetIsotope() == rgroup_id:
            return atom.GetIdx()
    raise StructureException(f"R-group R{rgroup_id} does not exist in the structure")


def get_rgroup_atom(mol: Chem.Mol, rgroup_id: int) -> int:
    """Alias for find_rgroup_atom for API compatibility."""
    return find_rgroup_atom(mol, rgroup_id)


def remove_rgroup(mol: Chem.RWMol, rgroup_id: int) -> int:
    """Remove the R-group dummy atom, return the index of its neighbor atom.

    The mol is modified in place. The returned index is the atom that was
    connected to the R-group and will be used for bonding during merging.
    """
    r_idx = find_rgroup_atom(mol, rgroup_id)
    return remove_rgroup_atom(mol, r_idx)


def remove_rgroup_atom(mol: Chem.RWMol, r_atom_idx: int) -> int:
    """Remove a specific R-group atom by index, return its neighbor's index."""
    atom = mol.GetAtomWithIdx(r_atom_idx)
    if atom.GetDegree() != 1:
        raise StructureException("Only one bond is allowed from Rgroup")

    neighbor = atom.GetNeighbors()[0]
    neighbor_idx = neighbor.GetIdx()

    mol.RemoveAtom(r_atom_idx)

    # Adjust index since removing an atom shifts indices
    if neighbor_idx > r_atom_idx:
        neighbor_idx -= 1

    return neighbor_idx


def merge_ignore_stereo(
    mol1: Chem.Mol,
    r_atom_idx1: int,
    mol2: Chem.Mol,
    r_atom_idx2: int,
) -> Chem.RWMol:
    """Merge mol2 into mol1 by connecting at R-group positions (ignoring stereo).

    If mol1 is mol2 (same object), performs cyclization.
    """
    mol1 = Chem.RWMol(mol1) if not isinstance(mol1, Chem.RWMol) else mol1
    mol2 = Chem.RWMol(mol2) if not isinstance(mol2, Chem.RWMol) else mol2
    if mol1 is mol2:
        # Cyclization
        Chem.Kekulize(mol1, clearAromaticFlags=True)
        atom1_idx = remove_rgroup_atom(mol1, r_atom_idx1)
        # After removal, r_atom_idx2 may have shifted
        adjusted_idx2 = r_atom_idx2 - 1 if r_atom_idx2 > r_atom_idx1 else r_atom_idx2
        atom2_idx = remove_rgroup_atom(mol1, adjusted_idx2)
        mol1.AddBond(atom1_idx, atom2_idx, Chem.BondType.SINGLE)
        return mol1

    Chem.Kekulize(mol1, clearAromaticFlags=True)
    Chem.Kekulize(mol2, clearAromaticFlags=True)

    atom1_idx = remove_rgroup_atom(mol1, r_atom_idx1)
    atom2_idx = remove_rgroup_atom(mol2, r_atom_idx2)

    # Combine molecules
    combined = Chem.RWMol(Chem.CombineMols(mol1, mol2))
    # atom2 index is offset by mol1's atom count
    new_atom2_idx = atom2_idx + mol1.GetNumAtoms()

    combined.AddBond(atom1_idx, new_atom2_idx, Chem.BondType.SINGLE)
    return combined


def is_single_stereo_bond(mol: Chem.Mol, r_atom_idx: int) -> bool:
    """Check if the bond connecting the R-atom has stereo information."""
    atom = mol.GetAtomWithIdx(r_atom_idx)
    if atom.GetDegree() != 1:
        raise StructureException(
            "RGroup is allowed to have single connection to other atom"
        )
    bond = atom.GetBonds()[0]
    bond_dir = bond.GetBondDir()
    return bond_dir in (
        Chem.BondDir.BEGINWEDGE,
        Chem.BondDir.BEGINDASH,
        Chem.BondDir.UNKNOWN,
    )


def merge(
    mol1: Chem.Mol,
    r_atom_idx1: int,
    mol2: Chem.Mol,
    r_atom_idx2: int,
) -> Chem.RWMol:
    """Merge mol2 into mol1 at R-group positions, preserving stereo bonds."""
    mol1 = Chem.RWMol(mol1) if not isinstance(mol1, Chem.RWMol) else mol1
    mol2 = Chem.RWMol(mol2) if not isinstance(mol2, Chem.RWMol) else mol2
    stereo1 = is_single_stereo_bond(mol1, r_atom_idx1)
    stereo2 = is_single_stereo_bond(mol2, r_atom_idx2)

    if stereo1 and stereo2:
        raise StructureException(
            "Both R atoms are connected to chiral centers"
        )

    Chem.Kekulize(mol1, clearAromaticFlags=True)
    Chem.Kekulize(mol2, clearAromaticFlags=True)

    if stereo2:
        # Keep the stereo bond direction from mol2 side
        bond2 = mol2.GetAtomWithIdx(r_atom_idx2).GetBonds()[0]
        bond_dir = bond2.GetBondDir()
        bond_type = bond2.GetBondType()

        neighbor2_idx = remove_rgroup_atom(mol2, r_atom_idx2)
        neighbor1_idx = remove_rgroup_atom(mol1, r_atom_idx1)

        combined = Chem.RWMol(Chem.CombineMols(mol1, mol2))
        new_neighbor2 = neighbor2_idx + mol1.GetNumAtoms()
        bond_idx = combined.AddBond(neighbor1_idx, new_neighbor2, bond_type) - 1
        combined.GetBondWithIdx(bond_idx).SetBondDir(bond_dir)
        return combined
    else:
        # Keep the bond from mol1 side (regardless of stereo type)
        bond1 = mol1.GetAtomWithIdx(r_atom_idx1).GetBonds()[0]
        bond_dir = bond1.GetBondDir()
        bond_type = bond1.GetBondType()

        neighbor1_idx = remove_rgroup_atom(mol1, r_atom_idx1)
        neighbor2_idx = remove_rgroup_atom(mol2, r_atom_idx2)

        if mol1 is mol2:
            mol1.AddBond(neighbor1_idx, neighbor2_idx, bond_type)
            return mol1

        combined = Chem.RWMol(Chem.CombineMols(mol1, mol2))
        new_neighbor2 = neighbor2_idx + mol1.GetNumAtoms()
        bond_idx = combined.AddBond(neighbor1_idx, new_neighbor2, bond_type) - 1
        if bond_dir != Chem.BondDir.NONE:
            combined.GetBondWithIdx(bond_idx).SetBondDir(bond_dir)
        return combined


def process_molecule_info(
    chunks: list[MoleculeInfo], caps: list[MoleculeInfo]
) -> MoleculeInfo:
    """Compute net MoleculeInfo by adding chunks and subtracting caps."""
    result = MoleculeInfo()
    mw = 0.0
    exact = 0.0
    atom_number_map: dict[str, int] = {}

    for mi in chunks:
        mw += mi.molecular_weight
        exact += mi.exact_mass
        _process_molecule_formula(atom_number_map, mi.molecular_formula, _ADDITION)

    for mi in caps:
        mw -= mi.molecular_weight
        exact -= mi.exact_mass
        _process_molecule_formula(atom_number_map, mi.molecular_formula, _SUBTRACTION)

    formula_parts = []
    for atom in sorted(atom_number_map.keys()):
        num = atom_number_map[atom]
        formula_parts.append(f"{atom}{num}")

    result.molecular_weight = mw
    result.exact_mass = exact
    result.molecular_formula = "".join(formula_parts)
    return result


def _process_molecule_formula(
    atom_number_map: dict[str, int], formula: str, operation_type: int
) -> None:
    """Parse a molecular formula and add/subtract atom counts."""
    atom = ""
    number = ""

    for ch in formula:
        if ch.isupper():
            if atom:
                _update_atom_number_map(atom_number_map, atom, number, operation_type)
                number = ""
            atom = ch
        elif ch.islower():
            atom += ch
        else:  # digit
            number += ch

    if atom:
        _update_atom_number_map(atom_number_map, atom, number, operation_type)


def _update_atom_number_map(
    atom_number_map: dict[str, int],
    atom: str,
    number: str,
    operation_type: int,
) -> None:
    """Update atom count map with addition or subtraction."""
    num = int(number) if number else 1

    if operation_type == _ADDITION:
        atom_number_map[atom] = atom_number_map.get(atom, 0) + num
    else:
        if atom not in atom_number_map:
            raise NotationException(
                f"Atom {atom} exists in cap group but not in main structure"
            )
        atom_number_map[atom] = atom_number_map[atom] - num


def get_rgroups_from_extended_smiles(extended_smiles: str) -> list[str]:
    """Extract R-group labels (e.g., R1, R2) from an extended SMILES string."""
    result: list[str] = []
    tokens = extended_smiles.split("R")
    if len(tokens) > 1:
        for token in tokens[1:]:
            numbers = ""
            for ch in token:
                if ch.isdigit():
                    numbers += ch
                else:
                    break
            if numbers:
                result.append(f"R{numbers}")
    return result


def decompose_mixture(mixture_extended_smiles: str) -> list[str]:
    """Decompose extended SMILES for a mixture into individual component SMILES."""
    simple_ext = get_simple_extended_smiles(mixture_extended_smiles)

    main_and_ext = re.split(SMILES_EXTENSION_SEPARATOR_REGEX, simple_ext)
    main = re.sub(WHITE_SPACE_REGEX, "", main_and_ext[0])
    extension = main_and_ext[1]

    all_smiles = re.split(SMILES_MIXTURE_DELIMITER_REGEX, main)
    token_atoms = re.split(ATOM_POSITION_DELIMITER_REGEX, extension)

    all_atom_list: list[str] = []
    for i, ta in enumerate(token_atoms):
        if i == 0 or i == len(token_atoms) - 1:
            if len(ta) != 0:
                all_atom_list.append(ta)
        else:
            all_atom_list.append(ta)

    pos = 0
    smiles_ext_list: list[str] = []

    for smiles in all_smiles:
        mol = get_molecule(smiles)
        if mol is None:
            raise StructureException(f"Cannot parse SMILES component: {smiles}")
        atom_count = mol.GetNumAtoms()

        parts = [SMILES_EXTENSION_SEPARATOR_SYMBOL + EXTENSION_COMPONENT_SEPARATOR_SYMBOL]
        for j in range(atom_count):
            parts.append(all_atom_list[pos])
            parts.append(ATOM_POSITION_DELIMITER_SYMBOL)
            pos += 1

        ext = "".join(parts)
        # Replace last delimiter with closing
        ext = ext[:-1] + EXTENSION_COMPONENT_SEPARATOR_SYMBOL + SMILES_EXTENSION_SEPARATOR_SYMBOL
        smiles_ext_list.append(f"{smiles} {ext}")

    return smiles_ext_list
