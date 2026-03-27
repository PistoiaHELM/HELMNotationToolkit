"""Parser for simple polymer notation strings.

This module provides functions that handle simple polymer notation. It is
unlikely that they will be called by clients directly.

Converted from org.helm.notation.tools.SimpleNotationParser (Java).

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

import copy
import re
from typing import TYPE_CHECKING

from rdkit import Chem

from helm_notation.constants import MONOMER_GROUP_COUNT_INTERVAL
from helm_notation.exceptions import (
    MonomerException,
    NotationException,
    StructureException,
)
from helm_notation.model.attachment import (
    BACKBONE_MONOMER_BRANCH_ATTACHMENT,
    BACKBONE_MONOMER_LEFT_ATTACHMENT,
    BACKBONE_MONOMER_RIGHT_ATTACHMENT,
    BRANCH_MONOMER_ATTACHMENT,
    Attachment,
)
from helm_notation.model.molecule_info import MoleculeInfo
from helm_notation.model.monomer import (
    BACKBONE_MONOMER_TYPE,
    BRANCH_MONOMER_TYPE,
    CHEMICAL_POLYMER_TYPE,
    NUCLEIC_ACID_POLYMER_TYPE,
    PEPTIDE_POLYMER_TYPE,
    UNDEFINED_MONOMER_TYPE,
    ID_R,
    Monomer,
)
from helm_notation.model.nucleotide import (
    ENDING_POSITION_TYPE,
    STARTING_POSITION_TYPE,
    Nucleotide,
)
from helm_notation.model.rgroup_structure import RgroupStructure
from helm_notation.tools.structure_parser import (
    find_rgroup_atom,
    get_molecule,
    get_molecule_info,
    get_unique_extended_smiles,
    merge,
    process_molecule_info,
    validate_smiles,
)

if TYPE_CHECKING:
    from helm_notation.monomer_store import MonomerStore

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

GROUP_LEVEL_DELIMITER_REGEX = r"\."
GROUP_LEVEL_DELIMITER = "."
MODIFICATION_START_SYMBOL = "["
MODIFICATION_END_SYMBOL = "]"

MODIFICATION_DELIMITER_REGEX = r"[\[|\]]"
BRANCH_START_SYMBOL = "("
BRANCH_END_SYMBOL = ")"

BRANCH_DELIMITER_REGEX = r"[(|)]"
RNA_SENSE_STRAND_ANNOTATION = "ss"
RNA_ANTISENSE_STRAND_ANNOTATION = "as"

# Module-level seed map for ad-hoc monomer ID generation
_seed_map: dict[str, int] = {}


# ---------------------------------------------------------------------------
# Helper: SimpleNotationGroupIterator (port of Java inner class)
# ---------------------------------------------------------------------------

class SimpleNotationGroupIterator:
    """Iterator class for groups within a simple notation string."""

    def __init__(self, notation_string: str) -> None:
        self.characters = list(notation_string)
        self.notation_string = notation_string
        self.position = 0

    def has_next_group(self) -> bool:
        return self.position < len(self.characters)

    def next_group(self) -> str:
        current_position = self.position

        while True:
            current_character = self.characters[current_position]
            if current_character == "[":
                current_position = get_matching_bracket_position(
                    self.characters, current_position, "[", "]"
                )
            elif current_character == "(":
                current_position = get_matching_bracket_position(
                    self.characters, current_position, "(", ")"
                )
            elif current_character != ".":
                current_position += 1

            if current_position < 0:
                current_position = len(self.characters)

            if current_position >= len(self.characters) or self.characters[current_position] == ".":
                break

        token = self.notation_string[self.position:current_position]
        self.position = current_position + 1
        return token


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------

def get_matching_bracket_position(
    characters: list[str],
    position: int,
    open_bracket: str,
    close_bracket: str,
) -> int:
    """Find the position of the matching closing bracket.

    Returns the index of the matching close bracket, or -1 if not found.
    """
    if position < (len(characters) - 1) and characters[position] == open_bracket:
        current_position = position
        opening_bracket_count = 1

        while True:
            current_position += 1
            current_character = characters[current_position]
            if current_character == open_bracket:
                opening_bracket_count += 1
            elif current_character == close_bracket:
                opening_bracket_count -= 1

            if opening_bracket_count <= 0 or current_position >= (len(characters) - 1):
                break

        if characters[current_position] == close_bracket:
            return current_position
        else:
            return -1
    else:
        return -1


# ---------------------------------------------------------------------------
# Internal helpers to get the default monomer store
# ---------------------------------------------------------------------------

def _get_default_monomer_store() -> MonomerStore:
    """Return the global MonomerFactory's monomer store."""
    from helm_notation.monomer_factory import get_instance as get_monomer_factory

    try:
        factory = get_monomer_factory()
    except Exception as ex:
        raise NotationException("Unable to initialize monomer factory") from ex
    return factory.monomer_store


# ---------------------------------------------------------------------------
# SMILES generation
# ---------------------------------------------------------------------------

def get_simple_polymer_smiles(
    polymer_notation: str,
    polymer_type: str,
    monomer_store: MonomerStore | None = None,
) -> str:
    """Generate the SMILES string for simple polymer notation.

    Returns extended SMILES for the simple polymer notation, still containing
    R-groups.
    """
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()

    struc = get_simple_polymer_structure(polymer_notation, polymer_type, monomer_store)
    if polymer_type == CHEMICAL_POLYMER_TYPE:
        if struc.molecule is None:
            raise NotationException(
                "Polymer notation contain non-specifc monomer structure"
            )
    mol = struc.molecule
    return get_unique_extended_smiles(mol)


# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------

def validate_simple_notation_for_chem(
    polymer_notation: str,
    monomer_store: MonomerStore | None = None,
) -> bool:
    """Validate chem simple notation. Returns True or raises exception."""
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()
    return validate_simple_notation(polymer_notation, CHEMICAL_POLYMER_TYPE, monomer_store)


def validate_simple_notation_for_peptide(
    polymer_notation: str,
    monomer_store: MonomerStore | None = None,
) -> bool:
    """Validate peptide simple notation. Returns True or raises exception."""
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()
    return validate_simple_notation(polymer_notation, PEPTIDE_POLYMER_TYPE, monomer_store)


def validate_simple_notation_for_rna(
    polymer_notation: str,
    monomer_store: MonomerStore | None = None,
) -> bool:
    """Validate RNA simple notation. Returns True or raises exception."""
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()
    return validate_simple_notation(polymer_notation, NUCLEIC_ACID_POLYMER_TYPE, monomer_store)


def validate_simple_notation(
    polymer_notation: str,
    polymer_type: str,
    monomer_store: MonomerStore | None = None,
) -> bool:
    """Check if simple polymer notation is valid. Returns True or raises exception."""
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()
    get_simple_polymer_structure(polymer_notation, polymer_type, monomer_store)
    return True


# ---------------------------------------------------------------------------
# Structure generation
# ---------------------------------------------------------------------------

def get_simple_polymer_structure(
    polymer_notation: str,
    polymer_type: str,
    monomer_store: MonomerStore | None = None,
) -> RgroupStructure:
    """Generate the RgroupStructure of simple polymer notation."""
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()

    monomer_list = get_monomer_list(polymer_notation, polymer_type, monomer_store)
    if not monomer_list:
        raise NotationException("Polymer notation contains no monomer")

    structure_list = _get_monomer_structure_list(monomer_list)

    if len(monomer_list) == len(structure_list):
        molecule: Chem.RWMol | None = None
        rmap: dict[str, int] = {}

        # Assign globally unique isotopes to all R-group dummy atoms so that
        # after CombineMols we can locate them by isotope without collisions.
        isotope_maps = _assign_unique_isotopes(structure_list)

        if polymer_type == CHEMICAL_POLYMER_TYPE:
            if monomer_list[0].can_smiles is not None:
                molecule = structure_list[0].molecule
                imap = isotope_maps[0]
                for label, uid in imap.items():
                    rmap[f"1:{label}"] = _find_rgroup_in_mol(molecule, uid)
        else:
            prev_monomer_index = -1
            # Track unique isotope IDs still available on the previous backbone
            prev_remaining: dict[str, int] = {}  # label -> unique_isotope

            for i in range(len(monomer_list)):
                m = monomer_list[i]
                ms = structure_list[i]
                mol = ms.molecule
                imap = isotope_maps[i]  # label -> unique_isotope for this monomer
                current_remaining = dict(imap)  # copy; we'll remove used ones

                if molecule is not None:
                    if m.monomer_type == BACKBONE_MONOMER_TYPE:
                        # Merge: prev's R2 connects to current's R1
                        prev_uid = prev_remaining[BACKBONE_MONOMER_RIGHT_ATTACHMENT]
                        cur_uid = current_remaining[BACKBONE_MONOMER_LEFT_ATTACHMENT]
                        prev_r_idx = _find_rgroup_in_mol(molecule, prev_uid)
                        cur_r_idx = _find_rgroup_in_mol(mol, cur_uid)
                        molecule = merge(molecule, prev_r_idx, mol, cur_r_idx)

                        del prev_remaining[BACKBONE_MONOMER_RIGHT_ATTACHMENT]
                        del current_remaining[BACKBONE_MONOMER_LEFT_ATTACHMENT]

                        # Record remaining R-groups from previous monomer
                        monomer_count = i
                        for label, uid in prev_remaining.items():
                            rmap[f"{monomer_count}:{label}"] = _find_rgroup_in_mol(molecule, uid)

                        prev_monomer_index = i
                        prev_remaining = current_remaining

                    elif m.monomer_type == BRANCH_MONOMER_TYPE:
                        prev_uid = prev_remaining[BACKBONE_MONOMER_BRANCH_ATTACHMENT]
                        cur_uid = current_remaining[BRANCH_MONOMER_ATTACHMENT]
                        prev_r_idx = _find_rgroup_in_mol(molecule, prev_uid)
                        cur_r_idx = _find_rgroup_in_mol(mol, cur_uid)
                        molecule = merge(molecule, prev_r_idx, mol, cur_r_idx)

                        del prev_remaining[BACKBONE_MONOMER_BRANCH_ATTACHMENT]
                        del current_remaining[BRANCH_MONOMER_ATTACHMENT]

                        # Record remaining R-groups from branch monomer
                        monomer_count = i + 1
                        for label, uid in current_remaining.items():
                            rmap[f"{monomer_count}:{label}"] = _find_rgroup_in_mol(molecule, uid)
                    else:
                        raise NotationException(
                            "Undefined Monomer Type is not supported in simple polymer"
                        )
                else:
                    # first monomer
                    molecule = mol
                    prev_monomer_index = i
                    if BACKBONE_MONOMER_LEFT_ATTACHMENT in current_remaining:
                        left_uid = current_remaining[BACKBONE_MONOMER_LEFT_ATTACHMENT]
                        rmap[f"1:{BACKBONE_MONOMER_LEFT_ATTACHMENT}"] = _find_rgroup_in_mol(
                            molecule, left_uid
                        )
                        del current_remaining[BACKBONE_MONOMER_LEFT_ATTACHMENT]
                    prev_remaining = current_remaining

            # Record remaining R-groups from the last backbone monomer
            monomer_count = prev_monomer_index + 1
            for label, uid in prev_remaining.items():
                rmap[f"{monomer_count}:{label}"] = _find_rgroup_in_mol(molecule, uid)

        structure = RgroupStructure()
        structure.molecule = molecule
        structure.rgroup_map = rmap
        return structure
    else:
        raise NotationException(
            "The number of monomers and structures do not match"
        )


# ---------------------------------------------------------------------------
# Monomer lookup
# ---------------------------------------------------------------------------

def get_monomer(
    monomer_id: str,
    polymer_type: str,
    monomer_store: MonomerStore | None = None,
) -> Monomer:
    """Return a Monomer object based on monomer ID and polymer type."""
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()

    monomers = monomer_store.get_monomers(polymer_type)
    if monomers is None or monomer_id not in monomers:
        raise NotationException(
            f"Unknow monomerID '{monomer_id}' for {polymer_type}"
        )
    return monomers[monomer_id]


# ---------------------------------------------------------------------------
# Monomer structure conversion
# ---------------------------------------------------------------------------

def _get_monomer_structure_list(monomer_list: list[Monomer]) -> list[RgroupStructure]:
    """Convert list of Monomer to list of RgroupStructure."""
    result: list[RgroupStructure] = []
    for monomer in monomer_list:
        ms = _get_monomer_structure(monomer)
        result.append(ms)
    return result


def _get_monomer_structure(monomer: Monomer) -> RgroupStructure:
    """Return the RgroupStructure for a single monomer."""
    ms: RgroupStructure | None = None
    structure_input: str | None = None

    # Prefer SMILES over molfile: molfiles from the monomer DB use ChemAxon
    # R# atom notation which RDKit cannot parse, whereas CXSMILES with R-group
    # labels are handled correctly by our _parse_cxsmiles function.
    if monomer.can_smiles:
        structure_input = monomer.can_smiles

    if structure_input is None and monomer.molfile:
        structure_input = monomer.molfile

    if structure_input is not None:
        ms = RgroupStructure()
        mol = get_molecule(structure_input)
        if mol is not None:
            # Dearomatize (Kekulize) the molecule
            try:
                Chem.Kekulize(mol, clearAromaticFlags=True)
            except Exception:
                pass

            # Build rgroup_map: R-group label -> atom index.
            # R-groups are represented as dummy atoms (atomic num 0) with
            # isotope = R-group number.
            rgroup_map: dict[str, int] = {}
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() == 0:
                    r_id = atom.GetIsotope()
                    if r_id > 0:
                        rgroup_map[f"R{r_id}"] = atom.GetIdx()

            ms.molecule = Chem.RWMol(mol) if not isinstance(mol, Chem.RWMol) else mol
            ms.rgroup_map = rgroup_map

    if ms is None:
        ms = RgroupStructure()
    return ms


def _find_rgroup_in_mol(mol: Chem.Mol, unique_isotope: int) -> int:
    """Find the atom index of an R-group dummy atom by its unique isotope.

    After CombineMols, atom indices shift, so we must search by isotope.
    """
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0 and atom.GetIsotope() == unique_isotope:
            return atom.GetIdx()
    raise NotationException(f"R-group with isotope {unique_isotope} not found in molecule")


def _assign_unique_isotopes(
    structure_list: list[RgroupStructure],
) -> list[dict[str, int]]:
    """Assign globally unique isotope IDs to R-group dummy atoms across all structures.

    Each monomer's R-groups (R1, R2, R3) get unique isotopes to avoid collisions
    when molecules are merged. Returns a list of maps: rgroup_label -> unique_isotope
    for each structure.
    """
    next_isotope = 100  # Start at 100 to avoid small isotope conflicts
    result: list[dict[str, int]] = []

    for ms in structure_list:
        label_to_unique: dict[str, int] = {}
        if ms.molecule is not None:
            rw = Chem.RWMol(ms.molecule)
            for atom in rw.GetAtoms():
                if atom.GetAtomicNum() == 0 and atom.GetIsotope() > 0:
                    original_label = f"R{atom.GetIsotope()}"
                    atom.SetIsotope(next_isotope)
                    label_to_unique[original_label] = next_isotope
                    next_isotope += 1
            ms.molecule = rw
        result.append(label_to_unique)

    return result


# ---------------------------------------------------------------------------
# Monomer list generation
# ---------------------------------------------------------------------------

def get_monomer_list(
    polymer_notation: str,
    polymer_type: str,
    monomer_store: MonomerStore | None = None,
) -> list[Monomer]:
    """Generate the list of Monomer from polymer notation."""
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()

    result: list[Monomer] = []
    ids = get_monomer_id_list(polymer_notation, polymer_type, monomer_store)
    for mid in ids:
        m = get_monomer(mid, polymer_type, monomer_store)
        result.append(m)
    return result


# ---------------------------------------------------------------------------
# Nucleotide list generation
# ---------------------------------------------------------------------------

def get_nucleotide_list(
    polymer_notation: str,
    validate: bool = True,
    monomer_store: MonomerStore | None = None,
    include_last: bool = True,
) -> list[Nucleotide]:
    """Return the nucleotide list for RNA polymer type.

    Validation of polymer notation could be slow; set *validate* to False
    to skip structure validation.
    """
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()

    if validate:
        validate_simple_notation_for_rna(polymer_notation, monomer_store)

    from helm_notation.nucleotide_factory import get_instance as get_nucleotide_factory

    ids: list[Nucleotide] = []
    reverse_nuc_map = get_nucleotide_factory().get_reverse_nucleotide_template_map()

    iterator = SimpleNotationGroupIterator(polymer_notation)
    while iterator.has_next_group():
        notation = iterator.next_group()
        symbol: str | None = None

        # last nucleotide will be handled differently
        tmp_notation = notation
        if not iterator.has_next_group() and notation.endswith(")"):
            tmp_notation = notation + "P"

        if tmp_notation in reverse_nuc_map:
            symbol = reverse_nuc_map[tmp_notation]
        else:
            chars = list(notation)
            base: str | None = None
            symbol = "X"

            # find base
            j = 0
            while j < len(chars):
                letter = chars[j]
                # skip modifications if not in branch
                if letter == MODIFICATION_START_SYMBOL:
                    matching_pos = get_matching_bracket_position(
                        chars, j, MODIFICATION_START_SYMBOL, MODIFICATION_END_SYMBOL
                    )
                    j += 1
                    if matching_pos == -1:
                        raise NotationException(
                            "Invalid Polymer Notation: Could not find matching bracket"
                        )
                    j = matching_pos
                # base is always a branch monomer
                elif letter == BRANCH_START_SYMBOL:
                    matching_pos = get_matching_bracket_position(
                        chars, j, BRANCH_START_SYMBOL, BRANCH_END_SYMBOL
                    )
                    j += 1
                    if matching_pos == -1:
                        raise NotationException(
                            "Invalid Polymer Notation: Could not find matching bracket"
                        )

                    base = notation[j:matching_pos]

                    if len(base) == 1:
                        symbol = base
                    else:
                        processed_id = process_node(
                            base, NUCLEIC_ACID_POLYMER_TYPE, "X", monomer_store
                        )
                        monomer_obj = get_monomer(
                            processed_id, NUCLEIC_ACID_POLYMER_TYPE, monomer_store
                        )
                        if monomer_obj.natural_analog is None:
                            symbol = "X"
                        else:
                            symbol = monomer_obj.natural_analog

                    j = matching_pos
                j += 1

        nuc = Nucleotide(symbol=symbol or "X", notation=notation)
        ids.append(nuc)

    if len(ids) > 0:
        ids[0].position_type = STARTING_POSITION_TYPE

    if len(ids) > 1:
        ids[-1].position_type = ENDING_POSITION_TYPE

    return ids


def get_strict_nucleotide_list(
    polymer_notation: str,
    validate: bool = True,
    monomer_store: MonomerStore | None = None,
    include_last: bool = True,
) -> list[Nucleotide]:
    """Return a strict nucleotide list for RNA polymer type.

    This is more strict than get_nucleotide_list in that it will utilize
    modX for a natural sequence X or include ? for unknown bases.
    """
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()

    if validate:
        validate_simple_notation_for_rna(polymer_notation, monomer_store)

    from helm_notation.nucleotide_factory import get_instance as get_nucleotide_factory

    ids: list[Nucleotide] = []
    reverse_nuc_map = get_nucleotide_factory().get_reverse_nucleotide_template_map()

    iterator = SimpleNotationGroupIterator(polymer_notation)
    while iterator.has_next_group():
        notation = iterator.next_group()
        symbol: str | None = None

        nuc: Nucleotide
        if notation in reverse_nuc_map:
            symbol = reverse_nuc_map[notation]
            nuc = Nucleotide(symbol=symbol, notation=notation)
        else:
            nuc = Nucleotide(symbol="temp", notation=notation)
            natural_analog = nuc.get_natural_analog(monomer_store)
            if nuc.is_modified:
                if not iterator.has_next_group() and nuc.unmodified_without_phosphate():
                    nuc.symbol = f"end{natural_analog}"
                else:
                    nuc.symbol = f"mod{natural_analog}"
            else:
                nuc.symbol = natural_analog
        ids.append(nuc)

    if len(ids) > 0:
        ids[0].position_type = STARTING_POSITION_TYPE

    if len(ids) > 1:
        ids[-1].position_type = ENDING_POSITION_TYPE

    return ids


# ---------------------------------------------------------------------------
# Nucleotide monomer list (private)
# ---------------------------------------------------------------------------

def _get_nucleotide_monomer_list(
    nucleotide_sequence: str,
) -> list[tuple[str, str]]:
    """Return list of (monomer_type, id) tuples from a nucleotide sequence."""
    ids: list[tuple[str, str]] = []
    chars = list(nucleotide_sequence)
    prev_letter = ""

    i = 0
    while i < len(chars):
        letter = chars[i]
        if letter == MODIFICATION_START_SYMBOL:
            matching_pos = get_matching_bracket_position(
                chars, i, MODIFICATION_START_SYMBOL, MODIFICATION_END_SYMBOL
            )
            i += 1

            if matching_pos == -1:
                raise NotationException(
                    "Invalid Polymer Notation: modified monomer must be enclosed by square brackets"
                )
            else:
                ids.append((BACKBONE_MONOMER_TYPE, nucleotide_sequence[i:matching_pos]))
            i = matching_pos

        elif letter == BRANCH_START_SYMBOL:
            if i == 0:
                raise NotationException(
                    "Invalid Polymer Notation: branch monomer is not allowed at the beginnig of notation"
                )

            if prev_letter == BRANCH_END_SYMBOL:
                raise NotationException(
                    "Invalid Polymer Notation: branch monomers cannot be connected with each other"
                )

            matching_pos = get_matching_bracket_position(
                chars, i, BRANCH_START_SYMBOL, BRANCH_END_SYMBOL
            )
            i += 1
            if matching_pos == -1:
                raise NotationException(
                    "Invalid Polymer Notation: modified monomer must be enclosed by brackets"
                )
            else:
                ids.append((BRANCH_MONOMER_TYPE, nucleotide_sequence[i:matching_pos]))
            i = matching_pos

        else:
            ids.append((BACKBONE_MONOMER_TYPE, nucleotide_sequence[i : i + 1]))

        prev_letter = letter
        i += 1

    return ids


# ---------------------------------------------------------------------------
# Monomer ID list generation
# ---------------------------------------------------------------------------

def get_monomer_id_list(
    polymer_notation: str,
    polymer_type: str,
    monomer_store: MonomerStore | None = None,
) -> list[str]:
    """Return the list of Monomer IDs from simple polymer notation."""
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()

    ids: list[str] = []

    group_iterator = SimpleNotationGroupIterator(polymer_notation)

    # CHEMICAL can have only one monomer
    if polymer_type == CHEMICAL_POLYMER_TYPE:
        mid = process_node(polymer_notation, polymer_type, "X", monomer_store)
        ids.append(mid)

    elif polymer_type == PEPTIDE_POLYMER_TYPE:
        while group_iterator.has_next_group():
            notation = group_iterator.next_group()
            chars = list(notation)

            cur_id: str | None = None

            i = 0
            while i < len(chars):
                letter = chars[i]
                if letter == MODIFICATION_START_SYMBOL:
                    matching_pos = get_matching_bracket_position(
                        chars, i, MODIFICATION_START_SYMBOL, MODIFICATION_END_SYMBOL
                    )
                    i += 1

                    if matching_pos == -1:
                        raise NotationException(
                            "Invalid Polymer Notation: modified monomer must be enclosed by square brackets"
                        )
                    else:
                        cur_id = notation[i:matching_pos]
                    i = matching_pos

                elif len(notation) == 1:
                    cur_id = notation
                else:
                    raise NotationException(f"Invalid Peptide Notation: {notation}")

                cur_id = process_node(cur_id, polymer_type, "X", monomer_store)
                ids.append(cur_id)
                i += 1

    elif polymer_type == NUCLEIC_ACID_POLYMER_TYPE:
        is_first_group = True
        while group_iterator.has_next_group():
            notation = group_iterator.next_group()
            is_last_group = not group_iterator.has_next_group()

            monomers = _get_nucleotide_monomer_list(notation)

            has_base = any(t[0] == BRANCH_MONOMER_TYPE for t in monomers)

            if len(monomers) == 1:
                if is_last_group:
                    ids.append(process_node(monomers[0][1], polymer_type, "R", monomer_store))
                elif is_first_group:
                    ids.append(process_node(monomers[0][1], polymer_type, "P", monomer_store))
                else:
                    raise NotationException(
                        f"Invalid Polymer Notation: Nucleotide sequence {notation} is incorrect"
                    )

            elif len(monomers) == 2:
                if is_last_group:
                    if has_base:
                        ids.append(process_node(monomers[0][1], polymer_type, "R", monomer_store))
                        ids.append(process_node(monomers[1][1], polymer_type, "X", monomer_store))
                    else:
                        ids.append(process_node(monomers[0][1], polymer_type, "R", monomer_store))
                        ids.append(process_node(monomers[1][1], polymer_type, "P", monomer_store))
                elif is_first_group:
                    if has_base:
                        ids.append(process_node(monomers[0][1], polymer_type, "X", monomer_store))
                        ids.append(process_node(monomers[1][1], polymer_type, "P", monomer_store))
                    else:
                        ids.append(process_node(monomers[0][1], polymer_type, "R", monomer_store))
                        ids.append(process_node(monomers[1][1], polymer_type, "P", monomer_store))
                elif not has_base:
                    ids.append(process_node(monomers[0][1], polymer_type, "R", monomer_store))
                    ids.append(process_node(monomers[1][1], polymer_type, "P", monomer_store))
                else:
                    raise NotationException(
                        f"Invalid Polymer Notation: Nucleotide sequence {notation} is incorrect"
                    )

            elif len(monomers) == 3:
                ids.append(process_node(monomers[0][1], polymer_type, "R", monomer_store))
                ids.append(process_node(monomers[1][1], polymer_type, "X", monomer_store))
                ids.append(process_node(monomers[2][1], polymer_type, "P", monomer_store))
            else:
                raise NotationException(
                    f"Invalid Polymer Notation: Nucleotide sequence {notation} is incorrect"
                )

            is_first_group = False

    return ids


# ---------------------------------------------------------------------------
# Node preprocessing (ad-hoc monomer creation)
# ---------------------------------------------------------------------------

def process_node(
    node_desc: str,
    polymer_type: str,
    natural_analog: str,
    monomer_store: MonomerStore,
) -> str:
    """Preprocess a monomer node.

    If the monomer is not contained in the monomer_store, the system will
    generate a temporary monomer and add it to the monomer store. The default
    attachment group will be Hydrogen, as R1-H, R2-H.

    Returns the monomer alternate_id.
    """
    # remove brackets
    if (
        len(node_desc) >= 2
        and node_desc[0] == MODIFICATION_START_SYMBOL
        and node_desc[-1] == MODIFICATION_END_SYMBOL
    ):
        node_desc = node_desc[1:-1]

    is_smiles_code = bool(re.search(r".*\$\|$", node_desc))

    if not is_smiles_code:
        monomers = monomer_store.get_monomers(polymer_type)
        if monomers is not None and node_desc in monomers:
            return node_desc

    smiles_db = monomer_store.smiles_monomer_db

    # validate smiles before adding temp. monomer (xhelm-63)
    try:
        if not validate_smiles(node_desc):
            raise NotationException(f"Unable to create ad hoc monomer for {node_desc}")
    except Exception:
        raise NotationException(f"Unable to create ad hoc monomer for {node_desc}")

    unique_smiles: str
    try:
        unique_smiles = get_unique_extended_smiles(node_desc)
    except Exception:
        unique_smiles = node_desc

    alternate_id: str | None = None
    if unique_smiles in smiles_db:
        alternate_id = smiles_db[unique_smiles].alternate_id
    else:
        from helm_notation.monomer_factory import get_instance as get_monomer_factory

        try:
            factory = get_monomer_factory()
        except Exception as ex:
            raise NotationException("Unable to initialize monomer factory") from ex

        alternate_id = generate_next_ad_hoc_monomer_id(polymer_type, monomer_store)

        attachment_db = factory.attachment_db
        r1h_att = attachment_db.get("R1-H")

        m: Monomer
        if polymer_type == CHEMICAL_POLYMER_TYPE:
            m = Monomer(
                polymer_type=polymer_type,
                monomer_type=UNDEFINED_MONOMER_TYPE,
                natural_analog=natural_analog,
                alternate_id=alternate_id,
            )
        elif polymer_type == NUCLEIC_ACID_POLYMER_TYPE:
            if natural_analog in ("P", "R"):
                m = Monomer(
                    polymer_type=polymer_type,
                    monomer_type=BACKBONE_MONOMER_TYPE,
                    natural_analog=natural_analog,
                    alternate_id=alternate_id,
                )
            else:
                m = Monomer(
                    polymer_type=polymer_type,
                    monomer_type=BRANCH_MONOMER_TYPE,
                    natural_analog=natural_analog,
                    alternate_id=alternate_id,
                )
        else:
            # Peptide
            m = Monomer(
                polymer_type=polymer_type,
                monomer_type=BACKBONE_MONOMER_TYPE,
                natural_analog=natural_analog,
                alternate_id=alternate_id,
            )

        m.ad_hoc_monomer = True
        m.can_smiles = unique_smiles
        m.name = "Dynamic"

        al: list[Attachment] = []
        start = 0
        pos = unique_smiles.find("R", start)
        number = ""
        while pos >= 0:
            pos += 1
            letter = unique_smiles[pos : pos + 1]
            while letter and re.match(r"\d", letter):
                number = number + letter
                pos += 1
                letter = unique_smiles[pos : pos + 1]

            try:
                if r1h_att is not None:
                    tmp_att = copy.deepcopy(r1h_att)
                    tmp_att.label = f"R{number}"
                    tmp_att.alternate_id = f"R{number}-H"
                    old_smi = tmp_att.cap_group_smiles
                    new_smi = old_smi.replace("R1", f"R{number}")
                    tmp_att.cap_group_smiles = new_smi
                    al.append(tmp_att)
            except Exception as ex:
                raise NotationException(
                    "Unable to create attachment by copying from attachment database"
                ) from ex

            start = pos
            pos = unique_smiles.find("R", start)
            number = ""

        m.attachment_list = al
        try:
            monomer_store.add_new_monomer(m)
        except Exception as ex:
            raise NotationException(
                "Unable to add adhoc new monomer into monomer databse"
            ) from ex

    return alternate_id


# ---------------------------------------------------------------------------
# Nucleotide sequence extraction
# ---------------------------------------------------------------------------

def get_trimmed_nucleotide_sequence(
    polymer_notation: str,
    monomer_store: MonomerStore | None = None,
) -> str:
    """Generate single letter (natural analog) nucleotide sequence for simple RNA
    notation, removing nucleotides without base from both ends.
    """
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()

    nuc_list = get_nucleotide_list(polymer_notation, True, monomer_store)

    start = 0
    na = nuc_list[start]
    while na.get_base_monomer(monomer_store) is None:
        start += 1
        na = nuc_list[start]

    end = len(nuc_list) - 1
    na = nuc_list[end]
    while na.get_base_monomer(monomer_store) is None:
        end -= 1
        na = nuc_list[end]

    sb: list[str] = []
    for i in range(start, end + 1):
        sb.append(nuc_list[i].get_natural_analog(monomer_store))
    return "".join(sb)


def get_nucleotide_sequence(
    polymer_notation: str | None = None,
    monomer_store: MonomerStore | None = None,
    nucleotide_list: list[Nucleotide] | None = None,
) -> str:
    """Generate single letter (natural analog) nucleotide sequence for simple RNA notation.

    Can be called with either polymer_notation or nucleotide_list.
    """
    if nucleotide_list is not None:
        return _nucleotide_sequence_from_list(nucleotide_list)

    if polymer_notation is None:
        raise NotationException("Either polymer_notation or nucleotide_list is required")

    if monomer_store is None:
        monomer_store = _get_default_monomer_store()

    nuc_list = get_nucleotide_list(polymer_notation, True, monomer_store)
    return _nucleotide_sequence_from_list(nuc_list)


def _nucleotide_sequence_from_list(nucleotide_list: list[Nucleotide]) -> str:
    """Generate single letter sequence from a nucleotide list."""
    sb: list[str] = []
    for nuc in nucleotide_list:
        sb.append(nuc.get_natural_analog())
    return "".join(sb)


def get_modified_nucleotide_sequence(
    polymer_notation: str | None = None,
    monomer_store: MonomerStore | None = None,
    nucleotide_list: list[Nucleotide] | None = None,
) -> str:
    """Convert RNA simple notation to modified nucleotide sequence.

    Can be called with either polymer_notation or nucleotide_list.
    """
    if nucleotide_list is not None:
        return _modified_nucleotide_sequence_from_list(nucleotide_list)

    if polymer_notation is None:
        raise NotationException("Either polymer_notation or nucleotide_list is required")

    if monomer_store is None:
        monomer_store = _get_default_monomer_store()

    nuc_list = get_nucleotide_list(polymer_notation, True, monomer_store)
    return _modified_nucleotide_sequence_from_list(nuc_list)


def _modified_nucleotide_sequence_from_list(nucleotide_list: list[Nucleotide]) -> str:
    """Generate modified nucleotide sequence from a nucleotide list."""
    sb: list[str] = []
    for nuc in nucleotide_list:
        sb.append(nuc.symbol)
    return "".join(sb)


# ---------------------------------------------------------------------------
# Complex notation conversion
# ---------------------------------------------------------------------------

def get_complex_notation_for_chem(
    simple_notation: str,
    monomer_store: MonomerStore | None = None,
) -> str:
    """Convert chemical simple notation to complex notation."""
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()
    return get_complex_notation(simple_notation, CHEMICAL_POLYMER_TYPE, monomer_store)


def get_complex_notation_for_peptide(
    simple_notation: str,
    monomer_store: MonomerStore | None = None,
) -> str:
    """Convert Peptide simple notation to complex notation."""
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()
    return get_complex_notation(simple_notation, PEPTIDE_POLYMER_TYPE, monomer_store)


def get_complex_notation_for_antisense_rna(
    simple_notation: str,
    monomer_store: MonomerStore | None = None,
) -> str:
    """Convert RNA simple notation to complex notation and annotate as antisense strand 'as'."""
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()
    notation = get_complex_notation_for_rna(simple_notation, monomer_store)
    return (
        notation[: len(notation) - 1]
        + NUCLEIC_ACID_POLYMER_TYPE
        + "1{"
        + RNA_ANTISENSE_STRAND_ANNOTATION
        + "}$"
    )


def get_complex_notation_for_sense_rna(
    simple_notation: str,
    monomer_store: MonomerStore | None = None,
) -> str:
    """Convert RNA simple notation to complex notation and annotate as sense strand 'ss'."""
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()
    notation = get_complex_notation_for_rna(simple_notation, monomer_store)
    return (
        notation[: len(notation) - 1]
        + NUCLEIC_ACID_POLYMER_TYPE
        + "1{"
        + RNA_SENSE_STRAND_ANNOTATION
        + "}$"
    )


def get_complex_notation_for_rna(
    simple_notation: str,
    monomer_store: MonomerStore | None = None,
) -> str:
    """Convert RNA simple notation to complex notation."""
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()
    return get_complex_notation(simple_notation, NUCLEIC_ACID_POLYMER_TYPE, monomer_store)


def get_complex_notation(
    simple_notation: str,
    polymer_type: str,
    monomer_store: MonomerStore | None = None,
) -> str:
    """Convert simple notation for a given polymer type to complex notation.

    Complex notation is needed to calculate molecule info such as MW and formula.
    """
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()
    validate_simple_notation(simple_notation, polymer_type, monomer_store)
    return f"{polymer_type}1{{{simple_notation}}}$$$$"


# ---------------------------------------------------------------------------
# Monomer count
# ---------------------------------------------------------------------------

def get_monomer_count_for_peptide(
    simple_notation: str,
    monomer_store: MonomerStore | None = None,
) -> int:
    """Return the total number of monomers in peptide simple notation."""
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()
    return get_monomer_count(simple_notation, PEPTIDE_POLYMER_TYPE, monomer_store)


def get_monomer_count_for_rna(
    simple_notation: str,
    monomer_store: MonomerStore | None = None,
) -> int:
    """Return the total number of monomers in RNA simple notation."""
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()
    return get_monomer_count(simple_notation, NUCLEIC_ACID_POLYMER_TYPE, monomer_store)


def get_monomer_count(
    simple_notation: str,
    polymer_type: str,
    monomer_store: MonomerStore | None = None,
) -> int:
    """Return the total number of monomers in a simple notation for a given polymer type."""
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()
    id_list = get_monomer_id_list(simple_notation, polymer_type, monomer_store)
    return len(id_list)


# ---------------------------------------------------------------------------
# Monomer replacement
# ---------------------------------------------------------------------------

def replace_monomer(
    simple_notation: str,
    polymer_type: str,
    existing_monomer_id: str,
    new_monomer_id: str,
    monomer_store: MonomerStore | None = None,
    validate: bool = True,
) -> str:
    """Replace existing monomer with new monomer in simple polymer notation."""
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()

    if validate:
        validate_monomer_replacement(
            polymer_type, existing_monomer_id, new_monomer_id, monomer_store
        )

    monomer_ids = get_monomer_id_list(simple_notation, polymer_type, monomer_store)
    for i in range(len(monomer_ids)):
        if monomer_ids[i] == existing_monomer_id:
            monomer_ids[i] = new_monomer_id

    return _get_simple_notation(monomer_ids, polymer_type, monomer_store)


def _get_simple_notation(
    monomer_id_list: list[str],
    polymer_type: str,
    monomer_store: MonomerStore | None = None,
) -> str:
    """Reconstruct simple notation from a list of monomer IDs."""
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()

    sb: list[str] = []

    if polymer_type == NUCLEIC_ACID_POLYMER_TYPE:
        for mid in monomer_id_list:
            monomer_obj = get_monomer(mid, polymer_type, monomer_store)
            monomer_type = monomer_obj.monomer_type
            natural_analog = monomer_obj.natural_analog

            if monomer_type == BACKBONE_MONOMER_TYPE:
                if natural_analog == ID_R:
                    if sb:
                        sb.append(GROUP_LEVEL_DELIMITER)
                if len(mid) > 1:
                    sb.append(MODIFICATION_START_SYMBOL)
                    sb.append(mid)
                    sb.append(MODIFICATION_END_SYMBOL)
                else:
                    sb.append(mid)
            else:
                sb.append(BRANCH_START_SYMBOL)
                if len(mid) > 1:
                    sb.append(MODIFICATION_START_SYMBOL)
                    sb.append(mid)
                    sb.append(MODIFICATION_END_SYMBOL)
                else:
                    sb.append(mid)
                sb.append(BRANCH_END_SYMBOL)

    elif polymer_type == PEPTIDE_POLYMER_TYPE:
        for mid in monomer_id_list:
            if sb:
                sb.append(GROUP_LEVEL_DELIMITER)
            if len(mid) > 1:
                sb.append(MODIFICATION_START_SYMBOL)
                sb.append(mid)
                sb.append(MODIFICATION_END_SYMBOL)
            else:
                sb.append(mid)

    elif polymer_type == CHEMICAL_POLYMER_TYPE:
        sb.append(monomer_id_list[0])

    return "".join(sb)


def get_simple_notation_from_nucleotide_list(
    nucleotide_list: list[Nucleotide],
) -> str:
    """Convert a nucleotide list back to simple notation."""
    sb: list[str] = []
    for nuc in nucleotide_list:
        if sb:
            sb.append(GROUP_LEVEL_DELIMITER)
        sb.append(nuc.notation)
    return "".join(sb)


# ---------------------------------------------------------------------------
# Monomer replacement validation
# ---------------------------------------------------------------------------

def validate_monomer_replacement(
    polymer_type: str,
    existing_monomer_id: str,
    new_monomer_id: str,
    monomer_store: MonomerStore | None = None,
) -> bool:
    """Validate that monomer replacement is valid."""
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()

    if not polymer_type:
        raise NotationException("Polymer type is required for monomer replacement")

    if not existing_monomer_id:
        raise NotationException("Existing monomer ID is required for monomer replacement")

    if not new_monomer_id:
        raise NotationException("New monomer ID is required for monomer replacement")

    monomers = monomer_store.get_monomers(polymer_type)
    if monomers is None or len(monomers) == 0:
        raise NotationException(f"Unknown polymer type [{polymer_type}] found")

    if existing_monomer_id not in monomers:
        raise NotationException(
            f"Existing monomer ID [{existing_monomer_id}] is invalid in polymer type {polymer_type}"
        )

    if new_monomer_id not in monomers:
        raise NotationException(
            f"New monomer ID [{new_monomer_id}] is invalid in polymer type {polymer_type}"
        )

    existing_monomer = monomers[existing_monomer_id]
    new_monomer = monomers[new_monomer_id]

    if polymer_type == NUCLEIC_ACID_POLYMER_TYPE:
        if existing_monomer.monomer_type == new_monomer.monomer_type:
            if existing_monomer.monomer_type == BACKBONE_MONOMER_TYPE:
                if existing_monomer.natural_analog != new_monomer.natural_analog:
                    raise NotationException(
                        f"Existing monomer natural analog [{existing_monomer.natural_analog}] "
                        f"and new monomer natural analog [{new_monomer.natural_analog}] are different"
                    )
        else:
            raise NotationException(
                f"Existing monomer type [{existing_monomer.monomer_type}] "
                f"and new monomer type [{new_monomer.monomer_type}] are different"
            )

    if not new_monomer.attachment_equals(existing_monomer):
        raise NotationException(
            f"Existing monomer attachment [{existing_monomer.attachment_list_string}] "
            f"and new monomer attachement [{new_monomer.attachment_list_string}] are different"
        )

    return True


# ---------------------------------------------------------------------------
# Peptide sequence extraction
# ---------------------------------------------------------------------------

def get_peptide_sequence(
    polymer_notation: str,
    monomer_store: MonomerStore | None = None,
) -> str:
    """Return the single letter natural analog peptide sequence."""
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()

    monomer_list = get_monomer_list(polymer_notation, PEPTIDE_POLYMER_TYPE, monomer_store)
    sb: list[str] = []
    for monomer_obj in monomer_list:
        single = monomer_obj.natural_analog
        if not single:
            single = "X"
        sb.append(single)
    return "".join(sb)


def get_modified_peptide_sequence(
    polymer_notation: str,
    delimiter: str | None = None,
    monomer_store: MonomerStore | None = None,
) -> str:
    """Return the modified peptide sequence with a specified delimiter."""
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()

    monomer_list = get_monomer_list(polymer_notation, PEPTIDE_POLYMER_TYPE, monomer_store)
    separator = delimiter if delimiter is not None else ""

    sb: list[str] = []
    for monomer_obj in monomer_list:
        if sb:
            sb.append(separator)
        sb.append(monomer_obj.alternate_id)
    return "".join(sb)


# ---------------------------------------------------------------------------
# Canonical notation
# ---------------------------------------------------------------------------

def get_simple_canonical_notation_for_rna(
    simple_notation: str,
    monomer_store: MonomerStore | None = None,
) -> str:
    """Return the canonical simple notation for an RNA polymer."""
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()

    validate_simple_notation(simple_notation, NUCLEIC_ACID_POLYMER_TYPE, monomer_store)

    # check for cycle info
    n_list = get_nucleotide_list(simple_notation, False, monomer_store)
    first_nuc = n_list[0]
    last_nuc = n_list[-1]
    standard_notation = simple_notation

    if first_nuc.get_sugar_monomer(monomer_store) is not None:
        if last_nuc.get_phosphate_monomer(monomer_store) is None:
            raise NotationException(
                "Non-cyclizable nucleic acid: last nucleotide contains no phosphate linker"
            )
    else:
        if last_nuc.get_phosphate_monomer(monomer_store) is not None:
            raise NotationException(
                "Non-cyclizable nucleic acid: duplicate phosphate linkers at both ends"
            )
        else:
            first_delimiter = simple_notation.index(GROUP_LEVEL_DELIMITER)
            back_part = simple_notation[first_delimiter + 1 :]
            front_part = simple_notation[:first_delimiter]
            standard_notation = back_part + front_part

    return _get_simple_canonical_notation(standard_notation)


def get_simple_canonical_notation_for_peptide(
    simple_notation: str,
    monomer_store: MonomerStore | None = None,
) -> str:
    """Return the canonical simple notation for a peptide polymer."""
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()
    validate_simple_notation(simple_notation, PEPTIDE_POLYMER_TYPE, monomer_store)
    return _get_simple_canonical_notation(simple_notation)


def get_simple_canonical_notation_for_chem(
    simple_notation: str,
    monomer_store: MonomerStore | None = None,
) -> str:
    """Return the canonical simple notation for a chemical polymer."""
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()
    validate_simple_notation(simple_notation, CHEMICAL_POLYMER_TYPE, monomer_store)
    return _get_simple_canonical_notation(simple_notation)


def get_simple_canonical_notation(
    simple_notation: str,
    polymer_type: str,
    monomer_store: MonomerStore | None = None,
) -> str:
    """Return the canonical simple notation for any polymer type."""
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()

    if polymer_type == NUCLEIC_ACID_POLYMER_TYPE:
        return get_simple_canonical_notation_for_rna(simple_notation, monomer_store)
    elif polymer_type == PEPTIDE_POLYMER_TYPE:
        return get_simple_canonical_notation_for_peptide(simple_notation, monomer_store)
    elif polymer_type == CHEMICAL_POLYMER_TYPE:
        return get_simple_canonical_notation_for_chem(simple_notation, monomer_store)
    raise NotationException(f"Unknown Polymer Type: {polymer_type}")


def _get_simple_canonical_notation(simple_notation: str) -> str:
    """Compute the canonical (lexicographically smallest rotation) of a simple notation."""
    rotation_list: list[str] = [simple_notation]
    start = 0
    current = simple_notation.find(GROUP_LEVEL_DELIMITER, start)

    while current > 0:
        back_part = simple_notation[current + 1 :]
        front_part = simple_notation[:current]
        new_notation = back_part + GROUP_LEVEL_DELIMITER + front_part
        rotation_list.append(new_notation)
        start = current + 1
        current = simple_notation.find(GROUP_LEVEL_DELIMITER, start)

    rotation_list.sort()
    return rotation_list[0]


def get_simple_canonical_notation_map_entry(
    simple_notation: str,
    polymer_type: str,
    monomer_store: MonomerStore | None = None,
) -> tuple[int, str]:
    """Return (offset, canonical_notation) for the canonical rotation.

    The offset indicates how many monomer positions the notation was rotated.
    """
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()

    entry_list: list[str] = []

    if polymer_type == NUCLEIC_ACID_POLYMER_TYPE:
        nuc_list = get_nucleotide_list(simple_notation, False, monomer_store)
        first_nuc = nuc_list[0]
        last_nuc = nuc_list[-1]
        offset = 0

        if first_nuc.get_sugar_monomer(monomer_store) is not None:
            if last_nuc.get_phosphate_monomer(monomer_store) is None:
                raise NotationException(
                    "Non-cyclizable nucleic acid: last nucleotide contains no phosphate linker"
                )
        else:
            if last_nuc.get_phosphate_monomer(monomer_store) is not None:
                raise NotationException(
                    "Non-cyclizable nucleic acid: duplicate phosphate linkers at both ends"
                )
            else:
                first_nuc_notation = first_nuc.notation
                last_nuc_notation = last_nuc.notation
                new_notation = last_nuc_notation + first_nuc_notation
                last_nuc.notation = new_notation
                nuc_list.pop(0)
                offset = 1

        for _ in range(len(nuc_list)):
            notation = get_simple_notation_from_nucleotide_list(nuc_list)
            entry_list.append(f"{notation} {offset}")
            nuc = nuc_list[0]
            if nuc.get_sugar_monomer(monomer_store) is not None:
                offset += 1
            if nuc.get_base_monomer(monomer_store) is not None:
                offset += 1
            if nuc.get_phosphate_monomer(monomer_store) is not None:
                offset += 1
            nuc_list.pop(0)
            nuc_list.append(nuc)

    elif polymer_type in (PEPTIDE_POLYMER_TYPE, CHEMICAL_POLYMER_TYPE):
        monomer_id_list = get_monomer_id_list(simple_notation, polymer_type, monomer_store)
        tmp_list = list(monomer_id_list)

        for i in range(len(monomer_id_list)):
            monomer_id = monomer_id_list[i]
            notation = _get_simple_notation(tmp_list, polymer_type, monomer_store)
            entry_list.append(f"{notation} {i}")
            tmp_list.pop(0)
            tmp_list.append(monomer_id)

    if entry_list:
        entry_list.sort()
        entry_string = entry_list[0]
        entries = entry_string.split()
        can_notation = entries[0]
        position = int(entries[1])
        return (position, can_notation)

    raise NotationException(f"Unknown Polymer Type: {polymer_type}")


# ---------------------------------------------------------------------------
# MoleculeInfo calculation
# ---------------------------------------------------------------------------

def get_molecule_info_for_polymer(
    notation: str,
    polymer_type: str,
    monomer_store: MonomerStore | None = None,
) -> MoleculeInfo:
    """Return the MoleculeInfo of a simple polymer using a divide and conquer approach.

    All R groups are capped.
    """
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()

    if polymer_type == CHEMICAL_POLYMER_TYPE:
        notation = process_node(notation, CHEMICAL_POLYMER_TYPE, "X", monomer_store)

    chunks: list[str] = []
    tmp_notation_parts: list[str] = []

    group_iterator = SimpleNotationGroupIterator(notation)
    groups: list[str] = []
    while group_iterator.has_next_group():
        groups.append(group_iterator.next_group())

    for i in range(1, len(groups) + 1):
        if tmp_notation_parts:
            tmp_notation_parts.append(GROUP_LEVEL_DELIMITER)
        tmp_notation_parts.append(groups[i - 1])

        if (i % MONOMER_GROUP_COUNT_INTERVAL) == 0:
            chunks.append("".join(tmp_notation_parts))
            tmp_notation_parts = []

    if tmp_notation_parts:
        chunks.append("".join(tmp_notation_parts))

    cap_mi_list: list[MoleculeInfo] = []

    # R2 caps
    for i in range(len(chunks) - 1):
        tmp_notation = chunks[i]
        monomer_id_list = get_monomer_id_list(tmp_notation, polymer_type, monomer_store)
        monomer_obj = get_monomer(monomer_id_list[-1], polymer_type, monomer_store)
        cap_mi = monomer_obj.get_cap_molecule_info("R2")
        if cap_mi is not None:
            cap_mi_list.append(cap_mi)

    # R1 caps
    for i in range(1, len(chunks)):
        tmp_notation = chunks[i]
        monomer_id_list = get_monomer_id_list(tmp_notation, polymer_type, monomer_store)
        monomer_obj = get_monomer(monomer_id_list[0], polymer_type, monomer_store)
        cap_mi = monomer_obj.get_cap_molecule_info("R1")
        if cap_mi is not None:
            cap_mi_list.append(cap_mi)

    if polymer_type == CHEMICAL_POLYMER_TYPE:
        tmp_notation = chunks[0]
        monomer_id_list = get_monomer_id_list(tmp_notation, polymer_type, monomer_store)
        monomer_obj = get_monomer(monomer_id_list[0], polymer_type, monomer_store)
        attachments = monomer_obj.attachment_list
        for att in attachments:
            label = att.label
            if label in ("R1", "R2"):
                continue
            cap_mi = monomer_obj.get_cap_molecule_info(label)
            if cap_mi is not None:
                cap_mi_list.append(cap_mi)

    chunk_mi_list: list[MoleculeInfo] = []
    for chunk in chunks:
        complex_notation = get_complex_notation(chunk, polymer_type, monomer_store)
        # Import here to avoid circular imports
        from helm_notation.tools.complex_notation_parser import (
            get_complex_polymer_smiles,
        )

        smiles = get_complex_polymer_smiles(complex_notation, monomer_store)
        tmp_mi = get_molecule_info(smiles)
        chunk_mi_list.append(tmp_mi)

    return process_molecule_info(chunk_mi_list, cap_mi_list)


# ---------------------------------------------------------------------------
# Ad-hoc monomer ID generation
# ---------------------------------------------------------------------------

def get_ad_hoc_monomer_id_prefix(polymer_type: str) -> str:
    """Return the prefix for ad-hoc monomer IDs."""
    if polymer_type == CHEMICAL_POLYMER_TYPE:
        return "CM#"
    elif polymer_type == PEPTIDE_POLYMER_TYPE:
        return "PM#"
    elif polymer_type == NUCLEIC_ACID_POLYMER_TYPE:
        return "NM#"
    else:
        return "AM#"


def generate_next_ad_hoc_monomer_id(
    polymer_type: str,
    store: MonomerStore,
) -> str:
    """Generate the next unique ad-hoc monomer ID."""
    global _seed_map

    from helm_notation.monomer_factory import get_instance as get_monomer_factory

    internal_monomers: dict[str, Monomer] | None = None
    try:
        internal_monomers = get_monomer_factory().monomer_db.get(polymer_type)
    except Exception:
        pass

    monomers = store.get_monomers(polymer_type)

    seed = _seed_map.get(polymer_type, 0)
    seed += 1
    _seed_map[polymer_type] = seed

    result = f"{get_ad_hoc_monomer_id_prefix(polymer_type)}{seed}"

    if (monomers is not None and result in monomers) or (
        internal_monomers is not None and result in internal_monomers
    ):
        return generate_next_ad_hoc_monomer_id(polymer_type, store)
    else:
        return result


def reset_seed() -> None:
    """Reset the ad-hoc monomer ID seed map."""
    global _seed_map
    _seed_map = {}


# ---------------------------------------------------------------------------
# Notation replacement: SMILES -> temporary IDs
# ---------------------------------------------------------------------------

def get_notation_by_replacing_smiles(
    simple_notation: str,
    polymer_type: str,
    monomer_store: MonomerStore | None = None,
) -> str:
    """Replace SMILES in simple notation with temporary IDs."""
    if monomer_store is None:
        monomer_store = _get_default_monomer_store()

    monomer_ids = get_monomer_id_list(simple_notation, polymer_type, monomer_store)
    return _get_simple_notation(monomer_ids, polymer_type, monomer_store)
