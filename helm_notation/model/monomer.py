"""Data model for monomers."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from helm_notation.model.attachment import Attachment
from helm_notation.model.molecule_info import MoleculeInfo

if TYPE_CHECKING:
    pass

# Polymer types
NUCLEIC_ACID_POLYMER_TYPE = "RNA"
PEPTIDE_POLYMER_TYPE = "PEPTIDE"
CHEMICAL_POLYMER_TYPE = "CHEM"
SUPPORTED_POLYMER_TYPES = [
    NUCLEIC_ACID_POLYMER_TYPE,
    PEPTIDE_POLYMER_TYPE,
    CHEMICAL_POLYMER_TYPE,
]

# Monomer types
BACKBONE_MONOMER_TYPE = "Backbone"
BRANCH_MONOMER_TYPE = "Branch"
UNDEFINED_MONOMER_TYPE = "Undefined"

STARTING_NAME = "5"
ATTACHMENT_LIST_DELIMITER = "$"

# Standard monomer IDs
ID_A = "A"
ID_G = "G"
ID_C = "C"
ID_U = "U"
ID_T = "T"
ID_R = "R"
ID_dR = "dR"
ID_P = "P"
ID_X = "X"
ID_ALA = "Ala"
ID_ARG = "Arg"
ID_ASN = "Asn"
ID_ASP = "Asp"
ID_CYS = "Cys"
ID_GLU = "Glu"
ID_GLN = "Gln"
ID_GLY = "Gly"
ID_HIS = "His"
ID_ILE = "Ile"
ID_LEU = "Leu"
ID_LYS = "Lys"
ID_MET = "Met"
ID_PHE = "Phe"
ID_PRO = "Pro"
ID_SER = "Ser"
ID_THR = "Thr"
ID_TRP = "Trp"
ID_TYR = "Tyr"
ID_VAL = "Val"
ID_CHEMICAL_STRUCTURE = "chemical structure"


@dataclass
class Monomer:
    """Data model for a monomer unit in HELM notation."""

    # Class-level constants (mirrors Java's static fields)
    NUCLIEC_ACID_POLYMER_TYPE = NUCLEIC_ACID_POLYMER_TYPE  # Java spelling
    PEPTIDE_POLYMER_TYPE = PEPTIDE_POLYMER_TYPE
    CHEMICAL_POLYMER_TYPE = CHEMICAL_POLYMER_TYPE
    BACKBONE_MOMONER_TYPE = BACKBONE_MONOMER_TYPE  # Java spelling
    BRANCH_MOMONER_TYPE = BRANCH_MONOMER_TYPE  # Java spelling
    UNDEFINED_MOMONER_TYPE = UNDEFINED_MONOMER_TYPE  # Java spelling

    id: int = 0
    alternate_id: str = ""
    natural_analog: str = ""
    name: str = ""
    can_smiles: str = ""
    molfile: str = ""
    monomer_type: str = ""
    polymer_type: str = ""
    attachment_list: list[Attachment] = field(default_factory=list)
    new_monomer: bool = False
    ad_hoc_monomer: bool = False

    def get_attachment(self, label: str) -> Attachment | None:
        """Get a specific attachment by label (case-insensitive)."""
        for att in self.attachment_list:
            if att.get_label().lower() == label.lower():
                return att
        return None

    def get_cap_molecule_info(self, label: str) -> MoleculeInfo | None:
        """Return MoleculeInfo for the cap group of the given R-group label."""
        from helm_notation.tools.structure_parser import get_molecule_info

        for att in self.attachment_list:
            if att.get_label().lower() == label.lower():
                cap_smi = att.cap_group_smiles
                return get_molecule_info(cap_smi)
        return None

    def add_attachment(self, attachment: Attachment) -> bool:
        """Add a new attachment if one with the same label doesn't already exist."""
        for a in self.attachment_list:
            if a.get_label().lower() == attachment.get_label().lower():
                return False
        self.attachment_list.append(attachment)
        return True

    @property
    def is_modified(self) -> bool:
        """Test if this monomer has been modified (alternate_id length > 1)."""
        return len(self.alternate_id) > 1

    def is_same_type(self, other: Monomer) -> bool:
        """Compare polymer type, monomer type, and natural analog."""
        analog1 = self.natural_analog or self.alternate_id
        analog2 = other.natural_analog or other.alternate_id
        return (
            self.monomer_type == other.monomer_type
            and self.polymer_type == other.polymer_type
            and analog1 == analog2
        )

    @property
    def attachment_list_string(self) -> str:
        """Return sorted, delimited string of attachment alternate IDs."""
        ids = sorted(a.alternate_id for a in self.attachment_list)
        return ATTACHMENT_LIST_DELIMITER.join(ids)

    def attachment_equals(self, other: Monomer) -> bool:
        return self.attachment_list_string == other.attachment_list_string

    def attachment_contains(self, other: Monomer) -> bool:
        return other.attachment_list_string in self.attachment_list_string

    def contain_any_atom(self) -> bool:
        """Check if this monomer's structure contains any 'A' (any) atom."""
        from helm_notation.tools.structure_parser import get_molecule

        smiles = self.can_smiles
        if smiles:
            mol = get_molecule(smiles)
            if mol is not None:
                for atom in mol.GetAtoms():
                    if atom.GetSymbol() == "A":
                        return True
        return False
