"""Data model for monomer attachment points."""

from dataclasses import dataclass, field


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
    """Attachment point on a monomer (R-group)."""

    id: int = 0
    alternate_id: str = ""
    label: str = ""
    cap_group_name: str = ""
    cap_group_smiles: str = ""

    def get_label(self) -> str:
        if not self.label:
            return " "
        return self.label

    def __str__(self) -> str:
        return f"{self.label}: {self.cap_group_name}"
