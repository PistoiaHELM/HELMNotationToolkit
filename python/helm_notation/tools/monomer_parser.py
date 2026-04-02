"""Port of org.helm.notation.tools.MonomerParser.

Reads/writes Monomer XML (JDOM elements → Python dataclasses) and validates
attachment / monomer structures.
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional

from helm_notation.exceptions import MonomerException
from helm_notation.model import Attachment, Monomer
from helm_notation.tools import molfile_encoder

# XML element / attribute names (mirror Java constants)
MONOMER_ELEMENT = "Monomer"
MONOMER_ID_ELEMENT = "MonomerID"
MONOMER_SMILES_ELEMENT = "MonomerSmiles"
MONOMER_MOL_FILE_ELEMENT = "MonomerMolFile"
MONOMER_TYPE_ELEMENT = "MonomerType"
POLYMER_TYPE_ELEMENT = "PolymerType"
NATURAL_ANALOG_ELEMENT = "NaturalAnalog"
MONOMER_NAME_ELEMENT = "MonomerName"
ATTACHMENTS_ELEMENT = "Attachments"
ATTACHMENT_ELEMENT = "Attachment"
ATTACHMENT_ID_ELEMENT = "AttachmentID"
ATTACHMENT_LABEL_ELEMENT = "AttachmentLabel"
CAP_GROUP_NAME_ELEMENT = "CapGroupName"
CAP_GROUP_SMILES_ELEMENT = "CapGroupSmiles"

_R_LABEL_RE = re.compile(r"R(\d+)")


# ------------------------------------------------------------------
# Attachment helpers
# ------------------------------------------------------------------

def parse_attachment(element: ET.Element, ns: str = "") -> Attachment:
    """Convert an <Attachment> XML element to an Attachment dataclass."""
    def _text(tag: str) -> str:
        child = element.find(f"{ns}{tag}")
        return child.text.strip() if child is not None and child.text else ""

    return Attachment(
        alternate_id=_text(ATTACHMENT_ID_ELEMENT),
        label=_text(ATTACHMENT_LABEL_ELEMENT),
        cap_group_name=_text(CAP_GROUP_NAME_ELEMENT),
        cap_group_smiles=_text(CAP_GROUP_SMILES_ELEMENT),
    )


def validate_attachment(att: Attachment) -> bool:
    """Validate an Attachment according to the same rules as the Java code."""
    if not att.alternate_id:
        raise MonomerException("Attachment must have unique ID")
    smi = att.cap_group_smiles
    if smi:
        from helm_notation.chem import ChemEngine

        if not ChemEngine.validate_smiles(smi):
            raise MonomerException("Attachment cap group SMILES is invalid")
        labels = get_attachment_labels(smi)
        if len(labels) != 1:
            raise MonomerException("Attachment must have one R group in SMILES")
        if labels[0] != att.label:
            raise MonomerException(
                "R group in monomer SMILES and R group label must match"
            )
    return True


def get_attachment_labels(smiles: str) -> List[str]:
    """Extract R-group labels from a SMILES string (e.g. ['R1', 'R2'])."""
    labels: List[str] = []
    tokens = smiles.split("R")
    for i in range(1, len(tokens)):
        digits = ""
        for ch in tokens[i]:
            if ch.isdigit():
                digits += ch
            else:
                break
        if digits:
            labels.append(f"R{digits}")
    return labels


# ------------------------------------------------------------------
# Monomer helpers
# ------------------------------------------------------------------

def parse_monomer(element: ET.Element, ns: str = "") -> Monomer:
    """Convert a <Monomer> XML element to a Monomer dataclass."""
    def _text(tag: str) -> str:
        child = element.find(f"{ns}{tag}")
        return child.text.strip() if child is not None and child.text else ""

    m = Monomer(
        alternate_id=_text(MONOMER_ID_ELEMENT),
        can_smiles=_text(MONOMER_SMILES_ELEMENT),
        monomer_type=_text(MONOMER_TYPE_ELEMENT),
        polymer_type=_text(POLYMER_TYPE_ELEMENT),
        natural_analog=_text(NATURAL_ANALOG_ELEMENT),
        name=_text(MONOMER_NAME_ELEMENT),
    )

    encoded_molfile = _text(MONOMER_MOL_FILE_ELEMENT)
    if encoded_molfile:
        try:
            m.molfile = molfile_encoder.decode(encoded_molfile)
        except Exception:
            raise MonomerException("Invalid monomer molfile")

    att_root = element.find(f"{ns}{ATTACHMENTS_ELEMENT}")
    if att_root is not None:
        for att_el in att_root.findall(f"{ns}{ATTACHMENT_ELEMENT}"):
            m.attachment_list.append(parse_attachment(att_el, ns))

    return m


def monomer_to_element(monomer: Monomer) -> ET.Element:
    """Serialize a Monomer dataclass back to an XML Element."""
    root = ET.Element(MONOMER_ELEMENT)

    def _add(tag: str, text: str) -> None:
        if text:
            el = ET.SubElement(root, tag)
            el.text = text

    _add(MONOMER_ID_ELEMENT, monomer.alternate_id)
    _add(MONOMER_SMILES_ELEMENT, monomer.can_smiles)
    if monomer.molfile:
        try:
            encoded = molfile_encoder.encode(monomer.molfile)
            _add(MONOMER_MOL_FILE_ELEMENT, encoded)
        except Exception:
            raise MonomerException("Invalid monomer molfile")
    _add(MONOMER_TYPE_ELEMENT, monomer.monomer_type)
    _add(POLYMER_TYPE_ELEMENT, monomer.polymer_type)
    _add(NATURAL_ANALOG_ELEMENT, monomer.natural_analog)
    _add(MONOMER_NAME_ELEMENT, monomer.name)

    if monomer.attachment_list:
        atts = ET.SubElement(root, ATTACHMENTS_ELEMENT)
        for att in monomer.attachment_list:
            att_el = ET.SubElement(atts, ATTACHMENT_ELEMENT)
            if att.alternate_id:
                ET.SubElement(att_el, ATTACHMENT_ID_ELEMENT).text = att.alternate_id
            if att.label:
                ET.SubElement(att_el, ATTACHMENT_LABEL_ELEMENT).text = att.label
            if att.cap_group_name:
                ET.SubElement(att_el, CAP_GROUP_NAME_ELEMENT).text = att.cap_group_name
            if att.cap_group_smiles:
                ET.SubElement(att_el, CAP_GROUP_SMILES_ELEMENT).text = att.cap_group_smiles

    return root
