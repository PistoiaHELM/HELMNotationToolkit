"""Parser for monomer XML definitions (replaces JDOM-based MonomerParser)."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

from lxml import etree

from helm_notation.exceptions import MonomerException, StructureException
from helm_notation.model.attachment import Attachment, PAIR_ATTACHMENT
from helm_notation.model.monomer import (
    BACKBONE_MONOMER_TYPE,
    BRANCH_MONOMER_TYPE,
    CHEMICAL_POLYMER_TYPE,
    NUCLEIC_ACID_POLYMER_TYPE,
    PEPTIDE_POLYMER_TYPE,
    SUPPORTED_POLYMER_TYPES,
    UNDEFINED_MONOMER_TYPE,
    Monomer,
)
from helm_notation.tools import molfile_encoder

if TYPE_CHECKING:
    pass

# XML element names
MONOMER_ELEMENT = "Monomer"
MONOMER_ID_ELEMENT = "MonomerID"
MONOMER_SMILES_ELEMENT = "MonomerSmiles"
MONOMER_MOL_FILE_ELEMENT = "MonomerMolFile"
MONOMER_TYPE_ELEMENT = "MonomerType"
POLYMER_TYPE_ELEMENT = "PolymerType"
NATURAL_ANALOG_ELEMENT = "NaturalAnalog"
MONOMER_NAME_ELEMENT = "MonomerName"
ATTACHMENTS_LIST_ELEMENT = "AttachmentList"
ATTACHMENTS_ELEMENT = "Attachments"
ATTACHMENT_ELEMENT = "Attachment"
ATTACHMENT_ID_ELEMENT = "AttachmentID"
ATTACHMENT_LABEL_ELEMENT = "AttachmentLabel"
CAP_GROUP_NAME_ELEMENT = "CapGroupName"
CAP_GROUP_SMILES_ELEMENT = "CapGroupSmiles"


def _get_child_text(element: etree._Element, tag: str) -> str | None:
    """Get text of a child element, handling namespaces."""
    ns = element.nsmap.get(None, "")
    if ns:
        child = element.find(f"{{{ns}}}{tag}")
    else:
        child = element.find(tag)
    if child is not None and child.text:
        return child.text
    return None


def get_attachment_from_element(element: etree._Element) -> Attachment:
    """Convert an Attachment XML element to an Attachment object."""
    att = Attachment()
    att.alternate_id = _get_child_text(element, ATTACHMENT_ID_ELEMENT) or ""
    att.label = _get_child_text(element, ATTACHMENT_LABEL_ELEMENT) or ""
    att.cap_group_name = _get_child_text(element, CAP_GROUP_NAME_ELEMENT) or ""
    att.cap_group_smiles = _get_child_text(element, CAP_GROUP_SMILES_ELEMENT) or ""
    return att


def get_attachment_element(att: Attachment) -> etree._Element:
    """Convert an Attachment object to an XML element."""
    element = etree.Element(ATTACHMENT_ELEMENT)

    if att.alternate_id:
        e = etree.SubElement(element, ATTACHMENT_ID_ELEMENT)
        e.text = att.alternate_id

    if att.label and att.label.strip():
        e = etree.SubElement(element, ATTACHMENT_LABEL_ELEMENT)
        e.text = att.label

    if att.cap_group_name:
        e = etree.SubElement(element, CAP_GROUP_NAME_ELEMENT)
        e.text = att.cap_group_name

    if att.cap_group_smiles:
        e = etree.SubElement(element, CAP_GROUP_SMILES_ELEMENT)
        e.text = att.cap_group_smiles

    return element


def get_monomer_from_element(element: etree._Element) -> Monomer:
    """Convert a Monomer XML element to a Monomer object."""
    m = Monomer()
    m.alternate_id = _get_child_text(element, MONOMER_ID_ELEMENT) or ""
    m.can_smiles = _get_child_text(element, MONOMER_SMILES_ELEMENT) or ""

    encoded_molfile = _get_child_text(element, MONOMER_MOL_FILE_ELEMENT)
    if encoded_molfile:
        try:
            m.molfile = molfile_encoder.decode(encoded_molfile) or ""
        except Exception:
            raise MonomerException("Invalid monomer molfile")
    else:
        m.molfile = ""

    m.monomer_type = _get_child_text(element, MONOMER_TYPE_ELEMENT) or ""
    m.polymer_type = _get_child_text(element, POLYMER_TYPE_ELEMENT) or ""
    m.natural_analog = _get_child_text(element, NATURAL_ANALOG_ELEMENT) or ""
    m.name = _get_child_text(element, MONOMER_NAME_ELEMENT) or ""

    # Parse attachments
    ns = element.nsmap.get(None, "")
    if ns:
        att_elem = element.find(f"{{{ns}}}{ATTACHMENTS_ELEMENT}")
    else:
        att_elem = element.find(ATTACHMENTS_ELEMENT)

    if att_elem is not None:
        if ns:
            att_children = att_elem.findall(f"{{{ns}}}{ATTACHMENT_ELEMENT}")
        else:
            att_children = att_elem.findall(ATTACHMENT_ELEMENT)
        m.attachment_list = [get_attachment_from_element(a) for a in att_children]

    return m


def get_monomer_element(monomer: Monomer) -> etree._Element:
    """Convert a Monomer object to an XML element."""
    element = etree.Element(MONOMER_ELEMENT)

    if monomer.alternate_id:
        e = etree.SubElement(element, MONOMER_ID_ELEMENT)
        e.text = monomer.alternate_id

    if monomer.can_smiles:
        e = etree.SubElement(element, MONOMER_SMILES_ELEMENT)
        e.text = monomer.can_smiles

    if monomer.molfile:
        e = etree.SubElement(element, MONOMER_MOL_FILE_ELEMENT)
        try:
            e.text = molfile_encoder.encode(monomer.molfile)
        except Exception:
            raise MonomerException("Invalid monomer molfile")

    if monomer.monomer_type:
        e = etree.SubElement(element, MONOMER_TYPE_ELEMENT)
        e.text = monomer.monomer_type

    if monomer.polymer_type:
        e = etree.SubElement(element, POLYMER_TYPE_ELEMENT)
        e.text = monomer.polymer_type

    if monomer.natural_analog:
        e = etree.SubElement(element, NATURAL_ANALOG_ELEMENT)
        e.text = monomer.natural_analog

    if monomer.name:
        e = etree.SubElement(element, MONOMER_NAME_ELEMENT)
        e.text = monomer.name

    if monomer.attachment_list:
        attachments = etree.SubElement(element, ATTACHMENTS_ELEMENT)
        for att in monomer.attachment_list:
            attachments.append(get_attachment_element(att))

    return element


def validate_attachment(attachment: Attachment) -> bool:
    """Validate an Attachment object."""
    from helm_notation.tools.structure_parser import validate_smiles

    if not attachment.alternate_id:
        raise MonomerException("Attachment must have unique ID")

    smiles = attachment.cap_group_smiles
    if smiles:
        if not validate_smiles(smiles):
            raise MonomerException("Attachment cap group SMILES is invalid")

        labels = _get_attachment_labels(smiles)
        if not labels or len(labels) != 1:
            raise MonomerException("Attachment must have one R group in SMILES")

        if labels[0] != attachment.label:
            raise MonomerException(
                "R group in monomer SMILES and R group label must match"
            )

    return True


def validate_monomer(monomer: Monomer) -> bool:
    """Validate a Monomer object."""
    from helm_notation.tools.structure_parser import (
        get_molecule,
        find_rgroup_atom,
        validate_smiles,
    )

    if monomer is None:
        raise MonomerException("Monomer is null")

    polymer_type = monomer.polymer_type
    if not polymer_type:
        raise MonomerException("Monomer has no polymer type defined")
    if polymer_type not in SUPPORTED_POLYMER_TYPES:
        raise MonomerException(f"Unknown polymer type '{polymer_type}'")

    monomer_type = monomer.monomer_type
    if not monomer_type:
        raise MonomerException("Monomer has no monomer type defined")

    if polymer_type == CHEMICAL_POLYMER_TYPE:
        if monomer_type != UNDEFINED_MONOMER_TYPE:
            raise MonomerException(
                f"Valid monomer type for chemical structures can only be '{UNDEFINED_MONOMER_TYPE}'"
            )
    else:
        if monomer_type not in (BACKBONE_MONOMER_TYPE, BRANCH_MONOMER_TYPE):
            raise MonomerException(
                f"Valid monomer type for simple polymer can only be "
                f"'{BACKBONE_MONOMER_TYPE}' or '{BRANCH_MONOMER_TYPE}'"
            )

    alternate_id = monomer.alternate_id
    if not alternate_id:
        raise MonomerException("Monomer has no monomerID defined")

    smiles = monomer.can_smiles
    molfile = monomer.molfile
    attachments = monomer.attachment_list
    error_note = f"{alternate_id} ({polymer_type})"

    if polymer_type != CHEMICAL_POLYMER_TYPE:
        if not smiles or not molfile or not attachments:
            raise MonomerException(
                "Monomers for specific polymer type must have structure info"
            )

    if smiles:
        if not validate_smiles(smiles):
            raise MonomerException(f"Monomer SMILES must be valid: {error_note}")

        attachment_labels = _get_attachment_labels(smiles)
        if len(attachment_labels) != len(set(attachment_labels)):
            raise MonomerException(
                f"Attachment labels on monomer must be unique: {error_note}"
            )

        if len(attachment_labels) != len(attachments):
            raise MonomerException(
                f"Attachment label number on monomer must match attachment number: {error_note}"
            )

        for att in attachments:
            validate_attachment(att)

        for label in attachment_labels:
            found = any(a.alternate_id.startswith(label) for a in attachments)
            if not found:
                raise MonomerException(
                    f"Attachment label in SMILES is not found in attachment list: {error_note}"
                )

    if alternate_id and polymer_type != CHEMICAL_POLYMER_TYPE:
        natural_analog = monomer.natural_analog
        if not natural_analog:
            raise MonomerException(
                f"Modified monomer must have natural analog defined: {error_note}"
            )
        if len(natural_analog) != 1:
            raise MonomerException(
                f"Natural analog must be single letter: {error_note}"
            )

    if not monomer.attachment_list:
        raise MonomerException(
            f"Monomer must have at least one attachment: {error_note}"
        )

    # Validate R-group connectivity
    if smiles:
        mol = get_molecule(smiles)
        if mol is not None:
            attachment_labels = _get_attachment_labels(smiles)
            for label in attachment_labels:
                rgroup_id = int(label[1:])
                try:
                    r_idx = find_rgroup_atom(mol, rgroup_id)
                except Exception:
                    raise MonomerException(
                        f"Unable to find R{rgroup_id} in monomer: {error_note}"
                    )
                atom = mol.GetAtomWithIdx(r_idx)
                if atom.GetDegree() != 1:
                    raise MonomerException(
                        f"R group can only connect with one atom in monomer: {error_note}"
                    )
                bond = atom.GetBonds()[0]
                if bond.GetBondType().name != "SINGLE":
                    raise MonomerException(
                        f"R group can only connect with another atom via single bond in monomer: {error_note}"
                    )

    return True


def validate_attachment_label(label: str) -> None:
    """Validate attachment label format (R# or 'pair')."""
    if label.lower() == PAIR_ATTACHMENT:
        return

    if not label or label[0] != "R":
        raise MonomerException("Invalid Attachment Label format")

    for ch in label[1:]:
        if not ch.isdigit():
            raise MonomerException("Invalid Attachment Label format")


def get_monomer_list(monomer_xml_string: str) -> list[Monomer]:
    """Parse XML string containing multiple monomers."""
    if not monomer_xml_string:
        return []

    root = etree.fromstring(monomer_xml_string.encode("utf-8"))
    result = []
    for child in root:
        m = get_monomer_from_element(child)
        if validate_monomer(m):
            result.append(m)
    return result


def get_monomer_from_xml(monomer_xml_string: str) -> Monomer | None:
    """Parse a single monomer from XML string."""
    if not monomer_xml_string:
        return None
    root = etree.fromstring(monomer_xml_string.encode("utf-8"))
    return get_monomer_from_element(root)


def fill_attachment_info(att: Attachment) -> None:
    """Fill attachment info from the global attachment database."""
    from helm_notation.monomer_factory import get_instance

    factory = get_instance()
    attachment_map = factory.attachment_db
    stored_att = attachment_map.get(att.alternate_id)
    if stored_att:
        att.label = stored_att.label
        att.cap_group_smiles = stored_att.cap_group_smiles
        att.cap_group_name = stored_att.cap_group_name


def _get_attachment_labels(extended_smiles: str) -> list[str]:
    """Extract R-group labels (R1, R2, etc.) from a SMILES string."""
    labels: list[str] = []
    i = 0
    while i < len(extended_smiles):
        r_pos = extended_smiles.find("R", i)
        if r_pos < 0:
            break
        r_pos += 1
        digits = ""
        while r_pos < len(extended_smiles) and extended_smiles[r_pos].isdigit():
            digits += extended_smiles[r_pos]
            r_pos += 1
        if digits:
            labels.append(f"R{digits}")
        i = r_pos
    return labels
