"""Parser for the XHELM XML format.

Extracts complex notation strings and monomer stores from XHELM documents.
"""

from __future__ import annotations

from lxml import etree

from helm_notation.exceptions import MonomerException
from helm_notation.monomer_store import MonomerStore


def parse_xhelm_string(xhelm_string: str) -> etree._Element:
    """Parse an XHELM XML string and return the root element.

    Args:
        xhelm_string: The XHELM XML string.

    Returns:
        The root XML element.
    """
    return etree.fromstring(xhelm_string.encode("utf-8"))


def get_complex_notation_string(root_element: etree._Element) -> str:
    """Extract the complex notation string from the root node of an XHELM document.

    Args:
        root_element: The root XML element of the XHELM document.

    Returns:
        The complex notation string.
    """
    helm_notation_element = root_element.find("HelmNotation")
    if helm_notation_element is None:
        # Try with namespace
        for child in root_element:
            tag = etree.QName(child).localname if isinstance(child.tag, str) else ""
            if tag == "HelmNotation":
                return child.text or ""
        return ""
    return helm_notation_element.text or ""


def get_monomer_store(root_element: etree._Element) -> MonomerStore:
    """Generate a monomer store from a given XHELM document.

    Args:
        root_element: The root XML element of the XHELM document.

    Returns:
        A MonomerStore populated with monomers from the document.
    """
    from helm_notation.tools.monomer_parser import get_monomer_from_element

    monomer_store = MonomerStore()
    monomer_list_element = root_element.find("Monomers")

    if monomer_list_element is None:
        # Try with namespace
        for child in root_element:
            tag = etree.QName(child).localname if isinstance(child.tag, str) else ""
            if tag == "Monomers":
                monomer_list_element = child
                break

    if monomer_list_element is not None:
        for monomer_element in monomer_list_element:
            tag = (
                etree.QName(monomer_element).localname
                if isinstance(monomer_element.tag, str)
                else ""
            )
            if tag == "Monomer" or monomer_element.tag == "Monomer":
                m = get_monomer_from_element(monomer_element)
                monomer_store.add_monomer(m)

    return monomer_store
