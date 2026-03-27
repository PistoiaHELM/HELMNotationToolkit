"""Exporter for the XHELM XML format.

Writes HELM notation and associated monomers to XHELM XML format.
"""

from __future__ import annotations

import logging

from lxml import etree

from helm_notation.exceptions import MonomerException
from helm_notation.monomer_store import MonomerStore

logger = logging.getLogger(__name__)

XHELM_ELEMENT = "Xhelm"
MONOMER_LIST_ELEMENT = "Monomers"
MONOMER_ELEMENT = "Monomer"
HELM_NOTATION_ELEMENT = "HelmNotation"


def write_xhelm(helm_notation: str, store: MonomerStore) -> str | None:
    """Write the XHELM XML format for the given notation and monomer store.

    Args:
        helm_notation: The HELM complex notation string.
        store: The MonomerStore containing monomer definitions.

    Returns:
        The XHELM XML string, or None on error.
    """
    root = etree.Element(XHELM_ELEMENT)

    helm_element = etree.SubElement(root, HELM_NOTATION_ELEMENT)
    helm_element.text = helm_notation

    monomer_list_element = etree.SubElement(root, MONOMER_LIST_ELEMENT)

    try:
        from helm_notation.tools.complex_notation_parser import get_polymer_node_list
        from helm_notation.tools.simple_notation_parser import get_monomer_list
        from helm_notation.tools.monomer_parser import get_monomer_element

        # Get all simple polymers
        seen: set[str] = set()
        monomer_set: list = []
        simple_polymers = get_polymer_node_list(helm_notation)

        for node in simple_polymers:
            monomers = get_monomer_list(node.label, node.type, store)
            for subnode in monomers:
                if not subnode.is_ad_hoc_monomer:
                    key = f"{subnode.polymer_type}:{subnode.alternate_id}"
                    if key not in seen:
                        seen.add(key)
                        monomer_set.append(subnode)

        # Distinct monomers
        for distinct_monomer in monomer_set:
            monomer_element = get_monomer_element(distinct_monomer)
            monomer_list_element.append(monomer_element)

    except Exception:
        logger.exception("Error generating XHELM")
        return None

    return etree.tostring(
        root, pretty_print=True, xml_declaration=True, encoding="unicode"
    )
