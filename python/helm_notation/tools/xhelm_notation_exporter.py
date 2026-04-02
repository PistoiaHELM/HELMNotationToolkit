from __future__ import annotations

import xml.etree.ElementTree as ET

from helm_notation.monomer_store import MonomerStore
from helm_notation.tools.monomer_parser import monomer_to_element


def write_xhelm(helm_notation: str, store: MonomerStore | None) -> str:
    root = ET.Element("Xhelm")
    helm_element = ET.SubElement(root, "HelmNotation")
    helm_element.text = helm_notation

    monomer_list_element = ET.SubElement(root, "Monomers")

    if store is not None:
        seen: set[tuple[str, str]] = set()
        for polymer_type in sorted(store.get_monomer_db().keys()):
            monomer_map = store.get_monomer_db().get(polymer_type, {})
            for monomer_id in sorted(monomer_map.keys()):
                monomer = monomer_map[monomer_id]
                if monomer.ad_hoc_monomer:
                    continue
                key = (monomer.polymer_type, monomer.alternate_id)
                if key in seen:
                    continue
                seen.add(key)
                monomer_list_element.append(monomer_to_element(monomer))

    _indent(root)
    return ET.tostring(root, encoding="unicode", xml_declaration=False)


def _indent(element: ET.Element, level: int = 0) -> None:
    indent = "   "
    i = "\n" + level * indent
    if len(element):
        if not element.text or not element.text.strip():
            element.text = i + indent
        for child in element:
            _indent(child, level + 1)
        if not element[-1].tail or not element[-1].tail.strip():
            element[-1].tail = i
    if level and (not element.tail or not element.tail.strip()):
        element.tail = i
