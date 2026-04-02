from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

from helm_notation.monomer_store import MonomerStore
from helm_notation.tools.monomer_parser import parse_monomer


def _to_root(xhelm: str | Path | ET.Element) -> ET.Element:
    if isinstance(xhelm, ET.Element):
        return xhelm
    if isinstance(xhelm, Path):
        return ET.parse(str(xhelm)).getroot()
    text = str(xhelm).strip()
    if text.startswith("<"):
        return ET.fromstring(text)
    return ET.parse(text).getroot()


def get_complex_notation_string(xhelm: str | Path | ET.Element) -> str:
    root = _to_root(xhelm)
    helm_element = root.find("HelmNotation")
    if helm_element is None or helm_element.text is None:
        return ""
    return helm_element.text.strip()


def get_monomer_store(xhelm: str | Path | ET.Element) -> MonomerStore:
    root = _to_root(xhelm)
    store = MonomerStore()
    monomer_list_element = root.find("Monomers")
    if monomer_list_element is None:
        return store

    for monomer_element in monomer_list_element.findall("Monomer"):
        monomer = parse_monomer(monomer_element)
        store.add_monomer(monomer)
    return store
