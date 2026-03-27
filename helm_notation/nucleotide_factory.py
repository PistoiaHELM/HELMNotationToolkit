"""Factory for loading nucleotide templates from XML."""

from __future__ import annotations

import importlib.resources
import logging
from pathlib import Path

from lxml import etree

from helm_notation.constants import NOTATION_DIRECTORY, NOTATION_SOURCE
from helm_notation.exceptions import NotationException

logger = logging.getLogger(__name__)

NUCLEOTIDE_TEMPLATE_FILE_NAME = "NucleotideTemplates.xml"
NUCLEOTIDE_TEMPLATE_FILE_PATH = NOTATION_DIRECTORY / NUCLEOTIDE_TEMPLATE_FILE_NAME


class NucleotideFactory:
    """Factory for nucleotide templates."""

    def __init__(self) -> None:
        self.nucleotide_templates: dict[str, dict[str, str]] = {}
        self.reverse_nucleotide_map: dict[str, str] = {}

    def get_nucleotide_templates(self) -> dict[str, dict[str, str]]:
        return self.nucleotide_templates

    def set_nucleotide_templates(
        self, templates: dict[str, dict[str, str]]
    ) -> None:
        self.nucleotide_templates = templates
        self.reverse_nucleotide_map = _build_reverse_map(
            templates, NOTATION_SOURCE
        )

    def get_reverse_nucleotide_template_map(self) -> dict[str, str]:
        return self.reverse_nucleotide_map


# Module-level singleton
_instance: NucleotideFactory | None = None


def get_instance() -> NucleotideFactory:
    """Get or create the singleton NucleotideFactory instance."""
    global _instance
    if _instance is None:
        _instance = NucleotideFactory()
        _initialize_nucleotide_templates(_instance)
    return _instance


def finalize() -> None:
    """Reset the factory singleton."""
    global _instance
    _instance = None


def _build_reverse_map(
    templates: dict[str, dict[str, str]], notation_source: str
) -> dict[str, str]:
    """Build reverse map: notation -> symbol."""
    result: dict[str, str] = {}
    normal_map = templates.get(notation_source)
    if normal_map:
        for symbol, notation in normal_map.items():
            result[notation] = symbol
    return result


def _initialize_nucleotide_templates(factory: NucleotideFactory) -> None:
    """Load nucleotide templates from XML."""
    xml_data: bytes | None = None

    local_path = NUCLEOTIDE_TEMPLATE_FILE_PATH
    if local_path.exists():
        try:
            xml_data = local_path.read_bytes()
            logger.info(f"{local_path} used for nucleotide templates initialization")
        except Exception:
            logger.info("Unable to use local nucleotide templates")

    if xml_data is None:
        resource_ref = importlib.resources.files("helm_notation.resources").joinpath(
            NUCLEOTIDE_TEMPLATE_FILE_NAME
        )
        xml_data = resource_ref.read_bytes()
        logger.info("Package resource used for nucleotide templates initialization")

    templates = _parse_nucleotide_templates_xml(xml_data)
    factory.nucleotide_templates = templates
    factory.reverse_nucleotide_map = _build_reverse_map(templates, NOTATION_SOURCE)


def _parse_nucleotide_templates_xml(
    xml_data: bytes,
) -> dict[str, dict[str, str]]:
    """Parse NucleotideTemplates.xml into a nested dict.

    Structure: {notation_source: {symbol: notation}}
    """
    root = etree.fromstring(xml_data)
    ns = root.nsmap.get(None, "")
    templates: dict[str, dict[str, str]] = {}

    # Find all TEMPLATE elements (uppercase in actual XML)
    if ns:
        template_list = root.findall(f".//{{{ns}}}TEMPLATE")
    else:
        template_list = root.findall(".//TEMPLATE")

    for template_elem in template_list:
        source = template_elem.get("notationSource", "")
        if source not in templates:
            templates[source] = {}

        if ns:
            nucleotides = template_elem.findall(f"{{{ns}}}NUCLEOTIDE")
        else:
            nucleotides = template_elem.findall("NUCLEOTIDE")

        for nuc_elem in nucleotides:
            if ns:
                symbol_elem = nuc_elem.find(f"{{{ns}}}SYMBOL")
                notation_elem = nuc_elem.find(f"{{{ns}}}MONOMER_NOTATION")
            else:
                symbol_elem = nuc_elem.find("SYMBOL")
                notation_elem = nuc_elem.find("MONOMER_NOTATION")

            if symbol_elem is not None and notation_elem is not None:
                symbol = symbol_elem.text or ""
                notation = notation_elem.text or ""
                templates[source][symbol] = notation

    return templates
