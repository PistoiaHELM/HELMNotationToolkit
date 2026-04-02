"""Port of org.helm.notation.NucleotideFactory.

Singleton that loads nucleotide templates from NucleotideTemplates.xml and
provides symbol↔notation mappings.
"""

from __future__ import annotations

import logging
import xml.etree.ElementTree as ET
from threading import Lock
from typing import Dict, Optional

from helm_notation import constants

logger = logging.getLogger(__name__)

_NS_PREFIX = "{lmr}"
_RESOURCE_XML = "NucleotideTemplates.xml"

_lock = Lock()
_instance: Optional["NucleotideFactory"] = None


class NucleotideFactory:
    """Singleton that provides nucleotide template lookups."""

    def __init__(self) -> None:
        # {notation_source: {symbol: monomer_notation}}
        self.nucleotide_templates: Dict[str, Dict[str, str]] = {}
        # reverse map {notation: symbol}
        self.reverse_nucleotide_map: Dict[str, str] = {}

    # ------------------------------------------------------------------
    # Public accessors
    # ------------------------------------------------------------------

    def get_nucleotide_templates(self) -> Dict[str, Dict[str, str]]:
        return self.nucleotide_templates

    def get_reverse_nucleotide_map(self) -> Dict[str, str]:
        return self.reverse_nucleotide_map

    def get_notation_for_symbol(self, symbol: str) -> Optional[str]:
        """Return the monomer notation for *symbol* using the default HELM source."""
        source_map = self.nucleotide_templates.get(constants.NOTATION_SOURCE, {})
        return source_map.get(symbol)

    def get_symbol_for_notation(self, notation: str) -> Optional[str]:
        return self.reverse_nucleotide_map.get(notation)

    # ------------------------------------------------------------------
    # Singleton access
    # ------------------------------------------------------------------

    @staticmethod
    def get_instance() -> "NucleotideFactory":
        global _instance
        with _lock:
            if _instance is None:
                factory = NucleotideFactory()
                factory._initialize()
                _instance = factory
            return _instance

    @staticmethod
    def reset() -> None:
        global _instance
        with _lock:
            _instance = None

    # ------------------------------------------------------------------
    # Initialisation
    # ------------------------------------------------------------------

    def _initialize(self) -> None:
        from helm_notation.resources import resource_path

        xml_path = resource_path(_RESOURCE_XML)
        tree = ET.parse(str(xml_path))
        root = tree.getroot()

        for template_el in root.findall(f"{_NS_PREFIX}TEMPLATE"):
            source = template_el.get("notationSource", constants.NOTATION_SOURCE)
            nuc_map: Dict[str, str] = {}
            for nuc_el in template_el.findall(f"{_NS_PREFIX}NUCLEOTIDE"):
                sym_el = nuc_el.find(f"{_NS_PREFIX}SYMBOL")
                not_el = nuc_el.find(f"{_NS_PREFIX}MONOMER_NOTATION")
                if sym_el is not None and not_el is not None:
                    sym = (sym_el.text or "").strip()
                    notation = (not_el.text or "").strip()
                    if sym and notation:
                        nuc_map[sym] = notation
            self.nucleotide_templates[source] = nuc_map

        # Build reverse map for default source
        helm_map = self.nucleotide_templates.get(constants.NOTATION_SOURCE, {})
        self.reverse_nucleotide_map = {v: k for k, v in helm_map.items()}

        logger.info(
            "NucleotideFactory initialized: %d templates, %d reverse entries",
            sum(len(v) for v in self.nucleotide_templates.values()),
            len(self.reverse_nucleotide_map),
        )
