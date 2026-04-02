"""Port of org.helm.notation.MonomerFactory.

Singleton that loads the monomer database from MonomerDBGZEncoded.xml and
exposes ``monomer_db``, ``smiles_monomer_db``, ``attachment_db``, and a
``MonomerStore`` view.
"""

from __future__ import annotations

import logging
import xml.etree.ElementTree as ET
from pathlib import Path
from threading import Lock
from typing import Dict, Optional

from helm_notation import constants
from helm_notation.exceptions import MonomerException
from helm_notation.model import Attachment, Monomer, MonomerCache
from helm_notation.monomer_store import MonomerStore
from helm_notation.tools.monomer_parser import parse_monomer

logger = logging.getLogger(__name__)

_RESOURCE_XML = "MonomerDBGZEncoded.xml"
_DEFAULT_NS = "lmr"
_NS_PREFIX = "{lmr}"

_lock = Lock()
_instance: Optional["MonomerFactory"] = None


class MonomerFactory:
    """Singleton factory that loads & caches the monomer database."""

    def __init__(self) -> None:
        self.monomer_db: Dict[str, Dict[str, Monomer]] = {}
        self.smiles_monomer_db: Dict[str, Monomer] = {}
        self.attachment_db: Dict[str, Attachment] = {}
        self._monomer_store: Optional[MonomerStore] = None

    # ------------------------------------------------------------------
    # Public accessors
    # ------------------------------------------------------------------

    def get_monomer_db(self) -> Dict[str, Dict[str, Monomer]]:
        return self.monomer_db

    def get_smiles_monomer_db(self) -> Dict[str, Monomer]:
        return self.smiles_monomer_db

    def get_attachment_db(self) -> Dict[str, Attachment]:
        return self.attachment_db

    def get_monomer_store(self) -> MonomerStore:
        if self._monomer_store is None:
            self._monomer_store = MonomerStore(self.monomer_db, self.smiles_monomer_db)
        return self._monomer_store

    # ------------------------------------------------------------------
    # Singleton access
    # ------------------------------------------------------------------

    @staticmethod
    def get_instance() -> "MonomerFactory":
        global _instance
        with _lock:
            if _instance is None:
                factory = MonomerFactory()
                factory._initialize()
                _instance = factory
            return _instance

    @staticmethod
    def finalize_monomer_cache() -> None:
        global _instance
        with _lock:
            _instance = None

    # ------------------------------------------------------------------
    # Initialisation
    # ------------------------------------------------------------------

    def _initialize(self) -> None:
        """Load monomer DB from the bundled XML resource."""
        from helm_notation.resources import resource_path

        xml_path = resource_path(_RESOURCE_XML)
        tree = ET.parse(str(xml_path))
        root = tree.getroot()

        # Iterate <PolymerList> → <Polymer polymerType="..."> → <Monomer>
        for poly_list in root.iter(f"{_NS_PREFIX}PolymerList"):
            for polymer_el in poly_list.findall(f"{_NS_PREFIX}Polymer"):
                ptype = polymer_el.get("polymerType", "")
                if ptype not in self.monomer_db:
                    self.monomer_db[ptype] = {}
                for monomer_el in polymer_el.findall(f"{_NS_PREFIX}Monomer"):
                    try:
                        m = parse_monomer(monomer_el, ns=_NS_PREFIX)
                        self.monomer_db[ptype][m.alternate_id] = m
                        _register_smiles(self.smiles_monomer_db, m)
                        _register_attachments(self.attachment_db, m)
                    except Exception as exc:
                        logger.warning("Skipping monomer: %s", exc)

        self._monomer_store = MonomerStore(self.monomer_db, self.smiles_monomer_db)
        logger.info(
            "MonomerFactory initialized: %d polymer types, %d SMILES entries",
            len(self.monomer_db),
            len(self.smiles_monomer_db),
        )


def _register_smiles(db: Dict[str, Monomer], monomer: Monomer) -> None:
    smi = monomer.can_smiles
    if smi:
        # Try to canonicalize – fall back silently
        try:
            from helm_notation.chem import ChemEngine

            smi = ChemEngine.get_unique_extended_smiles(smi)
        except Exception:
            pass
        monomer.can_smiles = smi
        if smi not in db:
            db[smi] = monomer


def _register_attachments(db: Dict[str, Attachment], monomer: Monomer) -> None:
    for att in monomer.attachment_list:
        if att.alternate_id and att.alternate_id not in db:
            db[att.alternate_id] = att
