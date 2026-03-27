"""Factory for loading and managing the monomer database (replaces singleton pattern)."""

from __future__ import annotations

import copy
import importlib.resources
import logging
from pathlib import Path
from typing import TYPE_CHECKING

from lxml import etree

from helm_notation.constants import NOTATION_DIRECTORY
from helm_notation.exceptions import MonomerException
from helm_notation.model.attachment import Attachment
from helm_notation.model.monomer import Monomer
from helm_notation.model.monomer_cache import MonomerCache
from helm_notation.monomer_store import MonomerStore

if TYPE_CHECKING:
    pass

logger = logging.getLogger(__name__)

MONOMER_DB_FILE_NAME = "MonomerDBGZEncoded.xml"
MONOMER_DB_FILE_PATH = NOTATION_DIRECTORY / MONOMER_DB_FILE_NAME

POLYMER_LIST_ELEMENT = "PolymerList"
POLYMER_ELEMENT = "Polymer"
POLYMER_TYPE_ATTRIBUTE = "polymerType"
ATTACHMENT_LIST_ELEMENT = "AttachmentList"


class MonomerFactory:
    """Factory that loads and manages the monomer database."""

    def __init__(self) -> None:
        self.monomer_db: dict[str, dict[str, Monomer]] = {}
        self.attachment_db: dict[str, Attachment] = {}
        self.smiles_monomer_db: dict[str, Monomer] = {}
        self._monomer_store: MonomerStore | None = None
        self._db_changed: bool = True

    @property
    def monomer_store(self) -> MonomerStore:
        if self._monomer_store is None:
            self._monomer_store = MonomerStore(
                self.monomer_db, self.smiles_monomer_db
            )
        return self._monomer_store

    def get_monomer_db(self, include_new: bool = True) -> dict[str, dict[str, Monomer]]:
        if include_new:
            return self.monomer_db
        return {
            pt: {
                mid: m
                for mid, m in mmap.items()
                if not m.new_monomer
            }
            for pt, mmap in self.monomer_db.items()
        }

    def get_smiles_monomer_db(self, include_new: bool = True) -> dict[str, Monomer]:
        if include_new:
            return self.smiles_monomer_db
        return {k: v for k, v in self.smiles_monomer_db.items() if not v.new_monomer}

    def get_attachment_db(self) -> dict[str, Attachment]:
        """Build and return a map of all unique attachments across all monomers."""
        result: dict[str, Attachment] = {}
        for _pt, monomers in self.monomer_db.items():
            for _mid, monomer in monomers.items():
                for att in monomer.attachment_list:
                    label = att.get_label()
                    if label and label not in result:
                        result[label] = att
        return result

    def get_polymer_types(self) -> list[str]:
        return sorted(self.monomer_db.keys())

    def add_new_monomer(self, monomer: Monomer) -> None:
        monomer.new_monomer = True
        self._add_monomer(monomer)
        self._db_changed = True

    def _add_monomer(self, monomer: Monomer) -> None:
        polymer_type = monomer.polymer_type
        if polymer_type not in self.monomer_db:
            self.monomer_db[polymer_type] = {}

        monomer_map = self.monomer_db[polymer_type]
        if monomer.alternate_id not in monomer_map:
            monomer_map[monomer.alternate_id] = copy.deepcopy(monomer)

        if monomer.can_smiles and monomer.can_smiles not in self.smiles_monomer_db:
            self.smiles_monomer_db[monomer.can_smiles] = monomer

        self._db_changed = True


# Module-level singleton
_instance: MonomerFactory | None = None
_db_changed: bool = True


def get_instance() -> MonomerFactory:
    """Get or create the singleton MonomerFactory instance."""
    global _instance
    if _instance is None:
        _instance = MonomerFactory()
        _initialize_monomer_cache(_instance)
    return _instance


def set_db_changed(changed: bool) -> None:
    global _db_changed
    _db_changed = changed


def finalize() -> None:
    """Reset the factory singleton."""
    global _instance, _db_changed
    _instance = None
    _db_changed = True


def _initialize_monomer_cache(factory: MonomerFactory) -> None:
    """Load monomer database from XML resource file."""
    from helm_notation.tools.monomer_parser import (
        get_attachment_from_element,
        get_monomer_from_element,
        validate_monomer,
    )
    from helm_notation.tools.structure_parser import get_unique_extended_smiles

    # Try local file first, then fall back to package resource
    xml_data: bytes | None = None

    local_path = MONOMER_DB_FILE_PATH
    if local_path.exists():
        try:
            xml_data = local_path.read_bytes()
            logger.info(f"{local_path} used for monomer cache initialization")
        except Exception:
            logger.info(f"Unable to use local monomer DB file: {MONOMER_DB_FILE_NAME}")

    if xml_data is None:
        resource_ref = importlib.resources.files("helm_notation.resources").joinpath(
            MONOMER_DB_FILE_NAME
        )
        xml_data = resource_ref.read_bytes()
        logger.info(f"Package resource used for monomer cache initialization")

    # Parse XML
    root = etree.fromstring(xml_data)
    ns = root.nsmap.get(None, "")

    # Build monomer DB
    if ns:
        polymer_list = root.find(f"{{{ns}}}{POLYMER_LIST_ELEMENT}")
    else:
        polymer_list = root.find(POLYMER_LIST_ELEMENT)

    if polymer_list is not None:
        if ns:
            polymers = polymer_list.findall(f"{{{ns}}}{POLYMER_ELEMENT}")
        else:
            polymers = polymer_list.findall(POLYMER_ELEMENT)

        for polymer_elem in polymers:
            polymer_type = polymer_elem.get(POLYMER_TYPE_ATTRIBUTE, "")
            monomer_map: dict[str, Monomer] = {}

            for monomer_elem in polymer_elem:
                try:
                    m = get_monomer_from_element(monomer_elem)
                    if validate_monomer(m):
                        monomer_map[m.alternate_id] = m
                except Exception as e:
                    logger.warning(f"Skipping invalid monomer: {e}")

            factory.monomer_db[polymer_type] = monomer_map

    # Build attachment DB
    if ns:
        attachment_list = root.find(f"{{{ns}}}{ATTACHMENT_LIST_ELEMENT}")
    else:
        attachment_list = root.find(ATTACHMENT_LIST_ELEMENT)

    if attachment_list is not None:
        for att_elem in attachment_list:
            try:
                att = get_attachment_from_element(att_elem)
                factory.attachment_db[att.alternate_id] = att
            except Exception as e:
                logger.warning(f"Skipping invalid attachment: {e}")

    # Build SMILES monomer DB
    for polymer_type, monomer_map in factory.monomer_db.items():
        for monomer_id, monomer in monomer_map.items():
            smiles = monomer.can_smiles
            try:
                smiles = get_unique_extended_smiles(smiles)
            except Exception:
                pass
            monomer.can_smiles = smiles
            factory.smiles_monomer_db[smiles] = monomer

    factory._db_changed = True
