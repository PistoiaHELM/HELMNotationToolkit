"""Port of org.helm.notation.MonomerStore."""

from __future__ import annotations

import copy
import logging
from typing import Dict, Optional

from helm_notation.model import Monomer

logger = logging.getLogger(__name__)


class MonomerStore:
    """Holds the monomer DB (by polymer-type + monomer-ID) and a SMILES→Monomer index."""

    def __init__(
        self,
        monomer_db: Optional[Dict[str, Dict[str, Monomer]]] = None,
        smiles_monomer_db: Optional[Dict[str, Monomer]] = None,
    ):
        self.monomer_db: Dict[str, Dict[str, Monomer]] = monomer_db or {}
        self.smiles_monomer_db: Dict[str, Monomer] = smiles_monomer_db or {}

    # ------------------------------------------------------------------
    def get_monomer_db(self) -> Dict[str, Dict[str, Monomer]]:
        return self.monomer_db

    def get_smiles_monomer_db(self) -> Dict[str, Monomer]:
        return self.smiles_monomer_db

    # ------------------------------------------------------------------
    def get_monomer(self, smiles: str) -> Optional[Monomer]:
        """Return monomer by SMILES key."""
        return self.smiles_monomer_db.get(smiles)

    def has_monomer(self, polymer_type: str, monomer_id: str) -> bool:
        return monomer_id in self.monomer_db.get(polymer_type, {})

    # ------------------------------------------------------------------
    def add_monomer(self, monomer: Monomer) -> None:
        """Add (or update) a monomer to both the ID-DB and the SMILES-DB."""
        pt = monomer.polymer_type
        if pt not in self.monomer_db:
            self.monomer_db[pt] = {}
        self.monomer_db[pt][monomer.alternate_id] = monomer

        smiles_str = monomer.can_smiles or ""
        if smiles_str:
            # try canonical normalisation, fall back to raw
            try:
                from helm_notation.chem import ChemEngine

                smiles_str = ChemEngine.get_unique_extended_smiles(smiles_str)
            except Exception:
                smiles_str = monomer.can_smiles

        if smiles_str and smiles_str not in self.smiles_monomer_db:
            cp = copy.copy(monomer)
            cp.can_smiles = smiles_str
            self.smiles_monomer_db[smiles_str] = cp

    # ------------------------------------------------------------------
    def clear(self) -> None:
        self.monomer_db.clear()
        self.smiles_monomer_db.clear()
