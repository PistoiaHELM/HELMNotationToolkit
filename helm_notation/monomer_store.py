"""Monomer store for managing monomers separately from the global factory."""

from __future__ import annotations

import copy
from typing import TYPE_CHECKING

from helm_notation.exceptions import MonomerException
from helm_notation.model.monomer import Monomer

if TYPE_CHECKING:
    pass


class MonomerStore:
    """Store for monomers, separating file-specific monomers from the global DB."""

    def __init__(
        self,
        monomer_db: dict[str, dict[str, Monomer]] | None = None,
        smiles_monomer_db: dict[str, Monomer] | None = None,
    ):
        self.monomer_db: dict[str, dict[str, Monomer]] = monomer_db or {}
        self.smiles_monomer_db: dict[str, Monomer] = smiles_monomer_db or {}

    def add_monomer(self, monomer: Monomer, db_changed: bool = False) -> None:
        """Add a monomer to the store."""
        from helm_notation.tools.structure_parser import get_unique_extended_smiles

        polymer_type = monomer.polymer_type
        alternate_id = monomer.alternate_id
        smiles_string = monomer.can_smiles

        try:
            smiles_string = get_unique_extended_smiles(smiles_string)
        except Exception:
            smiles_string = monomer.can_smiles

        has_smiles = bool(smiles_string)

        if polymer_type not in self.monomer_db:
            self.monomer_db[polymer_type] = {}

        monomer_copy = copy.deepcopy(monomer)
        if has_smiles:
            monomer_copy.can_smiles = smiles_string

        if alternate_id not in self.monomer_db[polymer_type]:
            self.monomer_db[polymer_type][alternate_id] = monomer_copy

            if has_smiles and smiles_string not in self.smiles_monomer_db:
                self.smiles_monomer_db[smiles_string] = monomer_copy

        if db_changed:
            from helm_notation.monomer_factory import set_db_changed
            set_db_changed(True)

    def has_monomer(self, polymer_type: str, alternate_id: str) -> bool:
        monomer_map = self.monomer_db.get(polymer_type)
        return monomer_map is not None and alternate_id in monomer_map

    def get_monomer(self, polymer_type: str, alternate_id: str) -> Monomer | None:
        monomer_map = self.monomer_db.get(polymer_type)
        if monomer_map:
            return monomer_map.get(alternate_id)
        return None

    def get_monomer_by_smiles(self, smiles: str) -> Monomer | None:
        return self.smiles_monomer_db.get(smiles)

    def get_monomers(self, polymer_type: str) -> dict[str, Monomer] | None:
        return self.monomer_db.get(polymer_type)

    def add_new_monomer(self, monomer: Monomer) -> None:
        monomer.new_monomer = True
        self.add_monomer(monomer, db_changed=True)

    def is_empty(self) -> bool:
        return not self.monomer_db

    def clear(self) -> None:
        self.monomer_db.clear()
        self.smiles_monomer_db.clear()

    @property
    def polymer_type_set(self) -> set[str]:
        return set(self.monomer_db.keys())

    def get_all_monomers_list(self) -> list[Monomer]:
        result: list[Monomer] = []
        for polymer_type in self.polymer_type_set:
            monomer_map = self.get_monomers(polymer_type)
            if monomer_map:
                result.extend(monomer_map.values())
        return result

    def __str__(self) -> str:
        parts: list[str] = []
        for monomer_map in self.monomer_db.values():
            line_parts = [
                f"{m.alternate_id}({m.polymer_type})"
                for m in monomer_map.values()
            ]
            parts.append("; ".join(line_parts))
        return "\n".join(parts)
