"""Data model for nucleotides."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from helm_notation.model.monomer import NUCLEIC_ACID_POLYMER_TYPE, Monomer

if TYPE_CHECKING:
    from helm_notation.monomer_store import MonomerStore

STARTING_POSITION_TYPE = 1
MIDDLE_POSITION_TYPE = 2
ENDING_POSITION_TYPE = 3


@dataclass
class Nucleotide:
    """Data model for a nucleotide unit."""

    STARTING_POSITION_TYPE = STARTING_POSITION_TYPE
    MIDDLE_POSITION_TYPE = MIDDLE_POSITION_TYPE
    ENDING_POSITION_TYPE = ENDING_POSITION_TYPE

    notation: str = ""
    position_type: int = MIDDLE_POSITION_TYPE
    symbol: str = ""

    @property
    def is_modified(self) -> bool:
        if self.notation:
            if "[" in self.notation:
                return True
            if "(X)" in self.notation:
                return True
            if self.notation.endswith(")"):
                return True
        return False

    def unmodified_without_phosphate(self) -> bool:
        if not self.notation:
            return True
        return "[" not in self.notation and "X" not in self.notation

    def get_natural_analog(self, monomer_store: MonomerStore | None = None) -> str:
        """Return the natural analog of this nucleotide."""
        if monomer_store is None:
            from helm_notation.monomer_factory import get_instance
            monomer_store = get_instance().monomer_store

        base_notation = self._get_base_symbol()
        if base_notation is None:
            return "X"

        # Strip brackets if present
        if (
            base_notation
            and base_notation.startswith("[")
            and base_notation.endswith("]")
        ):
            base_notation = base_notation[1:-1]

        try:
            monomers = monomer_store.get_monomers(NUCLEIC_ACID_POLYMER_TYPE)
            if monomers:
                m = monomers.get(base_notation)
                if m is None:
                    smiles_db = monomer_store.smiles_monomer_db
                    m = smiles_db.get(base_notation)
                if m:
                    return m.natural_analog
        except Exception:
            pass

        return "X"

    def get_phosphate_monomer(
        self, monomer_store: MonomerStore | None = None
    ) -> Monomer | None:
        if monomer_store is None:
            from helm_notation.monomer_factory import get_instance
            monomer_store = get_instance().monomer_store

        symbol = self._get_phosphate_symbol()
        if symbol:
            try:
                monomers = monomer_store.get_monomers(NUCLEIC_ACID_POLYMER_TYPE)
                if monomers:
                    return monomers.get(symbol)
            except Exception:
                pass
        return None

    def get_base_monomer(
        self, monomer_store: MonomerStore | None = None
    ) -> Monomer | None:
        if monomer_store is None:
            from helm_notation.monomer_factory import get_instance
            monomer_store = get_instance().monomer_store

        symbol = self._get_base_symbol()
        if symbol:
            # Strip brackets
            if symbol.startswith("[") and symbol.endswith("]"):
                symbol = symbol[1:-1]
            try:
                monomers = monomer_store.get_monomers(NUCLEIC_ACID_POLYMER_TYPE)
                if monomers:
                    return monomers.get(symbol)
            except Exception:
                pass
        return None

    def get_sugar_monomer(
        self, monomer_store: MonomerStore | None = None
    ) -> Monomer | None:
        if monomer_store is None:
            from helm_notation.monomer_factory import get_instance
            monomer_store = get_instance().monomer_store

        symbol = self._get_sugar_symbol()
        if symbol:
            try:
                monomers = monomer_store.get_monomers(NUCLEIC_ACID_POLYMER_TYPE)
                if monomers:
                    return monomers.get(symbol)
            except Exception:
                pass
        return None

    @property
    def linker_notation(self) -> str:
        p_symbol = self._get_phosphate_symbol()
        if not p_symbol:
            return ""
        if len(p_symbol) > 1:
            return f"[{p_symbol}]"
        return p_symbol

    @property
    def nucleoside_notation(self) -> str:
        linker_len = len(self.linker_notation)
        return self.notation[: len(self.notation) - linker_len]

    def _get_phosphate_symbol(self) -> str | None:
        nc = self.notation
        if not nc:
            return None

        if "(" in nc:
            start = nc.index(")")
            symbol = nc[start + 1 :]
            return re.sub(r"[\[\]]", "", symbol) or None

        if self.position_type == ENDING_POSITION_TYPE:
            if nc.startswith("["):
                idx = nc.index("]")
                symbol = nc[idx + 1 :]
                return re.sub(r"[\[\]]", "", symbol) or None
            else:
                return nc[1:] or None
        else:
            if nc.endswith("]"):
                # Find last '['
                for i in range(len(nc) - 1, -1, -1):
                    if nc[i] == "[":
                        symbol = nc[i:]
                        return re.sub(r"[\[\]]", "", symbol) or None
            else:
                return nc[-1] if len(nc) > 0 else None

        return None

    def _get_sugar_symbol(self) -> str | None:
        nc = self.notation
        if not nc:
            return None

        if "(" in nc:
            start = nc.index("(")
            symbol = nc[:start]
            return re.sub(r"[\[\]]", "", symbol) or None

        if self.position_type == STARTING_POSITION_TYPE:
            if nc.endswith("]"):
                for i in range(len(nc) - 1, -1, -1):
                    if nc[i] == "[":
                        symbol = nc[:i]
                        return re.sub(r"[\[\]]", "", symbol) or None
            else:
                symbol = nc[: len(nc) - 1]
                return re.sub(r"[\[\]]", "", symbol) or None
        else:
            if nc.startswith("["):
                idx = nc.index("]")
                symbol = nc[: idx + 1]
                return re.sub(r"[\[\]]", "", symbol) or None
            else:
                return nc[0] if nc else None

        return None

    def _get_base_symbol(self) -> str | None:
        nc = self.notation
        if not nc or "(" not in nc:
            return None
        start = nc.index("(")
        end = nc.index(")")
        symbol = nc[start + 1 : end]
        return re.sub(r"[\[\]]", "", symbol) or None
