"""Data model for molecular fragments."""

from dataclasses import dataclass, field
from typing import ClassVar


@dataclass
class Fragment:
    """A fragment containing one or more monomers (e.g., a nucleotide or linker)."""

    SUPPORTED_SOURCES: ClassVar[list[str]] = []

    source: str = ""
    type: str = ""
    symbol: str = ""
    notation: str = ""

    def is_supported_source(self, source_name: str | None = None) -> bool:
        name = source_name if source_name is not None else self.source
        if name is None:
            return False
        return any(s.lower() == name.lower() for s in self.SUPPORTED_SOURCES)
