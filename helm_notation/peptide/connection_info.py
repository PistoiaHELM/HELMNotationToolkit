"""Data model for connection information during peptide structure parsing."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class ConnectionInfo:
    """Tracks where a connection occurs in the monomer ID / connection map lists."""

    list_index: int = -1
    monomer_index: int = -1
    monomer_id: str = ""
    rgroup: str = ""
