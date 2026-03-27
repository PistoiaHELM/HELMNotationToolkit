"""Data model for molecular information."""

from dataclasses import dataclass, field


@dataclass
class MoleculeInfo:
    """Holds molecular information including weight, formula, and exact mass."""

    molecular_weight: float = 0.0
    molecular_formula: str = ""
    exact_mass: float = 0.0
