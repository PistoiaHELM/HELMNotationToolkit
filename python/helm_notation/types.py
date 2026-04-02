from dataclasses import dataclass


@dataclass(frozen=True)
class MoleculeInfo:
    molecular_weight: float
    exact_mass: float
    molecular_formula: str
