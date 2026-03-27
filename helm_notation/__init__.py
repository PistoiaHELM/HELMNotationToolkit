"""HELM Notation Toolkit - Python implementation using RDKit.

A Python toolkit for parsing, validating, and manipulating HELM
(Hierarchical Editing Language for Macromolecules) notation.

Quick start::

    import helm_notation as helm

    # Validate a notation
    helm.validate("PEPTIDE1{A.G.K}$$$$V2.0")

    # Get molecular info
    info = helm.get_molecule_info("PEPTIDE1{A.G.K}$$$$V2.0")
    print(info.molecular_weight)

    # Convert to SMILES
    smiles = helm.get_smiles("PEPTIDE1{A.G.K}$$$$V2.0")
"""

__version__ = "2.0.0"


def validate(notation: str) -> None:
    """Validate a HELM notation string. Raises on invalid notation."""
    from helm_notation.tools.complex_notation_parser import validate_complex_notation
    validate_complex_notation(notation)


def get_smiles(notation: str) -> str:
    """Convert HELM notation to canonical SMILES."""
    from helm_notation.tools.complex_notation_parser import get_complex_polymer_smiles
    return get_complex_polymer_smiles(notation)


def get_molecule_info(notation: str):
    """Get molecular info (MW, formula, exact mass) for a HELM notation."""
    from helm_notation.tools.complex_notation_parser import get_molecule_info as _get_info
    return _get_info(notation)


def canonicalize(notation: str) -> str:
    """Return the canonical form of a HELM notation string."""
    from helm_notation.tools.complex_notation_parser import get_canonical_notation
    return get_canonical_notation(notation)
