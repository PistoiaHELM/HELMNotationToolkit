"""Custom exceptions mirroring the Java toolkit's exception hierarchy."""


class NotationException(Exception):
    """Raised for HELM notation parsing / validation errors."""


class MonomerException(Exception):
    """Raised for monomer-related errors."""


class StructureException(Exception):
    """Raised for chemical structure errors."""


class CalculationException(Exception):
    """Raised for molecular-property calculation errors."""
