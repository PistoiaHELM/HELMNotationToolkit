"""Exception classes for HELM notation."""


class HelmError(Exception):
    """Base exception for all HELM errors."""


class NotationException(HelmError):
    """Exception for notation-related errors."""


class MonomerException(HelmError):
    """Exception for monomer-related errors."""


class StructureException(HelmError):
    """Exception for structure-related errors."""


class CalculationException(HelmError):
    """Exception for calculation-related errors."""
