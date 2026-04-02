"""Marker so resources/ is a package; also exposes helpers to locate data files."""

from importlib import resources as _res
from pathlib import Path

_PACKAGE = "helm_notation.resources"


def resource_path(filename: str) -> Path:
    """Return a concrete filesystem path to *filename* inside this sub-package."""
    ref = _res.files(_PACKAGE).joinpath(filename)
    # Python ≥3.9 Traversable
    return Path(str(ref))


def read_text(filename: str, encoding: str = "utf-8") -> str:
    return resource_path(filename).read_text(encoding=encoding)


def read_bytes(filename: str) -> bytes:
    return resource_path(filename).read_bytes()
