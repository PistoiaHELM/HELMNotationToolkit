"""Shared test fixtures for HELM Notation Toolkit tests."""

from __future__ import annotations

import os
from pathlib import Path

import pytest


TESTS_DIR = Path(__file__).parent
RESOURCES_DIR = TESTS_DIR / "resources"


@pytest.fixture(scope="session")
def monomer_factory():
    """Initialize and return the MonomerFactory singleton."""
    from helm_notation.monomer_factory import get_instance, finalize

    factory = get_instance()
    yield factory
    finalize()


@pytest.fixture(scope="session")
def nucleotide_factory():
    """Initialize and return the NucleotideFactory singleton."""
    from helm_notation.nucleotide_factory import get_instance, finalize

    factory = get_instance()
    yield factory
    finalize()


@pytest.fixture(scope="session")
def monomer_store(monomer_factory):
    """Return the MonomerStore from the factory."""
    return monomer_factory.monomer_store


def load_test_resource(filename: str) -> str:
    """Load a test resource file as string."""
    return (RESOURCES_DIR / filename).read_text(encoding="utf-8")
