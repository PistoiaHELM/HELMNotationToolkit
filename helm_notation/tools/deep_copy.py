"""Deep copy utility (replaces Java serialization-based DeepCopy)."""

from __future__ import annotations

import copy
from typing import TypeVar

T = TypeVar("T")


def copy_obj(obj: T) -> T:
    """Deep copy an object."""
    return copy.deepcopy(obj)
