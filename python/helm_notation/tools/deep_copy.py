"""Port of org.helm.notation.tools.DeepCopy."""

from __future__ import annotations

import copy
from typing import Any, List, TypeVar

from helm_notation.exceptions import NotationException

T = TypeVar("T")


def copy_list(src: List[T]) -> List[T]:
    """Deep copy a list.
    
    Args:
        src: Source list
    
    Returns:
        Deep copy of the list
    """
    return copy.deepcopy(src)


def copy_monomer(src: "Monomer") -> "Monomer":
    """Deep copy a Monomer object.
    
    Args:
        src: Source monomer
    
    Returns:
        Deep copy of the monomer
    """
    return copy.deepcopy(src)


def copy_attachment(src: "Attachment") -> "Attachment":
    """Deep copy an Attachment object.
    
    Args:
        src: Source attachment
    
    Returns:
        Deep copy of the attachment
    """
    return copy.deepcopy(src)


def copy_nucleotide(src: "Nucleotide") -> "Nucleotide":
    """Deep copy a Nucleotide object.
    
    Args:
        src: Source nucleotide
    
    Returns:
        Deep copy of the nucleotide
    """
    return copy.deepcopy(src)


def copy_serializable(src: Any) -> Any:
    """Deep copy any serializable object.
    
    Args:
        src: Source object
    
    Returns:
        Deep copy of the object
    """
    return copy.deepcopy(src)