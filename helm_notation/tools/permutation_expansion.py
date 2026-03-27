"""Permutation and expansion utilities for HELM notation processing.

Provides functions for generating permutations and expanding
combinatorial lists used in complex polymer notation handling.
"""

from __future__ import annotations

import copy


def linearize(
    original_list: list[list[list[str]]],
) -> list[list[str]]:
    """Flatten inner arrays into simple string lists.

    Args:
        original_list: List of lists of string arrays.

    Returns:
        List of flattened string lists.
    """
    result: list[list[str]] = []
    for inner_list in original_list:
        tmp: list[str] = []
        for sa in inner_list:
            for item in sa:
                tmp.append(item)
        result.append(tmp)
    return result


def expand(
    parent: list[list[list[str]]],
    child: list[list[str]],
) -> None:
    """Expand parent list by combining with child list elements (Cartesian product).

    Modifies parent in place.

    Args:
        parent: List of lists of string arrays to expand.
        child: List of string arrays to combine with each parent entry.
    """
    if len(parent) == 0:
        parent.append(child)
    else:
        remove: list[list[list[str]]] = []
        keep: list[list[list[str]]] = []
        for i in range(len(parent)):
            tmp = parent[i]
            remove.append(tmp)
            for j in range(len(child)):
                copied = copy.deepcopy(tmp)
                copied.append(child[j])
                keep.append(copied)
        for item in remove:
            parent.remove(item)
        parent.extend(keep)


def permutate_list(result: list[list[str]], items: list[str]) -> None:
    """Generate all permutations of items and append to result.

    Args:
        result: List to append permutations to.
        items: List of strings to permute.
    """
    a = list(items)
    n = len(items)
    _permutate(result, a, n)


def permutate(
    result: list[list[str]], a: list[str], n: int
) -> None:
    """Generate permutations using Heap's algorithm.

    Args:
        result: List to append permutations to.
        a: Array of strings being permuted.
        n: Current permutation size.
    """
    _permutate(result, a, n)


def _permutate(
    result: list[list[str]], a: list[str], n: int
) -> None:
    """Internal permutation via Heap's algorithm."""
    if n == 1:
        result.append(list(a))

    for i in range(n):
        _swap(a, i, n - 1)
        _permutate(result, a, n - 1)
        _swap(a, i, n - 1)


def _swap(a: list[str], i: int, j: int) -> None:
    """Swap elements at indices i and j."""
    a[i], a[j] = a[j], a[i]
