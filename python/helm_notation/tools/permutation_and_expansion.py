"""Port of org.helm.notation.tools.PermutationAndExpansion."""

from __future__ import annotations

from typing import List

from helm_notation.tools.deep_copy import copy_list


def linearize(original_list: List[List[str]]) -> List[List[str]]:
    """Linearize a list of lists.
    
    Args:
        original_list: List of lists of strings
    
    Returns:
        Flattened list of lists
    """
    result: List[List[str]] = []
    for inner_list in original_list:
        result.append(list(inner_list))
    return result


def expand(parent: List[List[str]], child: List[str]) -> None:
    """Expand parent list with child combinations.
    
    Modifies parent in place.
    
    Args:
        parent: Parent list of lists to expand
        child: Child list to combine with parent
    
    Raises:
        Exception: If expansion fails
    """
    if len(parent) == 0:
        parent.append(list(child))
    else:
        remove = [list(p) for p in parent]
        keep: List[List[str]] = []
        
        for i in range(len(parent)):
            tmp = parent[i]
            for j in range(len(child)):
                l = copy_list(tmp)
                l.append(child[j])
                keep.append(l)
        
        parent.clear()
        parent.extend(keep)


def permutate(l: List[str], ls: List[str]) -> None:
    """Generate all permutations of ls and add to l.
    
    Modifies l in place.
    
    Args:
        l: Output list to add permutations to
        ls: Source list to permute
    """
    a = list(ls)
    n = len(a)
    _permutate(l, a, n)


def _permutate(l: List[str], a: List[str], n: int) -> None:
    """Recursive helper for permutation generation."""
    if n == 1:
        result = list(a)
        l.append(result)
    else:
        for i in range(n):
            _swap(a, i, n - 1)
            _permutate(l, a, n - 1)
            _swap(a, i, n - 1)


def _swap(a: List[str], i: int, j: int) -> None:
    """Swap elements at indices i and j."""
    a[i], a[j] = a[j], a[i]