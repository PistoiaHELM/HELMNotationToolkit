"""Generic tree data structure for tracking peptide fragmentation hierarchy.

Originally based on a Java implementation by Yohann Coppel (Google, 2007),
licensed under the Apache License 2.0.
"""

from __future__ import annotations

from typing import Generic, TypeVar

T = TypeVar("T")

_INDENT = 2


class Tree(Generic[T]):
    """A simple n-ary tree that also maintains a global lookup map.

    The *locate* dict is shared across all nodes in the same tree so that
    any node can be found in O(1) time by its head value.
    """

    def __init__(self, head: T) -> None:
        self.head: T = head
        self.leafs: list[Tree[T]] = []
        self.parent: Tree[T] | None = None
        self._locate: dict[int, Tree[T]] = {id(head): self}

    # ------------------------------------------------------------------
    # Mutation helpers
    # ------------------------------------------------------------------

    def add_leaf_under(self, root: T, leaf: T) -> None:
        """Add *leaf* as a child of *root* (which must already be in the tree).

        If *root* is not found, it is first added as a child of ``self``
        and then *leaf* is added under it.
        """
        root_key = id(root)
        if root_key in self._locate:
            self._locate[root_key].add_leaf(leaf)
        else:
            self.add_leaf(root).add_leaf(leaf)

    def add_leaf(self, leaf: T) -> Tree[T]:
        """Add *leaf* as a direct child and return the new subtree node."""
        child: Tree[T] = Tree(leaf)
        self.leafs.append(child)
        child.parent = self
        child._locate = self._locate
        self._locate[id(leaf)] = child
        return child

    def set_as_parent(self, parent_root: T) -> Tree[T]:
        """Create a new root above this tree and return it."""
        new_root: Tree[T] = Tree(parent_root)
        new_root.leafs.append(self)
        self.parent = new_root
        new_root._locate = self._locate
        new_root._locate[id(self.head)] = self
        new_root._locate[id(parent_root)] = new_root
        return new_root

    # ------------------------------------------------------------------
    # Query helpers
    # ------------------------------------------------------------------

    def get_tree(self, element: T) -> Tree[T] | None:
        """Return the subtree whose head is *element*, or ``None``."""
        return self._locate.get(id(element))

    def get_successors(self, root: T) -> list[T]:
        """Return the immediate children values of *root*."""
        tree = self.get_tree(root)
        if tree is None:
            return []
        return [child.head for child in tree.leafs]

    def get_sub_trees(self) -> list[Tree[T]]:
        """Return the list of direct child subtrees."""
        return self.leafs

    @staticmethod
    def get_successors_from_collection(
        of: T, collection: list[Tree[T]]
    ) -> list[T]:
        """Search through a collection of trees for successors of *of*."""
        for tree in collection:
            if id(of) in tree._locate:
                return tree.get_successors(of)
        return []

    def get_depth(self) -> int:
        """Return the maximum depth of the tree (root is depth 1)."""
        steps: list[int] = []
        self._walk(self, steps, 1)
        return max(steps) if steps else 1

    # ------------------------------------------------------------------
    # Internal / display
    # ------------------------------------------------------------------

    def _walk(self, tree: Tree[T], steps: list[int], level: int) -> None:
        steps.append(level)
        if tree.leafs:
            next_level = level + 1
            for subtree in tree.leafs:
                self._walk(subtree, steps, next_level)

    def __str__(self) -> str:
        return self._print_tree(0)

    def _print_tree(self, increment: int) -> str:
        indent = " " * increment
        s = f"{indent}{self.head}"
        for child in self.leafs:
            s += "\n" + child._print_tree(increment + _INDENT)
        return s
