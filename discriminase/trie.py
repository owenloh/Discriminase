"""Prefix trie over 2-bit guides — the fast seed prefilter.

Each node is a ``(is_terminal, children_dict)`` tuple. Tuples beat a node class
here: less per-node overhead over millions of nodes. ``has_prefix`` answers
"does any protected guide share this guide's first k bases exactly?" in O(k),
which throws out the overwhelming majority of candidates before the BK-tree runs.
"""
import pickle

from bitarray import bitarray


class BitGuideTrie:
    def __init__(self):
        self.root = (False, {})

    def insert(self, bits: bitarray) -> None:
        # Single pass: descend once, remembering the leaf's parent dict + slot,
        # then flip the terminal flag in O(1). (Nodes are immutable tuples, so the
        # terminal bit is set by replacing the leaf entry in its parent's dict.)
        node = self.root
        parent_children = None
        idx = None
        for i in range(0, len(bits), 2):
            idx = (bits[i] << 1) | bits[i + 1]
            _, children = node
            if idx not in children:
                children[idx] = (False, {})
            parent_children = children
            node = children[idx]
        # Mark the end of a guide (skip empty inserts, which have no leaf)
        if parent_children is not None and not node[0]:
            parent_children[idx] = (True, node[1])

    def has_prefix(self, bits: bitarray, prefix_len: int) -> bool:
        node = self.root
        for i in range(0, 2 * prefix_len, 2):
            idx = (bits[i] << 1) | bits[i + 1]
            _, children = node
            if idx not in children:
                return False
            node = children[idx]
        return True

    def save(self, filename: str) -> None:
        with open(filename, 'wb') as f:
            pickle.dump(self, f)

    @staticmethod
    def load(filename: str) -> 'BitGuideTrie':
        with open(filename, 'rb') as f:
            return pickle.load(f)

    def count_guides(self) -> int:
        count = 0
        stack = [self.root]
        while stack:
            is_terminal, children = stack.pop()
            if is_terminal:
                count += 1
            stack.extend(children.values())
        return count
