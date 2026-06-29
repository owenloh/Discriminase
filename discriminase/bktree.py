"""BK-tree over 2-bit guides using Hamming distance.

A BK-tree indexes items by a metric distance so that an approximate-match query
only descends branches within ``[d - max_dist, d + max_dist]`` of each node —
the triangle inequality prunes the rest. Here it answers "is any protected guide
within `max_dist` mismatches of this candidate?" Each node stores
``(bits, genome, pos, is_rev)`` so a hit can be traced back to its source.
"""
import pickle

from bitarray import bitarray

try:
    # Optional C extension; ~10x faster. Build with: see README "Speed".
    from discriminase.bktree_cython import hamming_distance
except ImportError:
    def hamming_distance(bits1, bits2):
        assert len(bits1) == len(bits2), "Bitarrays must be of the same length"
        xor_result = bits1 ^ bits2
        distance = 0
        for i in range(0, len(xor_result), 2):
            if xor_result[i] or xor_result[i + 1]:
                distance += 1
        return distance


class BKTreeBitarray:
    def __init__(self, distance_fn):
        self.distance_fn = distance_fn
        self.tree = None  # (item, {distance: subtree}), item = (bits, genome, pos, is_rev)

    def insert(self, item):
        bits = item[0]
        if self.tree is None:
            self.tree = (item, {})
            return
        node = self.tree
        while True:
            node_bits = node[0][0]
            dist = self.distance_fn(bits, node_bits)
            children = node[1]
            if dist in children:
                node = children[dist]
            else:
                children[dist] = (item, {})
                break

    def search(self, bits, max_dist):
        """Return every node within `max_dist` of `bits`."""
        if self.tree is None:
            return []
        result = []
        nodes = [self.tree]
        while nodes:
            node = nodes.pop()
            dist = self.distance_fn(bits, node[0][0])
            if dist <= max_dist:
                result.append(node)
            for d in range(dist - max_dist, dist + max_dist + 1):
                child = node[1].get(d)
                if child:
                    nodes.append(child)
        return result

    def search_exists(self, bits, max_dist):
        """Fast existence check: True as soon as any node is within `max_dist`."""
        if self.tree is None:
            return False
        nodes = [self.tree]
        while nodes:
            node = nodes.pop()
            dist = self.distance_fn(bits, node[0][0])
            if dist <= max_dist:
                return True
            for d in range(dist - max_dist, dist + max_dist + 1):
                child = node[1].get(d)
                if child:
                    nodes.append(child)
        return False

    def nearest(self, bits, max_dist):
        """Return (distance, item) of the closest node within `max_dist`, else None.

        Useful for *explaining* a rejection: which protected guide collided, and
        how close it was. Stops early once distance 0 is found.
        """
        if self.tree is None:
            return None
        best = None
        nodes = [self.tree]
        while nodes:
            node = nodes.pop()
            dist = self.distance_fn(bits, node[0][0])
            if dist <= max_dist and (best is None or dist < best[0]):
                best = (dist, node[0])
                if dist == 0:
                    break
            for d in range(dist - max_dist, dist + max_dist + 1):
                child = node[1].get(d)
                if child:
                    nodes.append(child)
        return best

    def save(self, filename: str):
        with open(filename, 'wb') as f:
            pickle.dump(self, f)

    @staticmethod
    def load(filename: str) -> 'BKTreeBitarray':
        with open(filename, 'rb') as f:
            return pickle.load(f)

    def count_guides(self):
        def count(node):
            if node is None:
                return 0
            _, children = node
            return 1 + sum(count(c) for c in children.values())
        return count(self.tree)
