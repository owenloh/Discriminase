"""Pack a CRISPR guide into a single integer and compare guides cheaply.

A guide of length L <= 32 is two bits per base (A=0, C=1, G=2, T=3) packed into one
``uint64``. The guide is stored **PAM-proximal (seed) first, in the high bits**:

    value = sum_j  code(base[j]) << (2 * (L - 1 - j))

so base[0] (the seed start, next to the PAM) is the most significant pair. The point
of that convention: sorting the integers groups guides by their seed, which turns the
seed lookup into a binary search (see ``index.py``).

Hamming distance over the 2-bit fields is one branchless expression -- fold each
differing 2-bit pair down to a single set bit, then popcount:

    x = a ^ b
    folded = (x | (x >> 1)) & 0x5555...   # one bit per differing base
    distance = popcount(folded)

This is exact and vectorizes over a numpy ``uint64`` block (the candidate set the
binary search returns), which is the hot path of the search.
"""
from typing import List

import numpy as np

# A=0, C=1, G=2, T=3.  N folds to A (a deliberate simplification: keeps every base in
# 2 bits, at the cost of treating ambiguous N as A).
_CODE = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 0,
         "a": 0, "c": 1, "g": 2, "t": 3, "n": 0}
_BASE = "ACGT"

# Low bit of every 2-bit field: 0b0101...0101 across 64 bits.
_FOLD_MASK = np.uint64(0x5555555555555555)
_ONE = np.uint64(1)


def pack_guide(seq: str, guide_length: int) -> int:
    """Pack ``seq`` (length ``guide_length``, PAM-proximal first) into an int."""
    if len(seq) != guide_length:
        raise ValueError(f"expected length {guide_length}, got {len(seq)}: {seq!r}")
    value = 0
    for base in seq:
        try:
            value = (value << 2) | _CODE[base]
        except KeyError:
            raise ValueError(f"invalid base {base!r} in {seq!r}")
    return value


def unpack_guide(value: int, guide_length: int) -> str:
    """Inverse of :func:`pack_guide`."""
    bases = []
    for _ in range(guide_length):
        bases.append(_BASE[value & 0b11])
        value >>= 2
    return "".join(reversed(bases))


def seed_of(value: int, guide_length: int, seed_len: int) -> int:
    """The ``seed_len`` PAM-proximal bases (the high bits) as an integer."""
    if seed_len > guide_length:
        raise ValueError("seed_len cannot exceed guide_length")
    return value >> (2 * (guide_length - seed_len))


def hamming(a: int, b: int) -> int:
    """Number of differing bases between two packed guides (scalar)."""
    x = a ^ b
    folded = (x | (x >> 1)) & 0x5555555555555555
    return folded.bit_count()


def _popcount_u64(arr: np.ndarray) -> np.ndarray:
    """popcount of each ``uint64`` in ``arr`` -> ``uint8`` counts."""
    bitwise_count = getattr(np, "bitwise_count", None)
    if bitwise_count is not None:        # numpy >= 2.0: one vectorized call
        return bitwise_count(arr)
    # Fallback: byte-wise popcount via a 256-entry lookup table.
    lut = _POPCOUNT_LUT
    as_bytes = arr.view(np.uint8).reshape(-1, 8)
    return lut[as_bytes].sum(axis=1).astype(np.uint8)


_POPCOUNT_LUT = np.array([bin(i).count("1") for i in range(256)], dtype=np.uint8)


def hamming_array(g: int, arr: np.ndarray) -> np.ndarray:
    """Hamming distance from packed guide ``g`` to every packed guide in ``arr``.

    ``arr`` is a ``uint64`` numpy array (typically a contiguous block sliced out of
    the sorted index by a seed lookup). Returns an array of distances.
    """
    gv = np.uint64(g)
    x = arr ^ gv
    folded = (x | (x >> _ONE)) & _FOLD_MASK
    return _popcount_u64(folded)


def seed_neighbors(seed_value: int, seed_len: int, max_mismatch: int) -> List[int]:
    """All seed integers within ``max_mismatch`` substitutions of ``seed_value``.

    For ``max_mismatch <= 1`` this is small (1 + 3*seed_len). Higher values are
    supported but grow combinatorially -- the search only ever needs 0 or 1 in
    practice.
    """
    if max_mismatch <= 0:
        return [seed_value]

    results = {seed_value}
    frontier = {seed_value}
    for _ in range(max_mismatch):
        nxt = set()
        for s in frontier:
            for pos in range(seed_len):
                shift = 2 * pos
                cur = (s >> shift) & 0b11
                for code in range(4):
                    if code == cur:
                        continue
                    neighbor = (s & ~(0b11 << shift)) | (code << shift)
                    if neighbor not in results:
                        nxt.add(neighbor)
        results.update(nxt)
        frontier = nxt
    return list(results)
