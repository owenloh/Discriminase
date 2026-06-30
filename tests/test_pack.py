"""Pack/unpack and Hamming primitives."""
import numpy as np

from discriminase.pack import (
    hamming,
    hamming_array,
    pack_guide,
    seed_neighbors,
    seed_of,
    unpack_guide,
)


def _naive_hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def test_roundtrip():
    for seq in ["ACGT", "TTTTAAAACCCCGGGG", "ACGTACGTACGTACGTACGTACG"]:
        assert unpack_guide(pack_guide(seq, len(seq)), len(seq)) == seq


def test_n_folds_to_a():
    assert unpack_guide(pack_guide("NNNA", 4), 4) == "AAAA"


def test_seed_is_high_bits():
    # PAM-proximal first -> seed is the leading bases.
    L, sl = 8, 3
    seq = "ACGTACGT"
    v = pack_guide(seq, L)
    assert unpack_guide(seed_of(v, L, sl), sl) == seq[:sl]


def test_hamming_scalar_matches_naive():
    rng = np.random.default_rng(0)
    bases = "ACGT"
    for _ in range(500):
        a = "".join(rng.choice(list(bases), 12))
        b = "".join(rng.choice(list(bases), 12))
        assert hamming(pack_guide(a, 12), pack_guide(b, 12)) == _naive_hamming(a, b)


def test_hamming_array_matches_scalar():
    rng = np.random.default_rng(1)
    L = 16
    guides = [pack_guide("".join(rng.choice(list("ACGT"), L)), L) for _ in range(300)]
    g = guides[0]
    arr = np.array(guides, dtype=np.uint64)
    vec = hamming_array(g, arr)
    scal = np.array([hamming(g, x) for x in guides])
    assert np.array_equal(vec, scal)


def test_seed_neighbors_count_and_membership():
    L, sl = 10, 5
    v = pack_guide("ACGTACGTAC", L)
    s0 = seed_of(v, L, sl)
    assert seed_neighbors(s0, sl, 0) == [s0]
    nb1 = set(seed_neighbors(s0, sl, 1))
    assert len(nb1) == 1 + 3 * sl            # self + 3 substitutions per position
    assert s0 in nb1
    # every 1-neighbor is exactly one base off
    for nb in nb1:
        assert hamming(nb, s0) <= 1
