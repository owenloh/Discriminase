from bitarray import bitarray

from discriminase.trie import BitGuideTrie
from discriminase.bktree import BKTreeBitarray, hamming_distance


def bits(seq: str) -> bitarray:
    m = {'A': (0, 0), 'C': (0, 1), 'G': (1, 0), 'T': (1, 1)}
    out = bitarray()
    for base in seq:
        out.extend(m[base])
    return out


def test_trie_insert_and_prefix():
    trie = BitGuideTrie()
    guides = ['ACGT', 'ACGA', 'AC', 'GGT', 'T']
    for g in guides:
        trie.insert(bits(g))
    # every inserted guide is reachable as a prefix path
    for g in guides:
        assert trie.has_prefix(bits(g), len(g))
    # a path that was never inserted
    assert not trie.has_prefix(bits('C'), 1)
    # 'A' is an internal node (prefix of ACGT) even though 'A' alone wasn't inserted
    assert trie.has_prefix(bits('A'), 1)


def test_trie_terminal_count_and_dedup():
    trie = BitGuideTrie()
    for g in ['AC', 'ACGT', 'ACGA', 'GGT', 'T']:
        trie.insert(bits(g))
    assert trie.count_guides() == 5
    trie.insert(bits('ACGT'))  # duplicate must not change the count
    assert trie.count_guides() == 5


def test_hamming_distance():
    assert hamming_distance(bits('ACGT'), bits('ACGT')) == 0
    assert hamming_distance(bits('ACGT'), bits('ACGA')) == 1
    assert hamming_distance(bits('AAAA'), bits('TTTT')) == 4


def test_bktree_search_exists_and_nearest():
    bk = BKTreeBitarray(hamming_distance)
    for g in ['AAAAAA', 'AAAAAT', 'GGGGGG', 'ACGTAC']:
        bk.insert((bits(g), None, 0, False))
    # query within 1 mismatch of 'AAAAAA'
    assert bk.search_exists(bits('AAAAAC'), 1)
    # query far from everything
    assert not bk.search_exists(bits('TCTCTC'), 1)
    # nearest reports the actual minimum distance
    near = bk.nearest(bits('AAAAAC'), 2)
    assert near is not None and near[0] == 1
