"""GuideIndex: query correctness, the safety invariant, persistence."""
import numpy as np

from discriminase.index import GuideIndex
from discriminase.pack import hamming, pack_guide, seed_of, unpack_guide


def _brute(g, commensals, L, seed_len, d, s):
    """Reference: nearest collision under the seed-anchored model, or None."""
    sg = seed_of(g, L, seed_len)
    best = None
    for c in commensals:
        if hamming(seed_of(c, L, seed_len), sg) <= s and hamming(g, c) <= d:
            dd = hamming(g, c)
            best = dd if best is None else min(best, dd)
    return best


def _packed(seqs, L):
    return np.array([pack_guide(s, L) for s in seqs], dtype=np.uint64)


def test_exact_hit_is_distance_zero():
    L, sl = 8, 4
    commensals = _packed(["ACGTACGT", "TTTTAAAA"], L)
    idx = GuideIndex.from_packed(commensals, L, sl)
    g = pack_guide("ACGTACGT", L)
    hit = idx.query(g, total_max_mismatch=2, seed_max_mismatch=1)
    assert hit is not None and hit[0] == 0


def test_far_guide_is_spared():
    L, sl = 8, 4
    commensals = _packed(["ACGTACGT"], L)
    idx = GuideIndex.from_packed(commensals, L, sl)
    g = pack_guide("TGCATGCA", L)               # very different
    assert idx.is_spared(g, total_max_mismatch=2, seed_max_mismatch=1)


def test_seed_mismatch_beyond_s_not_flagged():
    # Commensal identical to g except in the seed; with s=0 it must NOT count.
    L, sl = 8, 4
    g_seq = "ACGTACGT"
    c_seq = "TCGTACGT"                            # differs only at seed position 0
    idx = GuideIndex.from_packed(_packed([c_seq], L), L, sl)
    g = pack_guide(g_seq, L)
    assert idx.is_spared(g, total_max_mismatch=4, seed_max_mismatch=0)   # spared
    assert not idx.is_spared(g, total_max_mismatch=4, seed_max_mismatch=1)  # flagged


def test_provenance_org_id():
    L, sl = 8, 4
    commensals = _packed(["ACGTACGT", "TTTTAAAA"], L)
    org = np.array([7, 3], dtype=np.uint8)
    idx = GuideIndex.from_packed(commensals, L, sl, org_ids=org)
    hit = idx.query(pack_guide("TTTTAAAA", L), 2, 1)
    assert hit is not None and hit[1] == 3


def test_dedup():
    L, sl = 8, 4
    commensals = _packed(["ACGTACGT", "ACGTACGT", "TTTTAAAA"], L)
    idx = GuideIndex.from_packed(commensals, L, sl)
    assert len(idx) == 2


def test_matches_bruteforce_random():
    """Safety invariant: the index never misses a collision the brute force finds."""
    rng = np.random.default_rng(42)
    L, sl, d, s = 12, 6, 3, 1
    commensal_seqs = ["".join(rng.choice(list("ACGT"), L)) for _ in range(800)]
    commensals = _packed(commensal_seqs, L)
    idx = GuideIndex.from_packed(commensals, L, sl)

    for _ in range(400):
        g_seq = "".join(rng.choice(list("ACGT"), L))
        g = pack_guide(g_seq, L)
        ref = _brute(g, commensals.tolist(), L, sl, d, s)
        hit = idx.query(g, d, s)
        if ref is None:
            assert hit is None, f"index flagged a spared guide {g_seq}"
        else:
            assert hit is not None, f"index MISSED a collision for {g_seq} (unsafe!)"
            assert hit[0] == ref


def test_export_web_roundtrip(tmp_path):
    L, sl = 10, 5
    commensals = _packed(["ACGTACGTAC", "TTTTAAAACC", "GGGGCCCCAA"], L)
    org = np.array([0, 1, 2], dtype=np.uint16)
    idx = GuideIndex.from_packed(commensals, L, sl, org_ids=org,
                                 organisms=[{"name": "a"}, {"name": "b"}, {"name": "c"}])
    prefix = str(tmp_path / "panel")
    manifest = idx.export_web(prefix)

    # the .guides.f64 is exactly the sorted guides, readable as plain float64
    back = np.fromfile(prefix + ".guides.f64", dtype="<f8")
    assert np.array_equal(back.astype(np.uint64), idx.guides)
    assert manifest["guide_length"] == L and manifest["seed_len"] == sl
    assert manifest["n_guides"] == len(idx)
    orgs = np.fromfile(prefix + ".orgs.u16", dtype="<u2")
    assert orgs.shape[0] == len(idx)


def test_save_load_roundtrip(tmp_path):
    L, sl = 10, 5
    commensals = _packed(["ACGTACGTAC", "TTTTAAAACC"], L)
    org = np.array([1, 2], dtype=np.uint8)
    organisms = [{"name": "org-a"}, {"name": "org-b"}, {"name": "org-c"}]
    idx = GuideIndex.from_packed(commensals, L, sl, org_ids=org, organisms=organisms)
    prefix = str(tmp_path / "panel")
    idx.save(prefix)

    loaded = GuideIndex.load(prefix)
    assert len(loaded) == 2
    assert loaded.guide_length == L
    assert loaded.seed_len == sl
    assert loaded.organisms == organisms
    g = pack_guide("ACGTACGTAC", L)
    assert loaded.query(g, 2, 1)[0] == 0
    assert unpack_guide(int(loaded.guides[0]), L)  # decodes without error
