"""Vectorized PAM extraction vs a slow reference; GC fraction."""
import numpy as np

from discriminase.genome import extract_packed_guides, gc_fraction
from discriminase.pack import pack_guide, unpack_guide


def _slow_extract(seq, L, pam, gap):
    """Reference: scan forward + reverse-complement the obvious slow way."""
    out = []

    def scan(s):
        for i in range(len(s) - len(pam) + 1):
            if s[i : i + len(pam)] == pam:
                start = i + len(pam) + gap
                if start + L <= len(s):
                    out.append(s[start : start + L])

    rc = seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    scan(seq)
    scan(rc)
    return out


def test_extraction_matches_slow_reference():
    rng = np.random.default_rng(7)
    L, pam, gap = 6, "TTT", 1
    seq = "".join(rng.choice(list("ACGT"), 400))
    expected = sorted(pack_guide(g, L) for g in _slow_extract(seq, L, pam, gap))
    got = sorted(int(x) for x in extract_packed_guides(seq, L, pam, gap))
    assert got == expected
    assert len(got) > 0  # the random seq should contain some TTT PAMs


def test_known_geometry():
    # PAM=TTT at index 0, gap=1, spacer starts at 0+3+1=4.
    seq = "TTTA" + "ACGTAC"          # spacer = ACGTAC
    guides = extract_packed_guides(seq, 6, "TTT", 1)
    decoded = {unpack_guide(int(x), 6) for x in guides}
    assert "ACGTAC" in decoded       # PAM-proximal first, forward strand


def _slow_skip(seq, L, pam, gap):
    """5' reference that drops any window containing a non-ACGT base, both strands."""
    out = []
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    rc = "".join(comp.get(b, b) for b in seq)[::-1]   # complement; N stays N

    def scan(s):
        for i in range(len(s) - len(pam) + 1):
            if s[i:i + len(pam)] == pam:
                st = i + len(pam) + gap
                if st + L <= len(s) and all(b in "ACGT" for b in s[st:st + L]):
                    out.append(s[st:st + L])

    scan(seq)
    scan(rc)
    return out


def test_ambiguous_bases_are_skipped_not_folded():
    rng = np.random.default_rng(5)
    seq = list("".join(rng.choice(list("ACGT"), 300)))
    for i in (50, 120, 200):
        seq[i] = "N"                      # inject ambiguous bases
    seq = "".join(seq)
    L, pam, gap = 6, "TTT", 1
    want = sorted(pack_guide(g, L) for g in _slow_skip(seq, L, pam, gap))
    got = sorted(int(x) for x in extract_packed_guides(seq, L, pam, gap, "5prime"))
    assert got == want                    # windows touching an N are dropped, not folded to A


def test_gc_fraction():
    g = pack_guide("CCCCGGGGAA", 10)  # 8 GC / 10 = 0.8
    frac = gc_fraction(np.array([g], dtype=np.uint64), 10)
    assert abs(frac[0] - 0.8) < 1e-9
    g2 = pack_guide("AAAATTTTAA", 10)  # 0 GC
    assert gc_fraction(np.array([g2], dtype=np.uint64), 10)[0] == 0.0
