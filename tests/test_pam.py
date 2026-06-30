"""IUPAC PAM matching, 5'/3' sides, and a SpCas9 (3' NGG) end-to-end run."""
import numpy as np

from discriminase.config import GuideFinderConfig, NUCLEASE_PRESETS
from discriminase.database import build_index
from discriminase.genome import extract_packed_guides
from discriminase.pack import pack_guide
from discriminase.search import find_sparing_guides

_IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG",
          "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT",
          "V": "ACG", "N": "ACGT"}


def _slow(seq, L, pam, side):
    """Reference extraction in 5'->3' guide order, both strands."""
    out = []

    def scan(s):
        for i in range(len(s) - len(pam) + 1):
            if all(s[i + j] in _IUPAC[pam[j]] for j in range(len(pam))):
                st = (i + len(pam)) if side == "5prime" else (i - L)  # protospacer adjacent to PAM
                if 0 <= st and st + L <= len(s):
                    out.append(s[st:st + L])

    rc = seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    scan(seq)
    scan(rc)
    return out


def _to_packed(guide_5to3, side):
    g = guide_5to3 if side == "5prime" else guide_5to3[::-1]   # stored PAM-proximal-first
    return pack_guide(g, len(g))


def test_iupac_both_sides_match_reference():
    rng = np.random.default_rng(11)
    seq = "".join(rng.choice(list("ACGT"), 400))
    for pam, side in [("TTT", "5prime"), ("NGG", "3prime"),
                      ("TNG", "5prime"), ("RYN", "3prime"),
                      ("TTTV", "5prime"), ("TTTN", "5prime")]:
        L = 6
        want = sorted(_to_packed(g, side) for g in _slow(seq, L, pam, side))
        got = sorted(int(x) for x in extract_packed_guides(seq, L, pam, side))
        assert got == want, f"mismatch for pam={pam} side={side}"


def test_presets_apply():
    cfg = GuideFinderConfig(**NUCLEASE_PRESETS["spcas9"])
    assert cfg.pam == "NGG" and cfg.pam_side == "3prime" and cfg.guide_length == 20


SHARED = "ACGTACGTACGTACGTACG"      # 19 nt, ~50% GC
UNIQUE = "TGCATGCATGCATGCATGC"      # 19 nt, ~50% GC, far from SHARED
PAD = "AT" * 12                     # no "GG", spawns no forward NGG


def test_spcas9_end_to_end(tmp_path):
    cfg = GuideFinderConfig(**NUCLEASE_PRESETS["spcas9"])
    cfg.guide_length = len(SHARED)
    cfg.db_dir = str(tmp_path)
    cfg.db_prefix = "cas9"
    cfg.email = ""
    cfg.max_guides = 0

    def cassette(spacer):              # 5'-[spacer][NGG]-3'  (PAM adjacent, on 3')
        return spacer + "AGG"

    commensal = PAD + cassette(SHARED) + PAD
    target = PAD + cassette(SHARED) + PAD + cassette(UNIQUE) + PAD

    index = build_index([{"name": "c", "seq": commensal}], cfg,
                        get_sequence=lambda src: (src["seq"], {}), progress=lambda *_: None)
    rows = find_sparing_guides(target, index, cfg, progress=False)
    kept = {r["guide_sequence"] for r in rows}

    assert UNIQUE in kept, "unique 3'-PAM guide must survive (and be reported 5'->3')"
    assert SHARED not in kept, "shared 3'-PAM guide must be dropped"
