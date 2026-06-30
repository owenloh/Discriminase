"""Vectorized genome encoding and PAM-anchored guide extraction.

A genome becomes a ``uint8`` array of base codes (A=0, C=1, G=2, T=3). Anything that
is not a definite A/C/G/T -- N, the IUPAC ambiguity codes (R,Y,S,W,K,M,B,D,H,V), gaps,
stray characters -- is tracked by :func:`valid_mask`, and any guide window that would
contain such a base is dropped (we never emit or index a guide with an ambiguous base).
PAM matching and guide packing are numpy operations over the array -- no Python
per-position loop, so a 5 Mbp genome scans in well under a second and uses only a few MB.

Both strands are covered by scanning the sequence and its reverse complement with the
*same* code path, so every extracted guide comes out PAM-proximal-first (seed in the
high bits) to match :mod:`pack`.
"""
import numpy as np

# ASCII -> 2-bit code. Default 0 (A); C/G/T set explicitly. N and any other byte -> A.
_ASCII_TO_CODE = np.zeros(256, dtype=np.uint8)
_ASCII_TO_CODE[ord("C")] = 1
_ASCII_TO_CODE[ord("G")] = 2
_ASCII_TO_CODE[ord("T")] = 3
_ASCII_TO_CODE[ord("c")] = 1
_ASCII_TO_CODE[ord("g")] = 2
_ASCII_TO_CODE[ord("t")] = 3


def encode(seq: str) -> np.ndarray:
    """DNA string -> ``uint8`` array of 2-bit base codes."""
    raw = np.frombuffer(seq.encode("ascii", "ignore"), dtype=np.uint8)
    return _ASCII_TO_CODE[raw]


def reverse_complement(codes: np.ndarray) -> np.ndarray:
    """Reverse complement in code space (A<->T, C<->G == ``3 - code``, reversed)."""
    return (3 - codes)[::-1]


# True only for A/C/G/T. Everything else -- N and the IUPAC ambiguity codes
# (R,Y,S,W,K,M,B,D,H,V), gaps, digits, anything -- is NOT a definite base.
_IS_ACGT = np.zeros(256, dtype=bool)
for _c in b"ACGTacgt":
    _IS_ACGT[_c] = True


def valid_mask(seq: str) -> np.ndarray:
    """Boolean array, True where the base is an unambiguous A/C/G/T."""
    raw = np.frombuffer(seq.encode("ascii", "ignore"), dtype=np.uint8)
    return _IS_ACGT[raw]


def _filter_valid_windows(starts: np.ndarray, valid: np.ndarray, L: int) -> np.ndarray:
    """Keep only windows ``[s, s+L)`` that contain no ambiguous/invalid base."""
    if starts.size == 0:
        return starts
    cum = np.concatenate(([0], np.cumsum((~valid).astype(np.int64))))
    return starts[(cum[starts + L] - cum[starts]) == 0]


# IUPAC nucleotide codes -> the set of base codes each one allows.
_IUPAC = {
    "A": (0,), "C": (1,), "G": (2,), "T": (3,),
    "R": (0, 2), "Y": (1, 3), "S": (1, 2), "W": (0, 3),
    "K": (2, 3), "M": (0, 1),
    "B": (1, 2, 3), "D": (0, 2, 3), "H": (0, 1, 3), "V": (0, 1, 2),
    "N": (0, 1, 2, 3),
}


def _pam_sets(pam: str):
    """Per-position allowed code sets; ``None`` for a fully-unconstrained (N) position."""
    out = []
    for b in pam.upper():
        if b not in _IUPAC:
            raise ValueError(f"invalid PAM letter {b!r} (use IUPAC codes)")
        allowed = _IUPAC[b]
        out.append(None if len(allowed) == 4 else allowed)
    return out


def _pam_positions(codes: np.ndarray, pam_sets) -> np.ndarray:
    """Start indices where the (possibly ambiguous) PAM matches ``codes``."""
    n, p = codes.shape[0], len(pam_sets)
    if n < p:
        return np.empty(0, dtype=np.int64)
    window = n - p + 1
    match = np.ones(window, dtype=bool)
    for j, allowed in enumerate(pam_sets):
        if allowed is None:
            continue                       # N -> any base, no constraint
        col = codes[j : j + window]
        hit = np.zeros(window, dtype=bool)
        for c in allowed:
            hit |= col == c
        match &= hit
    return np.nonzero(match)[0]


def _pack_windows(codes: np.ndarray, starts: np.ndarray, L: int, reverse: bool) -> np.ndarray:
    """Pack length-``L`` windows into ``uint64``, PAM-proximal base in the high bits.

    ``reverse`` flips the window so the PAM-proximal end leads (needed for a 3' PAM,
    where the guide sits *before* the PAM and its 3' end is PAM-proximal).
    """
    if starts.size == 0:
        return np.empty(0, dtype=np.uint64)
    windows = codes[starts[:, None] + np.arange(L)]
    if reverse:
        windows = windows[:, ::-1]
    windows = windows.astype(np.uint64)
    shifts = (2 * np.arange(L - 1, -1, -1)).astype(np.uint64)
    return (windows << shifts).sum(axis=1).astype(np.uint64)


def _strand_guides(codes, valid, L, pam_sets, side):
    """(window_starts, packed_guides) for one already-oriented strand.

    5' PAM (Cas12a): guide is *after* the PAM; PAM-proximal end leads naturally.
    3' PAM (SpCas9): guide is *before* the PAM; reverse so PAM-proximal end leads.
    The protospacer is adjacent to the PAM; any "gap" is expressed by padding the
    PAM with N (e.g. ``TTTN``). Windows containing an ambiguous/invalid base are dropped.
    """
    n, pamlen = codes.shape[0], len(pam_sets)
    if side not in ("5prime", "3prime"):
        raise ValueError(f"pam_side must be '5prime' or '3prime', got {side!r}")
    pos = _pam_positions(codes, pam_sets)
    starts = (pos + pamlen) if side == "5prime" else (pos - L)
    starts = starts[(starts >= 0) & (starts + L <= n)]
    starts = _filter_valid_windows(starts, valid, L)
    return starts, _pack_windows(codes, starts, L, reverse=(side == "3prime"))


def extract_packed_guides(
    seq: str, guide_length: int, pam: str = "TTT", side: str = "5prime"
) -> np.ndarray:
    """Every PAM-anchored guide in ``seq`` (both strands), packed to ``uint64``.

    Used for commensals, where only the guide identity matters (not its position).
    """
    codes = encode(seq)
    valid = valid_mask(seq)
    pam_sets = _pam_sets(pam)
    _, fwd = _strand_guides(codes, valid, guide_length, pam_sets, side)
    _, rev = _strand_guides(reverse_complement(codes), valid[::-1], guide_length, pam_sets, side)
    return np.concatenate([fwd, rev])


def extract_target_guides(
    seq: str, guide_length: int, pam: str = "TTT", side: str = "5prime"
):
    """Guides in ``seq`` with provenance: ``(packed, forward_start, strand)``.

    ``forward_start`` is the 0-based protospacer-window start on the forward strand;
    ``strand`` is 0 (+) or 1 (-). Reverse-strand coordinates map back to forward.
    """
    codes = encode(seq)
    n = codes.shape[0]
    valid = valid_mask(seq)
    pam_sets = _pam_sets(pam)

    f_starts, f_packed = _strand_guides(codes, valid, guide_length, pam_sets, side)
    r_starts, r_packed = _strand_guides(reverse_complement(codes), valid[::-1], guide_length, pam_sets, side)
    r_fwd = n - r_starts - guide_length        # rc window coord -> forward coord

    packed = np.concatenate([f_packed, r_packed])
    starts = np.concatenate([f_starts, r_fwd]).astype(np.int64)
    strands = np.concatenate(
        [np.zeros(f_starts.size, np.uint8), np.ones(r_starts.size, np.uint8)]
    )
    return packed, starts, strands


def gc_fraction(packed: np.ndarray, guide_length: int) -> np.ndarray:
    """GC fraction of each packed guide. C=01, G=10 are the two GC codes."""
    if packed.size == 0:
        return np.empty(0, dtype=np.float64)
    gc = np.zeros(packed.shape[0], dtype=np.int32)
    for k in range(guide_length):
        field = (packed >> np.uint64(2 * k)) & np.uint64(0b11)
        gc += ((field == 1) | (field == 2)).astype(np.int32)
    return gc / guide_length
