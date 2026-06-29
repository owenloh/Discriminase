"""The two-layer commensal-sparing search.

For each candidate guide from the target:
  1. Trie seed prefilter — reject if it shares its first `sig_cutoff` bases
     exactly with any protected (commensal) guide. O(k), kills most candidates.
  2. BK-tree approximate match — reject if any protected guide is within
     `max_mismatches` Hamming distance.
What survives hits the target but is far from everything in the commensal panel.
"""
import sys
from typing import List, Tuple

from bitarray import bitarray

from .config import GuideFinderConfig
from .encoding import TwoBitDNA
from .trie import BitGuideTrie
from .bktree import BKTreeBitarray


def decode_bits(bits: bitarray) -> str:
    decode_map = TwoBitDNA._decode_map
    return ''.join(decode_map[bits[i:i + 2].to01()] for i in range(0, len(bits), 2))


def _print_progress(done: int, total: int, found: int) -> None:
    if total == 0:
        return
    bar_len = 24
    filled = int(bar_len * done / total)
    bar = "#" * filled + "-" * (bar_len - filled)
    sys.stdout.write(f"\r  screening [{bar}] {done}/{total}  ({found} kept)")
    sys.stdout.flush()


def search_unique_guides(
    target_guides: List[Tuple[int, bool, bitarray]],
    protected_trie: BitGuideTrie,
    protected_bktree: BKTreeBitarray,
    config: GuideFinderConfig,
    verbose: bool = False,
    progress: bool = True,
) -> List[Tuple[int, bool, bitarray, str]]:
    prefix_len = config.sig_cutoff
    sim_thresh = config.similarity_threshold
    max_guides = config.max_guides
    total = len(target_guides)
    unique: List[Tuple[int, bool, bitarray, str]] = []

    try:
        for i, (idx, rev, bits) in enumerate(target_guides):
            if progress and (i % 200 == 0):
                _print_progress(i, total, len(unique))

            # 1. trie seed prefilter
            if config.if_sig_cutoff and protected_trie.has_prefix(bits, prefix_len):
                continue

            # 2. BK-tree approximate match (per-guide length, robust to mixed lengths)
            guide_len = len(bits) // 2
            max_dist = int((1 - sim_thresh) * guide_len)
            if protected_bktree.search_exists(bits, max_dist):
                continue

            decoded = decode_bits(bits)
            unique.append((idx, rev, bits, decoded))
            if verbose:
                print(f"\r  found {len(unique):>4}: pos={idx} {'-' if rev else '+'} {decoded}")
            if len(unique) >= max_guides:
                break
    except KeyboardInterrupt:
        print("\nInterrupted — returning guides found so far.")
        return unique

    if progress:
        _print_progress(total, total, len(unique))
        print()
    return unique
