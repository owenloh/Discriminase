"""Discriminase: commensal-sparing CRISPR guide discovery."""
from .config import GuideFinderConfig
from .encoding import TwoBitDNA
from .trie import BitGuideTrie
from .bktree import BKTreeBitarray, hamming_distance
from .guides import build_valid_target_guides, build_pam_valid_guides
from .search import search_unique_guides

__version__ = "0.1.0"

__all__ = [
    "GuideFinderConfig",
    "TwoBitDNA",
    "BitGuideTrie",
    "BKTreeBitarray",
    "hamming_distance",
    "build_valid_target_guides",
    "build_pam_valid_guides",
    "search_unique_guides",
]
