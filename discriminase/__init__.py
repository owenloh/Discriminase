"""Discriminase: commensal-sparing CRISPR guide discovery."""
from .config import GuideFinderConfig
from .index import GuideIndex
from .database import build_index, build_from_panel_csv
from .search import find_sparing_guides
from .genome import extract_packed_guides, extract_target_guides
from .pack import pack_guide, unpack_guide, hamming

__version__ = "0.3.0"

__all__ = [
    "GuideFinderConfig",
    "GuideIndex",
    "build_index",
    "build_from_panel_csv",
    "find_sparing_guides",
    "extract_packed_guides",
    "extract_target_guides",
    "pack_guide",
    "unpack_guide",
    "hamming",
]
