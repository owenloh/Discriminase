"""The commensal cut-site index: one sorted ``uint64`` array, memory-mapped.

This single structure replaces the old trie + BK-tree. Guides are stored
PAM-proximal-first (see :mod:`pack`), so the array sorted ascending groups guides by
seed and the seed lookup is a binary search over a contiguous block.

A target guide *g* "collides" with a commensal cut-site *c* when
``hamming(seed(g), seed(c)) <= seed_max_mismatch`` **and** ``hamming(g, c) <= d``.
:meth:`GuideIndex.query` answers that exactly: enumerate the (small) set of seeds
within ``seed_max_mismatch`` of *g*'s seed, binary-search each seed's block, verify
the full Hamming on that block, and return the nearest collision (or ``None`` if *g*
is commensal-sparing). No false negatives -> no commensal-cutting guide slips through.

On disk:
  {prefix}.idx.npy   sorted uint64[N]  -- the guides
  {prefix}.org.npy   uint8[N]          -- which commensal each guide came from
  {prefix}.meta.json                   -- params + organism table + counts
"""
import json
import os
from typing import List, Optional, Tuple

import numpy as np

from .pack import hamming_array, seed_neighbors, seed_of

FORMAT_VERSION = 1


class GuideIndex:
    def __init__(
        self,
        guides: np.ndarray,
        guide_length: int,
        seed_len: int,
        org_ids: Optional[np.ndarray] = None,
        organisms: Optional[List[dict]] = None,
        meta: Optional[dict] = None,
    ):
        if guide_length > 31:
            raise ValueError("guide_length must be <= 31 (uint64 seed-bound headroom)")
        if seed_len > guide_length:
            raise ValueError("seed_len cannot exceed guide_length")
        self.guides = guides                      # sorted uint64[N]
        self.org_ids = org_ids                    # uint8[N] or None
        self.guide_length = guide_length
        self.seed_len = seed_len
        self.organisms = organisms or []          # index -> {name, taxid, accession}
        self.meta = meta or {}
        self._shift = 2 * (guide_length - seed_len)

    def __len__(self) -> int:
        return int(self.guides.shape[0])

    # --- construction -------------------------------------------------------
    @classmethod
    def from_packed(
        cls,
        guides: np.ndarray,
        guide_length: int,
        seed_len: int,
        org_ids: Optional[np.ndarray] = None,
        organisms: Optional[List[dict]] = None,
        meta: Optional[dict] = None,
    ) -> "GuideIndex":
        """Sort + dedup packed guides into an index (keeps one org per unique guide)."""
        guides = np.asarray(guides, dtype=np.uint64)
        order = np.argsort(guides, kind="stable")
        guides = guides[order]
        if org_ids is not None:
            org_ids = np.asarray(org_ids, dtype=np.uint16)[order]

        if guides.size:
            keep = np.empty(guides.size, dtype=bool)
            keep[0] = True
            np.not_equal(guides[1:], guides[:-1], out=keep[1:])
            guides = guides[keep]
            if org_ids is not None:
                org_ids = org_ids[keep]

        return cls(guides, guide_length, seed_len, org_ids, organisms, meta)

    # --- query --------------------------------------------------------------
    def query(
        self, g: int, total_max_mismatch: int, seed_max_mismatch: int
    ) -> Optional[Tuple[int, int]]:
        """Return ``(distance, org_id)`` of the nearest colliding commensal, or None.

        ``None`` means *g* is commensal-sparing under the model (safe to keep).
        ``org_id`` is -1 when no provenance array is loaded.
        """
        seeds = seed_neighbors(
            seed_of(g, self.guide_length, self.seed_len),
            self.seed_len,
            seed_max_mismatch,
        )
        guides = self.guides
        best: Optional[Tuple[int, int]] = None
        for s in seeds:
            lo_key = np.uint64(s << self._shift)
            hi_key = np.uint64((s + 1) << self._shift)
            lo = int(np.searchsorted(guides, lo_key, side="left"))
            hi = int(np.searchsorted(guides, hi_key, side="left"))
            if lo == hi:
                continue
            dists = hamming_array(g, guides[lo:hi])
            within = np.nonzero(dists <= total_max_mismatch)[0]
            if within.size == 0:
                continue
            j = int(within[np.argmin(dists[within])])
            d = int(dists[j])
            if best is None or d < best[0]:
                org = int(self.org_ids[lo + j]) if self.org_ids is not None else -1
                best = (d, org)
                if d == 0:
                    return best  # cannot do better
        return best

    def is_spared(self, g: int, total_max_mismatch: int, seed_max_mismatch: int) -> bool:
        """True if *g* hits nothing in the panel (keep it)."""
        return self.query(g, total_max_mismatch, seed_max_mismatch) is None

    # --- persistence --------------------------------------------------------
    @staticmethod
    def paths(prefix: str):
        return f"{prefix}.idx.npy", f"{prefix}.org.npy", f"{prefix}.meta.json"

    def save(self, prefix: str) -> None:
        idx_path, org_path, meta_path = self.paths(prefix)
        os.makedirs(os.path.dirname(prefix) or ".", exist_ok=True)
        np.save(idx_path, self.guides)
        if self.org_ids is not None:
            np.save(org_path, self.org_ids)
        meta = dict(self.meta)
        meta.update(
            {
                "format_version": FORMAT_VERSION,
                "guide_length": self.guide_length,
                "seed_len": self.seed_len,
                "n_guides": len(self),
                "has_org_ids": self.org_ids is not None,
                "organisms": self.organisms,
            }
        )
        with open(meta_path, "w", encoding="utf-8") as fh:
            json.dump(meta, fh, indent=2)

    @classmethod
    def load(cls, prefix: str, mmap: bool = True) -> "GuideIndex":
        idx_path, org_path, meta_path = cls.paths(prefix)
        if not os.path.exists(idx_path):
            raise FileNotFoundError(f"no index at {idx_path}")
        mode = "r" if mmap else None
        guides = np.load(idx_path, mmap_mode=mode)
        org_ids = np.load(org_path, mmap_mode=mode) if os.path.exists(org_path) else None
        meta = {}
        if os.path.exists(meta_path):
            with open(meta_path, encoding="utf-8") as fh:
                meta = json.load(fh)
        return cls(
            guides,
            int(meta.get("guide_length", 0)) or _infer_len(guides),
            int(meta.get("seed_len", 10)),
            org_ids,
            meta.get("organisms", []),
            meta,
        )

    @staticmethod
    def exists(prefix: str) -> bool:
        return os.path.exists(f"{prefix}.idx.npy")

    # --- browser-loadable export -------------------------------------------
    def export_web(self, prefix: str) -> dict:
        """Write a browser-friendly index: raw little-endian binaries + manifest.

        Guides are stored as float64 (a 23-mer is <= 2^46, exact in float64, and
        JavaScript Numbers are float64 -- so the static site reads the file straight
        into a ``Float64Array`` with no BigInt or uint64 gymnastics). Org ids go to a
        parallel uint16 file. ``{prefix}.web.json`` carries the parameters.
        """
        os.makedirs(os.path.dirname(prefix) or ".", exist_ok=True)
        np.asarray(self.guides, dtype="<f8").tofile(f"{prefix}.guides.f64")
        if self.org_ids is not None:
            np.asarray(self.org_ids, dtype="<u2").tofile(f"{prefix}.orgs.u16")
        manifest = dict(self.meta)
        manifest.update(
            {
                "format_version": FORMAT_VERSION,
                "guide_length": self.guide_length,
                "seed_len": self.seed_len,
                "n_guides": len(self),
                "has_org_ids": self.org_ids is not None,
                "organisms": self.organisms,
                "guides_file": os.path.basename(f"{prefix}.guides.f64"),
                "orgs_file": os.path.basename(f"{prefix}.orgs.u16") if self.org_ids is not None else None,
                "encoding": "float64-le; guide = sum(code(base[j]) * 4**(L-1-j)), seed (PAM-proximal) first",
            }
        )
        with open(f"{prefix}.web.json", "w", encoding="utf-8") as fh:
            json.dump(manifest, fh, indent=2)
        return manifest


def _infer_len(guides: np.ndarray) -> int:
    """Best-effort guide length from the data if meta is missing."""
    if guides.size == 0:
        return 0
    top = int(guides.max())
    return max(1, (top.bit_length() + 1) // 2)
