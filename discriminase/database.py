"""Build the commensal cut-site index, streaming one genome at a time.

Peak memory is bounded by a single genome plus its (tiny) packed guides -- never the
whole panel -- which is what keeps the build off the OOM cliff the old trie/BK-tree
build fell off. Each genome's guides are written to a temporary shard on disk; only at
the end are the shards merged, sorted and de-duplicated into the final memmappable
index (see :mod:`index`).

``get_sequence`` is injected so the same builder serves local FASTAs (tests, custom
panels) and NCBI downloads (real commensal panels) without caring which.
"""
import os
import tempfile
import time
from datetime import datetime, timezone
from typing import Callable, List, Optional

import numpy as np

from .config import GuideFinderConfig
from .genome import extract_packed_guides
from .index import GuideIndex
from .loaders import read_sequence_file, read_taxid_csv

# A source describes one commensal: {"name", and one of "fasta"/"taxid"/"accession"}.
Source = dict
GetSequence = Callable[[Source], "tuple[str, dict]"]


def build_index(
    sources: List[Source],
    config: GuideFinderConfig,
    get_sequence: GetSequence,
    progress: Callable[[str], None] = print,
) -> GuideIndex:
    """Stream genomes -> shards -> one merged, sorted, de-duplicated GuideIndex."""
    os.makedirs(config.db_dir, exist_ok=True)
    organisms: List[dict] = []
    failed: List[dict] = []

    t0 = time.time()
    with tempfile.TemporaryDirectory(dir=config.db_dir) as tmp:
        shards: List[tuple] = []          # (guides_path, org_index, count)
        for src in sources:
            name = src.get("name", "?")
            try:
                seq, meta = get_sequence(src)
            except Exception as e:        # noqa: BLE001 - report and continue
                progress(f"  ! {name}: {e}")
                failed.append({**src, "error": str(e)})
                continue

            packed = extract_packed_guides(
                seq, config.guide_length, config.pam,
                config.pam_to_guide_gap, config.pam_side,
            )
            org_index = len(organisms)
            organisms.append(
                {
                    "name": name,
                    "taxid": meta.get("taxid"),
                    "accession": meta.get("accession"),
                    "length_bp": len(seq),
                    "n_guides": int(packed.size),
                }
            )
            shard = os.path.join(tmp, f"shard_{org_index}.npy")
            np.save(shard, packed)
            shards.append((shard, org_index, int(packed.size)))
            progress(
                f"  [{org_index + 1}/{len(sources)}] {name}: "
                f"{len(seq):,} bp -> {packed.size:,} guides"
            )
            del seq, packed               # free before the next genome

        if not shards:
            raise RuntimeError("no genomes ingested; cannot build index")

        total = sum(c for _, _, c in shards)
        progress(f"Merging {len(shards)} shards ({total:,} guides) ...")
        all_guides = np.empty(total, dtype=np.uint64)
        all_orgs = np.empty(total, dtype=np.uint16)
        at = 0
        for shard, org_index, count in shards:
            all_guides[at : at + count] = np.load(shard)
            all_orgs[at : at + count] = org_index
            at += count

    meta = {
        "pam": config.pam,
        "pam_side": config.pam_side,
        "pam_to_guide_gap": config.pam_to_guide_gap,
        "seed_max_mismatch": config.seed_max_mismatch,
        "total_max_mismatch": config.total_max_mismatch,
        "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S"),
        "failed": failed,
    }
    index = GuideIndex.from_packed(
        all_guides,
        config.guide_length,
        config.seed_len,
        org_ids=all_orgs,
        organisms=organisms,
        meta=meta,
    )
    index.save(config.index_prefix())
    progress(
        f"Index: {len(index):,} unique guides from {len(organisms)} organisms "
        f"-> {config.index_prefix()}.idx.npy  ({time.time() - t0:.1f}s)"
    )
    return index


# --- get_sequence implementations ------------------------------------------
def local_get_sequence(src: Source) -> "tuple[str, dict]":
    """Read a commensal genome from a local FASTA path (no network)."""
    seq = read_sequence_file(src["fasta"])
    return seq, {"accession": src.get("accession"), "taxid": src.get("taxid")}


def make_ncbi_get_sequence(config: GuideFinderConfig) -> GetSequence:
    """Fetch (and cache) a representative genome per commensal from NCBI."""
    from . import ncbi

    ncbi.setup_entrez(config.email)

    def _get(src: Source) -> "tuple[str, dict]":
        if src.get("fasta"):
            return local_get_sequence(src)
        if src.get("accession"):
            return ncbi.fetch_sequence_by_accession(src["accession"], config.cache_dir)
        if src.get("taxid"):
            return ncbi.fetch_sequence_by_taxid(src["taxid"], config.cache_dir, src.get("name"))
        raise ValueError(f"source has no fasta/accession/taxid: {src}")

    return _get


def build_from_panel_csv(csv_path: str, config: GuideFinderConfig,
                         progress: Callable[[str], None] = print) -> GuideIndex:
    """Build from a panel CSV (headers: taxid, organism_strain) via NCBI."""
    organisms = read_taxid_csv(csv_path)
    sources = [{"name": name, "taxid": taxid} for name, taxid in organisms]
    progress(f"Loaded {len(sources)} commensals from {csv_path}")
    return build_index(sources, config, make_ncbi_get_sequence(config), progress)
