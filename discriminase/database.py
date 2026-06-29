"""Build the protected-guide database (trie + BK-tree) from a commensal panel.

Fetches each commensal genome from NCBI, extracts every PAM-anchored guide,
and indexes them into a trie (seed prefix) and a BK-tree (Hamming). These two
files ARE the "things to spare" — the search rejects any target guide that
collides with them. The database is regenerable, so it is never committed; see
the README for how to (re)build or download it.
"""
import json
import os
import time
from datetime import datetime
from typing import List, Tuple

from .config import GuideFinderConfig
from .ncbi import setup_entrez, fetch_genome_by_taxid
from .guides import build_pam_valid_guides
from .loaders import read_taxid_csv
from .trie import BitGuideTrie
from .bktree import BKTreeBitarray, hamming_distance


def _build_trie(guides, guide_length: int) -> BitGuideTrie:
    trie = BitGuideTrie()
    for genome, pos, is_rev in guides:
        bits = genome.get_encoded_slice(pos, guide_length, is_rev)
        if len(bits) == 2 * guide_length:
            trie.insert(bits)
    return trie


def _build_bktree(guides, guide_length: int) -> BKTreeBitarray:
    bktree = BKTreeBitarray(hamming_distance)
    for genome, pos, is_rev in guides:
        bits = genome.get_encoded_slice(pos, guide_length, is_rev)
        if len(bits) == 2 * guide_length:
            bktree.insert((bits, genome, pos, is_rev))
    return bktree


def build_database(commensals_csv: str, config: GuideFinderConfig) -> dict:
    """Fetch commensal genomes and build + save the trie and BK-tree.

    Returns a summary dict with file paths and per-organism success/failure.
    """
    setup_entrez(config.email)
    os.makedirs(config.db_dir, exist_ok=True)

    organisms = read_taxid_csv(commensals_csv)
    print(f"Loaded {len(organisms)} commensal organisms from {commensals_csv}")

    genomes = []
    succeeded: List[Tuple[str, str]] = []
    failed: List[Tuple[str, str, str]] = []
    t0 = time.time()
    for name, taxid in organisms:
        try:
            print(f"  fetching {name} (txid{taxid}) ...", flush=True)
            genomes.append(fetch_genome_by_taxid(taxid, name))
            succeeded.append((name, taxid))
        except Exception as e:  # noqa: BLE001 - report and continue
            print(f"    ! failed: {e}")
            failed.append((name, taxid, str(e)))
    print(f"Fetched {len(genomes)}/{len(organisms)} genomes in {time.time() - t0:.1f}s")

    if not genomes:
        raise RuntimeError("No genomes fetched; cannot build database.")

    print("Extracting protected guides ...")
    t0 = time.time()
    guides = build_pam_valid_guides(genomes, config)
    print(f"  {len(guides)} guides in {time.time() - t0:.1f}s")

    print("Building trie ...")
    t0 = time.time()
    trie = _build_trie(guides, config.guide_length)
    trie.save(config.trie_path())
    print(f"  {trie.count_guides()} unique guides -> {config.trie_path()} ({time.time() - t0:.1f}s)")

    print("Building BK-tree ...")
    t0 = time.time()
    bktree = _build_bktree(guides, config.guide_length)
    bktree.save(config.bktree_path())
    print(f"  {bktree.count_guides()} guides -> {config.bktree_path()} ({time.time() - t0:.1f}s)")

    summary = {
        "guide_length": config.guide_length,
        "created": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "trie_file": config.trie_path(),
        "bktree_file": config.bktree_path(),
        "commensals_successful": [{"organism_strain": n, "taxid": t} for n, t in succeeded],
        "commensals_failed": [{"organism_strain": n, "taxid": t, "error": e} for n, t, e in failed],
    }
    meta_path = os.path.join(config.db_dir, f"{config.db_prefix}_len{config.guide_length}.meta.json")
    with open(meta_path, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    print(f"Wrote build manifest -> {meta_path}")
    return summary
