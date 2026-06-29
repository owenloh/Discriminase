"""Command-line interface: `discriminase build` and `discriminase find`.

Two phases, deliberately split:
  build  — fetch a commensal panel and index it (slow, once).
  find   — screen a target against that index (fast, repeatable).
"""
import argparse
import csv
import os
import re
import sys
import time

from .config import GuideFinderConfig
from . import ncbi
from .guides import build_valid_target_guides
from .search import search_unique_guides
from .trie import BitGuideTrie
from .bktree import BKTreeBitarray
from .loaders import read_sequence_file
from .database import build_database
from .encoding import TwoBitDNA


def _safe(name: str) -> str:
    return re.sub(r'[^A-Za-z0-9]+', '_', name).strip('_')


def _config_from_args(args) -> GuideFinderConfig:
    config = GuideFinderConfig()
    config.guide_length = args.guide_length
    config.db_dir = args.db_dir
    config.db_prefix = args.db_prefix
    if args.email:
        config.email = args.email
    if getattr(args, "processes", None):
        config.processes = args.processes
    # find-only knobs (absent on build)
    for attr, field in [
        ("similarity", "similarity_threshold"),
        ("min_gc", "min_gc"),
        ("max_gc", "max_gc"),
        ("max_guides", "max_guides"),
    ]:
        val = getattr(args, attr, None)
        if val is not None:
            setattr(config, field, val)
    if getattr(args, "no_seed_filter", False):
        config.if_sig_cutoff = False
    return config


def cmd_build(args) -> int:
    config = _config_from_args(args)
    build_database(args.commensals, config)
    return 0


def _resolve_target_genome(args) -> TwoBitDNA:
    if args.target_seq:
        seq = read_sequence_file(args.target_seq)
        print(f"      loaded sequence ({len(seq):,} bp) from {args.target_seq}")
        return TwoBitDNA("target_sequence", seq)
    if args.target_taxid:
        print(f"      fetching genome for txid{args.target_taxid} ...")
        return ncbi.fetch_genome_by_taxid(args.target_taxid)
    hit = ncbi.resolve_taxid(args.target)
    if hit:
        taxid, sci = hit
        print(f"      resolved: {sci}  (txid{taxid})")
        return ncbi.fetch_genome_by_taxid(taxid, sci)
    print(f"      no taxid match; trying direct genome search for '{args.target}'")
    return ncbi.fetch_genome_by_name(args.target)


def _ensure_database(config: GuideFinderConfig, commensals: str) -> None:
    if os.path.exists(config.trie_path()) and os.path.exists(config.bktree_path()):
        return
    if commensals:
        print(f"      database missing; building from {commensals} (one-time, slow) ...")
        build_database(commensals, config)
        return
    raise FileNotFoundError(
        f"No database at {config.trie_path()} / {config.bktree_path()}.\n"
        f"  Build it:   discriminase build --commensals <panel.csv>\n"
        f"  or pass --commensals to auto-build, or download a prebuilt release."
    )


def cmd_find(args) -> int:
    if not (args.target or args.target_taxid or args.target_seq):
        print("error: provide one of --target / --target-taxid / --target-seq", file=sys.stderr)
        return 2

    config = _config_from_args(args)
    t_all = time.time()

    print("[1/4] Resolving + fetching target ...")
    if args.target or args.target_taxid:
        ncbi.setup_entrez(config.email)
    genome = _resolve_target_genome(args)
    print(f"      target: {genome.name}  ({len(genome):,} bp)")

    print("[2/4] Extracting candidate guides ...")
    t0 = time.time()
    targets = build_valid_target_guides(genome, config)
    targets = targets[config.subsection[0]:config.subsection[1]]
    print(f"      {len(targets):,} PAM-valid, GC-filtered guides  ({time.time() - t0:.1f}s)")
    if not targets:
        print("      no candidate guides; nothing to screen.")
        return 0

    print("[3/4] Loading commensal database ...")
    _ensure_database(config, args.commensals)
    trie = BitGuideTrie.load(config.trie_path())
    bktree = BKTreeBitarray.load(config.bktree_path())

    print(f"[4/4] Screening vs commensals  (<= {config.max_mismatches} mismatches)")
    t0 = time.time()
    unique = search_unique_guides(targets, trie, bktree, config, verbose=args.verbose)
    print(f"      {len(unique)} commensal-sparing guides  ({time.time() - t0:.1f}s)")

    out = args.out or os.path.join("output", f"{_safe(genome.name)}_guides.csv")
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    with open(out, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["position", "strand", "guide_sequence"])
        for idx, rev, _bits, seq in unique:
            writer.writerow([idx, "-" if rev else "+", seq])

    print(f"\n  {len(unique)} guides -> {out}   (total {time.time() - t_all:.1f}s)")
    return 0


def _add_db_args(p):
    p.add_argument("--guide-length", type=int, default=23, help="Spacer length (default 23).")
    p.add_argument("--db-dir", default="database", help="Directory holding the trie/BK-tree.")
    p.add_argument("--db-prefix", default="protected_guides", help="Database filename prefix.")
    p.add_argument("--email", default=None, help="NCBI Entrez email (or set NCBI_EMAIL).")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="discriminase",
        description="Commensal-sparing CRISPR guide discovery: find guides that hit a "
                    "target bacterium while sparing a panel of commensal genomes.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    pb = sub.add_parser("build", help="Build the commensal (protected) database from a panel CSV.")
    pb.add_argument("--commensals", required=True, help="CSV with headers taxid,organism_strain.")
    pb.add_argument("--processes", type=int, default=4, help="Worker processes for genome scan.")
    _add_db_args(pb)
    pb.set_defaults(func=cmd_build)

    pf = sub.add_parser("find", help="Find commensal-sparing guides for a target.")
    g = pf.add_mutually_exclusive_group(required=True)
    g.add_argument("--target", help="Target organism NAME (resolved on NCBI).")
    g.add_argument("--target-taxid", help="Target NCBI taxonomy id.")
    g.add_argument("--target-seq", help="Target sequence file (.txt/.fa/.fasta/.fna).")
    pf.add_argument("--commensals", default=None, help="Panel CSV to auto-build the DB if missing.")
    pf.add_argument("--out", default=None, help="Output CSV path (default output/<target>_guides.csv).")
    pf.add_argument("--similarity", type=float, default=None, help="Similarity threshold (default 0.70).")
    pf.add_argument("--min-gc", type=float, default=None, help="Min GC fraction (default 0.40).")
    pf.add_argument("--max-gc", type=float, default=None, help="Max GC fraction (default 0.60).")
    pf.add_argument("--max-guides", type=int, default=None, help="Stop after this many hits (default 1000).")
    pf.add_argument("--no-seed-filter", action="store_true", help="Disable the trie seed prefilter.")
    pf.add_argument("--processes", type=int, default=4, help="Worker processes for genome scan.")
    pf.add_argument("-v", "--verbose", action="store_true", help="Print each guide as it is found.")
    _add_db_args(pf)
    pf.set_defaults(func=cmd_find)

    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
