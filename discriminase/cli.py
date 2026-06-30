"""Command-line interface.

  discriminase build          fetch a commensal panel and index it (slow, once)
  discriminase search-target  list candidate genomes for a name (pick by accession)
  discriminase find           screen a target against the index (fast, repeatable)

Target selection is accession/taxid-first on purpose: a species *name* maps to many
assemblies, so a bare name lists candidates and refuses to guess -- the accession is
the reproducible key.
"""
import argparse
import csv
import json
import os
import re
import sys
import time

from .config import GuideFinderConfig
from . import ncbi
from .database import build_from_panel_csv, build_index, make_ncbi_get_sequence
from .index import GuideIndex
from .loaders import read_sequence_file
from .search import find_sparing_guides


def _safe(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_") or "target"


def _config_from_args(args) -> GuideFinderConfig:
    cfg = GuideFinderConfig()
    for attr in ("guide_length", "db_dir", "db_prefix", "seed_len"):
        if getattr(args, attr, None) is not None:
            setattr(cfg, attr, getattr(args, attr))
    if getattr(args, "email", None):
        cfg.email = args.email
    if getattr(args, "processes", None):
        cfg.processes = args.processes
    for arg, field in [
        ("seed_mm", "seed_max_mismatch"),
        ("total_mm", "total_max_mismatch"),
        ("min_gc", "min_gc"),
        ("max_gc", "max_gc"),
        ("max_guides", "max_guides"),
    ]:
        val = getattr(args, arg, None)
        if val is not None:
            setattr(cfg, field, val)
    return cfg


# --- build ------------------------------------------------------------------
def cmd_web(args) -> int:
    """Serve the static client-side web app (no backend; runs in the browser)."""
    import functools
    import http.server
    import socketserver
    import webbrowser

    web_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "web"))
    if not os.path.isfile(os.path.join(web_dir, "index.html")):
        print(f"web app not found at {web_dir}", file=sys.stderr)
        return 1

    class Handler(http.server.SimpleHTTPRequestHandler):
        # ES modules require a JS MIME type; Python's default mislabels .mjs.
        extensions_map = {**http.server.SimpleHTTPRequestHandler.extensions_map,
                          ".mjs": "text/javascript", ".js": "text/javascript"}

        def end_headers(self):
            self.send_header("Cache-Control", "no-store")
            super().end_headers()

    handler = functools.partial(Handler, directory=web_dir)
    with socketserver.TCPServer(("", args.port), handler) as httpd:
        url = f"http://localhost:{args.port}"
        print(f"Discriminase web app: {url}   (Ctrl+C to stop)")
        try:
            webbrowser.open(url)
        except Exception:
            pass
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\nstopped")
    return 0


def cmd_build(args) -> int:
    cfg = _config_from_args(args)
    if not cfg.email:
        print("error: NCBI needs an email. Set NCBI_EMAIL or pass --email.", file=sys.stderr)
        return 2
    build_from_panel_csv(args.commensals, cfg)
    return 0


def cmd_export_web(args) -> int:
    """Build a panel and export it as a static prebuilt index for the web app."""
    cfg = _config_from_args(args)
    if not cfg.email:
        print("error: NCBI needs an email. Set NCBI_EMAIL or pass --email.", file=sys.stderr)
        return 2
    index = build_from_panel_csv(args.commensals, cfg)
    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)
    stem = f"{cfg.db_prefix}_len{cfg.guide_length}"
    index.export_web(os.path.join(out_dir, stem))

    listing_path = os.path.join(out_dir, "index.json")
    listing = []
    if os.path.exists(listing_path):
        with open(listing_path, encoding="utf-8") as fh:
            try:
                listing = json.load(fh)
            except json.JSONDecodeError:
                listing = []
    name = args.name or cfg.db_prefix.replace("_", " ")
    listing = [p for p in listing if p.get("prefix") != f"panels/{stem}"]
    listing.append({"name": f"{name} ({len(index.organisms)})", "prefix": f"panels/{stem}"})
    with open(listing_path, "w", encoding="utf-8") as fh:
        json.dump(listing, fh, indent=2)
    print(f"Exported prebuilt panel -> {out_dir}/{stem}.*  (listed in {listing_path})")
    return 0


# --- search-target ----------------------------------------------------------
def _print_candidates(cands) -> None:
    if not cands:
        print("  (no candidates found — try a broader name, or --all-assemblies)")
        return
    print(f"  {'#':>2}  {'accession':<16}  {'length':>12}  description")
    for i, c in enumerate(cands, 1):
        length = f"{c['length_bp']:,} bp" if c["length_bp"] else "?"
        print(f"  {i:>2}  {c['accession']:<16}  {length:>12}  {c['title'][:64]}")


def cmd_search_target(args) -> int:
    cfg = _config_from_args(args)
    if not cfg.email:
        print("error: NCBI needs an email. Set NCBI_EMAIL or pass --email.", file=sys.stderr)
        return 2
    ncbi.setup_entrez(cfg.email)
    cands = ncbi.search_target_candidates(
        args.query, retmax=args.limit, refseq_only=not args.all_assemblies
    )
    _print_candidates(cands)
    if cands:
        print(f"\n  Pick one:  discriminase find --target-accession {cands[0]['accession']} ...")
    return 0


# --- find -------------------------------------------------------------------
def _resolve_target(args, cfg):
    """Return (sequence, label, provenance dict). Lists candidates if a bare name."""
    if args.target_seq:
        seq = read_sequence_file(args.target_seq)
        return seq, _safe(os.path.basename(args.target_seq)), {"source": args.target_seq}
    ncbi.setup_entrez(cfg.email)
    if args.target_accession:
        seq, meta = ncbi.fetch_sequence_by_accession(args.target_accession, cfg.cache_dir)
        return seq, args.target_accession, {"accession": args.target_accession}
    if args.target_taxid:
        seq, meta = ncbi.fetch_sequence_by_taxid(args.target_taxid, cfg.cache_dir)
        return seq, f"txid{args.target_taxid}", {"taxid": args.target_taxid, **meta}
    # bare name: list candidates; only pick if the user explicitly chose an index
    cands = ncbi.search_target_candidates(args.target, retmax=args.limit)
    if args.pick:
        if not (1 <= args.pick <= len(cands)):
            print(f"error: --pick {args.pick} out of range (1..{len(cands)})", file=sys.stderr)
            return None
        chosen = cands[args.pick - 1]
        print(f"      picked #{args.pick}: {chosen['accession']}  {chosen['organism']}")
        seq, meta = ncbi.fetch_sequence_by_accession(chosen["accession"], cfg.cache_dir)
        return seq, chosen["accession"], {"accession": chosen["accession"], **chosen}
    print(f"Multiple assemblies match '{args.target}'. Pick one (accession is the key):\n")
    _print_candidates(cands)
    if cands:
        print(f"\n  Re-run:  discriminase find --target-accession {cands[0]['accession']} ...")
        print(f"  or:      discriminase find --target {args.target!r} --pick 1")
    return None


def cmd_find(args) -> int:
    if not (args.target or args.target_accession or args.target_taxid or args.target_seq):
        print("error: give one of --target / --target-accession / --target-taxid / --target-seq",
              file=sys.stderr)
        return 2

    cfg = _config_from_args(args)
    t_all = time.time()

    print("[1/3] Resolving target ...")
    resolved = _resolve_target(args, cfg)
    if resolved is None:
        return 0  # candidates listed; nothing to do until the user picks
    seq, label, provenance = resolved
    print(f"      target: {label}  ({len(seq):,} bp)")

    print("[2/3] Loading commensal index ...")
    if not GuideIndex.exists(cfg.index_prefix()):
        if args.commensals:
            print(f"      index missing; building from {args.commensals} (one-time) ...")
            build_from_panel_csv(args.commensals, cfg)
        else:
            print(f"error: no index at {cfg.index_prefix()}.idx.npy\n"
                  f"  build it:  discriminase build --commensals <panel.csv>\n"
                  f"  or pass --commensals to auto-build.", file=sys.stderr)
            return 1
    index = GuideIndex.load(cfg.index_prefix(), mmap=True)
    print(f"      {len(index):,} commensal guides from {len(index.organisms)} organisms")

    print(f"[3/3] Screening  (seed {cfg.seed_len} nt, <= {cfg.seed_max_mismatch} seed / "
          f"<= {cfg.total_max_mismatch} total mismatches)")
    t0 = time.time()
    rows = find_sparing_guides(seq, index, cfg)
    print(f"      {len(rows)} commensal-sparing guides  ({time.time() - t0:.1f}s)")

    out = args.out or os.path.join("output", f"{_safe(label)}_guides.csv")
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["position", "strand", "guide_sequence", "gc"])
        writer.writeheader()
        writer.writerows(rows)
    # reproducibility sidecar
    with open(out + ".meta.json", "w", encoding="utf-8") as fh:
        json.dump(
            {
                "target": provenance,
                "n_guides": len(rows),
                "params": {
                    "guide_length": cfg.guide_length, "pam": cfg.pam,
                    "seed_len": cfg.seed_len, "seed_max_mismatch": cfg.seed_max_mismatch,
                    "total_max_mismatch": cfg.total_max_mismatch,
                    "min_gc": cfg.min_gc, "max_gc": cfg.max_gc,
                },
                "index": cfg.index_prefix(),
            },
            fh, indent=2,
        )
    print(f"\n  {len(rows)} guides -> {out}   (total {time.time() - t_all:.1f}s)")
    return 0


# --- parser -----------------------------------------------------------------
def _add_db_args(p):
    p.add_argument("--guide-length", type=int, default=None, help="Spacer length (default 23).")
    p.add_argument("--db-dir", default=None, help="Directory holding the index.")
    p.add_argument("--db-prefix", default=None, help="Index filename prefix.")
    p.add_argument("--email", default=None, help="NCBI Entrez email (or set NCBI_EMAIL).")


def _add_model_args(p):
    p.add_argument("--seed-len", type=int, default=None, help="PAM-proximal seed length (default 10).")
    p.add_argument("--seed-mm", type=int, default=None, help="Max seed mismatches (default 1).")
    p.add_argument("--total-mm", type=int, default=None, help="Max total mismatches (default 4).")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="discriminase",
        description="Commensal-sparing CRISPR guide discovery: find guides that hit a "
                    "target while sparing a panel of commensal genomes.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    pw = sub.add_parser("web", help="Launch the point-and-click web UI.")
    pw.add_argument("--port", type=int, default=8501, help="Port (default 8501).")
    pw.set_defaults(func=cmd_web)

    pb = sub.add_parser("build", help="Build the commensal index from a panel CSV.")
    pb.add_argument("--commensals", required=True, help="CSV with headers taxid,organism_strain.")
    pb.add_argument("--processes", type=int, default=None, help="Workers (default: adapt to CPUs).")
    _add_db_args(pb)
    _add_model_args(pb)
    pb.set_defaults(func=cmd_build)

    pe = sub.add_parser("export-web", help="Build a panel and export a prebuilt index for the web app.")
    pe.add_argument("--commensals", required=True, help="CSV with headers taxid,organism_strain.")
    pe.add_argument("--out-dir", default="web/panels", help="Where to write the static index (default web/panels).")
    pe.add_argument("--name", default=None, help="Display name for the panel in the web app.")
    pe.add_argument("--processes", type=int, default=None, help="Workers (default: adapt to CPUs).")
    _add_db_args(pe)
    _add_model_args(pe)
    pe.set_defaults(func=cmd_export_web)

    ps = sub.add_parser("search-target", help="List candidate genomes for a name.")
    ps.add_argument("query", help="Organism name or search term.")
    ps.add_argument("--limit", type=int, default=25, help="Max candidates (default 25).")
    ps.add_argument("--all-assemblies", action="store_true", help="Include non-RefSeq assemblies.")
    _add_db_args(ps)
    ps.set_defaults(func=cmd_search_target)

    pf = sub.add_parser("find", help="Find commensal-sparing guides for a target.")
    g = pf.add_mutually_exclusive_group(required=True)
    g.add_argument("--target-accession", help="Target genome accession (reproducible; preferred).")
    g.add_argument("--target-taxid", help="Target NCBI taxonomy id (representative genome).")
    g.add_argument("--target-seq", help="Target sequence file (.txt/.fa/.fasta/.fna).")
    g.add_argument("--target", help="Target NAME -> lists candidates to pick from.")
    pf.add_argument("--pick", type=int, default=None, help="With --target: pick candidate #N.")
    pf.add_argument("--limit", type=int, default=25, help="Candidates to list for --target.")
    pf.add_argument("--commensals", default=None, help="Panel CSV to auto-build the index if missing.")
    pf.add_argument("--out", default=None, help="Output CSV (default output/<target>_guides.csv).")
    pf.add_argument("--min-gc", type=float, default=None, help="Min GC fraction (default 0.40).")
    pf.add_argument("--max-gc", type=float, default=None, help="Max GC fraction (default 0.60).")
    pf.add_argument("--max-guides", type=int, default=None, help="Stop after N hits (default 1000; 0=all).")
    pf.add_argument("--processes", type=int, default=None, help="Workers (default: adapt to CPUs).")
    _add_db_args(pf)
    _add_model_args(pf)
    pf.set_defaults(func=cmd_find)

    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
