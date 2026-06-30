# Discriminase architecture

**What it is.** Find CRISPR guide RNAs that hit a *target* bacterium while sparing a
panel of *commensal* genomes. Core operation: a fuzzy set-difference in DNA space —
keep target guides that are far (in a CRISPR-relevant sense) from every commensal
cut-site.

This document records *why* the engine is built the way it is. It supersedes the
original trie + BK-tree + pickle design, which was removed on purpose (see below).

## The redesign (2026-06-29)

The original design stored the "things to spare" as a **trie** (nested dicts) plus a
**BK-tree** (nested tuples, each node holding a pointer to a whole genome object),
serialized with `pickle`. Three fatal problems:

1. **Memory.** Both structures are millions of tiny Python objects. The trie pickle
   was 779 MB on disk → multiple GB live; the BK-tree similar. The build held all
   genomes + the full guide list + *both* trees in RAM at once (~5–10 GB peak) and
   the search `pickle.load`-ed both files fully. On a 15 GB box this OOMs and froze
   WSL. **This is what crashed the machine.**
2. **`pickle` is unsafe.** `pickle.load` of an untrusted DB = arbitrary code
   execution, and the DB is meant to be shareable.
3. **Wrong distance.** Uniform Hamming treats a PAM-proximal mismatch the same as a
   distal one; that is not how CRISPR specificity works.

### Measured on real genomes (validated 2026-06-30)

Building from two ~1.75 Mbp *C. coli* genomes (487k guides → 251,713 unique):

| | Old (trie + BK-tree, pickle) | New (sorted `uint64`, mmap) |
|---|---|---|
| Build time | — | **0.4 s** |
| Build peak RAM (Python) | ~5–10 GB → OOM | **54 MB** |
| Index on disk | hundreds of MB | **2.0 MB** |
| Load | unpickle → OOM / froze WSL | **0.001 s** (mmap) |
| Safety check (target = a commensal) | — | **0 guides kept** ✅ (every site collides at d=0) |

### The fix: a guide *is* a 64-bit integer

A guide of length L ≤ 32 packs into one `uint64` (2 bits/base). So the entire
commensal cut-site set is just a **sorted array of `uint64`** — not a tree of objects.

| | Old (trie + BK-tree, pickle) | New (sorted `uint64`, mmap) |
|---|---|---|
| On disk | ~1.27 GB | ~65 MB (8M guides × 8 B) |
| Load | unpickle 5–10 GB into RAM | `mmap` — instant, OS pages in only what's read |
| Build peak RAM | ~5–10 GB (OOM) | < ~500 MB (stream one genome at a time) |
| Security | `pickle` = code exec | plain integers |
| Provenance | genome ptr per BK node | parallel `uint8` org-id array (+8 MB) |

This is **less** code: trie + BK-tree + pickle collapse into one memmapped array + a
binary search.

### Packing convention (correctness-critical)

- Bases code as A=0, C=1, G=2, T=3 (2 bits).
- A guide is stored **PAM-proximal (seed) → distal**, with the seed in the **high**
  bits: `value = Σ code(base[j]) << (2·(L−1−j))`. So base[0] (seed start) is the most
  significant pair. **Sorting the integers groups guides by seed**, which is what
  makes the seed lookup a binary search.
- Hamming over 2-bit fields: `x = a ^ b; folded = (x | x>>1) & 0x5555…; popcount(folded)`.
  Exact; one branchless expression; vectorizes over a numpy `uint64` block.

### Distance model: seed-anchored (chosen 2026-06-29)

A target guide *g* collides with a commensal cut-site *c* iff
`hamming(seed(g), seed(c)) ≤ s` **and** `hamming(g, c) ≤ d`, where the seed is the
`seed_len` PAM-proximal bases. Rationale: CRISPR specificity is dominated by the
PAM-proximal region — a seed mismatch abolishes cutting, distal mismatches are
tolerated. This is faster and lower-memory than exact uniform Hamming (no weak
(d+1)-way pigeonhole index) **and** more biologically accurate.

Query for one *g*: enumerate seeds within `s` mismatches of `seed(g)` (for s≤1 this is
a tiny set), binary-search each seed's contiguous block in the sorted array, verify
full Hamming ≤ d on that small block, early-exit on the first collision. Exact w.r.t.
the model above (no false negatives → no commensal-cutting guide slips through).

**Defaults (tunable per nuclease, documented not dogma):** `seed_len=10`,
`seed_max_mismatch=1`, `total_max_mismatch=4`. Conservative-but-reasonable; the PAM is
Cas12a-style 5′ `TTT` so the seed is the spacer's PAM-proximal end.

### Storage format (`pickle`-free)

- `{prefix}.idx.npy` — sorted `uint64[N]`, the guides. `np.load(mmap_mode='r')`.
- `{prefix}.org.npy` — `uint8[N]` org index parallel to `.idx` (which commensal).
- `{prefix}.meta.json` — params, organism table (index → name/taxid/accession),
  counts, format version, build date.

### Build: streaming, low-memory

One genome at a time: fetch (cache the FASTA) → vectorized PAM scan + pack to `uint64`
→ append a per-genome shard to disk → free the genome. Then merge shards, sort, dedup
into `.idx`/`.org`. Peak RAM is bounded by one genome + the (tiny) packed array,
independent of panel size. Parallelism uses `os.cpu_count()` with memory-capped
concurrency (do not hold N genomes at once).

### Target selection (taxid-first)

Names are ambiguous (a species maps to thousands of assemblies); taxid/accession are
the reproducible keys. So: accept accession or taxid directly; a name search **lists
candidates** (accession · strain · assembly level · size) for the user to pick — it
never auto-takes `IdList[0]`. Custom/own sequences via FASTA upload or paste.

## Deferred (not now — revisit on a real signal)

- Position-weighted off-target scoring (CFD/MIT) instead of binary seed+Hamming.
- Indel/bulge off-targets.
- Hosted (vs local) web app.
- Compression beyond `uint64` (already ~20× smaller than the old DB).
