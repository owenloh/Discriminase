# How Discriminase compares

Honest positioning: Discriminase is **not** a new off-target algorithm and is **not**
experimentally validated. Its one real differentiator is the *framing* — designing
guides that **spare a panel of many genomes** (a microbiome), rather than avoiding
off-targets within a single reference. The engine is a clean, low-memory implementation
of a well-studied problem.

This is a **qualitative** comparison. A real head-to-head (run all tools on the same
target + panel; compare guide overlap and wall-clock) is deferred until the redesign
settles — see "What a fair benchmark must show".

## At a glance

| | **Discriminase** | Cas-OFFinder | GuideScan2 | CHOPCHOP / CRISPOR |
|---|---|---|---|---|
| Core question | guides that hit target, **spare N commensal genomes** | enumerate all off-target sites in given genome(s) | design + score guides, genome-wide | design guides for a gene, annotate off-targets |
| Multi-genome "spare" panel | **first-class** (the whole point) | possible (supply many FASTAs) but no spare/keep logic | single reference-centric | reference-centric |
| Off-target distance | seed-anchored: exact seed (±s) + Hamming ≤ d | exact, mismatches **+ DNA/RNA bulges** | mismatch enumeration + **CFD specificity score** | **CFD / MIT** scores + annotation |
| Scoring of guides | binary keep/drop + GC (no efficacy score) | none (enumerator) | validated specificity/efficiency | rich (efficiency, off-target, exonic, etc.) |
| PAM | configurable (default Cas12a 5′ `TTT`) | arbitrary | Cas9/Cas12a | several |
| Index / memory | one sorted `uint64` array, **mmap, ~65 MB for ~8M sites** | scans genomes (GPU/CPU), no persistent index | prebuilt genome index | server-side |
| Install / run | `pip install`, CLI **+ local web UI** | binary + OpenCL | CLI / web / DB | web (mostly) |
| Maturity | new, unvalidated | published, widely used | published, validated | published, widely used |

## Where each wins

- **Cas-OFFinder** — when you need *exhaustive, exact* off-target enumeration with
  bulges and an arbitrary PAM against specific genome(s). The gold standard for "find
  every near-match." No guide design or scoring; it's an enumerator.
- **GuideScan2** — when you want *validated specificity/efficiency scoring* and
  genome-wide guide design against a reference. Stronger biology than a binary cutoff.
- **CHOPCHOP / CRISPOR** — when you want a polished web tool that designs guides for a
  gene and annotates them richly against a reference genome.
- **Discriminase** — when the job is *commensal-sparing across many genomes at once*:
  give it a target and a microbiome panel, get guides predicted to hit the target and
  miss the whole panel. Plus a tiny, low-RAM, memmapped index and a click-to-use UI.

## Honest gaps (vs. the above)

1. **No position-weighted scoring.** Seed-anchored binary cutoff, not CFD/MIT. A
   distal mismatch and a near-seed mismatch outside the seed window are treated the
   same. Competitors model this; Discriminase deliberately approximates.
2. **No bulges/indels.** Substitution mismatches only. Cas-OFFinder handles bulges.
3. **Hamming is a model, not truth.** Real cutting depends on chromatin, expression,
   sequence context. None of this is modeled.
4. **Unvalidated.** No wet-lab confirmation; competitors are published and benchmarked.

## What a fair benchmark must show (the deferred real test)

- Same target + same panel for every tool (Discriminase panel = the FASTAs you give
  the others).
- **Guide overlap**: of Discriminase's "spared" guides, how many does Cas-OFFinder
  confirm have no panel hit within the tolerance? (Tests the no-false-negative claim
  on real data.)
- **Speed**: wall-clock build + query, like-for-like tolerance.
- **NOT** "2000× faster than BLAST" — BLAST is the wrong baseline for fixed-length
  k-mer off-target search. Compare to the tools above or don't claim a speedup.
