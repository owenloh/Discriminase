# Discriminase

**Commensal-sparing CRISPR guide discovery** — find guide RNAs that hit a target
bacterium while sparing the surrounding microbiome.

Given a target organism and a panel of commensal ("protected") genomes,
Discriminase returns guides that match the target but are far (in Hamming
distance) from everything in the panel. The search is a two-layer filter over
2-bit-packed DNA: a **trie** seed prefilter followed by a **BK-tree**
approximate-match check.

```
$ discriminase find --target "Salmonella enterica" --commensals data/commensals/gut_microbiome.csv

[1/4] Resolving + fetching target ...
      resolved: Salmonella enterica  (txid28901)
      target: Salmonella enterica  (4,857,432 bp)
[2/4] Extracting candidate guides ...
      38,402 PAM-valid, GC-filtered guides  (6.1s)
[3/4] Loading commensal database ...
[4/4] Screening vs commensals  (<= 6 mismatches)
  screening [########################] 38402/38402  (142 kept)
      142 commensal-sparing guides  (3.4s)

  142 guides -> output/Salmonella_enterica_guides.csv
```

## How it works

1. **Encode** every genome as 2 bits/base (A=00, C=01, G=10, T=11). The reverse
   complement is precomputed so both strands slice cheaply, and Hamming distance
   becomes an XOR over the packed bits.
2. **Build** (once, slow): fetch each commensal genome from NCBI, extract every
   PAM-anchored guide, and index them into a **trie** (exact seed prefixes) and a
   **BK-tree** (Hamming metric). These two files are the "things to spare."
3. **Find** (fast, repeatable): for each candidate guide in the target,
   * reject it if its first `sig_cutoff` (default 13) bases exactly match any
     commensal guide — the trie answers this in O(k) and removes most candidates;
   * otherwise reject it if any commensal guide is within `max_mismatches`
     Hamming distance — the BK-tree prunes via the triangle inequality.
   What survives hits the target and is distant from the whole commensal panel.

## Install

```bash
git clone https://github.com/owenloh/Discriminase.git
cd Discriminase
python -m venv .venv && source .venv/bin/activate
pip install -e .
```

NCBI asks callers to identify themselves. Set your email once:

```bash
export NCBI_EMAIL="you@example.com"   # or pass --email on any command
```

## Usage

### 1. Build a commensal database (one time)

```bash
discriminase build --commensals data/commensals/gut_microbiome.csv
# -> database/protected_guides_len23.trie
#    database/protected_guides_len23.bktree
#    database/protected_guides_len23.meta.json   (which organisms went in)
```

The panel CSV is just two columns:

```csv
taxid,organism_strain
511145,Escherichia coli MG1655
226186,Bacteroides thetaiotaomicron VPI-5482
...
```

### 2. Find commensal-sparing guides for a target

```bash
# by name (resolved to a taxid on NCBI)
discriminase find --target "Salmonella enterica"

# by taxid
discriminase find --target-taxid 28901

# by your own sequence file (.txt / .fa / .fasta / .fna)
discriminase find --target-seq my_target.fasta
```

Useful knobs: `--similarity 0.70` (off-target tolerance), `--min-gc/--max-gc`,
`--max-guides`, `--no-seed-filter`, `--out path.csv`. If the database does not
exist yet, pass `--commensals <panel.csv>` to `find` to auto-build it.

## The database is not in this repo (on purpose)

A built database is large (hundreds of MB to >1 GB) and fully regenerable from
the panel CSV, so it is **git-ignored**. Get one by either:

- running `discriminase build ...` yourself, or
- downloading a prebuilt artifact from the repo's **Releases** (files >100 MB are
  fine as release assets, unlike normal git objects).

> The database is serialized with `pickle`. Only load databases you built or
> trust — unpickling untrusted data can execute arbitrary code. A compact,
> safe-to-share binary format is on the roadmap.

## Speed (optional C extension)

The Hamming distance has a Cython implementation. It is optional — a pure-Python
fallback runs automatically — but ~10x faster:

```bash
pip install cython
cythonize -i discriminase/bktree_cython.pyx   # builds discriminase/bktree_cython*.so
```

## Honest notes

- **Nuclease / PAM.** The default PAM is a Cas12a/Cpf1-style 5′ `TTT` (set in
  `config.py`). It is matched exactly, so a following T is also accepted. This is
  *not* SpCas9 (which uses a 3′ `NGG`). Change `pam` / `pam_to_guide_gap` for a
  different system.
- **This is engineering, not a new algorithm.** Trie + BK-tree over packed DNA
  is a fast, clean implementation of a well-studied problem (CRISPR off-target
  search). The less-crowded angle is the *application*: multi-genome,
  commensal-sparing design rather than single-reference editing.
- **Benchmarks.** Speedups should be measured against purpose-built tools
  (Cas-OFFinder, GuideScan2, FlashFry), not BLAST — BLAST is the wrong baseline
  for fixed-length k-mer off-target search.

## Related tools

[Cas-OFFinder](https://github.com/snugel/cas-offinder) ·
[GuideScan2](https://guidescan.com) ·
[CHOPCHOP](https://chopchop.cbu.uib.no/) ·
[CRISPOR](http://crispor.tefor.net/) ·
[FlashFry](https://github.com/mckennalab/FlashFry)

## Roadmap

- 🔜 Streamlit web UI: search NCBI, click to pick target + commensal panel, run
  and download in the browser.
- ⏭️ Compact, portable, `pickle`-free database format.
- ⏭️ Honest benchmark suite vs Cas-OFFinder / GuideScan2.

## License

MIT — see [LICENSE](LICENSE).
