# Discriminase

**Commensal-sparing CRISPR guide design** — find guide RNAs that hit *one or more*
target bacteria while sparing every commensal around them.

🔗 **Live app (no install): https://discriminase-production.up.railway.app** — runs
entirely in your browser.

Use it to design **narrow-spectrum, sequence-specific** CRISPR tools: kill or edit a
single pathogen or strain without touching the rest of a microbial community. Give it a
target organism (or several) and a panel of commensal ("protected") genomes, and
Discriminase returns guides that match the target but stay clear — in a CRISPR-relevant,
seed-anchored sense — of every genome in the panel. With several targets, the web app
keeps only the guides shared by **every** target (exact match) before screening, so you
get guides that cut your whole target set and spare the panel.

**Who it's for:** iGEM teams and researchers building CRISPR antimicrobials,
microbiome-editing tools, strain-selective knockouts, or any guide that must hit one or
more organisms and spare their neighbours.

```
$ discriminase find --target-accession NZ_CP039503.1 \
                    --commensals data/commensals/gut_microbiome.csv

[1/3] Resolving target ...
      target: NZ_CP039503.1  (4,882,218 bp)
[2/3] Loading commensal index ...
      8,142,377 commensal guides from 52 organisms
[3/3] Screening  (seed 10 nt, <= 1 seed / <= 4 total mismatches)
  screening [########################] 38402/38402  (142 kept)
      142 commensal-sparing guides  (3.4s)

  142 guides -> output/NZ_CP039503_1_guides.csv
```

Prefer clicking to typing? Open the [live app](https://discriminase-production.up.railway.app),
or run it locally with `discriminase web`.

## How it works

A guide of up to 31 nt **is a single 64-bit integer** (2 bits/base, A=0 C=1 G=2 T=3),
stored PAM-proximal (seed) first. That one idea drives everything:

1. **Encode** each genome as 2-bit codes and slide a PAM-anchored window over both
   strands to pull out every candidate guide, packed to a `uint64`.
2. **Build** (once): stream the commensal panel one genome at a time, writing packed
   guides to disk, then merge into a single **sorted `uint64` array** — the index. It
   is memory-mapped, so loading is instant and uses almost no RAM.
3. **Find** (fast, repeatable): for each target guide, binary-search the index for
   commensals that share its **seed** (the PAM-proximal bases), then verify the full
   Hamming distance on that small block. A guide is kept only if **no** commensal is
   within the seed-anchored tolerance.

The distance model is **seed-anchored** because CRISPR specificity is dominated by the
PAM-proximal region: a guide collides with a commensal when their seeds match within
`seed_max_mismatch` *and* the full guides match within `total_max_mismatch`.

New to this? [docs/HOW_IT_WORKS.md](docs/HOW_IT_WORKS.md) is a no-skipping, visual
walkthrough from raw DNA to the keep/drop decision. Design rationale + measured
numbers: [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md).

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

## Use it

### Web app (easiest — runs in the browser, no backend)

The UI in `web/` is a **static client-side app**: all computation (genome download,
indexing, search) happens in the user's browser. Serve it locally with:

```bash
discriminase web        # serves web/ at http://localhost:8501
```

Or host the `web/` folder anywhere static (Netlify, GitHub Pages, Railway-static) and
share the link — there is no server to run. See [Host it](#host-it-no-backend).

In the app: set the nuclease/parameters (presets for Cas12a & SpCas9, or a custom PAM)
→ build a commensal panel (search NCBI and add, upload FASTA, or load a prebuilt panel)
→ pick a target (search-and-click, accession, paste, or upload) → **Run** → download
the CSV. Panels and parameters are saved in the browser, so a reload won't lose them.

### Command line

**1 · Build a commensal index (one time, slow — downloads genomes):**

```bash
discriminase build --commensals data/commensals/gut_microbiome.csv
# -> database/gut_microbiome_len23.idx.npy   (sorted uint64 guides)
#    database/gut_microbiome_len23.org.npy   (which organism each came from)
#    database/gut_microbiome_len23.meta.json (panel + parameters)
```

Panel CSV is two columns:

```csv
taxid,organism_strain
511145,Escherichia coli MG1655
226186,Bacteroides thetaiotaomicron VPI-5482
```

**2 · Pick a target.** A species name maps to many assemblies, so a name **lists
candidates** and never silently guesses — the accession is the reproducible key:

```bash
discriminase search-target "Salmonella enterica"
#  #  accession         length        description
#  1  NZ_CP191545.1     4,898,209 bp  Salmonella enterica strain PP14-31 chromosome...
#  ...
```

**3 · Find sparing guides** (pick one input — accession is preferred):

```bash
discriminase find --target-accession NZ_CP191545.1     # reproducible
discriminase find --target-taxid 28901                 # representative genome
discriminase find --target-seq my_target.fasta         # your own sequence
discriminase find --target "Salmonella enterica" --pick 1   # list + pick #N
```

Knobs: `--total-mm` (off-target tolerance), `--seed-len`, `--seed-mm`,
`--min-gc/--max-gc`, `--max-guides`, `--out`. Pass `--commensals <panel.csv>` to
`find` to auto-build the index if it is missing.

## The index is not in this repo (on purpose)

A built index is regenerable from the panel CSV, so it is **git-ignored**. Build it
yourself, or download a prebuilt artifact from **Releases** (large files are fine as
release assets). The format is plain NumPy arrays — **no `pickle`**, so loading a
shared index can't execute code.

## Host it (no backend)

The whole app is static files in `web/`. To put it online for others:

- **Netlify / GitHub Pages / Cloudflare Pages:** publish the `web/` folder. Done —
  share the URL. NCBI's genome API allows browser requests (CORS is open), so users'
  browsers fetch genomes directly; nothing runs on your server.
- **Railway:** the included `Dockerfile` + `Caddyfile` serve `web/` on Railway's
  `$PORT`. Step-by-step in [docs/DEPLOY.md](docs/DEPLOY.md).
- **Prebuilt panels (optional):** so users don't have to build a panel themselves,
  generate one offline and ship it as a static file:

  ```bash
  discriminase export-web --commensals data/commensals/gut_microbiome.csv \
                          --name "Gut microbiome" --out-dir web/panels
  ```

  It appears in the app's "prebuilt panel" dropdown. Users can still build their own
  panel (search NCBI in-browser) or upload FASTAs.

## Honest notes

- **Nuclease / PAM.** Default is a Cas12a/Cpf1-style 5′ `TTT` PAM (in `config.py`),
  matched exactly. *Not* SpCas9 (3′ `NGG`). Change `pam` / `pam_side` for another
  system; the guide sits right next to the PAM, so any spacing goes inside the PAM
  as `N` (e.g. `TTTN`).
- **Engineering, not a new algorithm.** Packed k-mers + a sorted index is a clean,
  fast take on a well-studied problem. The less-crowded angle is the *application*:
  multi-genome, commensal-sparing design.
- **Model, not truth.** Seed + Hamming is a deliberate approximation — no
  position-weighted (CFD/MIT) scoring, no bulges, no wet-lab validation. See how it
  stacks up against Cas-OFFinder / GuideScan2 / CRISPOR in
  [docs/COMPARISON.md](docs/COMPARISON.md).

## Related tools

[Cas-OFFinder](https://github.com/snugel/cas-offinder) ·
[GuideScan2](https://guidescan.com) ·
[CHOPCHOP](https://chopchop.cbu.uib.no/) ·
[CRISPOR](http://crispor.tefor.net/) ·
[FlashFry](https://github.com/mckennalab/FlashFry)

## License

MIT — see [LICENSE](LICENSE).
