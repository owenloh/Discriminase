# Deploying the web app on Railway

The app is **static files in `web/`** — there is no application server and no database.
Users' browsers download genomes (from NCBI) and do all the computation. Railway just
serves the files. This repo ships the two files that make that work:

- `Dockerfile` — a tiny [Caddy](https://caddyserver.com) image that serves `web/`.
- `Caddyfile` — listens on Railway's `$PORT`, gzips, and sets the correct
  `text/javascript` MIME type for the `.mjs` modules (browsers refuse modules served
  as anything else — this is the #1 thing that breaks naive static hosts).

## Deploy from the dashboard (~3 minutes)

1. Make sure the repo is on GitHub (it is: `owenloh/Discriminase`, branch `main`).
2. Go to **railway.app → New Project → Deploy from GitHub repo** and pick
   `owenloh/Discriminase`.
3. Railway detects the `Dockerfile` and builds it. **No environment variables needed**
   — `PORT` is injected automatically.
4. After the build goes green: **Settings → Networking → Generate Domain**. That URL is
   your public app. Share it.
5. (Optional) **Settings → Networking → Custom Domain** to use your own domain.

## Deploy from the CLI (alternative)

```bash
npm i -g @railway/cli
railway login
railway init          # run inside the repo; creates a project
railway up            # builds the Dockerfile and deploys
railway domain        # prints/creates a public URL
```

## Updating

Push to `main` → Railway auto-redeploys. The `Cache-Control: no-cache` header in the
`Caddyfile` means users get the new app code on their next load (no hard refresh).

## Shipping a prebuilt panel with the hosted app (optional)

Prebuilt panel binaries are git-ignored (they're large and regenerable). Without one,
users still build their own panel in-browser or upload FASTAs. To include a ready-made
panel so it shows up in the app's "prebuilt panel" dropdown:

```bash
discriminase export-web --commensals data/commensals/gut_microbiome.csv \
                        --name "Gut microbiome"
git add -f web/panels/gut_microbiome_len23.* web/panels/index.json
git commit -m "Ship prebuilt gut microbiome panel" && git push
```

(`-f` overrides the `.gitignore`. A ~52-organism panel is ~tens of MB — fine as a
normal file, well under limits, but don't commit hundreds of them.)

## Usage analytics ("how many people / what for")

The static app has no backend, so analytics is a deliberate add-on. Two routes — pick
in the chat and I'll wire it:

- **Self-hosted Umami on Railway (recommended).** Add Umami (open-source,
  privacy-friendly) as a second Railway service with a Postgres plugin — Railway has a
  one-click template. The app sends a tiny `umami.track("run", {...})` on each run/build.
  You get visitor counts **and** custom events ("what for"), you own the data, no custom
  backend code to maintain.
- **Custom backend.** A small server (serve `web/` + `POST /api/event` → Postgres +
  a protected `/api/stats`). Full control, more code and maintenance.

Either way, log **coarse, non-identifying** events by default (nuclease used, panel
size, whether the target came from a name/accession/upload, how many guides were found)
— not the actual sequences or specific target organisms, which can be sensitive in
research. We can dial that up if you want more detail.
