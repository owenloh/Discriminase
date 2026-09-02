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

The static app has no backend, so analytics is a separate Railway service:
self-hosted **Umami** (`ghcr.io/umami-software/umami:postgresql-latest`) backed by the
project's **Postgres**. `web/analytics.mjs` loads its `script.js` and sends coarse
`umami.track(...)` events; `web/config.js` holds the website id and script URL, and the
whole thing is a no-op if those are left empty.

Log **coarse, non-identifying** events only (nuclease used, panel size, whether the
target came from a name/accession/upload, how many guides were found) — never uploaded
FASTA contents or long sequences. `analytics.mjs` enforces this: uploads log
filename + length, and pasted DNA is only sent verbatim at or under 120 bases.

### Cost: Umami is billed for idle RAM, not traffic

Railway bills resident memory around the clock (~$10/GB-month) and CPU separately
(~$20/vCPU-month). At this site's traffic CPU is effectively zero, so **Umami's whole
bill is its idle footprint** — a Next.js server that sat at ~344 MB flat while serving
roughly 600 requests a month. It is not a leak: the high-water mark and the floor were
within 55 MB of each other over a week. It is just Node's baseline, billed 24/7.

These service variables cap that baseline (set on the **Umami** service in Railway):

| Variable | Value | Why |
| --- | --- | --- |
| `NODE_OPTIONS` | `--max-old-space-size=160 --max-semi-space-size=2` | V8 otherwise sizes its heap against the host's RAM and never feels pressure to collect. The second flag drops the scavenger's semi-spaces from the 16 MB default. |
| `UV_THREADPOOL_SIZE` | `2` | Four libuv worker threads is oversubscribed for this load. |
| `DISABLE_TELEMETRY` | `1` | Stops Umami's periodic phone-home to umami.is. |
| `NEXT_TELEMETRY_DISABLED` | `1` | Same for Next.js. |

Keep `--max-old-space-size` at 160 or above: the container runs `prisma migrate deploy`
before the server starts, and that CLI needs real headroom. If Umami ever crashloops on
boot after a version bump, raise it before assuming the image is broken.

**Serverless (app sleeping) does not work here.** Railway decides a service is idle from
its *outbound* packets, and Umami holds a Prisma connection pool open against Postgres,
so it would rarely sleep. Worse, Railway documents that the first request to a slept
service can return 502 — and for analytics that first request *is* the pageview you
wanted to record.

If the idle cost still isn't worth it, Umami Cloud's free tier (100k events/month, 3
websites, 6 months retention) speaks the same `umami.track` API: point `src` and
`websiteId` in `web/config.js` at it and delete both the Umami and Postgres services.
