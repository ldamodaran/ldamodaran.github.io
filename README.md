# ldamodaran.github.io

Personal site + interactive teaching apps, built with Quarto and
[shinylive](https://posit-dev.github.io/r-shinylive/) and served from
GitHub Pages out of `docs/`.

## Layout

| path | what it is |
|---|---|
| `index.qmd`, `research.qmd`, `teaching.qmd`, `news.qmd`, `cv.qmd` | site pages |
| `posts/` | News entries; copy `_template.qmd` to start one |
| `site/styles.scss` | the matte-green theme — six color variables at the top |
| `app/app.R` | Shiny app 1 — within-host norovirus evolution |
| `app/data/*.rds` | its precomputed inputs (~23 KB) |
| `app-h3n2/app.R` | Shiny app 2 — molecular epidemiology of H3N2 |
| `app-h3n2/data/*.rds` | its precomputed inputs (~56 KB) |
| `prep/01_prepare_data.R` | norovirus phylogenetics → `app/data/` |
| `prep/02_prepare_h3n2.R` | H3N2 phylogenetics → `app-h3n2/data/` |
| `site/data/` | alignments and case data offered as student downloads |
| `docs/` | **build output — served by GitHub Pages. Do not edit by hand.** |

## Adding a News post

```bash
cp posts/_template.qmd posts/my-post-slug.qmd
```

Set `title`, `description` and `date`, delete the `draft: true` line, then
`./build.sh --site` and push. Files starting with `_` are ignored by Quarto,
which is why the template never appears on the site.

Until the first post exists, `quarto render` prints
`WARN: The listing in 'news.qmd' ... doesn't match any files or folders`.
That is expected and harmless; it stops once `posts/` has a published entry.
Remove the "Nothing posted yet" panel from `news.qmd` at that point.

## Build

```bash
./build.sh           # render site + rebuild the app
./build.sh --site    # pages only (fast; use while editing text)
./build.sh --app     # app only
./build.sh --data    # re-run the R analysis first, then rebuild
```

`build.sh` renders Quarto *before* exporting the apps, because they live
inside the Quarto output directory. Both apps are exported into
`docs/teaching/` with `subdir=`, so they **share one copy** of the ~80 MB webR
runtime in `docs/teaching/shinylive/` rather than carrying one each.

## Why the analysis is precomputed

The app runs in the browser via WebAssembly, where only a subset of CRAN is
available. `ape`, `treedater`, `seqinr`, `ggtree` and `treeio` are not, so the scripts in
`prep/` do all of that on a real machine and save plain data frames —
including tree layouts converted to line segments. The apps themselves need
only `shiny`, `bslib`, `ggplot2` and `dplyr`, all of which are available in
WebAssembly (`shiny` and `bslib` ship in the shinylive base image).

Two corrections were made to the original R Markdown tutorials while porting
them, both worth knowing if you reuse the `.Rmd` files:

- treedater stores the strict-clock rate in `$mean.rate`. Both tutorials read
  `$clockRate`, which does not exist — R returns `NULL` silently, so every
  reported rate rendered blank. `prep/` uses `td_rate()`, which checks the
  real field names.
- the norovirus tutorial's closing paragraph describes purifying selection,
  but the data give dN/dS = 1.65 (9 non-synonymous changes against 3
  synonymous). The app describes positive selection and antigenic escape.

Re-run it only when the underlying data or analysis changes:

```bash
./build.sh --data
```

The scripts read their source data from `~/Dropbox/Norovirus_GII.4/` and
`~/Dropbox/scientific_reasources/teaching/BIO463/` (see the paths at the top of
each script). `01` uses `set.seed(463)` so the 350-sequence Sydney subsample is
reproducible.

## Publishing

Served at <https://ldamodaran.github.io/> from the `ldamodaran.github.io`
repo. **The repo must be public** for GitHub Pages to serve it on a free
account; pushing to a private repo works, but the site will not go live.

`docs/` is *not* committed -- it is a build artifact. GitHub Actions
(`.github/workflows/pages.yml`) installs Quarto, R and shinylive, runs
`./build.sh`, checks the output, and publishes it. So a deploy is just:

```bash
git add -A && git commit -m "update site" && git push
```

One-time setup: **Settings -> Pages -> Source: "GitHub Actions"**.

First CI build takes roughly 5-8 minutes, mostly downloading the ~80 MB webR
runtime; after that it is cached and builds run in a couple of minutes.

### Why the build runs in CI

The WebAssembly runtime under `docs/teaching/shinylive/` is ~82 MB. Committing
it meant pushing 68 MB over the wire and keeping it in git history forever,
with another ~40 MB added on every `shinylive` upgrade. Building on GitHub's
runners keeps this repo at a few megabytes.

What CI does *not* do is re-run the phylogenetics -- that needs `ape`,
`treedater`, `seqinr` and Bioconductor. Those results are precomputed into the
small `.rds` files under `app/data` and `app-h3n2/data`, which **are**
committed. Re-run them locally with `./build.sh --data` when the analysis
changes, then commit the updated `.rds` files.

### If CI ever gets in the way

`./build.sh` still works locally and produces exactly what CI produces. To go
back to committing the built site, drop `/docs/` from `.gitignore` and replace
the workflow's build job with a plain `actions/upload-pages-artifact` step --
but expect the large push to come back with it.
