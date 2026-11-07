#!/usr/bin/env bash
# Build the whole site into docs/ (which GitHub Pages serves).
#
#   ./build.sh            render site + rebuild the shinylive app
#   ./build.sh --site     render the site pages only (fast)
#   ./build.sh --app      rebuild the shinylive app only
#   ./build.sh --data     re-run the R analysis, then rebuild everything
#
# Order matters: quarto renders first, then the app is exported into
# docs/teaching/norovirus so a render can never delete it.
set -euo pipefail
cd "$(dirname "$0")"

QUARTO="${QUARTO:-$HOME/.local/bin/quarto}"
command -v "$QUARTO" >/dev/null 2>&1 || QUARTO=quarto

do_site=1; do_app=1; do_data=0
case "${1:-}" in
  --site) do_app=0 ;;
  --app)  do_site=0 ;;
  --data) do_data=1 ;;
  "")     ;;
  *) echo "unknown option: $1"; exit 1 ;;
esac

if [ "$do_data" = 1 ]; then
  echo "==> re-running analysis (needs full R: ape, treedater, seqinr)"
  Rscript prep/01_prepare_data.R
  Rscript prep/02_prepare_h3n2.R
fi

if [ "$do_site" = 1 ]; then
  echo "==> rendering site"
  "$QUARTO" render
fi

if [ "$do_app" = 1 ]; then
  echo "==> exporting shiny apps to WebAssembly"
  # Both apps export into docs/teaching with a subdir, so they SHARE one copy
  # of the ~80 MB webR runtime in docs/teaching/shinylive instead of two.
  rm -rf docs/teaching/shinylive docs/teaching/norovirus docs/teaching/h3n2
  mkdir -p docs/teaching
  Rscript -e 'shinylive::export("app",      "docs/teaching", subdir = "norovirus")'
  Rscript -e 'shinylive::export("app-h3n2", "docs/teaching", subdir = "h3n2")'
fi

touch docs/.nojekyll
echo "==> done. open docs/index.html, or: $QUARTO preview"
