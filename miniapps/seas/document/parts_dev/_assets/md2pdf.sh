#!/usr/bin/env bash
# md2pdf.sh — render a parts_dev markdown doc to a styled PDF (pandoc + xelatex).
#
#   ./_assets/md2pdf.sh <doc.md> [out.pdf]
#
# Style: serif body (Charter), sans colored headings with a hairline rule
# (Avenir Next), syntax-highlighted code in shaded/bordered boxes (Menlo),
# colored links, numbered sections + TOC.  See _assets/pdf-style.tex.
# Requires: pandoc + xelatex (TeX Live) with tcolorbox, titlesec, newunicodechar.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="${1:?usage: md2pdf.sh <doc.md> [out.pdf]}"
OUT="${2:-${SRC%.md}.pdf}"

pandoc "$SRC" -o "$OUT" \
  --pdf-engine=xelatex \
  --syntax-highlighting=kate \
  --toc --toc-depth=2 --number-sections \
  -V documentclass=article \
  -V papersize=letter \
  -V geometry:margin=1in \
  -V fontsize=11pt \
  -V mainfont="Charter" \
  -V sansfont="Avenir Next" \
  -V monofont="Menlo" \
  -V monofontoptions="Scale=0.84" \
  -V colorlinks=true \
  -H "$HERE/pdf-style.tex"

echo "wrote $OUT"
