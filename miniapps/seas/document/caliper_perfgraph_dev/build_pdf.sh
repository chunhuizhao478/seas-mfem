#!/usr/bin/env bash
# build_pdf.sh — render the ADER hot-path optimization plan to a house-style PDF.
#
#   ./build_pdf.sh [doc.md] [out.pdf]
#   (defaults: ader_hotpath_optimization_plan_2026-05-31.md -> same name .pdf)
#
# Uses the shared parts_dev house style (../parts_dev/_assets/pdf-style.tex +
# the same Charter/Avenir Next/Menlo fonts and pandoc flags as _assets/md2pdf.sh)
# PLUS a small doc-local supplement, because THIS doc has two things md2pdf.sh's
# defaults don't handle:
#   (1) the §1 perfgraph code block is 123 columns wide — at 1in margins it must
#       use a small code font (and fvextra break-safety) to keep rows aligned;
#   (2) a few math glyphs (partial, Sigma, varphi, 2/3, 1/2, leftarrow) are not
#       in the shared pdf-style.tex mapping list.
# The shared pdf-style.tex is left untouched.
#
# Requires: pandoc + xelatex (TeX Live) with tcolorbox, titlesec, newunicodechar,
#           fvextra, etoolbox, hyphenat.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ASSETS="$HERE/../parts_dev/_assets"
SRC="${1:-$HERE/ader_hotpath_optimization_plan_2026-05-31.md}"
OUT="${2:-${SRC%.md}.pdf}"

SUPP="$(mktemp -t ader_pdf_supp.XXXXXX.tex)"
trap 'rm -f "$SUPP"' EXIT
cat > "$SUPP" <<'TEX'
% Wide plain (no-language) code blocks — incl. the 123-col perfgraph — render as
% the standard `verbatim` env, which ignores \fvset.  Rebind it to fancyvrb's
% Verbatim so it takes a small font + line-break safety (rows stay aligned).
% Language-tagged blocks use pandoc's Highlighting (fancyvrb) -> handled by \fvset.
\usepackage{fvextra}
\fvset{breaklines=true,breakanywhere=true,fontsize=\fontsize{6.4}{7.7}\selectfont}
\RecustomVerbatimEnvironment{verbatim}{Verbatim}{%
  breaklines=true,breakanywhere=true,fontsize=\fontsize{6.4}{7.7}\selectfont}
% Wide pipe-tables: shrink so long \texttt cell content fits the columns.
\usepackage{etoolbox}
\AtBeginEnvironment{longtable}{\footnotesize}
\AtBeginEnvironment{tabular}{\footnotesize}
% Long inline \texttt{...} tokens: allow hyphenation + extra reflow vs. margin bleed.
\usepackage[htt]{hyphenat}
\setlength{\emergencystretch}{3.5em}
\AtBeginDocument{\sloppy}
% Glyphs this doc uses that are NOT in the shared pdf-style.tex list; map inside
% \AtBeginDocument so they win after the template's unicode-math claims them.
\AtBeginDocument{%
  \newunicodechar{∂}{\ensuremath{\partial}}%
  \newunicodechar{Σ}{\ensuremath{\Sigma}}%
  \newunicodechar{φ}{\ensuremath{\varphi}}%
  \newunicodechar{⅔}{\ensuremath{\tfrac{2}{3}}}%
  \newunicodechar{½}{\ensuremath{\tfrac{1}{2}}}%
  \newunicodechar{←}{\ensuremath{\leftarrow}}%
}
TEX

pandoc "$SRC" -o "$OUT" \
  --pdf-engine=xelatex \
  --syntax-highlighting=kate \
  --toc --toc-depth=2 --number-sections \
  -V documentclass=article -V papersize=letter \
  -V geometry:margin=1in -V fontsize=11pt \
  -V mainfont="Charter" -V sansfont="Avenir Next" \
  -V monofont="Menlo" -V monofontoptions="Scale=0.84" \
  -V colorlinks=true \
  -H "$ASSETS/pdf-style.tex" \
  -H "$SUPP"

echo "wrote $OUT"
