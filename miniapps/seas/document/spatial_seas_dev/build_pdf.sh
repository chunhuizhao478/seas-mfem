#!/usr/bin/env bash
# build_pdf.sh — render the spatial_seas (quasi-dynamic) driver plan to a
# house-style PDF.
#
#   ./build_pdf.sh [doc.md] [out.pdf]
#   (defaults: PLAN_spatial_seas_quasidynamic_driver_2026-05-31.md -> same name .pdf)
#
# Uses the shared parts_dev house style (../parts_dev/_assets/pdf-style.tex +
# the same Charter/Avenir Next/Menlo fonts and pandoc flags as _assets/md2pdf.sh)
# PLUS a small doc-local supplement, because THIS doc uses Greek/math glyphs
# (Γ Ω ε η λ ρ ψ ‖ ← ∫ and subscripts ₀ ₁ ₂ and …) that are NOT in the shared
# pdf-style.tex \newunicodechar list, and a few wide C++/call-graph code blocks.
# The shared pdf-style.tex is left untouched.
#
# Requires: pandoc + xelatex (TeX Live) with newunicodechar, fvextra, etoolbox,
#           hyphenat.
set -o pipefail   # NB: NOT set -u — conda activate trips on unbound vars

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ASSETS="$HERE/../parts_dev/_assets"
SRC="${1:-$HERE/PLAN_spatial_seas_quasidynamic_driver_2026-05-31.md}"
OUT="${2:-${SRC%.md}.pdf}"

SUPP="$(mktemp -t spatial_seas_pdf_supp.XXXXXX.tex)"
trap 'rm -f "$SUPP"' EXIT
cat > "$SUPP" <<'TEX'
% Wide plain code blocks (call graphs, C++ signatures) render as the standard
% verbatim env, which ignores \fvset.  Rebind it to fancyvrb's Verbatim so it
% takes a small font + line-break safety (rows stay aligned).
\usepackage{fvextra}
\fvset{breaklines=true,breakanywhere=true,fontsize=\small}
\RecustomVerbatimEnvironment{verbatim}{Verbatim}{%
  breaklines=true,breakanywhere=true,fontsize=\footnotesize}
% Wide pipe-tables: shrink so long \texttt cell content fits the columns.
\usepackage{etoolbox}
\AtBeginEnvironment{longtable}{\footnotesize}
\AtBeginEnvironment{tabular}{\footnotesize}
% Long inline \texttt{...} tokens: allow hyphenation + extra reflow vs margin bleed.
\usepackage[htt]{hyphenat}
\setlength{\emergencystretch}{3.5em}
\AtBeginDocument{\sloppy}
% Glyphs this doc uses that are NOT in the shared pdf-style.tex list.
\newunicodechar{Γ}{\ensuremath{\Gamma}}
\newunicodechar{Ω}{\ensuremath{\Omega}}
\newunicodechar{ε}{\ensuremath{\varepsilon}}
\newunicodechar{η}{\ensuremath{\eta}}
\newunicodechar{λ}{\ensuremath{\lambda}}
\newunicodechar{ρ}{\ensuremath{\rho}}
\newunicodechar{ψ}{\ensuremath{\psi}}
\newunicodechar{‖}{\ensuremath{\|}}
\newunicodechar{←}{\ensuremath{\leftarrow}}
\newunicodechar{∫}{\ensuremath{\int}}
\newunicodechar{₀}{\ensuremath{{}_{0}}}
\newunicodechar{₁}{\ensuremath{{}_{1}}}
\newunicodechar{₂}{\ensuremath{{}_{2}}}
\newunicodechar{…}{\ldots}
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
