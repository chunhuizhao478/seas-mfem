#!/bin/bash
# run_bench.sh <out_dir>
# Compile the ADER kernel microbenchmark against the repo's libmfem and run it
# single-core (uncontended) and, under SLURM, full-occupancy (one copy per core,
# memory-contended).  The contended A-vs-B ratio is the kernel plan's Phase-2
# go/no-go (proceed only if contended B beats A by >= 2x — locally, uncontended,
# B was only 1.25x, so the traffic case rests on this run).
#
# Toolchain-portable: the link line is DERIVED from an existing test target via
# `make -n`, so it inherits whatever BLAS/HDF5/MPI the current build uses (conda
# locally, from-source extern on Expanse).  Makefile is NOT modified.
set -o pipefail
OUT="${1:-.}"; mkdir -p "$OUT"
SEAS="$(cd "$(dirname "$0")/../.." && pwd)"          # .../miniapps/seas
BENCH="$SEAS/tests/bench/bench_ader_kernel_variants.cpp"
BIN="$OUT/bench_ader"
[[ -f "$BENCH" ]] || { echo "run_bench: missing $BENCH"; exit 1; }

# Grab the LINK command of a simple test target (it ends in -lmfem ...).  The
# Makefile recipe is backslash-continued, so JOIN continuations first (else the
# -lmfem/-L... tail is dropped and the link fails with undefined mfem symbols).
# awk continuation-join (portable across BSD/macOS + GNU/Linux sed differences).
LINK=$(cd "$SEAS" && make -n seas_test_lts_layout 2>/dev/null \
       | awk '{ if (sub(/\\[ \t]*$/,"")) { printf "%s ", $0 } else { print } }' \
       | grep -E -- '-o +seas_test_lts_layout' | tail -1)
[[ -n "$LINK" ]] || { echo "run_bench: could not derive a link line from 'make -n seas_test_lts_layout'"; exit 1; }

# Compile the bench IN PLACE of the test's own object (extra objects like
# lts_layout.o link harmlessly); force -O3; rename the output.
# Explicit -I<mfem-root> so the bench's #include "linalg/batched/batched.hpp"
# resolves even if a future config.mk strips the -I paths out of MFEM_LINK_FLAGS
# (review P0-4 — currently they carry through, but do not rely on it).
MFEMROOT="$(cd "$SEAS/../.." && pwd)"
CMD=$(printf '%s' "$LINK" \
      | sed -E "s#tests/unit/test_lts_layout\.o#${BENCH}#" \
      | sed -E "s#-o +seas_test_lts_layout#-O3 -I${MFEMROOT} -o ${BIN}#")
echo "run_bench compile:"; echo "  $CMD"
( cd "$SEAS" && eval "$CMD" ) || { echo "run_bench: compile FAILED"; exit 1; }
[[ -x "$BIN" ]] || { echo "run_bench: no binary produced"; exit 1; }

# macOS ONLY: dyld aborts on duplicate LC_RPATH load commands (the conda link
# line repeats -rpath); mirror the repo Makefile's install_name_tool dedup.
# No-op on Linux/Expanse (single rpath), so guarded by uname.
if [[ "$(uname)" == "Darwin" ]] && command -v install_name_tool >/dev/null 2>&1; then
    otool -l "$BIN" 2>/dev/null | awk '/ LC_RPATH$/{getline;getline; print $2}' \
      | sort | uniq -d | while read -r r; do install_name_tool -delete_rpath "$r" "$BIN" 2>/dev/null; done
fi

echo ""; echo "=== bench SINGLE-CORE (uncontended) ==="
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 "$BIN"

if [[ -n "${SLURM_NTASKS}" ]] && command -v srun >/dev/null 2>&1; then
    echo ""; echo "=== bench FULL-OCCUPANCY (${SLURM_NTASKS} independent copies, memory-contended) ==="
    echo "    (each rank runs the whole bench; compare its A/B/C/D us_per_elem to single-core:"
    echo "     the SLOWDOWN is the per-core bandwidth starvation; the A-vs-B ratio is the go/no-go)"
    OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 \
      srun --mpi=pmi2 -n "${SLURM_NTASKS}" "$BIN" 2>&1 \
      | grep -E 'variant|us[_ ]per[_ ]elem|GFLOP|speedup|checksum|Variant|===' | head -60
fi
