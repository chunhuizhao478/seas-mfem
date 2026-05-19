#!/usr/bin/env python3
"""
compare_tpv31_traces.py — A/B comparison harness for the TPV31
through-bi-material verification (plan §R.5 step 5,
PLAN_phase_R_exact_bimaterial_riemann_rev3.md).

Same RMS/peak structure as compare_tpv205_traces.py.  Pass criteria
per plan §R.5 step 5 / R.5.T-4: tolerances similar to what existing
SCEC SEAS / dynamic-rupture benchmarks accept — typically ~5-10%
peak relative error against the SCEC TPV31 trace bank (DRDG3D / SeisSol
reference codes — DIFFERENT solvers, so the tolerance is widely than
the byte-rounding band used for the TPV205 vs gold comparison).

This script is a thin wrapper around the TPV205 comparison utility —
the schema is identical (per-station .dat files; same fault-output
field set).  Documented separately because the wider tolerance band
and acquisition story differ.

Usage:
  python compare_tpv31_traces.py \\
      --new       tpv31/out/results \\
      --reference tpv31/benchmark_data/ \\
      --tol-rms   5e-2 \\
      --tol-peak  1e-1
"""
from __future__ import annotations

import os
import sys

# Re-use the TPV205 implementation; the per-station I/O contract is the
# same (SCEC on-fault station schema).  Path-add the tpv205 scripts dir
# so we can import without packaging.
_HERE = os.path.dirname(os.path.abspath(__file__))
_TPV205_SCRIPTS = os.path.normpath(
    os.path.join(_HERE, os.pardir, os.pardir, "tpv205", "scripts"))
sys.path.insert(0, _TPV205_SCRIPTS)

if __name__ == "__main__":
    # Argument forwarding: defer entirely to compare_tpv205_traces.main.
    # The pass-band defaults from compare_tpv205_traces are 5e-3 / 1e-2
    # (the "gold" comparison).  For TPV31 against the SCEC trace bank
    # the user should pass --tol-rms 5e-2 --tol-peak 1e-1 on the CLI.
    import compare_tpv205_traces
    sys.exit(compare_tpv205_traces.main())
