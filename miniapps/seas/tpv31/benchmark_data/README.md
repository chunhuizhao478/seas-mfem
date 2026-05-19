# TPV31 SCEC reference traces

This directory holds the SCEC TPV31 on-fault station traces used as the
reference for `compare_tpv31_traces.py` (plan §R.5 step 1,
`PLAN_phase_R_exact_bimaterial_riemann_rev3.md`).

## Source

Per the plan §R.5 Detailed Requirement 1, the reference traces should
be downloaded from the SCEC Code Verification Working Group submission
area:

  https://strike.scec.org/cvws/cgi-bin/cvws.cgi?problem=31

The TPV31 problem definition (see `tpv31/benchmark_document/TPV31_32_Description_v03.pdf`,
spec p. 12) requires submissions to provide:

  * **30 on-fault stations** (≥ 4 of which are required for the
    automated comparison)
  * Off-fault stations (depth and surface)
  * Rupture-time contour file

The reference set to download into this directory is the **full set**
of on-fault stations from a community-verified submission (e.g., the
SeisSol or DRDG3D entry).  At minimum, the 4 on-fault stations the
plan calls out (`faultst000dp075`, `faultst-080dp075`, etc.) must be
present for `R.5.T-4` to run.

## File layout

Each station produces a SCEC-format `.dat` file with the standard
header:

```
# t  slip-rate-1  slip-rate-2  slip-1  slip-2  traction-1  traction-2
```

The comparison harness (`compare_tpv31_traces.py`) parses the header
to map columns to fields; if the header is absent, it falls back to
the positional schema above.

## Acquisition gate (R.5.T-1)

If acquisition fails (login wall, network unavailable, or the SCEC
server is offline), the plan §R.5 Edge Case explicitly says **STOP
and ask the user for an alternate data source** rather than proceed
without reference data.

The local repository state at the time of Phase R rev-3 implementation
is: **this directory is otherwise empty** (no `.dat` files committed).
Acquisition is a pre-cluster-run step for the user; the script + plan
expect the user to populate this directory before running R.5.T-4.

## Acceptance tolerance

Per plan §R.5 step 5 / R.5.T-4, the pass-band against the SCEC trace
bank is documented in the comparison script's `--tol-rms` /
`--tol-peak` CLI flags.  Tolerances similar to what existing SCEC
dynamic-rupture benchmarks accept (~5–10% peak relative error)
because DRDG3D / SeisSol are different codes from this implementation.
