# Fix Report: REVIEW_phase3_code_2026-06-03.md — R-301 — 2026-06-04

## Summary
- Findings addressed: **R-301 (MODERATE) FIXED**. The four LOW items (R-302/303/304/305) are tracked below; not requested in this pass.
- Files modified: `drivers/spatial_seas_driver.cpp` (only — the QD driver; not in the [C2] no-touch set).
- Tests added: 0 (see "Test note" — the abort path is not cleanly unit-testable; the fix is verified by the np=4 collective no-regression).
- Build/run: clean; QD driver dry-runs np=1 AND np=4 on the constant + depth-profile configs all exit 0.

## Change Made
**R-301 (MODERATE) — per-DOF `b/f0/V0` uniformity guard was per-rank LOCAL, not MPI-global.**
The §4b `assert_uniform` lambda computed `Min/Max` over the **local** owned DOFs only (no reduction), so a globally non-uniform `b` that happened to be locally uniform on each rank (e.g. a region-partitioned `b` whose region boundary aligns with the mesh partition) bypassed the guard, and Phase 5 would then silently bake the scalar `b_default`.

Fix: reduce the local extrema to **global** via `MPIContext::GlobalMin/GlobalMax` (each an `MPI_Allreduce`) **before** the uniformity comparison. Details:
- Removed the per-rank `if (v.Size() <= 1) return;` early-out — it would let some ranks skip the collective reduction and **deadlock**. `assert_uniform` is called unconditionally on every rank inside the collective `if (cfg.rate_state.has_value())` block, so all ranks reach `GlobalMin/Max` in the same order.
- Empty local view contributes the reduction identity (`numeric_limits<real_t>::max()` for the min, `lowest()` for the max); if no rank has any fault DOF the post-reduce `hi < lo` and the check is skipped (no spurious abort).
- Used the project's `MPIContext` reduction wrappers (the `mpi` instance already constructed at §4b) rather than a raw `MPI_Allreduce`, matching the codebase idiom.

Note: the review observed that for `f0`/`V0` the resolver already rejects per-rule overrides, so only `b` can actually vary — the guard is now correct for all three regardless.

## Verification
- [x] **R-301** — uniformity guard is MPI-global. **np=4 is the decisive test** (it exercises the collective reduction across ranks): constant smoke and depth-profile configs both resolve+apply Phase-3 rate-state and exit 0 at np=4 (no deadlock; uniform `b/f0/V0` passes globally). np=1 also exit 0. Build clean (0 warnings).
- Driver-only change: shared headers untouched → elasticity / faultgeom-parity / bp5-analytic unaffected.

### Test note (why no new automated test)
A direct test of the abort path requires a `[friction.rate_state]` spatial rule that makes `b` vary across ranks but be uniform within each rank — i.e. a partition-dependent region placement (np≥2 with the region landing entirely on one rank). Constructing that deterministically across MPI partitions is impractical, and `MFEM_VERIFY` terminates the process (no `MFEM_USE_EXCEPTIONS`), so it cannot be caught in-process. The fix's correctness is established by construction (global reduce) + the np=4 collective no-regression (the uniform case now goes through the same collective path that would catch a global non-uniformity).

## Unresolved (Phase-3 LOW findings — not requested this pass)
- **R-302** (V_init dip set to 0 vs BP5 `V_zero=1e-20`) — forward note for Phase 5 (confirm the rate-state ODE uses `|V|`); no Phase-3 code change.
- **R-303** (`GetFaultDOFIntegrationPoints` values unverified in the parity test) — test-coverage hardening.
- **R-304** (`SetRateStatePerDOF` second-call σ_n staleness) — single-shot documentation/provenance.
- **R-305** (`SetRateStatePerDOF` unconstrained template parameter) — optional `static_assert` for diagnostics.
These are LOW; recommend addressing R-303 with the Phase-5 fault-operator work that first consumes the IP values.

## Ready for Re-Review: YES
