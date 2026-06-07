# Code Review: Part A implementation (central-flux contrast guard)

## Review Scope
- Plan: `document/bimaterial_dev/PLAN_unified_bimaterial_volume_and_fault_2026-06-06.md` (Part A)
- Files reviewed (the Part-A diff):
  `dynamic/godunov_flux_bimaterial.{hpp,cpp}`, `dynamic/wave_operator.hpp`,
  `dynamic/bimaterial_wave_operator.{hpp,inl}`, `spatial/code/spatial_friction.{hpp,cpp}`,
  `drivers/spatial_dyn_driver.cpp`, `tests/unit/test_bimaterial_contrast_guard.cpp`,
  `tests/unit/test_tpv_config_parse.cpp`, `Makefile`, the two TPV31 contrast-guard sbatch.
- Domain context: `CLAUDE.md`, `miniapps/seas/CLAUDE.md`, the implementation completion report.
- Adversarial pass over the IMPLEMENTED code (not the plan).  Verified correct (not findings):
  the `SetMixedFluxMode` override rebuilds the FULL set (base) then filters (mine), keeping
  the IMPL-8 size invariant and being idempotent on re-invocation; no iterator invalidation
  (erase is deferred to after both loops); the R-002 fallback assert; the histogram MPI_Reduce
  is collective-safe (gated on a rank-uniform `guard_on`, called by all ranks via the
  collective `SetMixedFluxMode`); `tol<0` adds no statements/collectives (byte-exact); the
  scalar path never calls the bimaterial build, so scalar + `tol>=0` is inert.

## Findings

### [R-001] MODERATE — test_bimaterial_contrast_guard: no DIRECT test that the guarded operator dispatches the reclassified face as upwind

**Category:** EDGE_CASE (test completeness for a critical-path behavior)

**Description:**
Phase 2 proves the reclassified face leaves `op_guard`'s `central_flux_face_set_`, and Phase 3
proves (on `op_none` vs `op_disabled`) that upwind != central with positive dissipation.  But
nothing asserts the actual RUNTIME consequence: that `op_guard` (the operator WITH the guard
enabled) dispatches the reclassified face through the bi-material UPWIND, i.e. that its
`InteriorFaceFlux_(removed, ...)` equals the upwind result and is per-side (two-valued), NOT
the single-valued central.  This is the whole point of the fix; it is currently only inferred.

**Trigger:** a regression that reclassifies the face in the SET but fails to change the
DISPATCH (e.g. a future change to the `mf_on_ && count()>0` branch) would pass all current
tests.

**Actual behavior:** the runtime dispatch effect is asserted only indirectly.

**Expected behavior:** assert `op_guard.InteriorFaceFlux_(removed)` == `op_none.InteriorFaceFlux_(removed)`
(both per-side upwind) and that it is two-valued (F_e1 != F_e2 under the contrast).

**Suggested fix (append to Phase 3 of `test_bimaterial_contrast_guard.cpp`):**
```cpp
   // Direct runtime check: the GUARDED operator dispatches the reclassified face as
   // the (per-side, two-valued) bi-material upwind — identical to the None operator.
   real_t Fg1[NUM_STATE], Fg2[NUM_STATE];
   op_guard.InteriorFaceFlux_(removed, Q_self, Q_nbr, nor, Fg1, Fg2);
   real_t worst_vs_none = 0.0, worst_two_valued = 0.0;
   for (int c = 0; c < NUM_STATE; ++c)
   {
      worst_vs_none    = std::max(worst_vs_none,    std::abs(Fg1[c] - F_up1[c]));
      worst_vs_none    = std::max(worst_vs_none,    std::abs(Fg2[c] - F_up2[c]));
      worst_two_valued = std::max(worst_two_valued, std::abs(Fg1[c] - Fg2[c]));
   }
   TEST_ASSERT(worst_vs_none <= 1.0e-9 * std::max(1.0, std::abs(F_up1[VX])),
               "guarded operator dispatches the reclassified face as the bi-material "
               "upwind (== None operator), not central");
   TEST_ASSERT(worst_two_valued > 1.0e-3,
               "the reclassified face's upwind dispatch is per-side (F_e1 != F_e2 "
               "under the contrast), i.e. the central single-valued path was abandoned");
```

**Test case:** the snippet above IS the test.

---

### [R-002] LOW — bimaterial_wave_operator.inl: `n_reclass_g` double-counts shared reclassified faces under MPI_Reduce

**Category:** QUALITY (diagnostic accuracy; currently inert)

**Description:**
`n_reclass = reclassified.size()` includes BOTH local and shared reclassified faces, and is
summed across ranks via `MPI_Reduce`.  A shared face appears in `central_flux_face_set_` on
both owning ranks, so a shared reclassification would be counted twice in the printed
`reclassified central -> upwind` total — the same double-count the histogram deliberately
avoids by excluding shared faces (R-005).  This is currently INERT (the local-side stub never
reclassifies a shared face), but the count is wrong the moment a real cross-rank exchange lands.

**Trigger:** np>1 with a real (non-stub) cross-rank material exchange that reclassifies a shared face.

**Actual behavior:** shared reclassified faces would be double-counted in `n_reclass_g`.

**Expected behavior:** the printed global count should reflect each face once.

**Suggested fix:** count only LOCAL reclassified faces for the global sum (disjoint across
ranks), mirroring the histogram's local-only convention.
```diff
+      // Count LOCAL reclassified faces only for the global sum (disjoint across
+      // ranks => no double-count; shared reclassification is inert under the stub).
-      long long n_reclass  = static_cast<long long>(reclassified.size());
+      long long n_reclass  = static_cast<long long>(n_reclass_local);
```
with `n_reclass_local` incremented next to each `reclassified.push_back(mesh_face)` in the
LOCAL loop only (not the shared loop).

**Test case:** not unit-testable without an np>1 real-exchange fixture (which does not exist);
LOW + inert, so a test is not required (rule 4) — the fix is a correctness-when-activated cleanup.

---

### [R-003] LOW — bimaterial_wave_operator.inl: matrix + mixed_flux=none + tol>=0 runs an extra collective + prints a zero histogram

**Category:** QUALITY

**Description:**
With `mixed_flux=none` (empty central set) on the matrix path, `BuildPerFaceCentralFluxMatrices_`
is still called; if `tol>=0` the `guard_on` block runs an `MPI_Reduce` and prints a
`reclassified 0/0` histogram.  Harmless (collective-safe, rank-uniform) and `tol>=0` is opt-in,
but it is noise for a no-central run.

**Trigger:** `--mixed-flux none --mixed-flux-contrast-tol 0.05` on the matrix path.

**Suggested fix:** skip the histogram block when the central set was empty AND nothing was
reclassified (no diagnostic value), keeping it collective-uniform:
```diff
-   if (guard_on)
+   if (guard_on && (n_local_hist > 0 || !reclassified.empty()))
```
(NOTE: this must remain rank-uniform.  `n_local_hist`/`reclassified` are per-rank, so a rank
with 0 corridor faces would skip while another prints — that is fine because only rank 0 prints
AND the MPI_Reduce would then be skipped on some ranks => DEADLOCK.  Therefore do NOT gate the
MPI_Reduce on per-rank emptiness; instead keep the MPI_Reduce unconditional under `guard_on`
and gate only the rank-0 PRINT on the GLOBAL total being > 0.)  Reworked safe fix:
```diff
-      if (hrank == 0)
+      if (hrank == 0 && (n_local_g > 0 || n_reclass_g > 0))
```

---

## Summary
- Critical issues: 0
- Moderate issues: 1 (R-001 — missing direct dispatch-as-upwind test)
- Low issues: 2 (R-002 shared reclass double-count [inert]; R-003 none+tol noise/print gating)
- Plan compliance: FULL (A1-A4 implemented; deviations documented in the completion report
  — `ContrastValue` helper and the L2 `0.5*cp` mechanism constant — are justified and correct)
- Verdict: PASS WITH FIXES — no critical bugs found.  The implementation is correct on the
  verified points above; apply R-001 (valuable: directly proves the runtime dispatch change),
  R-002 and R-003 (cleanups).

## Unreviewed Areas
- Parallel (np>1) execution of the guard + histogram MPI_Reduce — not runnable locally; reasoned
  to be collective-safe (rank-uniform `guard_on`, all ranks call via collective SetMixedFluxMode).
- The TOML-parse path under a real `make` (toml11 absent in worktree); verified via the main-repo
  toml11 `-I` override (98/98).
- TPV31 production effect of the guard — Frontera, user-submitted (no local full-mesh).
