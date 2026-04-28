# Code Review: TPV205 (SCEC TPV5) — fresh adversarial audit, 2026-04-27 (round 3)

## Review Scope
- Plan: no formal `*_plan_*.md` for TPV205 found in `document/`. Reviewed against
  the implementer's "Implementation Complete" report and against the SCEC TPV5
  spec (`tpv205/benchmark_document/TPV5_forwebsite.pdf`).
- Files reviewed:
  - `config/tpv205_params.hpp`
  - `dynamic/tpv205_friction.hpp`
  - `dynamic/tpv205_setup.hpp`
  - `dynamic/tpv205_substep_iterator.hpp`
  - `dynamic/tpv205_substep_iterator.cpp`
  - `drivers/tpv205_driver.cpp`
  - `Makefile` (TPV205 build rules)
  - `tpv205/mesh/tpv2053d_200m.geo` (mesh refinement)
  - `dynamic/wave_operator.inl` (read-only audit of fault-branch dispatch
    and the R-1600 shared-fault fallback)
- Domain context: `CLAUDE.md` (project root + miniapps/seas) and the
  TPV104 fix history that produced the R-1600/R-1601 fallback.

This is a fresh adversarial pass. Round 1 + 2 findings (R-001..R-015) have
been verified one-by-one against the current source; FIXED items are listed
in the *Round-1/2 verification* section below. **Two NEW critical/quality
findings (R-016, R-017) were discovered on the third pass**, including a
silent-wrong-physics MPI bug that the implementer's smoke test could not
catch.

---

## Findings

### [R-016] CRITICAL `dynamic/wave_operator.inl:~4518` × `drivers/tpv205_driver.cpp:1572` — shared-fault QPs at np > 1 silently run the rate-and-state Brent solver on LSW-repurposed DOFData, producing wrong physics

**Category:** BUG (NEW THIS PASS — class: cross-file MPI / dispatch contract violation)

**Description:**
The driver's banner at line 1572-1574 claims:
```
[tpv205_driver] MPI dispatch (nprocs=N): shared-fault QPs flow through the substep
side-channel via SetSubStepFaultImposedStates.
```
This is **false**. The wave operator's shared-fault path was patched for
TPV104 by R-1600/R-1601 to fall back to the inline solver for shared QPs
because the substep dispatch had a frame-convention divergence at MPI rank
boundaries. From `dynamic/wave_operator.inl` line ~4509-4524:

```cpp
// SHARED FALLBACK: this branch always runs the inline
// EvaluateADER regardless of substep_I_imp_*_flat_.
fault_flux_->EvaluateADER(fdata,
                          I_plus_local, I_minus_local,
                          dt,
                          I_imp_plus, I_imp_minus);
(void)substep_I_imp_plus_flat_;
(void)substep_I_imp_minus_flat_;
(void)substep_n_total_fault_qps_;
```

`fault_flux_->EvaluateADER` calls `Evaluate` → Brent on the **rate-and-
state** friction equation:
```
strength = a · σ_n · asinh(V · exp(ψ/a) / (2·V0))
```

For TPV205 the implementer **repurposed** `data.a = μ_s` and
`data.psi = μ_d`, so Brent is invoked with rate-and-state's `a, ψ` slots
holding LSW's `μ_s, μ_d`. The math then gives nonsensical strengths:

- Rupture-area QP (a=0.677, ψ=0.525, σ_n=120 MPa, V≈0.1 m/s):
  - `C = exp(0.525/0.677)/(2·1e-6) ≈ 1.08·10⁶ 1/m`
  - `f(V) = 0.677·asinh(0.1·1.08·10⁶) ≈ 0.677·12.3 ≈ 8.3`
  - `strength ≈ 120·10⁶·8.3 ≈ 1·10⁹ Pa` (vs LSW: `0.677·120·10⁶ ≈ 8.1·10⁷ Pa`)
  Brent therefore drives V → 0 to balance |τ| ≈ 70-90 MPa, while the
  correct LSW closed form gives V ≈ 0.08 – 1 m/s during rupture.

- Barrier QP (a=10000, ψ=0.525): strength ≈ 10¹³ Pa, so V → 0. (Coincides
  with LSW for barrier QPs — locked either way — but only by accident.)

In production at np > 1 the rupture front therefore **stalls at every MPI
rank boundary** that crosses the rupture area. The implementer's smoke
test did NOT catch this: 4 ADER steps × ~1 ms = 4 ms of physical time, in
which the rupture front has only travelled `cs · 4 ms ≈ 13 m` from the
nucleation centre — far less than any rank-boundary distance on the 200 m
or coarser meshes. At tfinal = 12 s the rupture sweeps the full 30 km × 15 km
area and **must** cross rank boundaries.

The TPV104 fallback was correct because TPV104 *is* rate-and-state; the
fallback just runs the same physics by a different code path. For TPV205
the fallback runs the **wrong physics**.

**Trigger:**
`mpirun -np 2 ./seas_tpv205_driver --tfinal 12.0 ...` — anything with
np > 1 and tfinal long enough for the rupture to reach a rank boundary.

**Actual behavior:**
Shared-fault QPs at rank boundaries report V ≈ 0 from the Brent solve
even when the LSW closed-form gives V ≈ 0.1 – 1 m/s. The rupture front
stalls at the first MPI partition boundary it reaches. The simulation
produces a silently wrong rupture history.

**Expected behavior:**
LSW physics applied identically to interior and shared fault QPs.
Rupture front propagates through MPI rank boundaries unaffected.

**Suggested fix:**
The wave operator's shared-fault fallback (`wave_operator.inl`) is on
the project-wide no-touch list per `CLAUDE.md`. The cleanest immediate
fix is to **abort at parse time when np > 1** so production runs cannot
silently produce wrong rupture physics:

In `drivers/tpv205_driver.cpp`, near the existing nprocs > 1 banner:

```diff
-   if (nprocs > 1 && rank == 0)
-   {
-      std::cout << "[tpv205_driver] MPI dispatch (nprocs=" << nprocs
-                << "): shared-fault QPs flow through the substep "
-                << "side-channel via SetSubStepFaultImposedStates.\n";
-   }
+   if (nprocs > 1)
+   {
+      // R-016: the wave operator's R-1600/R-1601 shared-fault fallback
+      // (wave_operator.inl L~4509-4524) runs `fault_flux_->EvaluateADER`
+      // — Brent on the rate-and-state law — for every shared-fault QP,
+      // ignoring `substep_I_imp_*_flat_`.  TPV205 repurposes `data.a` ←
+      // μ_s and `data.psi` ← μ_d, so Brent on these inputs produces
+      // strengths off by 1-3 orders of magnitude vs the LSW closed
+      // form.  At np > 1 the rupture front stalls at rank boundaries.
+      // Abort here until a TPV205-aware shared-fault dispatch lands in
+      // wave_operator.inl (currently on the no-touch list).  Override
+      // SEAS_FORCE_TPV205_MPI=1 only for diagnostic, single-step runs
+      // where the rupture front is known not to reach a rank boundary.
+      const char *force = std::getenv("SEAS_FORCE_TPV205_MPI");
+      if (!(force && force[0] == '1'))
+      {
+         if (rank == 0)
+         {
+            std::cerr << "[FATAL] TPV205 cannot run at nprocs=" << nprocs
+                      << " > 1: shared-fault QPs would silently use the "
+                      "rate-and-state Brent solver on LSW-repurposed "
+                      "DOFData (wave_operator.inl R-1600 fallback). "
+                      "Use SEAS_FORCE_TPV205_MPI=1 only for single-step "
+                      "diagnostic runs, or run at np=1.  See REVIEW R-016.\n";
+         }
+         MPI_Abort(comm, 16);
+      }
+      else if (rank == 0)
+      {
+         std::cerr << "[WARNING] SEAS_FORCE_TPV205_MPI=1 — running TPV205 "
+                   << "at np=" << nprocs << " is UNSAFE; rupture front "
+                   "will stall at rank boundaries.  Output is for "
+                   "diagnostic use only.\n";
+      }
+   }
```

The mid-term fix is to extend the wave operator's shared-fault dispatch
to consume `substep_I_imp_*_flat_` for non-rate-and-state laws (or to
have a separate `SetSubStepFaultImposedStatesAllQPs` that bypasses the
R-1600 frame-mismatch concern when the iterator's physics doesn't have
the rank-local frame asymmetry that motivated R-1600). That requires
touching `wave_operator.inl` and is out of scope for this fix round.

**Test case:**
```cpp
TEST(TPV205Driver, R016_AbortsAtNprocsGreaterThan1) {
   // Skipped under MPI_COMM_SELF; run with MPI launcher.
   //
   //   mpirun -np 2 ./seas_tpv205_driver --tfinal 12.0 ...
   //   Expected: process exits with error code 16 and stderr includes
   //             "TPV205 cannot run at nprocs=2".
   //
   // Pre-fix:  process completes (silently producing wrong rupture).
   // Post-fix: process aborts.
}

TEST(TPV205Driver, R016_OverrideAllowsNp2WithWarning) {
   //   SEAS_FORCE_TPV205_MPI=1 mpirun -np 2 ./seas_tpv205_driver --tfinal 0.001
   //   Expected: stderr includes "[WARNING] SEAS_FORCE_TPV205_MPI=1".
}
```

---

### [R-017] LOW `dynamic/tpv205_setup.hpp:TPV205StationWriter::WriteStep` — inline μ_eff formula doesn't use the canonical `LSWFrictionCoefficient_TPV205`

**Category:** QUALITY (POSSIBLE — flagged because the duplicate path bypasses
the R-002 strength-barrier short-circuit and would produce wrong outputs if
a future regression lets a barrier QP accumulate slip).

**Description:**
The station-writer's `WriteStep` computes the trace's "μ_eff" column by
inlining the LSW formula at lines 309-314:

```cpp
real_t mu_eff;
if (delta <= 0.0)        { mu_eff = mu_s_qp; }
else if (delta >= d_c_qp){ mu_eff = mu_d_qp; }
else                     { mu_eff = mu_s_qp - (mu_s_qp - mu_d_qp) * (delta / d_c_qp); }
```

This duplicates `LSWFrictionCoefficient_TPV205` (in `dynamic/tpv205_friction.hpp`)
*minus* the R-002 barrier short-circuit:
```cpp
if (mu_s >= 0.5 * TPV205Params::mu_s_barrier) { return mu_s; }
```
The ParaView writer in `drivers/tpv205_driver.cpp:1931` was updated in round 2
to call the canonical helper, but the station writer was overlooked — its
inline copy still has the original barrier-collapse semantics.

In production this is inert because R-003 (the SolveLSW barrier short-
circuit) keeps barrier QP slip ≡ 0, so the inline formula's first branch
returns `mu_s_qp` (correct). The 16 SCEC TPV205 stations in
`kStationsTPV205` are all inside the rupture area (|x| ≤ 12 km, depth ≤ 12 km)
— none are barrier QPs — so the divergence between the two μ_eff
implementations cannot affect any benchmark trace today.

But: any future regression that re-enables slip in a barrier QP, or any
new station added outside the rupture area, would trigger the divergence
silently. The two implementations should be the same one.

**Trigger:**
A barrier QP with δ > 0 (impossible while R-003 is in place) OR a future
station added outside the rupture area.

**Actual behavior:**
For a hypothetical barrier-zone station with δ ≥ d_c, the trace file
column 9 reads μ_d = 0.525 instead of the barrier sentinel μ_s = 10000.

**Expected behavior:**
Trace μ_eff matches the friction physics applied in the iterator
(`LSWFrictionCoefficient_TPV205`).

**Suggested fix:**
Replace the inline computation with a single call to the canonical helper,
matching the round-2 ParaView fix:

```diff
       const DOFData &d = dof_data[idx];

-      // LSW effective friction at current slip magnitude δ = |slip|.
-      const real_t delta = std::sqrt(d.slip1 * d.slip1
-                                     + d.slip2 * d.slip2);
-      const real_t mu_s_qp = d.a;     // μ_s stored in `a`
-      const real_t mu_d_qp = d.psi;   // μ_d stored in `psi`
-      const real_t d_c_qp  = d.Dc;
-      real_t mu_eff;
-      if (delta <= 0.0)        { mu_eff = mu_s_qp; }
-      else if (delta >= d_c_qp){ mu_eff = mu_d_qp; }
-      else                     { mu_eff = mu_s_qp
-                                          - (mu_s_qp - mu_d_qp)
-                                            * (delta / d_c_qp); }
+      // R-017: route through the canonical helper so the strength-
+      // barrier short-circuit (R-002) and any future μ-curve change
+      // applies uniformly to station traces, ParaView, and the
+      // friction solve.
+      const real_t delta = std::sqrt(d.slip1 * d.slip1
+                                     + d.slip2 * d.slip2);
+      const real_t mu_eff = LSWFrictionCoefficient_TPV205(delta,
+                                                           d.a, d.psi, d.Dc);
```

Add `#include "tpv205_friction.hpp"` to `tpv205_setup.hpp` if it is not
already on the include path.

**Test case:**
```cpp
TEST(TPV205Setup, R017_StationWriterMuEffMatchesCanonicalAtBarrier) {
   // Construct a single DOFData with d.a = mu_s_barrier (10000),
   //   d.psi = 0.525, d.Dc = 0.4, d.slip2 = 0.5  (post-d_c).
   // Both LSWFrictionCoefficient_TPV205(0.5, 10000, 0.525, 0.4) and
   // the station-writer's μ_eff column must equal 10000.
   // Pre-fix: station writer outputs 0.525.
}
```

---

## Round-1/2 verification (rerun against current source)

Each round-1/2 finding has been re-checked. Status:

| ID    | Status   | Evidence                                                                                  |
|-------|----------|-------------------------------------------------------------------------------------------|
| R-001 | FIXED    | `Makefile` lines 178-179, 372-383, 1194-1201, 1755-1759 add the TPV205 stanza.            |
| R-002 | FIXED    | `tpv205_friction.hpp:64` short-circuits `mu_s ≥ 0.5·mu_s_barrier`.                        |
| R-003 | FIXED    | `tpv205_friction.hpp:116-128` returns V=0, τ_corr=τ_trial in barrier zone.                |
| R-004 | FIXED    | `tpv205_driver.cpp:160` adds `LSWClosedForm`; banner emits "lsw-closed-form".             |
| R-005 | FIXED    | `tpv205_driver.cpp:1931` calls `LSWFrictionCoefficient_TPV205` for ParaView "state".      |
| R-006 | FIXED    | `tpv205_setup.hpp:397-400` returns empty default surface-station list.                    |
| R-007 | FIXED    | `tpv205_params.hpp:92-97` comment now matches behaviour ("V_ini = 0 is safe").            |
| R-008 | FIXED    | `tpv205_substep_iterator.cpp:94-95` populates `s.Theta`.                                  |
| R-009 | FIXED    | `tpv2053d_200m.geo` adds Surface(300/400) + Field[8..10] + Min field {2,5,6,10}.          |
| R-010 | FIXED    | `tpv205_setup.hpp:381` renames `.z` → `.fault_perp`.                                      |
| R-011 | OPEN     | `Tpv205SubStepIterator::Advance` is still exposed but unused. No regression observed.     |
| R-012 | FIXED    | `tpv205_driver.cpp:486` default mesh path now `tpv205/mesh/tpv2053d_200m.msh`.            |
| R-013 | FIXED    | `tpv205_driver.cpp:1190-1193` overrides `tau2_0 = tau_back` when `disable_nucleation`.    |
| R-014 | FIXED    | `tpv205_driver.cpp:1794` writes `pv_local_x3 = std::abs(fault_coords[i](2))`.             |
| R-015 | FIXED    | `tpv205_driver.cpp:415-420` adds `SubStepFaultImposedGuard` RAII.                         |

**Open from prior rounds:** R-011 (LOW; bit-rot risk on unused `Advance`).

---

## Summary
- Critical issues: **1 NEW** (R-016 silent wrong physics at np > 1)
- Moderate issues: **0 new**
- Low issues: **2** (R-011 prior-round, R-017 NEW μ_eff duplication)
- Plan compliance: **PARTIAL**. 14/15 prior-round findings fixed; the
  silent np > 1 bug (R-016) was not visible to the implementer's 4-step
  smoke test and constitutes a regression in production reliability.
- Verdict: **FAIL — must fix before proceeding to production**.
  R-016 silently produces wrong rupture physics on every np > 1 run that
  reaches the rupture-front rank-boundary crossing time. Either guard
  the driver against np > 1 (the no-touch-respecting fix above) or
  extend the wave operator's shared-fault dispatch to handle LSW (out
  of scope for this round; touches the no-touch list).

## Unreviewed Areas
- The actual numerical effect of the R-016 silent bug at np > 1 was
  not measured against a reference run; the diagnosis is from the
  arithmetic of the rate-and-state vs LSW strength formulas at the
  repurposed parameter values.
- `wave_operator.inl` was read for the R-1600 fallback branch only.
  Other code paths into `fault_flux_->Evaluate*` were not enumerated;
  if there is *any* additional path that calls into the fault flux for
  TPV205 (e.g., a future precomputed-flux wiring), R-016's class of
  bug would re-emerge there too.
- The MPI race in `TPV205SurfaceStationWriter::Open` (multiple ranks
  could open the same file if `mesh.FindPoints` returns hits on more
  than one rank) was not flagged because R-006 makes the default list
  empty, so the race is dormant. If a non-empty list is later added,
  this becomes a real concern.
- `Tpv205SubStepIterator::Advance` (the I/dt_macro variant) is still
  unused; not retested for arithmetic correctness this round.
