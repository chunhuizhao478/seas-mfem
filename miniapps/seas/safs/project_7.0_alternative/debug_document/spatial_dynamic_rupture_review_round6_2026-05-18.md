# Code Review (Round 6): Phase 4 driver + Phase H Stage 1 ctor (2026-05-18)

## Review Scope

- **Plan:** `safs/project_7.0_alternative/document/spatial_dynamic_rupture_plan.md` (rev-3)
- **Prior reviews / fix reports consumed (rounds 1-5):** see `debug_document/spatial_dynamic_rupture_review_round{2..5}_2026-05-18.md`
- **Files reviewed (new / modified since round-5 fix):**
  - `miniapps/seas/drivers/spatial_dyn_driver.cpp` (NEW — Phase 4, 1280 LOC)
  - `miniapps/seas/dynamic/wave_operator.hpp` (Phase H.1 Stage 1: new ctor + caches)
  - `miniapps/seas/dynamic/wave_operator.inl` (Phase H.1 ctor body + Phase H.3 heterogeneous CFL)
  - `miniapps/seas/dynamic/wave_operator.cpp` (3-line passthrough)
  - `miniapps/seas/dynamic/fault_face_flux.{hpp,cpp}` (touched — additive only per round-3 R-302)
  - `miniapps/seas/drivers/tpv205_driver.cpp` (modified — diff TBD; not relevant to Phase 4)
  - `miniapps/seas/fault/fault_geometry.hpp` (round-5 state)
  - `miniapps/seas/io/field_coefficient.hpp`, `io/tpv104_checkpoint.hpp` (round-5 state)
  - `miniapps/seas/tests/unit/test_phaseh_wave_operator_constant_parity.cpp` (NEW — 4 tests for the parity gate)
- **Domain context consulted:**
  - `miniapps/seas/CLAUDE.md` (incl. the new Gmsh-v2.2-only parser limitation)
  - `miniapps/seas/dynamic/tpv205_substep_iterator.{hpp,cpp}` (LSW closed-form per-substep solver)
  - `miniapps/seas/dynamic/tpv205_friction.hpp` (LSWFrictionCoefficient_TPV205 vs ForcedRupture)
  - `miniapps/seas/dynamic/spatial_setup.hpp` (round-5 state)
  - `miniapps/seas/spatial/code/spatial_friction.hpp` (LSWFrictionCoefficient_ForcedRupture)

## Findings

---

### [R-601] [CRITICAL] [drivers/spatial_dyn_driver.cpp + dynamic/tpv205_substep_iterator.cpp:100] — Substep iterator uses plain `LSWFrictionCoefficient_TPV205` and ignores `T_forced_rupture`; the SAFS driver's interior-fault dispatch silently runs plain LSW even with `[nucleation]` enabled

**Category:** BUG

**Description:**
The SAFS driver (Phase 4) routes interior-fault dispatch through the substep iterator pattern:

```cpp
// drivers/spatial_dyn_driver.cpp:1214-1216
AdvanceADERWithSubStep_Spatial(wave, substep_iterator, dof_data,
                               fault_coords, Q, dt_step,
                               cfg.numerics.ader_order, t, Q_new);
```

Inside `AdvanceADERWithSubStep_Spatial` (lines 334-428), the driver:
1. Calls `iterator.AdvanceWithSubStepStates(...)` which runs the LSW closed-form solve per substep and accumulates imposed states into `I_imp_*_flat`.
2. Installs `I_imp_*_flat` on the wave operator via `wave.SetSubStepFaultImposedStates(...)`.
3. Calls `wave.AdvanceADER(...)`, which inside `ComputeADERFaceFluxRHS` (wave_operator.inl:3594-3608) gates on `substep_I_imp_*_flat_ != nullptr` and CONSUMES the pre-computed buffer — bypassing the dispatch arm that calls `EvaluateADER_LSW_ForcedRupture`.

But `Tpv205SubStepIterator::StepOneQP_` (tpv205_substep_iterator.cpp:100) only knows about `LSWFrictionCoefficient_TPV205`:

```cpp
const real_t mu_eff = LSWFrictionCoefficient_TPV205(delta,
                                                    d.lsw_mu_s,
                                                    d.lsw_mu_d,
                                                    d.lsw_d_c);
```

The iterator does NOT consume `d.T_forced_rupture` or `d.t0_decay_forced`.  And it discards `t_macro_start` (it's marked `/*t_macro_start*/` at line 142).  So forced-rupture math (`f_2(t)`) NEVER fires in the iterator path.

The shared-fault fallback (wave_operator.inl:4772+) DOES use `EvaluateADER_LSW_ForcedRupture` — but only on shared QPs at MPI partition seams (rare; a few percent of fault DOFs at most).

**Net effect:** The SAFS driver sets `FaultFrictionLaw::LSW_ForcedRupture`, the dispatch arm is wired, the F-9 guard fires correctly, BUT the substep iterator silently ignores the forced-rupture fields and runs plain LSW on every interior-fault QP.  The rupture front near the hypocenter (always on interior faces) never receives the TPV26/27 forced-rupture friction reduction.  The simulation looks like it's running forced rupture (`LSW_ForcedRupture` set, guard happy) but produces PLAIN LSW physics everywhere except partition seams.

This is the worst kind of bug — silent, the dispatch arms LOOK connected, F-9 (round-5) passes, the F-7 unit test (round-3) passes (it tests the FaultFaceFlux method directly, not the iterator), but no rupture nucleates in a real SAFS run.

**Trigger:**
Any SAFS dynamic-rupture run with `[nucleation]` enabled and `--no-sidecar-material` (the only currently-supported production path).  The driver sets `LSW_ForcedRupture`, resolver populates `T_forced_rupture < 1e9` on hypocenter-zone DOFs, iterator's plain LSW ignores them, rupture never nucleates.

**Actual behavior:**
Forced-rupture mechanism wired but dead on the dominant (interior-fault, substep-buffer) path; runs plain LSW; rupture fails to nucleate.

**Expected behavior:**
Per plan §Phase H.6 — `EvaluateADER_LSW_ForcedRupture` (or the equivalent forced-rupture mu_eff computation) MUST fire on every interior fault QP when `LSW_ForcedRupture` is selected.

**Suggested fix:**
Two viable options:

**Option A (preferred — extend the iterator).**  Add forced-rupture awareness to `Tpv205SubStepIterator`.  Pass `T_forced_rupture` / `t0_decay_forced` + `t_macro_start` into `StepOneQP_` and switch `mu_eff` based on a friction-law tag:

```diff
@@ tpv205_substep_iterator.cpp::StepOneQP_ @@
-   const real_t mu_eff = LSWFrictionCoefficient_TPV205(delta,
-                                                       d.lsw_mu_s,
-                                                       d.lsw_mu_d,
-                                                       d.lsw_d_c);
+   real_t mu_eff;
+   if (friction_law_ == FaultFrictionLaw::LSW_ForcedRupture)
+   {
+      mu_eff = mfem::seas::spatial::LSWFrictionCoefficient_ForcedRupture(
+         delta, d.lsw_mu_s, d.lsw_mu_d, d.lsw_d_c,
+         t_sub_start, d.T_forced_rupture, d.t0_decay_forced);
+   }
+   else
+   {
+      mu_eff = LSWFrictionCoefficient_TPV205(delta,
+                                             d.lsw_mu_s, d.lsw_mu_d,
+                                             d.lsw_d_c);
+   }
```

Add a `void SetFrictionLaw(FaultFrictionLaw)` setter on `Tpv205SubStepIterator` and thread the per-substep `t_sub_start` (currently discarded) through.  The driver calls `substep_iterator.SetFrictionLaw(FaultFrictionLaw::LSW_ForcedRupture)` right after `SetSubSteps`.

**Option B (defensive — don't use iterator for forced rupture).**  Refuse to use the substep iterator when `FaultFrictionLaw::LSW_ForcedRupture` is selected; let the inline dispatch arm fire instead (loses per-substep precision but uses the correct physics).  In `AdvanceADERWithSubStep_Spatial`, branch on the friction-law tag: if forced-rupture, skip the iterator + substep-buffer pipeline and call `wave.AdvanceADER(Q, dt_step, ader_order, Q_new)` directly so the dispatch arm fires.

```diff
@@ spatial_dyn_driver.cpp::AdvanceADERWithSubStep_Spatial @@
+  if (wave.GetFaultFrictionLaw() == FaultFrictionLaw::LSW_ForcedRupture)
+  {
+     // Forced rupture requires the inline dispatch path (the substep
+     // iterator doesn't know about f_2(t)).  Call AdvanceADER directly
+     // with no substep buffer so the dispatch arm fires per DOF.
+     wave.AdvanceADER(Q, dt_step, ader_order, Q_new);
+     return;
+  }
   ... existing iterator + buffer pipeline ...
```

Option A is the right physics; Option B is the safety fallback that at least doesn't silently run wrong physics.  Land Option B NOW (one-line guard) and Option A as a follow-up.

**Test case:**
```cpp
// tests/unit/test_spatial_dyn_forced_rupture_dispatch.cpp (NEW)
static void T_SAFS_iterator_path_uses_forced_rupture()
{
   // Build a minimal SAFS-style fixture: 1 fault face with 1 DOF at
   // hypocenter (T_forced_rupture = 0, t0_decay_forced = 0.5).
   // Pre-stress / mu_s chosen so plain LSW is LOCKED but forced
   // rupture at f_2 = 0.5 UNLOCKS (same parameters as F-7).
   //
   // Run one macro step at t = 0.25 through AdvanceADERWithSubStep_Spatial.
   // Assert dof_data[0].slip_rate > 1e-3 (UNLOCKED via f_2).
   //
   // Pre-fix: iterator uses plain LSW, slip_rate stays 0 (LOCKED).
   // Post-fix (Option A): iterator uses forced rupture, slip_rate > 0.
   ...
}
```

---

### [R-602] [MODERATE] [drivers/spatial_dyn_driver.cpp:1178 — paraview_write helper] — ParaView "state" field uses plain `LSWFrictionCoefficient_TPV205`, masking forced rupture in the diagnostic output

**Category:** BUG (diagnostic)

**Description:**
In the per-step ParaView writer at lines 1167-1183, the per-DOF "state" field (the LSW friction coefficient mu_eff) is computed via:

```cpp
pv_local_state(i) = LSWFrictionCoefficient_TPV205(
                       delta_norm,
                       d.lsw_mu_s, d.lsw_mu_d, d.lsw_d_c);
```

For runs with `[nucleation]` enabled, the actual friction coefficient at this DOF is:

```cpp
LSWFrictionCoefficient_ForcedRupture(
   delta_norm, d.lsw_mu_s, d.lsw_mu_d, d.lsw_d_c,
   /*t_now=*/t, d.T_forced_rupture, d.t0_decay_forced)
```

Which can differ substantially (mu_s → mu_d transition during the t0_decay window).  ParaView output silently shows the WRONG mu_eff, masking forced rupture in any visualisation.

This is independent of R-601 — even when R-601 is fixed (iterator uses forced-rupture), this diagnostic line still emits plain LSW.

**Trigger:**
Any SAFS run with `[nucleation]` enabled; open ParaView output and inspect the "state" field at hypocenter DOFs.

**Actual behavior:**
"state" shows mu_eff per plain LSW formula; forced rupture invisible in viz.

**Expected behavior:**
"state" shows mu_eff per the actual friction law in use.

**Suggested fix:**
```diff
@@ drivers/spatial_dyn_driver.cpp:1178 @@
-      pv_local_state(i)             = LSWFrictionCoefficient_TPV205(
-                                         delta_norm,
-                                         d.lsw_mu_s, d.lsw_mu_d,
-                                         d.lsw_d_c);
+      pv_local_state(i) =
+         mfem::seas::spatial::LSWFrictionCoefficient_ForcedRupture(
+            delta_norm,
+            d.lsw_mu_s, d.lsw_mu_d, d.lsw_d_c,
+            time, d.T_forced_rupture, d.t0_decay_forced);
```

(`time` is the lambda capture of the outer `t` loop variable, passed in via the `paraview_write(step_num, time, V_max)` parameter at line 1161.)

**Test case:**
```cpp
// Verify pv_local_state at the hypocenter DOF (T_forced=0, t0=0.5) at
// t = 0.25 equals mu_s + (mu_d - mu_s) * 0.5, NOT mu_s.
```

---

### [R-603] [MODERATE] [drivers/spatial_dyn_driver.cpp:1262 — final checkpoint step number] — Final checkpoint records `step = nsteps` even when the time loop broke early

**Category:** EDGE_CASE

**Description:**
The time loop (lines 1208-1253) iterates `for (int step = step0; step < nsteps; ++step)` with an early-break on `if (dt_step <= 0.0) { break; }` (line 1211) when the clamped step would yield zero or negative dt.  After the loop, the final checkpoint at line 1262 is:

```cpp
WriteTpv104Checkpoint(prefix, t, dt_now, nsteps, Q, dof_data,
                      rank, nprocs, comm, "spatial_dyn");
```

Uses `nsteps` as the step number unconditionally.  If the loop broke early (e.g., `tfinal - t` reached zero with `step < nsteps - 1`), the checkpoint claims the run reached step `nsteps` when it actually reached some earlier step.  A subsequent restart with this checkpoint resumes at `step0 = nsteps`, the loop condition `step0 < nsteps` is false, and the restart immediately exits — losing the user's chance to extend the run.

Also: the in-loop checkpoint at line 1239 writes `step + 1` (correct).  If checkpoint_every divides `nsteps`, the final in-loop checkpoint writes `step + 1 = nsteps`, and the post-loop final checkpoint overwrites with the same `nsteps`.  Harmless duplication when the loop runs to completion; broken when the loop breaks early.

**Trigger:**
Run with a `dt` that doesn't divide `tfinal` evenly so the clamped final step is < `dt`; OR any future change that breaks out of the loop early.

**Actual behavior:**
Final checkpoint records `step = nsteps`; restart skips remaining work.

**Expected behavior:**
Track the actual last completed step number in a counter and record it.

**Suggested fix:**
```diff
@@ drivers/spatial_dyn_driver.cpp:1208-1253 @@
+  int last_step = step0;       // tracks the last completed step
   Vector Q_new(Q.Size());
   real_t V_max_global = 0.0;
   for (int step = step0; step < nsteps; ++step)
   {
      const real_t dt_step = std::min(dt_now, cfg.time.tfinal - t);
      if (dt_step <= 0.0) { break; }
      ...
+     last_step = step + 1;
   }

@@ final checkpoint @@
-      WriteTpv104Checkpoint(prefix, t, dt_now, nsteps, Q, dof_data,
+      WriteTpv104Checkpoint(prefix, t, dt_now, last_step, Q, dof_data,
                             rank, nprocs, comm, "spatial_dyn");
```

**Test case:**
```cpp
// Drive a run with tfinal = 0.25s, dt = 0.1s → last actual step = 3
// (steps 0, 1, 2 reach t = 0.3 > tfinal; on step 2 dt_step = tfinal - t
// = 0.05; on step 3 dt_step = 0, break).
// Run --tfinal=0.25 --dt-initial=0.1.  After completion, inspect
// the final checkpoint: step field must be 3, NOT ceil(0.25/0.1) = 3
// (coincidence — pick a better example like tfinal=0.27, dt=0.1).
```

---

### [R-604] [LOW] [drivers/spatial_dyn_driver.cpp:830-839] — Dead assignment to `ref_normal`; first initialisation `(0,0,1)` is immediately overwritten by `(0,1,0)`

**Category:** QUALITY

**Description:**
Lines 830-839 contain:

```cpp
Vector ref_normal(3);   ref_normal(0) = 0.0;  ref_normal(1) = 0.0;  ref_normal(2) = 1.0;
Vector up_vec(3);       up_vec(0)     = 0.0;  up_vec(1)     = 0.0;  up_vec(2)     = 1.0;
// (8 lines of comment block)
// For a SAFS curvilinear fault the canonical (BP5 / TPV205) choice is
// ref_normal = +y (pointing into the box's east half), up = +z.  ...
ref_normal(0) = 0.0; ref_normal(1) = 1.0; ref_normal(2) = 0.0;
```

The initial `ref_normal = (0, 0, 1)` is dead code — it's overwritten unconditionally to `(0, 1, 0)` 8 lines later.  A reader who scans the file top-down sees `ref_normal = (0, 0, 1)` and assumes that's the convention; only the secondary assignment changes it.  Easy to miss in a refactor.

**Suggested fix:**
```diff
-   Vector ref_normal(3);   ref_normal(0) = 0.0;  ref_normal(1) = 0.0;  ref_normal(2) = 1.0;
-   Vector up_vec(3);       up_vec(0)     = 0.0;  up_vec(1)     = 0.0;  up_vec(2)     = 1.0;
-   // For a SAFS curvilinear fault the canonical (BP5 / TPV205) choice is
-   // ref_normal = +y (pointing into the box's east half), up = +z.  The
-   // FaultBasis::Compute body orients each face's normal against
-   // ref_normal — non-orthogonal ref_normal still resolves a consistent
-   // sign per face.  We use (0,1,0) here; for runs that report sign-
-   // flipped basis vectors at every face the user should adjust via a
-   // TOML knob in a follow-up commit.
-   ref_normal(0) = 0.0; ref_normal(1) = 1.0; ref_normal(2) = 0.0;
+   // For a SAFS curvilinear fault the canonical (BP5 / TPV205) choice
+   // is ref_normal = +y (pointing into the box's east half), up = +z.
+   // The FaultBasis::Compute body orients each face's normal against
+   // ref_normal — non-orthogonal ref_normal still resolves a consistent
+   // sign per face.  Adjustable via a TOML knob in a follow-up commit
+   // if a SAFS run reports sign-flipped basis vectors at every face.
+   Vector ref_normal(3);
+   ref_normal(0) = 0.0;  ref_normal(1) = 1.0;  ref_normal(2) = 0.0;
+   Vector up_vec(3);
+   up_vec(0)     = 0.0;  up_vec(1)     = 0.0;  up_vec(2)     = 1.0;
```

---

### [R-605] [LOW] [dynamic/wave_operator.inl:5249-5294 — heterogeneous CFL bypass] — `ComputeMaxDt` heterogeneous branch is gated on `!per_elem_h_.empty()`; the scalar ctor never touches that vector, so the gate is correct — but in the SAFS driver's scalar path the byte-exact reduction property depends on the scalar ctor's `h_min_` being computed with the SAME per-element h formula, which it isn't audited against R-503 / R-310 round-3 fixes

**Category:** ASSUMPTION / POSSIBLE

**Description:**
The new heterogeneous CFL branch (lines 5264-5292) walks `per_elem_h_` and `per_elem_lmr_` to compute per-element CFL, MPI-reduces.  The scalar-fall-through branch (line 5294) uses `h_min_ / flux_.GetCp()`.

The comment claims (lines 5256-5263):

> For Mode::Constant input ... every element shares the same (lambda, mu, rho), so c_p,e is constant.  The local min over `cfl_factor * cfl * h_e / c_p` is exactly `cfl_factor * cfl * (min_e h_e) / c_p`, and the MPI MIN reduction over these per-rank locals matches the legacy path's reduction on `h_min_` ... Therefore the new path is BYTE-IDENTICAL to ... in Mode::Constant — proven by `test_phaseh_wave_operator_constant_parity` Test 3.

C-3 in `test_phaseh_wave_operator_constant_parity.cpp:285-290` does assert byte-identical `ComputeMaxDt` between scalar and hetero ctors.  But this relies on `h_min_` being computed with the EXACT same per-element h formula as `per_elem_h_` (inscribed-diameter for tets, vol^{1/dim} for non-tets).  The scalar ctor's `h_min_` formula is in an "extreme care" file (per CLAUDE.md) and was NOT touched in this round.  If the scalar ctor's `h_min_` calc subtly differs from `BuildGodunovFluxPool_`'s `per_elem_h_` calc (e.g., different normal direction convention for tet face area, different sign for negative-volume elements), the parity would fail.

C-3 passes today — so the formulas DO agree on the test fixture.  But a tetrahedron with negative-orientation vertices (rare but legal in MFEM) could expose a divergence.  POSSIBLE bug; downgrade to LOW because no concrete failure case is in hand.

**Suggested fix:**
Either (a) leave as-is and rely on C-3 to catch any future drift, or (b) refactor the scalar ctor's h_min_ loop to call into `BuildGodunovFluxPool_`'s per-element formula.  (b) touches the extreme-care file; defer.

---

### [R-606] [LOW] [dynamic/wave_operator.inl:587-616 — delegating ctor placeholder material] — `WaveOperator(MaterialField)` delegates to scalar ctor with `(1, 1, 1)` placeholder for Coefficient mode; MFEM_VERIFY catches it, but a future commit that comments out the verify would silently run with the placeholder material

**Category:** ASSUMPTION

**Description:**
The new `WaveOperator(MaterialField, BoundaryConfig)` ctor delegates to the scalar ctor:

```cpp
: WaveOperator(mesh, order,
               material.mode == MaterialField::Mode::Constant
                  ? material.lambda_const : real_t(1.0),
               material.mode == MaterialField::Mode::Constant
                  ? material.mu_const     : real_t(1.0),
               material.mode == MaterialField::Mode::Constant
                  ? material.rho_const    : real_t(1.0),
               bc)
{
   MFEM_VERIFY(material.mode == MaterialField::Mode::Constant,
               "Phase H Stage 1 gap: ... Coefficient mode requires ...");
   ...
}
```

The placeholder `(1.0, 1.0, 1.0)` is fed into the scalar ctor BEFORE the MFEM_VERIFY in the derived body fires.  The scalar ctor's `flux_` and `h_min_` are now built from `(1, 1, 1)` — physically nonsense but mathematically defined.  If a future commit removes the MFEM_VERIFY (e.g., to "test the hetero path"), Mult() runs on the (1, 1, 1) flux instead of the user's MaterialField — silent wrong physics.

The MFEM_VERIFY makes this safe today.  Defensive concern.

**Suggested fix:**
Construct the placeholder material from the heterogeneous material's element-0 evaluation (so even if Coefficient mode is later enabled and the verify removed, the placeholder is a sane representative value):

```diff
-  : WaveOperator(mesh, order,
-                 material.mode == MaterialField::Mode::Constant
-                    ? material.lambda_const : real_t(1.0),
-                 material.mode == MaterialField::Mode::Constant
-                    ? material.mu_const     : real_t(1.0),
-                 material.mode == MaterialField::Mode::Constant
-                    ? material.rho_const    : real_t(1.0),
-                 bc)
```

This is structurally awkward (the delegating ctor needs a mesh + material to compute the placeholder, before the body runs).  The cleaner long-term fix is to NOT delegate to the scalar ctor at all — refactor the scalar ctor's body into a private helper that both ctors call.  Defer.

---

### [R-607] [LOW] [drivers/spatial_dyn_driver.cpp:816-823 — ProbeNbfPerFace MPI_COMM_NULL guard awkward] — `#ifdef MFEM_USE_MPI` inside a function arg list

**Category:** QUALITY

**Description:**
Lines 818-822:

```cpp
const int nbf_per_face = ProbeNbfPerFace(pmesh, cfg.mesh.order,
                                          fault_int_faces,
                                          fault_shr_faces,
#ifdef MFEM_USE_MPI
                                          comm
#else
                                          MPI_COMM_NULL
#endif
                                         );
```

The `#else MPI_COMM_NULL` branch is unreachable — line 682 has `#error "spatial_dyn_driver requires MFEM_USE_MPI=YES."`.  Dead code that complicates the call site.  Either drop the `#else` branch entirely (the `#error` already guarantees MPI_USE_MPI is defined) or wrap the `ProbeNbfPerFace` declaration with `#ifdef MFEM_USE_MPI` so the function isn't visible in the non-MPI build.

**Suggested fix:**
```diff
-   const int nbf_per_face = ProbeNbfPerFace(pmesh, cfg.mesh.order,
-                                            fault_int_faces,
-                                            fault_shr_faces,
-#ifdef MFEM_USE_MPI
-                                            comm
-#else
-                                            MPI_COMM_NULL
-#endif
-                                           );
+   // MFEM_USE_MPI is unconditional in this driver (see #error at L682).
+   const int nbf_per_face = ProbeNbfPerFace(pmesh, cfg.mesh.order,
+                                            fault_int_faces,
+                                            fault_shr_faces,
+                                            comm);
```

---

### [R-608] [LOW] [drivers/spatial_dyn_driver.cpp:1149-1156 — substep quadrature is uniform-subdivision, not midpoint] — `Tpv205SubStepIterator::SetSubSteps` is called with uniform-subdivision dt and equal-weight quadrature

**Category:** EDGE_CASE / POSSIBLE

**Description:**
Lines 1150-1155:

```cpp
Tpv205SubStepIterator substep_iterator(fault_flux);
{
   const int O = std::max(1, cfg.numerics.ader_order);
   std::vector<real_t> deltaT(O, dt / static_cast<real_t>(O));
   std::vector<real_t> weights(O, 1.0 / static_cast<real_t>(O));
   substep_iterator.SetSubSteps(deltaT, weights);
}
```

Uses **equal subdivision** of the macro step (`deltaT[o] = dt/O`) and **equal weights** (`weights[o] = 1/O`).  The Tpv205 production driver and the plan §Phase 4 commentary suggest **midpoint** or **Gauss–Lobatto** nodes are typical for ADER quadrature.  Equal-weight subdivision is the O=1 limit's continuation but may not be the intended ADER-O production choice.

The driver does adjust per-step: `AdvanceADERWithSubStep_Spatial` re-scales `deltaT` by `dt_step / Σ deltaT_configured` (lines 365-371).  So the initial `deltaT = dt/O` and the per-step rescaling produce `deltaT_scaled = (dt/O) * (dt_step / dt) = dt_step / O` — fine, but still uniform.

POSSIBLE issue.  TPV205 driver may use a different quadrature and document a reason; an audit against `drivers/tpv205_driver.cpp` would resolve this.

**Suggested fix:**
Audit against tpv205_driver's iterator initialisation, then either (a) match the same quadrature for byte-similarity, or (b) document the equal-weight choice explicitly.

---

## Summary

- Critical issues: **1**  (R-601: forced-rupture math dead on interior-fault path)
- Moderate issues: **2**  (R-602 diag mu_eff plain LSW; R-603 final-checkpoint step count)
- Low issues: **5**  (R-604 dead ref_normal assignment; R-605 CFL parity reliance on scalar h_min_; R-606 placeholder material defensive; R-607 dead `#ifdef MFEM_USE_MPI` arg; R-608 uniform-weight substep quadrature audit)
- Plan compliance: **PARTIAL** — Phase 4 driver shipped and runs end-to-end; Phase H Stage 1 (new ctor + caches + heterogeneous CFL) shipped with byte-exact parity test; Phase H Stage 2 (per-element flux dispatch in hot loop) still deferred; Phase 6 sbatches + verify/compare scripts still mostly missing.
- Verdict: **FAIL — must fix R-601 before any SAFS dynamic-rupture run produces correct forced-rupture nucleation.**  The dispatch arm, the guard, the F-7/F-9 unit tests all pass — but the substep iterator that the driver actually uses runs plain LSW and ignores `T_forced_rupture`.  Same silent-failure shape as round-3 R-302 (forced-rupture dispatch arm missing) and round-4 R-401 (GetTime() returning 0); each round of fixes patches one layer and the next layer below reveals the same gap.  R-602/R-603 are independent bugs in the same driver that should land together.

## Unreviewed Areas

- **tpv205_driver.cpp modification** — file was touched but the diff was not audited (out of immediate scope for the SAFS Phase 4 review).
- **fault_face_flux.{hpp,cpp} touch** — round-3 R-302 added `EvaluateADER_LSW_ForcedRupture`; no further changes audited this round.
- **MFEM-v2.2 Gmsh constraint (CLAUDE.md)** — informational; the SAFS mesher emits v2.2 directly.  No code change required.
- **Round 1-5 carry-over findings** — assumed still resolved per the cumulative fix reports; not re-audited.
- **Full TPV byte-exact regression** (`make test-tpv104` / `make test-tpv205` / `make test-bp5-*`) — these gates SHOULD run before Phase H Stage 1 changes merge.  The constant-parity gate `C-1..C-4` is a strong sanity check on the new ctor, but the full per-driver byte-exact suite was not invoked in this round-6 review.
- **MPI np>1 behaviour of the SAFS driver** — the driver was not invoked under `mpirun -np 4`; the per-DOF table walk + MPIContext destructor sequencing + shared-fault dispatch under R-601 would all benefit from a parallel smoke run.
- **The new test `test_phaseh_wave_operator_constant_parity.cpp`** — verified by reading the source; not actually built + run in this review (the test was added by the implementer post-round-5 and should be added to the standard `make test-spatial-phase0-3b` umbrella or a `make test-phaseh` umbrella).
