# Code Review: TPV102 Dynamic Rupture Workflow (2026-04-14)

## Review Scope
- Plan: `document/system_dev/dynamic_rupture_plan_v4.md` (Phase 4a: TPV102 setup + driver)
- Files reviewed:
  - `config/tpv102_params.hpp` (144 lines)
  - `dynamic/tpv102_setup.hpp` (455 lines)
  - `drivers/tpv102_driver.cpp` (493 lines)
  - `dynamic/fault_face_flux.hpp` / `fault_face_flux.cpp` (159 lines impl)
  - `dynamic/wave_operator.hpp` / `wave_operator.inl` (840 lines impl)
  - `dynamic/wave_state.hpp` (72 lines)
  - `dynamic/godunov_flux.hpp` / `godunov_flux.cpp` (header + impl)
  - `dynamic/friction_solver.hpp` (96 lines)
  - `friction/state_evolution.hpp` (329 lines)
  - `friction/dieterich_ruina.hpp` (first 100 lines)
  - `tests/unit/test_tpv102_setup.cpp` (190 lines)
  - `tests/verification/test_tpv102_local.cpp` (547 lines)
- Domain context: `CLAUDE.md`, `CODEBASE_GUIDE.md`, SCEC TPV101/102 spec, Tandem reference at `/Users/chunhuizhao/projects/tandem/`

## Findings

### [R-001] [MODERATE] [tpv102_driver.cpp:414-424] — Station output mixes step-averaged and end-of-step quantities

**Category:** BUG

**Description:**
After the RK4 loop, the driver overwrites `V1`, `V2`, and `slip_rate` with RK4-weighted step averages (lines 416-421), but `tau1_corr`, `tau2_corr`, and `sigma_n_corr` retain stage-4 values written by `FaultFaceFlux::Evaluate` during `wave.Mult(Q_tmp, k4)`. The station writer then outputs these mixed quantities together. The friction balance equation `tau = sigma_n * f(V, psi) + eta * V` is **not** satisfied by the output fields at any single time point.

**Trigger:**
Every output step. The inconsistency grows with dt (for coarser meshes / lower CFL).

**Actual behavior:**
`d.V1` = RK4-weighted average over `[t, t+dt]`.  
`d.tau1_corr` = total corrected traction from stage 4 (approximately at time `t+dt`).  
`d.psi` = updated with averaged slip rate.  
These three fields are mutually inconsistent.

**Expected behavior:**
All output fields should correspond to the same time point (end of step). After the RK4 loop, re-evaluate friction at the final state to get consistent `V`, `tau`, `psi` for output. Alternatively, store stage-4 values for V1/V2/slip_rate alongside the stage-4 traction.

**Suggested fix:**
```diff
      // State update using RK4-weighted average slip rate.
      for (int i = 0; i < num_fault_total; i++)
      {
         real_t sr_avg = (sr_k1[i] + 2.0*sr_k2[i] + 2.0*sr_k3[i] + sr_k4[i]) / 6.0;
         dof_data[i].psi = UpdateStateAnalytic(psi_n[i], sr_avg, dof_data[i].Dc,
            dt_step, TPV102Params::f0, TPV102Params::b, TPV102Params::V0);
         dof_data[i].slip_rate = sr_avg;
-        dof_data[i].V1 = (V1_k1[i] + 2*V1_k2[i] + 2*V1_k3[i] + V1_k4[i]) / 6.0;
-        dof_data[i].V2 = (V2_k1[i] + 2*V2_k2[i] + 2*V2_k3[i] + V2_k4[i]) / 6.0;
+        // Use stage-4 values for output consistency with tau1_corr/sigma_n_corr.
+        // The state evolution still uses sr_avg (O(dt^2) coupling, see comment above).
+        dof_data[i].V1 = V1_k4[i];
+        dof_data[i].V2 = V2_k4[i];
         dof_data[i].slip1 += dof_data[i].V1 * dt_step;
         dof_data[i].slip2 += dof_data[i].V2 * dt_step;
      }
```

Note: this changes the slip accumulation to use stage-4 V, which is first-order for slip. The alternative (keep averaged V for slip, use stage-4 V only for output) requires splitting the output V from the integration V:
```diff
+        // For slip accumulation, use RK4-averaged V (correct O(dt^2) integration)
+        real_t V1_avg = (V1_k1[i] + 2*V1_k2[i] + 2*V1_k3[i] + V1_k4[i]) / 6.0;
+        real_t V2_avg = (V2_k1[i] + 2*V2_k2[i] + 2*V2_k3[i] + V2_k4[i]) / 6.0;
+        dof_data[i].slip1 += V1_avg * dt_step;
+        dof_data[i].slip2 += V2_avg * dt_step;
+        // For output, use stage-4 values (consistent with tau1_corr/sigma_n_corr)
+        dof_data[i].V1 = V1_k4[i];
+        dof_data[i].V2 = V2_k4[i];
```

**Test case:**
```cpp
// Verify output consistency: tau = sigma_n * f(V, psi) + eta * V
void test_R001_output_friction_balance()
{
   // After a time step, read station output
   // Check: |tau1_corr| ≈ sigma_n_corr * a * asinh(V1/(2*V0) * exp(psi/a)) + eta_s * V1
   // Tolerance: should be < 1% for dt < 0.01 s
   real_t f_V = d.a * std::asinh(d.slip_rate / (2.0 * V0) * std::exp(d.psi / d.a));
   real_t tau_expected = std::abs(d.sigma_n_corr) * f_V + d.eta_s * d.slip_rate;
   real_t tau_actual = std::sqrt(d.tau1_corr * d.tau1_corr + d.tau2_corr * d.tau2_corr);
   ASSERT(std::abs(tau_actual - tau_expected) / tau_expected < 0.01);
}
```

---

### [R-002] [MODERATE] [tpv102_driver.cpp:329,441] — Station output interval too coarse for SCEC benchmark comparison

**Category:** DEVIATION

**Description:**
Station output is written every `nsteps/100` steps (line 329), producing approximately 100 output points over the 12 s simulation. For a 100 m mesh with order 2, dt ≈ 0.001 s, giving output every ~0.12 s. SCEC reference solutions (PyLith, DR-DG3D) output at 0.007 s and 0.00012 s intervals respectively. The current ~100 points are 60-1000x coarser than reference data, making detailed benchmark comparison (rupture front arrival times, peak slip rate timing) unreliable.

**Trigger:**
Any benchmark comparison against SCEC reference solutions.

**Actual behavior:**
`output_interval = max(1, nsteps / 100)` — hardcoded to ~100 output points.

**Expected behavior:**
Output at a user-configurable time interval (e.g., `--output-dt 0.01`), defaulting to a resolution sufficient for benchmark comparison (at most every 0.01 s).

**Suggested fix:**
```diff
+   real_t output_dt = GetRealArg(argc, argv, "--output-dt", 0.01);
+   int output_interval = std::max(1, static_cast<int>(output_dt / dt));
-   int output_interval = std::max(1, nsteps / 100);
```

**Test case:**
```cpp
void test_R002_output_time_resolution()
{
   // Run 1 s simulation, verify station output has >= 100 points (dt_out <= 0.01 s)
   // Parse station file, check time spacing is uniform and <= 0.01 s
   std::vector<real_t> times = ParseStationTimes("station_flt_0_7.5.dat");
   ASSERT(times.size() >= 100);
   for (size_t i = 1; i < times.size(); i++) {
      ASSERT(times[i] - times[i-1] <= 0.011);  // allow 10% tolerance
   }
}
```

---

### [R-003] [MODERATE] [wave_operator.inl:605] — Ghost element DOF indexing assumes uniform ndof_per_el

**Category:** ASSUMPTION

**Description:**
In `ComputeSharedFaceFluxRHS`, ghost element data is indexed as `nbr_data[c][nbr_idx * ndof_per_el_ + i]` (line 605), where `ndof_per_el_` is computed from element 0. If the mesh contains mixed element types (e.g., hex and tet elements with different polynomial orders or DOF counts), ghost elements with a different type than element 0 would be indexed incorrectly, causing silent data corruption.

The ghost element's actual DOF count `ndof2` is correctly computed (line 559: `fe2->GetDof()`), but it's only used for the shape function evaluation loop bound — not for the data offset.

**Trigger:**
Mixed-element meshes where ghost elements have a different element type than element 0. Does NOT trigger for TPV102 (uniform tet mesh from Gmsh).

**Actual behavior:**
`Q_nbr[c] += shape2(i) * nbr_data[c][nbr_idx * ndof_per_el_ + i]` — uses global `ndof_per_el_` for offset.

**Expected behavior:**
Use the ghost element's actual DOF count or the ParGridFunction's proper DOF indexing.

**Suggested fix:**
```diff
            real_t Q_nbr[NUM_STATE];
            for (int c = 0; c < NUM_STATE; c++)
            {
               Q_nbr[c] = 0.0;
               for (int i = 0; i < ndof2; i++)
               {
-                 Q_nbr[c] += shape2(i) * nbr_data[c][nbr_idx * ndof_per_el_ + i];
+                 // For L2 spaces on uniform meshes, ndof2 == ndof_per_el_.
+                 // Assert this assumption explicitly.
+                 Q_nbr[c] += shape2(i) * nbr_data[c][nbr_idx * ndof_per_el_ + i];
               }
            }
+           MFEM_ASSERT(ndof2 == ndof_per_el_,
+                       "Mixed element types: ghost ndof " << ndof2
+                       << " != ndof_per_el " << ndof_per_el_);
```

**Test case:**
```cpp
void test_R003_ghost_dof_consistency()
{
   // Create a mixed mesh (hex + tet) with shared faces
   // Verify that ghost element DOF counts match ndof_per_el_
   // If they don't, the assertion fires
   // (Currently no mixed-mesh test exists for WaveOperator)
}
```

---

### [R-004] [MODERATE] [tpv102_setup.hpp:60-72] — Pre-stress direction assumes BuildFrame tangent alignment

**Category:** ASSUMPTION

**Description:**
`InitializeFaultDOFs` sets `tau1_0 = tau_ini` (along-strike shear pre-stress) and `tau2_0 = 0` as scalar values, without rotating them to the fault-local coordinate frame. The `FaultFaceFlux::Evaluate` function then adds these directly to the fault-local trial traction: `tau1_total = tau1_0 + tau1_trial`.

This is only correct if `GodunovFlux::BuildFrame` produces tangent vector `t1` aligned with the along-strike direction (x-axis) for the fault face normals. For a vertical planar fault at y=0 with normal = (0, ±1, 0), the standard `BuildFrame` implementation produces `t1 = (±1, 0, 0)`, making `tau1` correspond to the along-strike direction. The shear traction component `sigma_{nt1}` is invariant under simultaneous normal and tangent reversal, so the scalar pre-stress is correct for both face orientations.

However, for non-planar faults, tilted faults, or different BuildFrame implementations, the tangent directions would not align with the intended pre-stress direction, causing the pre-stress to be projected onto the wrong fault-local component.

**Trigger:**
Non-planar fault geometries, non-vertical faults, or a BuildFrame implementation that assigns t1 to a different direction than along-strike.

**Actual behavior:**
Pre-stress is a scalar added directly to the fault-local t1 component without rotation.

**Expected behavior:**
Pre-stress should be specified as a global-frame vector and rotated to the fault-local frame at each face QP, or the current approach should be documented as valid only for planar vertical strike-slip faults.

**Suggested fix (documentation-only for TPV102):**
```diff
 inline void InitializeFaultDOFs(std::vector<DOFData> &dof_data, int ndof,
                                 const std::vector<Vector> &fault_coords)
 {
+   // NOTE: tau1_0 and tau2_0 are specified in the fault-local coordinate frame.
+   // For TPV102's vertical planar fault at Y=0, BuildFrame aligns t1 with the
+   // x-axis (along-strike), so tau1_0 = tau_ini is the along-strike pre-stress.
+   // For non-planar faults, these would need to be rotated from global coordinates
+   // at each QP using the face normal and BuildFrame tangent vectors.
    MFEM_VERIFY(static_cast<int>(fault_coords.size()) >= ndof,
```

**Test case:**
```cpp
void test_R004_prestress_frame_consistency()
{
   // For a vertical fault at Y=0, verify that the fault-local t1 direction
   // from BuildFrame aligns with the along-strike (x) direction.
   real_t nor[3] = {0, 1, 0};  // fault normal
   real_t t1[3], t2[3];
   GodunovFlux::BuildFrame(nor, t1, t2);
   // t1 should be (1, 0, 0) or (-1, 0, 0)
   ASSERT(std::abs(std::abs(t1[0]) - 1.0) < 1e-12);
   ASSERT(std::abs(t1[1]) < 1e-12);
   ASSERT(std::abs(t1[2]) < 1e-12);

   // Also check opposite normal
   real_t nor2[3] = {0, -1, 0};
   GodunovFlux::BuildFrame(nor2, t1, t2);
   ASSERT(std::abs(std::abs(t1[0]) - 1.0) < 1e-12);
}
```

---

### [R-005] [LOW] [tpv102_driver.cpp:416-421] — slip_rate output inconsistent with V1, V2 components

**Category:** BUG

**Description:**
The written `slip_rate` is the average of per-stage magnitudes: `(sr_k1 + 2*sr_k2 + 2*sr_k3 + sr_k4) / 6`, where each `sr_ki = sqrt(V1_ki^2 + V2_ki^2)`. The written `V1` and `V2` are averages of per-stage components: `(V1_k1 + 2*V1_k2 + ...) / 6`.

By Jensen's inequality for the concave function sqrt: `avg(|V_ki|) >= |avg(V_ki)|`, so `slip_rate >= sqrt(V1^2 + V2^2)` in general. The output fields violate the identity `|V| = sqrt(V1^2 + V2^2)`.

For TPV102 (pure mode-II, V2 ≈ 0), the difference is negligible. For mixed-mode problems, the inconsistency would be noticeable.

**Trigger:**
Any simulation with significant V2 component (mixed-mode rupture). Negligible for TPV102.

**Actual behavior:**
`slip_rate != sqrt(V1^2 + V2^2)` in the output.

**Expected behavior:**
Either compute slip_rate from the averaged components (`sqrt(V1_avg^2 + V2_avg^2)`), or use the same source (stage-4 or averaged) for all three fields.

**Suggested fix:**
```diff
-        dof_data[i].slip_rate = sr_avg;
+        // Recompute slip_rate from averaged components for output consistency
+        dof_data[i].slip_rate = std::sqrt(dof_data[i].V1 * dof_data[i].V1
+                                        + dof_data[i].V2 * dof_data[i].V2);
```

Note: the `sr_avg` should still be used for the `UpdateStateAnalytic` call (state evolution uses the magnitude average, which is physically more meaningful for aging law).

---

### [R-006] [LOW] [tpv102_driver.cpp:427-431] — V_max diagnostic uses step-averaged slip rate

**Category:** BUG

**Description:**
After the final state update, `dof_data[i].slip_rate` is set to `sr_avg` (line 419). The V_max tracking (lines 427-431) then uses this averaged value. The instantaneous peak slip rate during the step (max of sr_k1..sr_k4) could be higher than the average. The V_max diagnostic slightly underestimates the true peak.

**Trigger:**
Every step. The underestimate is O(dt) for smoothly varying slip rates, larger during rapid acceleration (nucleation onset).

**Actual behavior:**
`V_max_local = max(dof_data[i].slip_rate)` where slip_rate = step average.

**Expected behavior:**
Track instantaneous V_max across all RK4 stages.

**Suggested fix:**
```diff
      real_t V_max_local = 0.0;
      for (int i = 0; i < num_fault_total; i++)
      {
-        V_max_local = std::max(V_max_local, dof_data[i].slip_rate);
+        // Track peak across all RK4 stages, not just the average
+        V_max_local = std::max(V_max_local,
+           std::max({sr_k1[i], sr_k2[i], sr_k3[i], sr_k4[i]}));
      }
```

---

### [R-007] [LOW] [wave_operator.inl:350-352] — Fault boundary face silently falls back to absorbing BC

**Category:** ASSUMPTION

**Description:**
In `ComputeFaceFluxRHS`, when a face has `bdr_attr == fault_attr` but is a pure boundary face (`e2 < 0`), the `ClassifyBoundaryFace` method returns `FaceBC::Fault`, which is then handled as absorbing BC:
```cpp
case FaceBC::Fault:
   flux_.Absorbing(nor, Q_self, F_h);
   break;
```

This silently absorbs waves at what should be a fault interface. For a correctly split BooleanFragments mesh, fault faces are always interior faces (`e2 >= 0`) and this path is never taken. But if the mesh is misconfigured (e.g., fault surface without a matching element on the other side), waves are silently absorbed instead of producing an error.

**Trigger:**
Misconfigured mesh where a fault boundary element has no corresponding element on the other side of the fault. Does not trigger with correct TPV102 meshes.

**Actual behavior:**
Absorbing BC applied silently at fault boundary face.

**Expected behavior:**
Either: (a) emit a warning the first time this path is hit, or (b) treat as an error.

**Suggested fix:**
```diff
            case FaceBC::Fault:
-              flux_.Absorbing(nor, Q_self, F_h);
+              // Fault boundary face with no element on the other side — mesh error.
+              // Fall back to absorbing BC but warn (first occurrence only).
+              {
+                 static bool warned = false;
+                 if (!warned) {
+                    mfem::out << "WARNING: Fault boundary face " << f
+                              << " has no neighbor element. Using absorbing BC."
+                              << " Check mesh configuration.\n";
+                    warned = true;
+                 }
+              }
+              flux_.Absorbing(nor, Q_self, F_h);
               break;
```

---

## Summary
- Critical issues: 0
- Moderate issues: 4 (R-001, R-002, R-003, R-004)
- Low issues: 3 (R-005, R-006, R-007)
- Plan compliance: **FULL** — all Phase 4a requirements implemented (parameters, setup, driver, RK4 time stepping, station output, nucleation, MPI support)
- Verdict: **PASS WITH FIXES** — the core physics and numerics are correct. The Boxcar function, `ComputeA` spatial distribution, nucleation perturbation, initial equilibrium, trial traction (Eq. 7), friction solver delegation, slip rate decomposition (Eq. 9), corrected traction (Eq. 10), imposed states (Eq. 11-12), and RK4 time stepping all match the plan and SCEC spec. Fix R-001 and R-002 before benchmark comparison; R-003 and R-004 before generalizing to non-TPV102 problems.

## Verified Correct (no findings)

The following aspects were audited and found to be correct:

- **SCEC parameters**: All physical constants (rho, cs, cp, mu, lambda, f0, b, V0, Dc, a_vw, a_vs, sigma_n, tau_ini) match the SCEC TPV101/102 spec.
- **Boxcar function** (`tpv102_params.hpp:73-81`): Matches SCEC Eq. (5) including tanh transition and continuity at boundaries.
- **ComputeA** (`tpv102_params.hpp:86-94`): Along-strike and down-dip Boxcar products correctly cover the full fault extent, with transition outside the fault boundary.
- **NucleationSpatial** and **NucleationTemporal** (`tpv102_params.hpp:98-115`): Match SCEC spec (C-infinity compact support in space, smooth ramp in time).
- **ComputeInitialPsi** (`tpv102_params.hpp:131-138`): Correct inversion of the regularized friction law for equilibrium psi.
- **UpdateStateAnalytic** (`state_evolution.hpp:298-324`): Correct exact solution of aging law in theta-space with numerically stable `expm1` for small V*dt/Dc.
- **Trial traction** (`fault_face_flux.cpp:37-63`): Correct Godunov-state traction from Eq. (7a-c). The +/- convention (Q_plus = self = Elem1, Q_minus = neighbor = Elem2) is internally consistent with the imposed state formulas and the flux accumulation.
- **Slip rate decomposition** (`fault_face_flux.cpp:106-107`): `V1 = V_abs * tau1_total / (strength + eta_s * V_abs)` correctly decomposes slip rate along traction direction, preserving `|V| = V_abs`.
- **Corrected traction** (`fault_face_flux.cpp:110-111`): `|tau_corr_total| = sigma_n * f(V_abs, psi)` — verified algebraically.
- **Imposed states** (`fault_face_flux.cpp:121-144`): Correctly derived from characteristic compatibility relations for both + and - sides.
- **RK4 state management**: psi saved/restored at each sub-stage, nucleation applied at correct stage times (t, t+dt/2, t+dt/2 [implicit from stage 2], t+dt).
- **MPI station ownership** (`tpv102_setup.hpp:210-261`): Global min-distance + rank tiebreaker prevents duplicate station files.
- **Shared fault face detection** (`wave_operator.inl:136-171`): Vertex-key matching correctly identifies shared fault faces without MPI_Allreduce on per-rank-varying arrays.
- **Surface station writer** (`tpv102_setup.hpp:358-450`): FindPoints-based element containment with proper ParMesh ownership.

## Unreviewed Areas
- `friction_solver.cpp` implementation of `SolveBrent` — delegates to the proven `DieterichRuinaFriction::SolveSlipRatePsi` (verified across 62 debug iterations per `CLAUDE.md`).
- `seas_dynamic_operator.hpp` — wrapper for coupled `WaveOperator + FaultFaceFlux`, not directly used by the TPV102 driver (driver does coupling inline).
- Gmsh `.geo` mesh files — **now audited in Round 3 below**.

---

## Second-Round Review: Deferred Performance Findings (2026-04-14)

### Scope

Fresh adversarial analysis of two deferred findings from a prior performance review, plus re-audit of the code paths involved (`ComputeSharedFaceFluxRHS` ghost exchange, fault face flux rotation pipeline, `GodunovFlux::Interior`, `GodunovFlux::BuildFrame/BuildRotation`).

### Deferred Finding: Ghost Exchange Batching (prior R-003)

**Claim:** 36 MPI `ExchangeFaceNbrData()` calls per step (9 components x 4 RK4 stages) add ~18 us on InfiniBand, negligible vs ~1 ms compute. Fix requires vdim=9 FE space with error-prone index remapping.

**Second-round analysis — AGREE with deferral. No correctness issues found.**

Detailed walkthrough of `ComputeSharedFaceFluxRHS` (wave_operator.inl:510-688):

1. **Ghost exchange correctness:** The loop at lines 535-543 copies one scalar component of Q into `ghost_gf_` (a persistent `ParGridFunction`), exchanges, and deep-copies the result into `nbr_data[c]`. The `Vector::operator=` at line 542 performs a deep copy, so subsequent iterations that overwrite `ghost_gf_`'s internal buffer do not corrupt earlier components. Verified: each `nbr_data[c]` owns independent memory. **Correct.**

2. **Early return safety:** Line 525 returns if `n_shared == 0`, skipping the exchange. `ExchangeFaceNbrData()` uses point-to-point MPI (Isend/Irecv to face neighbors), not collectives. A rank with zero shared faces has zero face neighbors, so skipping is safe — no other rank blocks on it. MFEM's ParMesh guarantees bilateral shared faces: if rank A shares with rank B, B also shares with A. **Correct.**

3. **Performance quantification:** 9 exchanges per `Mult()` x 4 stages = 36 exchanges/step. Each exchange involves ~O(ndof_per_el * num_ghost_faces) data per neighbor rank. On InfiniBand with eager protocol (messages < 12 KB), latency is ~0.5 us/call. Total: 18 us/step. For a 100m mesh with dt ~ 0.02 ms and ~50 us compute per `Mult()`, the exchange is ~9% of shared-face compute but <2% of total step compute. **Negligible for verification; marginal for production.**

4. **Batching risk assessment:** A vdim=9 `ParFiniteElementSpace` would change the DOF ordering from component-major (`Q[c * ndof_total + i]`) to MFEM's internal layout (byNODES or byVDIM). This requires changes to `ComputeVolumeRHS`, `ComputeFaceFluxRHS`, `ApplyMassInverse`, all fault flux dispatch code, and the driver's Q initialization. Estimated 15+ files touched with high regression risk. **Risk far exceeds benefit.**

5. **Hidden allocation cost (new observation):** `std::vector<Vector> nbr_data(NUM_STATE)` at line 533 allocates 9 empty Vectors per `Mult()` call. Each `nbr_data[c] = ...` triggers a heap allocation of ~O(num_ghost_elem * ndof_per_el * 8) bytes. With 4 `Mult()` calls/step, this is 36 heap alloc/free cycles per step. For typical ghost counts (~500-2000 elements), each allocation is 40-160 KB. Modern allocators handle this efficiently, and the cost is dominated by the memcpy, not malloc. **Not a bug, but could be pre-allocated as a member for marginal improvement.**

**Verdict: SAFE TO DEFER. No correctness issues. Performance impact <2% of total step.**

---

### Deferred Finding: Redundant Rotation in Fault Face Flux (prior R-005)

**Claim:** 8% per-fault-QP overhead from 2 extra 9x9 mat-vec products. Requires exposing private `ApplySplitFlux()` API.

**Second-round analysis — AGREE with deferral. No correctness issues found.**

Detailed walkthrough of the fault face rotation pipeline (wave_operator.inl:409-449 + godunov_flux.cpp:325-349):

**Current pipeline per fault QP:**

| Step | Operation | Location | FLOPs |
|------|-----------|----------|-------|
| 1 | `BuildFrame(nor, t1, t2)` | wave_operator.inl:411 | ~30 |
| 2 | `BuildRotation + BuildRotationInverse` | wave_operator.inl:412-414 | ~2×200 |
| 3 | `Q_plus_local = Tinv * Q_self` (+ Q_minus) | wave_operator.inl:417-426 | 2×81 = 162 |
| 4 | `FaultFaceFlux::Evaluate(...)` | wave_operator.inl:430 | ~500 (Brent) |
| 5 | `Q_imp_plus_g = T * Q_imp_plus` (+ minus) | wave_operator.inl:435-444 | 2×81 = 162 |
| 6a | `BuildFrame(nor, t1, t2)` **(redundant)** | godunov_flux.cpp:330 | ~30 |
| 6b | `BuildRotation + BuildRotationInverse` **(redundant)** | godunov_flux.cpp:335-336 | ~2×200 |
| 6c | `Tinv * Q_imp_plus_g` **(redundant: Tinv·T = I)** | godunov_flux.cpp:340-341 | 2×81 = 162 |
| 6d | `ApplySplitFlux(...)` | godunov_flux.cpp:345 | 2×81 = 162 |
| 6e | `T * F_rot` | godunov_flux.cpp:348 | 81 |

Steps 5 + 6a-6c rotate global → local → global → local, a round-trip that equals identity. The optimized path would skip steps 5, 6a, 6b, 6c entirely and call `ApplySplitFlux` directly on the fault-local imposed states from step 4, then rotate the flux back with one `T.Mult`:

**Optimized pipeline:** Steps 1-4 (unchanged) → `ApplySplitFlux(Q_imp_plus, Q_imp_minus, F_rot)` → `T.Mult(F_rot, F_h)`. Saves ~790 FLOPs/QP (steps 5 + 6a-6c) out of ~1990 total = **~40% of per-QP compute**, not 8%.

However, fault QP compute is a fraction of total step compute. For a 100m TPV102 mesh with ~15000 fault QPs and ~300000 volume QPs, the fault QP compute is ~5-10% of total. So the overall savings would be ~2-4% of step time.

**Correctness of the double rotation:** `BuildFrame` is deterministic for a given normal (lines 279-300 of godunov_flux.cpp: uses fixed `up = (0,0,1)` with fallback to `(1,0,0)`). Both calls in the pipeline receive the same `nor`, so they produce identical `t1, t2`. The round-trip `Tinv * T` introduces floating-point error of O(9*eps) ≈ 2e-15 per component. For stress values ~O(10^8 Pa), this gives ~O(10^-7) error — 8 orders of magnitude below the Brent solver tolerance. **Numerically benign.**

**Heap allocation overhead (new observation):** Each fault QP allocates 2 `DenseMatrix(9)` objects in wave_operator.inl (lines 412-413) and 2 more inside `Interior()` (godunov_flux.cpp:333-334). Each `DenseMatrix(9)` calls `SetSize(9,9)` which heap-allocates 648 bytes. That's 4 heap allocs per fault QP × 15000 QPs × 4 stages = 240,000 alloc/free cycles per step. This is a larger overhead than the FLOP cost of the redundant rotation. Pre-allocating the rotation matrices at the face level (outside the QP loop) would help both the redundant-rotation and heap-allocation issues. **Not a correctness bug, but the dominant overhead source in this path.**

**Verdict: SAFE TO DEFER for verification. The fix is simple (make `ApplySplitFlux` public, skip `Interior()` for fault faces) and low-risk, but the performance gain (~2-4% of step time) is not critical for TPV102 verification. Consider implementing before production Frontera runs, bundled with rotation matrix pre-allocation.**

---

### New Finding from Second-Round Audit

No new correctness bugs found in the ghost exchange or rotation paths. The following was verified during this round:

- **Shared fault face flux sign consistency:** Both ranks compute `F_h` with their own outward normal. Since `A_{-n}^+ = -A_n^-`, the flux from rank B's perspective is `F_h_B = -F_h_A`. Both ranks subtract: `rhs -= F_h`. Rank A gets `rhs_A -= F_h_A`, rank B gets `rhs_B -= (-F_h_A) = rhs_B += F_h_A`. This gives the correct DG accumulation (Elem1 -= F, Elem2 += F). **Correct.**
- **R-001 fix verification:** The driver now saves per-stage corrected tractions (lines 340-346) and computes RK4-weighted averages for output (lines 442-444). All output fields (V1, V2, slip_rate, tau1_corr, tau2_corr, sigma_n_corr) are now consistently step-averaged. The friction balance `|tau| = sigma_n * f(V, psi) + eta * V` still does not hold exactly for the averaged quantities (due to nonlinearity of f), but the error is O(dt^2), consistent with the acknowledged coupling accuracy. **Fix correctly applied.**
- **`BuildFrame` determinism:** Verified that the same normal always produces the same tangent frame (godunov_flux.cpp:279-300 uses a fixed `up` vector with a dot-product threshold). This guarantees the double rotation is numerically a round-trip (Tinv·T ≈ I to machine precision). **Correct.**
- **`BuildRotation`/`BuildRotationInverse` initialization:** Both methods call `SetSize(9,9)` then `= 0.0` (lines 199-200, 241-242), zeroing all 81 entries before filling. The Voigt-index stress rotation block and the velocity rotation block are filled completely. Off-diagonal blocks remain zero. **No uninitialized memory.**

---

## Round 3: Mesh and Boundary Attribute Audit (2026-04-14)

### Scope

Adversarial comparison of:
- Current meshes: `tpv102/mesh/tpv102_200m.geo`, `tpv102/mesh/tpv102_100m.geo`
- Reference mesh: `/Users/chunhuizhao/projects/farms_benchmark/meshgenerator/tpv1013d/tpv1013d_100m.geo` (SeisSol/Thomas Ulrich — same geometry for TPV101 and TPV102)
- Tandem reference: `/Users/chunhuizhao/projects/tandem/examples/tandem/3d/tpv102.geo`
- Driver boundary config: `drivers/tpv102_driver.cpp:147-150`
- Parameter file: `config/tpv102_params.hpp:48-49`

### Cross-Reference Table

| Property | Current TPV102 mesh | Tandem reference | farms_benchmark (SeisSol) |
|----------|-------------------|------------------|--------------------------|
| **Fault along-strike** | 30 km | **36 km** | **36 km (= 36e3 m)** |
| **Fault down-dip** | 15 km | **18 km** | **18 km (= 18e3 m)** |
| Domain size | 60×60×30 km | 100×100×50 km | 120×120×60 km |
| Coordinate units | km | km | **meters** |
| Free-surface tag | 1 | 1 | **101** |
| Fault tag | 3 | 3 | **103** |
| Absorbing tag | 5 | 5 | **105** |
| Bottom BC | absorbing (tag 5) | natural (tag 1) | absorbing (tag 105) |
| Mesh method | OpenCASCADE BooleanFragments | OpenCASCADE BooleanFragments | Extrude + Surface In Volume |

---

### [R-008] [CRITICAL] [tpv102_200m.geo, tpv102_100m.geo] — Fault surface only covers VW zone (30×15 km), missing 3 km transition zone

**Category:** BUG

**Description:**
Both current TPV102 meshes define the fault rectangle as 30 km along-strike × 15 km down-dip (lines 28-29: `l_f = 30; W_f = 15`). This exactly matches the velocity-weakening (VW) zone defined in the SCEC spec.

However, the fault surface must extend **beyond** the VW zone to include the VW→VS transition region. The SCEC spec defines a 3 km transition width (`w = 3 km`) via the Boxcar function. Both the Tandem reference mesh (`l_f = 36; W_f = 18`) and the farms_benchmark/SeisSol mesh (`Fault_length = 36e3; Fault_width = 18e3`) use a 36 km × 18 km fault — extending 3 km past the VW zone in each direction.

With the current 30×15 km fault:
- All fault DOFs have `a = a_vw = 0.008` (inside VW zone, Boxcar = 1)
- The fault terminates abruptly at the VW boundary with no VS arrest zone
- Rupture hits a free edge (welded → fault transition) instead of a physically correct VS barrier
- This creates artificial stress concentrations at fault tips and prevents proper rupture arrest

**Trigger:**
Every simulation. The rupture front reaches the fault edge at ~2-4 s and interacts with the artificial free edge for the remainder of the 12 s simulation.

**Actual behavior:**
Fault surface is 30 km × 15 km. Beyond the fault edge, faces are welded (regular interior Godunov flux). Rupture arrests by hitting a rigid wall, not by friction-controlled deceleration.

**Expected behavior:**
Fault surface is 36 km × 18 km. DOFs in the 3 km transition zone have `a` values between `a_vw` (0.008) and `a_vs` (0.016). Rupture decelerates gradually as it enters the VS zone, matching the SCEC benchmark physics.

**Suggested fix:**
```diff
 // Fault dimensions (km)
-l_f = 30;  // along-strike extent
-W_f = 15;  // depth extent
+l_f = 36;  // along-strike extent (VW 30 km + 3 km transition on each side)
+W_f = 18;  // depth extent (VW 15 km + 3 km transition beyond fault bottom)
```

Apply identically to both `tpv102_200m.geo` and `tpv102_100m.geo`.

Also update `tpv102_params.hpp` to reflect the actual fault surface size:
```diff
    // Fault geometry
-   static constexpr real_t fault_length = 30e3; ///< Along-strike [m]
-   static constexpr real_t fault_depth = 15e3;  ///< Down-dip [m]
+   static constexpr real_t fault_length = 36e3; ///< Along-strike [m] (VW 30 + 3 km transition each side)
+   static constexpr real_t fault_depth = 18e3;  ///< Down-dip [m] (VW 15 + 3 km transition at bottom)
```

**Test case:**
```cpp
void test_R008_fault_extends_into_transition_zone()
{
   // Verify that fault DOFs exist in the transition zone
   // At along-strike = 16 km (inside transition, outside VW zone):
   real_t a_16km = ComputeA(16e3, 7.5e3);
   // Should be between a_vw and a_vs (in transition)
   TEST_ASSERT(a_16km > TPV102Params::a_vw);
   TEST_ASSERT(a_16km < TPV102Params::a_vs);

   // At along-strike = 14 km (inside VW zone):
   real_t a_14km = ComputeA(14e3, 7.5e3);
   TEST_NEAR(a_14km, TPV102Params::a_vw, 1e-12);

   // Fault must have DOFs at 16 km for the transition to matter
   // With 30 km fault, max along-strike = 15 km → no DOFs at 16 km
   // With 36 km fault, max along-strike = 18 km → DOFs exist at 16 km
}
```

---

### [R-009] [CRITICAL] [tpv102_driver.cpp:147-150] — Hardcoded boundary attributes {1, 3, 5} incompatible with farms_benchmark mesh {101, 103, 105}

**Category:** BUG

**Description:**
The driver hardcodes boundary attributes:
```cpp
BoundaryConfig bc;
bc.natural_attrs = {1};     // free surface
bc.fault_attr = 3;          // fault
bc.absorbing_attrs = {5};   // absorbing sides + bottom
```

The farms_benchmark reference mesh (`tpv1013d_100m.geo`) uses:
- Physical Surface **101** = free surface (top)
- Physical Surface **103** = fault
- Physical Surface **105** = absorbing (sides + bottom)

When loading the farms_benchmark mesh with the current driver, **no faces match any boundary attribute**. The WaveOperator's `ClassifyBoundaryFace()` returns `Absorbing` for all unrecognized attributes, meaning ALL boundary faces (including the free surface and fault) are treated as absorbing. The fault coupling never activates. The free surface reflects no waves. The simulation produces completely wrong results with no error or warning.

**Trigger:**
Running the driver with the farms_benchmark mesh file without changing the boundary config.

**Actual behavior:**
All boundary faces treated as absorbing BC. Zero fault coupling. Free surface acts as absorbing boundary. No error message.

**Expected behavior:**
Driver should accept boundary attributes via CLI flags, or auto-detect common conventions (1/3/5 vs 101/103/105).

**Suggested fix:**
```diff
    real_t cfl_factor = GetRealArg(argc, argv, "--cfl", 0.5);
+   int bc_free = GetIntArg(argc, argv, "--bc-free", 1);
+   int bc_fault = GetIntArg(argc, argv, "--bc-fault", 3);
+   int bc_absorb = GetIntArg(argc, argv, "--bc-absorb", 5);

    // ...

    BoundaryConfig bc;
-   bc.natural_attrs = {1};     // free surface
-   bc.fault_attr = 3;          // fault
-   bc.absorbing_attrs = {5};   // absorbing sides + bottom
+   bc.natural_attrs = {bc_free};
+   bc.fault_attr = bc_fault;
+   bc.absorbing_attrs = {bc_absorb};
```

Usage with farms_benchmark mesh:
```
./seas_tpv102_driver --mesh tpv1013d_100m.msh --mesh-scale 1.0 \
    --bc-free 101 --bc-fault 103 --bc-absorb 105
```

**Test case:**
```cpp
void test_R009_boundary_attr_mismatch_detected()
{
   // Load a mesh with attr 101 for free surface
   // Set bc.natural_attrs = {1}
   // Verify that zero boundary faces are classified as FreeSurface
   int free_surface_count = 0;
   for (int b = 0; b < mesh.GetNBE(); b++) {
      if (bc.natural_attrs.count(mesh.GetBdrAttribute(b))) {
         free_surface_count++;
      }
   }
   // With mismatched attrs, this should be zero — detect and warn
   if (free_surface_count == 0) {
      WARN("No free surface faces found. Check --bc-free matches mesh tags.");
   }
}
```

---

### [R-010] [CRITICAL] [tpv102_driver.cpp:80] — Default mesh-scale 1000 wrong for meter-unit meshes

**Category:** BUG

**Description:**
The driver defaults to `--mesh-scale 1000.0` (line 80), designed for the current km-unit meshes. The farms_benchmark mesh uses meters as coordinates. Running with the default scale multiplies all coordinates by 1000, producing a domain of 120,000 km × 120,000 km × 60,000 km. The CFL time step scales with `h_min / cp`, producing `dt ≈ 100 / 6000 * cfl ≈ 0.0017 s` for the unscaled mesh but `dt ≈ 0.0017 * 1000 = 1.7 s` for the wrongly-scaled mesh. The simulation runs a few steps with a massive time step and produces garbage.

**Trigger:**
Running the driver with a meter-unit mesh file without passing `--mesh-scale 1.0`.

**Actual behavior:**
Mesh coordinates multiplied by 1000. Domain and time step are 1000× too large.

**Expected behavior:**
Either (a) default to `--mesh-scale 1.0` (safest — no silent scaling), or (b) add a validation check that the domain size matches `TPV102Params::domain_half` after scaling.

**Suggested fix:**
```diff
-   real_t mesh_scale = GetRealArg(argc, argv, "--mesh-scale", 1000.0);
+   real_t mesh_scale = GetRealArg(argc, argv, "--mesh-scale", 1.0);
```

Add a post-load validation:
```diff
+   // Validate domain size after scaling
+   Vector bbox_min, bbox_max;
+   serial_mesh.GetBoundingBox(bbox_min, bbox_max);
+   if (mesh_scale != 1.0) {
+      bbox_min *= mesh_scale;
+      bbox_max *= mesh_scale;
+   }
+   real_t domain_x = bbox_max(0) - bbox_min(0);
+   if (rank == 0 && (domain_x < 1e3 || domain_x > 1e6)) {
+      std::cerr << "WARNING: Domain X-extent = " << domain_x
+                << " m (after scale=" << mesh_scale
+                << "). Expected ~60000-120000 m for TPV102. "
+                << "Check --mesh-scale.\n";
+   }
```

**Test case:**
```cpp
void test_R010_mesh_scale_validation()
{
   // Load farms_benchmark mesh with default scale=1.0
   // Verify domain extent is ~120000 m
   Mesh m("tpv1013d_100m.msh");
   Vector bmin, bmax;
   m.GetBoundingBox(bmin, bmax);
   real_t Lx = bmax(0) - bmin(0);
   TEST_ASSERT(Lx > 50e3 && Lx < 200e3);  // reasonable range in meters
}
```

---

### [R-011] [MODERATE] [tpv102_200m.geo, tpv102_100m.geo] — Domain half-size 30 km too small for 12 s simulation

**Category:** DEVIATION

**Description:**
The current mesh domain is 60 km × 60 km × 30 km (half-size = 30 km). The P-wave speed is 6 km/s. P-wave travel from the fault to the nearest absorbing boundary is 30 km / 6 km/s = 5 s. Reflected energy returns to the fault at ~10 s, contaminating the last 2 s of the 12 s simulation.

Reference meshes use larger domains:
- Tandem: 100 km × 100 km × 50 km → reflections at ~16 s (safe for 12 s)
- farms_benchmark: 120 km × 120 km × 60 km → reflections at ~20 s (safe)

With first-order absorbing BC (which reflects ~30% of energy at oblique incidence), even a 50 km domain has detectable contamination. The current 30 km domain produces significant reflected-wave artifacts in the last 2 s.

**Trigger:**
Any simulation past ~10 s. Affects rupture arrest dynamics and post-rupture slip rate evolution.

**Actual behavior:**
Domain half-size = 30 km. Reflected P-waves reach the fault at ~10 s; reflected S-waves at ~17 s.

**Expected behavior:**
Domain half-size >= 50 km (matching Tandem). Using the farms_benchmark mesh (60 km half-size) solves this.

**Suggested fix:**
If keeping the current OpenCASCADE mesh format:
```diff
 // Domain half-sizes (km)
-Lx = 30;   // along-strike
-Ly = 30;   // fault-normal
-Lz = 30;   // depth
+Lx = 60;   // along-strike (matches farms_benchmark)
+Ly = 60;   // fault-normal
+Lz = 60;   // depth
```

Or simply use the farms_benchmark mesh with corrected attributes.

Also update `tpv102_params.hpp`:
```diff
    // Domain
-   static constexpr real_t domain_half = 30e3;  ///< Domain half-size [m] (each direction)
+   static constexpr real_t domain_half = 60e3;  ///< Domain half-size [m] (each direction)
```

**Test case:**
```cpp
void test_R011_domain_free_of_reflections()
{
   // P-wave travel time from fault to boundary:
   real_t t_reflect = TPV102Params::domain_half / TPV102Params::cp;
   // Round-trip (to boundary and back):
   real_t t_roundtrip = 2.0 * t_reflect;
   // Must exceed simulation time for clean results:
   TEST_ASSERT(t_roundtrip > TPV102Params::t_final);
   // With domain_half = 60 km: t_roundtrip = 120/6 = 20 s > 12 s ✓
   // With domain_half = 30 km: t_roundtrip = 60/6 = 10 s < 12 s ✗
}
```

---

### [R-012] [MODERATE] [tpv102_driver.cpp:147-150] — No validation that boundary attributes exist in the loaded mesh

**Category:** EDGE_CASE

**Description:**
After loading the mesh and setting `bc.natural_attrs = {1}`, `bc.fault_attr = 3`, `bc.absorbing_attrs = {5}`, the driver does not verify that these attributes actually exist in the mesh. If the mesh uses different tags (e.g., 101/103/105), the simulation proceeds with zero fault faces and zero free-surface faces, producing silently wrong results.

A diagnostic check after mesh loading would catch this immediately.

**Trigger:**
Any mesh with non-standard boundary attribute numbering.

**Actual behavior:**
Silent incorrect simulation.

**Expected behavior:**
Emit an error or warning if zero boundary faces match the configured attributes.

**Suggested fix:**
```diff
+   // Validate boundary attributes exist in the mesh
+   {
+      int n_free = 0, n_fault = 0, n_absorb = 0;
+      for (int b = 0; b < pmesh.GetNBE(); b++)
+      {
+         int attr = pmesh.GetBdrAttribute(b);
+         if (bc.natural_attrs.count(attr)) { n_free++; }
+         if (attr == bc.fault_attr) { n_fault++; }
+         if (bc.absorbing_attrs.count(attr)) { n_absorb++; }
+      }
+      int n_free_g = n_free, n_fault_g = n_fault, n_absorb_g = n_absorb;
+#ifdef MFEM_USE_MPI
+      MPI_Allreduce(&n_free, &n_free_g, 1, MPI_INT, MPI_SUM, comm);
+      MPI_Allreduce(&n_fault, &n_fault_g, 1, MPI_INT, MPI_SUM, comm);
+      MPI_Allreduce(&n_absorb, &n_absorb_g, 1, MPI_INT, MPI_SUM, comm);
+#endif
+      if (rank == 0)
+      {
+         std::cout << "BC faces — free: " << n_free_g
+                   << ", fault: " << n_fault_g
+                   << ", absorb: " << n_absorb_g << "\n";
+         if (n_fault_g == 0) {
+            std::cerr << "ERROR: No fault faces found (attr="
+                      << bc.fault_attr << "). Check mesh tags.\n";
+         }
+         if (n_free_g == 0) {
+            std::cerr << "WARNING: No free-surface faces found (attr="
+                      << *bc.natural_attrs.begin()
+                      << "). Check mesh tags.\n";
+         }
+      }
+      MFEM_VERIFY(n_fault_g > 0,
+                  "No fault faces with attr=" << bc.fault_attr
+                  << " found in mesh. Available attrs: check gmsh Physical Surface tags.");
+   }
```

**Test case:**
```cpp
void test_R012_missing_fault_attr_detected()
{
   // Create mesh with attr 103 for fault, set bc.fault_attr = 3
   // Verify the validation fires MFEM_VERIFY
   BoundaryConfig bc;
   bc.fault_attr = 3;
   int n_fault = 0;
   for (int b = 0; b < mesh.GetNBE(); b++) {
      if (mesh.GetBdrAttribute(b) == bc.fault_attr) n_fault++;
   }
   // Should be zero for mismatched mesh → error
   ASSERT(n_fault == 0); // demonstrates the bug
}
```

---

### Round 3 Summary

- Critical issues: 3 (R-008, R-009, R-010)
- Moderate issues: 2 (R-011, R-012)
- Mesh compliance: **FAIL** — the current meshes use wrong fault dimensions (30×15 instead of 36×18 km), and the driver cannot load the correct farms_benchmark mesh without code changes.

**Recommended action sequence:**
1. Fix R-009 first (add `--bc-free`, `--bc-fault`, `--bc-absorb` CLI flags) — enables using the farms_benchmark mesh immediately.
2. Fix R-010 (change default `--mesh-scale` to 1.0) — prevents silent wrong scaling.
3. Fix R-012 (add BC attribute validation) — catches future mismatches.
4. Fix R-008 (update fault dimensions in .geo files) — if generating new meshes rather than using the farms_benchmark mesh.
5. Fix R-011 (increase domain size) — if generating new meshes.

Alternatively: use the farms_benchmark mesh directly (`tpv1013d_100m.msh`) with `--mesh-scale 1.0 --bc-free 101 --bc-fault 103 --bc-absorb 105` after fixing R-009.
