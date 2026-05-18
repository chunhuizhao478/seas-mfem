# Code Review (Round 4): Post-fix audit of round-3 changes (2026-05-18)

## Review Scope

- **Plan:** `safs/project_7.0_alternative/document/spatial_dynamic_rupture_plan.md` (rev-3)
- **Prior reviews / fix reports consumed:**
  - Round-1: `debug_document/spatial_dynamic_rupture_review_2026-05-18.md`
  - Round-1 fix: `debug_document/spatial_dynamic_rupture_fix_2026-05-18.md`
  - Round-2: `debug_document/spatial_dynamic_rupture_review_round2_2026-05-18.md`
  - Round-2 fix: `debug_document/spatial_dynamic_rupture_fix_round2_2026-05-18.md`
  - Round-3: `debug_document/spatial_dynamic_rupture_review_round3_2026-05-18.md`
  - Round-3 fix: `debug_document/spatial_dynamic_rupture_fix_round3_2026-05-18.md`
- **Files reviewed (post-round-3-fix state):**
  - `miniapps/seas/dynamic/wave_operator.hpp`
  - `miniapps/seas/dynamic/wave_operator.inl` (lines 3617 and 4539 — new dispatch arms)
  - `miniapps/seas/dynamic/spatial_setup.hpp` (centroid fix + IP-aware overload)
  - `miniapps/seas/dynamic/godunov_flux_pool.cpp` (R-311 + R-312)
  - `miniapps/seas/fault/fault_geometry.hpp` (new BP5 ctor — R-304 + R-310 fixes)
  - `miniapps/seas/tests/unit/test_spatial_setup.cpp` (+S-7, +S-8)
  - `miniapps/seas/tests/unit/test_phaseh_godunov_flux_pool.cpp` (+P-5, +P-6)
  - `miniapps/seas/tests/unit/test_phaseh_lsw_forced_rupture.cpp` (+F-7)
- **Domain context consulted:**
  - `mfem/linalg/operator.hpp` (TimeDependentOperator::GetTime / SetTime contract)
  - `miniapps/seas/dynamic/seas_dynamic_operator.hpp` (SetTime fan-out)
  - `miniapps/seas/dynamic/tpv104_substep_iterator.hpp` (Advance(..., t_macro_start, ...))
  - `miniapps/seas/drivers/tpv104_driver.cpp`, `tpv205_driver.cpp` (no SetTime call sites)
  - `miniapps/seas/fault/fault_geometry.hpp:Print()`, `IsVelocityWeakening()`, `ComputeBP5Params()`
  - `miniapps/seas/CLAUDE.md`
- **Test sweep (round-4 baseline):**
  - All round-3 tests pass (515/515 from the round-3 fix report).

## Findings

---

### [R-401] [CRITICAL] [wave_operator.inl:3636 / 4549 — GetTime() source] — `GetTime()` returns the wave operator's stored `t`, which NO SEAS caller ever updates; forced rupture silently never fires

**Category:** BUG (introduced by R-302 fix)

**Description:**
The round-3 R-302 fix added two dispatch arms that route `FaultFrictionLaw::LSW_ForcedRupture` through `EvaluateADER_LSW_ForcedRupture(..., GetTime(), ...)`.  `GetTime()` is the MFEM `TimeDependentOperator` accessor that returns the operator's stored `t` member (initialised to 0 at ctor, only updated by `SetTime(t)`).

`grep "\.SetTime\b" miniapps/seas/{dynamic,drivers,solver}/` returns **zero call sites** inside the SEAS production code base (test_wave_operator.cpp's `wave.SetTime(t)` is the only caller).  Specifically:
- `Tpv104SubStepIterator::Advance(..., real_t t_macro_start, ...)` takes the absolute time as a parameter but never calls `wave.SetTime(t_macro_start)` on the wave operator — the time lives in `t_macro_start` and never reaches the WaveOperator's member.
- `drivers/tpv104_driver.cpp`, `tpv205_driver.cpp`, `tpv102_driver.cpp` never call `wave.SetTime(t)`.
- `seas_dynamic_operator.hpp:49` overrides `SetTime(t)` to forward to `wave_->SetTime(t_)`, but no caller invokes the overridden setter either (`seas_op->SetTime` is also not called anywhere in seas/).

So in the dispatch arms added by R-302, `GetTime()` returns 0.0 for every Mult call.  The forced-rupture helper then evaluates:

```cpp
if (t_now < T_forced)                  f2 = 0.0;        // t_now == 0
else if (t_now < T_forced + t0_decay)  f2 = ...;
else                                    f2 = 1.0;
```

For a DOF inside the nucleation zone with `T_forced(r) > 0` (every DOF except the exact hypocenter), `0 < T_forced` always — so `f2 = 0` forever, regardless of how much physical time has elapsed.  Forced rupture **never fires**.

For the hypocenter DOF (`T_forced = 0`), `0 < 0` is false and `0 < 0 + t0_decay` is true → `f2 = (0 - 0)/t0_decay = 0`.  Even at the hypocenter, `f_2` is identically 0 at every Mult call.

The dispatch arm is wired but the time it passes is stale.  R-302 LOOKS fixed (the round-3 F-7 unit test passes because the test calls `EvaluateADER_LSW_ForcedRupture` directly with an explicit `t_now=0.25`, bypassing `GetTime()`), but in the actual wave-operator dispatch the forced-rupture mechanism is a no-op.

This is the round-3 reviewer's exact warning text in the R-302 suggested fix:

> Note `GetTime()` (or whatever the existing ADER closure uses for time, e.g. `t_now`) needs to be the macro-step time the iterator passes in; verify against `Tpv104SubStepIterator::Step` to find the canonical name.

The round-3 fix did not verify, and the assumption that "GetTime() will be set by the driver" is false — there is no production caller that does.

**Trigger:**
SAFS driver (once Phase 4 lands) calls `wave.SetFaultFrictionLaw(LSW_ForcedRupture)`; runs `wave.Mult(Q, dQdt)` inside a time loop.  Every dispatch arm receives `GetTime() == 0` and `f_2 == 0` always.

**Actual behavior:**
Forced-rupture friction reduction never engages; the simulation runs as plain LSW with no nucleation, which produces no rupture for sub-critical pre-stress configurations.

**Expected behavior:**
The time passed to `EvaluateADER_LSW_ForcedRupture` MUST be the current simulation time, propagated from the driver's macro-step loop.

**Suggested fix:**
Two correct options:

**Option A (preferred, narrow):** plumb the macro-step time through the ADER call chain.  `Tpv104SubStepIterator::Advance` already receives `t_macro_start`; pass it to `wave.SetTime(t_macro_start)` at the top of `Advance` so the dispatch reads a current time.  Mirror the pattern in `Tpv205SubStepIterator::Step` (or whatever advance entry the LSW path uses).

```diff
@@ tpv104_substep_iterator.cpp::Advance @@
 void Tpv104SubStepIterator::Advance(...,
     real_t dt_macro,
     real_t t_macro_start,
     ...)
 {
+   // R-401: propagate the macro-step time to the wave operator so
+   // EvaluateADER_LSW_ForcedRupture (Phase H) reads a current time
+   // from WaveOperator::GetTime().  The forced-rupture path is a
+   // no-op without this — TimeDependentOperator::t defaults to 0
+   // and nothing else updates it.
+   wave_.SetTime(t_macro_start);
    ...
 }
```

(The iterator currently doesn't hold a `WaveOperator&`; if so, add a setter or plumb the operator into the iterator.)  Repeat the same change in `Tpv205SubStepIterator::Step` (or equivalent LSW path).

**Option B (defensive, broader-but-safe):** add an `MFEM_ASSERT` inside the new dispatch arms that GetTime() > 0 OR the friction-law is set up such that f_2 is intended to never fire (T_forced == 1e9 everywhere).  This at least crashes loudly instead of silently producing wrong physics:

```diff
@@ wave_operator.inl:3636 (and 4549) @@
                     else if (fault_friction_law_ ==
                              FaultFrictionLaw::LSW_ForcedRupture)
                     {
+                       MFEM_ASSERT(GetTime() > 0.0 ||
+                                   fdata.T_forced_rupture >= 1.0e8,
+                                   "WaveOperator: LSW_ForcedRupture "
+                                   "dispatch at GetTime()=0 with "
+                                   "T_forced_rupture < 1e8 — driver "
+                                   "forgot to call wave.SetTime(t).");
                        fault_flux_->EvaluateADER_LSW_ForcedRupture(
                           fdata,
                           I_plus_local, I_minus_local,
                           dt,
                           GetTime(),
                           I_imp_plus, I_imp_minus);
                     }
```

(Option A is the real fix; Option B is the diagnostic guard.)

**Test case:**
```cpp
// tests/unit/test_phaseh_dispatch_time.cpp (NEW)
static void T_PHASEH_dispatch_uses_current_sim_time()
{
   // Build minimal mesh + WaveOperator + FaultFaceFlux + 1 fault DOF
   // with T_forced=0, t0_decay=0.5.  Set SetFaultFrictionLaw(LSW_FR).
   // Call wave.SetTime(0.25) then wave.Mult(Q, dQdt).
   // Read dof_data[0].slip_rate (or some downstream signature).
   //
   // Then RESET dof_data, do NOT call SetTime, repeat.
   //
   // The two outputs MUST differ — the dispatch must consume the
   // current time, not the stored 0.  If they match, the dispatch is
   // still reading stale t.
   ...
}
```

---

### [R-402] [CRITICAL] [fault_geometry.hpp new BP5 ctor:241-270 + Print/IsVelocityWeakening readers] — Skipping `ComputeBP5Params()` leaves `a_values_/dc_values_/eta_values_/V_init_vec_` size 0; `Print()` and `IsVelocityWeakening()` crash or silently return wrong values

**Category:** BUG (introduced by R-310 fix)

**Description:**
The round-3 R-310 fix removed `ComputeBP5Params()` from the new BP5 ctor's body to avoid silent garbage from feeding world (x, z) into the BP5 (x2, x3) analytic functions.  The fix relied on the comment claim:

> Skip ComputeBP5Params() in the new ctor: the SAFS driver overwrites per-DOF analytic values via InitializeFaultDOFs_Spatial (Phase 5c) before any consumption.  ComputeSAFSParams (templated overload) will lazy-initialise the BP5 arrays if any caller reaches them — with the NaN coords above, that produces NaN-loud failure rather than silent garbage.

But ONLY `ComputeSAFSParams<StressSource>` and `ComputeSAFSParams(StressField3D&, ...)` lazy-call `ComputeBP5Params()` (via `if (a_values_.Size() != num_fault_dofs_)`).  These are the two SAFS-mode entry points.  Several other public FaultGeometry methods read `a_values_` / `eta_values_` DIRECTLY without lazy-init:

| Method                                              | Reads                                          | After new ctor (size 0) |
|-----------------------------------------------------|------------------------------------------------|-------------------------|
| `Print(os)` (line 703-740)                          | `a_values_(i)` for every i; `a_values_.Min()`; `a_values_.Max()`; `eta_values_(0)` | OUT-OF-BOUNDS abort, then NaN garbage |
| `IsVelocityWeakening(dof_idx)` (line 648-652)       | `a_values_(dof_idx)`                          | OUT-OF-BOUNDS abort     |
| `GetAValues()` (line 432)                           | `a_values_` directly                          | Empty Vector            |
| `GetEtaValues()` (line 435)                         | `eta_values_` directly                        | Empty Vector            |
| `GetTauPre()` (line 437)                            | `tau_pre_` directly                           | Empty Vector            |

If any driver or test or diagnostic calls `geom.Print(mfem::out)` on a FaultGeometry built via the new BP5 ctor before SAFS init has run (or even after, if SAFS init populated dof_data but not the FaultGeometry's analytic arrays), it crashes.  `Print()` is a standard SEAS diagnostic — every TPV/BP5 driver calls it as part of startup logging.

After `ComputeSAFSParams<StressSource>` lazy-triggers `ComputeBP5Params()` with NaN `coords_x2_/coords_x3_`:
- `a_values_(i) = bp5_params_.a_of_x2_x3(NaN, NaN)` → NaN
- `dc_values_(i) = bp5_params_.Dc_of_x2_x3(NaN, NaN)` → NaN
- `eta_values_(i) = bp5_params_.eta()` → constant, NOT NaN (eta doesn't depend on x2/x3)
- `tau_pre_(2*i), tau_pre_(2*i+1)` → NaN (then overwritten by SAFS projection)
- `V_init_vec_(2*i), V_init_vec_(2*i+1)` → NaN

`Print()` then runs without crash but emits `vw_count=0` (NaN < b is always false), `a_min=NaN, a_max=NaN` — silent garbage in the diagnostic.

`IsVelocityWeakening(i)` returns false for all i (NaN < b always false) — every DOF is classified as velocity-strengthening, which is the OPPOSITE of what a SAFS rupture run intends.

**Trigger:**
Any code that reads `a_values_` or `eta_values_` directly off a FaultGeometry built via the new BP5 ctor.  This includes:
1. `geom.Print(mfem::out)` — standard diagnostic (CRASH if called before ComputeSAFSParams; SILENT GARBAGE after).
2. `geom.IsVelocityWeakening(i)` — SAFS analyses may want to know VW vs VS distribution.
3. Direct read of `GetAValues()` / `GetEtaValues()` by any downstream postprocessor.

**Actual behavior:**
CRASH if read before lazy init; SILENT NaN garbage if read after.

**Expected behavior:**
Either (a) keep these readers consistent with the new ctor's "BP5 analytics not meaningful" intent by aborting with a clear message, OR (b) populate `a_values_` etc. with NaN sentinels eagerly in the new ctor so the failure mode is one consistent NaN, never an out-of-bounds crash.

**Suggested fix:**
Eagerly NaN-fill the BP5 analytic arrays in the new ctor (mirroring the coords_x2_/coords_x3_ NaN approach):

```diff
@@ fault_geometry.hpp new BP5 ctor body @@
       // Skip ComputeBP5Params() in the new ctor: the SAFS driver
       // overwrites per-DOF analytic values via InitializeFaultDOFs_
       // Spatial (Phase 5c) before any consumption.  ComputeSAFSParams
       // (templated overload) will lazy-initialise the BP5 arrays if
       // any caller reaches them — with the NaN coords above, that
       // produces NaN-loud failure rather than silent garbage.
+      // Eagerly size + NaN-fill the BP5 analytic per-DOF arrays so
+      // direct readers (Print, IsVelocityWeakening, GetAValues, ...)
+      // see a consistent NaN-loud state instead of an out-of-bounds
+      // crash on the empty vectors.
+      a_values_.SetSize(N);
+      eta_values_.SetSize(N);
+      dc_values_.SetSize(N);
+      tau_pre_.SetSize(2 * N);
+      V_init_vec_.SetSize(2 * N);
+      for (int i = 0; i < N; ++i)
+      {
+         a_values_(i)  = k_nan;
+         eta_values_(i) = k_nan;
+         dc_values_(i) = k_nan;
+         tau_pre_(2 * i + 0)     = k_nan;
+         tau_pre_(2 * i + 1)     = k_nan;
+         V_init_vec_(2 * i + 0)  = k_nan;
+         V_init_vec_(2 * i + 1)  = k_nan;
+      }
```

Then update `ComputeSAFSParams<StressSource>` so it does NOT lazy-call `ComputeBP5Params()` when the new ctor has run (it would overwrite the NaN with BP5(NaN, NaN) = NaN, which is the same value).  Actually the existing lazy check `if (a_values_.Size() != num_fault_dofs_)` is false after the eager NaN-fill, so the lazy init won't fire.  Good — keep that path as-is.

Additionally, document that `Print()` and `IsVelocityWeakening()` are not meaningful on a new-ctor-built FaultGeometry — emit a NaN-aware diagnostic line.  For belt-and-suspenders, optionally make `Print()` check `std::isnan(a_values_(0))` and emit "BP5 analytic state not populated (new SAFS ctor)" instead of garbage.

**Test case:**
```cpp
// tests/unit/test_spatial_setup.cpp — add S-9
static void S_9_print_does_not_crash_on_new_ctor()
{
   // Build a new-ctor FaultGeometry with 1 DOF, then call Print() into
   // a stringstream.  Without the eager NaN-fill, this throws (or
   // aborts) on a_values_(0) out-of-bounds.  With the fix, it emits a
   // line containing "nan" or "BP5 analytic not populated".
   Vector dof_coords_3d(3);  dof_coords_3d = 0.0;
   DenseMatrix dof_basis(9, 1);  dof_basis = 0.0;
   dof_basis(0, 0) = 1.0;  dof_basis(4, 0) = 1.0;  dof_basis(8, 0) = 1.0;
   Array<int> dof_to_elem(1);  dof_to_elem = 0;
   BP5Params params;
   FaultGeometry<mfem::Mesh> geom(params, dof_coords_3d, dof_basis,
                                  dof_to_elem, 1);

   std::ostringstream oss;
   geom.Print(oss);     // pre-fix: out-of-bounds abort; post-fix: emits NaN
   const std::string out = oss.str();
   TEST_ASSERT(!out.empty(), "Print() did not crash");
}
```

---

### [R-403] [MODERATE] [wave_operator.inl:3636-3643, 4549-4557] — `EvaluateADER_LSW_ForcedRupture` signature has 7 parameters; manual inspection of the new dispatch arms suggests the order matches the .hpp declaration, but no unit test exercises the actual dispatch path

**Category:** EDGE_CASE (covers a class of latent bugs)

**Description:**
The R-302 fix added two dispatch arms.  The call sites:

```cpp
fault_flux_->EvaluateADER_LSW_ForcedRupture(
   fdata,
   I_plus_local, I_minus_local,
   dt,
   GetTime(),
   I_imp_plus, I_imp_minus);
```

The declared signature in `fault_face_flux.hpp:392-398`:

```cpp
void EvaluateADER_LSW_ForcedRupture(DOFData &data,
                                    const real_t *I_plus,
                                    const real_t *I_minus,
                                    real_t dt,
                                    real_t t_now,
                                    real_t *I_imp_plus,
                                    real_t *I_imp_minus) const;
```

Argument count matches (7).  But — no unit test directly exercises the dispatch arm.  F-7 (added in round-3) calls `EvaluateADER_LSW_ForcedRupture` directly via `ffl.EvaluateADER_LSW_ForcedRupture(...)`.  The wave-operator dispatch path is not unit-tested; a future signature drift (e.g., reordering `dt` and `t_now`) would compile cleanly (both are `real_t`) and silently swap the values.

Compounding with R-401 (GetTime() returns 0), this means today's dispatch silently passes (..., dt, 0, ...) where one of those is intended as `t_now`.  A `dt` and `t_now` swap would mean the friction is evaluated at `t_now=dt`, which depending on the macro-step is similarly wrong.

**Trigger:**
Compile-time reorder of `dt` / `t_now` in either the signature or the call site.

**Actual behavior:**
Silent swap on real_t → real_t.

**Expected behavior:**
A test that ROUTES through the wave_operator dispatch (not the direct method call) would catch this.

**Suggested fix:**
Add a wave-operator-dispatch integration test in `test_phaseh_lsw_forced_rupture.cpp` paired with the R-401 fix (once SetTime is plumbed):

```cpp
static void F_8_wave_dispatch_forced_rupture_path()
{
   // Build minimal WaveOperator with 1 element, 1 fault face, 1 DOF.
   // Set FaultFrictionLaw::LSW_ForcedRupture.
   // wave.SetTime(0.25); wave.Mult(Q, dQdt);
   // Assert dof_data[0].slip_rate matches the value computed by a
   // direct EvaluateADER_LSW_ForcedRupture call with the same inputs
   // and t_now = 0.25.
}
```

Defensive sub-fix: add a static_assert or in-file marker that documents the parameter order.

**Test case:**
(See F-8 sketch above — full implementation requires R-401 fix first.)

---

### [R-404] [MODERATE] [fault_geometry.hpp new BP5 ctor — depths_ stays inconsistent with x3 convention] — `depths_(i) = world_z` in the new ctor; the legacy BP5 ctor populates `depths_(i) = coords_x3_(i)` which is the along-dip 2D fault coord

**Category:** EDGE_CASE / DEVIATION

**Description:**
Round-3 R-310 fix set `coords_x2_(i) = coords_x3_(i) = NaN` in the new BP5 ctor, but kept `depths_(i) = dof_coords_3d(3*i + 2)` = world z.  The comment:

```cpp
// Derive `depths_` from world z (used by ABSORBING / free-surface
// BC dispatch).  Leave `coords_x2_` / `coords_x3_` as NaN sentinels...
```

For a planar BP5 fault aligned with z-axis, the legacy ctor's `depths_(i) = coords_x3_(i)` (where coords_x3_ is the along-dip 2D fault coord) happens to equal world z because the fault IS planar in z.  For SAFS curvilinear faults, world z is NOT the along-dip coordinate — it's the depth in the box.

The semantic difference matters because `depths_` feeds into:
1. `ComputeDepthDependentParams` (BP2 only — never fires in the new BP5 ctor).
2. `FindClosestDOFToDepth(target_depth)` (line 614-640) — searches by `|depths_(i) - target_depth|`.  For SAFS, world z and along-dip are NOT the same; a caller probing for "the DOF closest to depth = -5000 m" would get the world-z-closest DOF, which is correct.  But if any other code assumes depths_ is along-dip, it would silently get the wrong DOF.
3. `Print()` reports `Depth range: [z_max/1000, z_min/1000] km`.  For SAFS this is the world-z range; reads correctly as depth km.

OK — depths_ = world z is consistent for SAFS (ABSORBING BC and depth searches both want world z).  But for callers that previously assumed depths_ = along-dip x3 (BP5 convention), the semantics has changed in the new ctor.  No abort guards this; latent risk.

**Suggested fix:**
Document the depths_ convention divergence in the new ctor's docstring, AND audit `GetDepths()` callers for any assumption that depths_ matches the BP5 along-dip x3:

```diff
@@ fault_geometry.hpp new BP5 ctor docstring @@
+/// `depths_(i)` is populated as WORLD z (matches the ABSORBING /
+/// free-surface BC dispatch convention).  For a planar BP5 fault
+/// aligned with the z-axis this is the same value the legacy ctor
+/// writes; for SAFS curvilinear faults this is the world depth, NOT
+/// the BP5 along-dip x3.  Any caller that assumed depths_ == x3
+/// will see different values here.
```

---

### [R-405] [MODERATE] [spatial_setup.hpp IP-aware overload] — Substantial body duplication vs the centroid overload; both will drift when DOFData fields are added

**Category:** QUALITY (drift risk)

**Description:**
The R-309 fix added an IP-aware `InitializeFaultDOFs_Spatial` overload.  The two overloads share ~30 lines of identical body (the loop that copies LSW fields, T_forced/t0_decay, zeroes RS slots).  Only the call to `seed_static_dof_fields` differs (IP-aware vs centroid).  Future additions to `DOFData` (e.g. when Phase H wiring lands) would require synchronised edits to BOTH overloads; a partial edit means one of the two leaves a field uninitialised.

**Suggested fix:**
Extract the per-DOF body into a helper that takes either the centroid IP or the supplied IP:

```diff
+/// Common per-DOF body for the centroid and IP-aware overloads.
+template <typename MeshT>
+inline void initialise_lsw_dof(DOFData& d,
+                               int elem,
+                               MeshT& mesh,
+                               const MaterialField& material,
+                               const mfem::IntegrationPoint& ip,
+                               real_t tau_pre_dip,
+                               real_t tau_pre_strike,
+                               real_t sigma_n_eff,
+                               const SlipWeakeningPerDOFParams& lsw,
+                               int i,
+                               real_t T_forced_s,
+                               real_t t0_decay_s)
+{
+   internal::seed_static_dof_fields<MeshT>(d, elem, mesh, material, ip,
+                                           tau_pre_dip, tau_pre_strike,
+                                           sigma_n_eff);
+   d.lsw_mu_s = lsw.mu_s(i);
+   d.lsw_mu_d = lsw.mu_d(i);
+   d.lsw_d_c  = lsw.d_c(i);
+   d.T_forced_rupture = T_forced_s;
+   d.t0_decay_forced  = t0_decay_s;
+   d.a = 0.0;  d.psi = 0.0;  d.Dc = 0.0;
+}
```

Then both overloads call `initialise_lsw_dof(...)` with their respective IP.  The centroid overload computes the centroid IP and passes it; the IP-aware overload passes `dof_ips[i]`.

---

### [R-406] [LOW] [godunov_flux_pool.cpp:round_sig] — R-312 guard at `|log_mag| > 290.0` returns `v` unrounded, breaking the dedup determinism for extreme inputs

**Category:** EDGE_CASE

**Description:**
The R-312 fix returns `v` (unrounded) when `|log_mag| > 290.0`.  This bypasses the rounding step, so two inputs that should logically dedup to the same bucket might land in different buckets:

```cpp
round_sig(1e-300, 6) == 1e-300   // unrounded
round_sig(1.000001e-300, 6) == 1.000001e-300  // also unrounded, different bucket
```

These would dedup to the same bucket under normal rounding.  In the unrounded path they DON'T.

This is a corner case that real SAFS material won't hit (Vp, Vs, rho are O(1) – O(1e10)).  But the guard converts an overflow-NaN failure mode into a determinism failure mode without comment.

**Suggested fix:**
Either (a) document the determinism trade-off in the guard comment, or (b) MFEM_VERIFY the input is in the safe range instead of silently bypassing the rounding:

```diff
-   if (!std::isfinite(log_mag) || std::abs(log_mag) > 290.0)
-   {
-      return v;
-   }
+   MFEM_VERIFY(std::isfinite(log_mag) && std::abs(log_mag) <= 290.0,
+               "GodunovFluxPool::round_sig: |log10(|v|)|=" << std::abs(log_mag)
+               << " is outside the safe rounding range; got v=" << v
+               << ".  Material values in SAFS are O(1)..O(1e10); this is "
+               "almost certainly a corrupted material input.");
```

Option (b) is strictly louder; (a) keeps the defensive bypass.  Either way, the silent dedup-bypass is the actual concern.

---

### [R-407] [LOW] [wave_operator.hpp:SetGodunovFluxPool — abort message] — Abort message says "follow-up edits" without pointing at where to land them

**Category:** QUALITY

**Description:**
The R-301 abort message:

```cpp
MFEM_ABORT("WaveOperator::SetGodunovFluxPool: per-element flux "
           "dispatch (Phase H.2), heterogeneous CFL (Phase H.3), "
           "and bi-material shared-face MPI exchange (Phase H.5) "
           "are NOT yet wired in wave_operator.inl.  Setting this "
           "pool would be a no-op that silently runs scalar-"
           "material physics on heterogeneous input.  Implement "
           "the dispatch follow-up (R-301) before re-enabling.");
```

The "R-301" reference points to an internal review document the user isn't going to read.  A better message would name the specific function bodies in wave_operator.inl that need a per-element flux dispatch added, so a future implementer doesn't have to re-derive the audit:

```diff
-              "Implement the dispatch follow-up (R-301) before re-enabling.");
+              "To enable: (1) add `flux_pool_->At(e)` branches to "
+              "every `flux_.Interior(...)` call in wave_operator.inl "
+              "(currently only `flux_` is used); (2) extend "
+              "ComputeMaxDt() to walk per-element c_p; (3) implement "
+              "ExchangeBiMaterialNeighbours_ for shared-face MPI.  "
+              "Then remove this abort.");
```

---

## Summary

- Critical issues: **2**  (R-401 GetTime() stale → forced rupture never fires; R-402 new BP5 ctor leaves a_values_ size 0 → Print/IsVelocityWeakening crash or NaN)
- Moderate issues: **3**  (R-403 dispatch never integration-tested; R-404 depths_ semantic divergence; R-405 IP-aware overload duplication risk)
- Low issues: **2**  (R-406 dedup determinism bypass; R-407 abort message clarity)
- Plan compliance: **PARTIAL** — Round-3 fixed R-301..R-313 as documented, but R-302's dispatch-arm fix is functionally inert because GetTime() returns 0 in production (R-401), and R-310's "skip ComputeBP5Params" fix breaks Print/IsVelocityWeakening on a new-ctor FaultGeometry (R-402).
- Verdict: **FAIL — must fix R-401 and R-402 before any SAFS dynamic-rupture run can produce correct physics or pass a Print() call.**  Both are silent at compile time AND silent at unit-test time (every existing test calls the new method directly or uses the centroid overload after sufficient init), so the round-3 test sweep (515/515) gives false confidence.  R-401 in particular is the worst kind of bug: the dispatch arm LOOKS connected, the unit test F-7 LOOKS pass, but the wave-operator-routed path is dead.

## Unreviewed Areas

- **Heterogeneous WaveOperator code** — R-301 abort means Phase H.2/H.3/H.5 dispatch + CFL + MPI exchange are still un-implemented; nothing new to review.
- **Phase 4 driver source** — still not present.
- **Phase 6 sbatches + verification scripts** — only `safs_smoke_8N_400r_dev.sbatch` + `verify_constant_tensor_projection.py` exist (Phase 3b helper); the plan-mandated `verify_spatial_dyn_smoke_safs.py` and `compare_spatial_dyn_velocity_models.py` are still missing.
- **MPI behaviour of the new BP5 ctor** — `ComputeGatherInfo()` is called when mpi_ctx is provided, but the rest of the code path was reviewed serial only.  An np=4 SAFS run could expose race conditions in the new dof_to_elem_new_ctor_ array distribution that the round-3 sweep didn't exercise.
- **Round-2 R-201 / R-202 / R-203 regressions** — confirmed passing in the round-3 fix report's test sweep (28/28 parser tests).  Not re-audited in round 4.
- **TPV byte-exact full regression** — only `seas_test_tpv104_checkpoint` (164/164) and `seas_test_compute_safs_params` (13/13) were run in round 3.  Full `make test-tpv104` / `make test-tpv205` / `make test-bp5-*` byte-exact gates should run before Phase H merges to confirm the wave_operator.inl dispatch arm additions did not subtly perturb the TPV205 LSW path (which uses the existing `LSW` arm — the new arm shouldn't fire, but the surrounding restructuring deserves a fresh diff).
