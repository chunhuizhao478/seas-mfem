# Code Review: TPV26/27 Phase 2 (forced rupture) + Phase 3 (depth cohesion) — 2026-07-09

> **File placement note.** `miniapps/seas/REVIEW.md` already exists with unrelated prior
> content, and this repo's convention is dated topical review files
> (`REVIEW_cross_rank_phase1_code_2026-06-07.md`, …).  This review is therefore written
> here rather than clobbering `REVIEW.md`.  **The `/code-fix` agent should read THIS file.**

## Review Scope
- **Plan:** `miniapps/seas/tpv26/PLAN_TPV26_27_spatial_dyn_driver_2026-06-28.md` §5 Phase 2 + Phase 3
- **Files reviewed:**
  - `dynamic/spatial_nucleation.{hpp,cpp}` (`ForcedRuptureSpec`, `ForcedRuptureTime`, `ResolveForcedRupture`)
  - `spatial/code/spatial_friction.{hpp,cpp}` (`NucleationKind::ForcedRupture`, `[nucleation.forced_rupture]` parser, `tpv2627_depth` validation)
  - `spatial/code/spatial_stress.{hpp,cpp}` (Phase 1 `Tpv2627DepthStressSource`, re-checked)
  - `dynamic/nucleation_factory.cpp`
  - `dynamic/friction_substep_iterator.{hpp,cpp}` (round-6 fix)
  - `dynamic/friction_iterator_factory.cpp`
  - `drivers/spatial_dyn_driver.cpp`
  - `Makefile`
  - `tests/unit/test_forced_rupture_resolver.cpp`, `test_forced_rupture_iterator_parity.cpp`, `test_tpv2627_cohesion_taper.cpp`
  - `tpv26/configs/tpv26_spatial_{forcedrupture,depthstress}_smoke.toml`
- **Domain context consulted:** `miniapps/seas/CLAUDE.md` (sign conventions, byte-exact contracts,
  Gmsh v2.2), `dynamic/wave_operator.inl` (seam fault dispatch, `FaultFrictionLaw`),
  `io/tpv104_checkpoint.hpp` (serialized DOFData fields), implementer completion report.

**Verified-not-broken (no finding):**
- `io/tpv104_checkpoint.hpp` serializes only `V1 V2 psi sigma_n_nuc slip1 slip2 slip_rate
  tau1_nuc tau2_nuc` — **not** `T_forced_rupture`/`t0_decay_forced`.  `InitializeFaultDOFs_Spatial`
  (driver `:2058`) runs *before* `ReadTpv104Checkpoint` (`:3016`), so `T(r)` is deterministically
  rebuilt from config on restart.  The plan's "no checkpoint change" claim holds.
- `StaticOverstress::ApplyIncrement/ApplyAbsolute` are true no-ops, so routing
  `ForcedRupture -> StaticOverstress` in `MakeNucleation` cannot perturb `tau_nuc`.
- The driver's post-`MakeNucleation` `dynamic_cast` chain has no trailing `else`-abort, so
  `StaticOverstress` falls through harmlessly.
- `seas_test_friction_substep_iterator_parity` (the pre-existing byte-exactness guard for the
  file I edited) still passes **36/36** — the plain-LSW `StepOneQP_` path is unchanged.

---

## Findings

### [R-001] CRITICAL [spatial_friction.cpp:parse_root, `tpv2627_depth` branch] — Pore-pressure gradient is never cross-checked against `water_density*g`, silently producing a fault that can never slip

**Category:** BUG

**Description:**
`Tpv2627DepthStressSource::Evaluate` builds `sigma11/sigma33/sigma13` from the **effective**
vertical stress `(sigma22 + Pf)` with `Pf = water_density * g * depth` (values taken from the
`[stress]` block).  The pore pressure that is actually *subtracted* from the projected normal
stress is a **separate, unrelated config knob**: `ComputeParams` subtracts
`P_p = [pore_pressure].P_p_grad_pa_per_m * depth`.

Nothing requires these to agree.  When they disagree, `sigma_n` is neither the total nor the
effective normal stress, and the on-fault ratio `tau/sigma_n` silently leaves the `(mu_d, mu_s)`
window that the entire TPV26/27 nucleation mechanism depends on.

**Trigger:**
A `tpv2627_depth` config that omits `[pore_pressure].P_p_grad_pa_per_m` (it defaults to `0.0`) or
sets it to any value != `water_density * g`.  The Phase-0 smoke config in this very repo has
`P_p_grad_pa_per_m = 0.0`; only `tpv26_spatial_depthstress_smoke.toml` happens to set `9800.0`.

**Actual behavior** (computed, at 10 km depth):

| `P_p_grad_pa_per_m` | `sigma_n_eff` | `tau_strike` | ratio  | outcome |
|---|---|---|---|---|
| `9800` (= wd*g)     | 175.64 MPa | 27.66 MPa | **0.1575** | `mu_d < r < mu_s` -> nucleates correctly |
| `0` (omitted)       | 273.64 MPa | 27.66 MPa | **0.1011** | `r < mu_d = 0.12` -> **fault can NEVER slip** |

No warning, no abort. The run completes with `V_max = 0` and looks "stable".

**Expected behavior:**
Abort at config-parse time with a message naming the required value.

**Suggested fix** — in `spatial/code/spatial_friction.cpp`, at the end of the
`else if (cfg.stress.kind == StressSourceKind::Tpv2627Depth)` branch (after the
`omega_bot_m > omega_top_m` verify).  `[pore_pressure]` is a required block and is parsed
earlier (`:894-903`), so `cfg.stress.pore_pressure` is populated here.

```diff
          MFEM_VERIFY(d.omega_bot_m > d.omega_top_m,
                      "[stress] kind=\"tpv2627_depth\" omega_bot_m ("
                      << d.omega_bot_m << ") must be > omega_top_m ("
                      << d.omega_top_m << ")");
+         // The depth profile builds sigma11/sigma33/sigma13 from the EFFECTIVE
+         // vertical stress (sigma22 + Pf), Pf = water_density*g*depth.
+         // ComputeParams then subtracts P_p = P_p_grad*depth.  If the two
+         // disagree, sigma_n is neither total nor effective and tau/sigma_n
+         // silently leaves (mu_d, mu_s) -> the fault can never slip.
+         {
+            const real_t pf_grad = d.water_density * d.g;
+            const real_t got = cfg.stress.pore_pressure.P_p_grad_pa_per_m;
+            const real_t tol = 1.0e-6 * std::max(static_cast<real_t>(1.0), pf_grad);
+            MFEM_VERIFY(std::abs(got - pf_grad) <= tol,
+                        "[stress] kind=\"tpv2627_depth\" requires "
+                        "[pore_pressure].P_p_grad_pa_per_m == water_density*g = "
+                        << pf_grad << " Pa/m; got " << got << ".  Otherwise the "
+                        "reported sigma_n is inconsistent with the depth profile "
+                        "and tau/sigma_n leaves (mu_d, mu_s) — the fault would "
+                        "never slip.");
+         }
```

**Test case** (extend `tests/unit/test_tpv2627_stress_source.cpp`):
```cpp
// Documents WHY R-001's guard exists: with P_p_grad = 0 the ratio leaves (mu_d, mu_s).
void Test_R001_Pf_Grad_Must_Match_WaterDensity_g()
{
   auto src = SpecSource();
   const real_t depth = 10000.0;
   const Proj p = project(src.Evaluate(-5000.0, 0.0, -depth));

   const real_t ratio_ok  = p.tau_strike / (p.sigma_n_total - 9800.0 * depth);
   const real_t ratio_bad = p.tau_strike / (p.sigma_n_total - 0.0 * depth);

   TEST_ASSERT(0.12 < ratio_ok && ratio_ok < 0.18,
               "P_p_grad = water_density*g -> ratio inside (mu_d, mu_s)");
   TEST_ASSERT(ratio_bad < 0.12,
               "P_p_grad = 0 -> ratio BELOW mu_d: fault can never slip "
               "(this is what the parser guard must reject)");
}
```

---

### [R-002] MODERATE [friction_substep_iterator.cpp:StepOneQP_ vs wave_operator.inl seam dispatch] — Interior and seam evaluate mu(delta,t) at different times, violating the plan's explicit seam-consistency requirement

**Category:** DEVIATION

**Description:**
Plan §5 Phase 2, requirement 2, bullet 4 states verbatim:

> *"Interior (iterator) and seam (`fault_face_flux.cpp:1088-1094`, which already calls the same
> helper) must use identical `μ(δ,t)` to avoid an MPI seam discontinuity."*

They do **not**:
- **Interior** (`StepOneQP_`) evaluates `mu` at `t_sub_end = t_macro_start + k*dt/O`, `k = 1..O`.
- **Seam** (`wave_operator.inl:4345` and `:5318`) calls
  `EvaluateADER_LSW_ForcedRupture(..., dt, GetTime(), ...)`, and `wave.SetTime(t)` is called
  **once per macro step** (`drivers/spatial_dyn_driver.cpp:3395`), so `GetTime()` is the
  macro-step **start** time, constant across all O sub-steps.

The two therefore disagree by up to one macro `dt`. Because `f_2(t) = (t - T)/t0` on the ramp,
the seam friction lags the interior by up to `Δmu ≈ (dt / t0_s) * (mu_s - mu_d)`.

**Trigger:**
Any forced-rupture run with fault faces split across MPI ranks (`shared > 0`). The smokes in
this branch all report `shared = 0` (the fault-locality partition keeps the tiny fault on one
rank), which is exactly why the implementer's np=1/2/4 runs agreed to the last digit — the seam
path was **never exercised**. Production TPV26/27 (100 m / 50 m, many ranks) will exercise it.

**Actual behavior:**
`mu_seam = mu(delta, t_macro_start)`, `mu_interior = mu(delta, t_macro_start + k*dt/O)`.
For the production 50 m mesh (`dt ≈ 4e-4 s`) and `t0_s = 0.5 s`:
`Δmu ≲ 4e-4/0.5 * 0.06 ≈ 4.8e-5` — small, but a systematic one-macro-step lag on every seam QP,
i.e. precisely the discontinuity the plan forbade.

**Expected behavior:**
Either both paths evaluate `mu` at the same reference time, or the residual mismatch is bounded
and documented.

**Suggested fix (mechanical, Option A — bound + document).** The seam ADER solve is structurally
*one* Riemann solve over the full macro `dt`, so it cannot be given sub-step resolution without
reworking the wave operator; forcing the interior to `t_macro_start` instead would discard the
sub-step fidelity the plan explicitly asked for ("thread the substep absolute end-time").
Therefore bound the mismatch at setup and document the convention.

In `drivers/spatial_dyn_driver.cpp`, immediately after `is_forced_rupture` is defined:
```diff
    const bool is_forced_rupture = wants_forced_rupture && is_lsw;
+   // R-002: seam QPs evaluate mu at GetTime() (macro-step start, wave op sets
+   // time once per step) while interior QPs evaluate at t_sub_end.  The
+   // resulting seam lag is bounded by (dt / t0_s) * (mu_s - mu_d).  Require a
+   // forced-rupture ramp that is long compared to the step so the lag stays
+   // negligible; a step ramp (t0_s == 0) would put a full mu_s->mu_d jump on
+   // the seam and is rejected for MPI runs.
+   if (is_forced_rupture && nprocs > 1)
+   {
+      MFEM_VERIFY(cfg.nucleation.forced_rupture.t0_s > 0.0,
+                  "[nucleation.forced_rupture].t0_s must be > 0 for MPI runs: "
+                  "seam QPs lag interior QPs by one macro step, so a step ramp "
+                  "(t0_s = 0) produces a full mu_s->mu_d seam discontinuity.");
+   }
```
and, at the top of `LinearSlipWeakeningIterator::StepOneQP_`'s forced-rupture branch, replace the
"no MPI seam discontinuity" claim in `friction_substep_iterator.hpp:WaveOpLaw()` doc-comment with
the accurate statement (the seam uses `GetTime()`; the lag is `O(dt/t0)`).

**Option B (needs a decision, do NOT apply mechanically):** pass `t_macro_start` instead of
`t_sub_end` for exact seam parity, sacrificing sub-step time resolution.

**Test case** (new, `tests/unit/test_forced_rupture_iterator_parity.cpp`):
```cpp
// R-002: bound the seam/interior mu mismatch.  Seam evaluates at t_macro_start,
// interior at t_sub_end; the gap must stay small relative to the ramp.
void Test_R002_Seam_Interior_Mu_Lag_Bounded()
{
   const real_t T = 1.0, t0 = 0.5, dt = 4.0e-4;
   const real_t mu_seam     = spatial::LSWFrictionCoefficient_ForcedRupture(
                                 0.0, kMuS, kMuD, kDc, /*t=*/1.2,      T, t0);
   const real_t mu_interior = spatial::LSWFrictionCoefficient_ForcedRupture(
                                 0.0, kMuS, kMuD, kDc, /*t=*/1.2 + dt, T, t0);
   const real_t bound = (dt / t0) * (kMuS - kMuD);
   TEST_ASSERT(std::abs(mu_interior - mu_seam) <= bound + 1e-15,
               "seam/interior mu lag bounded by (dt/t0)*(mu_s-mu_d)");
}
```

---

### [R-003] MODERATE [spatial_dyn_driver.cpp] — `forced_rupture` + a non-ADER integrator aborts deep in the time loop instead of at config time

**Category:** EDGE_CASE

**Description:**
`FaultFrictionLaw::LSW_ForcedRupture` has no instantaneous solve on the Mult/RK path.
`wave_operator.inl` correctly `MFEM_ABORT`s (interior at `:3031`, shared at `:3701`,
*"Use --time-integrator ader for forced-rupture configs."*), so there is **no silent wrong
physics**. But the abort fires only once the first RK stage evaluates the fault — after mesh
load, partition, stress projection, friction resolve and DOF init. On a production mesh that is
minutes of wasted work, and the error surfaces from deep inside the wave operator.

**Trigger:** `[nucleation] kind="forced_rupture"` together with `--time-integrator rk4` (or `rk45`).

**Actual behavior:** full setup, then `MFEM_ABORT` from `WaveOperator::Mult`.

**Expected behavior:** abort during config validation, next to the existing
`forced_rupture requires law="slip_weakening"` guard.

**Suggested fix** — `drivers/spatial_dyn_driver.cpp`, right after `is_forced_rupture`:
```diff
    const bool is_forced_rupture = wants_forced_rupture && is_lsw;
+   // The LSW_ForcedRupture flux dispatch exists only on the ADER path
+   // (EvaluateADER_LSW_ForcedRupture); the Mult/RK path aborts.  Fail fast.
+   MFEM_VERIFY(!is_forced_rupture
+               || cfg.numerics.time_integrator == spatial::TimeIntegratorKind::ADER,
+               "spatial_dyn_driver: [nucleation] kind=\"forced_rupture\" requires "
+               "the ADER time integrator — LSW_ForcedRupture has no instantaneous "
+               "solve on the Mult/RK path.  Use --time-integrator ader.");
```

**Test case** (extend `tests/unit/test_tpv2627_cohesion_taper.cpp` or a config test — parse-level,
since the driver guard is not unit-reachable):
```cpp
// R-003: a forced_rupture config must declare the ADER integrator.
void Test_R003_ForcedRupture_Requires_Ader()
{
   const SpatialFrictionConfig cfg = ParseSpatialFrictionConfigString(kForcedRuptureToml);
   TEST_ASSERT(cfg.nucleation.kind == NucleationKind::ForcedRupture, "kind parsed");
   TEST_ASSERT(cfg.numerics.time_integrator == TimeIntegratorKind::ADER,
               "forced_rupture configs must select ADER (driver enforces this)");
}
```

---

### [R-004] MODERATE [spatial_nucleation.cpp:ResolveForcedRupture] — `hypocenter_y_m` is parsed, stored, and silently ignored; `r` is wrong for any non-planar fault

**Category:** ASSUMPTION

**Description:**
`r = sqrt((x - x_hyp)^2 + (z - z_hyp)^2)` uses only `(x, z)`. `spec.hypocenter_y_m` is parsed,
documented as "stored for completeness", and **never read**. This is exact for the planar,
vertical `y = 0` TPV26/27 fault, but the *same driver* serves curved SAFS faults, and the three
sibling resolvers (`ResolveGradualOverstress`, `…CompactCircular`, `…InstantaneousCircular`) all
measure the in-fault-plane radius through the per-DOF `(dip, strike)` basis precisely because of
this. Nothing prevents `kind="forced_rupture"` on a curved fault, where the `(x, z)` measure
silently mis-places the entire forced-rupture front.

**Trigger:** `[nucleation] kind="forced_rupture"` on any fault that is not the planar `y = 0`
plane (e.g. a SAFS mesh), or any config that sets `hypocenter_y_m != 0` and expects it to matter.

**Actual behavior:** `r` ignores `y` entirely; a user-supplied `hypocenter_y_m` has zero effect.

**Expected behavior:** either reject non-planar use, or measure `r` in the fault plane via the
per-DOF basis like the sibling resolvers.

**Suggested fix (minimal, mechanical):** validate planarity inside `ResolveForcedRupture`, which
already receives the DOF coordinates.
```diff
    const int N = dof_coords_3d.Size() / 3;
    p.T_forced_s.SetSize(N);
    p.t0_decay_s.SetSize(N);
 
+   // r is measured in the (x, z) plane, which is exact ONLY for the planar,
+   // vertical y = 0 fault of TPV26/27.  Reject curved / offset faults rather
+   // than silently mis-placing the forced-rupture front.  (The sibling
+   // resolvers use the per-DOF dip/strike basis for the general case.)
+   for (int i = 0; i < N; ++i)
+   {
+      const real_t dy = dof_coords_3d(3 * i + 1) - spec.hypocenter_y_m;
+      MFEM_VERIFY(std::abs(dy) <= 1.0,
+                  "ResolveForcedRupture: fault DOF " << i << " is " << dy
+                  << " m off the hypocenter's fault-normal plane.  The (x, z) "
+                  "radius measure is valid only for the planar vertical y=0 "
+                  "TPV26/27 fault; use a basis-projected radius for curved faults.");
+   }
+
    for (int i = 0; i < N; ++i)
```

**Test case** (extend `tests/unit/test_forced_rupture_resolver.cpp`; replaces the current
`Test_Hypocenter_Offset_And_Plane` assertion that a y-offset is *ignored*, which currently
**enshrines the bug**):
```cpp
// R-004: a fault DOF far off the y=0 plane must be rejected, not silently
// projected.  (Today it is silently accepted with r computed from (x,z).)
void Test_R004_NonPlanar_Fault_Rejected()
{
   ForcedRuptureSpec spec = SpecDefaults();
   Vector coords(3);
   coords(0) = 0.0; coords(1) = 3000.0; coords(2) = 0.0;  // 3 km off-plane
   // EXPECT: ResolveForcedRupture aborts (MFEM_VERIFY).  Run as a death test,
   // or assert the planarity predicate directly once it is factored out.
}
```

---

### [R-005] MODERATE [tests] — The round-6 wiring itself is untested; a `dt_sub`/`t_sub_end` transposition would pass every existing test

**Category:** BUG (missing test coverage of the change under review)

**Description:**
`test_forced_rupture_iterator_parity.cpp` — despite its name — never constructs a
`LinearSlipWeakeningIterator`. It only exercises the two free friction-coefficient functions.
Consequently **none** of the three things the round-6 fix actually changed are covered:
1. that `StepOneQP_` receives `t_sub_end` and **not** `dt_sub` (the two are adjacent parameters
   of the same type `real_t`, and `dt_sub` is passed immediately before it — a classic
   transposition site);
2. that `MakeFrictionIterator` sets `forced_rupture_` from `cfg.nucleation.kind`;
3. that `WaveOpLaw()` flips to `LSW_ForcedRupture` (the driver asserts this at `:3034`, but no
   unit test does).

Confirmed: `grep -n ForcedRupture tests/unit/test_friction_iterator_factory.cpp
tests/unit/test_advance_interface_compiles.cpp` returns nothing.

A transposition (`StepOneQP_(..., dt_sub, dt_sub, ...)`) would make the forced front fire ~`1/dt`
times too late, and every current unit test would still pass. The end-to-end smoke happens to
catch it, but smokes are not run in `make test`.

**Trigger:** any future edit to `StepOneQP_`'s argument list or the factory.

**Suggested fix (a):** cover the factory + `WaveOpLaw()` — `tests/unit/test_friction_iterator_factory.cpp`:
```cpp
// R-005: forced_rupture cfg must produce an iterator that reports
// WaveOpLaw() == LSW_ForcedRupture (the driver asserts agreement).
void Test_R005_ForcedRupture_Factory_Flips_WaveOpLaw()
{
   spatial::SpatialFrictionConfig cfg = MakeLswConfig();     // existing helper
   cfg.nucleation.enabled = true;
   cfg.nucleation.kind    = spatial::NucleationKind::ForcedRupture;
   auto it = MakeFrictionIterator(cfg, flux, nullptr);
   TEST_ASSERT(it->WaveOpLaw() == FaultFrictionLaw::LSW_ForcedRupture,
               "forced_rupture -> WaveOpLaw() == LSW_ForcedRupture");

   cfg.nucleation.kind = spatial::NucleationKind::GradualOverstress;
   auto it2 = MakeFrictionIterator(cfg, flux, nullptr);
   TEST_ASSERT(it2->WaveOpLaw() == FaultFrictionLaw::LSW,
               "non-forced LSW -> WaveOpLaw() == LSW (byte-exact gate)");
}
```

**Suggested fix (b):** cover the `t_sub_end` plumbing by reusing the existing fixture in
`tests/unit/test_friction_substep_iterator_parity.cpp`: drive `Advance()` with
`t_macro_start = 0`, `dt_macro` and `O` sub-steps, one DOF with
`T_forced_rupture = 0.5 * dt_macro`, `t0_decay_forced = 0`. With `t_sub_end` the DOF weakens on
the *later* sub-steps; with `dt_sub` (each `= dt_macro/O < T_forced`) it never weakens. Assert
`d.slip2` advanced (or `mu` reached `mu_d`).

---

### [R-006] LOW [spatial_nucleation.cpp:ForcedRuptureTime] — Public API validated with `MFEM_ASSERT`, which compiles out in release builds

**Category:** BUG

**Description:**
`ForcedRuptureTime` is declared in the header and called directly by the unit test, but guards
`rcrit > 0`, `vs > 0`, `vr_factor > 0` with `MFEM_ASSERT` (a no-op in optimized builds).
`ResolveForcedRupture` re-checks with `MFEM_VERIFY`, so the *driver* path is safe; a direct
caller in a release build is not. (Impact is limited because the `!(r < rcrit)` early-out and
the `isfinite`/`> 1e9` clamp together absorb `rcrit = 0` and `vs = 0` into the `1e9` sentinel —
hence LOW, not MODERATE.)

**Suggested fix** — `dynamic/spatial_nucleation.cpp`:
```diff
-   MFEM_ASSERT(rcrit > 0.0, "ForcedRuptureTime: rcrit must be > 0");
-   MFEM_ASSERT(vs > 0.0, "ForcedRuptureTime: vs must be > 0");
-   MFEM_ASSERT(vr_factor > 0.0, "ForcedRuptureTime: vr_factor must be > 0");
+   MFEM_VERIFY(rcrit > 0.0, "ForcedRuptureTime: rcrit must be > 0; got " << rcrit);
+   MFEM_VERIFY(vs > 0.0, "ForcedRuptureTime: vs must be > 0; got " << vs);
+   MFEM_VERIFY(vr_factor > 0.0,
+               "ForcedRuptureTime: vr_factor must be > 0; got " << vr_factor);
```

---

### [R-007] LOW [spatial_stress.cpp:Tpv2627DepthStressSource ctor] — `water_density` is unvalidated; a negative value silently inverts the pore pressure

**Category:** EDGE_CASE

**Description:**
The ctor validates `rho > 0`, `g > 0`, `omega_bot > omega_top`, but not `water_density`. A
negative `water_density` yields `Pf < 0`, which flips the sign of the effective vertical stress
`(sigma22 + Pf)` and therefore of `sigma13` — producing a left-lateral background with no
diagnostic. (`b11/b33/b13` are legitimately free-signed, so they need no bound.)

**Suggested fix** — `spatial/code/spatial_stress.cpp`:
```diff
    MFEM_VERIFY(g_ > 0.0,
                "Tpv2627DepthStressSource: g must be > 0; got " << g_);
+   MFEM_VERIFY(wd_ >= 0.0,
+               "Tpv2627DepthStressSource: water_density must be >= 0 (0 disables "
+               "pore pressure); got " << wd_);
    MFEM_VERIFY(omega_bot_ > omega_top_,
```

**Test case** (extend `tests/unit/test_tpv2627_stress_source.cpp`): assert that with
`water_density = 0` the tensor is the *total*-stress profile (`Pf = 0`), i.e. `tau_strike` still
> 0 and `sigma_n_total = -sigma33` — pinning that `wd = 0` is the sanctioned "no pore pressure"
setting while negatives are rejected.

---

### [R-008] LOW [tpv27/configs/*] — Phase 3 was applied to only one of the two config families the plan names

**Category:** DEVIATION

**Description:**
Plan §5 Phase 3, *Files to modify*: **"The tpv26/tpv27 config TOMLs only."** The cohesion taper
rule was added to `tpv26/configs/tpv26_spatial_depthstress_smoke.toml` only.
`tpv27/configs/tpv27_spatial_smoke.toml` still carries no `[[friction.slip_weakening.spatial]]`
rule, so a TPV27 run has constant cohesion. Phase 3's acceptance criteria are met numerically by
the unit test, but the plan's file-level requirement is partially unmet.

**Suggested fix:** either add the identical taper rule to the tpv27 config, or record an explicit
deferral note in the plan (the production tpv26/tpv27 configs are Phase 4 deliverables).
```diff
+[[friction.slip_weakening.spatial]]
+kind                   = "depth"
+cohesion_floor_pa      = 0.40e6
+cohesion_grad_pa_per_m = 720.0
+cohesion_ref_depth_m   = 5000.0
+cohesion_taper_axis    = "z"
```

---

## Summary
- Critical issues: **1** (R-001)
- Moderate issues: **4** (R-002, R-003, R-004, R-005)
- Low issues: **3** (R-006, R-007, R-008)
- Plan compliance: **PARTIAL** — Phase 2 requirement "interior and seam must use identical
  mu(delta,t)" is violated (R-002); Phase 3 file-level scope partially unmet (R-008). All other
  plan requirements are implemented as specified, and the two documented `[CORR]` deviations
  (`LSW_ForcedRupture` already existed; taper keys live on `SpatialRule`) are correct and
  independently verified.
- **Verdict: PASS WITH FIXES.** R-001 must be fixed before any TPV26/27 production run — it
  silently produces a fault that cannot slip. R-002 must be resolved (or explicitly accepted)
  before any multi-rank production run, since local smokes never exercised the seam path.

## Fix Status (applied 2026-07-09, post-review)

| ID | Severity | Status | Where |
|---|---|---|---|
| R-001 | CRITICAL | **FIXED** + test | `spatial_friction.cpp` parser guard; verified it aborts on `P_p_grad=0` |
| R-002 | MODERATE | **RESOLVED (superseded)** | Dissolved by the unify plan (PLAN_unify_interior_shared_fault_substep_2026-07-09.md) Phases 2+4: shared fault QPs now consume the iterator's substep buffer, so seam and interior evaluate mu at the identical t_sub_end — one clock, no lag. The Option-A `t0_s > 0` MPI guard was removed (Phase 4); `Test_R002_Seam_Interior_Mu_Lag_Bounded` was replaced by `Test_Seam_Interior_Mu_Identical` (bit-exact). Neither Option A nor B is needed. |
| R-003 | MODERATE | **FIXED** + test | driver early ADER guard; verified it aborts *before* fault setup |
| R-004 | MODERATE | **FIXED** + test | `ResolveForcedRupture` planarity guard; the test that enshrined the y-ignored behaviour was rewritten into a death test |
| R-005 | MODERATE | **FIXED** (2 tests) | `F8_forced_rupture_flips_waveoplaw`; `R005b_forced_rupture_uses_t_sub_end` |
| R-006 | LOW | **FIXED** | `MFEM_ASSERT` -> `MFEM_VERIFY` |
| R-007 | LOW | **FIXED** + test | `water_density >= 0` |
| R-008 | LOW | **FIXED** | taper rule added to `tpv27_spatial_smoke.toml` (see behaviour note below) |

**Behaviour change (expected, R-008):** with the spec cohesion added, `tpv27_spatial_smoke.toml`
now reports `V_max = 0` — the 4 MPa surface cohesion locks the shallow 1.5 km *smoke* fault, just
as it does for `tpv26_spatial_depthstress_smoke.toml`. The config still loads and advances 166
steps, so the Phase-0 acceptance criterion (">= 5 steps") is met.

**Post-fix verification:** `test_friction_substep_iterator_parity` 39/39 (the 36 original
byte-exact parity assertions intact), `test_friction_iterator_factory` 21/21,
`test_tpv2627_stress_source` 17/17, `test_forced_rupture_resolver` 21/21,
`test_forced_rupture_iterator_parity` 31/31, `test_tpv2627_cohesion_taper` 17/17,
`test_spatial_friction_resolver` 127/127, `test_spatial_nucleation` 73/73,
`test_nucleation_factory` 22/22, `test_tpv31_nucleation_and_cohesion` 23/23.
Baseline smoke `V_max = 0.185543` unchanged (byte-exact gate holds); forced-rupture smoke
`V_max = 0.199649` identical at np = 1 / 2 / 4.

**R-002 is CLOSED (2026-07-10)** — superseded by the unify plan Phases 2+4 (see the Fix Status
table): the substep buffer is now consumed on shared fault QPs, giving seam and interior one clock.

## Unreviewed Areas
- `dynamic/wave_operator.inl` seam dispatch was **read** (to establish R-002/R-003) but its
  interior ADER fault path was not audited — it is outside the Phase 2/3 change set.
- MPI **seam** behaviour was never executed: every local smoke reported `shared = 0` because the
  fault-locality partition keeps the small smoke fault on a single rank, and no CLI flag disables
  that partition. R-002's discontinuity is therefore established by code reading, not measurement.
- The Phase-4 production configs / meshes (40x20 km, 100 m & 50 m) do not exist yet, so no
  convergence or SCEC-reference comparison was possible.
- `make test` was not run in full (each of ~230 targets links MFEM). The directly-affected tests
  were run: `test_friction_substep_iterator_parity` 36/36, `test_spatial_friction_resolver`
  127/127, `test_friction_iterator_factory` 17/17, plus the three new tests (20/20, 28/28, 17/17).
- Two pre-existing tests (`test_constant_tensor_sign`, `test_spatial_friction_config`) abort on a
  fixture (`spatial_friction_slip_weakening_safs_projected_stress.toml`) that is untracked in git
  and absent from the worktree — unrelated to this change set.
