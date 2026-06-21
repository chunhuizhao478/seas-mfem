# Code Review: PLAN_spatial_seas_quasidynamic_driver_2026-05-31.md (plan-stage review, 2026-06-01)

## Review Scope
- Plan: `miniapps/seas/document/spatial_seas_dev/PLAN_spatial_seas_quasidynamic_driver_2026-05-31.md`
- Implementation status: **none yet** — the worktree contains only the plan (`git status`: `?? miniapps/seas/document/spatial_seas_dev/`). This is therefore an audit of the plan's *normative contract* (the plan declares "Function signatures, file paths, and acceptance criteria are normative") against the actual in-tree code the plan instructs the implementing agent to call.
- Files cross-checked (main repo, identical to worktree for these paths):
  - `fault/rate_state_fault.hpp` (RateStateFaultOperator ctor + SetSAFSMode)
  - `fault/fault_geometry.hpp` (ctors, ComputeParamsFaultLocal, SetRateStatePerDOF target members, HasParams)
  - `domain/elasticity_operator.hpp` + `..._traction.inl` + `..._setup.inl` (ctor, getters, GetFaultDOFCoords3D, RestrictToOwnedFault, owned/local DOF counts, byNODES ordering)
  - `spatial/code/spatial_friction.hpp/.cpp` (ResolveRateState signature + sigma_n semantics, RateStatePerDOFParams)
  - `solver/seas_operator.hpp` (Mult/SetInitialCondition)
  - `solver/time_stepper.hpp` (DormandPrinceRK45 API)
  - `friction/dieterich_ruina.hpp`, `friction/state_evolution.hpp` (Constants, AgingLawPsi)
  - `linalg/hypre.cpp` (SetSystemsOptions/SetElasticityOptions ordering)
- Domain context: `CLAUDE.md` (root + `miniapps/seas/`), project memory (owned-vs-local seam, stale-`.o` hazard, SAFS sign rules).

## Findings

### [R-001] CRITICAL [drivers/spatial_seas_driver.cpp / Phase 3 + Phase 5 ordering] — Prescribed call order aborts at the `RateStateFaultOperator` constructor and creates a circular σ_n dependency

**Category:** BUG (plan deviation from a hard in-code precondition)

**Description:**
The plan splits the QD wiring across two phases with this temporal order:

- **Phase 3** (plan lines ~271–292): build `FaultGeometry geom(domain, seed, &mpi, false)` → call the **stress source** `geom.ComputeParams*(...)` / `geom.ComputeParamsFaultLocal(...)` → `resolver.ResolveRateState(...)` → `geom.SetRateStatePerDOF(rs)`.
- **Phase 5** (plan lines ~351–362): **construct** `RateStateFaultOperator<ParMesh,2> fault_op(&geom, &friction, &aging, seed, &mpi)` → `fault_op.SetSAFSMode(true, &geom.GetTauPre(), &geom.sigma_n_per_dof())`.

So `ComputeParams*` runs **before** the operator constructor. But the BP5 constructor enforces the **opposite** order with a hard verify:

```cpp
// fault/rate_state_fault.hpp:149
MFEM_VERIFY(geom_ == nullptr || !geom_->HasParams(),
            "RateStateFaultOperator: BP5 ctor must run BEFORE "
            "FaultGeometry::ComputeParams; otherwise BP5-mode "
            "tau_pre_ caches sidecar values silently.");
```

and `FaultGeometry::ComputeParamsFaultLocal` (and the sidecar `ComputeParams`) set `params_computed_ = true` (`fault/fault_geometry.hpp:695`), so `geom_->HasParams()` is `true` at the point Phase 5 constructs the operator. **The driver aborts at the constructor.**

There is a second, deeper conflict that prevents the trivial "just swap the two phases" fix:

1. The constructor **caches** per-DOF friction arrays at construction time:
   ```cpp
   // fault/rate_state_fault.hpp:160-171
   const Vector &V_init = geom_->GetVInit();          // asserts size == 2*num_nodes_
   Dc_values_     = geom_->GetDcValues();
   tau_pre_       = geom_->GetTauPre();
   V_init_values_ = geom_->GetVInit();
   ```
   `Dc_values_` and `V_init_values_` are **never re-sourced** by `SetSAFSMode` (only `tau_pre`/`sigma_n` are re-pointed). Therefore `SetRateStatePerDOF(rs)` (which fills `dc_values_`/`V_init_vec_`) **must run before** the constructor, or the operator caches NaN `Dc`/`V_init` (the `compute_bp5_params=false` ctor NaN-fills these — `fault_geometry.hpp:276-289`) and the friction RHS silently uses NaN.

2. `resolver.ResolveRateState(..., geom.sigma_n_per_dof())` consumes σ_n that is produced **only** by `ComputeParams*`. But (1) forces `ComputeParams*` to run **after** the constructor, while `ResolveRateState`→`SetRateStatePerDOF` must run **before** it. So the plan's data flow (σ_n from the stress source → resolver) is **circular** with the construction-order constraint.

**Trigger:** Any run that follows the plan literally — serial or parallel, BP5 or SAF. First failure is the `MFEM_VERIFY` at `rate_state_fault.hpp:149` during `seas_op` setup.

**Actual behavior (as planned):** `ComputeParams*` (Phase 3) sets `params_computed_=true`; Phase 5 constructor verify fails → `MFEM_VERIFY` abort. If a naive reviewer "fixes" it by moving the ctor earlier, `Dc_values_`/`V_init_values_` cache NaN → NaN friction RHS → NaN tripwire / equilibrium-init failure.

**Expected behavior:** The canonical order the operator documents (`rate_state_fault.hpp:140-148`: `ctor → ComputeParams → SetSAFSMode`), with `SetRateStatePerDOF` inserted **before** the ctor and σ_n for the resolver sourced **independently** of `geom.sigma_n_per_dof()`:

```
1. geom(domain, seed, &mpi, /*compute_bp5_params=*/false)
2. compute sigma_n_total independently (depth profile / sidecar / FaultLocalPrestress scalar)  →  Vector sigma_n_total
3. rs = resolver.ResolveRateState(..., sigma_n_total)          // NOT geom.sigma_n_per_dof()
4. geom.SetRateStatePerDOF(rs)                                 // fills a/Dc/eta/V_init; does NOT set params_computed_
5. RateStateFaultOperator fault_op(&geom, ...)                 // caches valid Dc/V_init; HasParams()==false → verify passes
6. geom.ComputeParams*/ComputeParamsFaultLocal(...)            // sets tau_pre/sigma_n, params_computed_=true
7. fault_op.SetSAFSMode(true, &geom.GetTauPre(), &geom.sigma_n_per_dof())
```

**Suggested fix (edit the plan's Phase 3/5 prose + the §"Critical-path call graph" ordering, and add a `SetRateStatePerDOF must NOT set params_computed_` note):**
```diff
 ### Phase 3 — driver wiring
-FaultGeometry<ParMesh> geom(domain, seed, &mpi, /*compute_bp5_params=*/false);
-// pre-stress + σ_n per DOF from the stress source (reuse spatial_dyn_driver.cpp:1337-1458 dispatch)
-//   FaultLocalPrestress → geom.ComputeParamsFaultLocal(tau_strike, tau_dip, sigma_n, P_p)
-...
-spatial::RateStatePerDOFParams rs = resolver.ResolveRateState(
-    *cfg.rate_state, dof_coords_3d, dof_to_elem, dof_to_attr,
-    material, pmesh, spatial::PorePressureSpec{}, geom.sigma_n_per_dof());
-geom.SetRateStatePerDOF(rs);
+FaultGeometry<ParMesh> geom(domain, seed, &mpi, /*compute_bp5_params=*/false);
+// Source TOTAL normal stress independently of geom (the stress source has
+// not run yet — it must run AFTER the fault-operator ctor; see Phase 5).
+Vector sigma_n_total;  /* from [stress] depth profile / FaultLocalPrestress scalar */
+spatial::RateStatePerDOFParams rs = resolver.ResolveRateState(
+    *cfg.rate_state, dof_coords_owned, dof_to_elem, dof_to_attr,
+    material, pmesh, pp_spec, sigma_n_total);   // see R-002 (owned coords) and R-003 (total σ_n)
+geom.SetRateStatePerDOF(rs);   // MUST NOT set params_computed_ (HasParams stays false)
+// Stress source (tau_pre/sigma_n) is deferred to Phase 5, AFTER the ctor.
```
```diff
 ### Phase 5 — construct operator, then stress source, then SAFS mode
 RateStateFaultOperator<ParMesh,2> fault_op(&geom, &friction, &aging, seed, &mpi);
+// Stress source runs HERE, after the ctor (rate_state_fault.hpp:149 requires
+// !geom.HasParams() at construction):
+//   FaultLocalPrestress → geom.ComputeParamsFaultLocal(tau_strike, tau_dip, sigma_n, P_p)
+//   sidecar/tensor     → geom.ComputeParams(...)
 fault_op.SetSAFSMode(true, &geom.GetTauPre(), &geom.sigma_n_per_dof());
```

**Test case (add to Phase-3/5 integration test; demonstrates the abort with the plan's order):**
```cpp
TEST(SpatialSeasWiring, R001_CtorMustPrecedeComputeParams) {
  // Reproduce the PLAN's order: ComputeParams BEFORE the operator ctor.
  FaultGeometry<ParMesh> geom(domain, seed, &mpi, /*compute_bp5_params=*/false);
  geom.ComputeParamsFaultLocal(/*tau_strike=*/ -2.0e6, /*tau_dip=*/ 0.0,
                               /*sigma_n=*/ 25.0e6, /*P_p=*/ 0.0);
  ASSERT_TRUE(geom.HasParams());                       // stress source ran
  // The plan now constructs the operator -> hard verify fails:
  EXPECT_DEATH(
      { RateStateFaultOperator<ParMesh,2> op(&geom, &fric, &aging, seed, &mpi); },
      "BP5 ctor must run BEFORE");
}
TEST(SpatialSeasWiring, R001_CorrectOrderSucceeds) {
  FaultGeometry<ParMesh> geom(domain, seed, &mpi, false);
  geom.SetRateStatePerDOF(rs_with_valid_Dc_Vinit);     // before ctor: Dc/V_init valid
  RateStateFaultOperator<ParMesh,2> op(&geom, &fric, &aging, seed, &mpi);  // HasParams()==false -> OK
  geom.ComputeParamsFaultLocal(-2.0e6, 0.0, 25.0e6, 0.0);                  // after ctor
  op.SetSAFSMode(true, &geom.GetTauPre(), &geom.sigma_n_per_dof());
  EXPECT_TRUE(op.IsSAFSMode());
  EXPECT_TRUE(std::isfinite(op.GetDcValuesCachedMin()));  // not NaN
}
```

---

### [R-002] CRITICAL [drivers/spatial_seas_driver.cpp / Phase 3 — fault-DOF coordinate sourcing] — Plan reuses the **local-ordered** `GetFaultDOFCoords3D` but pairs it with **owned-ordered** getters and an owned-sized assertion → parallel size/order mismatch

**Category:** BUG (the exact "owned-vs-local seam" the plan's own Risk table claims to have eliminated)

**Description:**
Phase 3 step 1 says the new getters "emit per-fault-DOF data in the canonical **owned** order (the order `RestrictToOwnedFault` produces)" and parenthetically claims "**(Coords/basis already exist via `GetFaultDOFCoords3D`/`GetFaultDOFBasis`.)**" (plan line ~264), then the driver snippet does:

```cpp
Vector dof_coords_3d; domain.GetFaultDOFCoords3D(dof_coords_3d);   // plan line ~285
spatial::RateStatePerDOFParams rs = resolver.ResolveRateState(
    *cfg.rate_state, dof_coords_3d, dof_to_elem, dof_to_attr, ...); // dof_to_elem/attr claimed "owned"
geom.SetRateStatePerDOF(rs);                                        // asserts rs.a.Size()==NumFaultDOFs() (owned)
```

But `GetFaultDOFCoords3D` is **local-ordered**, not owned-ordered:

```cpp
// domain/elasticity_operator_traction.inl:213-251
dof_coords_3d.SetSize(3 * num_fault_dofs_);          // LOCAL count
... iterate fault_interior_faces_ THEN ALL fault_shared_faces_ ...   // includes non-owned shared DOFs
```

and the two counts genuinely differ in parallel:
```cpp
// domain/elasticity_operator_setup.inl
num_fault_dofs_       = num_fault_faces_       * nbf_per_face_;   // :661  (LOCAL: interior + ALL shared)
num_owned_fault_dofs_ = num_owned_fault_faces_ * nbf_per_face_;  // :1055 (OWNED: only owned shared)
```
`FaultGeometry` and `RateStateFaultOperator` are built on the **owned** count (`fault_geometry.hpp:122` uses `GetNumOwnedFaultDOFs()` + `RestrictToOwnedFault`; `rate_state_fault.hpp:124` `num_nodes_ = geom->NumFaultDOFs()` = owned). And `ResolveRateState` derives `N = dof_coords_3d.Size()/3` then asserts:
```cpp
// spatial/code/spatial_friction.cpp:1879-1885
const int N = dof_coords_3d.Size() / 3;
MFEM_VERIFY(dof_to_elem.Size() == N, "...dof_to_elem.Size() != N");
MFEM_VERIFY(dof_to_attr.Size() == N, "...dof_to_attr.Size() != N");
```

So on `np>1`: `dof_coords_3d` has `N_local` entries, `dof_to_elem/attr` (owned getters) have `N_owned < N_local` → **`ResolveRateState` aborts** on the `dof_to_elem.Size() != N` verify. Even if the implementer instead makes the new getters local to match, `rs` then has `N_local` entries and **`SetRateStatePerDOF` aborts** on its `rs.a.Size()==NumFaultDOFs()` (owned) assert. The plan cannot be implemented as written for parallel runs (the Phase-3 acceptance criterion requires `mpirun -np 4`). In **serial** the two counts coincide, masking the bug — which is exactly how it would slip through a serial-only smoke test.

**Trigger:** `mpirun -np N` (N>1) with any fault partitioned across ranks such that at least one rank holds a shared fault face owned by a neighbor (the generic case).

**Actual behavior:** `MFEM_VERIFY` abort in `ResolveRateState` (or in `SetRateStatePerDOF`), or — if both verifies were removed — silently misaligned per-DOF `a/Dc/eta/V_init` (friction params from the wrong fault DOF).

**Expected behavior:** Coordinates fed to `ResolveRateState` must be in the **owned** order and length, consistent with the owned `dof_to_elem/attr` and `SetRateStatePerDOF`. Either restrict the existing local coords with `RestrictToOwnedFault(local_coords, owned_coords, 3)`, or have the new getter emit owned coords directly. The plan's "coords already exist via `GetFaultDOFCoords3D`" shortcut is false for the owned layout.

**Suggested fix (plan Phase 3 step 1 + the driver snippet):**
```diff
-   These walk the **same** `fault_interior_faces_` then `fault_shared_faces_` with `nbf_per_face_` QPs that
-   `SetupFaultInfo` already uses ... (Coords/basis already exist via
-   `GetFaultDOFCoords3D`/`GetFaultDOFBasis`.)
+   These walk the same `fault_interior_faces_` then `fault_shared_faces_` with `nbf_per_face_` QPs, then
+   **restrict to the owned set with `RestrictToOwnedFault(...)`** so they line up with `GetNumOwnedFaultDOFs()`.
+   NB: the existing `GetFaultDOFCoords3D`/`GetFaultDOFBasis` emit the LOCAL set (size `num_fault_dofs_`,
+   interior + ALL shared faces) — they are NOT owned-ordered and must be restricted before use here.
```
```diff
-Vector dof_coords_3d; domain.GetFaultDOFCoords3D(dof_coords_3d);
+Vector local_coords_3d; domain.GetFaultDOFCoords3D(local_coords_3d);   // LOCAL
+Vector dof_coords_3d;  domain.RestrictToOwnedFault(local_coords_3d, dof_coords_3d, /*comps=*/3);  // OWNED
+// (or add domain.GetOwnedFaultDOFCoords3D(...) that restricts internally)
```

**Test case:**
```cpp
TEST(SpatialSeasWiring, R002_OwnedVsLocalCoordsParallel) {  // run on np=4
  Vector local_coords; domain.GetFaultDOFCoords3D(local_coords);
  EXPECT_EQ(local_coords.Size(), 3 * domain.GetNumFaultDOFs());          // LOCAL
  // The PLAN feeds local coords with owned dof_to_elem -> mismatch:
  Array<int> dof_to_elem; domain.GetFaultDOFToElem(dof_to_elem);         // OWNED
  if (domain.GetNumOwnedFaultDOFs() != domain.GetNumFaultDOFs()) {
    EXPECT_NE(local_coords.Size() / 3, dof_to_elem.Size());              // sizes disagree
    EXPECT_DEATH(resolver.ResolveRateState(blk, local_coords, dof_to_elem,
                 dof_to_attr, material, pmesh, pp, sn), "dof_to_elem.Size\\(\\) != N");
  }
  // Correct path: restrict first
  Vector owned_coords; domain.RestrictToOwnedFault(local_coords, owned_coords, 3);
  EXPECT_EQ(owned_coords.Size() / 3, dof_to_elem.Size());                // now consistent
}
```

---

### [R-003] MODERATE [drivers/spatial_seas_driver.cpp / Phase 3 — `ResolveRateState` σ_n argument] — Effective σ_n passed into a parameter that expects **total** σ_n → pore pressure double-subtracted

**Category:** BUG (sign/semantics)

**Description:**
The resolver's last argument is **total** normal stress and it subtracts pore pressure internally:
```cpp
// spatial/code/spatial_friction.cpp:2007-2017
real_t sigma_n_total = sigma_n_total_per_dof(i);     // arg is TOTAL
const real_t P_p   = pp.P_p_pa + pp.P_p_grad_pa_per_m * depth_i;
real_t sigma_n_eff = sigma_n_total - P_p;            // subtracts pore pressure HERE
```
The plan (line ~290) passes `geom.sigma_n_per_dof()`. After `ComputeParamsFaultLocal(tau_strike, tau_dip, sigma_n, P_p)`, that member holds **effective** σ_n:
```cpp
// fault/fault_geometry.hpp:673-690
const real_t sigma_n_eff = sigma_n - P_p;            // already EFFECTIVE
sigma_n_per_dof_(i) = sigma_n_eff;
```
So the plan feeds **effective** σ_n where the resolver expects **total**. It "works" only when `PorePressureSpec{}` is empty (P_p = 0), which is what the BP5 snippet passes — but Phase 8 explicitly plans to "reuse `[stress.pore_pressure]`", i.e. a **non-empty** `PorePressureSpec`. In that case pore pressure is subtracted twice (once in `ComputeParamsFaultLocal`, once in `ResolveRateState`), yielding a too-low / possibly tensile `sigma_n_eff` and either an `MFEM_VERIFY(sigma_n_eff > 0)` abort (`spatial_friction.cpp:2064`) or a silently wrong friction regime. (This finding is partly subsumed by R-001's recommendation to source `sigma_n_total` independently — make that source **total**, not effective.)

**Trigger:** A SAF config (Phase 8) with a non-empty `[stress.pore_pressure]` and a `ComputeParams*` path that already applies pore pressure.

**Actual behavior:** σ_n_eff = σ_n − P_p − P_p → friction params computed at the wrong effective normal stress (or abort on the >0 verify).

**Expected behavior:** Pass **total** σ_n to `ResolveRateState` (and let it apply pore pressure once), or pass effective σ_n with an **empty** `PorePressureSpec` consistently and document the contract.

**Suggested fix (plan Phase 3 prose + driver snippet):**
```diff
-spatial::RateStatePerDOFParams rs = resolver.ResolveRateState(
-    ..., spatial::PorePressureSpec{}, geom.sigma_n_per_dof());
+// ResolveRateState's last arg is TOTAL normal stress (it subtracts pore
+// pressure internally via `pp`).  Pass a TOTAL-σ_n vector + the SAME pp spec;
+// do NOT pass geom.sigma_n_per_dof() (that is already EFFECTIVE after
+// ComputeParamsFaultLocal -> pore pressure would be subtracted twice).
+spatial::RateStatePerDOFParams rs = resolver.ResolveRateState(
+    ..., pp_spec, sigma_n_total_per_dof);
```

**Test case:**
```cpp
TEST(SpatialSeasWiring, R003_PorePressureNotDoubleSubtracted) {
  const double sn_total = 50e6, Pp = 10e6;
  PorePressureSpec pp; pp.P_p_pa = Pp;
  Vector sn_total_vec(N); sn_total_vec = sn_total;
  auto rs = resolver.ResolveRateState(blk, owned_coords, dte, dta, mat, mesh, pp, sn_total_vec);
  for (int i = 0; i < N; ++i)
    EXPECT_NEAR(rs.sigma_n_eff(i), sn_total - Pp, 1.0);     // 40 MPa, subtracted ONCE
  // The PLAN's path (effective fed as total + same pp) would give 30 MPa:
  Vector sn_eff_vec(N); sn_eff_vec = sn_total - Pp;          // what geom.sigma_n_per_dof() holds
  auto rs_bad = resolver.ResolveRateState(blk, owned_coords, dte, dta, mat, mesh, pp, sn_eff_vec);
  EXPECT_NEAR(rs_bad.sigma_n_eff(0), sn_total - 2*Pp, 1.0);  // 30 MPa — double-subtracted
}
```

---

### [R-004] MODERATE [drivers/spatial_seas_driver.cpp / Phase 5 — RK45 driver snippet] — Non-existent `SetTolerances` and wrong `Step(...)` signature in the normative time-loop

**Category:** BUG (API mismatch in a snippet the plan declares normative)

**Description:**
Phase 5 (plan lines ~364–367):
```cpp
DormandPrinceRK45 ode; ode.Init(seas_op);
ode.SetTolerances(cfg.time.rk45_atol, cfg.time.rk45_rtol);   // (1)
// loop: ode.Step(state, t, dt); clamp dt ≤ dt_max_years·yr; …   // (2)
```
The actual `DormandPrinceRK45` API (`solver/time_stepper.hpp`):
- (1) There is **no `SetTolerances`**. Tolerances are set with `SetAbsTol(real_t)` and `SetRelTol(real_t)` (`time_stepper.hpp:163-164`).
- (2) `Step` takes the **operator** as the first argument: `bool Step(TimeDependentOperator &op, Vector &state, real_t &t, real_t &dt)` (`time_stepper.hpp:227`). The plan's `Step(state, t, dt)` omits `op` and will not compile.
- The hand-rolled `clamp dt ≤ dt_max` is unnecessary — the stepper has `SetDtMax(real_t)` (`time_stepper.hpp:169`) and enforces the cap internally; the plan should set it once rather than clamp in the loop.

This contradicts the §Constraints item "do not alter existing signatures" — here the plan *misstates* an existing signature it intends to call.

**Trigger:** Compiling the Phase-5 driver as written.

**Actual behavior:** Compile errors: no member `SetTolerances`; no `Step(Vector&, real_t&, real_t&)` overload.

**Expected behavior:** Use the real API.

**Suggested fix (plan Phase 5):**
```diff
-DormandPrinceRK45 ode; ode.Init(seas_op);
-ode.SetTolerances(cfg.time.rk45_atol, cfg.time.rk45_rtol);
-double dt = (cfg.time.dt_init > 0) ? cfg.time.dt_init : 0.01*L_nuc/V_nuc;
-// loop: ode.Step(state, t, dt); clamp dt ≤ dt_max_years·yr; …
+DormandPrinceRK45 ode; ode.Init(seas_op);
+ode.SetAbsTol(cfg.time.rk45_atol);
+ode.SetRelTol(cfg.time.rk45_rtol);
+ode.SetDtMax(cfg.time.dt_max_years * seconds_per_year);   // internal cap; no manual clamp
+ode.SetMPIContext(&mpi);                                   // MPI_Allreduce(MAX) error reduction (CLAUDE.md)
+double dt = (cfg.time.dt_init > 0) ? cfg.time.dt_init : 0.01*L_nuc/V_nuc;
+// loop: ode.Step(seas_op, state, t, dt);   // operator is the FIRST argument
```

**Test case:**
```cpp
TEST(SpatialSeasTimeLoop, R004_RK45ApiCompilesAndSteps) {
  DormandPrinceRK45 ode; ode.Init(seas_op);
  ode.SetAbsTol(1e-7); ode.SetRelTol(1e-50); ode.SetDtMax(0.1*kSecPerYear);
  double t = 0, dt = 1e3;
  bool ok = ode.Step(seas_op, state, t, dt);   // 4-arg form
  EXPECT_TRUE(ok);
  EXPECT_LE(dt, 0.1*kSecPerYear + 1e-9);        // dt_max honored internally
}
```

---

### [R-005] MODERATE [Phase 5 — `DieterichRuinaFriction::Constants` from `rs`] — Undefined `rs_uniform` with wrong field names; per-DOF→scalar extraction is specified nowhere

**Category:** BUG / ASSUMPTION (normative snippet)

**Description:**
Phase 5 (plan lines ~352–356):
```cpp
DieterichRuinaFriction::Constants fc{ /*V0=*/rs_uniform.V0, /*f0=*/rs_uniform.f0,
                                      /*b=*/rs_uniform.b, /*Dc=*/seed.L0 };
AgingLawPsi aging(fc.b, fc.V0, fc.f0);
```
Three problems:
1. `rs_uniform` is never produced. The resolver returns `rs` (a `RateStatePerDOFParams`), whose `b`, `V_0`, `f_0` are **`Vector`s** (per-DOF), not scalars (`spatial_friction.hpp:687-692`). The Phase-3 edge case mandates asserting these are uniform but **never shows extracting the scalar** (e.g. `rs.b(0)`), leaving a gap the implementer must invent.
2. Field-name mismatch: the struct members are **`V_0`** and **`f_0`** (with underscores), not `V0`/`f0`. `rs.V0` / `rs.f0` do not compile.
3. `Dc = seed.L0`: the constructor's scalar `Dc` is overridden per-DOF by `geom`'s `Dc_values_` in the hot path (CLAUDE.md "Per-DOF Dc in BP5"), so seeding it from `seed.L0` is at best vestigial and at worst misleading; if any code path reads `friction.Constants().Dc` it would get the wrong value. (`Constants{V0,f0,b,Dc}` field order itself is correct — verified `dieterich_ruina.hpp:45-49`.)

**Trigger:** Compiling Phase 5 verbatim.

**Actual behavior:** Compile errors (`rs_uniform` undeclared; no member `V0`/`f0`).

**Expected behavior:** Extract uniform scalars from the asserted-uniform per-DOF vectors using the correct field names.

**Suggested fix (plan Phase 5):**
```diff
-DieterichRuinaFriction::Constants fc{ /*V0=*/rs_uniform.V0, /*f0=*/rs_uniform.f0,
-                                      /*b=*/rs_uniform.b, /*Dc=*/seed.L0 };
+// Phase-3 asserted rs.b / rs.f_0 / rs.V_0 are uniform; take element 0 as the scalar.
+MFEM_VERIFY(rs.b.Size() > 0, "empty rate-state");
+DieterichRuinaFriction::Constants fc{ /*V0=*/rs.V_0(0), /*f0=*/rs.f_0(0),
+                                      /*b=*/rs.b(0), /*Dc=*/rs.Dc(0) };
 AgingLawPsi aging(fc.b, fc.V0, fc.f0);
```
(Per-DOF `Dc` still flows through `geom`'s `Dc_values_`; `fc.Dc` is only the scalar fallback.)

**Test case:**
```cpp
TEST(SpatialSeasWiring, R005_UniformScalarExtraction) {
  auto rs = resolver.ResolveRateState(/*uniform b/f0/V0 block*/...);
  for (int i = 1; i < rs.b.Size(); ++i) {           // Phase-3 uniformity precondition
    EXPECT_NEAR(rs.b(i),   rs.b(0),   1e-12*std::abs(rs.b(0)));
    EXPECT_NEAR(rs.V_0(i), rs.V_0(0), 1e-12*std::abs(rs.V_0(0)));
    EXPECT_NEAR(rs.f_0(i), rs.f_0(0), 1e-12*std::abs(rs.f_0(0)));
  }
  DieterichRuinaFriction::Constants fc{ rs.V_0(0), rs.f_0(0), rs.b(0), rs.Dc(0) };
  EXPECT_GT(fc.V0, 0.0);   // compiles + sane
}
```

---

### [R-006] LOW [POSSIBLE] [Phase 4 — `SetElasticityOptions` on a DG space] — Plan applies the elasticity near-null-space on CG_AMG while simultaneously asserting "DG sparsity breaks the CFEM assumption" for GMRES_AMG

**Category:** ASSUMPTION (internal inconsistency; numerical-efficiency risk, not a crash)

**Description:**
Phase 4 sets `HypreBoomerAMG::SetElasticityOptions(fes_)` on the `CG_AMG` path but, for `GMRES_AMG`, says "no elasticity near-null-space ... DG sparsity breaks the CFEM assumption — keep `SetSystemsOptions` only." The DG-vs-CFEM concern applies equally to the CG path: `SetElasticityOptions` builds rigid-body modes by projecting coefficients onto the (DG) space (`linalg/hypre.cpp:5242+`) and BoomerAMG then coarsens the block-dense DG sparsity — this won't crash (verified: `SetSystemsOptions`/`SetElasticityOptions` handle `Ordering::byNODES` correctly, `hypre.cpp:5192-5212`), but the near-null-space/strength heuristics are tuned for CFEM elasticity and may give poor or non-`h`-independent convergence on SIPG/BR2 DG. MFEM's own DG-elasticity example (`ex17p`) does **not** use `SetElasticityOptions`. The plan already lists "AMG slow/non-convergent for DG elasticity" as a risk, so this is a consistency/expectation note rather than a new defect: the Phase-4 acceptance criterion "AMG iteration count bounded (h-independent within a factor)" may not hold, and the `residual_check` + MUMPS fallback are the right guards to keep mandatory rather than optional.

**Trigger:** SAF-scale DG mesh; observe AMG iteration count growth with refinement.

**Actual behavior:** Possibly slow/growing iteration counts; correctness preserved if `ksp_maxit` is large enough and `residual_check` catches non-convergence.

**Expected behavior:** Treat `SetElasticityOptions`-on-DG as unproven; make `residual_check` non-optional on the AMG paths and keep the MUMPS small-mesh cross-check as the correctness anchor.

**Suggested fix (plan Phase 4 prose):**
```diff
-  - Optional residual check (`cfg.solver.residual_check`): compute `||Kx-b||/||b||` and warn if `> 10*ksp_rtol`.
+  - **Mandatory** residual check on AMG paths: compute `||Kx-b||/||b||`; on `> 10*ksp_rtol` warn loudly
+    with iteration count (a silently non-converged solve corrupts traction→friction→blowup).  Note that
+    `SetElasticityOptions` rigid-body modes are a CFEM heuristic applied to a DG space (unproven; cf. the
+    GMRES note); if iteration counts grow with h, drop to `SetSystemsOptions`-only or MUMPS.
```

**Test case:** (covered by the Phase-4 `test_spatial_seas_iterative_vs_direct` equivalence test plus the multi-resolution iteration-count smoke check already specified; add a hard assert that the residual check fires on a deliberately under-iterated solve.)

---

## Summary
- Critical issues: 2  (R-001 ctor-ordering abort + circular σ_n dependency; R-002 owned-vs-local coords seam)
- Moderate issues: 3  (R-003 σ_n effective-as-total; R-004 RK45 API; R-005 `rs` field names / scalar extraction)
- Low issues: 1  (R-006 SetElasticityOptions-on-DG, POSSIBLE)
- Plan compliance (with the in-tree code it targets): **PARTIAL** — most signatures the plan cites are accurate (SolverType enum, `ElasticityDomainOperator` ctor `:106`, `ResolveRateState` arg list, `Constants{V0,f0,b,Dc}`, `AgingLawPsi(b,V0,f0)`, owned/local getters all verified present), but the **construction-order contract** (R-001) and the **owned-vs-local DOF layout** (R-002) — the two hardest parts and the plan's own top-2 risks — are specified incorrectly and will fail at runtime.
- Verdict: **FAIL — must fix before implementing.** R-001 and R-002 will abort every run (R-001 serial + parallel; R-002 parallel) and are not surface typos — they require reordering Phase 3/5 and adding an owned-restriction step. R-003/R-004/R-005 are concrete snippet fixes. Address R-001 and R-002 in the plan text before handing it to the implementing agent, since the plan is the normative contract.

## Unreviewed Areas
- **Physics correctness of BP5 parity numbers** (recurrence ~240 yr, dip slip ≈ 0): cannot be checked at plan stage — gated by Phase 7's golden comparison. Not reviewed.
- **SAF non-planar far-field loading function** (Phase 8 open question #1): flagged by the plan itself as an unresolved physics design item requiring user/SCEC input; out of scope for a code/plan audit.
- **`config/bp5_params.hpp` `seed` defaults** actually used by `RateStateFaultOperator` (`sigma_n_bp5_`, `Vp_bp5_`, `bp5_params_` cached at `rate_state_fault.hpp:129-131`): in SAFS mode σ_n/τ_pre route through `SetSAFSMode` pointers, but I did not exhaustively trace every read site to confirm none falls back to the BP5 scalar when `safs_mode_==true`. Recommend the implementer add an assert that no BP5-scalar path is reachable once `IsSAFSMode()`.
- **`Makefile` target wiring** (Phase 0): object-prerequisite list not audited against the real `seas_driver` rule; verify no wave OBJs leak in and no QD OBJ is missing at build time.
- **`tpv104_checkpoint.hpp` reuse for the QD state** (Phase 6): the schema fit (`[s_dip,s_strike,psi]*N_owned` vs the checkpoint's `DOFData`-oriented layout) was not verified; flagged as a Phase-6 implementation risk.
