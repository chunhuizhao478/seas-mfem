# Code Review (Round 3): Phase H + Phase 5 implementation (2026-05-18)

## Review Scope

- **Plan:** `safs/project_7.0_alternative/document/spatial_dynamic_rupture_plan.md` (rev-3)
- **Round-1 review:** `debug_document/spatial_dynamic_rupture_review_2026-05-18.md`
- **Round-1 fix report:** `debug_document/spatial_dynamic_rupture_fix_2026-05-18.md`
- **Round-2 review:** `debug_document/spatial_dynamic_rupture_review_round2_2026-05-18.md`
- **Round-2 fix report:** `debug_document/spatial_dynamic_rupture_fix_round2_2026-05-18.md`
- **Files reviewed (NEW since round-2):**
  - `miniapps/seas/dynamic/godunov_flux_pool.{hpp,cpp}` (NEW — Phase H.1)
  - `miniapps/seas/dynamic/wave_operator.hpp` (modified — Phase H additive)
  - `miniapps/seas/dynamic/wave_operator.inl` (modified — Phase H additive)
  - `miniapps/seas/dynamic/fault_face_flux.{hpp,cpp}` (modified — Phase H)
  - `miniapps/seas/dynamic/spatial_setup.hpp` (NEW — Phase 5c)
  - `miniapps/seas/fault/fault_geometry.hpp` (modified — Phase 5a)
  - `miniapps/seas/io/tpv104_checkpoint.hpp` (modified — Phase 5b)
  - `miniapps/seas/tests/unit/test_phaseh_godunov_flux_pool.cpp` (NEW)
  - `miniapps/seas/tests/unit/test_phaseh_lsw_forced_rupture.cpp` (NEW)
  - `miniapps/seas/tests/unit/test_spatial_setup.cpp` (NEW)
- **Domain context consulted:**
  - `miniapps/seas/CLAUDE.md`
  - `mfem/fem/intrules.hpp` (`IntegrationPoint::Init` semantics)
  - `mfem/fem/geom.hpp` (`Geometries.GetCenter`)
  - `miniapps/seas/dynamic/godunov_flux.hpp` (existing `GodunovFlux`)
  - `miniapps/seas/dynamic/tpv205_friction.hpp` (`LSWFrictionCoefficient_TPV205`)
  - `miniapps/seas/dynamic/tpv205_setup.hpp` (`InitializeFaultDOFs_TPV205`)
  - `miniapps/seas/spatial/code/spatial_friction.{hpp,cpp}` (round-1/2 baseline)
- **Test sweep (round-3):**
  - `seas_test_phaseh_godunov_flux_pool`   — 20 / 20
  - `seas_test_phaseh_lsw_forced_rupture`  — 32 / 32
  - `seas_test_spatial_setup`              — 63 / 63
  - `seas_test_spatial_friction_config`    — 28 / 28 (round-2 carry-over)
  - `seas_test_spatial_friction_resolver`  — 110 / 110 (round-2)
  - `seas_test_spatial_velocity_bundle`    — 6 / 6
  - `seas_test_spatial_stress_bundle`      — 6 / 6
  - `seas_test_spatial_constant_stress_source` — 45 / 45
  - `seas_test_compute_safs_params`        — 13 / 13 (byte-exact regression intact)
  - `seas_test_tpv104_checkpoint`          — 164 / 164 (DRIVER_TAG_V1 extension is back-compat)

## Findings

---

### [R-301] [CRITICAL] [wave_operator.hpp:SetGodunovFluxPool / wave_operator.inl hot loop] — `flux_pool_` is stored but NEVER read; heterogeneous-material runs silently use the scalar `flux_`

**Category:** DEVIATION / BUG (Phase H.2 / H.5 incomplete)

**Description:**
Phase H of the plan requires the WaveOperator hot loop to dispatch flux per element via the `GodunovFluxPool` so heterogeneous-material runs use the correct per-element Jacobian (plan §Phase H.2, §Phase H.5).  The implementer added `SetGodunovFluxPool(const GodunovFluxPool*)`, a non-owning `flux_pool_` member, and even built `GodunovFluxPool::Build` with full dedup logic — but the hot loop in `dynamic/wave_operator.inl` never references `flux_pool_`.  `grep flux_pool_` returns zero matches in `wave_operator.inl`.  The Phase H.3 heterogeneous-CFL extension is also missing: `ComputeMaxDt` (`wave_operator.inl:4964`) ends in `return cfl_mixed_flux_factor * cfl * h_min_ / flux_.GetCp();` — scalar-material `flux_` only.

The implementer documented the omission in the setter's docstring (`wave_operator.hpp:130-144`):

> Default `nullptr` ⇒ TPV/BP5 byte-exact: every hot-loop site that reads `flux_` continues to use the existing scalar-material `flux_` member with zero extra conditional cost. **The actual per-element dispatch ... is a Phase H follow-up edit to `wave_operator.inl` ... Until that follow-up lands, this setter only stores the pointer; the dispatch still uses `flux_` everywhere.**

This is a CRITICAL plan deviation that breaks the central guarantee of Phase H.  A driver that constructs a `GodunovFluxPool` from CVMH material and calls `wave.SetGodunovFluxPool(&pool)` gets a build-time success but a runtime simulation that propagates the wave equation with whatever scalar `(λ, μ, ρ)` the original ctor was given.  The whole point of Phase H — running SAFS with real material — silently does the wrong physics.

**Trigger:**
Any heterogeneous-material run: build a `MaterialField::MakeCoefficient(...)`, populate `GodunovFluxPool`, set on `WaveOperator`, run `Mult`.  The CFL is wrong (per H.3) and every flux evaluation uses the scalar seed (per H.2/H.5).

**Actual behavior:**
Wave propagation uses scalar-material flux; heterogeneity has no effect.

**Expected behavior:**
Either (a) implement per-element flux dispatch + heterogeneous-CFL as the plan specifies, or (b) abort in `SetGodunovFluxPool` until that work lands, with a "Phase H dispatch is not yet implemented; calling SetGodunovFluxPool is meaningless" message — so callers fail loud instead of silently producing wrong results.

**Suggested fix:**
Option (b) for the immediate safety; option (a) for the real work.  Apply option (b) NOW:

```diff
   void SetGodunovFluxPool(const GodunovFluxPool *pool)
-  { flux_pool_ = pool; }
+  {
+     MFEM_ABORT("WaveOperator::SetGodunovFluxPool: Phase H.2/H.3/H.5 "
+                "per-element flux dispatch + heterogeneous-CFL + "
+                "bi-material shared-face MPI exchange are NOT yet "
+                "implemented in wave_operator.inl.  Setting this pool "
+                "is a no-op that would silently produce scalar-material "
+                "physics on heterogeneous input.  Land the follow-up "
+                "edits before re-enabling this setter.");
+     flux_pool_ = pool;
+  }
```

Then track the dispatch-implementation work as its own follow-up.

**Test case:**
```cpp
// tests/unit/test_phaseh_wave_operator_heterogeneous_dispatch.cpp
// (this file does NOT exist; plan §Phase H wants
//  test_phaseh_wave_operator_layered.cpp).
static void T_PHASEH_layered_propagation_differs_from_scalar()
{
   // Build a 1-D layered Vs material (top 1000 m/s, bottom 3500 m/s).
   // Build TWO WaveOperators: (i) scalar Vs=1000, (ii) "scalar with
   // GodunovFluxPool of layered Vs".  Run Mult on a Gaussian initial
   // condition for 1 step.  Assert the dQdt outputs DIFFER by ≥ 1e-3
   // relative L_inf.  If they are byte-identical, the pool is being
   // ignored.
   ...
   const real_t rel_diff = max_relative_diff(dQdt_scalar, dQdt_with_pool);
   TEST_ASSERT(rel_diff > 1e-3,
               "layered material must change dQdt vs scalar Vs=1000");
}
```

---

### [R-302] [CRITICAL] [wave_operator.inl:3617, 4539] — `FaultFrictionLaw::LSW_ForcedRupture` enum value has NO dispatch arm; setting it silently routes through rate-and-state Brent

**Category:** BUG (Phase H wiring missing)

**Description:**
The enum value `FaultFrictionLaw::LSW_ForcedRupture = 2` was added to `wave_operator.hpp:73-78`, and the new method `FaultFaceFlux::EvaluateADER_LSW_ForcedRupture` was implemented in `fault_face_flux.cpp:838-932`.  But the wave operator's dispatch in `wave_operator.inl` at TWO sites — interior fault path (line 3617) and shared-fault fallback (line 4539) — uses only `if (... == LSW) { EvaluateADER_LSW(...) } else { EvaluateADER(...) }`.  There is **no arm for `LSW_ForcedRupture`**.  `grep "LSW_ForcedRupture" wave_operator.inl` returns ZERO matches.

So when a driver calls `wave.SetFaultFrictionLaw(FaultFrictionLaw::LSW_ForcedRupture)`, the dispatch evaluates to `else` and routes through `EvaluateADER` — the rate-and-state Brent solver — not `EvaluateADER_LSW_ForcedRupture`.  The forced-rupture code path is **dead**.

This is the SAFS dynamic-rupture driver's intended primary dispatch (plan §Phase 4 Detailed main() step 15):

```cpp
if (cfg.law == FrictionLawKind::SlipWeakening) {
   wave.SetFaultFrictionLaw(FaultFrictionLaw::LSW_ForcedRupture);
}
```

With R-302 unfixed, a SAFS LSW run with TPV26/27 forced rupture would dispatch on rate-and-state — using `DOFData::a` / `psi` / `Dc` which `InitializeFaultDOFs_Spatial` defensively ZEROES — and either crash inside Brent (singular bracket) or silently produce wrong slip.

**Trigger:**
`wave.SetFaultFrictionLaw(FaultFrictionLaw::LSW_ForcedRupture)` followed by ADER step.

**Actual behavior:**
`EvaluateADER` (rate-and-state) runs on LSW-init DOFData; Brent fails or returns junk.

**Expected behavior:**
The dispatch routes through `EvaluateADER_LSW_ForcedRupture`.

**Suggested fix:**
Add the third arm at BOTH dispatch sites in `wave_operator.inl`.  At line 3617 (interior fault path):

```diff
                     if (fault_friction_law_ == FaultFrictionLaw::LSW)
                     {
                        fault_flux_->EvaluateADER_LSW(
                           fdata,
                           I_plus_local, I_minus_local,
                           dt,
                           I_imp_plus, I_imp_minus);
                     }
+                    else if (fault_friction_law_ == FaultFrictionLaw::LSW_ForcedRupture)
+                    {
+                       fault_flux_->EvaluateADER_LSW_ForcedRupture(
+                          fdata,
+                          I_plus_local, I_minus_local,
+                          dt,
+                          GetTime(),
+                          I_imp_plus, I_imp_minus);
+                    }
                     else
                     {
                        fault_flux_->EvaluateADER(fdata,
                                                  I_plus_local, I_minus_local,
                                                  dt,
                                                  I_imp_plus, I_imp_minus);
                     }
```

Apply the same diff at line 4539 (shared-fault fallback).  Note `GetTime()` (or whatever the existing ADER closure uses for time, e.g. `t_now`) needs to be the macro-step time the iterator passes in; verify against `Tpv104SubStepIterator::Step` to find the canonical name.

**Test case:**
```cpp
// tests/unit/test_phaseh_lsw_forced_rupture_dispatch.cpp (NEW)
static void T_PHASEH_dispatch_routes_to_forced_rupture()
{
   // Build a minimal WaveOperator + FaultFaceFlux, populate one fault
   // DOFData with T_forced=0 and t0_decay=0.5 (forced rupture active),
   // set FaultFrictionLaw::LSW_ForcedRupture, run one ADER face step at
   // t = 0.25, assert that the resulting mu_eff used in the solve is
   // mu_s + (mu_d - mu_s) * 0.5 (the midway value) NOT mu_s.
   ...
   // With the bug present, the dispatch routes through EvaluateADER
   // (rate-and-state) and mu_eff is undefined; the test catches it
   // either by an explicit "didn't go through forced-rupture path"
   // assertion or by a sentinel value the new path writes that the
   // RS path does not.
}
```

---

### [R-303] [CRITICAL] [spatial_setup.hpp:seed_static_dof_fields] — Re-introduces R-002 / R-003 corner-vs-centroid bug; impedances evaluated at element CORNER

**Category:** BUG (regression of round-1 R-002 / R-003 fix)

**Description:**
Round-1 R-002 and R-003 identified that `ip.Init(0)` lands at the reference-element ORIGIN (a CORNER), not the centroid.  The round-1 fix replaced `ip.Init(0)` with `Geometries.GetCenter(mesh.GetElementBaseGeometry(e))` in `spatial_friction.cpp:resolve_forced_impl` and `resolve_rs_impl`.

`dynamic/spatial_setup.hpp:seed_static_dof_fields` (the new Phase 5c per-DOF impedance + LSW init) was authored AFTER the round-1 fix and **reintroduces the same bug**:

```cpp
// dynamic/spatial_setup.hpp lines 64-68
mfem::ElementTransformation* T = mesh.GetElementTransformation(elem);
IntegrationPoint ip;
ip.Init(0);    // ← CORNER, not centroid (R-002 bug returned)
real_t lam, mu, rho;
material.EvalAt(elem, *T, ip, lam, mu, rho);
```

The cascade is severe — these impedances seed every per-DOF flux:

```cpp
d.Zp_plus  = d.Zp_minus = rho * cp;     // ρ_corner * c_p,corner
d.Zs_plus  = d.Zs_minus = rho * cs;
d.eta_p    = 0.5 * rho * cp;             // radiation damping
d.eta_s    = 0.5 * rho * cs;
```

For `MaterialField::MakeCoefficient` (the SAFS production path with CVMH/CVMS sidecars), each impedance is biased by the heterogeneity scale of one element across the corner→DOF distance.  The downstream `EvaluateADER_LSW_ForcedRupture` and `EvaluateADER` consume these impedances directly.

The accompanying comment in `seed_static_dof_fields` even claims the centroid is being used:

> // Reference centroid IntegrationPoint —
> // when Phase 5a's `fault_dof_ip(i)` is wired in by the driver, the
> // caller should pass that IP directly via the eight-arg overload
> // below.  **The simple centroid path is exact for constant-mode
> // MaterialField and a good per-element average for the heterogeneous
> // coefficient mode.**

The claim is false: `ip.Init(0)` is the corner.  The comment is also misleading because the referenced "eight-arg overload below" does not exist — `seed_static_dof_fields` has only one signature (7 args).

**Trigger:**
Driver calls `InitializeFaultDOFs_Spatial` with a heterogeneous `MaterialField` (`MakeCoefficient` with any non-constant Coefficient).

**Actual behavior:**
Per-DOF impedances are biased by element-scale heterogeneity.

**Expected behavior:**
Evaluate at the element centroid (R-002 fix shape), OR — better — at the cached `fault_dof_ip(i)` if supplied via the new `FaultGeometry` ctor.

**Suggested fix:**
Apply the same R-002 fix to `seed_static_dof_fields`:

```diff
-   mfem::ElementTransformation* T = mesh.GetElementTransformation(elem);
-   IntegrationPoint ip;
-   ip.Init(0);
+   mfem::ElementTransformation* T = mesh.GetElementTransformation(elem);
+   // True reference centroid for the bulk element.  ip.Init(0) lands
+   // at the reference ORIGIN (a CORNER), which biases impedances by
+   // one element of heterogeneity in MaterialField::MakeCoefficient
+   // mode.  Centroid is the best per-element stop-gap until the
+   // driver threads FaultGeometry::fault_dof_ip(i) through here.
+   const mfem::Geometry::Type gtype = mesh.GetElementBaseGeometry(elem);
+   const mfem::IntegrationPoint& ip = mfem::Geometries.GetCenter(gtype);
    real_t lam, mu, rho;
    material.EvalAt(elem, *T, ip, lam, mu, rho);
```

Also fix the misleading "eight-arg overload" comment.  Long term, add an overload that accepts `const std::vector<IntegrationPoint>& dof_ips` and use it when the driver constructs `FaultGeometry` with non-empty `dof_ips`.

**Test case:**
```cpp
// tests/unit/test_spatial_setup.cpp — add S-7
static void S_7_layered_material_impedances()
{
   // Mirror F-3 / R-9 from test_spatial_friction_resolver: unit-cube
   // hex with Vs(z) = 1000 + 4000*z.  At z=0.5 (centroid) Vs=3000;
   // at z=0 (corner) Vs=1000.
   ...
   RhoConst rho_c;  MuFromZ mu_c;  LambdaFromZ lam_c;
   auto material = MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c);
   ...
   InitializeFaultDOFs_Spatial<mfem::Mesh>(dof_data, 1, dof_to_elem,
                                            material, mh.mesh(),
                                            lsw, tau_pre, sigma_n_eff,
                                            T_forced, t0_decay);
   const real_t Vs_centroid = 3000.0;
   const real_t Vs_corner   = 1000.0;
   const real_t rho = 2670.0;
   const real_t Zs_expected = rho * Vs_centroid;
   const real_t Zs_corner   = rho * Vs_corner;
   TEST_NEAR(dof_data[0].Zs_plus, Zs_expected, 1.0,
             "Zs uses Vs(centroid), not Vs(corner)");
   TEST_ASSERT(std::abs(dof_data[0].Zs_plus - Zs_corner)
               > 0.5 * std::abs(Zs_corner - Zs_expected),
               "Zs differs from buggy Vs(corner) value");
}
```

---

### [R-304] [MODERATE] [fault_geometry.hpp new BP5 ctor:280-282] — Unconditional `num_zero_normal_fallbacks_ = 0` defeats the R-001 SAFS-mode pre-flight guard

**Category:** ASSUMPTION / EDGE_CASE

**Description:**
The new BP5 `FaultGeometry` ctor (Phase 5a) explicitly sets the fallback counters to zero at the end of the body:

```cpp
// fault_geometry.hpp:280-282
num_dof_basis_fallbacks_   = 0;
num_zero_normal_fallbacks_ = 0;
num_t1_fallbacks_          = 0;
```

These counters are populated by the legacy ctor's `ComputePerDOFCoordsAndBasis_(domain_op)` to count fault DOFs whose basis was zero-normal or t1-degenerate.  The legacy path uses these as a sanity check that the fault basis is well-formed.

The SAFS R-001 pre-flight guard in `spatial_stress.cpp:apply_csm_impl` is:

```cpp
MFEM_VERIFY(geom.NumZeroNormalFallbacks() == 0,
            "ApplyCsmStressSidecar: FaultGeometry reports "
            << geom.NumZeroNormalFallbacks()
            << " zero-normal fallback DOF(s); the SAFS-mode "
            "projection would silently produce sigma_n = 0 "
            "at those DOFs.  Fix the fault-DOF basis (see "
            "ElasticityDomainOperator::GetFaultDOFBasis) before "
            "applying the CSM sidecar.");
```

When the new ctor is used (the WHOLE POINT of Phase 5a — avoid the throw-away `ElasticityDomainOperator`), this counter is forced to 0 regardless of whether the supplied `dof_basis` is well-formed.  The guard becomes a tautology.  Any caller that passes a degenerate basis (zero-vector normal, near-zero t1, etc.) silently passes the SAFS pre-flight and gets `n · σ · n = 0` everywhere.

The ctor's docstring claims:

> The caller (ComputePerDOFCoordsAndBasisFromWave) is expected to deliver an already-Gram-Schmidt-orthonormalised basis; verifying that here would be a noisy assertion in tests.

But (a) `ComputePerDOFCoordsAndBasisFromWave` does NOT exist yet (deferred per `spatial_setup.hpp:14-22`), and (b) checking that `||n_i|| > 0` and `|t1_i · n_i| < 1e-12` for each DOF is exactly the same cost as the legacy ctor's fallback counting — not "noisy" in any meaningful sense.

**Trigger:**
Driver uses the new BP5 ctor and supplies a `dof_basis` with any degenerate column.  Could come from a hand-rolled SAFS basis construction, a buggy `ComputePerDOFCoordsAndBasisFromWave` implementation when it lands, or a serialization round-trip that lost precision.

**Actual behavior:**
SAFS R-001 guard passes; projection silently yields zero σ_n at the degenerate DOFs.

**Expected behavior:**
Either (a) walk the supplied `dof_basis` columns and increment `num_zero_normal_fallbacks_` for each degenerate normal (norm < 1e-10), OR (b) hard-`MFEM_VERIFY` the basis is well-formed at ctor time.

**Suggested fix:**
Replace the unconditional zeroing with a walk over the supplied basis:

```diff
-     // No domain operator => no zero-normal / t1 fallback counters
-     // populated by ComputePerDOFCoordsAndBasis_.  The caller
-     // (ComputePerDOFCoordsAndBasisFromWave) is expected to deliver
-     // an already-Gram-Schmidt-orthonormalised basis; verifying that
-     // here would be a noisy assertion in tests.  SAFS-mode consumers
-     // gate on NumZeroNormalFallbacks() == 0; with no domain operator
-     // to count, the counters stay at their default 0.
-     num_dof_basis_fallbacks_   = 0;
-     num_zero_normal_fallbacks_ = 0;
-     num_t1_fallbacks_          = 0;
+     // Walk the supplied basis and count any degenerate normal /
+     // t1 columns.  The SAFS R-001 pre-flight guard in
+     // ApplyCsmStressSidecar gates on NumZeroNormalFallbacks() == 0,
+     // so silently zeroing this counter would let a malformed basis
+     // pass the guard and project to zero everywhere.
+     num_dof_basis_fallbacks_   = 0;
+     num_zero_normal_fallbacks_ = 0;
+     num_t1_fallbacks_          = 0;
+     constexpr real_t k_norm_tol = static_cast<real_t>(1e-10);
+     for (int i = 0; i < N; ++i)
+     {
+        const real_t nx = dof_basis(0, i), ny = dof_basis(1, i),
+                     nz = dof_basis(2, i);
+        const real_t t1x = dof_basis(3, i), t1y = dof_basis(4, i),
+                     t1z = dof_basis(5, i);
+        const real_t n_norm  = std::sqrt(nx*nx + ny*ny + nz*nz);
+        const real_t t1_norm = std::sqrt(t1x*t1x + t1y*t1y + t1z*t1z);
+        if (n_norm < k_norm_tol)
+        {
+           ++num_zero_normal_fallbacks_;
+           ++num_dof_basis_fallbacks_;
+        }
+        if (t1_norm < k_norm_tol)
+        {
+           ++num_t1_fallbacks_;
+        }
+     }
```

**Test case:**
```cpp
// tests/unit/test_spatial_setup.cpp — add S-8
static void S_8_new_ctor_counts_zero_normal()
{
   // Build a per-DOF basis where DOF 0 has a zero normal (badly built).
   Vector dof_coords_3d(3);  dof_coords_3d = 0.0;
   DenseMatrix dof_basis(9, 1);
   dof_basis = 0.0;
   // (n stays at 0,0,0; t1 = (0,1,0); t2 = (0,0,1) — degenerate)
   dof_basis(4, 0) = 1.0;
   dof_basis(8, 0) = 1.0;
   Array<int> dof_to_elem(1);  dof_to_elem = 0;
   BP5Params params;
   FaultGeometry<mfem::Mesh> geom(params, dof_coords_3d, dof_basis,
                                  dof_to_elem, 1);
   TEST_ASSERT(geom.NumZeroNormalFallbacks() == 1,
               "new ctor must count degenerate normals");
}
```

---

### [R-305] [MODERATE] [wave_operator.hpp:WaveOperator] — `WaveOperator(MaterialField, BoundaryConfig)` ctor MISSING (plan §Phase H.1)

**Category:** DEVIATION

**Description:**
Plan §Phase H §"Files to Modify" item 1 mandates:

> Add NEW constructor:
> ```cpp
> WaveOperator(MeshType& mesh, int order,
>              const MaterialField& material,
>              const BoundaryConfig& bc);
> ```
> The existing scalar-material ctor is preserved verbatim.

`wave_operator.hpp:101-104` shows only the existing scalar ctor:

```cpp
WaveOperator(MeshType &mesh, int order,
             real_t lambda, real_t mu, real_t rho,
             const BoundaryConfig &bc);
```

No `MaterialField`-taking overload exists.  The SAFS driver (per plan §Phase 4 step 9) would call:

```cpp
WaveOperator<ParMesh> wave(pmesh, cfg.mesh.order, material, bc);
```

This would fail to compile.  The driver can only get a `WaveOperator` via the scalar ctor with seed material — then must manually `SetGodunovFluxPool(&pool)` — which itself is a no-op per R-301.

**Trigger:**
SAFS driver implementation (Phase 4).

**Actual behavior:**
The `MaterialField` ctor signature does not exist; Phase 4 cannot follow the plan's recommended pattern.

**Expected behavior:**
The new ctor exists and internally walks the mesh, evaluates `MaterialField::EvalAt` per element to build a per-element `(λ, μ, ρ)` table, then builds the `GodunovFluxPool` internally.

**Suggested fix:**
Either add the ctor stub now and document its no-op semantics, OR explicitly defer it in the plan and reflect that in the file comment so the next implementer knows where to land it.  Stub:

```diff
@@ wave_operator.hpp public ctors @@
   WaveOperator(MeshType &mesh, int order,
                real_t lambda, real_t mu, real_t rho,
                const BoundaryConfig &bc);

+  /// Phase H.1 of spatial_dynamic_rupture_plan.md (rev-3): heterogeneous-
+  /// material WaveOperator ctor.  ABORTS until Phase H.2/H.3/H.5 land —
+  /// see SetGodunovFluxPool for the rationale.  When implemented, this
+  /// ctor will evaluate `material.EvalAt(e, T, ip_centroid)` per element,
+  /// build a GodunovFluxPool, populate the heterogeneous-CFL path, and
+  /// exchange bi-material shared-face neighbour material via MPI.
+  WaveOperator(MeshType &mesh, int order,
+               const MaterialField &material,
+               const BoundaryConfig &bc);
```

And in `.cpp`:

```cpp
template <typename MeshType>
WaveOperator<MeshType>::WaveOperator(MeshType& /*mesh*/, int /*order*/,
                                     const MaterialField& /*material*/,
                                     const BoundaryConfig& /*bc*/)
   : WaveOperator(/*delegating*/...)
{
   MFEM_ABORT("WaveOperator(MaterialField, BoundaryConfig): not yet "
              "implemented (Phase H.1 stub).  Use the scalar ctor + "
              "SetGodunovFluxPool path once Phase H.2/H.3/H.5 land.");
}
```

**Test case:**
```cpp
// tests/unit/test_phaseh_wave_operator_constant_parity.cpp (plan-mandated; MISSING)
static void T_PHASEH_constant_material_parity()
{
   // Build WaveOperator with MakeConstant material; assert Mult output
   // bit-identical to WaveOperator built with the scalar ctor on the
   // same (lambda, mu, rho).  When R-305 is fixed, this test passes.
}
```

---

### [R-306] [MODERATE] [wave_operator.inl:ComputeMaxDt] — Heterogeneous-CFL path missing (plan §Phase H.3)

**Category:** DEVIATION

**Description:**
Plan §Phase H.3 mandates extending `ComputeMaxDt` to compute per-element CFL and `MPI_Allreduce(MIN)` when the operator is in heterogeneous mode:

> Heterogeneous CFL: `ComputeMaxDt(cfl)` walks every element and `MPI_Allreduce(MIN)`s the per-element `cfl · h_e / c_p,e`.

`wave_operator.inl:4964-5002` shows only the scalar-material formula:

```cpp
return cfl_mixed_flux_factor * cfl * h_min_ / flux_.GetCp();
```

`flux_.GetCp()` returns the scalar c_p set at construction.  For heterogeneous material, the actual stability bound is `min_e (cfl · h_e / c_p,e)`, which can be ORDERS OF MAGNITUDE smaller than `cfl · h_min / c_p_scalar` if even one element has very high c_p (e.g., basement rock).  Running with the scalar formula will violate CFL in those elements and blow up.

**Trigger:**
Heterogeneous-material run (once R-305 + R-301 are fixed).  Currently latent because R-305 prevents construction.

**Suggested fix:**
Add a heterogeneous branch that walks all elements:

```diff
+   if (flux_pool_)
+   {
+      real_t dt_min = std::numeric_limits<real_t>::infinity();
+      for (int e = 0; e < ne_; ++e)
+      {
+         const real_t cp_e = flux_pool_->At(e).GetCp();
+         const real_t h_e  = ElementCharacteristicLength_(e);  // existing helper
+         dt_min = std::min(dt_min, cfl * h_e / cp_e);
+      }
+      real_t dt_global = dt_min;
+#ifdef MFEM_USE_MPI
+      if (auto* pmesh = dynamic_cast<ParMesh*>(&mesh_))
+      {
+         MPI_Allreduce(&dt_min, &dt_global, 1, MPITypeMap<real_t>::value(),
+                       MPI_MIN, pmesh->GetComm());
+      }
+#endif
+      return cfl_mixed_flux_factor * dt_global;
+   }
    return cfl_mixed_flux_factor * cfl * h_min_ / flux_.GetCp();
```

(Note: this references `ElementCharacteristicLength_` which may not exist yet; adapt to whatever per-element h accessor the codebase provides.  If none, this finding needs its own ticket.)

**Test case:**
```cpp
// tests/unit/test_phaseh_wave_operator_layered.cpp (plan-mandated; MISSING)
static void T_PHASEH_layered_cfl_min()
{
   // Two-layer Vs profile; check ComputeMaxDt returns
   // min(cfl * h / Vs_top, cfl * h / Vs_bot)
   // not   cfl * h / Vs_scalar.
}
```

---

### [R-307] [MODERATE] [tests/unit/test_phaseh_*] — Plan-mandated Phase H test files MISSING

**Category:** DEVIATION

**Description:**
Plan §Phase H "Files to Create" lists FOUR test files; only TWO exist:

| Plan-required test file | Status |
|--------------------------|--------|
| `test_phaseh_wave_operator_constant_parity.cpp` (4 tests) | **MISSING** |
| `test_phaseh_wave_operator_layered.cpp`           (3 tests) | **MISSING** |
| `test_phaseh_godunov_flux_pool.cpp`               (4 tests) | EXISTS (4 tests, but 2 plan-required ones missing) |
| `test_phaseh_lsw_forced_rupture.cpp`              (6 tests) | EXISTS (6 tests, matches plan) |

The existing `test_phaseh_godunov_flux_pool.cpp` covers:
- dedup correctness (P-1) ✓
- distinct triples (P-2) — not plan-required but useful
- rho = 0 aborts (P-3) ✓
- pool.At matches fresh GodunovFlux (P-4) — not plan-required

The plan asks specifically for:
1. dedup correctness ✓
2. boundary correctness (rho = 0) ✓
3. **MPI consistency (4 ranks see identical pool size given identical input)** — MISSING
4. **memory budget (NumUniqueTriples() ≤ expected_unique on the SAFS 1000 m cvmh fixture)** — MISSING

The TWO missing test files are the only direct verification that Phase H.1/H.2/H.3 actually work end-to-end.  Without them, the heterogeneous WaveOperator code path is entirely uncovered — exactly the situation that lets R-301/R-302 ship without anyone noticing the dispatch is dead.

**Trigger:**
N/A — coverage gap.

**Suggested fix:**
Create the two missing test files.  Skeletons:

```cpp
// tests/unit/test_phaseh_wave_operator_constant_parity.cpp
static void T_HET_const_material_parity_serial() { /* ... */ }
static void T_HET_const_material_parity_parallel_np2() { /* ... */ }
static void T_HET_const_material_parity_with_pml() { /* ... */ }
static void T_HET_const_material_byteflows_match_scalar() { /* ... */ }

// tests/unit/test_phaseh_wave_operator_layered.cpp
static void T_HET_two_layer_dQdt_differs_from_scalar() { /* ... */ }
static void T_HET_two_layer_arrival_time_differs() { /* ... */ }
static void T_HET_two_layer_cfl_min_is_smaller_layer() { /* ... */ }
```

These require R-305 (new ctor) and R-301 (live dispatch) to actually compile + pass, so this finding cannot be fully resolved until those land.  Tracking it separately keeps the gap visible.

Also extend `test_phaseh_godunov_flux_pool.cpp`:

```diff
+ // P-5  MPI consistency: 4 ranks with identical per-element input
+ //      build identical pools.
+ static void P_5_mpi_consistency() { ... }
+
+ // P-6  memory budget on SAFS 1000m cvmh fixture (skip if not present)
+ static void P_6_safs_memory_budget() { ... }
```

**Test case:**
(See above — the missing test files ARE the test cases.)

---

### [R-308] [MODERATE] [test_spatial_setup.cpp] — No layered-material test; S-1 / S-2 / S-3 all use `MakeConstant` which masks R-303

**Category:** EDGE_CASE / QUALITY (test coverage)

**Description:**
Plan §Phase 5 §Acceptance Test 2 mandates:

> 2. Layered material (`MakeCoefficient` with a depth-step Vs) produces distinct impedances at z = -1km vs z = -10km (relative diff > 1e-3).

`test_spatial_setup.cpp` has 6 tests (S-1..S-6).  Every test that calls `InitializeFaultDOFs_Spatial` uses `MaterialField::MakeConstant(...)`:

```cpp
// S-1 line 127
auto material = MaterialField::MakeConstant(lam, mu, rho);
// S-2 line 177
auto material = MaterialField::MakeConstant(32e9, 32e9, 2670.0);
// S-3 line 214
auto material = MaterialField::MakeConstant(32e9, 32e9, 2670.0);
```

`MakeConstant` materials short-circuit `EvalAt` and return the scalar triple regardless of IP location.  So the corner-vs-centroid bug (R-303) is INVISIBLE to every one of these tests — they pass byte-identically whether `ip.Init(0)` (corner) or `Geometries.GetCenter(...)` (centroid) is used.

The plan-mandated layered-material test that would catch R-303 is **missing**.

**Trigger:**
Reintroduce `ip.Init(0)` into `seed_static_dof_fields` (or fail to fix R-303) — every S-* test still passes.

**Suggested fix:**
Add the layered-material test (skeleton in the R-303 finding above as S-7).  Use `MaterialField::MakeCoefficient` with a z-varying Vs Coefficient (`MuFromZ` / `LambdaFromZ` / `RhoConst` mirroring the F-3 / R-9 pattern in `test_spatial_friction_resolver.cpp`).

---

### [R-309] [MODERATE] [spatial_setup.hpp:InitializeFaultDOFs_Spatial signature] — `fault_dof_ip(i)` cache from new BP5 ctor is NEVER consumed

**Category:** DEVIATION / EDGE_CASE

**Description:**
Phase 5a added `FaultGeometry::fault_dof_ip(int i)` and the new BP5 ctor accepts `const std::vector<IntegrationPoint>& dof_ips` so the driver can supply per-DOF reference IPs.  The plan's intent is that downstream consumers use these to evaluate material exactly at the fault DOF (R-111 / round-1 R-002 long-term fix).

`InitializeFaultDOFs_Spatial` does NOT accept an IP array.  Its signature only takes `dof_to_elem` + `material` + `mesh` — `seed_static_dof_fields` then uses `ip.Init(0)` (R-303) instead of consuming the cached IPs.  Even if a driver populates `fault_dof_ip` correctly via the new ctor, the impedance init silently ignores them.

The docstring in `spatial_setup.hpp:59-63` acknowledges this in passing:

> when Phase 5a's `fault_dof_ip(i)` is wired in by the driver, the caller should pass that IP directly via the eight-arg overload below.

But there IS no eight-arg overload.  The comment claims a feature that does not exist.

**Trigger:**
Driver builds `FaultGeometry` with non-empty `dof_ips`, expects `InitializeFaultDOFs_Spatial` to use them.

**Actual behavior:**
`dof_ips` data sits unused; impedances evaluated at `ip.Init(0)` (corner).

**Expected behavior:**
Add the documented overload that accepts `const std::vector<IntegrationPoint>* dof_ips` and uses each entry per DOF when non-null.

**Suggested fix:**
Add the missing overload and the corresponding `seed_static_dof_fields` overload:

```diff
@@ spatial_setup.hpp @@
+ /// Overload that consumes the per-DOF reference IPs cached by the
+ /// new FaultGeometry BP5 ctor (R-111).  When `dof_ips != nullptr`,
+ /// material is evaluated at `(*dof_ips)[i]` per DOF instead of at
+ /// the element centroid.
+ template <typename MeshT>
+ inline void InitializeFaultDOFs_Spatial(
+    std::vector<DOFData>&             dof_data,
+    int                                ndof,
+    const Array<int>&                  dof_to_elem,
+    const MaterialField&               material,
+    MeshT&                             mesh,
+    const SlipWeakeningPerDOFParams&   lsw,
+    const Vector&                      tau_pre,
+    const Vector&                      sigma_n_eff,
+    const Vector&                      T_forced_s,
+    const Vector&                      t0_decay_s,
+    const std::vector<IntegrationPoint>* dof_ips);
```

Also fix the misleading "eight-arg overload" comment in `seed_static_dof_fields`.

---

### [R-310] [MODERATE] [fault_geometry.hpp new BP5 ctor:257-271] — `coords_x2_` / `coords_x3_` populated with WORLD (x, z) confuses `ComputeBP5Params`

**Category:** BUG (latent)

**Description:**
The new BP5 ctor populates `coords_x2_` and `coords_x3_` with world coordinates:

```cpp
// fault_geometry.hpp:257-266
coords_x2_.SetSize(N);
coords_x3_.SetSize(N);
depths_.SetSize(N);
for (int i = 0; i < N; ++i)
{
   const real_t z = dof_coords_3d(3 * i + 2);
   coords_x2_(i) = dof_coords_3d(3 * i + 0);  // ← world x
   coords_x3_(i) = z;                          // ← world z
   depths_(i)    = z;
}
ComputeBP5Params();
```

The legacy BP5 ctor uses `coords_x2_` / `coords_x3_` as the **2-D fault-local coords** (along-strike and along-dip on the BP5 planar fault), and `ComputeBP5Params` reads them to populate the BP5 analytic per-DOF parameters `a_values_(x2, x3)`, `Dc_values_(x2, x3)`, `V_init_values_(x2, x3)`, plus the BP5 analytic `tau_pre_`.

Passing WORLD (x, z) of a SAFS curvilinear fault into BP5 analytic functions produces NONSENSE — e.g. a SAFS DOF at x = 480_000 m UTM gives `a_values_` lookup at "along-strike 480 km" which is far outside the BP5 domain `[-Lf/2, +Lf/2]`.

The docstring acknowledges the issue:

> These will be overwritten per-DOF by `InitializeFaultDOFs_Spatial` (Phase 5c) before the SAFS driver consumes them; the legacy BP5 ctor's 2-D (x2, x3) convention is intentionally NOT used in this path.

But the "will be overwritten" hand-wave assumes:
- The SAFS driver always calls `InitializeFaultDOFs_Spatial` BEFORE any consumer reads `a_values_` / `Dc_values_` / `V_init_values_` / `tau_pre_`.
- No diagnostic, no `--dry-run`, no checkpoint hand-back, and no logging touches these arrays in between.

This is fragile.  If `--dry-run` happens to log per-DOF BP5 params for sanity, the user sees garbage.  If a future test reads `a_values_` after the new ctor + before init, the test sees garbage.

**Trigger:**
Construct `FaultGeometry` via the new ctor; read `a_values_`, `Dc_values_`, `V_init_values_`, or `tau_pre_` before `InitializeFaultDOFs_Spatial` overwrites them.

**Actual behavior:**
Out-of-range BP5 analytic values.

**Expected behavior:**
Either (a) skip `ComputeBP5Params()` entirely in the new ctor (since the SAFS pipeline overwrites everything), OR (b) initialise the BP5 analytic arrays to NaN sentinels so any premature read fails loudly.

**Suggested fix:**
Option (a) is simpler and matches the comment's intent:

```diff
@@ fault_geometry.hpp new BP5 ctor body @@
       coords_x2_.SetSize(N);
       coords_x3_.SetSize(N);
       depths_.SetSize(N);
       for (int i = 0; i < N; ++i)
       {
          const real_t z = dof_coords_3d(3 * i + 2);
-         coords_x2_(i) = dof_coords_3d(3 * i + 0);
-         coords_x3_(i) = z;
+         // 2-D fault coords are NOT defined for SAFS curvilinear fault;
+         // leave coords_x2_ / coords_x3_ at NaN sentinels so any
+         // premature read of BP5 analytic params fails loudly.
+         coords_x2_(i) = std::numeric_limits<real_t>::quiet_NaN();
+         coords_x3_(i) = std::numeric_limits<real_t>::quiet_NaN();
          depths_(i)    = z;
       }
-
-      // Seed BP5 analytic per-DOF arrays (a, eta, Dc, V_init,
-      // tau_pre).  SAFS callers (Phase 5c) overwrite these per DOF
-      // before consumption.
-      ComputeBP5Params();
+      // Skip ComputeBP5Params: the SAFS driver overwrites every
+      // per-DOF analytic value via InitializeFaultDOFs_Spatial.
+      // Sentinel-NaN coords_x2_ / coords_x3_ above ensures any read
+      // of BP5 analytic params before the overwrite fails loudly.
```

---

### [R-311] [LOW] [godunov_flux_pool.cpp:Build] — "internal bug" guard at end is unreachable

**Category:** QUALITY

**Description:**
`godunov_flux_pool.cpp:102-104`:

```cpp
MFEM_VERIFY(ne == 0 || !unique_fluxes_.empty(),
            "GodunovFluxPool::Build: ne > 0 but produced 0 unique "
            "triples (internal bug)");
```

For `ne > 0`, the loop body iterates at least once.  Each iteration enters either the `if (it == key_to_idx.end())` branch (creates a new flux + entry) or the `else` branch (reuses an existing entry).  The first iteration ALWAYS enters the `if` branch (the map is empty), so `unique_fluxes_.size() >= 1` after iteration 0.  Therefore `ne > 0 && unique_fluxes_.empty()` cannot occur — the assertion is unreachable.

**Suggested fix:**
Remove the dead guard:

```diff
-   MFEM_VERIFY(ne == 0 || !unique_fluxes_.empty(),
-               "GodunovFluxPool::Build: ne > 0 but produced 0 unique "
-               "triples (internal bug)");
```

---

### [R-312] [LOW] [godunov_flux_pool.cpp:round_sig] — Unprotected `std::pow` for very large exponents

**Category:** EDGE_CASE

**Description:**
`round_sig`:

```cpp
const real_t scale = std::pow(static_cast<real_t>(10.0),
                              sig - 1 - std::floor(std::log10(mag)));
```

For very small `mag` like `1e-300` and `sig = 6`, the exponent is `6 - 1 - (-300) = 305`.  `scale = 1e305` is still representable in double.  But for `mag = 1e-310` (subnormal), the exponent is 315, and `scale = 1e315` overflows to +INF.  Then `v * scale = INF`, `round(INF) / INF = NaN`, and the key becomes `"nan|...|..."` — all such inputs collapse into one bucket.

Not a real-world issue for SAFS material (Vp, Vs, rho are all O(1) – O(1e10)), but worth guarding against.

**Suggested fix:**
Clamp the exponent to a safe range, or fall back to using the unrounded value when out of range:

```diff
 real_t round_sig(real_t v, int sig)
 {
    if (v == 0.0 || !std::isfinite(v)) { return v; }
    const real_t mag   = std::abs(v);
+   const real_t log_mag = std::log10(mag);
+   // Guard against catastrophic overflow in scale for subnormal mag.
+   if (!std::isfinite(log_mag) || std::abs(log_mag) > 290.0)
+   {
+      return v;  // mag at/near the dbl-precision range edge; bypass round.
+   }
    const real_t scale = std::pow(static_cast<real_t>(10.0),
-                                 sig - 1 - std::floor(std::log10(mag)));
+                                 sig - 1 - std::floor(log_mag));
    const real_t r = std::round(v * scale) / scale;
    return r;
 }
```

---

### [R-313] [LOW] [test_phaseh_lsw_forced_rupture.cpp F-1..F-6] — Tests exercise the helper directly, NOT the wave_operator dispatch

**Category:** EDGE_CASE / QUALITY (coverage)

**Description:**
All six tests in `test_phaseh_lsw_forced_rupture.cpp` call the helper / new method directly without going through `WaveOperator::Mult` → fault-dispatch.  With R-302 present (no LSW_ForcedRupture arm in wave_operator.inl), these tests would still pass byte-identically — the dead-code dispatch is invisible to the test suite.

**Suggested fix:**
Add one integration test that exercises the dispatch end-to-end (paired with the R-302 fix):

```cpp
static void F_7_dispatch_routes_through_forced_rupture()
{
   // Build minimal WaveOperator, FaultFaceFlux, FaultGeometry; populate
   // DOFData with T_forced = 0, t0_decay = 0.5, lsw_d_c = 0.4, slip = 0.
   // Set FaultFrictionLaw::LSW_ForcedRupture.
   // Call wave.AdvanceADER(t = 0.25, dt = 1e-3).
   // Read dof_data[0].tau1_corr; assert it matches the value computed
   // by EvaluateADER_LSW_ForcedRupture (not by EvaluateADER which would
   // crash or return garbage on LSW-init data).
}
```

---

### [R-314] [LOW] [drivers/] — Phase 4 `spatial_dyn_driver.cpp` still missing

**Category:** DEVIATION (scope-incomplete)

**Description:**
Plan §Phase 4 mandates `drivers/spatial_dyn_driver.cpp` (~1100 LOC) — the main executable.  `ls drivers/spatial*` returns no matches.  Phase 4 cannot complete until R-301/R-302/R-303/R-305/R-306 are resolved (the driver would compile but produce wrong physics).  Tracking informational only.

---

### [R-315] [LOW] [jobs/safs/, spatial/code/scripts/] — Phase 6 sbatches + verification scripts mostly missing

**Category:** DEVIATION (scope-incomplete)

**Description:**
Plan §Phase 6 mandates 3 sbatches + 2 verification scripts.  Present on disk:
- `jobs/safs/safs_smoke_8N_400r_dev.sbatch` (1 of 3)
- `spatial/code/scripts/verify_constant_tensor_projection.py` (Phase 3b helper; not the Phase 6 smoke verifier)

Missing:
- `jobs/safs/spatial_dyn_production_normal_48hr_safs.sbatch`
- `jobs/safs/spatial_dyn_cvm_compare_dev_2hr_safs.sbatch`
- `spatial/code/scripts/verify_spatial_dyn_smoke_safs.py`
- `spatial/code/scripts/compare_spatial_dyn_velocity_models.py`

Informational only; Phase 6 cannot complete until Phase 4 driver lands.

---

## Summary

- Critical issues: **3**  (R-301 flux_pool dispatch dead, R-302 LSW_ForcedRupture dispatch arm missing, R-303 ip.Init(0) regression in spatial_setup)
- Moderate issues: **7**  (R-304 zero-normal counter bypass, R-305 missing MaterialField ctor, R-306 missing heterogeneous CFL, R-307 missing Phase H test files, R-308 missing layered-material test, R-309 fault_dof_ip unused, R-310 BP5 ctor coords_x2 confusion)
- Low issues: **5**  (R-311 dead guard, R-312 round_sig overflow, R-313 helper-only tests, R-314 Phase 4 missing, R-315 Phase 6 missing)
- Plan compliance: **PARTIAL** — Phase H is 30% (helper + enum + method body exist; pool + dispatch + CFL + ctor all missing).  Phase 5 is 80% (checkpoint extension + new ctor + setup helpers exist; coords_x2 confusion + zero-normal bypass + fault_dof_ip unused).  Phase 4 and Phase 6 still not started.
- Verdict: **FAIL — must fix R-301, R-302, R-303 before any SAFS dynamic-rupture run can produce correct physics.**  R-301 turns heterogeneous-material runs into scalar-material runs.  R-302 routes the SAFS forced-rupture friction law through rate-and-state Brent.  R-303 re-introduces the round-1 corner-vs-centroid impedance bug.  All three are silent at compile and at unit-test time; they only manifest as wrong slip / wrong wave propagation in production.

## Unreviewed Areas

- **`wave_operator.inl`** beyond the two LSW dispatch sites — the file is ~5000 LOC; a full audit of every flux call site to verify no other `flux_` reference needs `flux_pool_` is required before R-301 can be cleared.
- **Phase H bi-material MPI exchange** — not implemented; nothing to review.
- **Phase 4 driver source** — not present; nothing to review.
- **TPV byte-exact full regression suite** (`make test-tpv104`, `make test-tpv205`, `make test-tpv-bp5-byte-exact`) — only `seas_test_tpv104_checkpoint` (164/164) and `seas_test_compute_safs_params` (13/13) were run.  The full per-driver byte-exact gates listed in §TPV Benchmark Isolation should be run before Phase H merges.
- **`seas_test_bp5_fault_operator`** — produces `Verification failed: Tagged attr-3 face key (6,7,15) was not recovered as a fault face` and exits 0 (silent failure).  Probably pre-existing on the SAFS branch (CLAUDE.md notes BP5 / TPV102 intermittent failures on this branch), but not confirmed against pre-Phase-H state in this review.
- **MPI behaviour** of all new code paths — every test ran serial or implicitly serial; the plan's `T-PHASEH-BIMATERIAL-SHARED-FACE-MPI` and `T-CHECKPOINT-DRIVER-TAG-BACK-COMPAT` regressions are not exercised.
