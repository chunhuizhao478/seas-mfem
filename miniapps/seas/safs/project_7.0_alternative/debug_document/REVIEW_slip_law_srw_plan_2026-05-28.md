# Code Review: Slip-Law + SRW plan (2026-05-28)

## Review Scope
- Plan: `safs/project_7.0_alternative/debug_document/spatial_dynamic_rupture_slip_law_srw_plan_2026-05-28.md`
- Files inspected (the code the plan would touch / depend on):
  - `friction/slip_law_srw_psi.hpp`
  - `dynamic/tpv104_substep_iterator.{hpp,cpp}`
  - `dynamic/friction_iterator.hpp`, `dynamic/friction_iterator_factory.cpp`
  - `dynamic/fault_face_flux.{hpp,cpp}` (`ComputeStageState`, `DOFData`)
  - `dynamic/spatial_setup.hpp` (`InitializeFaultDOFs_Spatial_RS`, `SeedEquilibriumPsi_RS`)
  - `spatial/code/spatial_friction.hpp`
  - `config/...depthprofile_Dc010_nuc8km_500m.toml`
- Domain context: `CLAUDE.md` (Brent-only friction solver; friction = extreme-care),
  `safs-nucleation-sigma-n-and-lengthscales` memory (depth profile → VW→VS arrest at 11 km).

This is a review of the **plan before implementation**. Findings are design
defects that would become code bugs if the plan is implemented verbatim.

## Findings

### [R-001] CRITICAL [plan §4.3/§4.4 + tpv104_substep_iterator.cpp:413–419,655–661] — Depth-varying b(z) is silently flattened to a scalar

**Category:** BUG (silent physics error)

**Description:**
Decision #3 keeps the depth-varying profile so that `b(z)` falls to give the
VW→VS transition (a−b=0) at ~11 km — the mechanism that **arrests down-dip
rupture** above the 16.6 km fault bottom. Decision #4 reuses
`Tpv104SubStepIterator`. But that iterator's ψ-update sources `b` from the
**scalar** `state_evo_.GetB()`, not from the per-DOF `d.b`:

```cpp
// tpv104_substep_iterator.cpp:413  (and the SeisSol variant at :655)
d.psi = UpdateStateAnalyticSlipLawSRW(d.psi, s.V_abs,
                                      d.Dc, dt_sub,      // per-QP L  ✓
                                      V_w[i], d.a,       // per-QP V_w, a ✓
                                      state_evo_.GetB(), // <-- SCALAR b  ✗
                                      state_evo_.GetV0(),
                                      state_evo_.GetF0(),
                                      state_evo_.GetMuW());
```

`SlipLawSRWPsi::GetB()` returns the single scalar passed to its constructor.
TPV104 is correct here (global b=0.014). For SAFS, `InitializeFaultDOFs_Spatial_RS`
(spatial_setup.hpp:340–342) sets per-DOF `d.a=rs.a`, **`d.b=rs.b`**, `d.Dc=rs.Dc`
— the per-DOF b is available and used by the aging path — but the SRW iterator
**ignores `d.b`** and uses the scalar instead. Result: every fault QP evolves ψ
with one constant b; the depth profile's b(z) (hence the VW→VS arrest) is
discarded. Worse, `f_LV = max(0, f0 − (b−a)·ln(V/V0))` then mixes **scalar b**
with **per-QP a (d.a)** — physically inconsistent across depth.

The plan's §4.3 ("construct `SlipLawSRWPsi` from … b …") and §4.4 ("byte-for-byte
identical except the nucleation call") both miss this: the b in the constructor
is *not* the per-QP b the physics needs, and the additive callback overload must
change **more than** the nucleation line.

**Trigger:** any SAFS SRW run with the depth profile enabled (the chosen config).
Manifests as: no down-dip rupture arrest; uniform velocity-weakening to the fault
bottom; wrong nucleation/propagation length scales.

**Actual behavior (if implemented as written):** ψ evolves with constant b
everywhere; b(z) silently ignored.

**Expected behavior:** ψ evolves with per-QP `d.b` (= the resolved depth-profile
b at each fault DOF).

**Suggested fix:**
In the new additive nuc_callback overload of `AdvanceWithSubStepStates`
(§4.4), source b from the per-DOF `d.b`, not the scalar accessor — i.e. the
overload differs from the original in **two** places (nucleation call *and* the b
argument):
```diff
  d.psi = UpdateStateAnalyticSlipLawSRW(d.psi, s.V_abs,
                                        d.Dc, dt_sub,
                                        V_w[i], d.a,
-                                       state_evo_.GetB(),
+                                       d.b,                 // per-QP depth-profile b
                                        state_evo_.GetV0(),
                                        state_evo_.GetF0(),
                                        state_evo_.GetMuW());
```
Add a precondition in the overload: `MFEM_VERIFY(std::isfinite(d.b) && d.b > 0.0, ...)`
(DOFData.b defaults to NaN). Update plan §4.3 to note the `SlipLawSRWPsi` b arg is
**unused in the SAFS path** (b comes from `d.b`); pass `b_default` only as a
placeholder, and update §4.4 to state the overload changes the b source too.
(The existing original `AdvanceWithSubStepStates`/`Advance` keep `GetB()` →
TPV104 byte-exact preserved.)

**Test case:**
```cpp
// tests/unit/test_slip_law_srw_depth_b.cpp
TEST(SlipLawSRWAdapter, PerDofBIsHonored) {
  // Two fault QPs with the SAME a, Dc, V, psi, V_w but DIFFERENT d.b
  // (e.g. d.b = 0.0261 at the patch vs d.b = 0.005 near the VW->VS edge).
  DOFData d_patch = make_dof(/*a=*/0.0127, /*b=*/0.0261, /*Dc=*/0.10, /*psi0=*/p0);
  DOFData d_edge  = make_dof(/*a=*/0.0127, /*b=*/0.0050, /*Dc=*/0.10, /*psi0=*/p0);
  // Drive one macro-step through SlipLawSRWFrictionIterator with identical Q/V_w.
  advance_one_step(d_patch); advance_one_step(d_edge);
  // If b were flattened to a scalar, psi would be identical. It must differ.
  EXPECT_GT(std::abs(d_patch.psi - d_edge.psi), 1e-9)
      << "per-DOF b(z) was flattened to a scalar (R-001)";
}
```

---

### [R-002] MODERATE [plan §4.7] — Parity test is ill-posed and will not catch R-001

**Category:** ASSUMPTION / inadequate test

**Description:**
§4.7 proposes a "bit-identical" parity test: drive `SlipLawSRWFrictionIterator`
and a bare `Tpv104SubStepIterator` with identical state + no-op nuc_callback and
assert identical DOFData. Two problems:
1. After the R-001 fix the adapter's overload uses `d.b`, while the bare
   iterator's original `Advance` uses `GetB()`. For a **depth-varying** b
   fixture they differ **by design** — the test would fail spuriously, or
   (if a constant b is used) it would pass while telling you nothing about the
   per-QP path.
2. The adapter passes `Method::Brent` (CLAUDE.md, §4.3), but `Advance`/
   `AdvanceWithSubStepStates` default to `Method::NewtonRaphsonStable`. Driving
   the bare iterator with the default → different solver → not bit-identical.

A "no-perturbation" parity test that uses constant b + matched Brent is fine as a
guard that the *adapter plumbing* adds nothing — but it cannot catch R-001
(constant b hides it). A **separate** depth-varying test (R-001's test case) is
required.

**Trigger:** running the §4.7 parity test as specified.

**Actual behavior:** false failure (depth-varying fixture) or false confidence
(constant fixture masks the b bug).

**Expected behavior:** (a) parity test uses **constant b** and **`Method::Brent`
on both sides** to verify the adapter adds no perturbation; (b) a distinct
depth-varying-b test (R-001) verifies per-QP b is honored.

**Suggested fix:** split §4.7 into two tests as above; pass `Method::Brent`
explicitly to the bare iterator in the parity test; never rely on the parity
test alone to cover the depth profile.

**Test case:** see R-001 test (depth-varying) + a constant-b/Brent parity test:
```cpp
TEST(SlipLawSRWAdapter, AdapterMatchesBareIteratorConstantB) {
  // constant b on all QPs, Method::Brent on BOTH sides, no-op nuc_callback.
  // Assert dof_data bit-identical after one macro-step.
}
```

---

### [R-003] MODERATE [POSSIBLE] [plan §4.3 + CLAUDE.md] — SRW friction solve uses Brent, but only Newton-stable is regression-validated for FL=103

**Category:** ASSUMPTION

**Description:**
§4.3 mandates `Method::Brent` (correct per CLAUDE.md and consistent with the
aging SAFS path). However, the TPV104 FL=103 path was validated against the
SeisSol reference traces with `NewtonRaphsonStable` (tpv104_substep_iterator.hpp:
106–142). The force solve's strength form is identical for aging and SRW, and the
aging SAFS path already runs Brent successfully on that form, so Brent is very
likely fine — but the plan *assumes* it without a verification step. A
non-converging or differently-converging Brent solve on the SRW ψ envelope would
surface as wrong V at high slip rate, only at run time.

**Trigger:** first SRW production run (no unit coverage of SRW+Brent today).

**Actual behavior:** unverified solver/law combination ships untested.

**Expected behavior:** an explicit SRW+Brent convergence smoke check before
production (e.g. a unit test solving a few representative (ψ, σ_n, a) points and
asserting the Brent root matches the strength balance to tolerance), and/or a
short local dry-run asserting bounded V.

**Suggested fix:** add to §8 (implementation order) an explicit "SRW+Brent
convergence check" step; add a unit test:
```cpp
TEST(SlipLawSRW, BrentRootMatchesStrengthBalance) {
  for (real_t psi : {0.3, 0.5, 0.6}) for (real_t Vstar : {1e-12, 1e-3, 1.0}) {
    // build stress from the analytic strength at Vstar, solve with Brent,
    // assert |V_solved - Vstar| / Vstar < 1e-6.
  }
}
```

---

### [R-004] LOW [POSSIBLE] [plan §2/§4.6 + spatial_setup.hpp:468–472] — "seed reused unchanged" is true only on the primary branch; the fallback branch is aging-specific

**Category:** BUG (latent; unreachable at the chosen V_init)

**Description:**
The plan states `SeedEquilibriumPsi_RS` is "law-agnostic … reused unchanged."
The **primary** branch (`InitialStatePsi`, tau_eff>0) is indeed the regularized
asinh inversion and is law-agnostic. The **fallback** branch (damping-dominated,
`tau_eff = tau0 − eta·V_init ≤ 0`) seeds the **aging-law** steady state
`ψ = f0 + b·ln(V0/V_init)`, which is *not* the SRW steady state
`ψ_ss = a·ln((2V0/V)·sinh(f_ss/a))`. For SAFS this branch is unreachable
(V_init=1e-12, eta·V_init ≈ 4.6e6·1e-12 ≈ 4.6e-6 Pa ≪ tau0 ~ tens of MPa), so it
never fires — but the plan's blanket "reused unchanged" hides the caveat.

**Trigger:** only if V_init or eta were raised until eta·V_init ≳ tau0 (not the
current config).

**Suggested fix:** annotate the plan: "seed reused unchanged; the tau_eff≤0
fallback is aging-specific but unreachable at V_init=1e-12 — do not raise V_init
into the damping-dominated regime without revisiting the SRW seed."

**Test case:**
```cpp
TEST(SeedEquilibriumPsiRS, SafsConfigUsesPrimaryBranch) {
  // V_init=1e-12, eta~4.6e6, tau0~3e7 -> tau_eff>0 -> InitialStatePsi branch.
  EXPECT_GT(tau0 - eta * 1e-12, 0.0);
}
```

---

### [R-005] LOW [plan §4.3] — `SetProductionMode()` is presented as required; it is a harmless no-op safety net here

**Category:** QUALITY (clarity, prevents an implementer misstep)

**Description:**
§4.3 says "call `state_evo_.SetProductionMode()` … production uses the per-QP
`_SRW` overloads." In fact `Tpv104SubStepIterator` evolves ψ via the **free
function** `UpdateStateAnalyticSlipLawSRW` and reads only the **non-virtual**
getters `GetB/GetV0/GetF0/GetMuW` — it never calls the base virtuals
(`Rate`/`SteadyState`/…) that production mode guards. So `SetProductionMode()` is
a defensive no-op on this path, not load-bearing. Stating it as required could
lead an implementer to expect the per-QP `_SRW` member overloads to be wired
(they are not used by this iterator), or to debug a non-existent dispatch.

**Suggested fix:** reword §4.3: "calling `SetProductionMode()` is a harmless
safety net (the iterator uses the free `UpdateStateAnalyticSlipLawSRW` + scalar
getters, never the base virtuals); set it anyway to fail loudly if a future edit
routes through a base virtual."

---

## Resolved during review (no finding)
- **Plan §4.4 `tau1_nuc` verification** — `ComputeStageState` (fault_face_flux.cpp:
  183–184, 320–321) folds **both** `data.tau1_nuc` and `data.tau2_nuc` into the
  total/corrected traction. The SAFS `gradual_overstress` dip+strike injection
  therefore works through the TPV104 iterator unchanged. (Moot for the current
  config anyway: `delta_tau_dip_pa = 0.0`.) The plan's verification step will
  pass.
- **Per-QP `a` and `Dc`** — `InitializeFaultDOFs_Spatial_RS` sets `d.a`, `d.Dc`
  per-DOF and the iterator consumes `d.a`/`d.Dc`; these are correctly per-QP.
- **`SlipLawSRWPsi` scalar `a` constructor arg** — unused in this path (iterator
  uses `d.a`); harmless placeholder. (Same now applies to `b` after R-001 fix.)

## Summary
- Critical issues: 1 (R-001)
- Moderate issues: 2 (R-002, R-003)
- Low issues: 2 (R-004, R-005)
- Plan compliance: PARTIAL — the chosen combination (decision #3 depth profile +
  decision #4 reuse iterator) is internally inconsistent as written; R-001 must
  be folded into §4.3/§4.4 before implementation.
- Verdict: **PASS WITH FIXES** — the plan is sound in structure; R-001 is a
  must-fix design correction (one extra one-line change in the new overload +
  two doc corrections). R-002/R-003 are must-do test/verification additions.
  R-004/R-005 are doc clarifications.

## Unreviewed Areas
- The exact line where the new nuc_callback overload is inserted in
  `tpv104_substep_iterator.cpp` (the .cpp body between the SEAS_DIAG blocks was
  read only via grep, not line-by-line) — confirm the overload is a faithful
  copy of the existing `AdvanceWithSubStepStates` body with exactly the two
  changes in R-001 (b source) and the nucleation call.
- The two sbatch edits (§5.2) were not line-diffed; the 250 m→500 m mesh-mode and
  CSV pre-flight changes should be re-reviewed once drafted.
