# Code Review (Round 2 — fresh audit of the corrected plan): fault-dealiasing, 2026-06-02

## Review Scope
- Plan: `miniapps/seas/document/fault_dealiasing_dev/seisol_overintegration_resample_speckle_2026-06-02.md` (post round-1 fixes)
- Reference code re-consulted this round:
  - SeisSol resample matrices: `codegen/matrices/dr_stroud_matrices_{3,4,5}.json` (re-verified: idempotent to 1e-14, **W-self-adjoint** `‖WR−(WR)ᵀ‖ ≈ 3e-18`, rank = `(N+1)(N+2)/2`). This **confirms** the round-1 fix Eq. 4.2 `R = V_N(V_NᵀWV_N)⁻¹V_NᵀW` is *exactly* SeisSol's operator — that finding is correctly resolved.
  - SEAS-MFEM: `dynamic/friction_substep_iterator.cpp:155-205` (LSW iterator), `dynamic/fault_state_channel.hpp:37-39`, `dynamic/tpv205_stations.hpp:204`, `dynamic/tpv31_stations.hpp:233`, `dynamic/fault_face_flux.hpp:27-140` (`DOFData`), `drivers/spatial_dyn_driver.cpp:217,925`, `tpv102/configs/tpv102_spatial.toml:49`.
- Domain context: `CLAUDE.md`, `miniapps/seas/CLAUDE.md` (gmsh `.msh` = tetrahedral; v2.2 parser note), round-1 fix report.

**Round-1 status:** all six round-1 findings (R-001…R-006) are correctly resolved in the current document — verified by re-reading the rewritten §3.3/3.5/4.1/4.2/4.3/5/6/7/8 and re-running the round-1 numeric checks (resample rank = degree N; worked example P₃ clean, P₄→c₂=−3/10). This round is a fresh adversarial pass and surfaces **new** issues that the corrected plan now *depends on*.

---

## Findings

### [R-001] CRITICAL [§4.2 / §6 Phase 3 sketch — LSW resample target] — For our LSW the resample is a silent no-op: μ reads `sqrt(slip1²+slip2²)` (un-resampled directional slip), and there is no slip-magnitude accumulator to resample

**Category:** ASSUMPTION / DEVIATION (plan-completeness gap that defeats the primary target)

**Description:**
The corrected Phase-3 sketch and §4.2 say, for LSW: build directional `slip1/slip2` from the **un-resampled** friction-solve rate (correct, matches SeisSol), and separately "resample the slip-rate magnitude and integrate it into the accumulated slip **(which drives the slip-weakening μ)**" — sketched as `GatherFaceState(f, dSlipMag); Apply(R, out, dSlipMag); ScatterAddAccumSlip(f, out);`.

This assumes a SeisSol-style **scalar accumulated-slip-magnitude state** that (a) exists and (b) is what the LSW strength reads. **Neither holds in our code:**
- Our LSW slip-weakening coefficient is evaluated from `delta = sqrt(d.slip1² + d.slip2²)` — the magnitude of the **directional** accumulated slip — in *four* places: `friction_substep_iterator.cpp:171`, `fault_state_channel.hpp:37`, `tpv205_stations.hpp:204`, `tpv31_stations.hpp:233`.
- `slip1/slip2` are accumulated from the **un-resampled** rate (`friction_substep_iterator.cpp:190-191`, `:125-126`; `tpv102_substep_iterator.cpp:182`; `tpv104_substep_iterator.cpp:410`).
- `DOFData` has **no** `accumulated_slip_mag` / `slip_magnitude` field (`grep` in `fault_face_flux.hpp` = 0 hits).

So `ScatterAddAccumSlip` would write a resampled magnitude into a field that **nothing reads**, while μ keeps using `sqrt(slip1²+slip2²)` built from the un-resampled directional slip. **The resample has zero effect on the LSW strength → zero effect on the radiated solution for LSW.** TPV31 — the document's *headline* target (§Purpose, §7.3) — is LSW. As written, the plan would implement a fix that is a no-op for TPV31.

**Trigger:** `/code-implement` Phase 3, LSW branch, on TPV31/TPV205.

**Actual behavior:** resampled slip-magnitude is computed and discarded; μ unchanged; σ_n unchanged. The Phase-0/§7 LSW acceptance ("σ_n flat") cannot be met by the resample for LSW.

**Expected behavior:** to port SeisSol's LSW resample faithfully, add a resample-able scalar accumulator and make the LSW strength read it (keeping directional `slip1/slip2` un-resampled for output/slip totals, as SeisSol does).

**Suggested fix (plan edit — call out the data-model change explicitly):**
```diff
 **Resample target depends on the law** (Section 4.2): **LSW** $\to$ resample the
 slip-rate magnitude and integrate it into the accumulated slip (which drives the
-slip-weakening $\mu$); **rate-state** $\to$ resample the state-variable increment
-$\Delta\psi$.
+slip-weakening $\mu$); **rate-state** $\to$ resample the state-variable increment
+$\Delta\psi$.
+
+> **Data-model note (LSW).** Our LSW strength currently reads
+> $\delta=\sqrt{\text{slip1}^2+\text{slip2}^2}$ (the magnitude of the
+> *directional* accumulated slip) at four sites
+> (`friction_substep_iterator.cpp:171`, `fault_state_channel.hpp:37`,
+> `tpv205_stations.hpp:204`, `tpv31_stations.hpp:233`). SeisSol instead keeps a
+> *separate* `accumulatedSlipMagnitude` fed by the **resampled** $|V|$, distinct
+> from the (un-resampled) directional `slip1/slip2`. To match it, **add
+> `DOFData::accumulated_slip_mag`** (fed `R(|V|)\,\Delta t$ per sub-step) and
+> **repoint all four LSW-strength sites to read it** when `--fault-resample` is
+> on. Without this, the LSW resample is a no-op (μ still reads the un-resampled
+> directional-slip magnitude) and TPV31 will not change.
```
The Phase-3 sketch's `ScatterAddAccumSlip(f, out)` must target this new field, and the four strength sites must switch from `sqrt(slip1²+slip2²)` to `accumulated_slip_mag` under the resample flag (default off ⇒ both paths read the directional magnitude ⇒ byte-exact).
(Rate-state is fine as written: `DOFData::psi` is a single per-QP scalar the friction solve reads, so resampling Δψ maps cleanly.)

**Test case:**
```python
def test_R001_lsw_resample_actually_drives_mu():
    # With over-integration + resample ON on an LSW face, the slip measure that
    # drives mu must be the integral of the RESAMPLED |V|, NOT sqrt(slip1^2+slip2^2)
    # built from the un-resampled directional slip.
    f = run_lsw_face(overint=2, resample=True)
    delta_mu_driver = f.slip_measure_used_in_mu          # what LSWFrictionCoefficient sees
    delta_directional = math.hypot(f.slip1, f.slip2)     # un-resampled directional magnitude
    assert delta_mu_driver == approx(f.accumulated_slip_mag)   # reads the resampled accumulator
    assert delta_mu_driver != approx(delta_directional)        # and it actually differs
    # Off => both equal (byte-exact):
    g = run_lsw_face(overint=1, resample=False)
    assert g.slip_measure_used_in_mu == approx(math.hypot(g.slip1, g.slip2))
```

---

### [R-002] MODERATE [§4.2 / §6 Phase 2 test (iv) / Acceptance — "R=I at the minimal rule"] — "resample alone ⇒ byte-exact" only holds where the fault rule is unisolvent (#GP=#DOF); the target meshes are tetrahedral (triangle faces), where this fails at N≥3

**Category:** ASSUMPTION (geometry-/order-dependent claim stated as universal)

**Description:**
The corrected plan states `R` is the identity "when #GP = #DOF (the minimal mass-matrix rule, e.g. tensor Gauss on a quad face: (N+1)² GPs = (N+1)² DOFs)", and builds the regression sub-claim on it: Acceptance "`--fault-resample` alone (no over-integration) ⇒ also byte-exact, because R=I (Phase 2 test iv)", and Phase-2 test (iv) "at the minimal rule (#GP=#DOF), R=I."

`R=I` requires `V_N` to be **square and invertible** (unisolvent point set), i.e. #GP = #DOF. The example given is a **quad** face — but the actual target meshes are **tetrahedral**: `tpv102/configs/tpv102_spatial.toml:49 path = "tpv102/mesh/tpv102_1000m.msh"`, loaded at `spatial_dyn_driver.cpp:925` (`Mesh smesh(cfg.mesh.path,...)`), and gmsh `.msh` for these benchmarks are tet (per `miniapps/seas/CLAUDE.md`). So the fault faces are **triangles**, and the minimal degree-`2N` triangle rule is unisolvent for degree `N` only at low order:
- N=1 (deg-2 rule, 3 pts) vs deg-1 DOFs (3): equal ⇒ R=I.
- N=2 (deg-4 rule, 6 pts) vs deg-2 DOFs (6): equal ⇒ R=I.
- **N=3 (deg-6 rule, ~12 pts) vs deg-3 DOFs (10): #GP > #DOF ⇒ R ≠ I.**
- N≥3 generally: the degree-`2N` simplex rule is over-determined ⇒ `R ≠ I` ⇒ **`--fault-resample` alone is NOT byte-exact**, and the "resample is a no-op without over-integration" framing (§3.3/§4.2/§4.3) is only strictly true at p≤2 on triangles (and at all orders on quads).

The production configs include p3 (`aderO3`/p2/p3 jobs), so N≥3 on triangle faces is in scope.

**Trigger:** Phase-2 unit test (iv) at order ≥3 on a triangular reference face; regression check `--fault-resample` alone at p3 on a tet mesh.

**Actual behavior:** test (iv) `R==I` **fails** at p3 on a triangle face (R≠I); `--fault-resample`-alone is not byte-exact there — contradicting the stated contract, and likely read as a regression by the implementer.

**Expected behavior:** condition the claim on unisolvence. The hard regression contract is **both knobs off ⇒ byte-exact** (always true). "Resample alone ⇒ byte-exact / R=I" holds **only when #GP=#DOF** (quads at all N; triangles at N≤2).

**Suggested fix (plan edits):**
```diff
-- **$R$ is the identity when \#GP $=$ \#DOF** (the minimal mass-matrix rule, e.g.
-  tensor Gauss on a quad face: $(N+1)^2$ GPs $=$ $(N+1)^2$ DOFs). It does nothing
-  until the flux is over-integrated (Section 4.1) — the two techniques are
-  **coupled**.
+- **$R$ is the identity only when the fault rule is *unisolvent* for degree $N$,
+  i.e. \#GP $=$ \#DOF** (tensor Gauss on a quad face: $(N+1)^2=(N+1)^2$ at all $N$;
+  on a **triangle** face only at $N\le 2$ — at $N\ge 3$ the minimal degree-$2N$
+  simplex rule has \#GP $>$ \#DOF, so $R\ne I$). Where $R=I$ it does nothing until
+  the flux is over-integrated; the techniques are **coupled**. The target `.msh`
+  meshes are tetrahedral (triangle faces), so this matters at $p\ge 3$.
```
```diff
-  (iv) at the minimal rule (\#GP $=$ \#DOF), $R=I$ (so Phase 2 alone is
-  byte-exact).
+  (iv) at a **unisolvent** rule (\#GP $=$ \#DOF — quad faces at all $N$, triangle
+  faces at $N\le2$), $R=I$; at an over-determined rule (triangle, $N\ge3$),
+  $R\ne I$ and resample alone is *not* byte-exact.
```
```diff
-- Both knobs off $\Rightarrow$ byte-exact vs current (regression contract
-  preserved). `--fault-resample` alone (no over-integration) $\Rightarrow$ also
-  byte-exact, because $R=I$ (Phase 2 test iv).
+- Both knobs off $\Rightarrow$ byte-exact vs current (the regression contract —
+  always holds). `--fault-resample` alone $\Rightarrow$ byte-exact **only where the
+  fault rule is unisolvent** ($R=I$); on triangle faces at $p\ge3$ it is not, so
+  the regression gate must use **both-off**, not resample-alone.
```

**Test case:**
```python
def test_R002_minimal_triangle_rule_is_over_determined_at_p3():
    # Unisolvent only at low order on simplices; over-determined (R != I) at N>=3.
    import numpy as np
    for N, ngp, ndof, expect_I in [(1,3,3,True),(2,6,6,True),(3,12,10,False)]:
        # ngp = #points of the degree-2N MFEM triangle rule; ndof = (N+1)(N+2)/2
        assert ngp == tri_rule_npoints(degree=2*N)
        assert ndof == (N+1)*(N+2)//2
        R = build_degreeN_resample(geom="triangle", order=N, overint_factor=1)
        assert np.allclose(R, np.eye(ngp), atol=1e-12) == expect_I
```

---

### [R-003] LOW [§6 Phase 4 — "disable resample on across-fault material contrast"] — No detection criterion given for "genuine across-fault material contrast"

**Category:** QUALITY (under-specified guard; implementable but unstated)

**Description:**
Phase 4 says "Add a guard that **disables resample on any genuine across-fault material contrast** (matching SeisSol's `BiMaterialFault`)", but gives no criterion for detecting one. The implementer needs a concrete test. `DOFData` already carries the per-side impedances `Zp_plus/Zp_minus`, `Zs_plus/Zs_minus` (`fault_face_flux.hpp:29-30`), so the contrast is detectable per-DOF; and this correctly keeps resample **on** for TPV31 (symmetric material ⇒ `Zp_plus==Zp_minus` at each depth) while turning it off for a true jump.

**Suggested fix:**
```diff
-- Resample is **on** for TPV31/TPV102/104/205 (symmetric material). Add a guard
-  that **disables resample on any genuine across-fault material contrast**
-  (matching SeisSol's `BiMaterialFault`, Section 4.2).
+- Resample is **on** for TPV31/TPV102/104/205 (symmetric material). Add a guard
+  that **disables resample on any genuine across-fault material contrast**,
+  detected per-DOF as `Zp_plus != Zp_minus` (or `Zs_plus != Zs_minus`) beyond a
+  relative tolerance (`fault_face_flux.hpp:29-30`). For TPV31 these are equal at
+  each depth (material identical on both sides), so the guard keeps resample on;
+  it fires only for a true fault-normal jump (matching SeisSol's `BiMaterialFault`,
+  Section 4.2).
```

**Test case:** n/a (LOW; covered by the R-001/R-002 tests plus a trivial `guard(Zp_plus==Zp_minus)==enabled` assertion).

---

## Summary
- Critical issues: 1 (R-001 — LSW resample is a no-op against our slip-magnitude model; defeats TPV31)
- Moderate issues: 1 (R-002 — "R=I / resample-alone byte-exact" only holds for unisolvent rules; tet meshes break it at p≥3)
- Low issues: 1 (R-003 — Phase-4 contrast-guard criterion unspecified)
- Plan compliance: N/A (plan document). **Round-1 findings: all resolved.** The `R` formula is confirmed exact against SeisSol (idempotent + W-self-adjoint + rank = degree-N).
- Verdict: **PASS WITH FIXES.** The mechanism is now correct and source-faithful. The two remaining substantive items are *data-model/geometry* gaps the plan must call out before `/code-implement`: (R-001) the LSW resample needs a new `accumulated_slip_mag` field + repointed strength sites or it no-ops for TPV31; (R-002) the "resample-alone byte-exact" claim must be conditioned on unisolvence because the targets are tet meshes. Neither changes the physics; both change what the implementer must build/test.

## Unreviewed Areas
- Exact MFEM triangle-rule point counts per requested degree were reasoned from standard simplex-cubature sizes (3/6/12/16 for deg 2/4/6/8), not enumerated from `mfem::IntegrationRules` at runtime; R-002's `ngp` values should be confirmed by querying `IntRules.Get(Geometry::TRIANGLE, 2N).GetNPoints()` in the Phase-0 harness. The *direction* of the claim (over-determined at N≥3 on simplices) is robust regardless of the exact counts.
- §4.4 ("N Gauss–Legendre time nodes") is unchanged context, not re-audited (explicitly "context, not the lever"; not on the implementation path).
- The PDF was re-rendered and spot-checked (pp. 6–11) but not re-diffed line-by-line against the `.md`; after these edits, rebuild via the `parts_dev/_assets` toolchain (without `--number-sections`).
```
