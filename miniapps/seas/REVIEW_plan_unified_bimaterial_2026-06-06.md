# Code Review: PLAN_unified_bimaterial_volume_and_fault_2026-06-06 (plan review)

## Review Scope
- Plan: `document/bimaterial_dev/PLAN_unified_bimaterial_volume_and_fault_2026-06-06.md`
- Cross-checked against: `dynamic/spatial_setup.hpp`, `dynamic/fault_face_flux.{hpp,cpp}`,
  `dynamic/bimaterial_wave_operator.{hpp,inl}`, `dynamic/wave_operator.{hpp,inl}`,
  `dynamic/heterogeneous_material.hpp`, `drivers/spatial_dyn_driver.cpp`,
  `tpv6/benchmark_document/2007RuthRalphletter2.pdf` (TPV6/7 spec).
- Domain context: `CLAUDE.md`, `miniapps/seas/CLAUDE.md`, the drdg3d comparison doc,
  the prior contrast-guard review (`REVIEW_plan_bimaterial_contrast_guard_2026-06-06.md`).
- Nature: PLAN review. "Actual/Expected" = what the plan specifies vs what it must;
  "fix" = an edit to the plan / a hard requirement for the implementer.
- Note: PART A's findings from the prior review (R-001..R-007 there) are folded into the
  plan (verified: two-sided dissipation test in A4, fallback assert in A2, shared-
  convention blocking step in A2, collective-safety in A2). This review focuses on the
  NEW material: PART B (fault) and PART C (TPV6).

## Findings

### [R-001] CRITICAL — Part B1: per-side material must be evaluated at the fault DOF (±ε off-fault), NOT at element centroids; centroid read breaks TPV31 and injects a spurious bimaterial fault

**Category:** BUG / DEVIATION

**Description:**
B1 says to populate `Zp_plus/Zp_minus` from "the +/- adjacent element materials" and the
edge-case note says "the per-side impedance is that element's CENTROID material — the
same approximation the volume pool uses." This is wrong for any fault-symmetric but
spatially-varying material (TPV31 `depth_profile_1d`, SAFS CVM):

- `depth_profile_1d` depends on physical depth z. On the ASYMMETRIC production mesh the
  `+y` and `-y` fault-adjacent elements have DIFFERENT centroid depths, so their pooled
  (centroid-evaluated) materials differ → `Zp_plus != Zp_minus` even though the fault is
  physically NOT bimaterial. Near a layer step (z = -5000) the two centroids can fall in
  DIFFERENT layers → a LARGE spurious contrast.
- Result: (a) TPV31 is no longer byte-exact (the plan's own acceptance criterion in B1
  fails), and (b) the fault becomes spuriously bimaterial → a spurious `[[v_n]]`/sigma_n
  coupling — i.e. it would INTRODUCE the very kind of leak this whole effort fixes.

Conversely, "evaluate at the fault-QP IP" (the current single-side approach, `spatial_setup.hpp:123`)
is wrong for the across-fault `halfspace` material: a `Mode::Coefficient` keyed on
`sign(y)` evaluated at the fault (y=0) is AMBIGUOUS, so both sides would read the same
value → NO contrast for TPV6.

Neither "centroid per side" nor "QP per side" is correct for BOTH material models.

**Trigger:** TPV31 (depth profile, asymmetric mesh) with B1's centroid read → spurious
`Zp_plus != Zp_minus`. TPV6 (halfspace) with a QP-at-fault read → `Zp_plus == Zp_minus`
(no contrast).

**Actual behavior (as planned):** centroid-per-side (breaks TPV31) — with the edge note
even endorsing it.

**Expected behavior:** evaluate each side's material at the fault-DOF physical point
displaced a small ε INTO that side's element along the fault normal, mapped to that
element's reference frame. For depth-only materials ±ε in y leaves z unchanged → both
sides equal (byte-exact). For across-fault materials +ε→near side, −ε→far side → correct
contrast.

**Suggested fix (rewrite B1's evaluation rule):**
```diff
- the setup must accept a material lookup keyed by (fault face, side) ... the per-side
- impedance is that element's centroid material — the same approximation the volume pool
- uses; document it.
+ For each fault QP with physical position x_f and fault unit normal n (pointing to the
+ "+"/near side), evaluate:
+   side(+):  pt_p = x_f + eps*n ; find pt_p's reference IP in elem_plus ;
+             material.EvalAt(elem_plus,  T_plus,  ip_p, lam,mu,rho) -> Zp_plus, Zs_plus
+   side(-):  pt_m = x_f - eps*n ; find pt_m's reference IP in elem_minus ;
+             material.EvalAt(elem_minus, T_minus, ip_m, lam,mu,rho) -> Zp_minus, Zs_minus
+ with eps a small fraction of the local element size (e.g. 1e-3 * h_elem).  This gives
+ the correct ONE-SIDED material limit for BOTH depth-symmetric (Zp_plus==Zp_minus) and
+ across-fault (Zp_plus!=Zp_minus) materials.  Do NOT use element centroids (asymmetric-
+ mesh depth difference -> spurious contrast) and do NOT evaluate exactly at x_f (y=0 is
+ the discontinuity for across-fault materials).
```

**Test case (intent):**
```cpp
// TPV31 depth_profile on an ASYMMETRIC mesh, fault DOF near z=-5000:
//   build per-side impedance via the eps-offset rule
//   ASSERT Zp_plus == Zp_minus (rel 1e-12)  // depth-symmetric => no spurious contrast
// TPV6 halfspace, any fault DOF:
//   ASSERT Zp_plus == Zp(near)=6000*2670 and Zp_minus == Zp(far)=3750*2225 (or swapped)
```

---

### [R-002] MODERATE — Part B0: the analytic R/T verification is only valid for a LOCKED fault and must pin the quantity + sign convention

**Category:** ASSUMPTION

**Description:**
B0's acceptance test asserts a normal-incidence P-wave gives `R=(Z2-Z1)/(Z2+Z1)`,
`T=2Z2/(Z1+Z2)`. Those are the WELDED-interface coefficients. The fault is a frictional
interface: the coefficients hold only in the LOCKED regime (V=0, tau_corr=tau_trial),
where the fault Riemann reduces to the welded bimaterial Riemann. On a SLIPPING fault the
relation differs. Also, reflection/transmission differ for the VELOCITY vs the STRESS
field and carry opposite signs depending on convention; the plan pins neither, so the
implementer can produce a "failing" test that is actually a sign/quantity mismatch.

**Trigger:** running the B0 test on a slipping fault, or comparing the wrong field/sign.

**Expected behavior:** test in the locked limit with an explicit field + sign convention.

**Suggested fix:**
```diff
- a normal-incidence P-wave produces reflection/transmission matching
-   R=(Z2-Z1)/(Z2+Z1), T=2Z2/(Z1+Z2) to 1e-9
+ With the fault LOCKED (set friction strength above |tau_trial| so V=0, tau_corr=
+ tau_trial), a normal-incidence P-wave's imposed-state flux reproduces the welded
+ bimaterial Riemann: assert the interface NORMAL VELOCITY and NORMAL STRESS match the
+ analytic welded star state, with the STRESS reflection coeff R_sigma=(Z2-Z1)/(Z2+Z1)
+ and velocity transmission T_v=2Z2/(Z1+Z2) (state the sign convention: compression
+ positive, n pointing to side 2).  The slipping regime is validated separately against
+ SeisSol, not against welded R/T.
```

**Test case (intent):**
```cpp
// locked fault, Z1!=Z2, incident normal-velocity step from side1:
//   sigma_n_corr / incident matches the welded star traction; sign per convention.
```

---

### [R-003] MODERATE — Part C2: TPV6 nucleation is a UNIFORM SQUARE overstress patch, not a circular cosine taper

**Category:** DEVIATION

**Description:**
The TPV6 spec (PDF p.4) is a 3000x3000 m SQUARE patch with UNIFORM initial shear 81.6
MPa inside (a step: 81.6 in-patch, 70 outside), TPV3-style instantaneous overstress.
C2 says "reuse `instantaneous_overstress_circular` with a square-equivalent or add a
square patch kind." Reusing the circular kind applies a COSINE-TAPERED CIRCULAR patch,
which does not match the spec (wrong shape + wrong taper) and will shift rupture
initiation.

**Trigger:** using `instantaneous_overstress_circular` for TPV6.

**Expected behavior:** a uniform square overstress patch.

**Suggested fix:**
```diff
- `[nucleation]` instantaneous overstress in the 3x3 km patch (reuse
-  `instantaneous_overstress_circular` with a square-equivalent or add a square patch kind)
+ Add a nucleation kind `instantaneous_overstress_square`: uniform delta_tau over a
+ square patch [|x-x0|<=L/2 and |z-z0|<=L/2] (L=3000 m), NO taper, applied to the strike
+ component at t=0 so in-patch tau_strike = 81.6 MPa (= 70 + 11.6).  Matches TPV3/TPV6.
```

**Test case (intent):**
```cpp
// square patch: a DOF at patch center and a DOF just inside a corner both get the FULL
// 11.6 MPa overstress (no taper); a DOF just outside gets 0.
```

---

### [R-004] MODERATE — Part C2: TPV6 station output requires DISPLACEMENT, which the velocity-stress solver does not carry

**Category:** EDGE_CASE / ASSUMPTION

**Description:**
The spec (PDF p.2) requires TPV6/7 stations to report DISPLACEMENT and VELOCITY on each
side. The wave operator state Q is (stress, velocity) — there is no displacement field.
C2 says the station writer emits "displacement & velocity on each side" without
specifying that displacement must be obtained by TIME-INTEGRATING the on-fault velocity
per station (a new per-station accumulator), with the correct initial condition (0) and
output cadence.

**Trigger:** implementing `tpv6_stations.hpp` to read a displacement field that does not
exist.

**Expected behavior:** accumulate displacement = ∫ v dt per side per station.

**Suggested fix:**
```diff
- per-SIDE station writer (displacement + velocity on `+` and `-` sides separately ...)
+ per-SIDE station writer.  Velocity is read from Q directly per side.  DISPLACEMENT is
+ NOT a state variable: maintain a per-station, per-side displacement accumulator
+ d += v * dt updated every step (trapezoidal over the macro step), initialized to 0,
+ written at the station output cadence.  Document the time-integration scheme.
```

**Test case (intent):**
```cpp
// constant velocity v0 for time T at a station => displacement output == v0*T (rel 1e-6).
```

---

### [R-005] MODERATE — Part B2: there are 5 distinct guarded Evaluate variants; a single global flag is only safe if B0 maps EACH variant to a per-side-A conversion site

**Category:** ASSUMPTION

**Description:**
The homogeneity guard appears in FIVE different entry points — `Evaluate` (`:359`),
`EvaluateTotal` (`:499`), `EvaluateADER_LSW` (`:778`), `EvaluateLSW` (`:912`),
`EvaluateADER_LSW_ForcedRupture` (`:1030`). B2 relaxes all of them with one
`SetPerSideFluxApplied(true)` flag set by the matrix operator. That is correct ONLY if
EVERY one of those variants, on the matrix path, has its imposed-state flux converted
with per-side A. B0 currently says "confirm per-side A at every conversion site"
generically; it must instead map EACH of the 5 variants to its conversion site and
confirm per-side A there. If any variant (e.g. an ADER corrector path) converts with a
single A, the global flag silently produces wrong flux on that path.

**Trigger:** a matrix-path run that enters a variant whose conversion site uses a single
A, with the flag set.

**Expected behavior:** per-variant verification gates the relaxation.

**Suggested fix:**
```diff
  B0 task 1: Trace EVERY fault-flux -> bulk conversion and confirm it applies PER-SIDE A
+   Produce a TABLE: {Evaluate, EvaluateTotal, EvaluateADER_LSW, EvaluateLSW,
+   EvaluateADER_LSW_ForcedRupture} x {its conversion site file:line} x {per-side A? Y/N}
+   x {used by which driver path (spatial RK / spatial ADER / TPV* native)}.
  B2: relax ONLY the guards whose variant is (a) used by a matrix-path driver AND (b)
+   confirmed per-side A in the B0 table.  For any variant NOT confirmed, KEEP the abort
+   (do not relax via the global flag).  Prefer a per-variant check over a single flag if
+   the table is mixed.
```

**Test case (intent):**
```cpp
// for each matrix-path variant used by the spatial driver: Zp_plus!=Zp_minus runs and
// matches the locked-fault welded Riemann (R-002); the others still abort.
```

---

### [R-006] LOW — Part B3: `halfspace_across_fault` sign-at-plane is undefined; rely on the ±ε eval, and pin the side->Q_plus/Q_minus mapping

**Category:** EDGE_CASE

**Description:**
B3 selects side by "sign of (x-x0).normal", which is 0/ambiguous exactly on the plane.
With R-001's ±ε evaluation the material is never queried exactly at the plane, so this is
benign — but B3 must say so, and must pin which half-space (`side_near`/`side_far`) maps
to `Q_plus`/`Q_minus` (the fault BP5 `ref_normal=(0,-1,0)` convention), or TPV6's fast/
slow sides land on the wrong y-sides (plausible-but-wrong rupture asymmetry).

**Suggested fix:**
```diff
+ B3: side selection uses sign((x-x0).n); the value exactly on the plane is never used
+ (R-001 evaluates at +/-eps).  Define side_near := (x-x0).n > 0 maps to Q_plus and
+ side_far := < 0 maps to Q_minus, consistent with ref_normal=(0,-1,0); assert this in a
+ parse test against a 2-element straddling fixture and document in the TPV6 config header.
```

---

### [R-007] LOW — Part C2: free-surface x bi-material corner is listed but has no acceptance test

**Category:** EDGE_CASE

**Description:**
C2 correctly flags the free-surface/bi-material-fault corner (TPV6 reaches z=0) as a
risk, but provides no concrete check. Without one, a spurious normal-stress source at the
surface-fault line would be invisible until the full Frontera run.

**Suggested fix:**
```diff
+ C2 acceptance: add a smoke check that, with the fault LOCKED and a quiescent initial
+ state, the on-fault sigma_n at the shallowest DOFs (z->0) shows NO drift over N steps
+ (the bi-material coupling is zero with no slip), confirming the free-surface + per-side
+ A composition introduces no spurious normal stress at the corner.
```

---

## Summary
- Critical issues: 1 (R-001 — per-side material evaluation breaks TPV31 byte-exactness
  and injects a spurious bimaterial fault; the ±ε-offset rule is required)
- Moderate issues: 4 (R-002 locked-fault R/T verification spec; R-003 square vs circular
  nucleation; R-004 displacement output not a state variable; R-005 per-variant guard
  mapping)
- Low issues: 2 (R-006 halfspace sign/mapping; R-007 free-surface corner test)
- Plan compliance: N/A (plan review vs codebase facts)
- Verdict: PASS WITH FIXES — the architecture is sound and the two foundational claims
  are verified (bimaterial operator inherits the base `Mult` per-side fault conversion;
  the fault-flux math is per-side). But R-001 is a must-fix BEFORE B1 implementation: as
  written, B1 would break TPV31 and create a spurious bimaterial fault. Fix R-001..R-005
  in the plan; R-006/R-007 are quick hardening.

## Unreviewed Areas
- Whether `ExchangeBiMaterialNeighbours_` already covers SHARED FAULT faces (B1/A2 both
  depend on it) — not traced to ground truth; B1 flags it as "reuse/extend", which is the
  right disposition, but the implementer must confirm before the np>1 fault path.
- The exact ADER-corrector fault conversion site (one of the 5 variants) — B0's job;
  flagged in R-005.
- PART A internals — covered by the prior review; only its fold-in was checked here.
- Frontera numerics (TPV31 A/B, TPV6 vs SCEC) — not evaluable locally (no full-mesh
  runs); acceptance thresholds reasonable but unverified until the runs return. TPV7 is
  ill-posed by design (grid-dependent) — validate against ensemble spread, per the plan.
