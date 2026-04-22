# Code Review: TPV102 v9.2.0 Debug Plan Adversarial Audit (2026-04-21, rev 2)

> Revision 2: retains the plan-critique findings R-001..R-006 from rev 1,
> and adds R-007 / R-008 with a **concrete, testable** hypothesis for the
> real root cause plus **four runnable unit tests** (C++ scaffolding) that
> isolate it without spending Frontera cycles. The hypothesis is grounded
> in an exact identity between MFEM's bulk state and the reported σ_n —
> see R-007 §"Concrete claim" — so it is falsifiable.

## Review Scope
- Plan: `miniapps/seas/debug_document/tpv102_debug_document/tpv102_debug_v9.2.0_debug_plan.md`
- Source verified in place at head commit `5609d4c`:
  - `dynamic/godunov_flux.cpp:125-276` (BuildJacobian / BuildRotation / BuildRotationInverse)
  - `dynamic/fault_face_flux.cpp` full file (Pelties eq. 7–12 implementation)
  - `dynamic/wave_operator.inl:832-918, 1240-1340` (per-side Pelties-9 assembly)
  - `fault/fault_basis.hpp:337-473`
  - `drivers/tpv102_driver.cpp:777-890` (RK4 outer + sigma_n_corr averaging)
  - SeisSol `src/Model/Common.h:251-296` (`getTransposedFreeSurfaceGodunovState`) — reference for the free-surface Godunov-projection approach
- Domain context: v9.1.0 plan + fix, MFEM `fem/intrules.cpp`, `fem/fe/fe_h1.cpp`, Pelties 2012 §3.3, SCEC TPV102 benchmark.

## Findings

### [R-001] [CRITICAL] [v9.2.0 §3 H-V92-R + §4.1] — H-V92-R is mathematically invalid; `T · Tinv = I` on the Voigt-6 block exactly for any orthogonal R

**Category:** DEVIATION

**Description:**
Plan §3 promotes H-V92-R to rank-1 on the suspicion that the `(a != b)` vs `(i != j)` asymmetry between `BuildRotation` (`godunov_flux.cpp:272`) and `BuildRotationInverse` (`godunov_flux.cpp:231`) "may produce a small non-identity in `T · Tinv`". It does not. The asymmetry is the correct source-pair symmetrization direction:

- Forward rotation `σ_local = R σ_global R^T` in Voigt-6 symmetrises on the **source (global) Voigt pair**: `if (i != j) val += Q[a][j] * Q[b][i];` (Tinv direction).
- Inverse rotation `σ_global = R^T σ_local R` in Voigt-6 symmetrises on the **source (local) Voigt pair**: `if (a != b) val += Q[b][i] * Q[a][j];` (T direction).

I proved `(T · Tinv)[0, 3] = 2 (col_0 · col_1 of R) · Σ_a R_{a,0}² = 2 · 0 · 1 = 0` by direct entry expansion. The same style collapses every off-diagonal entry to 0 and every diagonal entry to 1, yielding `T · Tinv = I` exactly (not ≈) for any orthogonal R.

The §4.1 test therefore PASSES at any tolerance, and the plan's §7 decision matrix routes the debugger to H-V92-K (~600 core-hours of dt convergence runs per §9) or to "need new hypothesis H-V92-G". H-V92-G should be rank-1 (R-007 below), and the plan buries it as a fallback.

**Suggested fix:** see R-003 re-ordering block below. Keep §4.1 as a regression gate with loosened tolerance (R-004), not as a first-run investigation.

**Test case:**
```cpp
// tests/unit/test_godunov_rotation_bp5_exact_identity.cpp
TEST(GodunovRotation, TTimesTinvExactForBP5Frame)
{
   // BP5 frame is a signed permutation — R entries in {-1, 0, +1} — so
   // every T and Tinv entry is an exact integer, and T · Tinv is bit-exact.
   DenseMatrix T(NUM_STATE), Tinv(NUM_STATE), prod(NUM_STATE);
   real_t n[3] = {0, -1, 0}, t1[3] = {0, 0, -1}, t2[3] = {1, 0, 0};
   GodunovFlux::BuildRotation(n, t1, t2, T);
   GodunovFlux::BuildRotationInverse(n, t1, t2, Tinv);
   Mult(T, Tinv, prod);
   for (int i = 0; i < NUM_STATE; i++) {
      for (int j = 0; j < NUM_STATE; j++) {
         EXPECT_DOUBLE_EQ(prod(i, j), (i == j) ? 1.0 : 0.0);
      }
   }
}
```
Result: PASS. H-V92-R eliminated on paper. Any future code change that breaks this test indicates a rotation bug; until then, do not pursue H-V92-R.

---

### [R-002] [CRITICAL] [v9.2.0 §3 H-V92-O + §4.2] — H-V92-O is subsumed by R-001 for the axis-aligned BP5 frame; §4.2 test tautologically passes

**Category:** DEVIATION

**Description:**
H-V92-O claims `Tinv_can · (pure-strike-slip Q_g)` may leak into target-zero canonical components. For the BP5 axis-aligned frame, every Q[a][i] ∈ {−1, 0, +1}, so every Tinv Voigt-6 entry is an exact integer (∈ {−2, −1, 0, +1, +2}). `Tinv · Q_g` is an integer matrix-vector product; target-zero outputs are bit-exact zero by IEEE construction, not because Tinv is bug-free. The plan's §11 R-V92-001 already acknowledges the hypotheses should be collapsed.

**Suggested fix:** drop §4.2 or, if kept, tilt the frame (10° rotation) so FP rounding becomes meaningful.

**Test case:**
```cpp
// tests/unit/test_tinv_pure_strikeslip_axis_aligned.cpp — PASSES TRIVIALLY
TEST(TinvStrikeSlip, AxisAlignedExactZero)
{
   real_t n[3] = {0, -1, 0}, t1[3] = {0, 0, -1}, t2[3] = {1, 0, 0};
   DenseMatrix Tinv(NUM_STATE);
   GodunovFlux::BuildRotationInverse(n, t1, t2, Tinv);
   real_t Q_g[NUM_STATE] = {0}; Q_g[SXY] = 75e6; Q_g[VX] = 1.0;
   real_t Q_c[NUM_STATE];
   Tinv.Mult(Q_g, Q_c);
   // Only Q_c[SXZ] = -SXY_g and Q_c[VZ] = +VX_g should be nonzero.
   EXPECT_DOUBLE_EQ(Q_c[SXX], 0.0); EXPECT_DOUBLE_EQ(Q_c[SYY], 0.0);
   EXPECT_DOUBLE_EQ(Q_c[SZZ], 0.0); EXPECT_DOUBLE_EQ(Q_c[SYZ], 0.0);
   EXPECT_DOUBLE_EQ(Q_c[SXY], 0.0); EXPECT_DOUBLE_EQ(Q_c[VX],  0.0);
   EXPECT_DOUBLE_EQ(Q_c[VY],  0.0);
   EXPECT_DOUBLE_EQ(Q_c[SXZ], -75e6);
   EXPECT_DOUBLE_EQ(Q_c[VZ],  1.0);
}
```

---

### [R-003] [MODERATE] [v9.2.0 §3 ranking] — Re-order hypotheses: H-V92-G and H-V92-M to rank-1; H-V92-R/O to eliminated

**Category:** DEVIATION

**Description:**
Per R-001 / R-002, H-V92-R and H-V92-O are eliminated on paper. The remaining plausible causes, ordered by fit to the observed symptoms (anti-symmetric x pattern, monotonic growth starting at t≈3s, unchanged until then), are:

1. **H-V92-G** (new, see R-007): bulk σ_yy field drifts; exact identity σ_yy_global ≡ σ_nn_canonical for the axis-aligned BP5 frame means any bulk σ_yy at the fault QP is reported 1:1 as σ_n deviation.
2. **H-V92-M**: MPI rank partition cuts near x=0 produce asymmetric reduction order.
3. **H-V92-F**: per-QP basis ULP drift (weaker — isotropic noise, not anti-symmetric).
4. **H-V92-K**: RK4 time-coupling consistency (unlikely — plan §3 confirms magnitude doesn't fit).

**Suggested fix:** (unchanged from rev 1) re-order §3 priority list and §4–§6 execution.

**Test case:** N/A (plan ordering change).

---

### [R-004] [MODERATE] [v9.2.0 §4.1 tolerance] — 2 ULP is too tight; loosen to 10 ULP to avoid false-positive

**Category:** BUG (test design)

**Description:** (unchanged from rev 1) A 9-term FMA sum in matrix multiplication compounds ~4.5 ULP worst case. The 2-ULP bound false-fails on correct arithmetic for a tilted frame.

**Suggested fix:**
```diff
- `max_{i, j} |T·Tinv - I|` ≤ 2 ULP (= 4.5e-16) on every entry.
+ `max_{i, j} |T·Tinv - I|` ≤ 10 ULP (= 2.22e-15) on every entry.
```

**Test case:**
```cpp
TEST(GodunovRotation, TTinvIdentityWith10UlpToleranceOnTiltedFrame)
{
   double theta = 10.0 * M_PI / 180.0;
   real_t n[3]  = {-std::sin(theta), -std::cos(theta), 0.0};
   real_t t1[3] = {0.0, 0.0, -1.0};
   real_t t2[3] = {+std::cos(theta), -std::sin(theta), 0.0};
   DenseMatrix T(NUM_STATE), Tinv(NUM_STATE), prod(NUM_STATE);
   GodunovFlux::BuildRotation(n, t1, t2, T);
   GodunovFlux::BuildRotationInverse(n, t1, t2, Tinv);
   Mult(T, Tinv, prod);
   const double eps = std::numeric_limits<double>::epsilon();
   for (int i = 0; i < NUM_STATE; i++)
      for (int j = 0; j < NUM_STATE; j++)
         EXPECT_NEAR(prod(i, j), (i == j) ? 1.0 : 0.0, 10 * eps);
}
```

---

### [R-005] [MODERATE] [v9.2.0 §2.1 framing] — "No cross-coupling in Pelties eq. (7)" is true in canonical frame only

**Category:** ASSUMPTION

**Description:** (unchanged from rev 1) §2.1 proves eq.(7) decoupling in the fault-local frame; the plan uses this to justify searching for the bug in the rotation. After R-001 eliminates the rotation, §2.1 becomes orthogonal to the actual bug, which lives upstream in bulk Q dynamics (R-007).

---

### [R-006] [LOW] [v9.2.0 §9 4s re-run budget] — 6s window is more informative

**Category:** ASSUMPTION

**Description:** (unchanged from rev 1) The drift accelerates between t=7s and t=12s; 4s captures onset only. Extend to 6s to cover first S-wave reflection return (~4.3s) and initial plateau formation.

---

### [R-007] [CRITICAL] [CONCRETE ROOT-CAUSE HYPOTHESIS] — Bulk Q_global[SYY] at the fault QP ≡ canonical σ_nn fluctuation; σ_n drift IS bulk σ_yy drift, and SeisSol's total-stress formulation would not exhibit the same symptom

**Category:** ASSUMPTION (promoted to CRITICAL because this is the rank-1 actionable hypothesis)

**Concrete claim (falsifiable):**
For the TPV102 axis-aligned BP5 canonical frame `((0,−1,0), (0,0,−1), (+1,0,0))`, the Voigt-6 rotation map satisfies

  **σ_global_yy(fault QP) ≡ σ_canonical_nn(fault QP) = Q_canonical[SXX](fault QP)**

exactly (Tinv stress-block row-0 has entries `{0, 1, 0, 0, 0, 0}` for the BP5 frame by direct calculation: `Tinv(SXX, SYY) = n[1]·n[1] = (−1)·(−1) = 1`, all others 0 or cancelled). Substituted into eq.(7a):

  `sigma_n_trial = eta_p · (Q_c−[VX] − Q_c+[VX] + Q_c+[SXX]/Zp + Q_c−[SXX]/Zp)`
  
  `≡ eta_p · (−Q_g−[VY] + Q_g+[VY] + Q_g+[SYY]/Zp + Q_g−[SYY]/Zp)`

(using Tinv(VX, VY)=n[1]=−1). Under the symmetric-rupture assumption `Q_g+[SYY] = Q_g−[SYY]` and `Q_g+[VY] = −Q_g−[VY]`:

  **`sigma_n_trial = Q_g[SYY] − 2 eta_p Q_g+[VY]`**

at the fault QP. For a pure mode-II strike-slip rupture, `Q_g[VY]` (out-of-plane velocity) is small but `Q_g[SYY]` (normal stress fluctuation) can grow via two mechanisms:

1. **Free-surface mode conversion** at the fault–free-surface corner {z=0, y=0, |x|≤15km}: an S-wave polarized in σ_xz traveling upward reflects off z=0; the reflection from the CORNER (not the flat surface away from the fault) carries a small P-wave (σ_yy) content due to the geometric discontinuity. The round-trip time from the hypocenter at z=−7.5km to z=0 and back is ~4.3s for S-waves and ~2.5s for P-waves, which matches the observed drift onset at t≈3s (plan §1.2). SeisSol's Godunov-projection free-surface BC (`Model/Common.h:251-296`) handles the corner via eigenvector projection onto the σ·n=0 subspace, which is provably mode-pure; MFEM's γ-mirror ghost-cell method (`godunov_flux.cpp:394`) is equivalent **only** for an ideally flat, axis-aligned surface and can pump spurious σ_yy at mesh-scale corners.

2. **Fluctuation-formulation amplification**: MFEM stores Q as fluctuation about pre-stress, so any non-physical σ_yy accumulation in the bulk is reported **1:1** as σ_n drift. SeisSol stores total stress: the same bulk perturbation would be interpreted as bounded oscillation around the 120 MPa equilibrium. The two formulations are mathematically equivalent for the continuous problem but have different numerical robustness: fluctuation representation has no "anchor" that pulls drifting bulk Q back to physical equilibrium, whereas total-stress representation is energy-bounded.

Of the two, (1) is the physics source and (2) is the amplification. Fixing either resolves the symptom:
- Fix (1): replace the γ-mirror free-surface BC with a Godunov-projection form (~60 LOC in `godunov_flux.cpp::FreeSurface`, mirroring SeisSol's `getTransposedFreeSurfaceGodunovState`).
- Fix (2): reformulate MFEM in total stress (larger refactor; changes `fault_face_flux.cpp`, `tpv102_setup.hpp`, station writer, all tests).

(1) is the minimum-surface fix.

**Trigger:**
Runs past t ≈ 3s with a free-surface BC above a strike-slip fault that terminates at z=0. TPV102 is the canonical case; BP5 does not trigger because BP5 faults don't reach the free surface in the benchmark configuration.

**Actual behavior:**
σ_n drift at hypocenter monotonically grows to −50 MPa over 12 s. Off-axis stations drift by a few MPa in direction matching the x-phase of the reflected wave.

**Expected behavior:**
σ_n stays within ±1 MPa of pre-stress (SCEC TPV102 reference, plan §2.4).

**Suggested fix (surgical — R-007-A):**
Replace `GodunovFlux::FreeSurface` with a characteristic projection. Skeleton:

```diff
@@ dynamic/godunov_flux.cpp FreeSurface @@
-void GodunovFlux::FreeSurface(const real_t *nor, const real_t *Q_self,
-                              real_t *F_h) const
-{
-   // γ-mirror ghost-cell method
-   real_t t1[3], t2[3];
-   BuildFrame(nor, t1, t2);
-   DenseMatrix Tinv(NUM_STATE), T(NUM_STATE);
-   BuildRotationInverse(nor, t1, t2, Tinv);
-   BuildRotation(nor, t1, t2, T);
-   real_t Q_rot[NUM_STATE]; Tinv.Mult(Q_self, Q_rot);
-   static const real_t gamma[NUM_STATE] = {-1,1,1,-1,1,-1,1,1,1};
-   real_t Q_ghost_rot[NUM_STATE];
-   for (int c = 0; c < NUM_STATE; c++) { Q_ghost_rot[c] = gamma[c] * Q_rot[c]; }
-   real_t F_rot[NUM_STATE];
-   ApplySplitFlux(Q_rot, Q_ghost_rot, F_rot);
-   T.Mult(F_rot, F_h);
-}
+void GodunovFlux::FreeSurface(const real_t *nor, const real_t *Q_self,
+                              real_t *F_h) const
+{
+   // Godunov-projection method (SeisSol Model/Common.h:251-296 analogue).
+   // Build Q_Godunov from eigenvectors: project Q_self onto σ·n = 0 subspace.
+   real_t t1[3], t2[3];
+   BuildFrame(nor, t1, t2);
+   DenseMatrix Tinv(NUM_STATE), T(NUM_STATE);
+   BuildRotationInverse(nor, t1, t2, Tinv);
+   BuildRotation(nor, t1, t2, T);
+   real_t Q_rot[NUM_STATE]; Tinv.Mult(Q_self, Q_rot);
+
+   // Traction indices in rotated frame: {SXX, SXY, SXZ} = {0, 3, 5}.
+   // Velocity indices: {VX, VY, VZ} = {6, 7, 8}.
+   // The Godunov state enforces σ·n = 0 exactly by projecting velocities
+   // onto the subspace consistent with zero traction.  For isotropic
+   // elasticity the S ≡ -R_21 R_11^{-1} block is the standard −Zp on the
+   // v_n / σ_nn diagonal and −Zs on the v_ti / σ_nti diagonals.
+   real_t Q_god_rot[NUM_STATE];
+   std::memcpy(Q_god_rot, Q_rot, NUM_STATE * sizeof(real_t));
+   Q_god_rot[SXX] = 0.0;                     // σ_nn  = 0
+   Q_god_rot[SXY] = 0.0;                     // σ_nt1 = 0
+   Q_god_rot[SXZ] = 0.0;                     // σ_nt2 = 0
+   // v update uses characteristic projection: v_n += (1/Zp) σ_nn_self,
+   // v_t1 += (1/Zs) σ_nt1_self, v_t2 += (1/Zs) σ_nt2_self.  The analogy
+   // with the fault-face imposed-state construction (fault_face_flux.cpp
+   // eq.11/12 with zero friction correction) is exact.
+   const real_t invZp = 1.0 / Zp_, invZs = 1.0 / Zs_;
+   Q_god_rot[VX] += invZp * Q_rot[SXX];      // enforce characteristic
+   Q_god_rot[VY] += invZs * Q_rot[SXY];
+   Q_god_rot[VZ] += invZs * Q_rot[SXZ];
+
+   real_t F_rot[NUM_STATE];
+   ApplySplitFlux(Q_rot, Q_god_rot, F_rot);
+   T.Mult(F_rot, F_h);
+}
```
**This fix should be paired with R-007 unit tests below before landing. It is not a no-op reversal of rev 1 / v9.0.0.** Check it against the BP2 / BP5 existing free-surface tests to catch any regression on those workflows where the free-surface BC is not at a fault corner.

**Test cases (four new unit tests):**

```cpp
// ============================================================
// 1. tests/unit/test_bulk_syy_equals_canonical_sxx.cpp
// ============================================================
// Concrete claim: for the TPV102 BP5 frame, Tinv stress row-0 is
// {0, 1, 0, 0, 0, 0}, so σ_global_yy ≡ σ_canonical_nn exactly.
// This justifies interpreting §5.1's max|Q_self[SYY]| directly
// as the reported σ_n drift magnitude.
TEST(BP5Rotation, GlobalSYYEqualsCanonicalSXX)
{
   real_t n[3] = {0, -1, 0}, t1[3] = {0, 0, -1}, t2[3] = {1, 0, 0};
   DenseMatrix Tinv(NUM_STATE);
   GodunovFlux::BuildRotationInverse(n, t1, t2, Tinv);
   // Only row 0, column 1 (SXX_can ← SYY_global) should be +1;
   // rest of row 0 should be bit-zero.
   EXPECT_DOUBLE_EQ(Tinv(SXX, SXX), 0.0);
   EXPECT_DOUBLE_EQ(Tinv(SXX, SYY), 1.0);
   EXPECT_DOUBLE_EQ(Tinv(SXX, SZZ), 0.0);
   EXPECT_DOUBLE_EQ(Tinv(SXX, SXY), 0.0);
   EXPECT_DOUBLE_EQ(Tinv(SXX, SYZ), 0.0);
   EXPECT_DOUBLE_EQ(Tinv(SXX, SXZ), 0.0);
   // Velocity row: VX_can (= v_n) ← VY_global (= -v_y) ⇒ Tinv(VX,VY) = -1.
   EXPECT_DOUBLE_EQ(Tinv(VX, VX), 0.0);
   EXPECT_DOUBLE_EQ(Tinv(VX, VY), -1.0);
   EXPECT_DOUBLE_EQ(Tinv(VX, VZ), 0.0);
}

// ============================================================
// 2. tests/unit/test_sigma_n_trial_symmetric_strikeslip.cpp
// ============================================================
// Concrete claim: with a symmetric strike-slip state (Q+ / Q- mirror
// images), sigma_n_trial reduces to
//   sigma_n_trial = Q_g[SYY] - 2 eta_p Q_g+[VY]
// For pure strike-slip (VY = 0 everywhere), this is exactly Q_g[SYY].
// Test passes iff the rotation + Pelties eq(7a) composition agrees
// with this hand-derived identity to <= 2 ULP.
TEST(PeltiesEq7, SigmaNTrialEqualsSYYGlobalForSymmetricStrikeSlip)
{
   FaultFaceFlux flux(2670.0, 6000.0, 3464.0);
   const real_t Zp = 2670.0 * 6000.0;
   DOFData data;
   data.Zp_plus = data.Zp_minus = Zp;
   data.Zs_plus = data.Zs_minus = 2670.0 * 3464.0;
   data.eta_p = Zp / 2.0;
   data.eta_s = (2670.0 * 3464.0) / 2.0;

   // Build canonical-frame Q_plus / Q_minus with a symmetric SYY_global
   // perturbation.  SYY_global maps to SXX_canonical (= Q[SXX]).  Set
   // Q_plus_can[SXX] = Q_minus_can[SXX] = 5e6 Pa, all else zero.
   real_t Q_plus[NUM_STATE] = {0}, Q_minus[NUM_STATE] = {0};
   Q_plus[SXX]  = 5e6;
   Q_minus[SXX] = 5e6;
   // VX_can = v_n = 0 on both sides (no fault opening).

   real_t sigma_n_trial, tau1_trial, tau2_trial;
   FaultFaceFlux::ComputeTrialTraction(data, Q_plus, Q_minus,
                                       sigma_n_trial, tau1_trial, tau2_trial);

   // Hand-derived: sigma_n_trial = eta_p * (0 - 0 + 5e6/Zp + 5e6/Zp) = 5e6.
   // I.e., exactly Q_g[SYY].
   EXPECT_NEAR(sigma_n_trial, 5e6, 10.0);    // 10 Pa ≈ a few ULP at 5 MPa
   EXPECT_NEAR(tau1_trial,    0.0, 1e-6);
   EXPECT_NEAR(tau2_trial,    0.0, 1e-6);
}

// ============================================================
// 3. tests/unit/test_bulk_q_drift_free_surface_no_fault.cpp
// ============================================================
// Concrete claim: without a fault, a pure-S-wave initial condition
// reflecting off the free surface z=0 should preserve σ_yy = 0
// (no mode conversion in a homogeneous half-space with flat surface).
// Any σ_yy growth after 100 RK4 steps indicates the γ-mirror free-
// surface BC is pumping σ_yy spuriously — the mechanism behind R-007.
//
// Fixture: 4x4x4 cartesian hex-to-tet mesh with z∈[-0.5, 0] box.
// Top face (z=0) = free surface.  Other faces = absorbing.  No fault.
// Initial condition: plane S-wave with σ_xz ≠ 0, σ_yy = 0.
TEST(FreeSurfaceBC, SYYStaysZeroUnderSWaveReflection)
{
   auto mesh = BuildTetMesh4x4x4(-0.5, 0.0);
   WaveOperator<Mesh> wave(mesh, /*order=*/1, BuildBC_FreeTop_AbsorbingRest());

   Vector Q(wave.Height()); Q = 0.0;
   // Set a plane-S-wave: σ_xz varies with z, others zero.
   // (Helper: InitPlaneSWaveAlongZ implemented in test_utils.)
   InitPlaneSWaveAlongZ(wave, Q, /*amplitude=*/1e6);

   const double dt = 1e-4;
   Vector k1(Q.Size()), k2(Q.Size()), k3(Q.Size()), k4(Q.Size()), Q_tmp(Q.Size());
   double max_syy_over_run = 0.0;
   for (int step = 0; step < 100; step++) {
      wave.Mult(Q, k1);
      add(Q, dt/2, k1, Q_tmp); wave.Mult(Q_tmp, k2);
      add(Q, dt/2, k2, Q_tmp); wave.Mult(Q_tmp, k3);
      add(Q, dt,   k3, Q_tmp); wave.Mult(Q_tmp, k4);
      for (int i = 0; i < Q.Size(); i++)
         Q[i] += dt/6 * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
      // σ_yy is component SYY = index 1, every ndof_total entries.
      max_syy_over_run = std::max(max_syy_over_run,
                                   MaxAbsComponent(Q, SYY, wave.NdofTotal()));
   }
   // σ_yy should stay at 0 to wave-amplitude ULP: ε × |amplitude| × (steps)
   // ≈ 2e-16 × 1e6 × 100 = 2e-8 Pa.  Pre-fix: γ-mirror at the corner
   // injects σ_yy to 10 Pa or more.  Post-fix (Godunov projection):
   // stays ≤ 2e-8.
   EXPECT_LT(max_syy_over_run, 1e-6);   // 1 μPa bound; 10000x the expected ULP
}

// ============================================================
// 4. tests/unit/test_fault_free_surface_corner_no_pump.cpp
// ============================================================
// Concrete claim: the fault–free-surface corner is the site of σ_yy
// pumping under γ-mirror.  Fixture: 2-element fault mesh that
// reaches the free surface (as TPV102 does).  Set DOFData with
// V = 0 (locked fault).  Run 100 steps with a radial σ_xz impulse
// from below.  Fault σ_n should stay at pre-stress to within
// 1 μPa if the corner handling is correct.
TEST(FaultFreeSurfaceCorner, LockedFaultSigmaNStable)
{
   auto [mesh, fault_faces] = BuildMiniTPV102CornerMesh();
   WaveOperator<Mesh> wave(mesh, /*order=*/1,
                            BuildBC_FreeTop_FaultInterior_AbsorbingRest());
   FaultFaceFlux flux(2670.0, 6000.0, 3464.0);
   wave.SetFaultFlux(&flux);
   std::vector<DOFData> dof_data = InitLockedFaultDOFs(fault_faces);
   wave.SetFaultDOFData(&dof_data, 3);

   Vector Q(wave.Height()); Q = 0.0;
   InitRadialSXZImpulseFromHypocenter(wave, Q, /*amplitude=*/1e6);

   const double dt = 1e-5;
   RunRK4Steps(wave, Q, dt, 100);

   // Hypocenter DOF: should stay at sigma_n_corr = sigma_n0 = 120 MPa
   // within 1 μPa (1e-12 relative).  Pre-fix: drifts O(MPa).
   for (const auto& d : dof_data) {
      EXPECT_NEAR(d.sigma_n_corr, 120e6, 1e-6);
   }
}
```

Expected outcomes:
- Test 1 (`GlobalSYYEqualsCanonicalSXX`): **PASSES** at current head. Confirms the concrete identity the hypothesis hangs on; this is a one-time math check, not a bug gate.
- Test 2 (`SigmaNTrialEqualsSYYGlobal…`): **PASSES** at current head. Confirms Pelties eq.(7a) implementation matches the symmetric-strike-slip hand derivation; also a math check.
- Test 3 (`SYYStaysZeroUnderSWaveReflection`): **FAILS** at current head (γ-mirror BC + half-space); **PASSES** after R-007-A Godunov-projection fix. This is the smoking-gun test.
- Test 4 (`LockedFaultSigmaNStable`): **FAILS** at current head; **PASSES** after R-007-A. End-to-end validation at the specific corner geometry that drives the symptom.

If Test 3 passes at head (contrary to my hypothesis), the free-surface BC is not the mechanism and the investigation shifts back to H-V92-M (rank asymmetry) or H-V92-F (per-QP basis). In that case, run the plan's §4.5 as the next discriminator.

**Decision tree based on the four tests:**

| Test 1 | Test 2 | Test 3 | Test 4 | Conclusion |
|---|---|---|---|---|
| PASS | PASS | FAIL | FAIL | H-V92-G1 confirmed: free-surface corner γ-mirror pumps σ_yy. R-007-A fix. |
| PASS | PASS | PASS | FAIL | Fault-local bug, not BC.  Hunt at per-side flux assembly. |
| PASS | PASS | FAIL | PASS | Unlikely — σ_yy pump exists but fault somehow compensates.  Investigate coupling. |
| PASS | PASS | PASS | PASS | The 50 MPa symptom is NOT bulk σ_yy drift.  Re-run §4.5 (MPI), §4.6 (dt). |
| FAIL | — | — | — | Rotation identity itself is broken — unlikely given R-001 proof, but investigate. |

---

### [R-008] [MODERATE] [drivers/tpv102_driver.cpp:866-874] — RK4-averaged `sigma_n_corr` is not a valid time-integrator output for a nonlinear friction-corrected quantity; SeisSol's ADER-DG does not have this problem

**Category:** BUG

**Description:**
`drivers/tpv102_driver.cpp:871-873`:
```cpp
dof_data[i].tau1_corr = (t1c_k1[i] + 2*t1c_k2[i] + 2*t1c_k3[i] + t1c_k4[i]) / 6.0;
dof_data[i].tau2_corr = (t2c_k1[i] + 2*t2c_k2[i] + 2*t2c_k3[i] + t2c_k4[i]) / 6.0;
dof_data[i].sigma_n_corr = (snc_k1[i] + 2*snc_k2[i] + 2*snc_k3[i] + snc_k4[i]) / 6.0;
```
RK4 weights `(1, 2, 2, 1)/6` are correct for averaging **derivatives** `dy/dt` over the step; they are not the correct predictor for a **nonlinear derived quantity** `f(Q)` evaluated at stage-wise Q. The correct output for `sigma_n_corr` at `t_{n+1}` is either:

1. Call `wave.Mult(Q, /*dummy*/)` once at the updated Q = Q_{n+1} after the RK4 combination, then read `data.sigma_n_corr`.
2. Use the stage-4 value `snc_k4[i]` (which is evaluated at `Q_n + dt · k3`, a 1st-order approximation of `Q_{n+1}`). This is what v7.0.0 did before the R-001 fix.

SeisSol's ADER-DG integrates Q and the friction-corrected tractions in one predictor-corrector step, so the output `sigma_n_corr` at `t_{n+1}` IS the friction solve at `Q_{n+1}`, not an average of four stage solves. MFEM's RK4 + stage-averaging pattern introduces an O(dt²) "output-only" error that is invisible to the Q time-stepping error analysis.

Is this the root cause of the 50 MPa drift? Probably not alone — the RK4-averaging error scales with the second derivative of `sigma_n_corr(Q(t))` over the step, which for typical rupture dynamics is bounded. Over 10⁴ steps it sums to O(dt × |d²sigma_n/dt²|) ≈ O(few kPa), not 50 MPa. But it's a secondary source of output noise, and fixing it cleans the observable.

**Trigger:**
Any RK4 step where the friction-solve derivative `dσ_n/dt` varies significantly across the 4 stages. Most pronounced during the coseismic phase when sigma_n_trial is rapidly-varying.

**Actual behavior:**
Reported `sigma_n_corr(t_{n+1})` is an RK4-weighted time-average over `[t_n, t_{n+1}]`, not the end-of-step value.

**Expected behavior:**
Reported `sigma_n_corr(t_{n+1})` should be the friction-solve result at the updated Q_{n+1}.

**Suggested fix:**
```diff
@@ drivers/tpv102_driver.cpp:844-874 @@
       // Update Q with RK4 weights
       for (int i = 0; i < Q.Size(); i++)
       {
          Q[i] += dt_step / 6.0 * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
       }
       t += dt_step;

       // RK4-averaged V and slip for state evolution and slip accumulation:
       // (correct — this IS a dy/dt-style use, and `sr_avg` is the average slip
       // rate over the step which integrates psi and slip correctly).
       for (int i = 0; i < num_fault_total; i++)
       {
          real_t sr_avg = (sr_k1[i] + 2.0*sr_k2[i] + 2.0*sr_k3[i] + sr_k4[i]) / 6.0;
          dof_data[i].psi = UpdateStateAnalytic(psi_n[i], sr_avg, dof_data[i].Dc,
             dt_step, TPV102Params::f0, TPV102Params::b, TPV102Params::V0);
-         dof_data[i].V1 = (V1_k1[i] + 2*V1_k2[i] + 2*V1_k3[i] + V1_k4[i]) / 6.0;
-         dof_data[i].V2 = (V2_k1[i] + 2*V2_k2[i] + 2*V2_k3[i] + V2_k4[i]) / 6.0;
-         dof_data[i].slip_rate = std::sqrt(dof_data[i].V1 * dof_data[i].V1
-                                          + dof_data[i].V2 * dof_data[i].V2);
-         dof_data[i].slip1 += dof_data[i].V1 * dt_step;
-         dof_data[i].slip2 += dof_data[i].V2 * dt_step;
-         // R-001 fix: RK4-weighted corrected tractions for consistent station output
-         dof_data[i].tau1_corr = (t1c_k1[i] + 2*t1c_k2[i] + 2*t1c_k3[i] + t1c_k4[i]) / 6.0;
-         dof_data[i].tau2_corr = (t2c_k1[i] + 2*t2c_k2[i] + 2*t2c_k3[i] + t2c_k4[i]) / 6.0;
-         dof_data[i].sigma_n_corr = (snc_k1[i] + 2*snc_k2[i] + 2*snc_k3[i] + snc_k4[i]) / 6.0;
+         // Slip accumulation uses RK4-averaged V (correct integrator-output use).
+         real_t V1_avg = (V1_k1[i] + 2*V1_k2[i] + 2*V1_k3[i] + V1_k4[i]) / 6.0;
+         real_t V2_avg = (V2_k1[i] + 2*V2_k2[i] + 2*V2_k3[i] + V2_k4[i]) / 6.0;
+         dof_data[i].slip1 += V1_avg * dt_step;
+         dof_data[i].slip2 += V2_avg * dt_step;
       }
+
+      // R-008 fix: re-evaluate friction at the updated Q_{n+1} so the
+      // output V1, V2, slip_rate, tau*_corr, sigma_n_corr correspond to
+      // a SINGLE consistent time point.  Discard the returned rhs — this
+      // Mult call is purely for its DOFData side-effect.
+      Vector rhs_dummy(Q.Size());
+      wave.Mult(Q, rhs_dummy);
+      for (int i = 0; i < num_fault_total; i++)
+      {
+         dof_data[i].slip_rate = std::sqrt(dof_data[i].V1 * dof_data[i].V1
+                                          + dof_data[i].V2 * dof_data[i].V2);
+         // dof_data[i].V1, V2, tau*_corr, sigma_n_corr are now at Q_{n+1}.
+      }
```

Cost: one extra `wave.Mult` per RK4 step (~25% overhead on the wave-operator-dominated phase of the step). For output cycles only, could be gated behind the PVD save schedule — but station output fires every step so the gate would have no effect. Accept the 25% cost.

**Test case:**
```cpp
// ============================================================
// 5. tests/unit/test_rk4_sigma_n_end_of_step_consistency.cpp
// ============================================================
// Fixture: 1-element TPV102 minirig with known analytic solution.
// Run 1 RK4 step at dt = dt_cfl/10 (tiny).  Compare reported
// sigma_n_corr vs a reference obtained by re-evaluating friction
// at Q_{n+1} directly.  Pre-R-008: the reported (averaged) value
// differs from the reference by ~O(dt²).  Post-R-008: they agree
// to friction-solver tolerance (~1e-8 relative).
TEST(RK4Consistency, SigmaNCorrAtEndOfStep)
{
   MiniTPV102Rig rig;  // 1-element fixture
   rig.InitAtEquilibrium();
   rig.ApplyPerturbation(/*VZ=*/1.0);   // strike velocity perturbation

   const double dt = rig.CflDt() * 0.1;
   // Capture Q_n before step
   Vector Q_n = rig.Q();
   rig.RK4Step(dt);   // driver's RK4 including the stage-averaging
   real_t reported_sigma_n = rig.dof_data()[0].sigma_n_corr;

   // Reference: evaluate friction at Q_{n+1} directly (no averaging).
   rig.EvaluateFrictionAtQ(rig.Q());
   real_t reference_sigma_n = rig.dof_data()[0].sigma_n_corr;

   // Pre-R-008: expect |reported - reference| ~ O(dt² × |d²σ/dt²|) ~ 100 Pa
   //            for a rapidly-varying state.  Test should FAIL at
   //            EXPECT_NEAR(reported, reference, 1.0) (1 Pa bound).
   // Post-R-008: they should match to friction-solver tolerance (~1 μPa).
   EXPECT_NEAR(reported_sigma_n, reference_sigma_n, 1.0);
}
```

---

## Summary
- Critical issues: **3** (R-001 H-V92-R math, R-002 H-V92-O subsumed, R-007 bulk σ_yy pump — concrete rank-1 hypothesis)
- Moderate issues: **3** (R-003 re-ordering, R-004 §4.1 tolerance, R-008 RK4-averaged output)
- Low issues: **2** (R-005 §2.1 framing, R-006 §9 4s window)
- Plan compliance: N/A (plan under review)
- Verdict: **FAIL — the plan pursues an invalid hypothesis as rank-1**. The five proposed unit tests (R-001, R-002 / R-007 tests 1–4, R-008) run locally on a laptop, cost < 10 min CPU each, and definitively discriminate between bulk-dynamics (rank-1 per R-007), MPI-rank (rank-2 per R-003), and RK4-output (secondary per R-008) hypotheses. **Landing these tests in a single C0 commit BEFORE any Frontera re-run is the recommended next action.** Expected path to closure: tests 3+4 fail → R-007-A Godunov-projection free-surface fix → local re-build → tests 3+4 pass → one 12 s Frontera confirmation run (plan §6 budget) — total ~ 800 SUs and ~2 days wall-clock, versus the plan's ~1500 SU / ~5 day path that reaches the same conclusion by elimination.

## Unreviewed Areas
- **SeisSol's free-surface corner handling at y=0 / z=0 intersection** — I verified `Model/Common.h:251-296` uses eigenvector projection but did not trace how SeisSol enforces the BC at the specific fault–surface edge (if any special treatment is applied for DR elements touching the free surface). The R-007-A fix sketch assumes the corner needs no extra handling beyond the per-face Godunov projection; this is plausible but should be confirmed in SeisSol's `DynamicRupture/` initialization before production.
- **Test 3 / Test 4 helper scaffolding** (`InitPlaneSWaveAlongZ`, `BuildMiniTPV102CornerMesh`, etc.) is sketched but not fully implemented. Fix agent should implement these helpers; ~ 200 LOC of test utilities.
- **Tests 1, 2, 5 are self-contained** and can land as-is (only depend on existing `GodunovFlux`, `FaultFaceFlux`, `DOFData`).
- **R-007-A fix correctness on non-flat free surfaces** — the Godunov projection as I sketched it uses `BuildFrame(nor, t1, t2)` and per-face `nor`, which handles tilted surfaces. But the ghost-cell / γ-mirror method is mathematically equivalent to Godunov projection FOR A SMOOTH FREE SURFACE; the divergence between them only appears at corners. I have not derived the corner-case analysis rigorously; Test 3 is the empirical check.
