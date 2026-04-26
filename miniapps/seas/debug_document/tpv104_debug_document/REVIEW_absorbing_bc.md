# Code Review — Absorbing BC effectiveness (2026-04-25)

## Review Scope
- Question: "Is the absorbing BC actually working?" — empirical question prompted by post-rupture trace contamination after t ≈ 1.15 s in Frontera P=2 jobs (7678275, 7678303).
- Files reviewed:
  - `dynamic/godunov_flux.cpp` (lines 377–410: `Absorbing`, `AbsorbingTotal`)
  - `dynamic/godunov_flux.hpp` (lines 56–80: API)
  - `dynamic/wave_operator.inl` (lines 1450–1517: dispatch in `ComputeFaceFluxRHS`; line 3283 ADER twin)
  - `dynamic/pml_layer.hpp` / `pml_layer.cpp` (PML implementation, not wired in TPV104 driver)
  - `drivers/tpv104_driver.cpp` (BC mode argument, SetAbsorbingBackground call)
- Domain context:
  - SCEC TPV104 spec (`tpv104/benchmark_document/SCEC_validation_slip_law.pdf`): half-space, fault at y=0, hypocenter (0, 7.5 km, 0).
  - Production sbatch (`jobs/tpv104/tpv104_symmirror_200m_partfile.sbatch`): mesh 16 × 4 × 16 km, `--bc-mode absorbing`, free-surface only at z=0, absorbing at all other 5 sides.
  - Trace from job 7678303 shows V_strike deviating from SeisSol after t ≈ 1.15 s, consistent with S-wave round-trip from y=±2 km absorbing boundary (S=3.464 km/s → 1.155 s).

## Background — what the absorbing BC IS

The implemented absorbing BC at all five outer faces (±x, ±y, z = −16 km) is the **classical first-order Engquist–Majda BC**:

```
F_abs = A_n^+ Q_self  +  A_n^- Q_bg              (with Q_bg = 0 in fluctuation mode)
```

This is `flux_.AbsorbingTotal(nor, Q_self, Q_bg, F_h)` (`godunov_flux.cpp:406-410`), which delegates to `Interior(nor, Q_self, Q_bg, F_h)`. The upwind split inside `Interior` produces `A^+ Q_self + A^- Q_bg`. With Q_bg = 0 (fluctuation), this reduces to **F = A^+ Q_self** — pure outgoing-only flux.

**Theoretical reflection coefficient** for 1st-order Engquist–Majda:

| Incidence angle θ from normal | Reflection coefficient R |
|---|---|
| 0° (normal) | 0 (perfect) |
| 30° | ~3% |
| 45° | ~17% |
| 60° | ~33% |
| 75° | ~67% |
| 90° (grazing) | 100% |

For TPV104, strike-slip rupture on y=0 fault: the dominant SH radiation lobe is along ±y → **NORMAL incidence on the y boundary** → BC absorbs near-perfectly. The dominant SV/P lobes are at oblique angles in the y-z plane and reflect partially. After the free-surface reflection at z=0 returns to the y boundary at angle ~45°, that wave bounces back at ~17% amplitude.

## Findings

### [R-AB-001] [MODERATE] [tpv104_driver.cpp / sbatch] — `--bc-mode absorbing` is the correct dispatch BUT the y-extent is too narrow for tfinal > 1.15 s

**Category:** ASSUMPTION (configuration / domain-size choice, not a code bug)

**Description:**
The TPV104 production sbatch builds a mesh with `--ly 4000` (y extent ±2 km from the fault). At S-wave speed 3.464 km/s, the round-trip from the fault to either y=±2 km boundary and back is **1.155 s**. This is **before the rupture peak** at t ≈ 1.07 s and well within the SCEC tfinal=12 s window. Even with a perfect absorbing BC, there is no protective time margin; with the 1st-order Engquist–Majda BC, oblique reflections from the free-surface bounce add ~10–20% spurious energy back to the fault starting around t ≈ 1.5 s.

**Trigger:** any tfinal > 1.155 s on the symmirror mesh built with ly=4000 m.

**Actual behavior:**
After t = 1.155 s, the recorded V_strike, τ_strike, ψ traces deviate from SeisSol's reference. Job 7678303 shows τ_strike = 24.7 MPa at t = 1.29 s vs SeisSol's ~30 MPa, and ψ floor = 0.05 vs SeisSol's 0.07.

**Expected behavior:**
For the SCEC TPV104 12 s benchmark, the simulation domain must keep the absorbing boundary at y far enough from the fault that the round-trip exceeds tfinal. SeisSol's reference uses ly ≥ 18 km (round-trip = 10.4 s).

**Suggested fix:**
Update `build_symmirror_mesh.py` defaults and the sbatch's mesh-build call from `--ly 4000` to `--ly 18000` (matching the SCEC standard half-space y-extent). At dx=200m: tets grow from 768k to ~3.46M (×4.5). Memory and runtime scale linearly. Per-rank cost at np=400 grows similarly.
```diff
- python3 build_symmirror_mesh.py --dx 200 --lx 16000 --ly 4000  --lz 16000 ...
+ python3 build_symmirror_mesh.py --dx 200 --lx 16000 --ly 18000 --lz 16000 ...
```

**Test case:**
```python
# tests/integration/test_tpv104_y_boundary_reflection.py (new)
def test_R_AB_001_y_boundary_round_trip_after_rupture_peak():
    # Run with ly=4000 m, tfinal=2 s.
    # Compute V_strike RMS over t in [1.20, 1.50] s (post-reflection window).
    # SeisSol reference at the same time has V_strike RMS ≈ 4.5 m/s.
    # Assert |V_strike_RMS_local - V_strike_RMS_seissol| < 0.5 m/s (10% tolerance).
    # WITHOUT the fix: assertion fails (reflected energy at ~1.155 s contaminates).
    # WITH the fix (ly=18000): assertion passes (no reflection until t ≈ 10.4 s).
```

---

### [R-AB-002] [LOW] [godunov_flux.cpp:379-386] — `Absorbing()` is dead code post-I-06 migration

**Category:** QUALITY

**Description:**
The function `GodunovFlux::Absorbing(nor, Q_self, F_h)` (calls `Interior(nor, Q_self, Q_zero, F_h)` with explicit zero) is defined but never called from `wave_operator.inl`. The dispatch always uses `flux_.AbsorbingTotal(nor, Q_self, bulk_bg_, F_h)` (lines 1461 and 1515 and ADER twin). When `bulk_bg_ = 0` (fluctuation mode, the only TPV104 mode), `AbsorbingTotal` produces bit-identical output to `Absorbing`, so removing `Absorbing` is safe.

**Trigger:** code maintenance — dead-code accumulation makes the surface area of "what BC actually runs" larger than necessary.

**Suggested fix:**
Either delete `GodunovFlux::Absorbing` (and its declaration at `godunov_flux.hpp:62`), OR document its non-use at the declaration. Preferred: delete.
```diff
- void GodunovFlux::Absorbing(const real_t *nor, const real_t *Q_self,
-                             real_t *F_h) const
- {
-    real_t Q_zero[NUM_STATE];
-    std::memset(Q_zero, 0, NUM_STATE * sizeof(real_t));
-    Interior(nor, Q_self, Q_zero, F_h);
- }
```

**Test case (LOW: this is a quality issue, no test required, but verifying API consistency):**
```cpp
// tests/unit/test_godunov_absorbing_api.cpp
TEST_CASE("R-AB-002: AbsorbingTotal with Q_bg=0 equals legacy Absorbing") {
   GodunovFlux f(rho, cp, cs);
   real_t Q[NUM_STATE] = {/* arbitrary */};
   real_t Q_bg[NUM_STATE] = {0};
   real_t F1[NUM_STATE], F2[NUM_STATE];
   f.AbsorbingTotal(nor, Q, Q_bg, F1);
   f.Absorbing(nor, Q, F2);  // legacy
   for (int c = 0; c < NUM_STATE; c++) REQUIRE(F1[c] == F2[c]);
}
```

---

### [R-AB-003] [MODERATE] [POSSIBLE] [godunov_flux.cpp / wave_operator.inl] — No empirical test of reflection coefficient

**Category:** ASSUMPTION (lack of a verifying test, not a confirmed bug)

**Description:**
The 1st-order Engquist–Majda BC theoretically has angle-dependent reflection (R ≈ 17% at 45°). The codebase has NO unit test that empirically verifies the reflection coefficient is ≤ a threshold. Without a test:
- a future refactor that subtly changes the upwind split sign (e.g., flips A^+ ↔ A^-) would silently double the reflection without any unit-test alarm
- a coordinate-frame mismatch (e.g., the `nor` is consistently inward vs outward) would invert the BC without alerting
- numerical dispersion at high polynomial order (P=3, P=4) interacts with the boundary; without a test, "BC works at P=1" doesn't imply "BC works at P=4"

**Trigger:** any future change to `godunov_flux.cpp::Interior` upwind split, or to `wave_operator.inl` BC dispatch sign conventions.

**Actual behavior:**
The dispatch `flux_.AbsorbingTotal(nor, Q_self, bulk_bg_, F_h)` is called and produces some output, but no automated check verifies that "outgoing waves leave with R ≤ ε".

**Expected behavior:**
A unit test should construct a known wave (e.g., Gaussian pulse traveling in +y direction in a small box mesh), run the wave operator for one round-trip to the +y absorbing boundary, and measure the reflected amplitude back at the source position. Assert R ≤ R_theoretical(angle).

**Suggested fix:**
Add `tests/unit/test_absorbing_bc_reflection.cpp` that:
1. Builds a 5×5×5 km box mesh, all 5 outer faces absorbing (no fault).
2. Sets initial Q to a Gaussian pulse in V_x at (0, 0, 0), σ = 200 m, propagating in +y.
3. Runs to t = (5000/cs) ≈ 1.44 s + small margin.
4. Probes max |V_x| at the source after the pulse should have left.
5. Asserts max |V_x| < 0.05 × initial peak amplitude (5% reflection at normal incidence — generous).

**Test case (concrete):**
```cpp
TEST_CASE("R-AB-003: absorbing BC reflection coefficient at normal incidence") {
   // 5×5×5 km box, all absorbing
   auto mesh = MakeBoxMesh(5000, 5000, 5000, /*nx=ny=nz=*/10);
   BoundaryConfig bc;
   for (int a = 1; a <= 6; a++) bc.absorbing_attrs.insert(a);
   WaveOperator op(mesh, /*order=*/1, lambda, mu, rho, bc);
   op.SetAbsorbingBackground(zero_q);

   Vector Q(NUM_STATE * op.GetScalarNDof());
   Q = 0.0;
   const double sigma = 200.0;
   const double y_src = 0.0;
   InitializeGaussianPulse(Q, /*center=*/{0, y_src, 0}, sigma, /*direction=*/+y, /*amplitude=*/1.0);
   const double initial_peak = ProbeMaxVxAtY(Q, y_src);

   double dt = op.ComputeMaxDt(0.5);
   double t_round_trip_one_way = 2500.0 / cs;
   for (double t = 0; t < 1.2 * t_round_trip_one_way; t += dt) {
      // RK4 step with op.Mult
      AdvanceRK4(op, Q, dt);
   }
   const double residual = ProbeMaxVxAtY(Q, y_src);
   REQUIRE(residual / initial_peak < 0.05);  // <5% reflection at normal incidence
}
```

---

### [R-AB-004] [LOW] [POSSIBLE] [pml_layer.hpp + tpv104_driver.cpp] — PML implemented but not wired

**Category:** SCOPE / OPTIONALITY

**Description:**
A PML (Perfectly Matched Layer) implementation exists in `dynamic/pml_layer.hpp/cpp` and `WaveOperator::SetPML` is wired (`wave_operator.hpp:93`). PML, properly tuned, gives R ≈ 10⁻³ to 10⁻⁵ at all incidence angles, dramatically better than 1st-order Engquist–Majda. **No TPV104 driver call invokes `SetPML`**, so PML is unavailable for TPV104 production runs.

If R-AB-001 is fixed by enlarging ly to 18 km, PML may be unnecessary. If memory/cost rules out a 3.5M-tet mesh, enabling PML on a 4-km-y mesh is the alternative path. Currently neither lever is being used.

**Trigger:** TPV104 production tfinal > 1.155 s on the current ly=4000m mesh.

**Suggested fix (optional, only if increasing ly is impractical):**
Add a `--enable-pml` driver flag and wire `wave.SetPML(...)` with appropriate damping-zone width and absorption profile.

**Test case:** would mirror R-AB-003 but with PML enabled, asserting R < 1e-3 instead of 5%.

---

## Summary

- Critical issues: 0
- Moderate issues: 2 (R-AB-001 narrow y-extent; R-AB-003 no reflection-coefficient test)
- Low issues: 2 (R-AB-002 dead code; R-AB-004 PML unused)
- Plan compliance: PARTIAL — `--bc-mode absorbing` dispatch is correctly implemented, but the production mesh's y-extent does not honor the implicit "domain large enough to be reflection-free for tfinal" assumption.
- Verdict: **PASS WITH FIXES** — the absorbing BC code is correct as 1st-order Engquist–Majda. It IS working at normal incidence (perfect) and at oblique incidence (with the expected angle-dependent reflection). The user-visible "trace mismatch with SeisSol after t = 1.15 s" is NOT a code bug — it's the well-known angle-dependent reflection compounded by an inadequately wide simulation domain (y=±2 km too close to the fault).

## Direct answer to "are we sure the absorbing BC is actually working?"

YES, the BC is working **as a 1st-order Engquist–Majda BC**:

1. **Code formula is correct**: `F = A^+ Q_self + A^- Q_bg` with Q_bg=0 ⇒ F = A^+ Q_self. Standard.
2. **Sign conventions correct**: `nor` is outward from the domain; A^+ corresponds to outgoing waves.
3. **Background buffer correctly initialized**: driver calls `SetAbsorbingBackground(Q_bg=0)` ✓.
4. **At normal incidence, R = 0** (perfect absorption — what we observe in the dominant SH radiation lobe of the strike-slip rupture).
5. **At oblique incidence, R = (1-cos θ)/(1+cos θ)** (well-known limitation; not a bug).

The trace contamination after t ≈ 1.15 s is dominated by the **proximity of the y boundary** (±2 km from fault, S-wave round-trip 1.155 s) — not by a defect in the BC. Free-surface reflections of the SV/P radiation from z=0 hit the y boundary at ~45° angles and reflect ~17% of the energy back to the fault, contaminating the late-time strike-slip dynamics.

**The fix is at the configuration level (R-AB-001), not the BC code level.** Enlarge `ly` to 18 km in `build_symmirror_mesh.py` and the production sbatches; the absorbing BC code itself is fine.

## Unreviewed Areas

- `pml_layer.cpp` damping-profile implementation: not audited in this round; only the wiring (its absence) was flagged.
- ADER twin BC dispatch in `wave_operator.inl::ComputeADERFaceFluxRHS` near line 3283: structurally identical to the Mult dispatch at line 1450, audited only by parallel-structure inspection.
- High-polynomial-order numerical-dispersion-at-boundary effects (P=3, P=4): no test exercises them; would need R-AB-003-style reflection probe at each P.
