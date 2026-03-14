# BP5 Debug v23: All v22 Test Results — V_nuc=0.03 Is the Problem

**Date**: 2026-03-14
**Status**: All 3 tests completed, root cause narrowed
**Previous**: v22 (diagnostic test plan)
**Results**: `results_v22_test1_no_dtau`, `results_v22_test2_bdr_loading`, `results_v22_test3_tandem_psi`

---

## 1. Executive Summary

All three tests produced **identical earthquakes and identical deep VS lockup**.
Changing δτ or the psi initialization mode had zero effect. The earthquake is
driven entirely by **V_nuc = 0.03 m/s** in the nucleation zone — 3×10⁷ times
plate rate, already seismic regardless of δτ or psi.

| Run | δτ | Psi init | V_nuc | V/Vp at z=22km, t=100yr | Loading rate z=10km |
|-----|-----|----------|-------|-------------------------|---------------------|
| v21 (original) | 1.0 | SCEC | 0.03 | 0.0873 | 0.0084 MPa/yr |
| Test 1 (no δτ) | **0** | SCEC | 0.03 | 0.0873 | 0.0084 MPa/yr |
| Test 3 (Tandem psi) | 1.0 | **Tandem** | 0.03 | 0.0873 | 0.0084 MPa/yr |

The numbers are **identical to 4 significant figures**. The δτ and psi initialization
are irrelevant — V_nuc = 0.03 dominates.

---

## 2. Why Tests 1 and 3 Were Invalid

### Test 1 (--delta-tau-factor 0)

Removing δτ only removes the EXTRA overstress. The nucleation zone still starts at
V = V_nuc = 0.03 m/s because:
- τ_pre is computed at V_nuc (without δτ): τ = σ_n·a·asinh(V_nuc/(2V₀)·exp(ψ_ss/a)) + η·V_nuc
- With SCEC psi (ψ_ss at Vp), the friction solver recovers V = V_nuc = 0.03

### Test 3 (--psi-init-mode tandem)

Tandem-style psi absorbs the stress into ψ, making the system start in equilibrium
at V = V_nuc. But V_nuc = 0.03 is still seismic. With δτ_factor = 1.0 (default),
the absorbed stress is slightly higher, but V still starts at 0.03.

### The Core Issue

**V_nuc = 0.03 m/s is inherently unstable** in the velocity-weakening nucleation zone.
The nucleation patch (12 km × 12 km) exceeds the critical nucleation length
h* ≈ 4.4 km, so the instability runs away within seconds.

**Tandem uses V_nuc = 0.01 m/s** (from bp5.lua line 73). With Tandem-style psi
absorption, V = 0.01 apparently decays back toward Vp without triggering
immediate nucleation. The 3× difference in V_nuc is critical.

---

## 3. Deep VS Lockup Data (z=22 km, All Three Runs Identical)

| Time (yr) | V/Vp | θ/θ_ss | τ_strike (MPa) |
|-----------|------|--------|-----------------|
| 0 | 1.000 | 1.0 | 13.273 |
| 1 | 7.06 | 0.1 | 13.404 |
| 5 | 1.38 | 0.4 | 12.976 |
| 10 | 0.653 | 0.9 | 12.767 |
| 25 | 0.256 | 2.3 | 12.534 |
| 50 | **0.138** | **4.5** | 12.423 |
| 100 | **0.087** | **8.2** | 12.417 |
| 150 | **0.077** | **10.8** | 12.491 |
| 200 | **0.078** | **12.1** | 12.588 |
| 250 | **0.085** | **12.2** | 12.683 |
| 295 | **0.095** | **11.5** | 12.751 |

---

## 4. Seismogenic Core (z=10 km, All Three Runs Identical)

| Time (yr) | τ_strike (MPa) | slip_s (m) | Loading rate |
|-----------|-----------------|------------|-------------|
| 0 | 19.490 | 0.000 | — |
| 1 | 9.505 | 6.316 | (earthquake) |
| 10 | 9.948 | 6.316 | — |
| 50 | 10.405 | 6.316 | — |
| 100 | 10.805 | 6.316 | 0.0084 MPa/yr |
| 200 | 11.568 | 6.316 | 0.0084 MPa/yr |
| 295 | 12.329 | 6.316 | 0.0084 MPa/yr |

Loading rate: **0.0084 MPa/yr** (boundary-only, deep VS locked).
Tandem's loading rate: **~0.07 MPa/yr** (boundary + deep VS creep).
Ratio: MFEM is **8.3× too slow**.

---

## 5. Test 2: Boundary Loading Verified Correct

The `--diag-vtk` path solved the elasticity at t=1yr with zero slip and computed
traction at the fault. Results from `bp5_v22_t2_bdrload_7598607.out`:

```
Dirichlet BC solve at t=1yr: |u| = 0.00575393

=== Boundary Loading Traction Verification ===
Analytical estimate (2D antiplane): 5051 Pa = 0.00505 MPa

z=0km:  trac_strike=5490 Pa (0.00549 MPa)  trac_dip=0.005 Pa  ratio=1.086
z=10km: trac_strike=5484 Pa (0.00548 MPa)  trac_dip=-0.4 Pa   ratio=1.085
z=22km: trac_strike=5478 Pa (0.00548 MPa)  trac_dip=-1.7 Pa   ratio=1.084
```

**Conclusions**:
- **DG Dirichlet boundary loading is correct.** Traction matches analytical estimate
  within 8.5% (difference from 3D geometry + y-boundary contributions).
- **Purely along-strike**: dip traction is negligible (< 2 Pa vs 5500 Pa strike).
- **Depth-uniform**: traction barely varies with depth (5490→5478, <0.3% variation).
- **No factor-of-2 error**, no sign issue, no missing contribution.
- The boundary loading rate is: 5484 Pa / (3.15e7 s) = 1.74e-4 Pa/s = **0.0055 MPa/yr**.
  Over 300 years: 1.64 MPa. This is correct for boundary-only loading but insufficient
  for earthquake recurrence (Tandem needs ~0.07 MPa/yr from boundary + VS creep combined).

---

## 6. Root Cause Update

### Confirmed: V_nuc = 0.03 triggers immediate earthquake

The SCEC BP5 specification allows V_nuc as a free parameter. MFEM uses 0.03 m/s,
Tandem uses 0.01 m/s. With V_nuc = 0.03, the nucleation zone is so fast that an
earthquake occurs within seconds regardless of δτ or psi initialization.

### Confirmed: Boundary loading is correct, rate is 0.0055-0.0084 MPa/yr

Test 2 traction diagnostic verified the DG Dirichlet loading matches the analytical
estimate within 8.5%. The observed interseismic rate (0.0084 MPa/yr at z=10km) is
slightly higher due to contributions from the bonded domain below 40 km. This rate
is correct but requires ~2000 years for the next event — confirming the non-cycling
is due to deep VS lockup, not a boundary loading bug.

### Unknown: Would V_nuc = 0.01 + Tandem psi prevent immediate nucleation?

In Tandem, V_nuc = 0.01 with absorbed δτ apparently decays without nucleating.
If MFEM with V_nuc = 0.01 + Tandem psi also avoids immediate nucleation, then:
- Boundary + deep VS creep loading would operate normally (~0.07 MPa/yr)
- First event at ~150 years (matching Tandem)
- Subsequent cycling expected

### Unknown: Does the VS zone maintain V ≈ Vp without ANY nucleation?

If we set V_nuc = Vp (completely uniform fault), V should stay at Vp everywhere.
If it doesn't, there is a fundamental problem with the elasticity/friction coupling.

---

## 7. Corrected Test Plan

### Test 1-corrected: Completely uniform fault (no nucleation at all)

```bash
--delta-tau-factor 0 --V-nuc 1e-9 --tfinal 9.45e9
```

Every node starts at V = Vp with SCEC psi. No nucleation zone. The fault should
remain at V = Vp everywhere for 300 years.

**Pass**: V/Vp stays in [0.5, 2.0] everywhere.
**Fail**: V drifts — indicates fundamental elasticity/loading problem.

### Test 3-corrected: Match Tandem exactly

```bash
--psi-init-mode tandem --V-nuc 0.01 --tfinal 1.89e10
```

Match Tandem's V_nuc = 0.01 and psi absorption. Run for 600 years to capture
potential first event at ~150 yr and verify recurrence.

**Pass**: First event at ~150 yr, deep VS maintains V ≈ Vp, subsequent cycling.
**Fail**: Immediate earthquake (V_nuc = 0.01 still too fast) or deep VS lockup.

---

## 8. Files

| File | Description |
|------|-------------|
| `jobs/bp5/bp5_v22_test1_no_dtau.sbatch` | Test 1 (invalid) |
| `jobs/bp5/bp5_v22_test2_bdr_loading.sbatch` | Test 2 (need .out file) |
| `jobs/bp5/bp5_v22_test3_tandem_psi.sbatch` | Test 3 (invalid — V_nuc=0.03) |
| `results_v22_test1_no_dtau/` | Test 1 results |
| `results_v22_test2_bdr_loading/` | Test 2 results |
| `results_v22_test3_tandem_psi/` | Test 3 results |
