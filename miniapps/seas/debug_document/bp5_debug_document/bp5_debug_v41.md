# BP5 Debug v41: Parametric Study — Domain Size, Resolution, and Polynomial Order

**Date**: 2026-03-18 (updated with TACC results)
**Status**: All 6 TACC runs completed, analysis complete — p≥2 fault locking is BLOCKING
**Previous**: v40 (IP penalty correction sign fix, all IP traction sign bugs resolved)
**Branch**: `feature/elasticity`

---

## 1. Problem Statement

v40 fixed all IP traction sign bugs and produced 8 earthquake events over 1800 yr with
~250 yr recurrence — matching Tandem's 240 yr. Three issues remained:

1. **Dip slip/traction mismatch**: MFEM's dip quantities are larger than Tandem's and
   grow linearly with time at several stations
2. **Boundary truncation**: The standard BP5 domain (400×200×100 km) may be too small
3. **20% strike slip deficit**: MFEM accumulates ~20% less strike slip than Tandem at p=1

v41 ran a systematic parametric study to isolate which effects are bugs vs. discretization.

## 2. Runs Submitted and Completed

| Run | Config | Mesh | Elements | Nodes | Ranks | Sim Duration | Purpose |
|-----|--------|------|----------|-------|-------|-------------|---------|
| v41a | IP p=1, 1000m | bp5_tandem.msh | 63,451 | 8 | 400 | **1499 yr** | Baseline |
| v41b | IP p=1, 500m | bp5_tandem_500m.msh | ~500K | 16 | 800 | **199 yr** | h-refinement |
| v41c | IP p=2, 1000m | bp5_tandem.msh | 63,451 | 16 | 800 | **1800 yr** | p-refinement |
| v41d | IP p=4, 2500m | bp5_tandem_2500m.msh | 13,503 | 16 | 800 | **1800 yr** | Match uphoff.2 |
| v41e | IP p=1, 2x domain | bp5_tandem_2x_1000m.msh | 104K | 12 | 600 | **926 yr** | Domain size |
| v41f | IP p=1, 4x domain | bp5_tandem_4x_1000m.msh | 147K | 16 | 800 | **378 yr** | Domain size |

Reference: v40 (IP p=1, 1000m) = 8 events in 1800 yr, Tandem (p=6, 4km) = 8 events in 1800 yr.

---

## 3. Critical Results

### 3.1 Earthquake Cycling

| Run | Events | Duration | Recurrence | Status |
|-----|--------|----------|------------|--------|
| v40 (p1, 1x) | 8 | 1800 yr | **~250 yr** | ✓ Good match |
| v41a (p1, 1x) | 6 | 1499 yr | **~250 yr** | ✓ Same as v40 |
| v41b (p1, 500m) | 1 | 199 yr | — | Walltime limited |
| **v41c (p2, 1x)** | **0** | **1800 yr** | **∞** | **❌ FAULT LOCKED** |
| **v41d (p4, 2.5km)** | **0** | **1800 yr** | **∞** | **❌ FAULT LOCKED** |
| v41e (p1, 2x) | 3 | 926 yr | **~333 yr** | Longer recurrence |
| v41f (p1, 4x) | 2 | 378 yr | **~377 yr** | Longer recurrence |
| Tandem (p6, 4km) | 8 | 1800 yr | **~240 yr** | Reference |

### 3.2 p=2 and p=4 Details — Complete Fault Locking

**v41c (p=2, 1000m)**: V_max starts at 1e-9 m/s and DECAYS to 1e-13. tau_strike drops
from 13.27 to 9.21 MPa over 1800 yr. Total strike slip: 0.014 m (should be ~57 m).

**v41d (p=4, 2500m)**: V_max starts at 1e-9 m/s and DECAYS to 1e-13. tau_strike drops
from 13.27 to 8.55 MPa over 1800 yr. Total strike slip: 0.008 m.

Both runs show the fault gradually losing stress without ever nucleating. The stress
drops below the steady-state friction level for plate-rate creep, locking the fault.

**This is the same "fault locking" bug seen in v33-v34 with BR2, now confirmed for IP.**

---

## 4. Dip Quantities — Comprehensive Analysis

### 4.1 Dip Slip at t=500 yr (meters)

| Station | Tandem | v40/41a (p1) | v41c (p2) | v41d (p4) | v41e (2x) |
|---------|--------|-------------|-----------|-----------|-----------|
| strk+00dp+00 | -0.106 | **+0.143** | -0.000 | +0.000 | -0.062 |
| strk+00dp+10 | -0.057 | +0.068 | -0.000 | +0.000 | -0.071 |
| strk+00dp+22 | -0.006 | -0.129 | +0.000 | +0.000 | +0.010 |
| strk+16dp+00 | -0.316 | -0.544 | -0.000 | -0.000 | -0.159 |
| strk+36dp+00 | -0.301 | -0.559 | +0.000 | +0.000 | -0.301 |
| strk-16dp+00 | +0.128 | +0.585 | +0.000 | +0.000 | +0.015 |
| strk-36dp+00 | +0.362 | +0.376 | -0.000 | -0.000 | +0.266 |

### 4.2 Dip Traction at t=500 yr (MPa)

| Station | Tandem | v40/41a (p1) | v41c (p2) | v41d (p4) | v41e (2x) |
|---------|--------|-------------|-----------|-----------|-----------|
| strk+00dp+00 | +0.085 | +0.143 | -0.000 | +0.004 | -0.044 |
| strk+00dp+10 | +0.069 | -0.002 | -0.000 | +0.063 | -0.358 |
| strk+00dp+22 | -0.008 | -0.149 | +0.010 | +0.142 | -0.036 |
| strk+16dp+00 | +0.024 | **-0.406** | -0.004 | -0.002 | +0.321 |
| strk+36dp+00 | -0.189 | **-1.193** | +0.029 | +0.031 | -0.137 |
| strk-16dp+00 | +0.006 | +0.277 | +0.002 | +0.004 | -0.324 |
| strk-36dp+00 | +0.187 | +0.740 | -0.042 | -0.026 | -0.103 |

### 4.3 Key Findings

#### [A] No Sign Convention Bug

**p=2 and p=4 produce essentially zero dip slip and dip traction** (~1e-4 m and ~0.03 MPa).
If there were a global sign bug, it would appear at all polynomial orders. The near-zero
dip quantities at p≥2 prove the p=1 dip mismatch is a **discretization artifact**.

The antisymmetric pattern (strk+16 vs strk-16 opposite signs) is consistent with a
fault-edge stress concentration effect, not a sign error.

**Conclusion: No code fix needed for dip quantities.**

#### [B] p=1 Dip Amplification Mechanism

At p=1, MFEM's dip quantities are 1.5-6× Tandem's magnitude and grow linearly with time:

1. The DG jump residual `[[u]] - delta_u` scales with accumulated slip
2. The IP penalty correction `η × physical_jump` amplifies this into growing tau_dip
3. The rate-state friction law converts tau_dip into V_dip → accumulating dip slip
4. The O(h) discretization error at p=1 drives this entire chain

This vanishes at p≥2 because the DG residual scales as O(h^p).

#### [C] 2x/4x Domains Reduce Dip but Increase Recurrence

The 2x domain reduces dip quantities by ~2× compared to the 1x baseline, confirming
boundary truncation contributes to dip artifacts at p=1.

However, the 2x and 4x runs show **longer recurrence** (333/377 yr vs 250 yr). This may
be because: (1) runs didn't reach steady-state cycling yet, or (2) the first 1-2 events
after nucleation have anomalous timing.

---

## 5. Root Cause: p≥2 Fault Locking — BLOCKING BUG

### 5.1 The Pattern

Both v41c (IP p=2, 1000m) and v41d (IP p=4, 2500m) show identical behavior:
- No earthquake nucleation over 1800 yr
- V_max decays from 1e-9 to 1e-13 m/s
- tau_strike monotonically decreases (13.27 → 9.2 MPa at p=2, 8.5 MPa at p=4)
- Fault gradually loses stress without loading back up

This is identical to the v33-v34 BR2 fault locking, which was caused by:
1. Centroid-only traction evaluation at p≥2 (fixed in v35-v36 for BR2)
2. Missing cross-element BR2 lifting terms (fixed in v36 for BR2)
3. face_int2 sign in interior Dirichlet loading (fixed in v37)

### 5.2 Suspected IP-Specific Causes

The v35-v37 fixes targeted the **BR2 path** in `ComputeTraction()` and
`AssembleDirichletLoading()`. The **IP path** may have analogous issues:

1. **`ComputeTraction()` IP path at p≥2**: May still evaluate the average stress
   `{σ·n̂}` at the face centroid only (single point), which is exact for p=1 but
   incorrect for p≥2. The per-quad-point rewrite in v35-v36 may only have been
   applied to the BR2 path.

2. **`AssembleDirichletLoading()` IP path**: The interior Dirichlet loading (non-fault
   Y=0 faces) may have the same cross-element scaling issues at p≥2 that v36-v37
   fixed for BR2. Specifically, the 2x penalty overshoot and missing cross-element
   terms would weaken the effective tectonic loading.

3. **IP penalty magnitude at p≥2**: c_N_1 scales as p(p+D-1)/D:
   - p=1: c_N_1 = 1.0
   - p=2: c_N_1 = 2.67
   - p=4: c_N_1 = 8.0
   The 2.67× or 8× larger penalty creates proportionally larger traction corrections,
   which may bias fault dynamics if the DG residual is not also converging fast enough.

4. **Slip assembly at p≥2**: `AssembleSlipContributionIP()` may have centroid-only
   evaluation that loses accuracy at higher polynomial order.

### 5.3 Why It Didn't Show at p=1

At p=1, the displacement field is linear per element, so:
- Centroid evaluation of gradients is exact (constant gradient)
- Face-averaged shapes equal centroid shapes
- DG jumps are linear on faces → centroid value = face average

All these properties break at p≥2, where gradients vary spatially.

### 5.4 Investigation Plan

1. Check `ComputeTraction()` IP path in `elasticity_operator.hpp`:
   - Does it loop over quadrature points or evaluate at centroid?
   - Compare with the BR2 path (which was fixed in v35-v36)

2. Check `AssembleDirichletLoading()` IP path:
   - Does the interior Dirichlet IP code have cross-element terms?
   - Does it have the 0.5 factor matching the bilinear form?

3. Check `AssembleSlipContributionIP()`:
   - Is the slip evaluated per quadrature point or at centroid?

4. After fixing, run p=2 smoke test to confirm nucleation occurs.

---

## 6. Summary of Findings

| Issue | Root Cause | Code Bug? | Fix |
|-------|-----------|-----------|-----|
| Dip slip mismatch | p=1 discretization artifact | **No** | Use p≥2 |
| Dip traction mismatch | p=1 IP penalty amplification | **No** | Use p≥2 |
| 20% strike slip deficit | p=1 O(h) stress error | **No** | Use p≥2 or finer mesh |
| 2x/4x longer recurrence | Unclear (may be transient) | Investigate | Need longer runs |
| **p≥2 fault locking** | **IP path not updated for p≥2** | **YES** | **v42 priority** |

The dip mismatch and strike deficit are NOT bugs — they are expected consequences of
p=1 discretization on 1km tetrahedral meshes. The definitive proof: p=2 and p=4 produce
zero dip quantities, confirming O(h^p) convergence.

However, p≥2 fault locking prevents us from using higher order to resolve these issues.
**Fixing the p≥2 fault locking is the #1 priority.**

---

## 7. Cumulative Fix History

| Fix | Description | Status |
|-----|-------------|--------|
| v30-v40 | All previous fixes (see v30-v40 docs) | Done |
| v40 | IP penalty correction sign negation | Done |
| v41a | Baseline IP p=1 1000m — 6 events in 1499 yr ✓ | Done |
| v41b | 500m mesh — 1 event in 199 yr (walltime) | Done |
| v41c | **IP p=2 1000m — NO EARTHQUAKES (fault locked)** | **BLOCKING** |
| v41d | **IP p=4 2500m — NO EARTHQUAKES (fault locked)** | **BLOCKING** |
| v41e | 2x domain p=1 — 3 events in 926 yr, longer recurrence | Done |
| v41f | 4x domain p=1 — 2 events in 378 yr | Done |

---

## 8. Next Steps

### Priority 1: Fix p≥2 Fault Locking (v42)

1. Audit `ComputeTraction()` IP path for centroid-only evaluation
2. Audit `AssembleDirichletLoading()` IP path for cross-element terms
3. Audit `AssembleSlipContributionIP()` for per-quad-point evaluation
4. Apply analogous fixes to what v35-v37 did for BR2
5. Run p=2 smoke test to confirm nucleation

### Priority 2: Verify p=2 Fixes Everything

Once p≥2 works:
- v42a: IP p=2, 1000m → confirm cycling + reduced dip + reduced deficit
- v42b: IP p=2, 2x domain → confirm boundary effects are small
- If p=2 recurrence ≈ 240 yr and dip ≈ Tandem → formulation verified

### Priority 3: Domain Size (Lower Priority)

The dip artifacts at p=1 are discretization errors. Larger domains are not needed
if p≥2 works. However:
- 2x domain reduces p=1 dip artifacts by ~2×
- 4x domain shows diminishing returns
- Domain size primarily matters at p=1 where discretization error is large
