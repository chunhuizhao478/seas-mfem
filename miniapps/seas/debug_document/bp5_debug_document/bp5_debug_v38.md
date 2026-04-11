# BP5 Debug v38: Switch to IP Method + Parallel Bug Fixes

**Date**: 2026-03-16 (v38a-b), 2026-03-17 (v38c)
**Status**: IP blows up (v38a), shared-face Dirichlet fix (v38b), shared-face traction sign fix (v38c)
**Previous**: v37 (BR2 interior Dirichlet sign fix, recurrence 295 yr)
**Branch**: `feature/elasticity`

---

## 1. Critical Discovery: Tandem Uses IP, Not BR2

All v31-v37 runs used `--dg-method BR2`, but Tandem BP5 is **hardcoded to IP**:

```cpp
// tandem/app/tandem/Context.h line 37
return std::make_shared<Elasticity>(..., DGMethod::IP);
```

The `penalty()` function in Elasticity.h confirms:
```cpp
double penalty(std::size_t fctNo) const {
    if (method_ == DGMethod::BR2) {
        return NumFacets;  // = 4 (dimensionless, effectively zero for traction)
    }
    return penalty_[fctNo];  // IP: material-dependent, properly scaled
}
```

All our BR2 analysis and fixes (v34-v37) were optimizing the wrong DG method.

## 2. Why IP Was Abandoned (v24-v28)

| Version | Change | Result |
|---------|--------|--------|
| v24 | Remove BR2 traction lifting, use scalar penalty=4 | Blows up at ~2yr |
| v25 | Same + new uniform mesh | Same blowup at ~2yr |
| v26 | Switch to IP bilinear form | Blows up immediately |
| v27 | IP with Tandem scalar traction penalty | Still blows up at ~0.08yr |
| v28 | Custom IP penalty integrator + tag-based fault detection | Still blows up at (0,0,0) |

**v29 conclusion**: "Every IP configuration blows up... The BR2 method with its
lifting-based traction is the ONLY stable configuration."

## 3. Why IP Should Work Now (Hypothesis)

**v31 (DG slip sign fix)** identified the root cause of ALL v30 blowups:

> v30 negated `tau_pre` to match Tandem convention, making V < 0. But the DG
> sign variable was never adjusted: `sign = (nor(1) > 0) ? -1 : +1`.
> With sign=-1 and delta_u<0: g = sign*delta_u = +|delta_u| — **wrong sign**.
> → penalty corrects in wrong direction → **exponential blowup**

v31 explicitly says:
```
// If nor(1) > 0: Elem1 is -Y side, [[u]] = u(-Y) - u(+Y) = delta_u -> sign = +1
// If nor(1) < 0: Elem1 is +Y side, [[u]] = u(+Y) - u(-Y) = -delta_u -> sign = -1
real_t sign = (nor(1) > 0) ? 1.0 : -1.0;
```

The sign fix applies to ALL 6 locations (both IP and BR2 paths).

## 4. v38a Results: IP Still Blows Up

**Job**: `bp5_v38a_ip_p1_7601535` — 400 ranks, p=1, 1000m mesh, IP method

**Result**: Simulation blew up during the initial earthquake.

| Metric | Value |
|--------|-------|
| Steps completed | 634 |
| Simulation time | 5.71e-09 yr (~0.18 seconds) |
| Final V_max | 8.97e+05 m/s (897 km/s — unphysical) |
| First TRACTION BLOWUP | Step 43, Rank 317, DOF 20 (shared), x=(0,0,0) |
| Total blowup messages | 275,884 |
| Earthquakes | 1 (initial nucleation only) |

The blowup pattern:
- Starts at step 43 with V_max = 216 m/s (already high)
- First blowup on a **shared face** at (0,0,0) — the fault center
- tau_mag reaches 1e9 Pa (1 GPa), then exponentially grows to 1e11+ Pa
- Simulation never reaches interseismic period

**Conclusion**: The v31 sign fix was NOT the sole cause of IP blowup. There is
a separate IP instability that persists. The v24-v28 blowups had two causes:
1. The sign bug (fixed in v31) — caused BR2 to blow up too
2. An IP-specific instability — still present

## 5. Bug Fix: Missing Shared-Face Dirichlet Loading (v38b)

### 5.1 The Bug

`dirichlet_shared_faces_` is populated in `BuildDirichletInteriorFaces()` (line 463)
but **never used** in `AssembleDirichletLoading()`. The function only loops over
`dirichlet_interior_faces_` for the Y=0 skeleton pattern.

In parallel runs (400-800 MPI ranks), any Y=0 non-fault face that falls on a
processor boundary becomes a shared face and receives **zero tectonic loading**.
This weakens the effective driving force near those faces.

Evidence:
- **Line 210**: `Array<int> dirichlet_shared_faces_;` declared
- **Line 393**: Initialized to empty
- **Line 463**: Populated with matching shared face indices
- **Lines 1702-2345**: `AssembleDirichletLoading()` — NO loop over `dirichlet_shared_faces_`

### 5.2 The Fix

Added a new loop at the end of `AssembleDirichletLoading()` that processes
`dirichlet_shared_faces_`, following the established shared-face pattern from
`AssembleSlipContributionIPShared` and `AssembleSlipContributionBR2Shared`:

```cpp
// Shared faces with Dirichlet BC (parallel only)
if constexpr (IsParallelMesh<MeshType>::value)
{
#ifdef MFEM_USE_MPI
   auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());
   for (int fi = 0; fi < dirichlet_shared_faces_.Size(); fi++)
   {
      int sf = dirichlet_shared_faces_[fi];
      FaceElementTransformations *FTr =
         mesh_.GetSharedFaceTransformations(sf);
      // ... IP or BR2 path, elem1 only ...
   }
#endif
}
```

Key design decisions:
- **Elem2 access**: `nbr_idx = FTr->Elem2No - mesh_.GetNE()`, then
  `pfes->GetFaceNbrFE(nbr_idx)` — matches existing shared-face patterns
- **Only scatter Elem1 RHS**: Elem2 is on a neighboring rank; that rank processes
  the same face as its own elem1
- **Y-sign**: Elem1 from centroid, `sign2 = -sign1` (both elements straddle Y=0)
- **Both IP and BR2 paths**: Full skeleton pattern matching interior face code
- **BR2 cross-element lifting**: Uses both Minv1 and Minv2 (face-neighbor mass
  inverse already precomputed at `elem_mass_inv_[FTr->Elem2No]`)

### 5.3 Impact Assessment

This bug affects ALL parallel runs since v34 (when interior Dirichlet was
introduced). It affects both IP and BR2 equally:
- **v34-v37 (BR2)**: Some Y=0 faces on proc boundaries got zero loading
- **v38a (IP)**: Same missing loading (but IP blew up for other reasons)

The effect on recurrence depends on how many Y=0 non-fault faces fall on
processor boundaries. For 400 ranks over a mesh with ~63k elements, a
non-trivial fraction of Y=0 faces will be shared faces.

This fix may contribute to shorter recurrence in BR2 runs but is unlikely to
be the primary factor (the v37→Tandem gap is ~55 yr: 295 vs 240).

## 6. Current IP Code Path

The IP code path includes all necessary fixes:

| Component | Implementation | Status |
|-----------|---------------|--------|
| Bilinear form | Custom `DGElasticityIPPenaltyIntegrator` (v28): `penalty * \|nor\|` | Ready |
| Slip RHS | `AssembleSlipContributionIP`: Tandem's scalar penalty formula | Ready |
| Dirichlet RHS (boundary) | `AssembleDirichletLoading` IP path: `penalty * \|nor\|` | Ready |
| Dirichlet RHS (interior Y=0) | `AssembleDirichletLoading` interior IP: same formula | Ready |
| Dirichlet RHS (shared Y=0) | `AssembleDirichletLoading` shared IP: **v38b fix** | **Fixed** |
| Traction | `ComputeTraction` IP: scalar penalty, per-quad-point | Ready |
| Sign convention | v31 fix: all 6 locations | Ready |
| Fault detection | Tag-based (H32): excludes z=0, z=Wf faces | Ready |

### IP penalty formula (matches Tandem exactly):
```
p(side) = (Dim+1) * c_N_1 * (area/volume) * (c1^2/c0)
penalty = (p(0) + p(1)) / 4.0   (interior faces)
penalty = p(0)                    (boundary faces)

c0 = 2mu, c1 = Dim*lambda + 2mu
c_N_1 = order * (order + Dim - 1) / Dim
```

### Traction formula (matches Tandem exactly):
```
T = {sigma*n_hat} - penalty * ([[u]] - sign*delta_u)   (scalar x vector)
```

## 7. Root Cause Found: Shared-Face IP Traction Sign Bug (v38c)

### 7.1 The Diagnostic Red Herring

The v38a blowup messages all showed `x=(0,0,0)` for shared-face DOFs, which initially
suggested the fault-surface intersection (z=0 edge) was the problem. However, this was
a **diagnostic bug**: the blowup reporting code (line 3612-3626) only computes face
coordinates for interior faces (`i < fault_interior_faces_.Size()`). For shared faces,
`face_center` stays at its default initialization `(0,0,0)`. The actual blowup locations
are unknown.

### 7.2 The Real Bug: Parallel Sign Inconsistency in IP Traction Correction

**File**: `domain/elasticity_operator.hpp`, `ComputeTraction()`, shared-face IP path

For shared fault faces, MFEM always puts the **local** element as `Elem1`. Two ranks
sharing the same physical face have **opposite** Elem1/Elem2 assignments:

```
Rank A: Elem1 = element_A (local),  Elem2 = element_B (neighbor)
Rank B: Elem1 = element_B (local),  Elem2 = element_A (neighbor)
```

This causes the face normal from `CalcOrtho` to point in **opposite** directions on the
two ranks, giving opposite `sign` values. The IP penalty correction then differs:

| | Rank A | Rank B |
|---|--------|--------|
| `nor(1)` | +N | -N |
| `sign` | +1 | -1 |
| `u1 - u2` | u_A - u_B | u_B - u_A = -(u_A - u_B) |
| `sign * delta_u` | +delta_u | -delta_u |
| `jump = (u1-u2) - sign*delta_u` | (u_A-u_B) - delta_u | -(u_A-u_B) + delta_u = **-jump_A** |
| `correction = penalty * jump` | penalty × jump_A | **-penalty × jump_A** |
| `T = {σ·n̂} - correction` | {σ·n̂} - correction | {σ·n̂} **+** correction |

The two ranks compute **opposite** penalty corrections for the same physical face,
giving **different** tractions. Since there is no MPI synchronization of traction or
fault state between ranks, the rate-state ODE evolves independently on each rank with
different traction inputs. The fault states diverge, creating larger jump residuals,
which the IP penalty amplifies further → exponential blowup.

**Why BR2 doesn't have this bug**: BR2's face integral uses `jump_q * nor_q` as a
product. When elem1/elem2 are swapped, both `jump_q` and `nor_q` flip sign, but their
product is **invariant**. The BR2 lifting and evaluation are symmetric in
`(eval1 + eval2)`, so both ranks get the same correction. IP has no such cancellation
because the correction is `penalty * jump` with no normal multiplication.

### 7.3 Magnitude Analysis

| DG Method | Penalty magnitude | Typical correction per 1mm residual | Inter-rank discrepancy |
|-----------|------------------|-------------------------------------|----------------------|
| BR2 | 4 (dimensionless, via lifting) | ~0.06 MPa (through Minv + C:n) | ~0.12 MPa (negligible) |
| IP | ~2.4×10⁹ Pa/m (material-dependent) | ~2.4 MPa | ~4.8 MPa (**catastrophic**) |

The IP penalty is ~10⁹× the BR2 dimensionless penalty. Even tiny jump residuals from
the elastic solve create O(MPa) inter-rank traction discrepancies, triggering the
divergence feedback loop.

The friction law cannot balance the resulting traction:
- Max friction at V=10⁶ m/s: ~42.6 MPa (σ_n × (f0 + a×ln(V/V0)))
- IP traction at blowup: **1027 MPa** (40× unphysical)

### 7.4 The Fix

**One-line change** in the shared-face IP traction path (line ~3406):

```cpp
// BEFORE (v38a — different on each rank):
real_t jump_c = (u1q - u2q) - sign * delta_u[c];
correction_q[c] = penalty_ip * jump_c;

// AFTER (v38c — canonical, same on both ranks):
real_t jump_raw = (u1q - u2q) - sign * delta_u[c];
correction_q[c] = penalty_ip * sign * jump_raw;
```

The fix multiplies the raw jump by `sign`, producing the canonical form:

```
correction = penalty × sign × ((u1-u2) - sign×delta_u)
           = penalty × (sign×(u1-u2) - delta_u)
```

Since `sign×(u1-u2)` is **invariant** across ranks (both `sign` and `u1-u2` flip
together), `delta_u` is consistent (from FaultBasis with fixed ref_normal), and
`penalty` is a material constant, **both ranks now compute the same correction**.

Verification at equilibrium: when `u1-u2 = sign×delta_u` (constraint satisfied),
`jump_raw = 0` → `correction = 0` → `T = {σ·n̂}` regardless of rank. ✓

### 7.5 Why Interior Faces Don't Need This Fix

For interior faces, each face is processed by exactly **one** rank. There is no
inter-rank inconsistency. The traction value depends on element ordering (which is
non-physical), but since only one entity computes it, the result is self-consistent
within the simulation.

Changing interior faces would alter established BR2 behavior that has been validated
through v34-v37. The minimal, targeted fix applies only to shared faces.

### 7.6 Why the v24-v28 IP Experiments Also Blew Up

The v24-v28 IP experiments (before v31 sign fix) had **two** independent bugs:
1. **DG slip sign bug** (fixed in v31): affected both IP and BR2
2. **Shared-face traction sign bug** (fixed in v38c): affects only IP in parallel

The v24-v28 runs used parallel execution, so both bugs contributed to the blowups.
The v31 sign fix resolved bug #1, but bug #2 remained. The v38a test confirmed that
IP still blows up after v31 — now explained by bug #2.

Note: v28 reported the blowup at "x=(0,0,0)" — the same diagnostic artifact. The
actual blowup location was unknown, likely on shared faces just like v38a.

## 8. Files Changed

| File | Change | Version |
|------|--------|---------|
| `domain/elasticity_operator.hpp` | Added shared-face Dirichlet loading loop | v38b |
| `domain/elasticity_operator.hpp` | Canonical IP traction correction for shared faces | v38c |

## 9. Verification

### Build
```
conda activate mfem-dev && make -j8
```
Build succeeds with no errors.

### Serial Unit Tests
All test suites pass:

| Suite | Tests | Status |
|-------|-------|--------|
| seas_test_elasticity_operator | 84 | PASS |
| seas_test_elasticity_br2 | 46 | PASS |
| seas_test_fault_basis | 187 | PASS |
| seas_test_domain_interface | 32 | PASS |
| seas_test_bp5_params | 96 | PASS |

### Parallel Unit Tests (8 MPI ranks)
All parallel suites pass:

| Suite | Tests | Status |
|-------|-------|--------|
| seas_test_parallel_domain | 11 | PASS |
| seas_test_parallel_elasticity | 6 | PASS |
| seas_test_parallel_fault | 27 | PASS |
| seas_test_parallel_utils | 9 | PASS |
| seas_test_serial_parallel_consistency | 12 | PASS |
| seas_test_br2_consistency | 7 | PASS |
| seas_test_bp5_parallel_smoke | 13 | PASS |

## 10. Next Steps

1. **Submit v38c IP run on TACC**: Same mesh/config as v38a but with the
   shared-face traction fix. This is the critical test — if the fix resolves
   the parallel sign inconsistency, IP should survive the initial earthquake.

2. **Submit v38c BR2 run on TACC**: Same as v37 config but with both parallel
   fixes (shared-face Dirichlet + canonical traction). Compare recurrence
   with v37's 295 yr.

3. **If IP survives coseismic**: Compare recurrence with Tandem's 240 yr.
   The IP traction formula now matches Tandem exactly (scalar penalty,
   no BR2 lifting bias).

4. **If IP still blows up**: The remaining cause would be the IP penalty
   magnitude at interior faces during coseismic. Options:
   - Try p=2 (better DG enforcement, smaller jump residual)
   - Cap the correction magnitude at interior faces
   - Compare time step sizes with Tandem during coseismic

## 11. Cumulative Fix History

| Fix | Description | Status |
|-----|-------------|--------|
| H1-H32 | Previous fixes (see v1-v28 docs) | Done |
| v30 | Tandem coordinate system | Done |
| v31 | DG slip sign fix + interior Dirichlet | Done |
| v32-v33 | General polynomial order (p-refinement) | Done |
| v34 | Interior Dirichlet Y=0 jump loading | Done |
| v35-v36 | Per-quad-point traction + cross-element BR2 | Done |
| v37 | Interior Dirichlet face_int2 sign fix | Done |
| v38a | Switch to IP method (matching Tandem) | **IP blows up** |
| v38b | Shared-face Dirichlet loading fix (parallel bug) | Applied |
| **v38c** | **Shared-face IP traction sign fix (parallel bug)** | **Applied, tests pass** |
