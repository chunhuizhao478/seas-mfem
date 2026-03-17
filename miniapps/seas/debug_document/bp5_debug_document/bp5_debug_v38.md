# BP5 Debug v38: Switch to IP Method + Shared-Face Dirichlet Fix

**Date**: 2026-03-16
**Status**: IP blows up (v38a), shared-face Dirichlet fix applied (v38b)
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

## 7. Open Question: What Causes IP Blowup?

The v38a blowup occurs during the initial coseismic phase (V_max > 200 m/s)
at the fault center (0,0,0). Possible causes:

1. **IP penalty too large during coseismic**: The penalty scales as
   `(area/volume) * c1^2/c0`. During rapid slip, the DG jump `[[u]] - delta_u`
   is large, and the scalar penalty may overcorrect, creating oscillations.
   BR2 avoids this because its traction correction uses the anisotropic
   elasticity tensor, which distributes the correction more naturally.

2. **IP bilinear form mismatch**: The custom `DGElasticityIPPenaltyIntegrator`
   uses `penalty * |nor|` scaling. If this doesn't exactly match the traction
   computation's penalty formula, there could be a consistency gap that
   manifests as instability during high slip rates.

3. **Time step control**: During coseismic, the adaptive time stepper may not
   reduce dt enough for IP's stiffer penalty. BR2's softer penalty (sigma=4,
   dimensionless) may be more forgiving of large time steps.

## 8. Files Changed

| File | Change |
|------|--------|
| `domain/elasticity_operator.hpp` | Added shared-face Dirichlet loading loop (v38b) |

## 9. Verification

### Build
```
conda activate mfem-dev && make -j8
```
Build succeeds with no errors.

### Serial Unit Tests
All test suites pass (identical to v37).

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
| seas_test_mpi_context | 11 | PASS |

## 10. Next Steps

1. **Submit v38b BR2 run on TACC**: Same as v37 config but with shared-face
   Dirichlet fix. Compare recurrence with v37's 295 yr.

2. **Investigate IP blowup**: The IP instability during coseismic needs deeper
   analysis. Options:
   - Compare IP vs BR2 penalty magnitudes during first earthquake
   - Check if IP bilinear form penalty matches traction penalty exactly
   - Try reducing IP penalty by a factor (e.g., 0.5x) to test sensitivity
   - Check Tandem's time step control during coseismic

3. **Eliminate IP bilinear form as cause**: Run with `--dg-method IP` but
   using BR2 bilinear form + IP traction (hybrid approach) to isolate
   whether the blowup comes from the bilinear form or traction computation.

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
| **v38b** | **Shared-face Dirichlet loading fix (parallel bug)** | **Applied, tests pass** |
