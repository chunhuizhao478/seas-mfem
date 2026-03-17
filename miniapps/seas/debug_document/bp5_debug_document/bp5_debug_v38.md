# BP5 Debug v38: Switch to IP Method (Matching Tandem)

**Date**: 2026-03-16
**Status**: Code ready, TACC submission pending
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

## 3. Why IP Should Work Now

**v31 (DG slip sign fix)** identified the root cause of ALL v30 blowups:

> v30 negated `tau_pre` to match Tandem convention, making V < 0. But the DG
> sign variable was never adjusted: `sign = (nor(1) > 0) ? -1 : +1`.
> With sign=-1 and delta_u<0: g = sign*delta_u = +|delta_u| — **wrong sign**.
> → penalty corrects in wrong direction → **exponential blowup**

This sign bug existed in v24-v28 (it was introduced in v30's tau_pre negation,
but the coordinate system change in v30 was BEFORE the IP experiments in
v26-v28). After v31's fix, BR2 works perfectly. **IP was never retested.**

v31 explicitly says:
```
// If nor(1) > 0: Elem1 is -Y side, [[u]] = u(-Y) - u(+Y) = delta_u -> sign = +1
// If nor(1) < 0: Elem1 is +Y side, [[u]] = u(+Y) - u(-Y) = -delta_u -> sign = -1
real_t sign = (nor(1) > 0) ? 1.0 : -1.0;
```

The sign fix applies to ALL 6 locations (both IP and BR2 paths).

## 4. Current IP Code Path (No Changes Needed)

The IP code path already includes all necessary fixes:

| Component | Implementation | Status |
|-----------|---------------|--------|
| Bilinear form | Custom `DGElasticityIPPenaltyIntegrator` (v28): `penalty * \|nor\|` | Ready |
| Slip RHS | `AssembleSlipContributionIP`: Tandem's scalar penalty formula | Ready |
| Dirichlet RHS (boundary) | `AssembleDirichletLoading` IP path: `penalty * \|nor\|` | Ready |
| Dirichlet RHS (interior Y=0) | `AssembleDirichletLoading` interior IP: same formula | Ready |
| Traction | `ComputeTraction` IP: scalar penalty, per-quad-point (v35-era) | Ready |
| Sign convention | v31 fix: all 6 locations | Ready |
| Fault detection | Tag-based (H32): excludes z=0, z=Wf faces | Ready |
| Interior Dirichlet | v34: Y=0 non-fault faces with jump loading | Ready |

### IP penalty formula (matches Tandem exactly):
```
p(side) = (Dim+1) * c_N_1 * (area/volume) * (c1²/c0)
penalty = (p(0) + p(1)) / 4.0   (interior faces)
penalty = p(0)                    (boundary faces)

c0 = 2μ, c1 = Dim·λ + 2μ
c_N_1 = order * (order + Dim - 1) / Dim
```

### Traction formula (matches Tandem exactly):
```
T = {σ·n̂} - penalty * ([[u]] - sign·δu)   (scalar × vector)
```

## 5. v38a Test Plan

### 5.1 Full Runs (1800 yr each)

| Job | Config | Nodes | Key Question |
|-----|--------|-------|-------------|
| `bp5_v38a_1000m_ip_p1_full` | IP, p=1, 1000m | 8 | Does IP work after v31 sign fix? Recurrence? |
| `bp5_v38a_1000m_br2_p1_full` | BR2, p=1, 1000m | 8 | Baseline comparison (same code, different method) |

### 5.2 Expected Outcomes

**IP p=1**: If IP now works (no blowup):
- Traction should use Tandem's exact formula (scalar penalty)
- No BR2 lifting bias → may match Tandem's 240yr recurrence more closely
- The O(h^p) traction correction bias from BR2 lifting is eliminated

**BR2 p=1 (baseline)**: Should reproduce v37 results (~295 yr) for reference.

### 5.3 Diagnostic Checks

1. **First 10 steps**: If traction stays O(15 MPa), IP is stable
2. **t = 2 yr**: This is where v24-v28 IP blew up — critical checkpoint
3. **V_max time series**: earthquake nucleation and cycling
4. **tau_dip**: should be small and stable (no growing antisymmetric pattern)

## 6. Why IP Should Give Better Results Than BR2

From v32_v33 analysis:
- BR2 traction includes `σ × {{C:L([[u]]-δ)}}·n̂` — a non-trivial O(h^p) correction
- At p=1, h=1000m: correction ~0.06 MPa = 0.4% of fault stress
- Rate-and-state amplifies: `δV/V ≈ (1/a) × δτ/σ_n ≈ 40%` slip-rate bias
- This accumulates over centuries → 1.82× longer recurrence (435yr vs 240yr)

IP traction: `T = {σ·n̂} - penalty × jump` where penalty is properly scaled:
- During interseismic: jump ≈ 0 → correction ≈ 0
- The traction is essentially `{σ·n̂}` (average stress only)
- No O(h^p) bias from BR2 lifting

This is the same approach as Tandem, which achieves 240yr recurrence at p=6.

## 7. Files Changed

**None.** v38a uses the existing IP code path with `--dg-method IP`.

The only change is the sbatch submission scripts.

## 8. Cumulative Fix History

| Fix | Description | Status |
|-----|-------------|--------|
| H1-H32 | Previous fixes (see v1-v28 docs) | Done |
| v30 | Tandem coordinate system | Done |
| v31 | DG slip sign fix + interior Dirichlet | Done |
| v32-v33 | General polynomial order (p-refinement) | Done |
| v34 | Interior Dirichlet Y=0 jump loading | Done |
| v35-v36 | Per-quad-point traction + cross-element BR2 | Done |
| v37 | Interior Dirichlet face_int2 sign fix | Done |
| **v38a** | **Switch to IP method (matching Tandem)** | **TACC pending** |
