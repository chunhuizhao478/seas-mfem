# BP5 Debug v5: Persistent Blowup After Shared Face Fix

## Problem

After adding shared face support to `ElasticityDomainOperator` (bp5_debug_v4), BP5 still blows up on Frontera (56 MPI ranks):
```
zeroIn: F(a) and F(b) must have different signs.
Rank 25: F(a)=2.1642750689123562e+10  b=4680
Rank 34: F(a)=3.3719270683915901e+10  b=7291
```
The error DECREASED from v3 (F(a)=2.8e10) to v4 (F(a)=2.16e10), confirming the shared face fix helped but didn't fully resolve the issue.

F(a) = tau_abs ≈ 21.6 GPa. Pre-stress ≈ 15 MPa. So elastic traction ≈ 21.6 GPa. With mu=32 GPa, this requires strain ~0.34 — impossible physically. The displacement solution is wildly wrong at certain DOFs.

## Analysis Summary

### What was ruled out:
1. **`scalar_fes_` being serial** — Only used for `GetFE(Elem1No)` (local elements); not a bug
2. **Mass inverse indexing** — `elem_mass_inv_[FTr->Elem2No]` is correct because `PrecomputeMassInverse` stores at `ne + i` and `FTr->Elem2No = ne + i` for shared faces
3. **Missing ExchangeFaceNbrData** — Called in `PrecomputeMassInverse` and `AssembleStiffness` before any shared face access; not needed again in slip assembly (only mass inverses and shapes needed, not DOF values)
4. **Double-counting shared faces** — Interior faces and shared faces are disjoint sets in MFEM
5. **Sign convention in shared face assembly** — Matches interior face pattern and antiplane reference

### Bug found: Normal convention inconsistency in ComputeTraction BR2 correction

In `ComputeTraction()` BR2 penalty correction (lines 1909-1960 for interior, lines 2143-2180 for shared):

- **face_int** uses `nor` (unnormalized, from `CalcOrtho`) — CORRECT for face integral
- **Elasticity tensor coupling** uses `basis.normal` (UNIT normal) — INCONSISTENT

Compare with:
- **BR2 stiffness matrix** (`dg_elasticity_br2_integrator.hpp`): `TestNormal` uses unnormalized `n_q` in BOTH face_int AND tensor coupling
- **BR2 RHS slip assembly** (`AssembleSlipContributionBR2`): uses unnormalized `n_q` in BOTH face_int AND tensor coupling (lines 854-856)

The traction correction should use the same convention for consistency. The correct formula for point evaluation of the BR2 traction correction is `T -= η * C : r_e([[u]] - δ) · n_hat` where:
- r_e uses the full face integral (unnormalized normal)
- The final traction uses the unit normal n_hat

However, the face_int in the traction code also OMITS the quadrature weight `w_q` that the matrix and RHS include. This means the face_int represents an "integrand" not an "integral". For consistency with the unit normal in the tensor coupling, this is actually a self-consistent (but different) convention — it computes the BR2 correction density rather than the integral.

**Impact estimate**: For hex elements with h=1000m, the correction magnitude is ~288 MPa per meter of slip — significant but NOT 21 GPa. This bug alone cannot explain the blowup.

### Most likely remaining cause: Issue only reproducible at scale

The error occurs only on Frontera (56 ranks) but not locally (2-4 ranks). Possible explanations:
1. **MUMPS solver issue at scale** — factorization failure returning garbage silently
2. **Mesh partitioning edge case** — specific partitioning configurations that don't appear with few ranks
3. **Undetected fault faces** — faces at fault zone boundaries with floating-point coordinate issues
4. **CG solver (non-MUMPS build)** — if Frontera build uses CG+ILU instead of MUMPS, convergence may fail silently

## Fix

### Phase 1: Add diagnostic output and local reproduction test

#### 1.1 Add diagnostic traction bounds check in `ComputeTraction()`
**File**: `miniapps/seas/domain/elasticity_operator.hpp`

After computing all traction values, add bounds check that prints rank, DOF index, face type (interior/shared), and traction magnitude for any DOF with |tau| > 1 GPa.

#### 1.2 Add solver convergence check
**File**: `miniapps/seas/domain/elasticity_operator.hpp`, in `Solve()`

After `solver_->Mult(B_, X_)`, check `||u||_inf` and `||RHS||_inf` for blowup indicators.

#### 1.3 Add RHS assembly diagnostics
Check max RHS norm before and after slip contribution to identify if the blowup comes from slip assembly or Dirichlet loading.

### Phase 2: Fix BR2 traction correction normal convention

**File**: `miniapps/seas/domain/elasticity_operator.hpp`

In `ComputeTraction()` BR2 correction, change the elasticity tensor coupling to use unnormalized `nor` instead of `basis.normal`, then normalize the entire correction by `||nor||^2`:

For **interior faces** (line ~1937-1960):
```cpp
// BEFORE (inconsistent):
real_t tn = lambda_val_ * (u == s ? 1.0 : 0.0) * basis.normal[i]
   + mu_val_ * ((i == u ? 1.0 : 0.0) * basis.normal[s]
               + (i == s ? 1.0 : 0.0) * basis.normal[u]);

// AFTER (consistent with matrix and RHS):
real_t tn = lambda_val_ * (u == s ? 1.0 : 0.0) * nor(i)
   + mu_val_ * ((i == u ? 1.0 : 0.0) * nor(s)
               + (i == s ? 1.0 : 0.0) * nor(u));
```

Then normalize the correction:
```cpp
real_t nor_sq = nor * nor;
for (int c = 0; c < dim; c++) {
   correction[c] = br2_penalty * 0.5 * sum_c / nor_sq;
}
```

Same change for **shared faces** (line ~2128-2180).

**Rationale**: The BR2 lifting face_int uses `nor` (unnormalized). The tensor coupling should also use `nor` for consistency (as in the matrix and RHS). Then divide by `||nor||^2` to convert from force to traction (force/area).

### Phase 3: Write BP5 parallel smoke test

**New file**: `miniapps/seas/tests/parallel/test_bp5_parallel_smoke.cpp`

This test runs the ACTUAL BP5 simulation for 2-3 time steps with a coarse mesh (dx~4000m) and 2-4 ranks. It should detect the blowup that unit tests miss.

Key checks:
- All traction values bounded (|tau| < 1 GPa)
- No NaN/Inf in displacement
- zeroIn converges for all fault DOFs
- Slip rate bounded (|V| < 1e4 m/s)

### Phase 4: Run and diagnose

1. Build and run the smoke test with `mpirun -np 2` and `mpirun -np 4`
2. If error reproduces locally: use diagnostics from Phase 1 to identify problematic DOFs
3. If error doesn't reproduce: increase to 8-16 ranks or use a finer mesh
4. Analyze diagnostic output to determine root cause

## Files Modified

1. **`miniapps/seas/domain/elasticity_operator.hpp`** — Add diagnostics (1.1-1.3), fix BR2 correction normals (Phase 2)
2. **`miniapps/seas/tests/parallel/test_bp5_parallel_smoke.cpp`** — New parallel smoke test (Phase 3)
3. **`miniapps/seas/Makefile`** — Add build target for smoke test

## Verification

1. `conda activate mfem-dev && make -j4` — build
2. `./seas_test_elasticity_operator` — all serial tests still pass
3. `mpirun -np 2 ./seas_test_parallel_elasticity` — all 6 parallel tests still pass
4. `mpirun -np 4 ./seas_test_parallel_elasticity` — all parallel tests still pass
5. `mpirun -np 2 ./seas_test_bp5_parallel_smoke` — BP5 smoke test passes (no blowup)
6. `mpirun -np 4 ./seas_test_bp5_parallel_smoke` — BP5 smoke test passes
7. If smoke test passes locally, push to Frontera and test with 56 ranks
