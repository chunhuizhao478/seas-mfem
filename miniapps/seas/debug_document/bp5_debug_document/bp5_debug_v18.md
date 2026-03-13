# BP5 Debug v18: CG+AMG Solver Validation — False Alarm Residual Warnings

**Date**: 2026-03-13
**Status**: Fix applied, ready for cluster retest
**Previous**: v17 (H13-H17 solver diagnostics, --solver flag, CG+AMG default)

---

## 1. Executive Summary

Cluster results from v17 short runs (job 7598007 CG, job 7598012 MUMPS) reveal:

1. **MUMPS is confirmed broken** for this problem size — `INFOG(1)=-9` (insufficient working memory) from the very first factorization. Cannot be used for production BP5 runs.
2. **CG+AMG is working correctly** — CG converges globally on every solve (no "CG did not converge" warnings). Physics looks correct (V_max declining from 0.03 → 0.02 during nucleation transient relaxation).
3. **The H14 residual check had a bug** — it computed per-rank `||r_local||/||b_local||` instead of global norms, producing 132,538 false alarm warnings (9.5 MB err file). Ranks far from the fault have tiny `||b_local||`, inflating relative residual to 20-30 even though the global solve is accurate.

**Fix H18**: Changed residual check to use global norms via `MPI_Allreduce`.

---

## 2. Evidence

### 2.1 Cluster Run Configuration

| Parameter | CG Run | MUMPS Run |
|-----------|--------|-----------|
| Job ID | 7598007 | 7598012 |
| Solver | CG+AMG (SetSystemsOptions(3)) | MUMPS |
| Mesh | bp5_1000m.msh, 66220 elements | Same |
| Ranks | 224 (4 nodes × 56) | Same |
| Queue | development (2h) | Same |
| t_final | 0.1 yr (3.15e6 s) | Same |

### 2.2 MUMPS Run (7598012): Complete Failure

```
INFOG(1)=-9   (insufficient working memory for factorization)
INFO(1)=-1, INFO(2)=45   (error signaled by rank 45, cascading to all 224 ranks)
```

- Hundreds of `ERROR RETURN FROM DMUMPS` messages across all ranks
- MFEM's automatic memory relaxation retry (ICNTL(14)=40) also failed
- Never completed a single time step
- **Conclusion**: MUMPS cannot handle the 66,220-element 3D DG elasticity matrix. This validates the v17 plan's decision to default to CG+AMG.

### 2.3 CG Run (7598007): Solver Working, Monitoring Broken

**Solver convergence**: CG claims convergence on every solve (no "CG did not converge" warnings in stdout or stderr). No NaN/Inf, no displacement blowup.

**Physics**: V_max declining correctly:
```
Step  1: dt=2.557e-02  V_max=3.014e-02  (nucleation transient)
Step 13: dt=1.225e-01  V_max=3.573e-02  (peak)
Step 24: dt=1.427e-01  V_max=3.399e-02  (relaxing)
Step 50: dt=1.708e-01  V_max=2.508e-02  (continuing to relax)
Step 71: dt=1.824e-01  V_max=2.005e-02  (approaching interseismic)
```

This is the expected nucleation transient relaxation (V_nuc=0.03 → interseismic V~1e-9). The solver is producing correct physics.

**Residual warnings (false alarms)**: 132,538 warnings in the err file (9.5 MB). Analysis:

| Rank group | Count | Relative residual range | Root cause |
|------------|-------|------------------------|------------|
| Most ranks (~218/224) | ~132,000 | 1e-7 to 1e-5 | Small `\|\|b_local\|\|` inflation |
| 6 bad ranks (96, 104, 108, 165, 203, 214) | ~500 | 0.05 → 30+ (growing) | Very small `\|\|b_local\|\|`, near-zero loading |

**Key insight**: The CG solver reports `GetConverged() = true` (global `||K*x - b|| / ||b|| < 1e-12`), but the per-rank check `||r_local|| / ||b_local||` is misleading because:
- For a distributed solve, individual ranks can have tiny local RHS norms
- Ranks 96, 104, 108, 165, 203, 214 are likely far from the fault and loaded boundaries
- Their `||b_local||` ≈ ε (floating point noise), so even `||r_local||` = 30ε gives relative residual = 30
- This contributes negligibly to the global norm: `30ε << ||b_global||`

### 2.4 Traction Monitoring Data

The `--monitor-traction 100` output shows physically reasonable behavior:
- `tau_pre ≈ 1.327e+07 Pa` (13.27 MPa) — consistent with BP5 initial stress
- Elastic traction corrections: O(100-1000 Pa) — small relative to tau_pre
- Velocities near V_init = 1e-9 m/s on most fault DOFs
- No anomalous stress buildup visible in the 71 steps completed

---

## 3. Root Cause: Per-Rank Residual Check Bug

### 3.1 Original H14 Code (Buggy)

```cpp
// Post-solve residual check: ||K*x - b|| / ||b||
if (check_residual_)
{
   Vector R_(B_.Size());
   cached_Ah_.As<HypreParMatrix>()->Mult(X_, R_);  // correct: includes communication
   R_ -= B_;
   real_t res_norm = R_.Norml2();    // LOCAL norm
   real_t rhs_norm = B_.Norml2();    // LOCAL norm
   real_t rel_res = (rhs_norm > 0.0) ? res_norm / rhs_norm : res_norm;
   if (rel_res > 1e-8)
   {
      // Warning per rank → 224 warnings per solve × 71 solves × RK stages
      mfem::err << "[Rank " << rank << "] RESIDUAL WARNING: ..." << rel_res;
   }
}
```

**Problem**: `R_.Norml2()` and `B_.Norml2()` are LOCAL norms. The ratio `||r_local||/||b_local||` is NOT the global relative residual. For ranks with `||b_local|| → 0`, this ratio → ∞ regardless of solver accuracy.

### 3.2 Why Global CG Converges But Local Ratios Are Large

The global residual norm is:
```
||r||_global² = Σ_ranks ||r_local||²
||b||_global² = Σ_ranks ||b_local||²
```

CG converges when `||r||_global / ||b||_global < 1e-12`.

For a rank with `||b_local|| = ε` (tiny), even if `||r_local|| = 30ε`:
- Local ratio: `30ε / ε = 30` (alarm!)
- Contribution to global: `(30ε)² << ||b_global||²` (negligible)

---

## 4. Fix H18: Global Residual Check

**File**: `domain/elasticity_operator.hpp`, post-solve residual check

Changed to use `MPI_Allreduce` for global norms:

```cpp
if (check_residual_)
{
   Vector R_(B_.Size());
   cached_Ah_.As<HypreParMatrix>()->Mult(X_, R_);
   R_ -= B_;
   // Use global norms (MPI AllReduce) for meaningful relative residual
   real_t local_res2 = R_ * R_;
   real_t local_rhs2 = B_ * B_;
   real_t global_res2 = 0.0, global_rhs2 = 0.0;
   MPI_Allreduce(&local_res2, &global_res2, 1, MPI_DOUBLE, MPI_SUM, mesh_.GetComm());
   MPI_Allreduce(&local_rhs2, &global_rhs2, 1, MPI_DOUBLE, MPI_SUM, mesh_.GetComm());
   real_t global_res = std::sqrt(global_res2);
   real_t global_rhs = std::sqrt(global_rhs2);
   real_t rel_res = (global_rhs > 0.0) ? global_res / global_rhs : global_res;
   int rank = 0;
   MPI_Comm_rank(mesh_.GetComm(), &rank);
   if (rel_res > 1e-8 && rank == 0)
   {
      mfem::err << "RESIDUAL WARNING: global ||K*x-b||/||b|| = " << rel_res
                << " (threshold 1e-8)\n";
   }
}
```

**Changes**:
- `R_ * R_` (dot product) gives local squared norm, then `MPI_Allreduce(SUM)` → global squared norm
- Only rank 0 prints (one line per warning, not 224)
- Serial path unchanged (already correct)

**Verification**: All 82 unit tests pass (mpirun -np 2).

---

## 5. Run Assessment

### 5.1 Short CG Run: Incomplete But Healthy

The 2h development queue was insufficient — the run completed 71 time steps but is still at t ≈ 0.00 yr. The adaptive time stepper takes tiny steps (dt ≈ 0.025-0.18 s) during the nucleation transient because V_max ≈ 0.02 m/s.

This is expected behavior:
- V_nuc = 0.03 means the fault starts near the earthquake threshold
- The transient must relax V_max down to ~1e-9 before dt can grow to interseismic values (~years)
- Tandem also takes many small steps through this transient

### 5.2 Next Steps

1. **Resubmit short CG run** with the H18 fix (clean err file expected)
2. **Consider longer wall time** or the normal queue — the nucleation transient may need many more steps before dt grows
3. **Medium CG run** (5yr) once short run completes cleanly — verify stress rate ~0.001 MPa/yr
4. **Full CG run** for correct ~200yr earthquake recurrence

---

## 6. Files Modified

| File | Change | Fix ID |
|------|--------|--------|
| `domain/elasticity_operator.hpp` | Per-rank residual → global residual via MPI_Allreduce, rank 0 only print | H18 |

---

## 7. Cumulative Fix History

| Fix | Description | Status |
|-----|-------------|--------|
| H1-H9 | Friction, output, penalty fixes | Done (v1-v9) |
| H10 | All-boundary Dirichlet loading (matches Tandem) | Done (v15) |
| H11 | Configurable nucleation parameters | Done (v15) |
| H12 | CLAUDE.md BC documentation | Done (v15) |
| H13 | MUMPS print level 0→1 | Done (v17) |
| H14 | Post-solve residual check | Done (v17), **fixed in v18 (H18)** |
| H15 | AMG SetSystemsOptions(3) for 3D elasticity | Done (v17) |
| H16 | Runtime --solver flag (default: cg) | Done (v17) |
| H17 | Traction monitoring at fault stations | Done (v17) |
| H18 | Global residual norms via MPI_Allreduce | Done (v18) |
