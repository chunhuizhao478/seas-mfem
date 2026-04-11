# BP5 Debug v17: Spurious Early Second Seismic Event — Solver Diagnostics

**Date**: 2026-03-13
**Status**: Analysis complete, fixes designed
**Previous**: v15 (H10 all-boundary Dirichlet, H11 configurable nucleation)

---

## 1. Executive Summary

Post-v15 results (all-boundary Dirichlet loading matching Tandem) on the cluster (224 ranks, MUMPS, 1000m mesh, 66220 elements) show a **spurious second seismic event at ~2.8 years**. Expected earthquake recurrence is ~200 years per SCEC BP5-QD benchmarks.

Station data reveals interseismic stress rebuilds at **~6 MPa/yr** vs expected **~0.001 MPa/yr** — a 10,000x discrepancy. The boundary loading has been verified correct (matches Tandem exactly). The issue is in the linear solver.

**Root cause**: MUMPS direct solver reports `INFOG(1)=-9` (insufficient working memory) during factorization. The stiffness matrix factorization may be inaccurate, producing wrong displacements → wrong traction → anomalous stress evolution.

---

## 2. Evidence

### 2.1 Cluster Run Configuration

- **Cluster output**: `/Users/chunhuizhao/Downloads/seas-mfem/results_1000m_bp5/bp5_1000m_7597783.out`
- **224 MPI ranks**, MUMPS direct solver
- **66220 tetrahedral elements**, 1000m mesh resolution
- MUMPS error: `INFOG(1)=-9, INFOG(2)=116701`

### 2.2 Station Data Analysis

| Station | After EQ#1 τ (MPa) | τ at ~2.8yr (MPa) | Loading Rate | Expected |
|---------|--------------------|--------------------|-------------|----------|
| (0km, 10km) VW center | 9.14 | 18.75 | ~6 MPa/yr | ~0.001 MPa/yr |
| (-24km, 10km) nuc zone | 9.0 | 17.0 | ~6 MPa/yr | ~0.001 MPa/yr |
| (0km, 0km) surface | ~9 | ~18 | ~6 MPa/yr | ~0.001 MPa/yr |

**Key observations:**
- EQ#1 occurs at t=0 (expected — nucleation zone starts at V_nuc=0.03)
- After EQ#1, fault locks properly (V ≈ 10^-33 m/s)
- Stress rebuilds 10,000x too fast → triggers spurious EQ#2 at ~2.8yr
- All stations show the same anomalous loading rate → systematic solver error

### 2.3 MUMPS Error Details

From the cluster log:
```
MUMPS: INFOG(1)=-9, INFOG(2)=116701
```

- `INFOG(1)=-9`: Insufficient working memory for numerical factorization
- `INFOG(2)=116701`: Additional memory (in MB) that would be needed
- MFEM's MUMPS wrapper has auto-retry logic that increases `ICNTL(14)` (working memory percentage), but our `SetPrintLevel(0)` suppresses all retry/warning messages — we cannot confirm whether retries succeeded

### 2.4 Why the Solver Matters

In our SEAS code, the stiffness matrix K is:
1. Assembled once (`AssembleStiffness()`, `elasticity_operator.hpp:487-563`)
2. Factored once (MUMPS) or preconditioned once (AMG)
3. **Reused for every RK45 stage solve** — hundreds of thousands of calls to `Solve()`

If the factorization is inaccurate:
- Every `solver_->Mult(B_, X_)` produces slightly wrong displacement
- `ComputeTraction()` extracts wrong stress on the fault
- The error accumulates coherently (always in the same direction) because the same bad factor is reused
- Result: systematic stress loading bias → 10,000x too fast

---

## 3. Root Cause Analysis

### 3.1 Current Solver Setup

**File: `domain/elasticity_operator.hpp` lines 526-551**

```cpp
// MUMPS path (current default when compiled with MFEM_USE_MUMPS)
auto *mumps = new MUMPSSolver(mesh_.GetComm());
mumps->SetMatrixSymType(MUMPSSolver::MatType::SYMMETRIC_POSITIVE_DEFINITE);
mumps->SetPrintLevel(0);   // ← PROBLEM: suppresses all INFOG warnings
mumps->SetOperator(*cached_Ah_.As<HypreParMatrix>());

// CG+AMG path (fallback)
auto *cg = new CGSolver(mesh_.GetComm());
cg->SetRelTol(1e-12);
cg->SetMaxIter(10000);
cg->SetPrintLevel(0);
auto *amg = new HypreBoomerAMG(*cached_Ah_.As<HypreParMatrix>());
amg->SetPrintLevel(0);
```

### 3.2 Solve Method

**File: `domain/elasticity_operator.hpp` lines 1699-1808**

```cpp
X_ = 0.0;           // zero initial guess
B_ = rhs;
solver_->Mult(B_, X_);  // solve K*X = B
```

The solve already checks for NaN/Inf and CG convergence, but does **not** check:
- MUMPS factorization quality (INFOG status)
- Post-solve residual `||K*x - b|| / ||b||`

### 3.3 Why CG+AMG Wasn't Used

The code defaults to MUMPS when `MFEM_USE_MUMPS` is defined (compile-time). The cluster build has MUMPS enabled. There is no runtime `--solver` flag to select CG+AMG.

### 3.4 CG+AMG Missing Elasticity-Aware AMG

The current AMG setup uses default BoomerAMG settings. For 3D elasticity, `SetSystemsOptions(dim)` enables rigid-body-mode-aware coarsening, which significantly improves convergence. Without it, CG+AMG may require many more iterations for 3D elasticity.

---

## 4. Proposed Fixes

### H13: MUMPS Print Level (Quick Win)
**File**: `domain/elasticity_operator.hpp` line 531

Change:
```cpp
mumps->SetPrintLevel(0);  →  mumps->SetPrintLevel(1);
```

This reveals MUMPS retry messages, INFOG warnings, and factorization statistics. No performance impact.

### H14: Post-Solve Residual Check
**File**: `domain/elasticity_operator.hpp`, after `solver_->Mult(B_, X_)` (line 1753)

Add:
```cpp
// Post-solve residual check
if (check_residual_)
{
   Vector r(B_.Size());
   cached_Ah_.As<HypreParMatrix>()->Mult(X_, r);  // r = K*x
   r -= B_;                                         // r = K*x - b
   double res_norm = r.Norml2();
   double rhs_norm = B_.Norml2();
   double rel_residual = (rhs_norm > 0) ? res_norm / rhs_norm : res_norm;
   if (rel_residual > 1e-8)
   {
      mfem::out << "WARNING: Solve relative residual = " << rel_residual
                << " (threshold 1e-8)" << std::endl;
   }
}
```

This catches both MUMPS factorization errors and CG convergence failures. The `check_residual_` flag is controlled by a CLI option.

### H15: CG+AMG with Elasticity-Aware AMG Settings
**File**: `domain/elasticity_operator.hpp` line 544-546

Add `SetSystemsOptions(dim)` for rigid-body-mode-aware AMG:
```cpp
auto *amg = new HypreBoomerAMG(*cached_Ah_.As<HypreParMatrix>());
amg->SetSystemsOptions(3);   // 3D elasticity: near-null-space = 6 rigid body modes
amg->SetPrintLevel(0);
```

This is critical for 3D elasticity — default BoomerAMG treats the system as scalar, missing the block structure.

### H16: Runtime Solver Selection (`--solver` CLI Flag)
**File**: `tests/verification/bp5_verification_full.cpp`

Add:
```
--solver mumps|cg    Select linear solver (default: cg)
```

Default to CG+AMG for robustness. MUMPS remains available for small problems or debugging.

Wire to `ElasticityDomainOperator` constructor via a `use_mumps` parameter (already exists as `use_mumps_` member).

### H17: Traction Monitoring at Key Stations
**File**: `fault/rate_state_fault.hpp`, `ComputeRHS()` method

Log `tau_pre`, elastic `traction`, and total `tau` at 2-3 fault DOFs nearest benchmark stations every N steps:
```
Step 100: station(0,10km): tau_pre=25.2 traction=-16.1 tau=9.1 V=1e-33
```

Controlled by `--monitor-traction` CLI flag. Low overhead (3 DOF lookups per step).

---

## 5. Implementation Order

| Priority | Fix | Rationale |
|----------|-----|-----------|
| 1 | H14 (residual check) | Most diagnostic value — confirms solver accuracy |
| 2 | H15 (elasticity-aware AMG) | Resolves the root cause for CG path |
| 3 | H16 (--solver flag, default cg) | Makes CG+AMG the default, avoids MUMPS memory issue |
| 4 | H13 (MUMPS print level) | Quick, helps debug future MUMPS runs |
| 5 | H17 (traction monitoring) | Lower priority, useful for long runs |

---

## 6. Files to Modify

| File | Changes |
|------|---------|
| `domain/elasticity_operator.hpp` | H13: MUMPS print level 0→1; H14: post-solve residual check; H15: AMG `SetSystemsOptions(3)` |
| `tests/verification/bp5_verification_full.cpp` | H16: `--solver` and `--check-residual` CLI flags |
| `fault/rate_state_fault.hpp` | H17: traction monitoring (optional) |

---

## 7. Verification

### Step 1: Local test (8 ranks, 1000m mesh, 1 step)
```bash
conda activate mfem-dev
cd /Users/chunhuizhao/projects/seas-mfem/miniapps/seas
make -j seas_bp5_full
mpirun -np 8 ./seas_bp5_full --mesh bp5/mesh/bp5_v2_1000m.msh \
   --mesh-scale 1000 --tfinal 1e-9 --solver cg --check-residual \
   --output-dir bp5/diag_output/v17_test
```
**Expected**: Residual < 1e-8 with both MUMPS and CG solvers.

### Step 2: Short run (~0.1 yr)
```bash
mpirun -np 8 ./seas_bp5_full --mesh bp5/mesh/bp5_v2_1000m.msh \
   --mesh-scale 1000 --tfinal 3.15e6 --solver cg --check-residual \
   --output-dir bp5/diag_output/v17_short
```
**Expected**: No residual warnings; traction at (0,10km) ≈ tau_pre + small elastic correction.

### Step 3: Medium run (~5 yr) with CG+AMG
```bash
mpirun -np 8 ./seas_bp5_full --mesh bp5/mesh/bp5_v2_1000m.msh \
   --mesh-scale 1000 --tfinal 1.58e8 --solver cg --check-residual \
   --output-dir bp5/diag_output/v17_5yr
```
**Expected**: Stress loading rate at (0, 10km) ≈ 0.001 MPa/yr (not 6 MPa/yr). No spurious second event.

### Step 4: Cluster rerun
Submit with CG+AMG default + residual check on full 1000m mesh. Verify no spurious second event through at least 10yr.

---

## 8. Previously Implemented Fixes

- H1-H7: Various friction, output, penalty fixes
- H8: Bottom-only Dirichlet BC (superseded by H10)
- H9: Output floor removed
- H10: All-boundary Dirichlet loading (matches Tandem)
- H11: Configurable nucleation parameters
- H12: Updated CLAUDE.md BC documentation
