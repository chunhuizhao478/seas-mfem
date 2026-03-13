# BP5 Debug v20: Fix Initial State Bug & Medium Run

**Date**: 2026-03-13
**Status**: Implementation
**Previous**: v19 (MUMPS BLR, SuperLU, STRUMPACK solver options)

---

## 1. Problem: Nucleation Delay (192s vs Tandem 48.6s)

MUMPS BLR short run (job 7598114, 1000m mesh) nucleates at **192s** — 4x slower than Tandem (48.6s) despite having a **stronger** perturbation (V_nuc=0.03 + delta_tau vs V_nuc=0.01).

### Root Cause: InitialStatePsi Absorbs delta_tau

**SCEC BP5 Specification**:
- **Eq. 18**: θ(x₂,x₃,0) = L/V_init — initial state at plate rate **everywhere**, including nucleation zone
- **Eq. 23**: τ₀ᵢ = σₙ·a·asinh[Vᵢ/(2V₀)·exp(ψ_ss/a)] + ηVᵢ + δτ, "while keeping the initial state variable θ(x₂,x₃,0) **unchanged**"
- δτ is a **genuine overstress** above friction equilibrium — the state variable is NOT adjusted

**Bug in `fault/rate_state_fault.hpp` (BP5 vector path, lines 323-326)**:
```cpp
// BUG: InitialStatePsi inverts tau_abs (which includes delta_tau)
// to compute psi0, absorbing delta_tau into a higher psi
real_t psi0 = dr_friction_->InitialStatePsi(
   tau_abs, V_abs_init, sigma_n_bp5_, eta, a);
state(i * StatePerNode + PsiIndex) = psi0;
```

**What happens**: `InitialStatePsi()` computes `tau_eff = tau_abs - eta*V_init`, then inverts friction to find psi. Since `tau_abs` includes `delta_tau`, the resulting `psi0` is HIGHER than SCEC-specified, making friction stronger. The system starts in equilibrium at V=V_nuc rather than being overstressed. V_nuc decays → slow natural nucleation at 192s.

---

## 2. Fix: Use SCEC-Specified psi(0) Directly (H24)

**File**: `fault/rate_state_fault.hpp`, BP5 vector path

**Before** (buggy):
```cpp
real_t psi0 = dr_friction_->InitialStatePsi(
   tau_abs, V_abs_init, sigma_n_bp5_, eta, a);
state(i * StatePerNode + PsiIndex) = psi0;
```

**After** (fixed):
```cpp
// SCEC Eq. 18: psi(0) = f0 + b*ln(V0/V_init) everywhere
// delta_tau is genuine overstress, not absorbed into state
real_t psi0 = bp5_params_.psi_init();
state(i * StatePerNode + PsiIndex) = psi0;
```

`bp5_params_.psi_init()` returns `f0 + b*log(V0/V_init)` — exactly SCEC Eq. 18.

### Effect of fix:
- **Outside nucleation zone**: tau0 has no delta_tau. `InitialStatePsi(tau0, V_init, ...)` would recover psi_init anyway. Using psi_init directly is equivalent. **No change in behavior.**
- **Inside nucleation zone**: tau0 includes delta_tau. With psi = psi_init (unchanged), friction is at steady-state for plate rate, but stress is higher by delta_tau. SolveSlipRateVectorPsi will solve for V > V_nuc. **Genuine overstress drives immediate acceleration** instead of equilibrium decay.

---

## 3. Local Verification

8-rank local test confirms the fix works:
- **Initial V_max = 0.0493** — above V_nuc (0.03), overstress drives immediate acceleration
- **V_max increases** in first ~12 steps (0.0495 → 0.0580) — earthquake starts immediately
- **Earthquake #1 detected at t=0** — no 192s delay

Before fix: V_max started at V_nuc=0.03 and decayed (overstress absorbed into psi).
After fix: V_max starts above V_nuc and accelerates (genuine overstress preserved).

---

## 4. Cluster Run: Medium MUMPS BLR

**Sbatch**: `jobs/bp5/bp5_v20_med_mumps_blr.sbatch`
- Solver: MUMPS BLR (fastest for this problem)
- tfinal: 9.45e9 s (~300 years) — goal: capture second seismic event
- 8 nodes, 400 ranks, 1000m mesh

---

## 5. Cumulative Fix History

| Fix | Description | Status |
|-----|-------------|--------|
| H1-H12 | Friction, output, penalty, BC, nucleation fixes | Done |
| H13 | MUMPS print level 0→1 | Done (v17) |
| H14 | Post-solve residual check | Done (v17), fixed in v18 (H18) |
| H15 | AMG SetSystemsOptions(3) → SetElasticityOptions | Done (v17), upgraded (v18) |
| H16 | Runtime --solver flag (default: cg) | Done (v17) |
| H17 | Traction monitoring at fault stations | Done (v17) |
| H18 | Global residual norms via MPI_Allreduce | Done (v18) |
| H19 | CG RelTol 1e-12 → 1e-10 | Done (v18) |
| H20 | GMRES+BlockILU solver option | Done (v18) |
| H21 | MUMPS BLR, SuperLU, STRUMPACK solver options | Done (v19) |
| H22 | GMRES+BlockILU solver option, faster than CG+AMG | Done (v19) |
| H23 | MUMPS BLR set as default solver | Done (v19) |
| **H24** | **Fix InitialStatePsi absorbing delta_tau → use psi_init()** | **Done (v20)** |
