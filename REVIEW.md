# Plan Review: Dynamic Rupture Implementation Plan v4 (2026-04-12)

## Review Scope
- Plan: `miniapps/seas/document/system_dev/dynamic_rupture_plan_v4.pdf` (22 pages)
- Codebase analyzed: `miniapps/seas/` (all modules: domain, fault, friction, solver, config, io, common)
- Key files read: `domain/domain_operator.hpp`, `domain/elasticity_operator.hpp`, `solver/seas_operator.hpp`, `fault/rate_state_fault.hpp`, `fault/fault_basis.hpp`, `friction/dieterich_ruina.hpp`, `constitutive/constitutive_model.hpp`, `domain/boundary_config.hpp`, `config/seas_config.hpp`, `solver/time_stepper.hpp`, `Makefile`
- Domain context: CLAUDE.md, CODEBASE_GUIDE.md, ARCHITECTURE.md, Tandem bp5.geo, SeisSol reference
- Cross-reference: `bp5/mesh/bp5_v2.geo` (existing distance field grading), existing sbatch patterns

## Findings

### [R-001] [CRITICAL] [Plan §1.2 + §3.1] — WaveOperator inheritance model conflicts with existing architecture

**Category:** DEVIATION

**Description:**
The plan states: "The WaveOperator inherits `DomainOperator<MeshType>` and implements `TimeDependentOperator::Mult()`." This implies:

```cpp
class WaveOperator : public DomainOperator<MeshType>, public TimeDependentOperator
```

But the existing architecture uses **composition, not dual inheritance**. `SEASQuasiDynamicOperator` inherits `TimeDependentOperator` and OWNS a `DomainOperator*` via a pointer:

```cpp
// solver/seas_operator.hpp:53
class SEASQuasiDynamicOperator : public TimeDependentOperator {
   DomainOperator<ParMesh>* domain_;  // composition
   RateStateFaultOperator<ParMesh,2>* fault_;
};
```

The plan's dual inheritance creates three problems:
1. `DomainOperator::Solve()` is **pure virtual** (designed for quasi-static equilibrium). WaveOperator must implement it, but during FD time-stepping `Solve()` has no meaning — it would be a no-op pure virtual, which is a code smell.
2. `DomainOperator::ComputeTraction()` extracts traction from a displacement-based DG stiffness solve. WaveOperator's traction comes from the Riemann solver at fault faces — a fundamentally different mechanism. Forcing it through the same interface is misleading.
3. The Phase 5 hybrid operator (`seas_hybrid_operator.hpp`) needs to switch between a QD `ElasticityDomainOperator` and FD `WaveOperator`. With composition, this is trivial (swap the pointer). With dual inheritance, the hybrid operator must manage two different type hierarchies.

**Trigger:**
Attempting to implement `WaveOperator : public DomainOperator<MeshType>, public TimeDependentOperator` and discovering that `Solve()` and `ComputeTraction()` semantics don't map to explicit FD time-stepping.

**Actual behavior:**
Plan prescribes dual inheritance.

**Expected behavior:**
Follow the existing composition pattern:
```cpp
// Option A: Follow SEASQuasiDynamicOperator pattern
class WaveOperator : public TimeDependentOperator {
   // OWNS mesh, FE spaces, mass matrix, etc. — no DomainOperator base
   void Mult(const Vector &Q, Vector &dQdt) const override;  // RHS
};

class SEASDynamicOperator : public TimeDependentOperator {
   WaveOperator* wave_;
   RateStateFaultOperator<ParMesh,2>* fault_;
   void Mult(const Vector &state, Vector &dstate_dt) const override;
};
```

```cpp
// Option B: If DomainOperator reuse is desired, add virtual defaults
class DomainOperator {
   virtual void Solve(...) { MFEM_ABORT("Not supported in this mode"); }
   virtual void ComputeTraction(...) { MFEM_ABORT("Not supported"); }
   // New methods for FD:
   virtual void ComputeRHS(const Vector &Q, Vector &dQdt) {}
   virtual real_t GetCFL() const { return 1e30; }
};
```

**Suggested fix:**
Rewrite Plan §1.2 Overview and §3.1 Phase 1 to use composition. The `WaveOperator` should inherit ONLY `TimeDependentOperator` (like `SEASQuasiDynamicOperator` does). Shared fault-coupling code (`FaultBasis`, `FaultGeometry`, etc.) is accessed through composition, not inheritance. Add a brief justification for diverging from the `DomainOperator` base class.

**Test case:**
Verify that the hybrid operator (Phase 5) can switch between QD and FD modes by swapping composed operators, not by type-casting through inheritance.

---

### [R-002] [CRITICAL] [Plan §2.4.7 + §3.6.3] — Brent bracket sign analysis is inverted in documentation

**Category:** BUG

**Description:**
Section 2.4.7 states the Brent bracket guarantee:

> "Bracket: [V_lo = 0, V_hi = Θ/η_s]. Guaranteed sign change because **g(0) = Θ > 0** and **g(Θ/η_s) = −|σ_n|f < 0**."

The residual function is g(V̂) = |σ_n| f(V̂, ψ) + η_s V̂ − Θ (from Eq. 8 rearranged). Evaluating:

- g(0) = |σ_n| f(0, ψ) + 0 − Θ = 0 − Θ = **−Θ < 0** (not +Θ)
- g(Θ/η_s) = |σ_n| f(Θ/η_s, ψ) + Θ − Θ = |σ_n| f(Θ/η_s, ψ) = **+|σ_n|f > 0** (not −|σ_n|f)

The signs are **swapped**. The bracket [0, Θ/η_s] IS correct (sign change exists), but the documented polarity is wrong. If the implementer uses Brent's method with the documented signs (e.g., asserting `g(V_lo) > 0` as a precondition), the assertion will fire and crash.

**Trigger:**
Implementing the Brent solver with precondition checks based on the documented signs.

**Actual behavior:**
Plan states g(0) > 0 and g(V_hi) < 0.

**Expected behavior:**
g(0) = −Θ < 0 and g(V_hi) = |σ_n|f > 0. The lower bracket is negative, upper is positive.

**Suggested fix:**
In §2.4.7, replace:
```diff
- Bracket: [V_lo = 0, V_hi = Θ/η_s]. Guaranteed sign change because g(0) = Θ > 0
- and g(Θ/η_s) = −|σ_n|f < 0.
+ Bracket: [V_lo = 0, V_hi = Θ/η_s]. Guaranteed sign change because g(0) = −Θ < 0
+ and g(Θ/η_s) = +|σ_n|f(Θ/η_s, ψ) > 0.
```

Also add: "This is the opposite polarity from the QD Brent solver (where F(V_lo) > 0 and F(V_hi) < 0) — ensure the Brent implementation handles both orderings."

**Test case:**
Test 34 (`TestBrentNREquivalence`): verify that Brent and NR produce the same result AND that the bracket g(0) < 0 < g(V_hi) holds for 1000 random parameter sets.

---

### [R-003] [MODERATE] [Plan §3.2.4] — SBI DtN kernel formula drops from 2D to 1D without justification

**Category:** BUG

**Description:**
The SBI section presents two inconsistent formulas:

- Eq. (18): τ̂(k_x, k_z) = −μ |k| δ̂(k_x, k_z), where |k| = √(k_x² + k_z²) — **2D Fourier space** (correct for 3D problem)
- Eq. (19): T(z, t) = −μ · F⁻¹[|k_z| û(x_b, k_z, t)] — **1D Fourier space** (only z-wavenumber)

For a 3D simulation with a 2D boundary surface (e.g., the y-z plane at x = x_b), the DtN map requires a **2D** Fourier transform in (y, z) with wavenumber magnitude |k| = √(k_y² + k_z²). Eq. (19) uses only k_z, which would be correct for a 2D simulation but not 3D.

The implementation files (`dynamic/sbi_kernel.hpp/.cpp`) will need 2D FFTs (FFTW r2c 2D plans), not 1D. The ~200 LOC estimate may undercount if this wasn't anticipated.

**Trigger:**
Implementing the SBI kernel for the 3D BP5 domain (boundary is a 2D surface).

**Actual behavior:**
Eq. (19) shows a 1D Fourier transform, implying a 1D FFT implementation.

**Expected behavior:**
Eq. (19) should use 2D Fourier transform:

T(y, z, t) = −μ · F₂D⁻¹[|k| û(x_b, k_y, k_z, t)]

where |k| = √(k_y² + k_z²) and F₂D is the 2D discrete Fourier transform.

**Suggested fix:**
Replace Eq. (19) with the 2D formulation. Update the implementation note to specify FFTW 2D r2c plans (`fftw_plan_dft_r2c_2d`). Add a note that for TPV102 (which is a 3D problem with a 2D fault), the 2D transform is required. Update LOC estimate for `sbi_kernel.cpp` from ~200 to ~300 to account for 2D FFT management.

**Test case:**
Test 15c (`TestSBIZeroReflection`): verify E_res/E_in < 10⁻¹⁰ using the 2D kernel on a 3D mesh (not a 2D slice).

---

### [R-004] [MODERATE] [Plan §3.2.3] — PML corner treatment in code sketch is ambiguous

**Category:** ASSUMPTION

**Description:**
The plan correctly states that at corners where two PML regions overlap: "d(x) **D** = d_x(x) **D**_x + d_y(x) **D**_y + d_z(x) **D**_z." But the implementation code sketch (p. 14) computes the scalar `d = pml_layer_->ComputeDamping(x_q)` and then applies it uniformly to all damped components:

```cpp
real_t d = pml_layer_->ComputeDamping(x_q);
if (d > 0.0) {
    for (int c = 0; c < 9; c++)
        if (D[c] > 0)
            rhs_e[c*ndof+i] -= w * d * shape(i) * Q_qp[c];
}
```

At a corner (e.g., where x-PML and z-PML overlap), `σ_xx` (index 0) should be damped by d_x, `σ_xz` (index 5) should be damped by d_x + d_z (it has both an x-index and z-index), and `v_x` (index 6) should be damped by d_x only. But the code uses a single scalar `d` for all components. The `D[c]` array would need to be a per-component damping value, not a binary flag.

The single-scalar approach is correct when only ONE PML direction is active (non-corner). At corners, each state component needs its own damping coefficient.

**Trigger:**
Simulating a domain with corner PML regions (e.g., where x-boundary and z-boundary meet).

**Actual behavior:**
All damped components get the same d value. At corners, this over-damps some components (e.g., σ_yy should only be damped by d_x for x-PML, not by d_z, but would get d_x + d_z at a corner).

**Expected behavior:**
Per-component damping: `d_c = d_x * D_x[c] + d_y * D_y[c] + d_z * D_z[c]` computed for each component c.

**Suggested fix:**
Replace the code sketch with:
```cpp
real_t dx = pml_layer_->ComputeDamping_x(x_q);
real_t dy = pml_layer_->ComputeDamping_y(x_q);
real_t dz = pml_layer_->ComputeDamping_z(x_q);
// D_total[c] = dx*Dx[c] + dy*Dy[c] + dz*Dz[c]
static const int Dx[] = {1,0,0,1,0,1,1,0,0};
static const int Dy[] = {0,1,0,1,1,0,0,1,0};
static const int Dz[] = {0,0,1,0,1,1,0,0,1};
for (int c = 0; c < 9; c++) {
    real_t d_c = dx*Dx[c] + dy*Dy[c] + dz*Dz[c];
    if (d_c > 0.0)
        for (int i = 0; i < ndof; i++)
            rhs_e[c*ndof+i] -= w * d_c * shape(i) * Q_qp[c];
}
```

Also update `ComputeDamping()` to return 3 directional values, not a single scalar.

**Test case:**
Test 18b (`TestPMLCorner`): specifically checks that a corner PML with two overlapping layers is stable for 500 steps. Add an additional check that the per-component damping is correct by verifying energy decay rates match the theoretical prediction.

---

### [R-005] [MODERATE] [Plan §2.4.4 + §3.3] — Imposed state derivation omits fault-local → global rotation step

**Category:** ASSUMPTION

**Description:**
Section 2.4.4 derives the imposed state (Eqs. 11-12) in **fault-local** coordinates (normal n, tangent t_1, tangent t_2). The code dictionary (p. 11) maps these directly to global-frame indices:

| Math | Code | Meaning |
|------|------|---------|
| v_{t1}^{+,imp} | `imposed_plus[VY]` | Imposed tangent-1 velocity, plus side |
| τ_1^{corr} | `t1_corr` | Corrected tangent-1 traction |

This mapping (`VY` = tangent-1) is only valid when the fault normal aligns with the x-axis (i.e., the fault is a YZ plane). For a general fault orientation in 3D (including BP5's fault at Y=0), the tangent directions don't align with coordinate axes.

The existing codebase handles this via `FaultBasis::ProjectTraction()` (fault-global → fault-local) and `FaultBasis::EmbedSlip()` (fault-local → global). The imposed state must follow the same pattern:

1. Rotate Q^± to fault-local using `FaultBasis`
2. Compute trial traction + friction solve in fault-local frame
3. Construct imposed state Q^{±,imp} in fault-local frame (Eqs. 11-12)
4. **Rotate Q^{±,imp} back to global frame** using `FaultBasis` inverse
5. Feed global-frame Q^{±,imp} to the Godunov flux (Eq. 4)

Step 4 is not mentioned in the plan. The code dictionary's direct mapping to `VY`, `SXY` etc. skips the rotation.

**Trigger:**
Implementing the fault-face flux for BP5, where the fault normal is in the Y-direction, not X.

**Actual behavior:**
Plan maps fault-local variables directly to global indices without rotation.

**Expected behavior:**
Plan should explicitly state that the imposed state computation requires the `FaultBasis` rotation before and after the Riemann solve, consistent with how the existing QD traction extraction works.

**Suggested fix:**
Add a subsection "2.4.6 Implementation: Fault-Local to Global Rotation" that documents:
1. The rotation T (from §2.1.4) must be applied to Q^± before entering the trial traction computation
2. After constructing Q^{±,imp} in fault-local coordinates, apply T⁻¹ to get global-frame values
3. Reference `FaultBasis::ProjectTraction()` and `FaultBasis::EmbedSlip()` as the existing implementation of T and T⁻¹
4. Update the code dictionary to show fault-local indices (0-8) rather than global enums (SXX, VX, etc.)

**Test case:**
Test in Phase 3 (`test_fault_face_flux.cpp`): verify that for a fault at Y=0 (BP5 orientation), the imposed state produces the correct slip rate when the fault normal is (0,1,0), not (1,0,0).

---

### [R-006] [MODERATE] [Plan §3.4.2] — Frontera sbatch template uses conda instead of module loads

**Category:** DEVIATION

**Description:**
The Phase 4b Frontera job script (p. 17) uses:
```bash
module load intel/19 impi/19 phdf5/1.10.4
conda activate mfem-dev
```

But ALL existing production sbatch files in `jobs/bp5/` use explicit module loads without conda:
```bash
module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8
export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:..."
```

`conda activate` in a non-interactive SLURM batch script often fails because conda's shell initialization (`conda init`) hasn't been sourced. The script would need `source ~/.bashrc` or `eval "$(conda shell.bash hook)"` first, which is fragile. The existing pattern with explicit module loads and `LD_LIBRARY_PATH` is reliable.

Additionally, the template is missing: MUMPS, HYPRE, ParMETIS, PETSc module loads that the existing production pattern requires, and the `LD_LIBRARY_PATH` export.

**Trigger:**
Submitting the TPV102 job on Frontera — conda activation fails silently, linking against wrong libraries.

**Actual behavior:**
Template uses `conda activate mfem-dev` which may fail in batch mode.

**Expected behavior:**
Follow the established pattern from `jobs/bp5/bp5_phase7_new_driver_production.sbatch`.

**Suggested fix:**
Replace the Frontera template (p. 17) with a pattern matching the existing sbatch files:
```bash
#!/bin/bash
#SBATCH -J tpv102_seas_mfem
#SBATCH -o tpv102_%j.out
#SBATCH -e tpv102_%j.err
#SBATCH -p normal
#SBATCH -N 4
#SBATCH -n 224
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

export LC_ALL=C
export LANG=C

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8

export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${LD_LIBRARY_PATH}"

cd /scratch2/10024/zhaochun/seas-project/seas-mfem
cd miniapps/seas

test -x ./seas_tpv102_driver || { echo "ERROR: driver not built"; exit 1; }

ibrun ./seas_tpv102_driver \
    --config tpv102/config/tpv102.toml \
    --override output.output_dir="tpv102/results/${SLURM_JOB_ID}"
```

---

### [R-007] [MODERATE] [Plan §3.4.2] — TPV102 driver CLI interface diverges from seas_driver convention

**Category:** DEVIATION

**Description:**
The Frontera template shows:
```bash
ibrun ./seas_tpv102_driver \
    --config tpv102/config/tpv102.toml \
    --mesh tpv102/mesh/tpv102_fine.msh \
    --output-dir tpv102/results/${SLURM_JOBID}
```

This uses `--config`, `--mesh`, `--output-dir` flags. But the existing `seas_driver` uses positional TOML path + `--override` syntax:
```bash
ibrun ./seas_driver config/bp5_production.toml \
    --override output.output_dir="results_dir" \
    --override output.output_prefix="prefix"
```

The mesh is specified INSIDE the TOML file (`[mesh] file = "..."`) and the output directory is overridden via `--override output.output_dir=...`. Having two different CLI conventions for drivers in the same project creates confusion.

**Trigger:**
Users familiar with `seas_driver` trying to run TPV102 with the same syntax, or vice versa.

**Actual behavior:**
TPV102 driver uses a different CLI convention than seas_driver.

**Expected behavior:**
TPV102 driver should use the same TOML + `--override` pattern as `seas_driver`. The TOML file contains all parameters (mesh path, BCs, solver settings, etc.) and CLI overrides only adjust output paths.

**Suggested fix:**
In §3.4.2, change the driver invocation to match the existing convention:
```bash
ibrun ./seas_tpv102_driver tpv102/config/tpv102.toml \
    --override output.output_dir="tpv102/results/${SLURM_JOB_ID}"
```

And ensure `tpv102.toml` includes `[mesh] file = "tpv102/mesh/tpv102_fine.msh"` (consistent with `bp5_production_new_driver.toml` pattern).

---

### [R-008] [MODERATE] [Plan §3.5] — Phase 5 QD→FD velocity initialization creates discontinuity at fault tips

**Category:** EDGE_CASE

**Description:**
Section 5.2 describes warm-start velocity initialization:

> "At the fault, set v_{t1}^+ = +V_{qd,1}/2 and v_{t1}^- = −V_{qd,1}/2. In the bulk, set **v = 0**."

This creates a velocity discontinuity at the fault TIPS (where fault elements meet non-fault elements). Fault-adjacent elements have initialized velocity ±V_{qd}/2, while neighboring bulk elements have v = 0. This discontinuity will radiate artificial P-waves and S-waves at the fault tips during the first few FD time steps.

The plan's 5-step damped ramp (ramping nucleation perturbation from 0 to full strength over 5 steps) mitigates the nucleation transient, but it does NOT address the fault-tip discontinuity.

For BP5 with V_qd ~ 10⁻⁹ m/s (interseismic), this is negligible. But during a QD→FD transfer triggered by high slip rate (V_qd ~ 0.1 m/s at nucleation), the fault-tip discontinuity could produce significant artifacts.

**Trigger:**
QD→FD transfer during nucleation when slip rate is high near fault tips.

**Actual behavior:**
Sharp velocity jump from ±V_qd/2 to 0 at fault tips.

**Expected behavior:**
Smooth velocity taper near fault tips. Options:
1. Taper the initialized velocity over a few elements near fault tips: v(x) = V_qd(x)/2 × taper(dist_to_tip)
2. Extend the 5-step damped ramp to also apply to fault-tip velocities
3. Accept the artifact but increase the equilibrium correction step (Section 5.2, Problem 3) damping to absorb it

**Suggested fix:**
Add to §5.2 "Problem 4: Fault-tip velocity discontinuity" with the taper solution. Use the existing `FaultGeometry::GetFaultCoords2D()` to identify DOFs near fault tips and apply a Gaussian taper over ~3 elements.

**Test case:**
Test 41 (`TestTransientSuppression`): should verify that the transient energy at fault tips is < 1% of signal, not just at the fault center.

---

### [R-009] [LOW] [Plan §3.5.2] — FD→QD transfer interface is underspecified

**Category:** ASSUMPTION

**Description:**
Section 5.2 describes FD→QD transfer as:
1. Accumulated slip: δ_new = δ_frozen + ∫V dt
2. QD elasticity solve: **Ku = f(δ_new)**

But the existing QD solver (`SEASQuasiDynamicOperator`) doesn't have a `Solve(slip → displacement)` as a standalone entry point. The QD operator's `Mult()` computes the full RHS (domain solve + traction + friction). The transfer would need to:
1. Set the fault state vector (slip, psi) from FD accumulated values
2. Call `DomainOperator::Solve(time, slip_bc, displacement)` to get the new QD displacement field
3. Recompute traction to verify equilibrium

Step 2 exists (`ElasticityDomainOperator::Solve()`), but the plan doesn't specify:
- How `psi` (state variable) is transferred back — is it just the final FD value?
- Whether the initialization sequence (4-phase from CLAUDE.md) needs to be partially re-executed
- How to handle the fact that the QD displacement was computed with a DIFFERENT slip distribution than what the FD phase produced

**Trigger:**
Implementing `regime_transfer.cpp` and discovering that the QD re-initialization after FD isn't a simple "set slip and solve."

**Actual behavior:**
Plan says "one MUMPS solve, O(n_e^1.5)" but doesn't describe the full re-initialization sequence.

**Expected behavior:**
Add a subsection "5.1.2 FD→QD Transfer: Detailed Steps" with:
1. Set fault state: slip = δ_new, psi = psi_final_FD
2. Call `domain_->Solve(t, slip_bc, displacement)` — full quasi-static solve
3. Call `domain_->ComputeTraction(displacement, slip, traction)` — get new QD traction
4. Verify friction equilibrium: check that `|σ_n f(V, ψ) + η V - τ| < tol` at all fault DOFs
5. If not satisfied, run 1-2 Init-style correction steps (from the 4-phase initialization)

---

### [R-010] [LOW] [Plan §4.3] — Test count discrepancy: testing summary says 65, but table sums to 65 only if test_godunov_flux.cpp has 8 tests (not 14)

**Category:** QUALITY

**Description:**
The Testing Strategy Summary (§4.1) lists:

| Phase | File | # Tests |
|-------|------|---------|
| 1 | test_godunov_flux.cpp, test_wave_operator.cpp | 14 |

But the project layout (§4.3) lists:
- `test_godunov_flux.cpp` — Phase 1: **8 tests**
- `test_wave_operator.cpp` — Phase 1: **6 tests**

8 + 6 = 14. This is internally consistent (14 total for Phase 1). BUT: the Phase 1 description (§3.1) says "14 unit tests" with a reference to "see v3 Phase 1 for full details." If v3 had a different test breakdown, the reader can't verify the count.

More importantly, the test table doesn't include the test counts for:
- Phase 4a `test_tpv102_local.cpp` — listed as "6 local integration tests" in §3.4.1
- But §4.1 says "test_tpv102_setup.cpp, test_tpv102_local.cpp → 11 tests" (5 + 6 = 11 ✓)

Minor inconsistency but worth cleaning up to avoid confusion during implementation tracking.

**Suggested fix:**
In §4.1, split the Phase 1 row to show individual file counts:
```
Phase 1: test_godunov_flux.cpp (8), test_wave_operator.cpp (6)  → 14
```

---

## Summary
- Critical issues: 2 (R-001, R-002)
- Moderate issues: 6 (R-003, R-004, R-005, R-006, R-007, R-008)
- Low issues: 2 (R-009, R-010)
- Plan compliance: N/A (this is a plan review, not implementation review)
- Verdict: **PASS WITH FIXES** — R-001 (architecture) and R-002 (bracket signs) must be resolved before implementation begins. R-003 through R-008 should be addressed before the relevant phase.

### What the Plan Gets Right

The mathematical derivations (trial traction Eq. 7, imposed state Eqs. 11-12, friction balance Eq. 8) were verified against the eigenstructure of the plan's A-matrix and are **correct**. The equation-to-code dictionaries are a major improvement over v3. The reuse map is **verified accurate** by codebase inspection — FaultBasis, FaultGeometry, DieterichRuinaFriction, ConstitutiveModel, and ProbeOutput have stable interfaces that genuinely transfer to dynamic without modification. The phased testing strategy (65 tests across 6 phases) is thorough.

## Unreviewed Areas
- v3 plan content referenced by "see v3" in Sections 1.4, 2.1.3, 2.1.4, 3.1 — not independently verified
- SeisSol source code cross-check (plan references `ElasticSetup.h:28-77`) — not read during this review
- GPU kernel (Phase 6) correctness of MFEM_FORALL parallelism — requires MFEM GPU expertise
- Detailed LOC estimates per file — not independently validated
