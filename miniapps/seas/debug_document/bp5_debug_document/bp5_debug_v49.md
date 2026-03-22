# BP5 Debug v49: Multi-Direction Investigation Plan for 1000m Mesh Instability

**Date**: 2026-03-20
**Status**: PLANNING. Comprehensive investigation plan with prioritized directions.
**Previous**: v48 (1000m serial crash confirmed — eliminates parallel/BLR/shared face hypotheses)
**Branch**: `feature/elasticity`

---

## 1. Summary of What We Know

### 1.1 The Definitive v48 Finding

The 1000m mesh at p=2 with ×3 penalty crashes **even in serial** (v48i). The crash
signature is deterministic — identical tau_mag and slip to 6 digits across serial,
8-rank, and 400-rank runs with both exact MUMPS and BLR solvers.

| Hypothesis | Status | Evidence |
|-----------|--------|---------|
| Shared face sign bug | DISPROVED | Normals flip correctly (v48 Section 10) |
| Shared face data consistency | DISPROVED | Serial crashes too (v48i) |
| MUMPS-BLR accuracy | DISPROVED for 1000m | Crashes with exact MUMPS (v48f) |
| Parallel rank count | DISPROVED | Serial crashes (v48i) |
| **Formulation instability: 1000m + p≥2 + ×3 penalty** | **CONFIRMED** | All evidence converges |

### 1.2 What Triggers the Instability

The combination requires ALL THREE:
1. **1000m mesh resolution** (not 2500m) — higher penalty ∝ 1/h
2. **p ≥ 2** (multi-DOF, nbf=6) — p=1 works on 1000m
3. **×3 penalty correction** — 1/3 penalty avoids it (compensating error)

### 1.3 The Crash Mechanism

From v48i serial crash: `slip=(0.00131165, -0.236655)`. The strike component
(0.0013) matches V_nuc × dt = 0.01 × 0.131. The dip component (-0.237) is
182× too large and negative. This develops within the RK stages of the **first**
time step through a feedback cascade:

```
Stage 1: slip=0 → solve → traction perturbation δτ (the SEED)
Stage 2: δτ → velocity → slip → solve → amplified δτ'
Stage 3: δτ' → velocity → larger slip → ... → GPa blowup
```

The perturbation is **deterministic** (same values regardless of solver/ranks).
This means the seed δτ comes from the formulation itself, not numerical noise.

### 1.4 The Central Question

**Why does the 1000m mesh produce a non-zero seed perturbation δτ when slip=0,
while the 2500m mesh does not?**

Everything downstream (cascade amplification, penalty scaling) is secondary.
If the seed is zero, there is no cascade regardless of penalty magnitude.

### 1.5 Why Tandem Works

Tandem uses the same penalty formula (Uphoff et al.) at p=4 on similar meshes.
v47 Section 18 analysis showed the data flow is mathematically equivalent at
every step for flat faces. The difference must be in a subtle implementation
detail. Identifying this difference is the ultimate goal.

---

## 2. Root Cause Candidates

### Candidate 1: Quadrature Order Mismatch

K uses order 2p (6 quad points on triangle for p=2).
f uses order 2p+1 (7 quad points).

**Theoretical analysis**: Both rules are exact for degree 2p integrands on flat
faces. Should NOT cause mismatch. **But**: if any face in the 1000m mesh is
not perfectly flat (tet face vertices off-plane by floating point error), the
rules could give slightly different results. This difference would scale with
penalty magnitude.

**Probability**: Low. Linear tets have mathematically flat faces.
**Test difficulty**: Trivial — change one number.

### Candidate 2: Element Ordering Across Code Paths

All three code paths (K, f, traction) call `GetInteriorFaceTransformations`
with the same face index → same FTr → same Elem1/Elem2.

**Probability**: Very low.
**Test difficulty**: Part of the K/f diagnostic.

### Candidate 3: K/f Penalty Per-DOF Consistency

K penalty assembled through MFEM's integrator. f penalty through our code.
Same formula, same inputs. For flat faces with polynomial integrands, both
quadrature rules give exact results.

**Probability**: Very low for the penalty term alone. The formulas are provably
identical for flat faces.
**Test difficulty**: Moderate — requires building a diagnostic.

### Candidate 4: Traction Normal Inconsistency

ComputeTraction stress uses `basis.normal` (fixed ≈ (0,-1,0)), while the penalty
correction uses `sign` from `CalcOrtho`. If they disagree on a specific face,
the combined traction is wrong.

**Probability**: Low (both derive from CalcOrtho, just at different times).
**Test difficulty**: Trivial — print both normals.

### Candidate 5: Sharp Nucleation Boundary + Multi-DOF Interpolation (NEW)

The nucleation zone has a **sharp boundary**: V jumps from 0.01 to 1e-9 m/s
(7 orders of magnitude) between adjacent DOFs. At p≥2 with nbf=6 DOFs per
face, a face straddling this boundary has DOFs with wildly different slip values.

The multi-DOF slip is interpolated to quadrature points via polynomial basis
functions. For a sharp discontinuity, the polynomial interpolation creates
**Gibbs-like oscillations** — the interpolated slip at some quad points
OVERSHOOTS the actual DOF values. These overshoots create artificial traction
perturbations that seed the cascade.

At p=1 (nbf=1), each face has a single averaged slip → no interpolation
oscillation → stable. At 2500m, there are fewer faces at the boundary and
the penalty is lower → oscillations are below the cascade threshold.

**Probability**: HIGH. This directly explains all three trigger conditions:
- 1000m: more boundary faces, higher penalty amplifies oscillations
- p≥2: multi-DOF interpolation creates oscillations (p=1 doesn't interpolate)
- ×3 penalty: amplifies the oscillation-induced perturbation

**Test difficulty**: Very easy — smooth the nucleation boundary.

### Candidate 6: ODE Coupling Eigenvalue Exceeds RK Stability Limit

The coupled fault-elasticity system has an effective eigenvalue
`λ ≈ penalty / η`. For 1000m: `λ ≈ 4e10 / 4.62e6 ≈ 8665 s⁻¹`.
With dt=0.13s: `dt × λ ≈ 1127`, far beyond RK45's stability limit (~6).

**BUT**: if the solve is exact, the penalty cancels in the traction
computation (K×u=f → jump matches prescribed → penalty correction = 0).
The effective eigenvalue is then the stress coupling, NOT penalty/η.

The eigenvalue only involves penalty/η if there's a non-zero residual in
the jump. This residual can come from:
- Solver error (disproved — exact MUMPS also crashes)
- Interpolation oscillation (Candidate 5)
- K/f mismatch (Candidate 3)

**Probability**: Medium — this is the AMPLIFICATION mechanism, not the SEED.
The seed must come from another candidate.
**Test difficulty**: Hard to directly test eigenvalues.

---

## 3. Prioritized Investigation Plan

### Priority 0: Quick Experiments (< 10 lines of code each)

These are fast tests that can immediately confirm or eliminate hypotheses.
Each can be run in a single SLURM job and analyzed in minutes.

#### Direction 1: Smoothed Nucleation Boundary (Tests Candidate 5)

**Rationale**: This is the highest-probability candidate. If the sharp V
boundary is the trigger, smoothing it eliminates the Gibbs oscillation seed.

**Code change**: In the BP5 initialization function where V_init is set per DOF,
replace the sharp box boundary with a smooth transition:

```cpp
// BEFORE (sharp boundary):
real_t V = (inside_nucleation_zone) ? V_nuc : V_init;

// AFTER (smooth Gaussian taper):
real_t dist = distance_from_nucleation_center;
real_t r_nuc = nucleation_half_width;
real_t taper = 0.5 * (1.0 + std::tanh((r_nuc - dist) / (2.0 * h_element)));
real_t V = V_init + (V_nuc - V_init) * taper;
```

The taper width should be ~2× the element size so the transition spans at
least 2 elements (no single face straddles the full jump).

**Run on**: 1000m mesh, serial, p=2, ×3 penalty.
**Expected outcomes**:
- Crash disappears → Candidate 5 CONFIRMED as root cause.
  The sharp boundary + multi-DOF interpolation is the issue.
- Still crashes → Candidate 5 eliminated. The seed is elsewhere.

**Effort**: ~10 lines. **Priority**: HIGHEST.

#### Direction 2: Match Quadrature Order (Tests Candidate 1)

**Rationale**: The cheapest possible code change. Even though theoretical
analysis says it shouldn't matter, testing costs nothing.

**Code change**: In `AssembleSlipContributionIP` and
`AssembleSlipContributionIPShared`, change quadrature order from 2p+1 to 2p:

```cpp
// BEFORE:
const IntegrationRule &ir = IntRules.Get(FTr->FaceGeom, 2 * face_order + 1);

// AFTER:
const IntegrationRule &ir = IntRules.Get(FTr->FaceGeom, 2 * face_order);
```

Also match in `ComputeTraction` (which also uses 2p+1).

**Run on**: 1000m mesh, serial, p=2, ×3 penalty.
**Expected outcomes**:
- Crash disappears → quadrature mismatch IS the issue (surprising but
  possible if faces aren't perfectly flat in floating point).
- Still crashes → Candidate 1 eliminated.

**Effort**: Change 3 numbers. **Priority**: HIGH (trivial to test).

#### Direction 3: Normal Consistency Check (Tests Candidate 4)

**Rationale**: If `CalcOrtho` normal and `basis.normal` disagree on any face,
the combined traction (stress + penalty) is inconsistent.

**Code change**: In `ComputeTraction`, add a one-time print comparing both
normals for all fault faces:

```cpp
static bool nor_check_done = false;
if (!nor_check_done)
{
   for (int fi = 0; fi < fault_interior_faces_.Size(); fi++)
   {
      Vector nor(dim);
      CalcOrtho(FTr->Jacobian(), nor);
      const auto &basis = fault_basis_.GetBasis(fi);
      real_t dot = 0;
      for (int d = 0; d < dim; d++)
         dot += nor(d) / nor.Norml2() * basis.normal[d];
      if (std::abs(std::abs(dot) - 1.0) > 1e-10)
         mfem::out << "[NOR-CHECK] fi=" << fi << " MISMATCH dot=" << dot
                    << " nor_unit=(...) basis.normal=(...)\n";
   }
   nor_check_done = true;
}
```

**Run on**: Both 2500m and 1000m.
**Expected outcomes**:
- Some faces have |dot| ≠ 1 → normals are misaligned → likely bug.
- All faces have |dot| = 1 → Candidate 4 eliminated.

**Effort**: ~15 lines. **Priority**: HIGH (quick, eliminates a candidate).

---

### Priority 1: Trace the Instability Seed

If Priority 0 experiments don't resolve the issue, we need to find exactly
WHERE the initial perturbation comes from.

#### Direction 4: First-Stage Traction Dump

**Rationale**: At the first RK evaluation (stage 1), slip=0 for all DOFs. The
elasticity solve with zero slip should produce traction that exactly matches
the equilibrium stress. Any deviation is the SEED of the cascade.

**Implementation**: After the first ComputeTraction call (detected by checking
if all slip values are zero), dump per-DOF traction for all fault DOFs in the
nucleation zone:

```cpp
if (first_call && max_slip < 1e-20)
{
   for (int i = 0; i < num_fault_dofs_; i++)
   {
      real_t tau_dip = traction(2*i);
      real_t tau_strike = traction(2*i+1);
      if (std::abs(tau_dip) > 1.0 || std::abs(tau_strike) > 1.0)
      {
         // Traction should be ~0 when slip=0 (equilibrium initialized in τ_pre)
         mfem::out << "[SEED] DOF=" << i << " tau=(" << tau_dip << ","
                    << tau_strike << ") at x=(...)\n";
      }
   }
}
```

**What this reveals**: Which DOFs have non-zero initial traction perturbation.
If these DOFs are at the nucleation boundary faces → supports Candidate 5.
If they're uniformly distributed → suggests solver or formulation issue.

**Run on**: 1000m serial AND 2500m serial (compare seed magnitude).
**Effort**: ~30 lines. **Priority**: HIGH.

#### Direction 5: RK Stage-by-Stage State Dump

**Rationale**: Traces exactly how the cascade grows through RK stages.
Identifies which stage, which DOF, and which component (strike vs dip)
first exceeds physical bounds.

**Implementation**: In the DormandPrinceRK45 stepper, after each stage
evaluation, print summary stats:

```cpp
// After computing k_[stage]:
real_t max_V_dip = 0, max_V_strike = 0, max_slip_dip = 0;
int max_V_dip_dof = -1;
for (int i = 0; i < num_dofs; i++)
{
   real_t V_dip = std::abs(k_[stage](i * 3 + 0));
   real_t V_str = std::abs(k_[stage](i * 3 + 1));
   if (V_dip > max_V_dip) { max_V_dip = V_dip; max_V_dip_dof = i; }
   if (V_str > max_V_strike) { max_V_strike = V_str; }
}
mfem::out << "[RK-STAGE " << stage << "] max_V_dip=" << max_V_dip
          << " (DOF " << max_V_dip_dof << ")"
          << " max_V_strike=" << max_V_strike
          << " max_slip_dip=" << max_slip_dip << "\n";
```

**What this reveals**: The exact amplification factor per stage.
If V_dip grows 1000× per stage → eigenvalue analysis is correct.
If V_dip grows slowly for 5 stages then explodes → non-linear trigger.

**Run on**: 1000m serial (will crash, but prints stages 1-3+ before crash).
**Effort**: ~40 lines. **Priority**: MEDIUM-HIGH.

---

### Priority 2: Rule Out / Confirm Candidates

These are more involved diagnostics that methodically eliminate hypotheses.

#### Direction 6: K/f Consistency Diagnostic (Tests Candidate 3)

**Rationale**: Even though theoretically K×u should equal f for flat faces,
this diagnostic provides DEFINITIVE proof. If K×u = f to machine precision,
Candidate 3 is eliminated beyond doubt. If not, we found the bug.

**Implementation**: The detailed diagnostic from the original v49 Sections 3-6:
- Call `DGElasticityIPPenaltyIntegrator::AssembleFaceMatrix` for a fault face
- Construct u_test with prescribed slip jump
- Compare K×u_test with f from our assembly
- Report per-DOF differences

**Faces to test**: 3 faces — nucleation interior, nucleation boundary, far field.

**Run on**: Both 2500m (validation) and 1000m (actual test).
**Expected**: K×u = f on both meshes to ~1e-14 relative error.
**Effort**: ~100 lines. **Priority**: MEDIUM.

#### Direction 7: Penalty Magnitude and Mesh Quality Analysis

**Rationale**: Understanding the QUANTITATIVE differences between meshes.

**Implementation**: Dump for all fault faces:
```
[MESH-DIAG] fi=<idx> penalty=<val> nl_q=<val> detJ1=<val> detJ2=<val>
            A_face=<val> V_elem1=<val> V_elem2=<val> aspect_ratio=<val>
            centroid=(<x>,<y>,<z>) zone=<nuc_interior|nuc_boundary|far_field>
```

Compute summary statistics:
- Mean/max/min penalty on 2500m vs 1000m
- Aspect ratio distribution
- Number of nucleation-boundary faces on each mesh
- Penalty × nbf product (total DOF-weighted coupling strength)

**What this reveals**: Whether the 1000m mesh has specific faces with
anomalous penalty (e.g., very thin elements) that the 2500m mesh doesn't.

**Run on**: Both meshes, serial.
**Effort**: ~50 lines. **Priority**: MEDIUM.

#### Direction 8: Residual Check After First Solve

**Rationale**: After the first solve (slip=0), check ||K×u - f|| / ||f||.
If the residual is large, the solver is inaccurate on this system.

**Implementation**: Already partially present (`--check-residual` flag).
Enhance to also compute the displacement jump at fault face DOFs:
```
max |Σ_k shape1(k,q)*u1(k) - Σ_k shape2(k,q)*u2(k)| for each fault face
```
With slip=0, this should be near-zero. Any non-zero jump is solver error.

**Run on**: 1000m serial.
**Effort**: ~40 lines. **Priority**: MEDIUM.

---

### Priority 3: Deep Analysis

If Priorities 0-2 don't identify the root cause, these provide deeper insight.

#### Direction 9: Tandem on Same Mesh

**Rationale**: Verify that Tandem actually works on the bp5_tandem.msh mesh
at 1000m scale with p=2. If Tandem also fails, the issue is fundamental to
this mesh, not our code.

**Implementation**: Configure Tandem with the same mesh and p=2.
**Effort**: Configuration + SLURM job. **Priority**: LOW (time-consuming).

#### Direction 10: Eigenvalue Analysis of Coupled System

**Rationale**: Compute the largest eigenvalue of the linearized ODE system
dF/dy around the initial state. If |λ_max × dt| > 6 (RK stability limit),
the system is formally unstable with explicit RK at this dt.

**Implementation**: Form the Jacobian numerically (perturb each state variable,
measure the response) for a reduced system (one face). Very expensive for the
full system.

**Alternative**: Use the MUMPS determinant (INFOG(12,13)) to estimate the
condition number of K, then bound the coupled eigenvalue.

**Effort**: High. **Priority**: LOW (theoretical, not actionable).

#### Direction 11: 1500m Mesh Test (Intermediate Resolution)

**Rationale**: Find the critical mesh resolution where the instability first
appears. If 2500m works and 1000m crashes, where's the threshold? Testing at
1500m or 2000m narrows the search and reveals whether the transition is
gradual (conditioning) or sharp (specific mesh feature).

**Implementation**: Generate a 1500m mesh and run serial p=2 with ×3 penalty.
**Effort**: Mesh generation + SLURM job. **Priority**: LOW.

---

## 4. Execution Sequence

### Phase 1: Quick Experiments (can run in parallel)

All three Priority 0 experiments can be implemented in a single code version
and run as separate SLURM jobs:

| Job | Direction | Code change | Mesh | What it tests |
|-----|-----------|-------------|------|--------------|
| v49a | 1 (smooth nuc) | Gaussian taper in V_init | 1000m serial | Candidate 5 |
| v49b | 2 (quad match) | 2p+1 → 2p in RHS | 1000m serial | Candidate 1 |
| v49c | 3 (normal check) | Print normals | 1000m serial | Candidate 4 |
| v49d | 4 (traction seed) | Print first-stage τ | 1000m serial | Seed location |
| v49e | 4 (traction seed) | Print first-stage τ | 2500m serial | Seed comparison |

**Decision point after Phase 1:**

```
v49a (smooth nuc) works?
├── YES → Candidate 5 confirmed.
│   The sharp boundary is the trigger. Fix: smooth the nucleation
│   initialization (match Tandem's approach). Done.
└── NO → Candidate 5 eliminated. Continue Phase 2.

v49b (quad match) works?
├── YES → Candidate 1 confirmed (surprising).
│   Fix: use 2p quadrature everywhere. Done.
└── NO → Candidate 1 eliminated.

v49c (normals) show mismatch?
├── YES → Candidate 4 confirmed.
│   Fix: align traction normals. Done.
└── NO → Candidate 4 eliminated.

v49d/e (traction seed) — WHERE is the seed?
├── Nucleation boundary DOFs → supports Candidate 5 mechanism
├── Uniformly distributed → suggests eigenvalue/conditioning issue
└── No seed on 2500m, seed on 1000m → confirms mesh-dependent source
```

### Phase 2: Detailed Diagnostics (only if Phase 1 inconclusive)

| Job | Direction | What it tests |
|-----|-----------|--------------|
| v49f | 5 (RK stages) | Cascade growth rate per stage |
| v49g | 6 (K/f dump) | Candidate 3 |
| v49h | 7 (mesh quality) | Element quality differences |
| v49i | 8 (residual) | Solver accuracy |

### Phase 3: Deep Analysis (only if Phase 2 inconclusive)

| Job | Direction | What it tests |
|-----|-----------|--------------|
| v49j | 9 (Tandem) | Reference comparison |
| v49k | 10 (eigenvalue) | Stability theory |
| v49l | 11 (1500m mesh) | Resolution threshold |

---

## 5. Implementation Notes

### 5.1 Smoothed Nucleation (Direction 1)

The BP5 nucleation zone is defined by the rectangular region:
- Along strike: |x| < l_vw/2 = 30 km
- Along dip: depth between H - Wf + (Wf - w_nuc)/2 and H - (Wf - w_nuc)/2

The initial V is V_nuc inside, V_init outside. The smoothing should use a
function of the distance from the nucleation zone boundary.

Check Tandem's initialization: does Tandem use a sharp or smooth boundary?
From `tandem/examples/tandem/3d/bp5.lua`:
```lua
bp5_outside = BP5.new({eps=1e-3})
```
The `eps=1e-3` parameter slightly enlarges the nucleation zone by 1mm.
This is effectively sharp. **But Tandem uses WarpAndBlend nodes for p≥3
which may provide natural smoothing through the interpolation basis.**

At p=2 with GaussLobatto nodes, WarpAndBlend and GaussLobatto are identical
(v47 Section 18.3). So the node placement is not the difference.

**Key question**: Does Tandem evaluate the nucleation zone membership per-node
or per-face? If per-face (using face centroid), all DOFs on a face get the
same V → no intra-face variation → no interpolation oscillation. If per-node
(using DOF coordinates), boundary faces have mixed V values.

Look at Tandem's `bp5.lua` Vinit function — it takes a node coordinate and
returns V. So it IS evaluated per-node. But then Tandem should have the same
sharp boundary issue...

Unless Tandem's solver (iterative, not direct) naturally damps the oscillation
through the residual tolerance.

### 5.2 First-Stage Traction Dump (Direction 4)

The traction from ComputeTraction with slip=0 should be purely from the
stress term `{σ(u)·n̂}` (the penalty correction is zero when [[u]] = 0
and δu = 0). The stress comes from the equilibrium displacement under
far-field BCs.

The total traction seen by the friction law is τ_total = τ_pre + τ_computed.
At initialization, τ_pre is set so that τ_total = σ_n × f(V_init, ψ_init) + η × V_init.

If the traction computation returns exactly the expected value, τ_total is
at equilibrium and V doesn't change. Any deviation creates a velocity
perturbation.

The threshold for cascade onset: if `δτ × dt / η > ε_critical`, the
perturbation survives one RK stage and grows. With penalty/η ≈ 8665 for
1000m, even δτ = 1 Pa creates δV ≈ 1/4.62e6 ≈ 2e-7 m/s, which over dt=0.13s
gives δslip ≈ 3e-8 m, which at stage 2 gives δτ' ≈ penalty × 3e-8 ≈ 1200 Pa.
Amplification: 1200× per stage. After 3 stages: 1.7e9 → GPa. **So even a
1 Pa initial perturbation blows up on the 1000m mesh.**

On the 2500m mesh: penalty/η ≈ 3466, amplification ≈ 450× per stage. After
3 stages: 9e7 → 90 MPa. This is comparable to σ_n × f ≈ 15 MPa, so the
friction law non-linearity kicks in and may stabilize. **The 2500m mesh is
RIGHT AT the stability boundary.**

This analysis predicts:
- Any mesh with penalty/η > ~1000 will be unstable at dt=0.13s
- The critical element size is h_crit where penalty(h_crit)/η ≈ threshold
- Smoothing the nucleation boundary reduces the initial perturbation from
  ~1-10 Pa to ~0.001 Pa, which may keep the cascade below GPa

### 5.3 RK Stage Dump (Direction 5)

The RK45 Dormand-Prince has 7 stages. The state vector is
(s_dip, s_strike, ψ) per DOF. The rate vector k_[stage] contains
(ds_dip/dt, ds_strike/dt, dψ/dt) = (V_dip, V_strike, dψ/dt).

At stage 1: V should be the equilibrium value (V_init or V_nuc).
At stage 2: V should be nearly the same (small perturbation).
If V_dip at any DOF exceeds 1 m/s by stage 3, the cascade is growing.

Print: stage number, max |V_dip|, DOF index of max, max |slip_dip|.

---

## 6. Decision Tree (Complete)

```
Phase 1 (Quick Experiments)
│
├── v49a: Smoothed nucleation on 1000m serial
│   ├── STABLE → ROOT CAUSE IS THE SHARP NUCLEATION BOUNDARY
│   │   Action: Investigate how Tandem handles this.
│   │   - Does Tandem use a smooth boundary? (check bp5.lua more carefully)
│   │   - Does Tandem's iterative solver naturally damp oscillations?
│   │   - Implement proper smooth initialization. Run benchmark.
│   │   **INVESTIGATION COMPLETE.**
│   └── CRASH → Not the sharp boundary. Continue.
│
├── v49b: Quadrature order match on 1000m serial
│   ├── STABLE → ROOT CAUSE IS QUADRATURE ORDER MISMATCH
│   │   Action: Use 2p everywhere. Investigate why flat faces show this.
│   │   **INVESTIGATION COMPLETE.**
│   └── CRASH → Not quadrature order. Continue.
│
├── v49c: Normal consistency check
│   ├── MISMATCH found → ROOT CAUSE IS NORMAL INCONSISTENCY
│   │   Action: Align traction normals. Test fix.
│   │   **INVESTIGATION COMPLETE.**
│   └── All normals consistent → Not normals. Continue.
│
├── v49d/e: First-stage traction dump (1000m + 2500m)
│   ├── Seed found at nucleation boundary DOFs
│   │   → Confirms oscillation hypothesis even if smooth test didn't fix it
│   │   → Investigate slip interpolation in detail
│   ├── Seed found uniformly across DOFs
│   │   → Solver accuracy or conditioning issue
│   │   → Run Direction 8 (residual check)
│   └── No seed > 1 Pa on either mesh
│       → The perturbation comes from WITHIN the RK stages, not stage 1
│       → Run Direction 5 (RK stage dump)
│
└── Phase 1 decision:
    ├── Root cause found → Fix and verify
    └── No root cause yet → Phase 2

Phase 2 (Detailed Diagnostics)
│
├── v49f: RK stage-by-stage dump
│   → Identifies exact amplification factor and trigger stage
│
├── v49g: K/f consistency dump
│   → Definitively confirms or eliminates Candidate 3
│
├── v49h: Mesh quality analysis
│   → Identifies anomalous elements on 1000m mesh
│
└── v49i: Solver residual check
    → Measures actual solve accuracy on 1000m

Phase 3 (Deep Analysis, if needed)
│
├── v49j: Tandem on same 1000m mesh → reference comparison
├── v49k: Eigenvalue estimation → stability theory
└── v49l: 1500m mesh → find resolution threshold
```

---

## 7. Revision History

| Version | Change | Status |
|---------|--------|--------|
| v48 | Sign hypothesis disproved; BLR partially valid for 2500m | Done |
| v48i | 1000m serial crash confirmed | DEFINITIVE |
| **v49** | **Multi-direction investigation plan** | **PLANNING** |
| v49a (planned) | Smoothed nucleation test | P0 — HIGHEST |
| v49b (planned) | Quadrature order match test | P0 |
| v49c (planned) | Normal consistency check | P0 |
| v49d (planned) | First-stage traction dump (1000m) | P1 |
| v49e (planned) | First-stage traction dump (2500m) | P1 |
| v49f (planned) | RK stage-by-stage dump | P1 |
| v49g (planned) | K/f consistency diagnostic | P2 |
| v49h (planned) | Mesh quality analysis | P2 |
| v49i (planned) | Solver residual check | P2 |
