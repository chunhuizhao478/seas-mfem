# SEAS Miniapp - Claude Code Instructions

MFEM-based SIPG/BR2 Discontinuous Galerkin code for SCEC SEAS benchmark problems (BP1, BP2, BP5). Couples a DG elasticity domain solver with regularized Dieterich-Ruina rate-and-state friction on an embedded fault interface via a quasi-dynamic approximation. Time integration uses adaptive Dormand-Prince RK45.

## Before ANY Code Change

1. Read `ARCHITECTURE.md` for class hierarchy, data flow, and execution model
2. Read debug documents in `debug_document/bp5_debug_document/` for history of what has been tried
3. Run `make test` (unit tests) to confirm baseline passes
4. For BP5 changes: check Tandem reference at `/Users/chunhuizhao/projects/tandem/`

## Reference Materials

| Resource | Location |
|----------|----------|
| Tandem source code | `/Users/chunhuizhao/projects/tandem/` |
| Porting reports | `/Users/chunhuizhao/projects/mfem/miniapps/seas/document/` |
| Debug history | `debug_document/bp5_debug_document/bp5_debug_v*.md` |
| SCEC BP5 spec | https://strike.scec.org/cvws/seas/download/SEAS_BP5_QD.pdf |
| Architecture | `ARCHITECTURE.md` (this directory) |

## Critical Numerical Details

These lessons were learned through extensive debugging (v1-v62). Violating any of them causes simulation failure.

### Sign Conventions

- **Slip rate direction**: PARALLEL to traction tau (not antiparallel). `V_vec = (V_abs / tau_abs) * tau_vec`. Antiparallel sign creates positive feedback -> unbounded growth (debug v8).
- **Pre-stress direction**: PARALLEL to initial velocity. `tau_pre = tau0_scalar * Vi / |Vi|` (debug v8).
- **Dip direction**: (0, 0, +1) = downward into earth. Was wrongly (0, 0, -1) = upward (debug v10).
- **DG face sign**: `sign = (nor(0) > 0) ? -1.0 : 1.0` applied consistently to slip embedding and traction extraction.
- **Normal stress**: sigma_n > 0 = compression (geology convention).
- **Depth coordinate**: Z=0 at surface, Z<0 is depth. `Wf` is positive, fault extends from z=0 to z=-Wf.
- **Fault-local tangent frame (BP5 / TPV102, project-wide)**: Uses the
  `FaultBasis` (Tandem) convention — `tangent1 = dip, tangent2 = strike`.
  For TPV102's vertical y=0 fault with `ref_normal=(0,-1,0)` and `up=(0,0,1)`:
  `can_t1 = (0, 0, -1)` (down-dip), `can_t2 = (+1, 0, 0)` (along strike).
  Therefore in `DOFData`: `V1/slip1/tau1_0/tau1_corr` are the **dip** components
  and `V2/slip2/tau2_0/tau2_corr` are the **strike** components.  TPV102 is
  pure strike-slip, so `tau2_0 = tau_ini`, `V2 = V_ini`, and `tau1_0 = 0`,
  `V1 = 0` (debug v7.0.0 R-801).  Pre-R-801 the interior-fault branch used
  `GodunovFlux::BuildFrame` (t1=strike) while the shared-fault branch used
  BP5 (t1=dip); the two have been unified on the BP5 convention.

### Friction Solver

- **Must use Brent's method**, not Newton. With large psi/a, Newton fails because F(V_lo=1e-30) > 0 while true solution is V << 1e-30 (debug v1).
- **Degenerate bracket guard**: When psi << 0, check `if (Fb >= 0.0) return V_hi = tau/eta` (frictionless limit). Allows RK45 to reject gracefully (debug v7).
- **No artificial tau floor**: Use `if (tau_abs <= 0.0) return 0`, not `if (tau_abs < 1e-30)`. The latter creates discontinuity in ODE RHS (debug v13 H7).
- **Vector slip (BP5)**: Solve for |V| via Brent in log10(V) space, then decompose direction separately.

### DG Method

- **IP penalty must include elasticity tensor coupling** for 3D. Scalar-only penalty causes blowup during nucleation (debug v12-v13 H3).
- **BR2 face_int uses unnormalized normal** from `CalcOrtho`, while `TestNormal` uses unit normal. These are different by |J_F| (debug v11-v12).
- **Shared faces**: Must use `GetSharedFaceTransformations()`. Only Elem1 is local; Elem2 via `pfes->GetFaceNbrFE()`. Slip assembly: only add Elem1 RHS contributions.
- **Stiffness matrix**: Assembled once and reused across all time steps and RK stages.

### Boundary Conditions (BP5)

- **All non-free-surface boundaries must get Dirichlet loading**: `u_X = sgn(Y) * Vp * t / 2`. With only bottom-boundary Dirichlet, nucleation zone starves (6.3 vs 50 kPa/yr loading rate), causing 3.7x longer recurrence (debug v9, v14-v15).
- **BCMode::FarField** is correct: Dirichlet on attrs 1-4 (vertical faces), Natural on 5-6 (top/bottom).
- Tandem mesh tags: Physical Surface 1 = Natural, 3 = Fault, 5 = Dirichlet.

### Time Stepping

- **Initial dt**: `0.01 * L_nuc / V_nuc` (not 0.5*). Too large -> RK45 stage amplification during nucleation (debug v7).
- **Max dt**: `0.1 * seconds_per_year` for SCEC output compliance (debug v13 H4).
- **RK45 tolerances**: atol=1e-7, rtol=1e-50 (pure absolute, matches PETSc TS default).
- **Parallel error reduction**: `MPI_Allreduce(MPI_MAX)` prevents deadlock from rank divergence.

### Initialization

- **4-phase**: (1) PreInit: slip=0, psi=steady-state. (2) Solve domain -> traction. (3) Init: compute psi from equilibrium. (4) Verify equilibrium error < 1e-6.
- **Per-DOF Dc in BP5**: `PsiToTheta(psi)` must use per-DOF Dc (0.13m in nucleation, 0.14m elsewhere), not global constant (debug v13 H1).

### Tandem vs SCEC Spec

- Tandem uses V_nuc=0.01, no delta_tau. SCEC spec: V_nuc=0.03, delta_tau=eta*V_nuc.
- Made configurable: `--V-nuc` and `--delta-tau-factor` CLI flags.
- **Follow SCEC spec unless explicitly matching Tandem** for cross-verification.

## Files Requiring Extreme Care

Any change to these files requires running full verification tests:

| File | Risk | Why |
|------|------|-----|
| `friction/dieterich_ruina.hpp` | Friction solver, Brent method | Wrong bracket/tolerance -> silent wrong V |
| `domain/elasticity_operator.hpp` | DG assembly, traction, BCs | Changes affect every simulation output |
| `fault/fault_basis.hpp` | Coordinate transforms | Sign error -> positive feedback -> blowup |
| `fault/rate_state_fault.hpp` | State evolution, slip assembly | Wrong interleaving -> corrupted state |
| `solver/seas_operator.hpp` | Coupling logic, init sequence | Wrong phase order -> bad equilibrium |
| `solver/time_stepper.hpp` | RK45, error control | Parallel error reduction critical for MPI |
| `config/bp5_params.hpp` | Spatial parameter functions | Wrong a(z) or Dc -> wrong friction regime |

## Building and Testing

### Environment

```bash
conda activate mfem-dev       # For building (mpicxx, MPI, MUMPS, etc.)
conda activate pythonenv       # For Gmsh mesh generation
```

### Build

```bash
cd miniapps/seas
make all                       # Build everything (tests + drivers)
make seas_pseas                # Build parallel BP2 driver only
make seas_bp5_full             # Build BP5 full verification driver
```

### Testing

```bash
# Unit tests (fast, ~2 min total)
make test

# Individual test targets
make test-friction             # Friction law validation
make test-elasticity-operator  # DG elasticity assembly
make test-bp5-integration      # BP5 integration test
make test-bp5-smoke            # BP5 parallel smoke test (4 ranks)

# Parallel tests
mpirun -np 4 seas_test_parallel_elasticity
mpirun -np 8 seas_test_bp5_parallel_smoke

# Full verification (long-running)
mpirun -np 8 seas_bp5_full --mesh bp5/mesh/bp5_1000m.msh --tfinal 56844000000
```

### What Constitutes a Regression

1. V_max monotonically increasing (should peak then decrease in events)
2. Traction > 1 GPa (physical limit ~25 MPa for sigma_n = 25 MPa)
3. dt going to zero or NaN
4. zeroIn bracket failure (F(a) and F(b) same sign)
5. Recurrence >> 300 years or << 100 years (expected ~240 years)
6. Event slip << 1m or >> 10m
7. Significant dip slip in BP5 (should be ~0 for pure strike-slip)

## Do NOT

- Change sign conventions without understanding the full chain (slip -> traction -> friction -> state)
- Hardcode numerical constants; derive from parameters/mesh for generalizability
- Use Newton-Raphson for the friction solver (Brent is required)
- Apply scalar-only DG penalty for 3D elasticity (must include elasticity tensor)
- Mix structural refactoring with numerical changes in one commit
- Revert a previous fix without citing the debug document and getting explicit approval
- Assume Tandem is always correct; follow SCEC spec when they differ
