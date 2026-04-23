# TPV102 Nucleation Fix — Persistent Fault-Prestress Channel (2026-04-22)

## Companion to
- `tpv102_nucleation_code_review_2026-04-22.md` (the bug analysis).
- The user's verification on 2026-04-22 confirming Finding 1 (CRITICAL),
  Finding 2 (MODERATE), Finding 3 (LOW), and the SeisSol
  `initialStressInFaultCS`-analog as the right abstraction.

## Root cause (one-liner)
`ApplyNucleationTotal` injects nucleation as a single-DOF point source into
bulk `Q[SXY]`. The wave operator (`wave_operator.inl:905-913`) radiates that
poke outward in `O(h/cp) ≈ 3×10⁻⁵ s` — far shorter than the 1 s nucleation
ramp, so the friction solver sees a residual ~kPa instead of the requested
25 MPa.

## Fix direction (implemented here)
Option (a) from the user's verdict: **add a separate per-DOF nucleation
prestress field that `EvaluateTotal` reads at every call** — the SeisSol
`initialStressInFaultCS` analog. The static background prestress contract
(bulk Q carries it; `tau*_0 = sigma_n0 = 0` enforced by assert) is
preserved; the new `tau*_nuc` channel is orthogonal to that assert and
carries the time-varying nucleation perturbation.

This mirrors the working fluctuation-Q `ApplyNucleation` at
`tpv102_setup.hpp:135-148`, which OVERWRITES `tau2_0 = tau_ini + dtau` per
call. We just write to a different DOFData slot so the total-Q assert keeps
catching the original double-counting bug.

## Changes

### 1. `dynamic/fault_face_flux.hpp` — DOFData
Add three nucleation prestress fields (zero-initialized):
```cpp
real_t sigma_n_nuc = 0, tau1_nuc = 0, tau2_nuc = 0;
```
TPV102 only writes `tau2_nuc` (pure strike-slip). The other two are kept
zero and are present only for symmetry with the existing
`{sigma_n0, tau1_0, tau2_0}` triple (no extra cost — same struct layout
expansion).

### 2. `dynamic/fault_face_flux.cpp` — EvaluateTotal
After `ComputeTrialTraction(...)`, build separate friction-input totals
WITHOUT mutating `tau*_trial`:
```cpp
const real_t sigma_n_fric = sigma_n_trial + data.sigma_n_nuc;
const real_t tau1_fric    = tau1_trial    + data.tau1_nuc;
const real_t tau2_fric    = tau2_trial    + data.tau2_nuc;
```
- Friction solve and slip-rate decomposition use `*_fric` (TOTAL,
  parallels SeisSol `totalTraction*` at `RateAndState.h:222-240`).
- Riemann imposed-state path uses `tau*_corr = tau*_trial - eta_s*V`
  (TRIAL scale — parallels SeisSol `tractionResults.traction*`).
- Output for stations: `data.tau*_corr = tau*_corr + tau*_nuc` (TOTAL).

**Round-2 correction (R-001 CRITICAL fix from REVIEW.md)**: the first
draft mutated `tau*_trial += data.tau*_nuc` directly, which leaked the
nucleation amplitude into the Riemann velocity-jump
`(2/Zs)·(tau_corr - Q_bulk[SXY])`.  At full nucleation
(tau2_nuc=25 MPa, Zs≈9.25 MPa·s/m, V_abs≈0.22 m/s) this overshot by
~24× and would have radiated a spurious shear pulse into bulk Q every
step — the very wave-radiation pathology the fix was supposed to
avoid.  The `*_fric` / `*_corr` separation above is bit-for-bit the
SeisSol semantics.

The `data.psi`, `{sigma_n0, tau1_0, tau2_0} == 0` asserts are unaffected.

### 3. `dynamic/tpv102_setup_total.hpp` — new ApplyNucleationTotalPrestress
A new overwrite-style function (parity with the fluctuation-Q
`ApplyNucleation` at `tpv102_setup.hpp:135-148`):
```cpp
inline void ApplyNucleationTotalPrestress(std::vector<DOFData> &dof_data,
                                          const std::vector<Vector> &fault_coords,
                                          real_t t);
```
Per fault DOF: `dof_data[i].tau2_nuc = NucleationPerturbation(x, z, t)`.
The legacy bulk-Q `ApplyNucleationTotal(Vector &Q, ...)` and its
`FaultQPNodalMap` / `FaultQPNucleationState` machinery REMAIN in place so
the existing R-002 / R-009 / R-I06-003 unit tests continue to pass — those
tests verify per-call arithmetic, not the multi-step radiation problem.
The driver no longer calls them.

### 4. `drivers/tpv102_driver.cpp`
- Remove the `FaultQPNodalMap fault_qp_map;`,
  `FaultQPNucleationState nuc_state;`, and `BuildFaultQPNodalMap(...)`
  setup block.
- Replace each of the 5 `ApplyNucleationTotal(Q, fault_qp_map, nuc_state,
  fault_coords, ...)` call sites with
  `ApplyNucleationTotalPrestress(dof_data, fault_coords, t_arg)`.
- The driver still calls `ZeroDOFDataPreStressTotal(dof_data, ...)` so the
  static-prestress assert in `EvaluateTotal` continues to fire on
  contract violation; `tau2_nuc` is unaffected by that assert.

### 5. `tests/unit/test_persistent_nuc_prestress_channel.cpp` (new)
Unit-test the fix without any production-mesh run.  20 subtests, all
green, runtime <1 s on the local laptop, no MPI:
- **Test A — additivity + Riemann decoupling**: build a minimal
  DOFData with symmetric fault-local Q_± (`SXX = +sigma_n`,
  `SXZ = +tau_ini`; sign per `InitializeStateTotal` rotation), call
  `EvaluateTotal` with and without `tau2_nuc = nuc_dtau`. Subtests A1-A7
  check the additive shift and orthogonality of `tau1_nuc` /
  `sigma_n_nuc` channels.  **Subtest A8 (R-002 fix)**: assert
  `|Q_imp_plus[VZ] - Q_imp_minus[VZ]| == V_abs` to within 5%.  Under
  the round-1 `tau*_trial`-mutation bug this would have been ~24×
  V_abs; under the corrected `*_fric` / `*_corr` separation it is
  exactly V_abs (the friction-reaction contribution only).
- **Test B — overwrite semantics**: 3 fault QPs (hypo, off-hypo,
  far-field), pre-poison `tau2_nuc`, call `ApplyNucleationTotalPrestress`,
  verify each QP overwritten to `NucleationPerturbation(x, z, t)`.
  Idempotency + temporal-ramp tracking checks.
- **Test C — default no-op**: `tau*_nuc=0` produces bit-exact identical
  `Q_imp` and `data.*_corr` as a fixture without the new channel.
  Locks in that the round-1 fix does not break callers that don't opt in.

### 6. `dynamic/tpv102_setup_total.hpp` — ZeroDOFDataPreStressTotal (R-004)
Extended to also reset `{sigma_n_nuc, tau1_nuc, tau2_nuc} = 0` per DOF.
Today this is a no-op (the new fields default-construct to 0), but a
future re-init path that calls `InitializeFaultDOFs` after a non-zero
nucleation history (checkpoint restart, test fixture rebuild) would
otherwise carry stale nuc state across the reset.

## What this fix does NOT change
- The bulk-Q background prestress (in `Q[SXY]` via
  `InitializeStateTotal`).
- The `EvaluateTotal` zero-prestress assert on
  `{sigma_n0, tau1_0, tau2_0}`.
- The fluctuation-Q `ApplyNucleation` at `tpv102_setup.hpp:135-148`
  (which has been working correctly all along — the user's review noted
  this is the "prior-art fix" already in the same codebase).
- Friction solver, slip-rate decomposition, ADER predictor, RK4 endpoint
  re-eval — all untouched.

## Frontera plan (NOT submitted yet — awaits user approval)
After local unit tests pass, generate
`jobs/tpv102/tpv102_200m_p1_0.01s_400rank_ader_persistent_nuc.sbatch`
copying the working sbatch's full `module load` + `LD_LIBRARY_PATH` from
the `tpv102_200m_p1_0.01s_400rank_ader_total_init.sbatch` reference.
Stop and ask the user to review before submission.
