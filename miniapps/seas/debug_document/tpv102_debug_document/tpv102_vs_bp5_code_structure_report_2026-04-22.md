# TPV102 vs BP5 Code-Structure Comparison: Unnecessary Complexity and Reusable Stability

Date: 2026-04-22
Companion to: `tpv102_mfem_vs_seissol_codeflow_pepper_report_2026-04-22.md`

## Scope

BP5 (quasi-dynamic) is proven working; TPV102 (dynamic) exhibits the triangle-scale
"pepper" pattern identified in the SeisSol-vs-MFEM comparison. This document answers:

1. Which TPV102 code paths are *structurally unnecessary* (i.e. could be deleted
   or replaced by an existing BP5 construct) and therefore add bug surface?
2. Which BP5 abstractions could be reused directly or with minor extension so
   that the TPV102 path becomes as debug-stable as BP5?

The pepper report already diagnosed the likely failure *location* (fault-face
bookkeeping, interior vs shared branches, canonical-frame reconstruction,
`DOFData` side-channel). This document gives the *mechanism*: TPV102 re-builds
from scratch what BP5 already owns, and the duplication is exactly where pepper
is most naturally generated.

## Headline Numbers

| Artifact | BP5 | TPV102 |
|---|---|---|
| Driver LOC | 680 (`drivers/seas_driver.cpp`) | 1,643 (`drivers/tpv102_driver.cpp`) |
| Fault-coupling module | `domain/elasticity_operator*.inl` + `fault/rate_state_fault.hpp` | `dynamic/wave_operator.inl` (3,234 LOC) + `dynamic/fault_face_flux.cpp` (464 LOC) + `dynamic/tpv102_setup_total.hpp` (697 LOC) |
| Fault-coupling routines | 1 per-face traction integrator + 1 `ComputeRHS` | 4 near-duplicate face-flux routines (`Compute{,Shared,ADER,ADERShared}FaceFluxRHS`) |
| Per-QP state channel | state vector owned by ODE stepper | `DOFData` side-channel mutated inside flux routine |
| Cross-rank fault infrastructure | `SharedFaultCommBlock` + `FaultScatter` (owner-resolved, 197 LOC) | both ranks keep independent `DOFData` + R-101 verifier gates the drift (adds `shared_fault_elem1_on_plus_`, `shared_fault_key`, cross-rank gather/pair logic) |
| Time integrator | one `DormandPrinceRK45` (+ optional PETSc TS) | hand-rolled RK4 with 32 stage buffers + endpoint-refresh + alternate ADER-O(2/3/4) branch |
| Nucleation path | N/A (rate-state only) | `FaultQPNodalMap` + `FaultQPNucleationState` + `ApplyNucleationTotal` (writes into bulk `Q[SXY]` via shape_max-weighted nearest-nodal DOF on each side) — 200+ LOC with uniqueness asserts |
| Diagnostics compiled into production | trace/face_trace_logger, off by default | `SEAS_DIAG_FAULT_FLUX`, `--debug-qnorm`, `pv_local_*_k4` stage-4 ParaView CellData, `--paraview-bulk`, `psi_stability_warned` tripwire — all active |

## Part A — Unnecessary TPV102 Code (Highest Bug Surface First)

### A1. Four near-duplicate fault-face flux routines in `wave_operator.inl`

Line numbers: `ComputeFaceFluxRHS` (831–1322, 491 LOC), `ComputeSharedFaceFluxRHS`
(1323–1773, 450 LOC), `ComputeADERFaceFluxRHS` (1774–2106, 332 LOC),
`ComputeADERSharedFaceFluxRHS` (2107–2393, 286 LOC).

The code itself documents this (wave_operator.inl:1769):

> The duplication (vs. refactoring ComputeFaceFluxRHS into a shared core) is
> intentional: `wave_operator.inl` is on the CLAUDE.md "Files Requiring Extreme
> Care" list and the RK4 path must remain byte-identical.

What is duplicated inside each:

- Fault-QP canonical-frame reconstruction from `FaultBasisQPData::sign_flipped`
  (`can_n`, `can_t1`, `can_t2`, `T_can`, `Tinv_can`) — wave_operator.inl:1092–1110
  and 1588–1606 and 1961–1971 and 2278–2293.
- `elem1_on_plus` decision (different mechanism interior vs shared: interior
  uses `!qpd.sign_flipped`; shared uses the geometry-derived
  `shared_fault_elem1_on_plus_[sf_idx_in_fault]`).
- `Q_plus_local / Q_minus_local` selection (4 places, same formula).
- Dual-mode `EvaluateTotal` vs `Evaluate` dispatch (4 places, same `has_bulk_bg_`
  branch).
- Back-rotation `Q_imp_{plus,minus}_g` via `T_can` (4 places, same loop).
- Per-side Godunov flux + DG accumulation with `elem1_on_plus ? -1 : +1`
  (4 places, same sign algebra).

Why it matters for pepper: the pepper report's hypothesis H1 and H2 are
*literally* "interior/shared branch mismatch" and "canonical-basis /
side-ordering inconsistency." Four variants of the same 80-line block is
exactly the pattern that produces a single-branch bug manifesting as
triangle-scale speckle.

**What BP5 does instead**: one `ComputeTractionImpl` (elasticity_operator_traction.inl:488)
that has a `mode` parameter discriminating {plain, components, diagnostics}
and a single shared-face extension path. No ADER path because BP5 is
quasi-dynamic. No canonical-frame rebuild per step because the traction is
projected via the integrator's precomputed face basis.

### A2. `DOFData` as mutable per-QP side-channel populated inside `FaultFaceFlux::Evaluate`

`fault_face_flux.cpp:206–211`:

```cpp
data.slip_rate    = V_abs;
data.V1           = V1;
data.V2           = V2;
data.tau1_corr    = data.tau1_0 + tau1_corr;
data.tau2_corr    = data.tau2_0 + tau2_corr;
data.sigma_n_corr = data.sigma_n0 + sigma_n_corr;
```

`Evaluate` is called once per face QP per RK4 stage (4× per step). Because
these writes are destructive, the driver captures stage-wise snapshots into
`sr_k1..k4`, `V1_k1..k4`, `V2_k1..k4`, `t1c_k1..k4`, `t2c_k1..k4`, `snc_k1..k4`
(tpv102_driver.cpp:1000–1014, 32 `std::vector<real_t>` buffers), weighted-averages
them after stage 4, and then calls `wave.Mult(Q, k_endpoint)` a FIFTH time to
refresh `DOFData` with the endpoint values (tpv102_driver.cpp:1411–1428).

Observations:

- `psi` must explicitly be guarded `psi-pure` (assertions at fault_face_flux.cpp:83
  and 217) because the driver's coupled RK4 on `psi` integrates it *outside*
  `Evaluate`. This invariant is purely a consequence of mixing per-QP state
  into the flux routine.
- The "capture, average, endpoint-refresh" dance is R-V92-K01 debug scar
  tissue (fault-observables-at-t vs t+dt/2 time-labelling bug). It wouldn't
  exist if fault observables lived in the ODE state.

**What BP5 does instead**: fault observables (`slip_rate_`, `V_max_`) live as
private members of `RateStateFaultOperator`, written *only* by `ComputeRHS`
(the ODE RHS function). Output and the time stepper read them via
`GetSlipRate()`, `GetGlobalMaxSlipRate()`. There is no mutable per-QP state
written by the domain solve.

Why it matters for pepper: each of the 4 face-flux routines writes into
`DOFData`, and the output path reads `DOFData`. Any per-stage or per-branch
inconsistency in what gets written lands *directly* on every fault output
triangle — which is precisely the observed symptom.

### A3. Nucleation via `FaultQPNodalMap` injection into bulk `Q[SXY]`

Files: `dynamic/tpv102_setup_total.hpp` (697 LOC total), principal API
`BuildFaultQPNodalMap`, `ApplyNucleationTotal`, `FaultQPNucleationState`,
`VerifyFaultQPNodalMapUnique`, `PickNearestNodalDOF`.

Semantics: for each fault QP on each side, find the element Lagrange node
whose shape function peaks there, write `-delta/shape_max` into
`Q[SXY, elem*ndof + local_dof]` so the shape-weighted QP sum increments by
`-delta`, and maintain a per-QP `dtau_applied_{plus,minus}` tracker so
subsequent calls apply deltas instead of overwrites. Uniqueness of (elem,
dof) pairs across QPs is asserted at setup.

Why it is bug-prone:

- Both sides of each QP get an independent write; the + and − side deltas
  use *different* `shape_max` scalings. A shape_max mismatch ≡ triangle-scale
  amplitude mismatch across the fault.
- Partition-seam QPs write only on the locally-owned side; the peer rank
  writes the other side. There is no direct verification that the two
  additions actually sum to the intended QP-level `delta`.
- Called 5× per step (4 stages + endpoint-refresh) with
  `update_state=true`/`false` semantics, and each call iterates all QPs.
- Requires the uniqueness check (`VerifyFaultQPNodalMapUnique`) because the
  shape_max scaling diverges otherwise.

**What exists in the fluctuation path** (`tpv102_setup.hpp::ApplyNucleation`,
14 LOC): `dof_data[i].tau2_0 = TPV102Params::tau_ini + dtau`. One write per
fault QP, no shape_max, no MPI side-split, no uniqueness theorem. The only
reason TPV102 uses `ApplyNucleationTotal` today is the v9.3.0 total-Q
migration (I-06), which pushed pre-stress into bulk `Q` so the wave field
propagates against a nonzero background (for the absorbing-BC total-Q
semantics). The migration can be *scoped to the BC path* without dragging
nucleation out of `DOFData`.

### A4. Hand-rolled RK4 + alternate ADER branch + endpoint re-evaluation

`tpv102_driver.cpp:994–1619` implements the time loop by hand:

- 32 stage buffers (`sr_k1..4`, `V1_k1..4`, `V2_k1..4`, `t1c_k1..4`,
  `t2c_k1..4`, `snc_k1..4`, `psi_k1..4`).
- Nucleation injection at each stage time (4 `ApplyNucleationTotal` calls).
- `Q_tmp = Q + dt/2·k1`, `Q_tmp = Q + dt/2·k2`, `Q_tmp = Q + dt·k3` — explicit
  RK4 Butcher.
- Butcher-weighted averaging of `DOFData.{V1,V2,slip_rate,tau1_corr,tau2_corr,
  sigma_n_corr}` into themselves.
- A 5th `wave.Mult(Q, k_endpoint)` call to refresh `DOFData` to the endpoint.
- ParaView `pv_local_*_k4` stage-4-only snapshots shipped into production
  VTU output as a discriminator for an already-resolved debug (R-V92-E02).
- A separate 150-line ADER branch that replaces all the above with
  `wave.AdvanceADER(Q, dt, ader_order, Q_new)` and a forward-Euler `psi`
  update, preserved as "the alternative" by flag `--time-integrator rk4|ader`.
- Psi stability tripwire (`psi_stability_warned`) active in both RK4 and ADER
  branches.

**What BP5 does instead** (`seas_driver.cpp:525–656`):

```cpp
while (t < t_final && step < max_steps) {
   real_t dt;
   bool accepted = ode_solver.Step(seas_op, state, t, dt);
   if (!accepted) { continue; }
   step++;
   V_max = seas_op.GetMaxSlipRate();
   // output, eq detection, ...
}
```

One call into `DormandPrinceRK45::Step`. MPI error reduction, V-guard, PI
controller, FSAL endpoint re-use — all inside the stepper. No driver-visible
stage buffers, no endpoint re-evaluation, no alternate integrator.

### A5. Duplicated cross-rank fault infrastructure

TPV102 carries:

- `shared_fault_elem1_on_plus_` — per-shared-fault-face Boolean derived from
  comparing Elem1 centroid to face centroid along `ref_normal`
  (wave_operator.inl:389–466, ~78 LOC of ctor logic).
- `shared_fault_dof_offset_` — per-shared-fault-face → DOFData offset map
  (wave_operator.hpp:365).
- `VerifySharedFaultDOFDataConsistency` — gather-and-pair-by-FaceVertexKey
  verifier, run once after step 0, tolerance 1e-10 (wave_operator.inl in
  the 2650–3230 range).
- `shared_fault_key.hpp` (`FaceVertexKey` struct) deliberately held bit-identical
  to BP5's nested copy in `elasticity_operator.hpp` (per the shared_fault_key.hpp
  header comment).

The symmetric BP5 infrastructure:

- `SharedFaultCommBlock` + `FaultScatter` (`common/fault_scatter.hpp`, 197 LOC)
  — canonical owner resolution via `FaceVertexKey`, followed by one
  `MPI_Alltoall`-style data scatter. Only the owner rank holds authoritative
  fault state.

Why TPV102's version is higher bug surface:
1. Both ranks hold independent `DOFData` entries for each shared QP (by
   design — each rank runs its own `Evaluate`), so *by construction* they
   can drift. R-101 was added to detect drift; the drift doesn't actually
   happen under the current canonical frame, but the invariant now has to
   be policed rather than held by construction.
2. `shared_fault_elem1_on_plus_` is an ad-hoc geometric classification that
   depends on MFEM's `CalcOrtho` convention for shared faces (the
   wave_operator.inl comment documents the pitfall). `SharedFaultCommBlock`
   doesn't need this at all — it uses vertex ID comparison, which is
   rotation-invariant and mesh-topology-agnostic.

## Part B — BP5 Abstractions Reusable for TPV102 Stability

### B1. `FaultBasis` as the *sole* source of per-QP rotations (not a supplement)

Already included by TPV102, but only as input to a per-stage reconstruction
of `can_n`, `can_t1`, `can_t2`, `T_can`, `Tinv_can` (wave_operator.inl:1095–1110
and 1592–1606, recomputed every stage on every face QP).

Change: precompute the canonical rotation matrices per QP at setup time
(one `DenseMatrix T_can[nq_total]`, one `DenseMatrix Tinv_can[nq_total]`,
stored inside `WaveOperator`). The 4 face-flux routines then read from the
cache instead of rebuilding. This collapses ~300 LOC of per-stage rotation
arithmetic and eliminates the principal H2 pepper mechanism (sign-flip /
tangent swap in canonical-frame rebuild).

### B2. `FaceQuadrature` as the fault-side L2 space

`fault/face_quadrature.hpp` gives a proper H1 face basis + reference mass
matrix inverse, already validated at p=0..4 on triangles and in use by
`elasticity_operator_traction.inl` for projecting QP traction to fault
DOFs. Nucleation is by far the most natural use:

> Proper nucleation: project `delta_tau_face(x)` onto the face L2 basis,
> store the per-face coefficient vector, consume it in the flux routine
> as a *smooth field on the fault surface*.

That removes:

- `FaultQPNodalMap`, `PickNearestNodalDOF`, the `1/shape_max` scaling.
- `VerifyFaultQPNodalMapUnique` (no Voronoi-cell collision concern because
  projection is defined for any basis).
- `FaultQPNucleationState` / `dtau_applied_{plus,minus}` / the
  `update_state=true/false` `ApplyNucleationTotal` dance.
- The per-side-independent writes that are a prime A3 pepper mechanism.

`FaceQuadrature` is ~280 LOC and already in `fault/`.

### B3. `SharedFaultCommBlock` + `FaultScatter` (drop R-101 and `shared_fault_elem1_on_plus_`)

`common/fault_scatter.hpp` resolves shared-fault cross-rank ownership via
the vertex-triple key that is *already shared* (the shared_fault_key.hpp
header explicitly plans this merge). With FaultScatter in the TPV102 path:

- Only the owner rank stores / writes `DOFData` for a shared QP.
- The non-owner rank receives the evaluated `Q_imp_±` from the owner.
- `shared_fault_elem1_on_plus_` becomes redundant (ownership is canonical,
  not geometry-derived).
- R-101 `VerifySharedFaultDOFDataConsistency` becomes redundant (there is
  only ever one `DOFData` entry per shared QP).
- The second of the four face-flux routines (`ComputeSharedFaceFluxRHS`,
  450 LOC) collapses into a single dispatch helper that reads scattered
  imposed states — the interior path handles the flux math.

### B4. `RateStateFaultOperator<MeshType, 2>` as the dynamic friction operator

BP5 already has a validated vector friction solver with `SolveSlipRateVectorPsi`,
state evolution, psi-space initialization, NaN guards, per-DOF `a`, `Dc`,
`tau_pre`, `eta`. Its `ComputeRHS(traction, state, rate)` is exactly what
TPV102 needs, minus the `tau_pre_` add (which for dynamic becomes the
endpoint of the Riemann trial traction).

Proposed refactor (follows the CLAUDE.md rule "duplicate into a new file
under `dynamic/`, don't modify BP5 source"):

```
dynamic/dynamic_rate_state_fault.hpp    ← duplicate of rate_state_fault.hpp
                                          with ComputeRHS signature taking a
                                          Riemann trial traction instead of
                                          tau_pre+elastic, and a pre-stress
                                          discriminator for the total-Q path.
```

Gains:

- Single implementation of the friction solve across BP5 and TPV102.
- `FaultFaceFlux::Evaluate` reduces to the Riemann + flux-construction steps
  (Pelties Eq. 7 + Eq. 11-12) — it stops owning the friction solve. That
  eliminates A2 (`DOFData` mutable side-channel), which is the principal
  path by which stage-wise bookkeeping drives pepper-level cellwise
  artefacts into the output.
- `psi` integration moves into the RK45 state vector (as BP5 does), and
  the driver stops weaving `psi_k1..4` with `V_k1..4`.

### B5. `SEASQuasiDynamicOperator`-style top-level operator

BP5's top-level `Mult` is a 6-step physics/numerics contract that is *not*
modified by the driver (solver/seas_operator.hpp:286–320). A dynamic analogue
would be:

```
SEASDynamicOperator::Mult(const Vector &state, Vector &rate)
  1. Extract (Q_bulk, slip, psi) from state
  2. wave_op_->Mult(Q_bulk, rate_Q)       // pure elastodynamics, no DOFData
  3. wave_op_->SamplePlusMinusAtFault(Q_bulk, Q_plus, Q_minus)
  4. fault_flux_->ComputeRiemann(Q_plus, Q_minus, traction_trial)
  5. dynamic_fault_op_->ComputeRHS(traction_trial, psi, rate_slip, rate_psi)
  6. Assemble fault-flux contribution into rate_Q using imposed states
```

Gains:

- The driver's RK4/ADER loop becomes a 3-line call to
  `DormandPrinceRK45::Step(seas_dyn_op, state, t, dt)`.
- The 32 stage buffers (A4) are gone; RK stepper owns all stage arithmetic.
- Output is a pure reader of the ODE state (as in BP5), so there is no
  "endpoint-refresh" step.
- Nucleation becomes an explicit modification of `rate` (e.g.
  `rate_psi -= delta_tau_rate`) or of the `traction_trial` before the
  friction solve — no bulk-Q injection, no FaultQPNodalMap.

### B6. `DormandPrinceRK45` (with adaptive dt and MPI error reduction)

Already handles the MPI rank-consistent accept/reject
(`MPI_Allreduce(err_norm, MPI_MAX)` — CLAUDE.md calls out this as critical),
V-guard, PI controller, FSAL endpoint re-use. Drop-in replacement for
TPV102's `for (step=0; step<nsteps; step++)` fixed-dt loop, yielding:

- Automatic time-step control around nucleation (which is presently
  enforced manually by CFL + `psi_stability_warned`).
- FSAL gives the endpoint re-evaluation for free — the R-V92-K01 refresh
  becomes the next step's stage-1.

## Part C — Recommended Consolidation Order

Each item in Part B is independent; they can be landed in any order. The
order below minimizes compounded risk and tests each change against
pepper independently.

| # | Change | Pepper mechanism eliminated | Risk |
|---|---|---|---|
| 1 | B1 — cache canonical `T_can`/`Tinv_can` per QP at setup | H2 (canonical-basis/side-ordering) | Low (cache of already-computed values) |
| 2 | B3 — move shared faults to `SharedFaultCommBlock` + `FaultScatter`; delete `shared_fault_elem1_on_plus_`, R-101 verifier, `ComputeSharedFaceFluxRHS` core logic | H1 (interior/shared branch mismatch) + H4 (`DOFData` ownership) | Medium (MPI behaviour change; covered by `seas_test_bp5_parallel_smoke` pattern duplicated for dynamic) |
| 3 | B4 — dynamic `RateStateFaultOperator` in `dynamic/`; `FaultFaceFlux::Evaluate` stops writing `DOFData.{V,slip_rate,tau_corr,sigma_n_corr}` | H3 (per-side Q trace mismatch) + A2 `DOFData` side-channel | Medium-high (friction split from flux; needs unit test against EvaluateTotal per-QP) |
| 4 | B5 — top-level `SEASDynamicOperator`; driver RK loop becomes `ode_solver.Step(...)` | A4 hand-rolled RK + endpoint-refresh scar tissue | Medium (full driver restructure; ADER branch becomes an alternative stepper, not a 150-LOC in-driver branch) |
| 5 | B2 — `FaceQuadrature`-based nucleation as proper L2 projection | A3 nucleation-via-Q injection bookkeeping | Low (nucleation amplitude is small compared to wave field; quick validation) |
| 6 | B6 — `DormandPrinceRK45` as the default stepper | A4 residual (already gutted by 4) | Low (stepper is well-tested on BP5) |

After #1 and #2, the pepper report's three top hypotheses (H1, H2, H4) are
structurally impossible, regardless of what else is or isn't wrong with the
physics. That is the highest-leverage first pair of changes for the current
bug.

## Part D — Items That Look Like Unnecessary Code But Aren't

For completeness, some TPV102 complexity is load-bearing and should *not*
be consolidated:

- `GodunovFlux` and `FaultFaceFlux::ComputeTrialTraction`: dynamic-specific,
  no BP5 analogue. Keep.
- The ADER predictor-corrector (`ComputeADERTimeIntegrated`, `AdvanceADER`):
  genuine alternative integrator with O(dt^4) accuracy for the wave equation.
  Keep — but *as* an `ODESolver` under the same top-level `Mult`, not as a
  150-line branch inside the driver.
- PML (`PMLLayer`): dynamic-specific. Keep.
- Absorbing / free-surface total-Q BC variants
  (`AbsorbingTotal`, `FreeSurfaceTotal`, `FreeSurfaceGodunovTotal`): needed
  for TPV102 with background pre-stress. Keep.
- `TPV102StationWriter`, `TPV102SurfaceStationWriter`: SCEC benchmark output
  format; dynamic-specific. Keep.

## Bottom Line

The TPV102 fault-coupling path carries roughly 2× the LOC and 4× the
fault-face code paths of the BP5 equivalent, and the duplication is
concentrated in exactly the places the pepper report flagged as the
highest-risk triangle-scale artefact generators. BP5 already owns all the
abstractions needed — `FaultBasis`, `FaceQuadrature`, `FaultScatter`,
`RateStateFaultOperator`, `SEASQuasiDynamicOperator`, `DormandPrinceRK45` —
and using them (with dynamic-specific duplicates under `dynamic/` per the
project's no-modify-BP5 rule) would remove four of the five named
mechanisms that can produce the pepper pattern by construction rather than
by invariant policing.

The highest-leverage first step is B1 + B3 (precompute `T_can` per QP,
scatter shared-fault state via `SharedFaultCommBlock`). Those two together
collapse the 450-LOC shared-fault flux routine and eliminate the per-stage
canonical-frame rebuild that the pepper report identifies as H1 and H2 —
the top two ranked hypotheses.
