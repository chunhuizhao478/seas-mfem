# TPV104 Frontera launch checklist

Four sbatch scripts in strict dependency order. **Each `sbatch` requires
your explicit approval per `feedback_frontera_approval.md`.** Do not
submit any of these from the agent loop — paste the `sbatch` command
yourself after reviewing the script and the state of the repo on
`/scratch2/`.

| Order | Script | Purpose | Wall | SU |
|-------|--------|---------|------|----|
| 1 | `tpv104_build_banner_sanity.sbatch` | Build + banner + `--dry-run` | ~2 min | ~0.04 |
| 2 | `tpv104_mesh_build.sbatch` | Gmsh → `tpv104_{200m,500m,1000m}.msh` | ~5 min | ~0.08 |
| 3 | `tpv104_init_stations.sbatch` | `tfinal=0` init-only, ψ_ini check at 9 stations (P2_D) | ~5 min | ~0.5 |
| 4 | `tpv104_t3_mesh_coupled.sbatch` | Mesh-coupled `tfinal=3.0 s` ADER-O2 run | ~5 h | ~2000 |

## Local verification already performed

On the developer laptop (serial, single-rank, TPV102 1000 m mesh used as
a stand-in since TPV104 geometry is byte-identical modulo header):

```
$ ./seas_tpv104_driver --mesh tpv102/mesh/tpv102_1000m.msh --tfinal 0.0 \
                       --output-dir /tmp/tpv104_init --output-prefix tpv104
...
Mesh: 215143 elements, 1 ranks
BC faces — free: 5786, fault: 1646, absorb: 5312
Fault QPs (global): 4938 (local: 4938, shared: 0)
Background: τ_strike = 40 MPa, σ_n = 120 MPa, ψ_ini(a_in) ≈ 0.5636, V_ini = 1e-16 m/s
[tpv104_driver] tfinal = 0 — init-only run. Station t=0 row written, exiting.

--- P2_D check ---
stations_ok = 9 / 9
P2_D: PASS
```

All 9 TPV104 fault stations report `ψ_init = 5.6359184263e-01`
(absolute error 2.6e-9, well under the 1e-8 threshold).

## Driver invariants

**R7-001 disclosure.** The production driver path runs one-shot
`wave.AdvanceADER` with Brent friction under the hood. The
`--friction-solver newton-stable` and `--fault-iterator substep` flags
are accepted for banner / smoke-test parity but do NOT reach the
solver dispatch — that would require exposing per-sub-step I± from
`wave_operator.inl` which is on the extreme-care no-touch list
(CLAUDE.md [C2]). For TPV104's background ψ/a ≈ 56, Brent is robust
and numerically equivalent to the stable-asinh Newton variant to
tolerance, so the physics result is unaffected. Banner + sbatch
headers + `--verify-dispatch` all consistently describe Brent +
one-shot as the runtime path.

- **Friction solver (runtime)**: Brent (log10-V, Tandem-verified),
  hard-coded via `FaultFaceFlux::EvaluateADERTotal`'s default argument.
  `SolveSlipRateNewtonStable` on `FrictionCoefficientStable` (plan §4.10
  Step 5) is compiled + unit-tested but NOT reached through the one-shot
  path.
- **Fault iterator (runtime)**: one-shot `wave.AdvanceADER`.
  `Tpv104SubStepIterator::Advance` is compiled + unit-tested (plan §4.10
  Step 7) but NOT wired into this driver.
- **ψ update (runtime)**: macro-step `UpdateStateAnalyticSlipLawSRW` (FVW
  analytic exponential relaxation, free-function dispatch) — plan §3.3.
  Plan §3.12 per-SUB-step ψ cadence is blocked on R7-001 option (a) /
  [C2] relaxation.
- **Nucleation (runtime)**: macro-step `ApplyNucleationIncremental_TPV104`
  (cumulative smoothStepIncrement) — plan §3.9 telescoping identity
  makes the total perturbation at macro-step boundaries correct.
- **V_w side-channel**: per-QP `PopulateVwSideChannel_TPV104` at init;
  read inside the ψ update.
- **DOFData init**: `InitializeFaultDOFs_TPV104` with per-QP `a(x,z)`
  and `ψ_ini` from FVW steady-state inversion.
- **Q representation**: fluctuation-Q (plan §3.10).  Pre-stress in
  `DOFData.sigma_n0 / tau*_0`; bulk Q starts at 0.
- **Station output**: `TPV104StationWriter` emits the nine SCEC trace
  columns (t, h-slip, h-rate, h-stress, v-slip, v-rate, v-stress, σ_n,
  ψ).

## Approval protocol

Before each `sbatch` the agent:
1. Shows script content + expected wall / SU.
2. Waits for your explicit "approved" or "submit".
3. Proposes the `sbatch` command — you run it.

No `ssh`, no `allocation check`, no `sbatch` from the agent loop.

## Toolchain

MFEM-side: `intel/19.1.1` + `impi/19.0.9` + `hypre/2.31.0` + `mumps/5.3`
+ `parmetis` + `petsc/3.15` + `fftw3/3.3.8` (copied verbatim from the
working `tpv102_200m_p1_0.01s_50rank_init.sbatch` per
`feedback_sbatch_modules.md`).

## Open items (non-blocking for t = 3.0 s run)

- **R7-001: Tpv104SubStepIterator integration.** Currently the t=3.0 s
  driver uses `wave.AdvanceADER` for the bulk + Riemann step, then
  updates ψ analytically via `UpdateStateAnalyticSlipLawSRW` at the
  macro-step boundary.  `Tpv104SubStepIterator` (§4.10 Step 7,
  per-sub-step accumulation) is available and unit-tested but NOT
  spliced into `wave.AdvanceADER`'s internal friction solve.  Wiring
  requires exposing per-sub-step I± from `wave_operator.inl` which is
  on the extreme-care no-touch list.  Per-sub-step cadence (plan §3.9,
  §3.12) is a Phase-3 tightening if Probe 2 / Probe 5 diffs against
  SeisSol exceed the target 1e-6 rel in §5.6.
- **R7-001 disclosure.** The `--friction-solver` and `--fault-iterator`
  CLI flags are parsed for banner parity + CLI-typo rejection but have
  NO effect on the dispatched solver / iterator until the R7-001
  wiring lands.  Banner + sbatch headers disclose this explicitly.
  Run `./seas_tpv104_driver --dry-run --verify-dispatch` to see the
  [dispatch] lines confirming Brent / one-shot / slip-SRW.
- **ParaView output.** Intentionally omitted to minimise failure
  surface on the first run.  The TPV102 driver's pv_out block is the
  reference for adding this once t=3.0 s runs are stable.
