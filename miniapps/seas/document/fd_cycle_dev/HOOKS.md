# Hooks of the fully dynamic cycle branch

Branch `fully-dynamic-cycle`, plan "Fully dynamic BP5 cycle, phase plan rev 3"
(rules M3, M5, M9, with M3 changed as below).

A hook is the smallest additive public method in a shared file. Hooks live
on this branch. Chunhui Zhao, 16 Oct 2026: "all your work should be in a
separate branch, and only commit to nitsche-dynamic-rupture when I ask you
to". So a hook is committed here, by additions only (`make -f cycle/cycle.mk
cycle-owned-check` enforces both), and reaches `nitsche-dynamic-rupture` only
when he asks, together with the rest of the branch. A hook taken from another
branch is copied byte for byte. When a hook lands, add its file to
`CYCLE_HOOK_FILES` in `cycle/cycle.mk` if it is new there, add its row below,
and add its pin to `tests/cycle/test_cycle_api_pin.cpp`.

## Landed

The branch point was `nitsche-dynamic-rupture` at `010b56e`. The branch has since
merged it at `4c97e32` (merge `7baf0a4`, 20 Oct 2026), at `b1f5446`, the QD
speedup (22 Oct 2026), and at `bc91ac91` (merge `e9a35a05`, 26 Oct 2026), to
which `nitsche-dynamic-rupture` then fast-forwarded, so H1a, H1b and H2 are on
the trunk; `cycle-owned-check` diffs against the latest merge base.

| Hook | File and method | Commit | Pinned |
|------|-----------------|--------|--------|
| H1a | `fault/fault_geometry.hpp`: `FaultGeometry` 4-argument constructor, `SetRateStatePerDOF` and `GetVwValues`, copied byte for byte from spatial-seas-driver `f6cf3f5` (`fault/fault_geometry.hpp:183`, `:678`; 154 lines, additions only) | "Fully dynamic BP5 cycle P1: hooks H1a and H1b on this branch" | yes |
| H1b | `fault/fault_geometry.hpp`: `FaultGeometry::SetInitialFieldsPerDOF(tau0, V_init, sigma_n)`, appended after H1a (44 lines, additions only) | same commit | yes |
| H2 | `domain/nitsche_wave_operator.hpp`: `NitscheWaveOperator::SetBackgroundTraction(S0)`, after `NucleationFactor` (12 lines, additions only): replaces the projected prestress by per-DOF global tractions. M9: a trial merge with `nitsche-benchmarks` (H-B1, `81759da`) is clean, and H2 sits outside the lines H-B1 changes | "Fully dynamic BP5 cycle P2: first dynamic event" | yes |
| H-2D1 | `domain/nitsche_wave_operator.hpp`: `NitscheWaveOperator::SetNormalStressModifier(m)` and `NormalStressModifier()` after `SetBackgroundTraction`, the reduction in `CloseDOF` and one `Advance` per fault step in `AdvanceFault_` (57 lines, additions only), with the interface `fault/fault_normal_stress_modifier.hpp` (new): the law sees sigma_eff less the modifier's reduction (thermal pressurization, TPV105). Branch `2d-inplane` (2D in-plane design rev 3, section 5.9; `document/2d_inplane/SHARED_CHANGES.md`, phase D part 2); null by default, every bit unchanged | "2D in-plane D2: normal-stress hook, TPV105-2D deck" (branch `2d-inplane`) | no: the cycle does not call it |
| H3 | `domain/nitsche_wave_operator.hpp`: `NitscheWaveOperator::SetConstantLoad(b)`, `ConstantLoad()` and `StiffnessMult(x, Kx)`, after `SetBackgroundTraction`; one line in `Acceleration_` (`F_ += b`) and two members (additions only): the frozen load of the absorbing side walls, `[cycle.absorbing] faces = "outer"` (plan 11.2). Global stepping only: `[cycle.lts]` with absorbing faces is refused | "Fully dynamic BP5 cycle P5: absorbing outer faces in the dynamic regime" | yes |

## Planned

| Hook | What it adds | Lands with |
|------|--------------|------------|
| H4 | optional | |

## Trunk functions the p=2 run calls (no new hooks)

The p=2 run's speed-ups and benchmark output need no shared change: they call
trunk functions as the trunk drivers do, all pinned in
`tests/cycle/test_cycle_api_pin.cpp` (section "The p=2 run").

| Deck option | Trunk functions | Trunk driver it follows |
|-------------|-----------------|-------------------------|
| `[cycle.qd].operators = "matrix"` | `SEASQuasiDynamicOperator::SetQDOperators`, `QDStageOperators::PrintSummary` | `seas_bp5_full --qd-operators matrix` |
| `[cycle.qd.green]` | `SetQDGreen`, `QDGreenConfig`, `RefreshDisplacement`, `QDGreenActive`, `PrintQDGreenCheckSummary` | `seas_bp5_full --qd-green*` |
| the FSAL guard (always) | `TSRKGetTableau` | `seas_bp5_full` (QD speed-up v2 item A) |
| `[cycle.lts] mode = "rate2"` | `BuildNitscheLtsAssignment`, `BuildLtsAwarePartition`, `ComputeLtsPartitionImbalance`, `NitscheLtsRun` | `seas_nitsche_dyn_driver --lts rate2` |
| `[cycle.benchmark]` | `ElasticityDomainOperator::SolveStage` / `Solve` (QD velocity), `GetFaceQuadrature`, `FaceQuadrature` rules, `ParMesh::FindPoints` | none (new output) |

## Shared changes merged from other tracks

A change taken whole from another branch is cherry-picked with `-x` and its
commit is listed in `CYCLE_IMPORTS` in `cycle/cycle.mk`; `cycle-owned-check`
then requires the branch to change that commit's files by exactly its patch.
Once `nitsche-dynamic-rupture` carries the change and the branch merges it,
the commit leaves `CYCLE_IMPORTS`.

| Change | Home branch, commit | Merged here | Gates re-run |
|--------|---------------------|-------------|--------------|
| Station sigma_n fix: station V and tau re-solved with the normal traction the solver used (`SEASQuasiDynamicOperator::GetFrictionNormalTraction`, BP5 station writers, `seas_driver`, `seas_bp5_full`, unit test). Chunhui Zhao, 16 Oct 2026: use the same QD path as the Vista run (job 1074099) for fully dynamic development | `gpu-qd-bp5` `ef70517` | cherry-picked; `QDPhase::FillStationFields` passes the same normal traction. Since the merge of `nitsche-dynamic-rupture` `b1f5446` the commit itself is in the trunk, so it left `CYCLE_IMPORTS` | API pin, parity, 50 steps vs `gpu-qd-bp5` `seas_bp5_full`, restart; P0 job checks time, slip, state, global vs the older golden |
| Track 2's H-B1 (`FaultPointState` 3 -> 8 scalars, `CloseDOF` trailing `slip_abs`) | `nitsche-benchmarks` `bb89a65`; reached `nitsche-dynamic-rupture` with `e4182a7` | with the trunk merge `7baf0a4` (trunk `4c97e32`) | L-10, T12, P2, P3, P4, snapshots, API pin |
