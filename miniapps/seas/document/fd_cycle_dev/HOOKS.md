# Hooks of the fully dynamic cycle branch

Branch `fully-dynamic-cycle`, plan "Fully dynamic BP5 cycle, phase plan rev 3"
(rules M3, M5, M9).

A hook is the smallest additive public method in a shared file. It is
committed first on its home branch with Chunhui Zhao's typed go, then merged
here, so the same commit is in both histories. A hook taken from another
branch is copied byte for byte. When a hook lands, add its row below and its
pin to `tests/cycle/test_cycle_api_pin.cpp`.

## Landed

None. The branch point is `nitsche-dynamic-rupture` at `010b56e`.

| Hook | File and method | Home branch, commit | Merged here | Pinned |
|------|-----------------|---------------------|-------------|--------|

## Planned

| Hook | What it adds | Lands with |
|------|--------------|------------|
| H1a | `FaultGeometry` 4-argument constructor and `SetRateStatePerDOF`, copied from spatial-seas-driver `f6cf3f5` (`fault/fault_geometry.hpp:183`, `:678`) | P1 |
| H1b | `FaultGeometry::SetInitialFieldsPerDOF(tau0, V_init, sigma_n)` | P1 |
| H2 | `NitscheWaveOperator::SetBackgroundTraction(S0)`, after track 2's H-B1 (M9) | P2 |
| H3 | `NitscheWaveOperator::SetConstantLoad` (stage B absorbing side walls only) | P5, if stage B |
| H4 | optional | |

## Shared changes merged from other tracks

| Change | Home branch, commit | Merged here | Gates re-run |
|--------|---------------------|-------------|--------------|
| Track 2's H-B1 (`FaultPointState` 3 -> 8 scalars, `CloseDOF` trailing `slip_abs`) | `nitsche-dynamic-rupture`, not landed yet | | L-00, T12, P0/P1 BP5 gates, API pin |
