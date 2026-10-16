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

The branch point is `nitsche-dynamic-rupture` at `010b56e`.

| Hook | File and method | Commit | Pinned |
|------|-----------------|--------|--------|
| H1a | `fault/fault_geometry.hpp`: `FaultGeometry` 4-argument constructor, `SetRateStatePerDOF` and `GetVwValues`, copied byte for byte from spatial-seas-driver `f6cf3f5` (`fault/fault_geometry.hpp:183`, `:678`; 154 lines, additions only) | "Fully dynamic BP5 cycle P1: hooks H1a and H1b on this branch" | yes |
| H1b | `fault/fault_geometry.hpp`: `FaultGeometry::SetInitialFieldsPerDOF(tau0, V_init, sigma_n)`, appended after H1a (44 lines, additions only) | same commit | yes |

## Planned

| Hook | What it adds | Lands with |
|------|--------------|------------|
| H2 | `NitscheWaveOperator::SetBackgroundTraction(S0)`, after track 2's H-B1 (M9) | P2 |
| H3 | `NitscheWaveOperator::SetConstantLoad` (stage B absorbing side walls only) | P5, if stage B |
| H4 | optional | |

## Shared changes merged from other tracks

| Change | Home branch, commit | Merged here | Gates re-run |
|--------|---------------------|-------------|--------------|
| Track 2's H-B1 (`FaultPointState` 3 -> 8 scalars, `CloseDOF` trailing `slip_abs`) | `nitsche-benchmarks` `bb89a65`; reaches `nitsche-dynamic-rupture` only when track 2 finishes and Chunhui Zhao asks | | L-00, T12, P0/P1 BP5 gates, API pin |
