---
title: "Implementation Plan — Strong Rate Weakening (FL=103) in the Quasi-Dynamic SEAS Driver"
subtitle: "Wire SlipLawSRWPsi + per-DOF V_w into RateStateFaultOperator so spatial_seas_driver can run the SAF QD cycle model with TPV104 strong-velocity-weakening friction"
author: "SEAS-MFEM-SAFS"
date: "2026-06-20"
---

# Plan: Strong Rate Weakening in the QD driver (`spatial_seas_driver`)

## STATUS — IMPLEMENTED 2026-06-20 (S1–S3 + the R-001 SAF loading function)

Phases S1–S3 are implemented to the one-path design and locally verified:
- **S1** `fault_geometry.hpp`: `V_w_values_` + `GetVwValues()` + the `SetRateStatePerDOF`
  fill (size-0-or-N guard).
- **S2** `rate_state_fault.hpp`: `srw_evo_` (dynamic_cast in the BP5 ctor + per-DOF
  V_w/a size verify) and the SINGLE dispatch chokepoint `StateRate_`/`StateSteady_`
  at the two vector-path sites.  **Zero edits to any friction-law file.**
- **S3** `spatial_seas_driver.cpp`: builds `SlipLawSRWPsi` (production mode) when
  `state_evolution="slip_law_strong_rate_weakening"`; the abort guard is removed.
- **R-001 SAF loading**: `MakeSAFDirichletFunc(Vp, y0)` (boundary_config.hpp) +
  `[boundary].plate_loading` ("bp5" default | "saf_recenter_y"); the driver computes
  y0 = box mid-plane from the mesh.

Verification (local): new `seas_test_spatial_seas_srw_dispatch` (serial + np=4) —
`vw_err=0` (S1), `max_srw_err=0` (operator psi-rate == `Rate_SRW` with per-DOF
V_w/a), `max_aging_err=0` (aging byte-identical), and a DECISIVE per-DOF proof:
at a V_w-sensitive operating point (psi<<0 -> frictionless-limit V in the
V_w band) `n_sensitive=27160/27936` DOFs have per-DOF != DOF-0-scalar, so the
operator's match to the per-DOF reference could NOT come from a scalar impl
(strengthened after review finding R-006 — the original psi=0.6 point gave V~0,
`Rate_SRW`~0, a vacuous test).  Regression: faultgeom parity (279372/0) +
qd-coupling-smoke (R-001 order) PASS with fresh binaries.  Review findings fixed:
R-S1 (driver clears resolver V_w for non-SRW), bp5_analytic+SRW guard, y0 via
GetBoundingBox, dead-literal.

PRE-EXISTING QD-driver crash — ROOT-CAUSED + FIXED (was NOT the SRW work).
Symptom: a full driver time-loop segfaulted at np=1 in `ComputeRHS` (via Mult ->
RK45 Step), identically for the AGING config — so pre-existing, not the SRW
dispatch.  ROOT CAUSE (lldb on a -g -O0 build: `stp d9,d14,[x8,#-0x10]`, x8=0x10
=> null `rate` data, just before `bl StateRate_`; a size-dump showed every INPUT
vector correctly sized): the driver never called `DormandPrinceRK45::Init(seas_op)`
before the loop, so the 7 RK stage vectors (`k_[i].SetSize(op.Width())`) were
size-0/null-data, and the first `Step -> Mult -> ComputeRHS` wrote `rate(0)` to
0x0.  Init succeeded because it uses a freshly-allocated `rate_temp`; only the
RK45 path used the unallocated `k_[i]` (hence law-independent).  FIX (1 line,
driver): add `ode.Init(seas_op);` before the loop.  VERIFIED: with the fix the
np=1 run completes SetInitialCondition and ENTERS the RK45 step (ComputeRHS +
friction solver run); no segfault (was a ~30 s crash).  Unblocks every QD
time-loop run (aging, SRW, BP5-parity-via-spatial_seas).
Two SEPARATE pre-existing issues remain (not the RK45 bug, not the SRW work):
(1) cg_amg is slow on DG elasticity locally (np=1 OK; np=2 init stalls) — the
plan's known AMG-on-DG perf risk; (2) the MUMPS direct-solve path segfaults
INSIDE the MUMPS library during the init Solve at np=1 (use cg_amg locally).
A loose `ksp_rtol=1e-2` produces a non-physical traction (~93 GPa) -> a clean
FRIC-GUARD NaN abort; rtol=1e-8 (the config default) avoids it.
Pending: clean full-step verification at rtol=1e-8 (slow locally) + BP5 parity +
the SAF run (Frontera).

## Why this plan exists

The SAF QD cycle config `config/safs_qd/safs_qd_500m_rssrw_v1.toml` adopts the
SeisSol `safs_seisol_v1_0_0_RSSRW` friction (SCEC FL=103 = strong rate
weakening). The friction LAW class already exists
(`friction/slip_law_srw_psi.hpp::SlipLawSRWPsi`, ported + unit-tested for the
TPV104 *dynamic* driver), and the shared config/resolver already produce
everything it needs (`state_evolution="slip_law_srw"` parses; the resolver fills
per-DOF `rs.V_w` from the `boxcar_taper` rule). **But the QD fault operator does
not drive it.**

Verified gap (read, not assumed):

- `fault/rate_state_fault.hpp` calls ONLY the base-virtual state law:
  `evolution_->SteadyState(...)` (init, lines ~307, ~323) and
  `evolution_->Rate(V, psi, Dc)` (RHS, lines ~513, ~616). It never calls
  `SlipLawSRWPsi::Rate_SRW(V,psi,L,V_w[i],a[i])` / `SteadyState_SRW(V,V_w[i],a[i])`.
- `SlipLawSRWPsi::SetProductionMode()` makes the base virtuals THROW (guard
  against exactly this silent fallthrough). So in production mode the QD operator
  would abort; without production mode it would silently use the scalar
  `V_w_default`/`a_` instead of the per-DOF depth profile — wrong physics either
  way.
- `drivers/spatial_seas_driver.cpp` builds `AgingLawPsi` unconditionally
  (~line 841); it never reads `cfg.rate_state->state_evolution`.

A guard now ABORTS any `slip_law_srw` QD config with a pointer to this plan
(`spatial_seas_driver.cpp`, just after the R-702 bp5_analytic guard), so there
is no silent aging-law fallback (CLAUDE.md). This plan removes that guard by
implementing the real wiring.

## The key simplification that makes this small + safe

In the regularized rate-and-state formulation this repo uses (psi-space, Tandem
convention), the **instantaneous fault strength is identical** for aging-law and
strong rate weakening:

```
  tau(V, psi) = sigma_n_eff * a * asinh( V/(2*V0) * exp(psi/a) )  +  eta*V
```

The strong-rate-weakening physics lives ENTIRELY in the **state evolution** (the
steady state that psi relaxes toward), not in `tau(V,psi)`:

```
  aging law :  dpsi/dt = (b*V0/Dc) * ( exp((f0 - psi)/b) - V/V0 )
  slip+SRW  :  dpsi/dt = -(V/Dc) * ( psi - psi_ss_srw(V) )
               psi_ss_srw(V) = a * ln( (2*V0/V) * sinh( f_ss(V)/a ) )
               f_ss(V)  = f_w + ( f_LV(V) - f_w ) / ( 1 + (V/V_w)^8 )^(1/8)
               f_LV(V)  = max( 0, f0 - (b - a)*ln(V/V0) )
```

Consequences (all verified against the source):

1. `friction/dieterich_ruina.hpp` (the Brent force-balance solver) is
   **UNCHANGED** — it inverts `tau(V,psi)` for V, which is law-independent. No
   touch to that extreme-care file.
2. The 4-phase init's psi-from-equilibrium step (solve `tau_traction =
   tau(V_init, psi)` for psi) is **law-independent** — also unchanged.
3. The ONLY behavioral change is: (a) which `StateEvolution` object is built, and
   (b) the operator calling the per-DOF `_SRW` overloads (which need per-DOF
   `V_w` and per-DOF `a`) instead of the base virtuals — at the init
   `SteadyState` sites and the RHS `Rate` sites.

So the surface area is: per-DOF `V_w` plumbing (resolver -> geom -> fault op,
mirroring the existing per-DOF `Dc`/`V_init` plumbing) + a state-law dispatch in
4 call sites + the driver's law selection. The aging-law path stays
byte-for-byte (the dispatch is gated on an SRW flag), which is the regression
guard for the BP5 benchmark.

---

## Phase S1 — per-DOF V_w plumbing (`fault_geometry.hpp`)

**Goal.** `FaultGeometry` carries a per-owned-DOF `V_w` array, filled by the
resolver, readable by the fault operator — exactly like `Dc_values_`,
`V_init_vec_`, `a_values_`.

**Changes (additive).**
- `fault/fault_geometry.hpp`: add `mfem::Vector V_w_values_;` (size `NumFaultDOFs()`
  = owned count) + a const getter `const mfem::Vector& V_w_per_dof() const`.
- `SetRateStatePerDOF(const spatial::RateStatePerDOFParams& rs, const mfem::Vector& init_vel_dir)`:
  also copy `rs.V_w` -> `V_w_values_` when `rs.V_w.Size() == NumFaultDOFs()`
  (the resolver fills it for SRW; for aging-law `rs.V_w` is empty/NaN — leave
  `V_w_values_` size 0). Assert `rs.V_w.Size()` is 0 or `NumFaultDOFs()`.
  MUST NOT set `params_computed_` (unchanged — R-001 order preserved).

**Acceptance.** `test_spatial_seas_faultgeom_parity` extended: with an
`slip_law_srw` rate-state block, `geom.V_w_per_dof()` matches the resolver's
`rs.V_w` element-wise to 1e-12; with `aging_law`, `V_w_per_dof().Size()==0`.

---

## Phase S2 — SRW dispatch in the QD fault operator (`rate_state_fault.hpp`, EXTREME CARE)

**Goal.** When the state law is `SlipLawSRWPsi` (production mode), the operator
evolves psi by the per-DOF SRW law; the aging-law path is unchanged.

**DESIGN DIRECTIVE (user, 2026-06-20): minimal friction-law change + ONE SRW path.**
- **Reuse the SINGLE existing `SlipLawSRWPsi`** — make **ZERO edits** to the
  friction-law files (`friction/slip_law_srw_psi.hpp`, `dieterich_ruina.hpp`,
  `state_evolution.hpp`). No new SRW law, no second SRW implementation.
- The QD path uses **ONLY the per-QP `_SRW` route** (`Rate_SRW` /
  `SteadyState_SRW`). The class's base-virtual scalar route
  (`Rate(V,psi,L)` -> `Rate_SRW(...,V_w_default_,a_)`) is **never** used in the
  QD operator: we call `SetProductionMode()` so the base virtuals THROW. So there
  is exactly **one** SRW evaluation route in the path (no scalar-vs-per-QP
  duality to trace).
- **One dispatch chokepoint:** route ALL four call sites through a SINGLE private
  helper, not four scattered `if (srw_)` branches. One place to read, one place
  to breakpoint.

**Changes (additive; aging-law byte-identical).**
- Detect SRW once at ctor: `srw_evo_ = dynamic_cast<const SlipLawSRWPsi*>(evolution_);`
  (non-null iff SRW). If non-null, cache the per-DOF `V_w` (`geom->V_w_per_dof()`)
  and per-DOF `a` (`geom->a_per_dof()`); assert both size `num_nodes_` (clear
  abort if `V_w` is empty — "slip_law_srw selected but no per-DOF V_w; check
  [friction.rate_state] V_w_default / boxcar_taper").
- Add TWO private helpers — the **only** place the aging-vs-SRW choice is made:
  ```cpp
  inline real_t StateRate(int i, real_t V, real_t psi) const {
     if (srw_evo_) return srw_evo_->Rate_SRW(V, psi, Dc_[i], Vw_[i], a_[i]);
     return evolution_->Rate(V, psi, Dc_[i]);          // unchanged aging path
  }
  inline real_t StateSteady(int i, real_t V) const {
     if (srw_evo_) return srw_evo_->SteadyState_SRW(V, Vw_[i], a_[i]);
     return evolution_->SteadyState(V, Dc_[i]);         // unchanged aging path
  }
  ```
- Replace the 2 init `SteadyState` sites + 2 RHS `Rate` sites with
  `StateSteady(i, V)` / `StateRate(i, V, psi)`. (These are the ONLY behavioural
  lines that change; everything else in the operator is untouched.)
- The Jacobian/`RateDerivative*` virtuals: the QD RK45 path uses `ComputeRHS`
  (not a Jacobian), so the derivative virtuals are not on the QD critical path;
  confirm by grep before relying on it. If any QD path calls
  `RateDerivativeV/Theta`, route them through a third helper the same way.

**Why this is safe in an extreme-care file.** The branch lives in ONE helper
whose `else` arm is the existing call verbatim, so every non-SRW driver (BP5,
TPV102, all existing rate-state) is bit-for-bit unchanged; the SRW arm only
activates when `state_evolution="slip_law_strong_rate_weakening"` is selected.
No friction-law file is edited at all.

**Acceptance.**
- New unit test `test_spatial_seas_srw_dispatch`: a 2-DOF fault with two
  different `V_w[i]`; assert the operator's psi-rate at each DOF equals
  `SlipLawSRWPsi::Rate_SRW(V, psi, Dc[i], V_w[i], a[i])` to 1e-14 (proves the
  per-DOF wiring, not just a scalar).
- Regression: ALL existing `make test` rate-state tests pass unchanged
  (aging-law path untouched); rebuild after editing the header (no header deps,
  project memory).

---

## Phase S3 — driver law selection (`spatial_seas_driver.cpp`)

**Goal.** Build the correct state law from config; remove the SRW abort guard.

**Changes.**
- Replace the unconditional `AgingLawPsi aging(...)` with a dispatch:
  ```
    StateEvolution* evo = nullptr;
    AgingLawPsi    aging(fc.b, fc.V0, fc.f0);
    SlipLawSRWPsi  srw(/*a*/rs.a(0), /*b*/fc.b, /*V0*/fc.V0, /*f0*/fc.f0,
                       /*muW=*/cfg.rate_state->f_w_default,
                       /*V_w_default=*/cfg.rate_state->V_w_default);
    if (cfg.rate_state->state_evolution == SlipLawStrongRateWeakening) {
       srw.SetProductionMode();         // per-DOF V_w/a via _SRW; base virtuals throw
       evo = &srw;
    } else {
       evo = &aging;
    }
    RateStateFaultOperator<ParMesh,2> fault_op(&geom, &friction, evo, seed, &mpi);
  ```
  (a/b/f0/V0 are the Phase-3 asserted-uniform scalars for b/f0/V0; a is per-DOF
  and reaches the operator via geom, the `rs.a(0)` here is only the class's scalar
  fallback. f_w/V_w come from the rate-state block.)
- DELETE the `state_evolution==SlipLawStrongRateWeakening` abort guard added near
  the R-702 check (it is now implemented).
- Keep `SetProductionMode()` so any accidental base-virtual call aborts loudly
  rather than silently using `V_w_default`.

**Acceptance.**
- `./seas_spatial_seas_driver --config config/safs_qd/safs_qd_500m_rssrw_v1.toml
  --dry-run` no longer aborts at the guard; prints the SRW law in the banner.
- BP5 parity (Phase 7) configs (aging-law) still pass — the SRW branch is not
  taken for them.

---

## Phase S4 — SAF QD smoke + cycle run (Frontera)

**Goal.** Run `safs_qd_500m_rssrw_v1.toml` end to end.

- `jobs/safs_qd/safs_qd_500m_rssrw_v1_dev_2hr.sbatch` (staged): build, init
  (equilibrium < 1e-6), advance >= 100 steps, finite V_max, AMG converging.
- Then `..._production.sbatch`: many cycles; verify the SAF cycle is physically
  plausible (recurrence, event slip, V_max peaks-then-decays per CLAUDE.md
  regression list).
- Caliper solver-scaling report (build `MFEM_USE_CALIPER=YES`).

---

## Acceptance gates (summary)

| Gate | What it proves | Where |
|------|----------------|-------|
| All existing `make test` pass unchanged | aging-law / BP5 / TPV path byte-identical (regression) | local |
| `test_spatial_seas_faultgeom_parity` (+V_w) | per-DOF V_w plumbed resolver->geom | local |
| `test_spatial_seas_srw_dispatch` | operator calls `_SRW` with the right per-DOF V_w/a | local |
| BP5 parity (Phase 7, aging-law) | SRW branch did not regress the benchmark | Frontera |
| SAF dev smoke | SRW QD config builds + advances stably | Frontera |

## Risk + mitigation

- **Extreme-care file `rate_state_fault.hpp`.** Mitigate: the SRW arm is purely
  additive and gated; the aging-law `else` arm is the existing code verbatim
  (regression-proved by `make test`). Single-purpose commit (no refactor mixed
  in, CLAUDE.md).
- **Per-DOF a/V_w ordering.** The owned-DOF order is the same `FaultGeometry`
  `i=0..N-1` the resolver/geom/operator already share (Phase 3); the V_w array
  rides the exact path `Dc_values_`/`V_init_vec_` already take, so it inherits
  the verified ordering.
- **Spatially-varying b/f0/V0.** Out of scope (Phase 3 asserts they are uniform);
  v1.0.0 has uniform b/f0/V0 (only a, V_w vary). If a later SAF model needs
  per-DOF b/f0/V0, that is the separate Phase-3b change.
- **No local full-mesh validation** (project memory). Mitigate: the per-DOF
  dispatch is unit-tested locally with a 2-DOF fault; production correctness is
  gated by the Frontera smoke + the unchanged BP5 benchmark.

## Estimated size

~25 lines in `rate_state_fault.hpp` (2 cached arrays + 2 dispatch helpers + 4
one-line call-site swaps; NO friction-law file edited), ~10 in
`fault_geometry.hpp` (array + getter + fill), ~12 in the driver (law dispatch,
minus the deleted guard), 1 new ~120-line unit test. Three single-purpose
commits (S1 geom, S2 operator, S3 driver+test). Reuses the single existing
`SlipLawSRWPsi` via its per-QP `_SRW` route only.
