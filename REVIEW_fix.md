# Fix Report — Phase 1 (RS primitives) round-3 review

Companion to `REVIEW.md` round-3 (Phase 1). The three round-3 findings
were all LOW hardening gaps; all are now fixed and covered by tests.
Prior findings (R-001 / R-008 / R-009) were already verified RESOLVED by
the reviewer and are unchanged here.

## Fixes applied

- **R-012 (LOW) — FIXED.** `SeedEquilibriumPsi_RS`
  (`dynamic/spatial_setup.hpp`) now validates the `rs` per-DOF vectors up
  front, before the loop:
  `MFEM_VERIFY(rs.a.Size() >= nd && rs.Dc.Size() >= nd && rs.V_init.Size() >= nd, ...)`
  (the three vectors the seed reads), matching the sibling
  `InitializeFaultDOFs_Spatial_RS`. A too-short `rs` now aborts with a
  clear message instead of an OOB read.

- **R-013 (LOW) — FIXED.** The `nuc_callback` overload of
  `Tpv102SubStepIterator::AdvanceWithSubStepStates`
  (`dynamic/tpv102_substep_iterator.cpp`) now rejects an empty
  `std::function` up front (next to the other arg checks) with a
  `std::runtime_error` ("pass a no-op `[](real_t,real_t){}` to opt out"),
  so the hot loop can no longer throw an opaque `std::bad_function_call`.
  Chose the throw form (over `if (nuc_callback) ...` skip) to match the
  overload's existing `throw std::runtime_error(...)` validation style and
  to avoid silently swallowing a mis-wired callback (CLAUDE.md).

- **R-014 [POSSIBLE] (LOW) — FIXED.** `SeedEquilibriumPsi_RS` now branches
  on the *effective* traction `tau_eff = tau0 − d.eta_s · V_init > 0`
  instead of `tau0 > 0`. When the radiation damping dominates
  (`η·V_init ≥ tau0`) the locked steady state `f0 + b·ln(V0/V_init)`
  absorbs the case, so `InitialStatePsi`'s `log(negative)` → NaN path is
  unreachable. The priority SAFS run (`V_init = 1e-9`) and the existing
  fixture (`V_init = 0.010`) keep `tau_eff > 0`, so their seeds are
  byte-unchanged.

## New / extended tests

- `seas_test_seed_equilibrium_psi_rs` — now **14/14**:
  - **S5 (R-012):** `dof_data.size()=2` with `rs` vectors size 1 →
    `RunInChild_` confirms abort; control with matched sizes → finite ψ.
  - **S6 (R-014):** `tau0=1 Pa`, `V_init=1 m/s` (η·V_init ≫ tau0) →
    `std::isfinite(ψ)` and `ψ == f0 + b·ln(V0/V_init)` (locked branch).
    Pre-fix this DOF seeded NaN.
- `seas_test_tpv102_nuc_callback_parity` — now **15/15**:
  - **P3 (R-013):** a default-constructed `std::function` callback throws
    `std::runtime_error` (NOT `std::bad_function_call`).

## Regression

- The plain `AdvanceWithSubStepStates` overload and the `Advance` method
  are untouched (standalone TPV102 unaffected). The R-013 guard lives only
  in the new callback overload.
- S1–S4 / P1a / P1b / P2 from round-2 still pass unchanged; the R-014
  branch change is a no-op for every `tau_eff > 0` DOF (all existing
  cases).

## Ready for Re-Review: YES (Phase 1 hardening complete; proceeding to Phase 2)

---

# Fix Report — Round-4 review (Phase 1 + Phase 2, R-015…R-018)

All four Round-4 findings were LOW and latent (off the wired path); all
are now fixed and covered by tests.

## Fixes applied

- **R-015 (LOW) — FIXED.** `RateStateAgingFrictionIterator`
  (`dynamic/friction_iterator.hpp`) now `= delete`s copy/move ctors +
  assignment.  `it_` binds `const AgingLawPsi&` to the member `law_`, so a
  member-wise copy/move would dangle; the adapter is only ever held via
  `unique_ptr` (in-place construct), so forbidding copy/move turns the
  footgun into a compile error.  Compile-time guard: `static_assert`s in
  the factory test.

- **R-016 (LOW) — FIXED.** `LswFrictionIterator::Advance` now rejects an
  empty `nuc_callback` with `std::runtime_error` before forwarding —
  symmetric with the Tpv102 (R-013) guard.  The owned
  `Tpv205SubStepIterator` (byte-exact oracle) is untouched; the guard
  lives in the adapter, so the wired LSW path no longer throws an opaque
  `std::bad_function_call`.

- **R-017 (LOW) — FIXED.** `SeedEquilibriumPsi_RS` now guards
  `d.sigma_n0 > 0` (alongside the V_init / size / V0 guards).
  `InitialStatePsi` divides by `sigma_n0`, so `==0 -> inf`, `<0 -> NaN` ψ.
  Verified `ResolveRateState` already enforces `sigma_n_eff > 0`
  (`spatial_friction.cpp:1073`), so this stays LOW — a defensive tripwire,
  not a wired-path bug.  Also updated the `@brief` to narrate the
  `tau_eff ≤ 0` (damping-dominated) branch (the doc nit the reviewer
  flagged).

- **R-018 (LOW / QUALITY) — ADDRESSED via option (c).** The plan mandated
  the verbatim copy of the plain `AdvanceWithSubStepStates` body (keep it
  as the byte-exact oracle), and the reviewer noted options (a)/(b) are
  what Phase 5's `SubStepIteratorBase` will subsume — so restructuring now
  would conflict with the plan and double Phase-5 work.  Chose the
  reviewer's option (c): broadened the parity fixture.  New **P4** drives
  O=3 sub-steps over 3 macro-steps with BOTH dip and strike pre-stress +
  trial traction, asserting the callback overload stays bit-identical to
  the plain overload at every step — so a future divergence outside the
  narrow P1a fixture is caught.

## New / extended tests

- `seas_test_seed_equilibrium_psi_rs` — now **16/16**:
  **S7 (R-017):** `sigma_n0 == 0` and `< 0` → `RunInChild_` confirms abort.
- `seas_test_tpv102_nuc_callback_parity` — now **18/18**:
  **P4 (R-018):** broadened bit-identity fixture (3 macro-steps, O=3,
  dip+strike, nonzero pre-stress) + meaningfulness checks (dip slip
  accumulates, nucleation fires).
- `seas_test_friction_iterator_factory` — now **10/10** + 4 compile-time
  `static_assert`s: **F5 (R-016)** LSW empty callback → `runtime_error`;
  **static_assert (R-015)** adapter is non-copy/move-constructible/assignable.

## Regression

- Driver recompiles after the `friction_iterator.hpp` changes and full
  driver links (R-015/R-016 are header-only, additive guards/deletions).
- All prior Round-2/3 test cases still pass unchanged.

## Ready for Re-Review: YES (R-015…R-018 resolved; proceeding to Phase 3)

---

# Fix Report — Phase 3 (SAFS + rate-and-state run) — verification + RS print-derived

This pass **verified** the Phase-3 implementation (driver RS branch, D3.1,
`spatial_friction.cpp` guards, the new RS config + sbatch, and the 3 Phase-3
unit tests) by building and running it for the first time — the round-4
`REVIEW.md` audit built/ran only Phase 1+2; the Phase-3 source was written but
never compiled or executed. It then closed the one gap that surfaced: the
SAFS-RS sbatch passes `--print-derived`, and `PrintDerivedAndCheck` is LSW-only.

## Phase-3 plan requirements — status (all PRESENT in source; now verified)

| Req | What | Where | Verified |
|-----|------|-------|----------|
| 1 | remove the `is_lsw` reject | driver:715-720 (reject gone; law verify + `is_lsw` kept) | dispatch run reaches `law: rate_state` |
| 2 | friction-law tag from `is_lsw` | driver:975 | `[verify-dispatch] friction law: RateAndState` |
| 3 | RS resolver branch (R-002 gate + R-003 zero `PorePressureSpec{}`) | driver:1167-1197 | `test_resolve_rate_state_guards` R-003 8/8; dry-run no abort |
| 4 | DOF-init branch + `SeedEquilibriumPsi_RS` | driver:1265-1285 | dry-run constructs `dof_data`, seeds ψ, exits clean |
| 5 | `MakeFrictionIterator` dispatch | driver:1803-1815 | `[verify-dispatch] friction iter: RateStateAgingFrictionIterator (Tpv102 aging, Brent)` |
| 6 | D3.1 factory negation `-sigma_xy_pa` | driver:1144-1149 | `test_constant_tensor_sign` 27/27 |
| 7 | SAFS-RS config (a<b, V_0=1e-6, gradual_overstress) | config/...rate_state...Dc2.toml | dry-run + print-derived (VW, V_0 guard passes) |
| 8 | `test_constant_tensor_sign` (Dc2/Dc8/Dc10+RS, parametrized) | tests/unit | 27/27 (base+Dc2+Dc8+Dc10+rate_state) |
| 9 | RS sbatch (8N/400r dev) | jobs/safs/...ratestate...sbatch | present |
| 10 | ParaView RS state channel = ψ | driver:1881-1890 | `test_paraview_state_channel_rs` 7/7 |
| 11 | reject per-DOF b/f_0/V_0 (R-006) | spatial_friction.cpp:1013-1024 | `test_resolve_rate_state_guards` R-006 |
| 12 | relax a<b validator (R-011) | spatial_friction.cpp:1077-1083 | `test_resolve_rate_state_guards` R-011 |

**Build + run evidence (`conda activate mfem-dev`):**
- `seas_test_resolve_rate_state_guards` → **8/8** (R-003 no-double-PP, R-006 per-DOF b reject, R-011 a>b allowed + a≤0 abort).
- `seas_test_constant_tensor_sign` → **27/27** (D3.1 golden bit-identity over base/Dc2/Dc8/Dc10/rate_state).
- `seas_test_paraview_state_channel_rs` → **7/7** (R-005 RS=ψ, LSW byte-identical).
- `seas_spatial_dyn_driver` **compiles + links**; `mpirun -np 8 … --dry-run --verify-dispatch` on the RS config reaches `law: rate_state`, reports `RateAndState (aging)` + `RateStateAgingFrictionIterator (Tpv102 aging, Brent)` + `gradual_overstress` + scalar flux, then exits clean — the **plan's Phase-3 acceptance criterion is met**.

## DEVIATION FROM PLAN — RS `--print-derived` path (user-approved)

**What:** The Phase-3 sbatch `jobs/safs/spatial_dyn_resolution_Dc2_ratestate_8N_400r_dev_2hr_safs.sbatch:178` and the new RS config's documented run command both pass `--print-derived`. The driver's `--print-derived` block called `PrintDerivedAndCheck(pd_cfg, lsw, …)` unconditionally, and that function is **LSW-only** — it dereferences `lsw.mu_s/mu_d/d_c` (L_nuc, barrier detection, the TRIGGER/LOCKED/STRESS-DROP gates). For an RS run `lsw` is empty, so it aborts at `spatial_print_derived.cpp:132` (`MFEM_VERIFY(lsw.mu_s.Size()==N)` → `MPI_ABORT`).

**Why this is a deviation:** the Phase-3 plan (rev 6) **never mentions `--print-derived`**. Its acceptance criterion is `--dry-run --verify-dispatch` only — which passes without any change. There is therefore no plan spec for RS-derived diagnostics; implementing them is *added scope*. The user was asked and chose **"Full RS-aware gate"** (vs. gating it off or dropping the flag), and asked that the deviation be reported here and the RS print path be handed to the reviewer for correctness checking.

**What was implemented (all RS physics grounded in this repo's own kernels — no Tandem, per `feedback_safs_dynamic_no_tandem`):**
- New overload `spatial::PrintDerivedAndCheckRS(cfg, rs, …)` in `dynamic/spatial_print_derived.{hpp,cpp}` (LSW overload untouched/byte-identical).
- RS nucleation length **`L_nuc = μ·Dc / ((b−a)·σ_n_eff)`** — the direct rate-and-state analog of the LSW `μ·d_c/((μ_s−μ_d)·σ_n)` at `spatial_print_derived.cpp:153` (raw, no Day/Andrews prefactor, per `feedback_lnuc_no_prefactor`), over velocity-WEAKENING (b>a) DOFs only.
- Steady-state friction **`f_ss(V) = a·asinh[(V/2V_0)·exp(ψ_ss/a)]`, `ψ_ss = f_0 + b·ln(V_0/V)`** — the repo's regularized Dieterich-Ruina coefficient (`friction/dieterich_ruina.hpp:272-295`) at the seeded steady-state ψ (`friction/state_evolution.hpp:196-200`), with the identical large-(ψ/a) overflow guard. Verified: f_ss(V_init)=0.634539 matches both the hand formula f_0+(a−b)ln(V_init/V_0) and the live 8-rank run.
- RS gate: hard **FAIL** (abort iff `abort_on_failure`) on the three physically unambiguous cases — (a) all DOFs velocity-strengthening (b≤a; RS analog of "all barriers"), (b) nucleation patch empty, (c) in-patch entirely velocity-strengthening. **WARN-only** for resolution (L_nuc/h_min<10), the overstress trigger (|τ_pre+δτ|−f_ss·σ_n<0), and the outside-asperity ratio (|τ_pre|/(f_ss·σ_n)≥1) — these are heuristics, not hard equilibrium violations, and RS nucleation is more forgiving than the LSW static-yield gate.
- Driver `--print-derived` block now branches `is_lsw ? PrintDerivedAndCheck : PrintDerivedAndCheckRS` (driver:1357-1402).

**Verification:**
- `seas_test_spatial_print_derived` → **33/33** (24 original LSW assertions unchanged + 5 new RS cases T-RS01…T-RS05: f_ss formula, L_nuc formula, all-VS gate, VW-nucleation PASS, VS-patch diagnosis).
- Live `mpirun -np 8 … --dry-run --print-derived` on the RS config: prints RS-grounded `[derived]` block (VW=127113/127113, f_ss=0.634539, L_nuc 132–421 km with L_nuc/h_min 670–2139, overshoot +16.8 MPa "drives acceleration"), gate **PASS**, clean exit — was aborting before.
- LSW no-regression: `--print-derived` on the LSW Dc2 config still prints the LSW L_nuc (1200–3831 m) and `PASS` (the LSW overload is byte-unchanged; only wrapped in `if (is_lsw)`).

**Hand-off:** the RS print path (`PrintDerivedAndCheckRS` + the driver branch + T-RS01…T-RS05) is the part requiring an independent correctness review — particularly the choice of L_nuc formula, the WARN-vs-FAIL gate severities, and the outside-asperity-ratio semantics under the equilibrium ψ seed.

## Files changed (Phase 3 verification pass)

- [modified] `dynamic/spatial_print_derived.hpp` — declare `PrintDerivedAndCheckRS` (RS overload).
- [modified] `dynamic/spatial_print_derived.cpp` — add `SteadyStateFrictionRS` helper + `PrintDerivedAndCheckRS`.
- [modified] `drivers/spatial_dyn_driver.cpp` — branch the `--print-derived` block (`is_lsw` → LSW; else RS).
- [modified] `tests/unit/test_spatial_print_derived.cpp` — add `RunPrinterRS` + T-RS01…T-RS05 (uses the existing Makefile target, no new target).

## Ready for Re-Review: YES (Phase 3 implemented + verified; RS print-derived deviation flagged for correctness review)

---

# Fix Report — Round 5 review (Phase 3, R-019…R-024)

The Round-5 audit ran the wired RS path and found my R-006/R-011 edits had
broken a pre-existing baseline test (`make test` red — R-019, CRITICAL), plus
one MODERATE and four LOW items. All six are addressed below.

## R-019 CRITICAL — resolver baseline test broke (`make test` SIGABRT) — FIXED

`tests/unit/test_spatial_friction_resolver.cpp` pinned the *old* RS contract
that R-006/R-011 changed; it is in `make test` and aborted at R-2 (exit 134).
Updated it to the new contract — the only edit is to stop it contradicting the
new positive coverage in `test_resolve_rate_state_guards.cpp`:
- **R-2** ("RS box override"): dropped the per-DOF `r.b = 0.020` (R-006 now
  rejects it → would SIGABRT outside a fork) and kept the `r.a = 0.005`
  override; the in-slab assertion now expects `b == 0.015` (the default, since
  per-DOF `b` is rejected).
- **R-5** ("validator"): inverted from "a≥b aborts" to **"a≥b does NOT abort"**
  (R-011), and added a second forked case asserting **`a≤0` still aborts**
  (positivity preserved).

**Verified:** `./seas_test_spatial_friction_resolver` → **90/90, exit 0**
(was exit 134). Grep confirms only two test files call `ResolveRateState`
(this one + the new guards test, whose `rule.b=0.02` is inside a fork that
*expects* the abort) — no other test depends on the changed behaviour.

## R-020 MODERATE — 8-rank RS parallel ψ correctness unvalidated/unguarded — runtime guard LANDED; test status documented

The reviewer's gate is satisfied by **either** the np=2 ψ-consistency test
**or** a runtime warning. Landed the warning (`drivers/spatial_dyn_driver.cpp`,
right after `SetFaultFrictionLaw`): at `law=rate_state && nprocs>1`, rank 0
prints a WARNING that shared (rank-seam) fault QPs use end-of-step ψ (1st-order)
and the R-004 gate is unvalidated → "Treat parallel RS results as PRELIMINARY."
**Verified:** fires on an `mpirun -np 2 … --dry-run` of the RS config.

On the np=2 ψ-consistency *test* (plan R-004, two assertions):
- **Part (a) — cross-rank ψ bit-identity post-reconcile — ALREADY COVERED.** The
  existing `tests/unit/test_shared_fault_reconcile_cross_rank.cpp` runs the np=2
  shared-fault **rate-and-state** leg, injects a cross-rank seed, and asserts
  `VerifySharedFaultDOFDataConsistency` `worst_rel == 0` over a 9-field payload
  that **includes ψ** (`FIELD_NAMES[5]="psi"`, `wave_operator.inl:6102`). That
  is exactly part (a), already green.
- **Part (b) — interior-vs-shared ψ temporal drift — DEFERRED (with reason).** A
  one-macro-step version would be **misleading**: the sub-step iterator advances
  `d.psi` identically for shared and interior QPs in a single step
  (`tpv102_substep_iterator.cpp:348`); the 1st-order degradation is in the
  shared-QP **flux** (end-of-step ψ fed to the macro-step traction), a
  *multi-step flux-feedback* effect, not a single-step ψ difference. A faithful
  part (b) needs the multi-step time loop (np=1 interior twin vs np=2 shared),
  which the plan itself keeps "open until validated." The runtime warning is the
  reviewer's accepted interim close-out; the multi-step test remains the real
  close-out before any 8-rank RS result is presented as validated.

## R-021 LOW — sub-critical (L_nuc ≫ patch) RS setup passed the gate silently — FIXED

Added a WARN-only block to `PrintDerivedAndCheckRS` (mirrors the existing
`L_nuc/h_min` resolution WARN): when `nuc.enabled` and the in-patch `L_nuc_min`
exceeds `max(radius_dip, radius_strike)`, it warns the patch is sub-critical and
may not nucleate. **Verified:** the shipped Dc=2.0 smoke config (L_nuc 132 km vs
4 km patch) now prints `WARNING: in-patch L_nuc (min 132008 m) exceeds the
nucleation patch radius (4000 m) — … sub-critical …`. WARN-only per the
reviewer's "at minimum a WARN" (the smoke Dc is a documented cluster-tuning
knob), so the structural PASS still prints.

## R-022 LOW — seed test V_init at 0.010, not production 1e-9 — FIXED

`tests/unit/test_seed_equilibrium_psi_rs.cpp`: `kVInit` 0.010 → **1.0e-9**
(production quasi-static rate) + updated the fixture comment. **Verified:**
**16/16**, including S4(b) "t=0 solve sits at V_init (|V−V_init| ≤ 1e-6·V_init)"
— the seed→Brent round-trip recovers 1e-9 at the production ψ/a regime, so no
bracket-widening finding to surface.

## R-023 LOW — `IFrictionIterator::WaveOpLaw()` dead on the driver path — FIXED

Moving the iterator construction ~900 lines earlier is unsafe (it needs
`fault_flux`/`rs`, resolved late) and dropping `WaveOpLaw()` breaks the factory
test that asserts it. Instead added an `MFEM_VERIFY` right after the factory
builds the iterator: `substep_iterator.WaveOpLaw()` must equal
`is_lsw ? LSW : RateAndState`. This makes the accessor **live on the production
path** as a single-source-of-truth cross-check rather than dead surface area; it
trips only on a factory/`is_lsw` mismatch. (Runs in a full execution; the
`--dry-run` exits before iterator construction.)

## R-024 LOW — ParaView-state-channel test mirrored the driver instead of exercising it — FIXED

Extracted the selection into `dynamic/fault_state_channel.hpp`
(`inline real_t FaultStateChannelValue(bool is_lsw, const DOFData&)` — LSW
friction coefficient vs RS ψ). The driver's ParaView writer and the R-005 test
now **both call it**, so the test guards the real code. Added the header to
`DYNAMIC_HEADERS` (it was an explicit list, not a glob) so incremental builds
track it. **Verified:** `seas_test_paraview_state_channel_rs` **7/7**.

## Verification summary (Round-5 fixes)

| Test | Result |
|------|--------|
| `seas_test_spatial_friction_resolver` (R-019) | **90/90, exit 0** (was 134/SIGABRT) |
| `seas_test_seed_equilibrium_psi_rs` (R-022) | **16/16** (S4(b) round-trip at V_init=1e-9) |
| `seas_test_paraview_state_channel_rs` (R-024) | **7/7** (via shared helper) |
| `seas_test_spatial_print_derived` (R-021) | **33/33** |
| `mpirun -np 2 … --dry-run --print-derived` (RS) | R-020 + R-021 warnings fire; clean exit |

## Files changed (Round-5)

- [modified] `tests/unit/test_spatial_friction_resolver.cpp` — R-019: R-2/R-5 to the new R-006/R-011 contract.
- [modified] `drivers/spatial_dyn_driver.cpp` — R-020 np>1 RS warning; R-023 `WaveOpLaw()` consistency assert; R-024 use `FaultStateChannelValue` (+ include).
- [modified] `dynamic/spatial_print_derived.cpp` — R-021 sub-critical-patch WARN in `PrintDerivedAndCheckRS`.
- [modified] `tests/unit/test_seed_equilibrium_psi_rs.cpp` — R-022 `kVInit` 0.010 → 1e-9.
- [new] `dynamic/fault_state_channel.hpp` — R-024 shared `FaultStateChannelValue`.
- [modified] `tests/unit/test_paraview_state_channel_rs.cpp` — R-024 call the shared helper.
- [modified] `miniapps/seas/Makefile` — R-024 add `fault_state_channel.hpp` to `DYNAMIC_HEADERS`.

## Unresolved / deferred
- **R-020 part (b)** (interior-vs-shared ψ temporal-drift test): deferred — a
  one-step version is misleading (see above); the multi-step close-out is gated
  behind the runtime warning, which is landed. The full np=2 ψ test file was
  **not** created; part (a) is covered by `test_shared_fault_reconcile_cross_rank.cpp`.

## Ready for Re-Review: YES (R-019 critical fixed — baseline green; R-020 guard landed; R-021…R-024 resolved)

---

# Fix Report — Round 6 review (R-025 + carried R-020(b)/R-021) + RS+SAFS run deliverable

Round 6 verified the Round-5 fixes (all RESOLVED), ran the RS time loop at np=2
(V_max=1e-9 at step 0 — the equilibrium-ψ seed is correct at runtime), and
raised one new readiness finding (R-025) plus two carried gates. This pass
addresses them and ships the RS+SAFS run deliverable.

## R-025 MODERATE — RS run window vs reflections — FIXED (correctly; reviewer's `use_pml` suggestion is a no-op — see deviation)

**DEVIATION FROM THE SUGGESTED FIX (skill rule 5).** The review suggested
`use_pml=true`. I verified that **`use_pml` is a no-op in `seas_spatial_dyn_driver`**:
`WaveOperator` has `SetPML()`/`pml_layer_` (`wave_operator.hpp:246,821`), but the
spatial driver constructs `WaveOperator(pmesh, order, λ, μ, ρ, bc)` (`:933`) and
**never calls `SetPML`** — `cfg.numerics.use_pml` is consumed only by the CLI
override (`:624`), the banner (`:743`), and the reflection-warning gate
(`:958`). So `use_pml=true` would **silence the legitimate warning without
absorbing anything** — strictly worse than the status quo. The box exterior does
use first-order **absorbing BCs** (`bc.absorbing_attrs={103,104}`, `:828`), which
partially (not perfectly) damp reflections.

**Correct fix applied:** bound `tfinal` to the box reflection-free window
`min_box_dim/cp_max ≈ 6.94 s` (for the 500 m mesh + cp=5996). Set
`tfinal = "6s"` in `spatial_friction_rate_state_safs_projected_stress_resolution_Dc2.toml`
(was `"100s"`) with a comment documenting the no-op-PML caveat and the absorbing
BCs. **Verified:** the RS dry-run now reports `tfinal: 6 s` with **no** reflection
warning. Wiring a real PML layer into the spatial driver is a recommended
follow-up (out of scope — feature addition).

## R-020(b) MODERATE — psi-explicit np=2 consistency test — CREATED (and the overstated claim corrected)

The Round-5 report claimed part (a) was "covered by
`test_shared_fault_reconcile_cross_rank.cpp`." The review correctly flagged this
as **overstated** — that test asserts an aggregate `worst_rel` and has zero
literal `psi` references. **Fix:** created
`tests/unit/test_shared_fault_rs_psi_consistency_np2.cpp` (MPI np=2), wired into
the Makefile (`seas_test_shared_fault_rs_psi_consistency_np2` +
`test-shared-fault-rs-psi-consistency-np2`; excluded from the default `test:`
umbrella like the other np=2 tests). It builds the 2-tet fault partitioned
across the seam (fault face SHARED), sets up TPV102 aging RS, **perturbs
`d.psi` on rank 0 only**, runs one `AdvanceADER` step (whose
`ComputeADERSharedFaceFluxRHS` reconcile broadcasts payload[5]=ψ), then **gathers
`d.psi` from both ranks and asserts bit-identity** — psi-EXPLICIT.

**Verified (`mpirun -np 2`): 3/3** — (1) PRE-step ranks DIFFER (perturbation
took; non-vacuous), (2) POST-step ψ BIT-IDENTICAL across ranks (R-020 part a),
(3) `VerifySharedFaultDOFDataConsistency` worst_rel==0. This closes part (a).
**Part (b)** (interior-vs-shared temporal-drift bound) remains the deferred
close-out — it needs the multi-step time loop + an np=1 interior twin; the
driver's "Treat parallel RS results as PRELIMINARY" warning (R-020) gates it.

## R-021 physics — sub-critical patch — ADDRESSED via the new TPV102-default config

The Dc=2.0 m smoke config is left as a documented smoke (sub-critical WARN now
prints). The **new deliverable config uses default TPV102 Dc=0.02 m**, giving
`L_nuc = μ·Dc/((b−a)·σ_n) ≈ 1.65–5.27 km` vs the 4 km patch — i.e. the patch is
(barely) **super-critical**, and the sub-critical WARN correctly does **not**
fire (verified on the np=2 dry-run). An advisory `L_nuc/h_min min ≈ 8.4 < 10`
(marginal process-zone resolution at the most-stressed DOF) prints — informative;
the user tunes Dc/mesh later.

## RS+SAFS run deliverable (default TPV102 friction params)

Created, per the request (default TPV102 params, to be tuned later):
- `safs/project_7.0_alternative/config/spatial_friction_rate_state_safs_projected_stress_tpv102_defaults.toml`
  — SAFS mesh/stress/material/nucleation + `[friction.rate_state]` = TPV102
  defaults (`a=0.008<b=0.012`, `Dc=0.02`, `V_0=1e-6`, `f_0=0.6`, `V_init=1e-12`,
  `eta=auto`); `use_pml=false` (no-op note); `tfinal="6s"` (R-025).
- `jobs/safs/spatial_dyn_ratestate_tpv102def_8N_400r_dev_2hr_safs.sbatch`
  — 8N/400r dev sbatch defaulting to the new config, `SAFS_TFINAL=6s`, the
  R-020/R-025 caveats in the header, `--print-derived` (RS-aware path).

**Verified (np=2 `--dry-run --print-derived --verify-dispatch`):** `law:
rate_state`, no reflection warning, fully velocity-weakening (127074/127074),
`f_ss(V_init)=0.655262`, `L_nuc 1.65–5.27 km` (not sub-critical), overshoot
+15.8 MPa "drives acceleration", outside ratio 0.978<1, dispatch =
`RateStateAgingFrictionIterator (Tpv102 aging, Brent)` + scalar flux, `PASS`,
clean exit.

## Verification summary (Round-6)

| Check | Result |
|-------|--------|
| `seas_test_shared_fault_rs_psi_consistency_np2` (np=2) | **3/3** (psi-explicit cross-rank bit-identity) |
| New TPV102-defaults config `--dry-run --print-derived` (np=2) | RS dispatch OK; super-critical; PASS; clean exit |
| Existing Dc2 RS config `--dry-run` (tfinal=6s) | no reflection warning; clean exit |

## Files changed / created (Round-6)

- [modified] `safs/.../config/spatial_friction_rate_state_safs_projected_stress_resolution_Dc2.toml` — R-025: `tfinal` 100s→6s + no-op-PML note.
- [new] `tests/unit/test_shared_fault_rs_psi_consistency_np2.cpp` — R-020(b) psi-explicit np=2 test.
- [modified] `miniapps/seas/Makefile` — R-020(b): SRC/OBJ vars, link target, obj rule, `test-shared-fault-rs-psi-consistency-np2` runner.
- [new] `safs/.../config/spatial_friction_rate_state_safs_projected_stress_tpv102_defaults.toml` — RS+SAFS run config (TPV102 defaults).
- [new] `jobs/safs/spatial_dyn_ratestate_tpv102def_8N_400r_dev_2hr_safs.sbatch` — RS+SAFS run sbatch.

## Unresolved / deferred
- **R-020 part (b)** (interior-vs-shared ψ temporal-drift) — deferred; needs the
  multi-step time loop + np=1 twin. The np>1 RS path stays PRELIMINARY (runtime
  warning + the new part-(a) test gate it).
- **PML not wired** in the spatial driver — recommended follow-up; meanwhile
  `tfinal` is bounded by the reflection window and absorbing BCs (103/104) damp
  partially.
- **R-021 fine physics** — the TPV102-default config is super-critical but the
  prestress/patch/Dc are the user's to tune for a specific rupture; `L_nuc/h_min`
  is marginal (~8.4) at the peak DOF.

## Ready for Re-Review: YES (R-025 fixed via tfinal bound — use_pml is a no-op, documented; R-020(b) psi-explicit np=2 test passes 3/3; RS+SAFS TOML + sbatch shipped and validated)

---

# Round 8 fix report — depth-dependent RS a(z)/b(z) audit (R-028, R-029)

## Summary
- Findings addressed: 2 of 2 (R-028 MODERATE, R-029 LOW).
- Files modified: `dynamic/fault_face_flux.hpp`, `dynamic/spatial_setup.hpp`,
  `drivers/spatial_dyn_driver.cpp`, `tests/unit/test_seed_equilibrium_psi_rs.cpp`.
- Tests added: 1 (`S9_unset_b_aborts`, +2 assertions).
- Suite: PASS — seed 21/21, depth_profile 23/23, resolver 98/98, guards 10/10,
  parity 18/18 (byte-exact), ader_tpv102_smoke 4/4 (byte-exact), np2 reconcile 3/3;
  shipped-config dry-run PASS, R-029 WARN correctly silent.

## Changes
1. **[R-028]** `DOFData.b` default `0.0` → `std::numeric_limits<real_t>::quiet_NaN()`
   (`fault_face_flux.hpp`), with the comment corrected: a forgotten `b` now
   propagates NaN through the aging-law `UpdateStateAnalytic` (psi→NaN, loud),
   whereas `0.0` was SILENT-WRONG (psi snaps to f0 for psi<f0, no NaN). Added a
   per-DOF tripwire `MFEM_VERIFY(std::isfinite(d.b) && d.b > 0.0, …)` in
   `SeedEquilibriumPsi_RS` (`spatial_setup.hpp`). Updated the four seed-test
   fixtures (`MakeFixture`, S5 control, S7, S8) to set `d.b` (post-init state) so
   the intended guard fires; added `S9_unset_b_aborts`. All production RS init
   paths already set `d.b` (TPV102=`TPV102Params::b`, SAFS-RS=`rs.b(i)`), so the
   byte-exact gates are unchanged.
2. **[R-029]** Added a depth-profile WARN in the driver `--print-derived` summary.
   **Deviation (documented):** the reviewer suggested "WARN when the two CSV depth
   ranges differ"; I made it **mesh-aware** — WARN only when the ranges differ AND
   `fault_max_depth > min(a_max, amb_max)` (the actual condition where `b=a−(a−b)`
   mixes a flat-clamped curve with a varying one). The coarse "ranges differ" WARN
   would false-alarm the shipped config (CSVs to 52.7/60 km but the 16.5 km fault
   never reaches the clamped region), training users to ignore WARNs. Verified
   silent on the shipped config; fires when a (future deeper) fault exceeds the
   shallower CSV.

## New tests
- `S9_unset_b_aborts` (test_seed_equilibrium_psi_rs) — covers R-028: `DOFData.b`
  default-constructs as NaN (not 0.0); an unset `d.b` aborts the seed.

## Unresolved / deferred (unchanged from Round 6/7)
- R-020 part (b) interior-vs-shared ψ temporal drift; PML not wired; L_nuc/h_min
  marginal in the shallow VW zone (WARN-only, user's to tune).

## Ready for Re-Review: YES
