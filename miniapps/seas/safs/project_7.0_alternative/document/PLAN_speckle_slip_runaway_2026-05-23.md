# Implementation Plan: SAFS slip "speckle" runaway — catch, localize, bound (single-rank, post-reconcile)

> **Revision 2 (2026-05-23).** Incorporates the plan audit
> `REVIEW_PLAN_speckle_slip_runaway_2026-05-23.md` (verdict: PASS WITH FIXES). The
> headline tripwire in rev-1 was logically dead (it compared slip against its own
> kinematic ceiling); rev-2 re-anchors it to an independent physical `c_s·t` ceiling.
> Two byte-exact premises were corrected and three secondary findings resolved. See
> **§Audit resolution** for the finding-by-finding mapping.

## Overview
The cross-rank reconcile is done (`worst_rel == 0`), but `max_slip` still grows
unphysically (62875 → 105366 m) while sampled `V_max ≈ 4.67 m/s`. The companion review
`REVIEW_speckle_slip_runaway_debug_2026-05-23.md` identified the mechanism (code-confirmed
below): an intermediate sub-step velocity spike is integrated into slip but is invisible to
`V_max`, ratchets a DOF past `d_c`, and the now-permanently weakened DOF keeps slipping.
This plan adds, in dependency order, (1) a per-sub-step diagnostic trace that discriminates
the spike source, (2) an **honest** sub-step `V_max` plus a slip-plausibility tripwire
anchored to the physical shear-wave-speed ceiling so the runaway is *visible and localized*,
and (3) a finiteness/`c_s` sanity bound that flags the unphysical slip rate at its source.
Per the 2026-05-23 scoping decision, the **tension fix (R-005d), within-step non-relaxation
(R-006d), and the resolution / physical-`D_c` study (R-007d) are deferred** to a follow-up
plan gated on what the Phase-1 trace shows. (Deferred-finding IDs carry a `d` suffix to
avoid collision with the plan-audit finding IDs.)

**Decisions locked (user, 2026-05-23):**
- **Guard policy:** finiteness (`NaN`/`Inf`) always aborts; the `c_s`-bound and
  slip-plausibility checks are **env-gated and NON-FATAL by default** (log + localize the
  first offender, run proceeds so the full speckle field still reproduces). An env flag
  flips them to abort. Mirrors `SEAS_R101_NONFATAL` / `SEAS_DIAG_BLOWUP`.
- **Scope:** observability + bounds now (Phase 1/2/3); physics deferred.

## The mechanism (verified against source, 2026-05-23)
Radiation-damping LSW balance, per sub-step node `o`, per fault QP `i`:

- trial+total traction → `tau_abs = |τ_total| = sqrt(τ1_total² + τ2_total²)`
  (`tpv205_friction.hpp:106-107`);
- `μ(δ) = μ_s − (μ_s−μ_d)·min(δ/d_c, 1)`, monotone non-increasing, `δ = sqrt(slip1²+slip2²)`
  (`LSWFrictionCoefficient_TPV205`, `:53-68`); **once `δ ≥ d_c` ⇒ `μ_d` forever** (`:66`);
- `τ_str = μ(δ)·max(σ_n_total, 0)` (`:133-134`);
- `V_abs = (tau_abs − τ_str)/η_s` if `tau_abs > τ_str` else 0 — **no upper bound** (`:137-144`);
- slip integrates **every** sub-step: `d.slip1 += s.V1·dt_sub; d.slip2 += s.V2·dt_sub`
  (`tpv205_substep_iterator.cpp:122-123`) — **this is the ONLY slip-accumulation site in the
  whole code** (verified: no `slip1/slip2 +=` anywhere in `wave_operator.inl`; the substep
  iterators are the sole accumulators);
- `slip_rate = V_abs` is written **only on the last sub-step**: `if (last_sub_step)
  WriteBackState(d,s)` (`:131-134`) → `data.slip_rate = s.V_abs` (`fault_face_flux.cpp:307`),
  and **additionally overwritten on shared faces** by the macro-dt shared-fault solve
  (`wave_operator.inl:4222`) and the cross-rank reconcile (`:5157`);
- the driver's `V_max` reads `dof_data[i].slip_rate` (`spatial_dyn_driver.cpp:1891`),
  `n_rupturing` uses `slip_rate>0.5` (`:1916`), and `[DIAG-ONSET]` argmaxes by `slip_rate`
  (`:1919`).

⇒ A transient overshoot at an **intermediate** node (`o < O−1`) makes a large `V_abs`,
integrated into `slip` but unseen by `slip_rate`/`V_max` (only `o = O−1` is written) ⇒
`max_slip ≫ V_max·t`. The spike ratchets `δ` past `d_c` ⇒ the DOF sits at `μ_d` and keeps
slipping ⇒ `max_slip` grows monotonically. The barrier short-circuit (`:116-128`,
`μ_eff ≥ ½·μ_s_barrier ⇒ V=0`) means runaway DOFs are necessarily **inside** the
30 km × 15 km rupture area.

**Why the rev-1 tripwire was dead (audit R-001):** slip is the time-integral of velocity,
so `slip_i ≤ (max|V|)·t ≤ V_substep_peak·t` is a *kinematic identity*. A bound of the form
`k·V_substep_peak·t` with `k ≥ 1` is therefore satisfied by construction and can never trip
— it compares slip against its own ceiling. This is structurally the **same** failure as
the original `slip_rate`-reconcile bug: a check defined in terms of the very quantity it is
trying to validate is blind. The fix (Phase 2) anchors the bound to an **independent**
reference — the physical shear-wave speed `c_s`, which does not depend on the observed
velocities.

## Constraints
- **Byte-exact regression contract (HARD).** Every change in Phases 1–3 MUST be inert on a
  normal run: `test_tpv205_evaluate_ader_lsw_parity`, `test_tpv205_friction`,
  `test_tpv104_substep_iterator*`, `test_tpv102_*`, and the TPV102/104/205 driver regressions
  must stay byte-identical. **The invariant is "the edits are inert," NOT "the modified code
  path is unused by TPV*."** (Audit R-003: TPV205 *does* route through the modified
  `AdvanceWithSubStepStates` — `tpv205_driver.cpp` → `AdvanceADERWithSubStep` →
  `iterator.AdvanceWithSubStepStates`, `tpv205_substep_iterator.cpp:293`.) Inertness is
  guaranteed by: (a) the new `DOFData` field is numeric, zero-initialized, and read **only**
  by the SAFS monitor — never by any TPV*/BP5 output, imposed-state accumulation, or
  checkpoint; (b) all traces/bounds are env-gated (default off) or pure observability that
  does not feed back into state; (c) the only edit inside the shared helper
  `SolveLSW_TPV205` is an always-on `MFEM_VERIFY(isfinite(...))` that cannot change
  control flow on finite input. The parity regression is the gate that proves inertness.
- **Checkpoint compatibility (audit R-004, verified).** The TPV104 checkpoint
  (`io/tpv104_checkpoint.hpp`, used by the spatial driver) is **plain text, field-enumerated**
  (writes/reads exactly: `psi, slip_rate, V1, V2, slip1, slip2, tau1_nuc, tau2_nuc,
  sigma_n_nuc`) — there is no `sizeof`/`memcpy(DOFData)` serialization. Adding a field to
  `DOFData` therefore does NOT affect checkpoint/restart, and the new field MUST NOT be added
  to the checkpoint schema (it is a per-macro-step transient, recomputed every step;
  serializing it would be meaningless and would break the byte-identical V1 file format).
- **No magic constants** (CLAUDE.md): the velocity ceiling derives from the material
  shear-wave speed `c_s = sqrt(μ/ρ)` (`spatial_dyn_driver.cpp:1199` `cs_seed`; iterator via a
  new `FaultFaceFlux::cs()`), and the slip floor from the per-DOF `d_c` (`lsw_d_c`).
- **Report, don't silently work around** (CLAUDE.md): no clamp-and-continue. Default is
  log-and-localize (non-fatal); finiteness aborts; an env flag opts into abort for the bounds.
- Friction solver stays closed-form LSW (no Newton/Brent change); this plan does not touch
  the rate-and-state path.

---

## Phase 1: per-sub-step trace (discriminator) + honest sub-step `V_abs` field (R-004 substrate)

### Goal
A single env-gated run produces a per-sub-step trace at every spiking node, and each DOF
carries the max `V_abs` it saw over the current macro step's writes — across BOTH the
iterator sub-steps and the shared-fault macro solve — so the monitor can stop being blind.

### Files to Modify
- `dynamic/fault_face_flux.hpp` — add one numeric field to `DOFData`, adjacent to `slip_rate`
  (`:50`): `real_t slip_rate_substep_max = 0;` (max `|V|` over the current macro step's
  writes; zero-init; read only by the SAFS monitor; **not** checkpointed). Update the struct
  doc-comment to state it is transient and intentionally excluded from the checkpoint schema.
- `dynamic/tpv205_substep_iterator.cpp` — in `AdvanceWithSubStepStates` (the predictor path,
  `:262-...`), inside the sub-step loop after the `StepOneQP_` call (`:379-380`):
  1. running max: `d.slip_rate_substep_max = std::max(d.slip_rate_substep_max, s.V_abs);`
     on **every** sub-step `o` (covers the intermediate spike that `slip_rate` misses). The
     per-macro-step reset is done by the driver (below), so no `o==0` special-case.
  2. the env-gated trace (below). Do **not** touch the no-callback `Advance` overload
     (`:137`) or `StepOneQP_`.
- `dynamic/wave_operator.inl` — at the shared-fault `WriteBackState` site (`:4222`), append:
  `fdata_qq.slip_rate_substep_max = std::max(fdata_qq.slip_rate_substep_max, states[qq].V_abs);`
  so shared QPs (whose `slip_rate` is the macro-dt solve, not the iterator's) also contribute
  to the honest max. (Audit R-005.)
- `drivers/spatial_dyn_driver.cpp` — in the per-step loop, **before**
  `AdvanceADERWithSubStep_Spatial` (`:1880`), reset the field for all DOFs:
  `for (int i=0;i<num_fault_total;++i) dof_data[i].slip_rate_substep_max = 0.0;`
  This makes it a clean per-macro-step max regardless of which path (iterator or shared)
  writes a given DOF. (Audit R-005.)

### Detailed Requirements
1. **Honest sub-step max** (R-004 substrate). The driver reset + iterator/shared max-updates
   yield a per-macro-step max over every `V_abs` written this step. Cost: one compare per QP
   per write; negligible. It does not change any imposed-state accumulation.
2. **Trace** (the discriminator), self-selecting on the spike, **runtime-gated only** (no
   preprocessor wrapper — audit R-007). `s`, `d`, `i`, `o`, `O`, `t_substep_end` are all in
   scope; `mu_eff`/`tau_strength` are reconstructed here so `StepOneQP_` stays untouched:
   ```cpp
   // tpv205_substep_iterator.cpp, in AdvanceWithSubStepStates, just after StepOneQP_.
   static const bool   slip_diag  = [] {
      const char *e = std::getenv("SEAS_DIAG_SLIP"); return e && e[0] && e[0] != '0'; }();
   static const double slip_v_thr = [] {
      const char *e = std::getenv("SEAS_DIAG_SLIP_VTHR"); return e ? std::atof(e) : 10.0; }();
   if (slip_diag && s.V_abs > slip_v_thr)
   {
      const real_t delta  = std::sqrt(d.slip1*d.slip1 + d.slip2*d.slip2);
      const real_t mu_eff = LSWFrictionCoefficient_TPV205(delta, d.lsw_mu_s,
                                                          d.lsw_mu_d, d.lsw_d_c);  // mirrors StepOneQP_:100
      const real_t sn_pos = std::max<real_t>(s.sigma_n_total, 0.0);
      std::fprintf(stderr,
         "[SLIP] t=%.6e o=%d/%d qp=%d last=%d V_abs=%+.6e tau_abs=%+.6e "
         "tau_str=%+.6e sigma_n_tot=%+.6e sigma_n_pos=%+.6e delta=%+.6e "
         "mu_eff=%+.6e d_c=%+.6e\n",
         t_substep_end, o, O, i, (o == O-1) ? 1 : 0,
         s.V_abs, s.Theta, mu_eff*sn_pos, s.sigma_n_total, sn_pos,
         delta, mu_eff, d.lsw_d_c);
      std::fflush(stderr);
   }
   ```
   `s.Theta == tau_abs` (set at `StepOneQP_:94`). The `last=` field flags whether the spike is
   on the node `slip_rate` samples (the R-004 blindness). The `mu_eff` reconstruction uses the
   identical function and arguments as `StepOneQP_:100`, so it is faithful to the iterator's
   solve **by construction** (audit R-006 — see the forced-rupture note in §Edge Cases).
3. **Discrimination key** (record in the follow-up doc after the run):
   - `tau_abs` spikes at intermediate `o`, `sigma_n_tot` stays > 0 ⇒ **predictor overshoot /
     under-resolution** → follow-up R-006d/R-007d.
   - `sigma_n_tot < 0` (tensile), `sigma_n_pos → 0`, `tau_str → 0` ⇒ **tension free-slide**
     → follow-up R-005d.
   - neither, but `delta` still grows ⇒ accumulation/`dt_sub` bug → re-open Phase 3.

### Edge Cases to Handle
- `O == 1` (single sub-step): the driver reset + single max gives `slip_rate_substep_max ==
  slip_rate`; honest `V_max` equals the current `V_max` (no regression).
- **Shared QPs (audit R-005):** the iterator accumulates their slip and maxes the field; the
  shared-fault macro solve overwrites their `slip_rate` and (via the `:4222` edit) also maxes
  the field; the cross-rank reconcile may set `slip_rate` to the boss's value. Because the
  monitor (Phase 2) reports `max(slip_rate, slip_rate_substep_max)`, under-reporting is
  impossible regardless of ordering or reconcile.
- **Forced rupture (audit R-006):** SAFS uses `gradual_overstress`, NOT forced rupture —
  `T_forced_rupture` stays at the `1.0e9` sentinel (`spatial_dyn_driver.cpp:1210` supplies the
  "never forced" sentinel; the wave-op dispatch resolves to plain `EvaluateADER_LSW`), and the
  iterator's `StepOneQP_` has no forced-rupture branch. So `LSWFrictionCoefficient_TPV205` is
  the correct coefficient and the trace mirrors the iterator. The implementer MUST verify
  (AC below) that the traced config has `T_forced_rupture == 1.0e9`; if a config ever activates
  forced rupture, the iterator/wave-op `μ_eff` consistency is a separate pre-existing question,
  out of scope here.
- Log flooding: the `> slip_v_thr` gate (default 10 m/s ≫ physical ~5 m/s) limits prints to
  genuine spikes. If still excessive, add a per-step print cap (note it; do not silently drop).

### Acceptance Criteria
- [ ] `make seas_spatial_dyn_driver` builds; `SEAS_DIAG_SLIP=1` Dc2 run to `t≈1 s` emits at
      least one `[SLIP] ... last=0 ...` line with `V_abs > 10·V_max_macro` (proves the spike is
      on an intermediate node, invisible to `slip_rate`).
- [ ] With `SEAS_DIAG_SLIP` unset: byte-identical output to pre-change on
      `test_tpv205_evaluate_ader_lsw_parity` and the TPV205 driver regression.
- [ ] On a serial/interior fixture, `dof_data[i].slip_rate_substep_max ≥ dof_data[i].slip_rate`
      for every interior DOF every step. (On shared QPs the monitor uses `max(...)`, so the
      invariant is asserted as `max(slip_rate, slip_rate_substep_max) ≥ slip_rate` — trivially.)
- [ ] Checkpoint round-trip test (`test_*tpv104_checkpoint*` or the spatial restart smoke)
      still passes: the new field is absent from the serialized schema and restart is
      unaffected.
- [ ] Verified: the SAFS config under trace has `T_forced_rupture == 1.0e9` (grep the resolved
      `T_forced_s` vector / `spatial_setup.hpp` path), confirming the plain-LSW dispatch.

### Dependencies
- Depends on: nothing. Required by: Phase 2 (reads `slip_rate_substep_max`).

---

## Phase 2: honest blow-up monitor + slip-plausibility tripwire anchored to `c_s` (R-004 + R-003-of-review)

### Goal
The driver's blow-up monitor reports the sub-step-aware `V_max` (so it can no longer be fooled
by a last-node-only sample), and a per-DOF tripwire — anchored to the **physical `c_s·t`
ceiling**, independent of the observed velocities — names the DOF whose accumulated slip is
physically impossible. Non-fatal by default.

### Files to Modify
- `drivers/spatial_dyn_driver.cpp`, the per-step monitor (`:1888-1953`):
  1. **R-004 honest `V_max`:** reduce over `V_honest_i = max(dof_data[i].slip_rate,
     dof_data[i].slip_rate_substep_max)` to get `V_substep_max_step`; track
     `V_substep_peak_global`. Print it alongside the existing `V_max` in `[DIAG]` so the gap is
     visible. Use `V_substep_max_step` for the `>10` onset trigger and for `[DIAG-ONSET]`, and
     add a **slip-magnitude** argmax (`sqrt(slip1²+slip2²)`) distinct from the existing
     `slip_rate` argmax — the runaway DOF is the max-slip one, not the max-`slip_rate` one.
  2. **Slip-plausibility tripwire (audit R-001 fix — anchored to `c_s`, NOT to peak V):**
     ```cpp
     // env gates (parse once): SEAS_DIAG_SLIP enables; SEAS_DIAG_SLIP_ABORT=1 makes fatal.
     // PHYSICAL ceiling: cumulative slip cannot exceed (a multiple of) the shear-wave
     // speed integrated over elapsed time. c_s = cs_seed = sqrt(mu/rho) (driver scope, :1199).
     // This reference is INDEPENDENT of the observed velocities (the rev-1 bug was anchoring
     // to V_substep_peak, which is slip's own kinematic ceiling -> never trips).
     const real_t k_slip   = [] { const char*e=std::getenv("SEAS_SLIP_PLAUS_KCS");
                                  return e?std::atof(e):3.0; }();
     const real_t k_floor  = [] { const char*e=std::getenv("SEAS_SLIP_PLAUS_KFLOOR");
                                  return e?std::atof(e):10.0; }();
     long long n_implausible_local = 0; int worst_dof = -1; real_t worst_slip = 0;
     for (int i = 0; i < num_fault_total; ++i) {
        const DOFData &d = dof_data[i];
        const real_t slip_i = std::sqrt(d.slip1*d.slip1 + d.slip2*d.slip2);
        // per-DOF physical ceiling; d.lsw_d_c is the LSW critical slip distance.
        const real_t slip_bound = k_slip * cs_seed * t + k_floor * d.lsw_d_c;
        if (slip_i > slip_bound) { n_implausible_local++;
           if (slip_i > worst_slip) { worst_slip = slip_i; worst_dof = i; } }
     }
     // reduce n_implausible (SUM) + MAXLOC(worst_slip) across ranks; rank 0 logs the count;
     // the owning rank logs worst_dof xyz + full DOFData state (reuse [DIAG-ONSET] fields).
     // if SEAS_DIAG_SLIP_ABORT: MFEM_VERIFY(n_implausible_g == 0, "...names worst_dof, slip,
     //                                       bound, cs_seed, t...");
     ```
     - **Why this fires on the real runaway:** with `cs_seed ≈ 3000 m/s`, `k_slip=3`, `t≈1 s`,
       the ceiling is `≈ 9e3 m` (plus `~10·d_c`). The observed `1e5 m` runaway exceeds it by an
       order of magnitude → trips. Physical event slip `~4 m ≪ 9e3 m` → never false-trips.
       (Contrast rev-1: `4·V_peak·t` with the velocities that integrated to `1e5 m` makes the
       bound `≫ 1e5 m`, so it stayed silent on the exact failure it was meant to catch.)
     - **Floor (audit R-002):** `k_floor·d_c` is a physical absolute (`d_c ~ 0.4 m` ⇒ floor
       `~4 m`), NOT the rev-1 `V_plate·tfinal ≈ 1e-9 m`. The formula has no division by `t`, so
       the floor's only role is the `t→0` baseline; the `c_s·t` term dominates almost
       immediately. `k_floor` chosen so the floor sits just above physical event slip yet far
       below runaway.
     - **Heterogeneous-`c_s` caveat:** `cs_seed` is the homogeneous reference
       (`sqrt(μ_const/ρ_const)`); SAFS varies `c_s` per region via `heterogeneous_material`. For
       an order-of-magnitude ceiling this is robust (the runaway exceeds any reasonable
       `c_s·t`). A per-DOF refinement (`Zs/ρ`) is possible but deferred; note it inline.
     - The implausible-DOF **count** answers the review's open question (handful of ratcheted
       DOFs vs spreading region).

### Edge Cases to Handle
- Early time (`t` small): `slip_bound = k_floor·d_c > 0` prevents false trips; assert
  `slip_bound > 0`.
- `cs_seed <= 0` (mis-config): `MFEM_VERIFY(cs_seed > 0, ...)` once before the loop — never use
  a non-positive ceiling.
- The reductions are collective — every rank calls `MPI_Allreduce` unconditionally even with
  `num_fault_total == 0` (follow the existing `V_max` reduction at `:1894-1897`).

### Acceptance Criteria
- [ ] On the pre-fix Dc2 run, `[DIAG]` shows `V_substep_max ≫ V_max` at the steps where
      `max_slip` jumps (quantifies R-004 blindness).
- [ ] On the pre-fix Dc2 run, the tripwire reports a **non-zero** implausible-DOF count and
      names the worst DOF's xyz once `max_slip` exceeds `k_slip·cs_seed·t` — i.e. it actually
      fires on the runaway (the rev-1 bound did not). Unit-level proof: a synthetic DOF with
      `slip = 1e5 m` at `t = 1 s`, `cs_seed = 3000`, `k_slip = 3` ⇒ `slip_bound ≈ 9e3 m` ⇒
      tripped; a DOF with `slip = 4 m` ⇒ not tripped.
- [ ] With `SEAS_DIAG_SLIP_ABORT=1` the tripwire aborts naming that DOF; default (unset) it
      logs and the run completes.
- [ ] With `SEAS_DIAG_SLIP` unset: TPV102/104/205 driver regressions byte-identical (the new
      reductions/branches are gated or pure observability).

### Dependencies
- Depends on: Phase 1 (`slip_rate_substep_max`). Required by: Phase 3 (shares the env gate).

---

## Phase 3: bound the unbounded `V_abs` at its source (review R-002)

### Goal
A non-finite slip rate aborts always; an unphysical (`≫ c_s`) slip rate is detected and
localized at the friction solve rather than silently integrated.

### Files to Modify
- `dynamic/fault_face_flux.hpp` — add a public accessor to `FaultFaceFlux`:
  `real_t cs() const { return cs_; }` (next to the ctor; `cs_` exists at `:447`).
- `dynamic/tpv205_friction.hpp` `SolveLSW_TPV205` — add **only** an always-on finiteness guard
  right after `V_abs` is set (`:139-144`); no signature change, no control-flow change on
  finite input:
  ```cpp
  MFEM_VERIFY(std::isfinite(V_abs),
              "SolveLSW_TPV205: non-finite V_abs (tau_abs=" << tau_abs
              << " tau_strength=" << tau_strength << " eta_s=" << eta_s << ")");
  ```
- `dynamic/tpv205_substep_iterator.cpp` — in `AdvanceWithSubStepStates`, co-located with the
  Phase-1 trace, the env-gated **non-fatal** `c_s` sanity check (default off):
  ```cpp
  static const double vsane_kcs = [] {
     const char *e = std::getenv("SEAS_SLIP_VSANE_KCS"); return e ? std::atof(e) : 3.0; }();
  static const bool   vsane_abort = [] {
     const char *e = std::getenv("SEAS_DIAG_SLIP_ABORT"); return e && e[0] && e[0]!='0'; }();
  const real_t v_sane_max = vsane_kcs * flux_.cs();      // c_s-derived, NOT a magic constant
  if (slip_diag && s.V_abs > v_sane_max) {
     std::fprintf(stderr, "[SLIP-VSANE] t=%.6e o=%d qp=%d V_abs=%+.6e > %g*c_s=%+.6e\n",
                  t_substep_end, o, i, s.V_abs, vsane_kcs, v_sane_max);
     std::fflush(stderr);
     if (vsane_abort) {
        MFEM_ABORT("SolveLSW_TPV205 V_abs=" << s.V_abs << " m/s exceeds " << vsane_kcs
                   << "*c_s=" << v_sane_max << " at qp " << i << " — transient overshoot.");
     }
  }
  ```

### Detailed Requirements
1. `v_sane_max = k·c_s` with `k` env-tunable (default 3.0). A slip RATE approaching `c_s` is
   unphysical even for supershear; the observed runaway implies a time-averaged `V ≫ c_s` at
   the worst DOF, so `k∈[1,5]` all catch it. (Note: this is the per-sub-step `V_abs` rate bound;
   Phase 2's `c_s·t` bound is the cumulative-slip ceiling — same physical anchor, two scales.)
2. **Heterogeneous caveat (documented, optional refinement):** a per-DOF `c_s` could be
   recovered as `Zs/ρ`, but `ρ` is not on `DOFData`; `flux_.cs()` is the homogeneous reference.
   Defer per-DOF refinement; note it inline.
3. The finiteness `MFEM_VERIFY` is always on because a `NaN`/`Inf` `V_abs` is never acceptable
   and cannot occur on a correct run (so the regression is unaffected).

### Edge Cases to Handle
- `flux_.cs() <= 0` (mis-initialized): `MFEM_ASSERT(flux_.cs() > 0)` once at iterator entry.
- Legitimate near-`c_s` rupture (supershear test, if any): default `k=3` avoids tripping; the
  bound is a tripwire, not a clamp, so even a trip only logs.

### Acceptance Criteria
- [ ] `test_tpv205_friction`: a `SolveLSW_TPV205` call with `tau_abs = 50·tau_strength` returns
      a finite (large) `V_abs` and does NOT abort (the in-solver guard is finiteness-only); and
      a `σ_n_total > 0` byte-identity check proves the in-solver edit is inert.
- [ ] A unit test feeding a manufactured non-finite input trips the `MFEM_VERIFY` (if
      non-finite cannot be manufactured through the public API, document that and rely on the
      caller-side `c_s` test).
- [ ] Iterator-level test: a two-sub-step fixture with an injected `tau_abs` overshoot on
      sub-step 0 → `[SLIP-VSANE]` logs once; with `SEAS_DIAG_SLIP_ABORT=1` it aborts naming the
      qp; with the env unset, no log and byte-identical imposed-state output.
- [ ] TPV205 parity regression byte-identical with all `SEAS_*` env unset.

### Dependencies
- Depends on: Phase 1 (trace gate + caller hook). Required by: nothing.

---

## Audit resolution (REVIEW_PLAN_speckle_slip_runaway_2026-05-23.md)
- **R-001 CRITICAL** (tripwire can never fire): FIXED — Phase 2 anchors the bound to the
  independent `k·c_s·t` ceiling instead of `k·V_substep_peak·t`. The §Overview now explains the
  structural-blindness equivalence to the original reconcile bug.
- **R-002 MODERATE** (1e-9 floor): FIXED — floor is `k_floor·d_c` (physical), and the bound has
  no `t`-division so the floor is only the `t→0` baseline.
- **R-003 MODERATE** (false byte-exact premise): FIXED — Constraints now state the invariant as
  "edits are inert," verified that TPV205 routes through the modified overload
  (`tpv205_substep_iterator.cpp:293`), with the parity regression as the gate.
- **R-004 MODERATE** (layout invariant / checkpoint): RESOLVED — verified the TPV104 checkpoint
  is text/field-enumerated (no `memcpy(DOFData)`); the new field is safe, **not** serialized,
  and a checkpoint round-trip AC is added.
- **R-005 MODERATE [POSSIBLE]** (shared-QP blind spot): FIXED — monitor reports
  `max(slip_rate, slip_rate_substep_max)`; the field is also max-updated at the shared
  `WriteBackState` (`wave_operator.inl:4222`) and reset per macro step in the driver. AC adjusted
  for shared QPs.
- **R-006 [POSSIBLE]** (forced-rupture `μ_eff`): RESOLVED — verified SAFS uses plain LSW
  (`T_forced_rupture == 1e9` sentinel, `spatial_dyn_driver.cpp:1210`; iterator `StepOneQP_` has
  no forced branch). The trace mirrors `StepOneQP_:100` exactly; a verification AC is added.
- **R-007 LOW** (dead `#ifndef NDEBUG_SEAS_SLIP`): FIXED — removed; the trace is runtime
  env-gated only.

---

## Testing Strategy
- **Unit (per phase, in existing files):**
  - `tests/unit/test_tpv205_friction.cpp` — review-R-002: 50× overshoot returns finite `V_abs`;
    (if constructible) non-finite input trips `MFEM_VERIFY`; **and** a `σ_n_total > 0`
    byte-identity check proving the in-solver edit is inert.
  - A new iterator/driver fixture (model on `test_tpv104_substep_iterator.cpp`) covering
    R-004/R-005/Phase-2/Phase-3: two sub-steps, sub-step 0 `V=100`, sub-step 1 `V=1`; assert
    `slip_rate == 1` (unchanged), `slip_rate_substep_max == 100`, slip `≈ (100+1)·dt_sub`;
    `max(slip_rate, slip_rate_substep_max) == 100`; with the env gate on, `[SLIP]`/`[SLIP-VSANE]`
    fire; the `c_s·t` tripwire trips for a synthetic `slip=1e5 @ t=1` and does not for `slip=4`.
  - Wire the new fixture into a `make` target and into `test`/`test-parallel` (match the
    reconcile test's R-004 wiring precedent).
- **Byte-exact regression (gate for every phase):** `test_tpv205_evaluate_ader_lsw_parity`,
  `test_tpv205_friction`, `test_tpv104_substep_iterator*`, `test_tpv102_*`, TPV102/104/205 driver
  smoke + the checkpoint round-trip — all byte-identical / passing with `SEAS_*` env unset. This
  is the primary guard that the observability/bounds are inert in production.
- **Frontera (user runs):** the Dc2 sbatch with `SEAS_DIAG_SLIP=1` to `t≈1 s`: capture the
  `[SLIP]` trace (discriminator), the `V_substep_max ≫ V_max` gap (R-004), the implausible-DOF
  count + worst-DOF xyz (Phase-2 tripwire — now actually fires), and the first `[SLIP-VSANE]`
  (Phase-3). One run classifies the spike source and feeds the deferred follow-up.

## Risk Assessment
- **Byte-exact breakage** is the dominant risk. Mitigated by keeping `SolveLSW_TPV205` edits to
  a finiteness-only `MFEM_VERIFY`, keeping `StepOneQP_` untouched, putting all
  observability/bounds in callers behind env gates, the new field unread by any TPV* output,
  and gating every phase on the parity + checkpoint regressions.
- **Tripwire still mis-anchored** — mitigated by the explicit unit AC that a `1e5 m @ t=1`
  synthetic DOF trips and a `4 m` one does not, and by anchoring to `c_s` (independent of the
  observed velocity), per the audit meta-note.
- **Shared-QP/reconcile interaction** — mitigated by reporting `max(slip_rate,
  slip_rate_substep_max)` (under-reporting impossible) plus the shared-path max-update and the
  driver per-step reset.
- **Heterogeneous `c_s`** makes a single global ceiling approximate — acceptable for a tripwire
  (runaway exceeds any reasonable `c_s·t` by an order of magnitude); per-DOF refinement deferred.
- **Log flooding** under `SEAS_DIAG_SLIP` — mitigated by the `V_abs > vthr` self-selection; add
  a per-step cap if needed (note it, don't silently drop).
- **Misattribution** — this plan makes the runaway observable and bounded; it does NOT fix the
  root excitation. The Phase-1 trace MUST classify the source before the deferred follow-up
  touches physics (R-005d) or mesh (R-007d). Do not ship "refine the mesh" until the trace rules
  out tension (R-005d) and the accumulation path is bounded (Phases 2–3).

---

## Deferred follow-up (separate plan, gated on the Phase-1 trace verdict)
Not implemented here by decision (2026-05-23). Open once the `[SLIP]` trace classifies the spike:
- **R-005d tension free-slide** (`tpv205_friction.hpp:133-134`): if the trace shows
  `σ_n_total < 0` at the spike, evaluate a no-opening response (zero tangential `V` under
  tension). **Must** first check the SCEC TPV205 spec and Tandem/SeisSol convention (CLAUDE.md:
  always refer to Tandem); the branch must be inert for `σ_n > 0` to preserve the TPV205
  byte-exact regression. Changes physics — requires explicit sign-off.
- **R-006d within-step bulk non-relaxation** (`AdvanceWithSubStepStates`): all sub-step `Q̃[o]`
  are predicted from start-of-macro-step `Q` and do not see the slip accumulated by earlier
  sub-steps. If the trace shows sustained high `tau_abs` across nodes while `δ` grows, the
  remedy is smaller macro `dt` or a corrected predictor — investigation, not a one-liner.
- **R-007d resolution + physical `D_c`** study: `{h, h/2, h/4}` on the nucleation patch + the
  physical `D_c = 1.0` config (Dc2 widens `L_nuc` but leaves `h` fixed, so it is not a
  substitute for refining `h`). Decide refine-vs-ill-posed from whether `max_slip(t)` and the
  speckle count converge as `h → 0`.
