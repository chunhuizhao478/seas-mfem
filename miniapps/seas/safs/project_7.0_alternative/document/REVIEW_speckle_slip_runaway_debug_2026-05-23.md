# Code Review / Debugging Plan: SAFS slip "speckle" runaway (single-rank, post-reconcile) — 2026-05-23

## Review Scope
- **Symptom:** post-reconcile (`worst_rel = 0`, `DIVERGED = 0`), `max_slip = 62875 → 105366 m
  and growing` while sampled `V_max ≈ 4.67 m/s` (physical slip should be ~4 m). Both ranks hold
  the *same* unphysical slip ⇒ single-rank cause (issue B is closed; this is the B.2.5
  discriminator firing).
- **Plan:** `PLAN_shared_fault_reconcile_fix_2026-05-23.md` (B.2.5 predicted this), CLAUDE.md
  ("No artificial tau floor"; sign/normal-stress conventions).
- **Files reviewed:** `dynamic/tpv205_substep_iterator.cpp` (`StepOneQP_`, `AdvanceWithSubStepStates`),
  `dynamic/tpv205_friction.hpp` (`SolveLSW_TPV205`, `LSWFrictionCoefficient_TPV205`),
  `dynamic/fault_face_flux.cpp` (`WriteBackState`), `drivers/spatial_dyn_driver.cpp`
  (`V_max`/`max_slip` diagnostics).
- **Goal of this doc:** turn the runaway into a *bounded, ordered debugging plan* — one
  instrumentation step that discriminates the mechanism, then the candidate code-level amplifiers
  (each with a concrete guard), then the root-excitation track.

## The mechanism, read off the code

Three code facts combine into the observed signature:

1. **Slip is accumulated on EVERY sub-step node**, magnitude-monotone:
   `d.slip1 += s.V1 * dt_sub; d.slip2 += s.V2 * dt_sub;` (`tpv205_substep_iterator.cpp:122-123`).
2. **`slip_rate` (what `V_max` samples) is written ONLY on the last sub-step**:
   `if (last_sub_step) { flux_.WriteBackState(d, s); }` (`:131-134`), and `WriteBackState` sets
   `data.slip_rate = s.V_abs` (`fault_face_flux.cpp:307`). The driver's `V_max` reads
   `dof_data[i].slip_rate` (`spatial_dyn_driver.cpp:1891`).
3. **`V_abs` has no upper bound** (`tpv205_friction.hpp:139`:
   `V_abs = (tau_abs - tau_strength) / eta_s`), and once `delta = |slip| ≥ d_c` the DOF is
   **permanently** at `mu_d` (`LSWFrictionCoefficient_TPV205:66`, monotone non-increasing — LSW
   never re-locks).

⇒ A single intermediate sub-step node with a transient stress overshoot produces a large `V_abs`,
which is integrated into `slip` (fact 1) but **not** seen by `V_max` (fact 2, only the last node
is sampled) — this is *exactly* `max_slip ≫ V_max·t`. That spike ratchets `delta` past `d_c`
(fact 3) ⇒ the DOF is stuck fully-weakened ⇒ it keeps slipping ⇒ `max_slip` grows monotonically.
The "speckle" is the set of DOFs that took one spurious spike and are now permanently weakened.

The root *excitation* (why a node overshoots) is most likely the ADER sub-step predictor at an
under-resolved cohesive zone (`L_nuc/h_min < 10`), but the **amplifiers that turn a transient
overshoot into permanent garbage are code-level and fixable independent of resolution.** Confirm
the excitation before refining the mesh.

---

## Findings

### [R-001] ACTION (do this first) [tpv205_substep_iterator.cpp:StepOneQP_] — per-sub-step trace at the argmax-slip DOF discriminates the three candidate mechanisms

**Category:** instrumentation (the enabling diagnostic)

**Description:**
The global scalars (`max_slip`, `V_max`) cannot tell you *why* a DOF runs away. One env-gated
trace at the worst DOF, printed **per sub-step** (not per macro-step), settles it in one run:
log `(o, V_abs, tau_abs, tau_strength, sigma_n_total, sigma_n_pos, delta, mu_eff)`. The pattern
tells you which finding below is firing:
- `tau_abs` spikes at an intermediate `o`, `sigma_n_total` stays compressive ⇒ **predictor
  overshoot / under-resolution** (R-006, R-007).
- `sigma_n_total < 0` (tensile) at the spike, `sigma_n_pos → 0`, `tau_strength → 0` ⇒ **tension
  free-slide** (R-005).
- neither spikes but `slip` still grows ⇒ **accumulation/`dt_sub` bug** (look again at R-002 path).

**Suggested fix (add the trace; zero overhead when off):**
```cpp
// In StepOneQP_, after Step 4 (SolveLSW_TPV205), behind an env gate.
#ifdef SEAS_DIAG_SLIP
   static const bool diag = [] {
      const char *e = std::getenv("SEAS_DIAG_SLIP"); return e && e[0] && e[0] != '0';
   }();
   if (diag && d.diag_print) {   // set diag_print=true on the argmax-slip DOF at init/onset
      const real_t sn_pos = std::max<real_t>(s.sigma_n_total, 0.0);
      std::fprintf(stderr,
         "[SLIP] o=%d V_abs=%+.6e tau_abs=%+.6e tau_str=%+.6e "
         "sigma_n_tot=%+.6e sigma_n_pos=%+.6e delta=%+.6e mu_eff=%+.6e\n",
         /*o*/-1, s.V_abs, s.Theta, mu_eff*sn_pos,
         s.sigma_n_total, sn_pos,
         std::sqrt(d.slip1*d.slip1 + d.slip2*d.slip2), mu_eff);
   }
#endif
```
(Pass the sub-step index `o` through `StepOneQP_`'s signature so the print shows *which* node
spikes. Pick the witness DOF by `argmax(sqrt(slip1²+slip2²))` at the first step `max_slip` exceeds,
say, 100 m — the driver already finds an argmax for `[DIAG-ONSET]`,
`spatial_dyn_driver.cpp:1937-1952`; reuse that to set `diag_print`.)

**Test case:** run the Dc2 sbatch with `SEAS_DIAG_SLIP=1` to `t≈1 s`; assert the trace captures at
least one sub-step with `V_abs > 10·V_max_macro`. That single line tells you which of R-002/R-005/
R-006 to fix.

---

### [R-002] CRITICAL [tpv205_friction.hpp:SolveLSW_TPV205] — `V_abs` is unbounded; a transient trial-traction overshoot becomes an arbitrarily large slip increment with no finite/sanity guard

**Category:** BUG (the amplifier)

**Description:**
`V_abs = (tau_abs - tau_strength) / eta_s` (`:139`) has no upper bound and no finiteness check.
`tau_abs` is the *trial* traction at a sub-step node; the ADER predictor can overshoot it
transiently at an under-resolved front. The result is integrated straight into `slip`
(`StepOneQP_:122`) with no opportunity to reject it, and ratchets the DOF past `d_c` (R-003). A
slip rate of, e.g., 1e4 m/s is ~3 orders above the shear-wave speed `c_s` and is physically
impossible; the solver accepts it silently.

This is consistent with CLAUDE.md's "No artificial tau floor" rule (don't perturb the ODE RHS),
but a **sanity tripwire that flags/aborts on an unphysical `V_abs`** is not a tau floor — it
converts a silent garbage-accumulation into a loud, localized failure you can debug.

**Trigger:** a sub-step node where `tau_abs ≫ tau_strength` (predictor overshoot, or
`tau_strength → 0` under tension, R-005).

**Actual behavior:** `V_abs` unbounded → garbage slip increment → permanent weakening (R-003) →
`max_slip` runaway, invisible to `V_max` (R-004).

**Expected behavior:** an unphysical slip rate is detected at the source, not silently integrated.

**Suggested fix (sanity guard, derived from material `c_s` — no magic constant):**
```cpp
   // Closed-form V_abs from radiation damping balance.
   if (tau_abs > tau_strength) { V_abs = (tau_abs - tau_strength) / eta_s; }
   else                        { V_abs = 0.0; }

+  // Sanity tripwire: a slip rate far above the shear-wave speed is
+  // unphysical (radiation damping should preclude it). c_s = 2*eta_s/rho;
+  // pass a characteristic velocity in or derive it from the DOFData impedance.
+  // Flag rather than silently integrate a garbage increment.
+  MFEM_VERIFY(std::isfinite(V_abs),
+              "SolveLSW_TPV205: non-finite V_abs (tau_abs=" << tau_abs
+              << " tau_strength=" << tau_strength << " eta_s=" << eta_s << ")");
+  // (debug build) localise the first unphysical spike; v_cap is c_s-derived,
+  // NOT a hardcoded clamp — it aborts so you can SEE the DOF, not silence it.
+  MFEM_ASSERT(V_abs <= v_sane_max,
+              "SolveLSW_TPV205: V_abs=" << V_abs << " m/s exceeds the sane bound "
+              << v_sane_max << " (~k*c_s); transient stress overshoot at this QP.");
```
where `v_sane_max` is threaded in from the caller as a multiple of `c_s` (e.g. `c_s` itself —
derivable from `eta_s` and `rho`, both already on `DOFData`). **Decision needed:** whether to
*abort* (find the bug) or *clamp* (limp past it) — for debugging, abort; for production resilience,
clamp + count. Recommend abort first to localise, per CLAUDE.md "REPORT the exact error … ASK
before changing the approach."

**Test case:**
```cpp
// Unit test on SolveLSW_TPV205 directly:
// feed tau_abs = 50*tau_strength (a 50x overshoot), assert the guard trips
// (or, in clamp mode, V_abs <= v_sane_max). Confirms the unbounded path is closed.
```

---

### [R-003] CRITICAL [tpv205_substep_iterator.cpp:StepOneQP_ + tpv205_friction.hpp:LSWFrictionCoefficient_TPV205] — slip is a one-way ratchet with no unphysical-growth tripwire; one spurious spike permanently corrupts a DOF and `max_slip` grows with no abort

**Category:** BUG (the persistence mechanism)

**Description:**
`delta = |slip|` only grows (slip accumulates each sub-step, `:122-123`), and
`LSWFrictionCoefficient_TPV205` is monotone non-increasing — once `delta ≥ d_c` the DOF sits at
`mu_d` **forever** (`:66`). That is correct LSW physics for *real* slip, but it means a single
**spurious** increment (R-002) permanently weakens the DOF, which then keeps slipping (low
strength ⇒ `V_abs = (tau_abs − μ_d·σ_n)/η_s > 0` indefinitely) ⇒ `max_slip` ratchets up with **no
detection and no abort**. The run "looks bounded" via `V_max` (R-004) while silently corrupting
state.

The fix is not to make LSW reversible (it shouldn't be); it is (a) stop the spurious increment at
the source (R-002), and (b) add a **per-DOF physical-plausibility tripwire**: total slip cannot
exceed the time-integral of a physical slip rate. `slip_i ≤ ∫|V| dt` is automatic, but
`slip_i ≫ V_max_seen · t` flags exactly the decoupling you observed.

**Trigger:** any DOF that takes one R-002 spike inside the rupture area (`mu_eff < barrier`).

**Actual behavior:** permanent weakening, unbounded `max_slip`, no abort.

**Expected behavior:** a DOF whose accumulated slip is physically impossible for the slip rates
seen aborts (or is flagged with full state) instead of silently growing.

**Suggested fix (driver-side tripwire — cheap, catches the class):**
```cpp
// spatial_dyn_driver.cpp, in the per-step monitor (near :1888-1898):
// track the running max sub-step-aware slip rate (see R-004) and the elapsed time;
// a DOF slip that exceeds, e.g., a few * (V_peak_seen * t) is unphysical.
const real_t slip_i = std::sqrt(d.slip1*d.slip1 + d.slip2*d.slip2);
MFEM_VERIFY(slip_i <= slip_plausibility_factor * V_peak_seen_global * t + slip_floor,
            "spatial_dyn_driver: DOF " << i << " slip=" << slip_i
            << " m is unphysical vs V_peak=" << V_peak_seen_global
            << " m/s over t=" << t << " s — slip ratchet (see R-002/R-003).");
```
(`slip_plausibility_factor` ~ a few; `slip_floor` a small absolute m to avoid early-time
division-by-tiny-`t`. Both are *plausibility* knobs, not physics — and the abort message names the
DOF so you can trace it.)

**Test case:**
```cpp
// Drive a single DOF with one injected V spike that pushes delta past d_c,
// continue stepping; assert the driver tripwire aborts (pre-fix it runs away silently),
// and that with R-002's guard in place the spike never enters so the tripwire never fires.
```

---

### [R-004] MODERATE [tpv205_substep_iterator.cpp + spatial_dyn_driver.cpp] — `V_max` samples only the last sub-step's `slip_rate`, so the blow-up monitor is blind to the intermediate spikes that drive the runaway

**Category:** BUG (monitoring blind spot — the analogue of R-001/R-002 from the reconcile review)

**Description:**
`WriteBackState` runs only on `last_sub_step` (`:131-134`), so `dof_data.slip_rate` — and hence
the driver's `V_max` (`:1891`) and `n_rupturing` (`:1916`) — reflect only the **final** sub-step
node. The slip integral (`:122-123`) sees every node. So a DOF can integrate a huge intermediate
`V` (inflating `slip`) while reporting a small `slip_rate`. This is precisely why the run reported
`V_max ≈ 4.67` while `max_slip → 1e5`: the monitor literally could not see the spike. The
blow-up guard is being fooled the same way `slip_rate` fooled the reconcile diagnostics earlier.

**Trigger:** any intermediate sub-step `V` larger than the final sub-step `V`.

**Suggested fix:** track the max `V_abs` *over sub-steps* and expose it, so the monitor is honest:
```cpp
// In AdvanceWithSubStepStates, accumulate a per-DOF (or global) sub-step V max:
//   v_substep_max[i] = std::max(v_substep_max[i], s.V_abs);  // inside StepOneQP_ via out-param
// Return it to the driver and use it for the blow-up monitor / R-003 tripwire
// (V_peak_seen_global) instead of, or alongside, dof_data.slip_rate.
```

**Test case:**
```cpp
// Two-substep fixture where substep 0 has V=100 and substep 1 has V=1:
// assert dof_data.slip_rate==1 (current) BUT the new v_substep_max==100,
// and that slip accumulated ~ (100+1)*dt_sub — proving the monitor now sees the spike.
```

---

### [R-005] MODERATE [POSSIBLE] [tpv205_friction.hpp:SolveLSW_TPV205] — tension clamp lets the fault slide freely (`τ_strength → 0`) under a transient tensile excursion, a candidate spike source

**Category:** ASSUMPTION (physics modeling)

**Description:**
`sigma_n_pos = max(sigma_n_total, 0)` (`:133`) ⇒ under a transient tensile `sigma_n_total < 0`,
`tau_strength = 0` and `V_abs = tau_abs / eta_s` with **no** frictional resistance (`:134-139`).
The comment calls this "fault opens under tension … freely slides." But the normal channel is the
*unconstrained* one (plan B.2 fact 3: slip does not change `σ_n` inside the solve, and there is no
no-opening limit — this is what reached `σ_n = −124 GPa` in the original blow-up). Combined with
R-002's missing `V` cap, a tensile transient is a direct spike source. Whether the physically
correct response under tension is "free slide" or "open + zero V" (no-opening lock) should be
checked against the SCEC TPV205 spec and SeisSol's `FrictionSolverCommon` — if SeisSol zeros V (or
the slip increment) under tension, this is a deviation that manufactures slip.

**Trigger:** `sigma_n_total < 0` at a sub-step node (confirm via R-001 trace).

**Suggested fix (only if R-001 shows tension at the spike, and only after checking the spec):**
```cpp
   const real_t sigma_n_pos = std::max<real_t>(sigma_n_total, 0.0);
-  const real_t tau_strength = mu_eff * sigma_n_pos;
+  const real_t tau_strength = mu_eff * sigma_n_pos;
+  // Under tension the fault opens; the SCEC/SeisSol convention zeros the
+  // tangential slip response rather than letting it slide frictionlessly.
+  if (sigma_n_total <= 0.0) { V_abs = 0.0; V1 = 0.0; V2 = 0.0;
+     tau1_corr = tau1_trial; tau2_corr = tau2_trial; return; }
```
**Do NOT apply blind** — verify against the spec first; this changes TPV205 physics and must not
break the byte-exact TPV205 regression on runs where `σ_n` never goes tensile (the branch is inert
when `sigma_n_total > 0`, so the regression should be unaffected — confirm).

**Test case:**
```cpp
// SolveLSW_TPV205 with sigma_n_total = -1e6 (tension), tau_abs = 1e7:
// current -> V_abs = 1e7/eta_s (large); fixed -> V_abs == 0 (no-opening).
// Plus: assert with sigma_n_total > 0 the result is byte-identical to before.
```

---

### [R-006] MODERATE [POSSIBLE] [tpv205_substep_iterator.cpp:AdvanceWithSubStepStates] — within a macro step the bulk does not relax as the fault slips, so trial traction can stay high across sub-steps and over-accumulate slip

**Category:** ASSUMPTION (operator-split sub-step scheme)

**Description:**
All `Q_pointwise_*[o]` are predicted from the **start-of-macro-step** bulk via
`ComputeADERSubStepStates` (`spatial_dyn_driver.cpp:416-429`); they do **not** see the slip
accumulated by earlier sub-steps. The friction solve weakens `mu_eff` across sub-steps (delta
grows) but the *trial traction* it is handed does not relax in response to the slip. So at a DOF
that should be shedding stress as it slips, the trial `tau_abs` can remain high at later nodes ⇒
`V_abs` stays high ⇒ slip over-accumulates within the macro step (up to ~O sub-steps' worth). For
normal rupture rates this is negligible; in the spike regime it compounds the runaway. This is a
property of the predictor/sub-step split, not a one-line bug — flag it as a candidate contributor
to confirm via R-001 (does `tau_abs` stay high across nodes while `delta` grows?).

**Trigger:** large `V_abs` sustained across multiple sub-step nodes within one macro step.

**Suggested fix:** none mechanical — if R-001 confirms this dominates, the remedy is a smaller
macro `dt` (so the bulk relaxes between steps) or a corrected sub-step scheme that feeds the
imposed state back into the predictor. Investigation, not a patch.

**Test case:** n/a (design-level; the R-001 trace is the test — look for sustained high `tau_abs`
across `o` with growing `delta`).

---

### [R-007] INVESTIGATION [ComputeADERSubStepStates / mesh] — the root excitation is most likely an under-resolved cohesive zone; confirm with a resolution + physical-`D_c` study before refining blindly

**Category:** root-cause track (the "speckle" proper)

**Description:**
B.2.5 predicted this exact outcome (`worst_rel → 0` but speckle survives ⇒ under-resolution
`L_nuc/h_min ≈ 3–10 < 10` or SSO). If R-001 shows `tau_abs` overshoot with compressive `σ_n`
(no tension), the excitation is the predictor at a front the mesh can't resolve. The test is a
**mesh-refinement convergence study**: if `max_slip`/the speckle *converges* as `h → 0`, it's
numerical under-resolution (refine the nucleation patch so `L_nuc/h ≥ ~5–10`); if it does **not**
converge, the nucleation is genuinely ill-posed and you report ensemble behavior, not a single
trajectory. Also run the physical `D_c = 1.0` config (the plan's Phase-4 step 3) — the inflated
`Dc2` widens `L_nuc` but the absolute `h` is unchanged, so `Dc2` is *not* a substitute for
refining `h`.

**Suggested action:** (1) confirm mechanism via R-001; (2) if overshoot-without-tension, run
`{h, h/2, h/4}` on the nucleation patch and plot `max_slip(t)` + the speckle count; (3) decide
refine-vs-ill-posed from convergence. Independent of this, R-002/R-003/R-004 should land first so
the runaway is *caught and localized* rather than silently growing during the study.

---

## Recommended order of operations

1. **R-001** (per-sub-step trace at the argmax DOF) — one run, tells you which branch you're in.
2. **R-004** (honest sub-step `V_max`) + **R-003** (slip-plausibility tripwire) — so the runaway
   is *visible and localized*, not silent. Cheap, no physics change.
3. **R-002** (V sanity guard) — close the unbounded path; abort-to-localize first.
4. Branch on R-001's verdict: tension ⇒ **R-005** (check spec, then no-opening); sustained
   overshoot ⇒ **R-006/R-007** (dt / mesh-refinement study).

## Summary
- Critical issues: 2 (R-002 unbounded `V_abs`; R-003 ratchet with no tripwire)
- Moderate issues: 3 (R-004 monitor blind spot; R-005 tension free-slide [POSSIBLE];
  R-006 within-step bulk non-relaxation [POSSIBLE])
- Investigation/instrumentation: 2 (R-001 trace; R-007 resolution study)
- Verdict: **the runaway is real, not cosmetic, and is currently silent.** Land R-004/R-003/R-002
  (catch + localize + bound) before the R-007 resolution study, then fix the confirmed root.
  Do **not** ship "it's under-resolution, refine the mesh" until R-001 rules out the tension
  (R-005) and accumulation (R-002) amplifiers — those are code bugs a finer mesh would only
  postpone.

## Unreviewed Areas
- `ComputeADERSubStepStates` / `ComputeADERTimeIntegrated` (the CK predictor) — not opened here;
  R-006/R-007 hinge on whether it overshoots at the front. Trace via R-001 before diving in.
- The gradual_overstress `nuc_callback` writing `tau*_nuc` per sub-step — not re-audited for an
  unbounded ramp; if R-001 shows `tau_abs` driven by `tau*_nuc` rather than the trial, audit the
  nucleation accumulator next.
- Whether `max_slip` is one DOF or many — the driver's `maxslip_g` is a global max
  (`spatial_dyn_driver.cpp:1926`); add a count of DOFs with `slip > 100 m` to know if the speckle
  is a handful of ratcheted DOFs or a spreading region.
