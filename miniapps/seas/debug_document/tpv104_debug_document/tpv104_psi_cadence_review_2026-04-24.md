# Code Review: TPV104 ψ-Cadence Analysis (2026-04-24)

## Review Scope
- **Subject**: User's side-by-side analysis of MFEM vs SeisSol ψ-update cadence and the proposed next-step options (wire `Tpv104SubStepIterator` after `wave.AdvanceADER`, vs run a `--dt 5e-5` halving experiment).
- **Files referenced**:
  - `miniapps/seas/friction/slip_law_srw_psi.hpp:92-124` (`UpdateStateAnalyticSlipLawSRW`).
  - `miniapps/seas/drivers/tpv104_driver.cpp:993-1054` (production time-loop body).
  - `SeisSol/src/DynamicRupture/FrictionLaws/CpuImpl/FastVelocityWeakeningLaw.h:43-78`.
  - `SeisSol/src/DynamicRupture/FrictionLaws/CpuImpl/RateAndState.h:30-78` and `:151-195`.
  - `SeisSol/src/DynamicRupture/FrictionLaws/RateAndStateCommon.h:32-43`.
  - `SeisSol/src/DynamicRupture/Misc.h` (`TimeSteps = ConvergenceOrder`).
- **Domain context**: `miniapps/seas/CLAUDE.md`, `tpv104_debug_plan_2026-04-24.md`, `research_seissol_tpv104_2026-04-24.md` (revised).

## Bottom-line up front

**Structurally agree on the cadence claim** — the formula is byte-identical, neither side has a ψ floor, and the driver loop is what differs (1 call/macro-step in MFEM vs O × N_kaneko + 1 calls/macro-step in SeisSol). The R7-001(b) disclosure remains accurate.

**Disagree on the magnitude argument** — the analysis claims that at the rupture front the macro-step path "jumps directly to ψ_ss in one step" with `exp1m ≈ 0.91`, but the cited numbers (`V = 9.45 m/s, dt_macro = 0.28 ms, L = 0.4 m`) give `V·dt/L = 6.73e-3` and `exp1m = 0.67%` — not 91%. The "one-jump" regime would require `dt_macro ≈ 100 ms`, which is two orders of magnitude above the CFL-stable dt for TPV104. Verified by `python3 -c "..."` arithmetic below. **This means the cadence gap is much smaller per macro step than the analysis claims, and the Phase-3 probe-2 drift over the rupture phase is bounded by O(dt²) accumulated over the pulse, not by a single-step jump.**

**Recommendation on next steps**: **dt-halving first; iterator-wiring second.** The analysis's qualitative argument for cadence-being-the-bug rests on the "one-jump" magnitude claim that does not hold under realistic dt. Halving dt is a cheap test that distinguishes "cadence is the dominant error" (peak collapses) from "something else is wrong" (peak persists or grows). Iterator wiring around the [C2] no-touch list is a much larger, riskier change; it should be undertaken only if dt-halving confirms cadence is the root cause.

---

## Findings

### [R-001] [MODERATE] [BUG] arithmetic-error-in-magnitude-claim — `V·dt/L ≈ 2.4 → exp1m = 0.91` is wrong by ≈400×

**Category:** ASSUMPTION / BUG (in the analysis, not the code).

**Description:**
The analysis writes:
```
For our Frontera run at the hypocenter:
- macro V·dt/L = 9.45 · 1e-4 / 0.4 ≈ 2.4 → exp1m = 0.91 → essentially full jump to ψ_ss
- per-sub-step V·dt_sub/L would be 0.012 → exp1m = 0.012 → 1.2% step, gradual
```

The first computation is internally inconsistent. `9.45 × 1e-4 / 0.4 = 2.36e-3`, not `2.4`. The arithmetic appears to substitute `1e-1` for `1e-4` somewhere. With the user's stated `V = 9.45 m/s, dt_macro = 0.28 ms = 2.85e-4 s, L = 0.4 m`, the actual values are:

```
V·dt_macro/L = 9.45 × 2.85e-4 / 0.4 = 6.73e-3
exp1m_macro  = 1 − exp(−6.73e-3)    = 6.71e-3 ≈ 0.67%

V·dt_sub/L (O=5) = 9.45 × 5.7e-5 / 0.4 = 1.35e-3
exp1m_sub        = 1 − exp(−1.35e-3) = 1.35e-3 ≈ 0.135%
```

Verified with `python3` (logged with this review).

**Trigger:**
Any reader who relies on the analysis's claim that `exp1m_macro ≈ 0.91` to justify "ψ jumps directly to ψ_ss in one macro step" will reach a wrong conclusion about the structural difference between the two codes.

**Actual behavior of analysis:**
Reports `exp1m_macro = 0.91` (full jump per macro step). Concludes that MFEM's ψ commits to ψ_ss in one step at the rupture front.

**Expected behavior:**
Reports `exp1m_macro ≈ 0.67%` per macro step at `dt = 0.28 ms`. The "full jump" regime requires `V·dt/L ≳ 2`, i.e. `dt ≈ 100 ms` at `V = 10 m/s, L = 0.4 m` — **two orders of magnitude above the CFL-stable dt for TPV104**. Under any stable run, neither MFEM nor SeisSol jumps to ψ_ss in one step.

**Suggested fix (in the analysis, not in code):**
Replace the table row with the correct numbers and reframe the conclusion:

```diff
- For our Frontera run at the hypocenter:
- - macro V·dt/L = 9.45 · 1e-4 / 0.4 ≈ 2.4 → exp1m = 0.91 → essentially full jump to ψ_ss
- - per-sub-step V·dt_sub/L would be 0.012 → exp1m = 0.012 → 1.2% step, gradual
+ For the Frontera run at the hypocenter (V = 9.45 m/s, dt = 0.28 ms, L = 0.4 m):
+ - macro V·dt/L     = 6.73e-3 → exp1m_macro = 0.67% per macro step
+ - per-sub-step at O=5: V·dt_sub/L = 1.35e-3 → exp1m_sub  = 0.135% per sub-step
+ Both paths advance ψ smoothly per step; the difference is which V they
+ use (macro-averaged vs per-sub-step).  Over the full rupture pulse
+ (~500 macro steps in the rapid-V phase) the cumulative deviation in
+ ψ is bounded by O(dt²), not by a one-step jump.
```

**Test case** (Python, no MFEM build needed):
```python
import math
def exp1m(V, dt, L): return 1 - math.exp(-V*dt/L)
assert abs(exp1m(9.45, 2.85e-4, 0.4) - 0.00671) < 1e-5
assert abs(exp1m(9.45, 5.7e-5,  0.4) - 0.00135) < 1e-5
# one-jump regime requires dt ≈ 100 ms:
assert abs(exp1m(9.45, 0.10,    0.4) - 0.906)   < 1e-2
```

---

### [R-002] [MODERATE] [BUG] Kaneko inner-iteration count claimed `5`, actual is `2` (plus `+1` final call)

**Category:** ASSUMPTION (in the analysis).

**Description:**
The table row reads:
```
Calls per macro step: 1   |   O × 5 (O sub-steps × 5 Kaneko inner iterations)
```

The actual SeisSol settings are at `RateAndStateCommon.h:32-43`:
```cpp
struct Settings {
  const uint32_t maxNumberSlipRateUpdates{60};   // Newton iterations on V (NOT updateStateVariable calls)
  const uint32_t numberStateVariableUpdates{2};   // Kaneko outer iterations
  const double newtonTolerance{1e-8};
};
```

So per ADER sub-step:
1. `updateStateVariableIterative` (`RateAndState.h:151-195`) does **2** Kaneko outer iterations, each calling `updateStateVariable` once and then `invertSlipRateIterative` (Newton on V — does NOT call `updateStateVariable`).
2. `calcSlipRateAndTraction` (`RateAndState.h:197-256`) calls `updateStateVariable` **once more** at line 211 with the final mean V.

**Total `updateStateVariable` calls per sub-step: 3.** Per macro step at O=5 (`TimeSteps = ConvergenceOrder`): `3 × 5 = 15`, **not** `O × 5 = 25`.

**Trigger:**
The "5 Kaneko" misstatement does not affect the qualitative cadence-gap argument, but it overstates the per-macro-step work imbalance and could mislead the agent that wires the iterator (it would size buffers for 25 calls, not 15).

**Suggested fix:**
```diff
- Calls per macro step: 1   |   O × 5 (O sub-steps × 5 Kaneko inner iterations)
+ Calls per macro step: 1   |   O × (numberStateVariableUpdates + 1) = O × 3 = 15 at O=5
+                            (2 Kaneko outer iters per sub-step + 1 final call in
+                             calcSlipRateAndTraction; numberStateVariableUpdates is
+                             const at 2 in RateAndStateCommon.h:41)
```

**Test case:**
Static-grep against the SeisSol source:
```bash
grep -A2 'numberStateVariableUpdates' \
  SeisSol/src/DynamicRupture/FrictionLaws/RateAndStateCommon.h | grep '{2}'
```

---

### [R-003] [LOW] [QUALITY] stale line citation `tpv104_driver.cpp:1024-1037`

**Category:** QUALITY.

**Description:**
The analysis cites `drivers/tpv104_driver.cpp:1024-1037` for the call site of `UpdateStateAnalyticSlipLawSRW`. Lines 1024-1037 are the **R7-007 disclosure comment block**, not the call. The actual call is at `tpv104_driver.cpp:1040-1051`:

```cpp
1040:         dof_data[i].psi = UpdateStateAnalyticSlipLawSRW(
1041:            psi_n[i],
1042:            dof_data[i].slip_rate,
1043:            dof_data[i].Dc,
1044:            dt_step,
1045:            V_w[i],
1046:            dof_data[i].a,
1047:            TPV104Params::b,
1048:            TPV104Params::V0,
1049:            TPV104Params::f0,
1050:            TPV104Params::f_w);
```

with `psi_n[i] = dof_data[i].psi` saved at line 1001 before the macro step body. The `slip1/slip2` accumulation is at lines 1052-1053.

**Suggested fix:**
```diff
- Called from drivers/tpv104_driver.cpp:1024-1037 once per macro step
+ Called from drivers/tpv104_driver.cpp:1040-1053 once per macro step
+ (psi_n captured at :1001 before wave.AdvanceADER at :1019)
```

---

### [R-004] [LOW] [CONFIRMED] no ψ floor on either side

**Category:** ASSUMPTION (correctly identified).

**Description:**
The analysis correctly notes that **neither** `UpdateStateAnalyticSlipLawSRW` nor SeisSol's `updateStateVariable` clamps ψ to a non-negative range. ψ is allowed to go negative — and physically does, because for TPV104 inside the VW core at `V ≈ 10 m/s`:

```
f_LV(10)  = max(0, 0.6 − (0.014−0.01)·ln(10/1e-6))
          = max(0, 0.6 − 0.0644) = 0.5355
f_ss(10)  = 0.1 + (0.5355 − 0.1) / (1 + (10/0.1)^8)^(1/8)
          ≈ 0.1 + 0.4355 / 100 = 0.1044
ψ_ss(10)  = 0.01·log((2e-6/10) · sinh(0.1044/0.01))
          ≈ −0.0568
```

(Verified with `python3`.)

So ψ_ss < 0 at high V is the **physically correct** asymptote. Adding a `psi = max(0, psi)` floor would suppress the FVW weakening signature at the rupture front, which is exactly what TPV104 is supposed to test. The analysis's recommendation to **not** add the floor is correct.

**No action.**

---

### [R-005] [MODERATE] [POSSIBLE BUG IN PROPOSAL] iterator wiring as the primary fix path is premature

**Category:** ASSUMPTION (in the analysis's recommendation).

**Description:**
The analysis offers two follow-ups: (a) wire the iterator after `wave.AdvanceADER` using `I_±` reconstructed from `Q_n` and `Q_{n+1}`, or (b) run dt-halving (`--dt 5e-5`) as a cheap experiment. Implicit in the framing is that (a) is the "principled fix" and (b) is a sanity check.

This ordering is **inverted from what the data justify**.

* The cadence-gap argument's central numerical claim (`exp1m_macro = 0.91`, "one-step jump to ψ_ss") is unsupported under realistic dt (R-001). The actual per-step jump is ~0.7%, which means the cumulative deviation between macro-step and per-sub-step ψ trajectories is O(dt²) per step — small.
* If the cadence gap is the dominant error, halving dt (or equivalently running with `--dt 5e-5`) should ~quadruple-collapse the peak-V error (O(dt²) scaling). This is a 1-line CLI change, no source edit.
* If halving dt does **not** collapse the peak-V error, the cadence gap is **not** the dominant error and iterator-wiring won't fix it either — saving days of work on a [C2] workaround.
* Reconstructing per-sub-step `I_±` from only the macro-step endpoints `(Q_n, Q_{n+1})` requires assumptions about the predictor trajectory that ADER does not expose. A literal "linear interpolation" reconstruction is **not** what SeisSol's ADER predictor produces; the ADER predictor inside `wave.AdvanceADER` evaluates a Taylor-expanded space-time polynomial at the sub-step quadrature nodes, not a linear interpolant. So option (a) as framed would deliver `I_±` values that disagree with SeisSol's per-sub-step `qInterpolated[o]` by O(dt²) — the same error magnitude that the iterator was supposed to remove. **The proposed (a) cannot recover the per-sub-step cadence without exposing the ADER predictor's internal sub-step states, which requires editing `wave_operator.inl`.**

**Suggested ordering**:

```diff
- Want me to dig further into how to wire the iterator without touching
- the [C2] no-touch list ... Or do you want to first investigate dt-halving?
+ Run dt-halving first.  It is a 1-line CLI change with no source edit:
+   --dt 5e-5   (or smaller, until O(dt²) scaling is verified)
+ Compare the resulting peak-V at the hypocenter to the dt = 2.85e-4
+ baseline.
+   - If peak-V error drops by ≥ 4× (consistent with O(dt²)), the cadence
+     gap is the dominant error, and iterator wiring is the right next
+     step.
+   - If peak-V error is unchanged or grows, the cadence gap is NOT
+     dominant; iterator wiring will not fix it.  Look elsewhere
+     (nucleation cadence I-04, friction-solver dispatch I-15, or
+     tangent-frame mapping I-11).
+ Only after dt-halving confirms cadence is the bug, propose iterator
+ wiring — and DO NOT attempt to reconstruct per-sub-step I± from
+ (Q_n, Q_{n+1}); that approach trades the macro-step error for an
+ ADER-predictor-mismatch error of the same O(dt²) magnitude.
```

**Test case:**
The dt-halving experiment itself is the test. Acceptance criterion: at `--dt 1e-5` vs `--dt 2.85e-4`, the peak-V relative error against the SeisSol reference at the hypocenter should drop by `≥ (2.85e-4 / 1e-5)² ≈ 800×` if cadence is the dominant error. If it drops by less than 4×, look elsewhere.

---

### [R-006] [LOW] [CONFIRMED] formula byte-equivalence — no action

**Category:** ASSUMPTION (correctly identified).

**Description:**
The analysis claims:
```
The analytic formula is byte-identical to MFEM's UpdateStateAnalyticSlipLawSRW.
Same psi_ss * exp1m + exp1v * stateVarReference. No floor here either.
```

Confirmed by line-for-line read:

| MFEM `slip_law_srw_psi.hpp:120-123` | SeisSol `FastVelocityWeakeningLaw.h:71-74` |
|---|---|
| `preexp1 = -V * (dt / L);` | `preexp1 = -localSlipRate * (timeIncrement / localSl0);` |
| `exp1v = std::exp(preexp1);` | `exp1v = std::exp(preexp1);` |
| `exp1m = -std::expm1(preexp1);` | `exp1m = -std::expm1(preexp1);` |
| `return psi_ss * exp1m + exp1v * psi_old;` | `return steadyStateStateVariable * exp1m + exp1v * stateVarReference;` |

**No action.** This is the part of the analysis that is rock-solid.

---

## Summary

- **Critical issues**: 0
- **Moderate issues**: 3 (R-001 magnitude error, R-002 Kaneko count, R-005 ordering of next steps)
- **Low issues**: 3 (R-003 stale line numbers, R-004 confirmation, R-006 confirmation)
- **Plan compliance**: PARTIAL — the analysis correctly identifies the cadence gap (matches the R7-007 disclosure already in the driver source and the I-17 entry in the revised research doc) but quantifies it incorrectly.
- **Verdict**: **PASS WITH FIXES** for the structural finding (cadence gap is real, formula is identical, no floor either side, no `max(0, ψ)` hack). **FAIL** on the magnitude argument and the proposed ordering of next steps (R-001, R-002, R-005). Take the dt-halving path before attempting any iterator wiring.

## Direct answer to "do you agree on this analysis?"

**Partially.**

| Claim | Verdict |
|-------|---------|
| Formula identical between codes | ✅ Confirmed |
| No ψ floor on either side | ✅ Confirmed |
| MFEM updates ψ once per macro step with time-averaged V | ✅ Confirmed |
| SeisSol updates ψ multiple times per macro step (per sub-step + Kaneko averaging) | ✅ Confirmed (3 calls × O sub-steps; analysis's "5 Kaneko" is wrong) |
| `max(0, ψ)` floor would be a hack, not a fix | ✅ Confirmed |
| `exp1m_macro ≈ 0.91`, "ψ jumps directly to ψ_ss in one macro step" | ❌ Wrong (actual ≈ 0.67% per step at the cited dt; 91% requires dt ≈ 100 ms) |
| Cadence gap is the dominant cause of the peak-V mismatch | ⚠️ Unverified (analysis's magnitude argument is wrong; dt-halving will verify or refute) |
| Iterator wiring is the principled fix | ⚠️ Possibly, but premature; option (a) as framed (reconstruct I± from `Q_n, Q_{n+1}`) cannot match SeisSol's per-sub-step predictor without editing `wave_operator.inl` |
| dt-halving is the right next experiment | ✅ Strongly endorsed; do this first |

## Unreviewed Areas
- The actual peak-V plot from the Frontera run that motivated this investigation — not provided in the analysis text. The recommendation to dt-halve assumes there IS a measurable peak-V mismatch in the production run; if not, the cadence question is moot.
- The full ADER-predictor internals in `wave_operator.inl` — flagged in R-005 as the technical reason option (a) is harder than the analysis suggests, but not audited line-by-line here. A separate review is needed before any iterator-wiring patch goes in.
