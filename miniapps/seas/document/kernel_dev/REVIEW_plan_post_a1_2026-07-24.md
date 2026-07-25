# Code Review: PLAN_performance_program_2026-07-23.md — post-A1 consistency audit (2026-07-24)

## Review Scope
- Plan: `document/kernel_dev/PLAN_performance_program_2026-07-23.md` (215 lines, as of the A1-results commit)
- Files reviewed: the plan; `document/comm_dev/RESULTS_a1_merge_2026-07-24.md`; `document/comm_dev/RESULTS_a0_wiggle_2026-07-24.md`; `document/comm_dev/DESIGN_a1_per_tick_merge_2026-07-24.md`; `document/comm_dev/PLAN_lts_comm_reduction_2026-07-22.md`
- Domain context: measured artifacts from jobs 52422891 (A0) and 52472765 (A1); the plan's own governing rule (measured body, quarantined projections)
- Method: single-reviewer adversarial pass. Every quantitative claim recomputed by script; every causal claim checked against the A1 experiment, which is the newest and most discriminating evidence.

The A1 result (bitwise-correct merge, **zero** speedup) falsified the mechanism behind three of Track A's four remaining phases. The A1 banner and the work-list strikes were added, but the plan was **not swept for downstream statements built on the dead mechanism**. That sweep is this review.

## Findings

### [R-001] [CRITICAL] [PLAN:154 — Known-unknowns #1] — A falsified causal claim is still marked "RESOLVED"

**Category:** BUG (falsified claim carried as settled)

**Description:**
Line 154: *"~~Whether exposed wait scales with exchange-round rate.~~ **RESOLVED by A0: yes, super-proportionally (1.807× per 1.587×).**"* This is now known to be **false as stated**. A0 varied sync count with rounds/sync fixed (271→171 syncs, 125 rounds/sync in both legs), so round rate and sync rate moved together — that experiment could not distinguish them. A1 (job 52472765) is the deconfounding experiment: it varied round rate alone (52 syncs in both legs, 125→64 rounds/sync) and wait did **not** move (0.971×). The resolved answer is: wait scales with **sync rate**, NOT round rate.

**Trigger:** Anyone reading the known-unknowns list to learn what is settled — the section's entire purpose.

**Actual behavior:** The plan asserts, as a resolved fact, the hypothesis its own newest experiment refuted.

**Expected behavior:** The entry states the corrected resolution and cites both experiments as the pair that separated the variables.

**Suggested fix:**
```diff
-1. ~~Whether exposed wait scales with exchange-round rate.~~ **RESOLVED by A0: yes, super-proportionally (1.807× per 1.587×).**
+1. ~~Whether exposed wait scales with exchange-round rate.~~ **RESOLVED by A0+A1 jointly: it does NOT —
+   wait scales with SYNC rate.** A0 (rounds/sync fixed, syncs ×0.63) → wait ×0.55; A1 (syncs fixed,
+   rounds/sync ×0.51) → wait ×1.03. A0 alone confounded the two variables; A1 deconfounded them
+   (job 52472765).
```

**Test case:**
```bash
# Fails while the falsified resolution stands:
! grep -q "RESOLVED by A0: yes, super-proportionally" \
    document/kernel_dev/PLAN_performance_program_2026-07-23.md
```

---

### [R-002] [CRITICAL] [PLAN:102-111] — Two adjacent sections assert opposite conclusions

**Category:** BUG (internal contradiction)

**Description:**
Lines 91–100 (the A1 banner) correctly state the merge bought nothing and A2/A3 are devalued. Ten lines later, the un-updated paragraph at 102–107 says *"exposed wait does scale with the exchange-round rate… That converts A1–A3 from a fitted model into measurement-backed work… A1–A3 pay off by removing sync points rather than bytes."* Every clause is falsified — and A1–A3 did **not** remove sync points, which is precisely why they failed. Lines 109–111 then give operational guidance for A3 (*"A3 carries a liveness gate…"*) as if A3 were live, though the table above strikes it.

**Trigger:** Reading the Track-A section top to bottom.

**Actual behavior:** The plan simultaneously says A1–A3 are measurement-backed and that A1 measured zero.

**Expected behavior:** The stale paragraph becomes a short historical note; the A3 guidance is conditioned on A3 ever being revived.

**Suggested fix:**
```diff
-**A0 is DONE (2026-07-24) and it PASSED** — see the A0 RESULT section above. It confirmed the
-assumption the entire track rested on: exposed wait does scale with the exchange-round rate, and
-super-proportionally (1.807× per 1.587×). That converts A1–A3 from a fitted model into
-measurement-backed work. It also delivered a 1.370× end-to-end speedup for one config line, and
-redirected the *mechanism*: the wait is synchronisation exposure, not bandwidth (wire time is
-0.15 % of it), so A1–A3 pay off by removing sync points rather than bytes.
-
-**A3 carries a liveness gate, not just a correctness gate.** Bitwise identity cannot detect this
-code's actual recorded failure mode for exchange changes — the R-1600 unmatched-collective hang
-(`wave_operator.inl:3661-3672`). Deadlock produces no wrong bits; it produces no bits.
+*(Historical note: after A0, this section read A0's result as "wait scales with round rate" and
+declared A1–A3 measurement-backed. A1 falsified that reading — see the banner above. If A3 is ever
+revived, it needs a liveness gate in addition to bitwise identity: the recorded R-1600
+unmatched-collective hang, `wave_operator.inl:3661-3672`, produces no wrong bits — it produces no
+bits.)*
```

**Test case:**
```bash
! grep -q "converts A1–A3 from a fitted model into measurement-backed work" \
    document/kernel_dev/PLAN_performance_program_2026-07-23.md
```

---

### [R-003] [CRITICAL] [PLAN:168-214 — Appendix inventory] — The endpoint is overstated ~1.9× by savings A1 just zeroed

**Category:** BUG (falsified numbers presented as current projections)

**Description:**
The inventory still credits A1 = 611, A2 = 320, A3 = 72 s/sim-s and concludes *"Endpoint if everything lands: ~1107 s/sim-s ≈ 4.6× SeisSol. Read it as '~5×, maybe'."* Those rows are measured-dead (A1) or devalued by the same mechanism (A2/A3). Recomputed: **~2110 s/sim-s ≈ 8.8×** on kernel items alone, or a band **[6.2, 8.8]×** folding in A6's unsized [0, ~640]. Also falsified with them: structural-fact #1 (*"A1 + A2 + B1/B2 carry ~85 %"*; *"A2 is the second-largest comm item, not a rounding error"*) and the trailing bullets (*"Track A, if wait scales with rounds: 2341 → ~300"*; *"Composed: ~3.8–5.6×"*). The appendix is declared non-load-bearing, but the endpoint is the number that gets quoted upward — a ~2× error there produces wrong decisions.

**Trigger:** Anyone quoting the program's expected endpoint.

**Actual behavior:** Endpoint 4.6× / "~5×, maybe" — unreachable on current evidence.

**Expected behavior:** A1/A2/A3 rows struck with the measured verdict; A6 added with its bounds; endpoint restated ~8.8× (kernels only) with the conditional A6 band; structural-fact #1 and both trailing bullets corrected.

**Suggested fix (pattern; apply the work-list strike-through style):**
```diff
-| **A1** | merge 125→**64** rounds/sync | wait (1251) | 611 | 2512 | 10.5× | projected from A0's law |
-| **A2** | single-round tick, 64→**32** | residual wait | 320 | 2192 | 9.2× | same law; **NOT a minor item** |
-| A3 | split-post overlap | residual wait | 72 | 2120 | 8.9× | low; anti-synergistic with B |
+| ~~A1~~ | merge 125→64 | wait | ~~611~~ **0 — MEASURED** (52472765: wait ×0.971) | 3123 | 13.1× | dead |
+| ~~A2~~ | 64→32 | wait | ~~320~~ ~0 | 3123 | 13.1× | devalued — the mechanism A1 falsified |
+| ~~A3~~ | overlap | wait | ~~72~~ ~0 | 3123 | 13.1× | devalued — wire is 0.15 % of the wait |
+| **A6** | rank rebalancing | wait skew | **0–640, UNSIZED** | 2483–3123 | 10.4–13.1× | bounds only; size first |
...
-**Endpoint if everything lands: ~1107 s/sim-s ≈ 4.6× SeisSol.** Read it as "~5×, maybe".
+**Endpoint: ~2110 s/sim-s ≈ 8.8× SeisSol on kernel items alone; [6.2, 8.8]× if A6 delivers its
+unsized upper bound.** The pre-A1 "~5×, maybe" is withdrawn.
```

**Test case:**
```python
def test_R003_endpoint_recompute():
    claimed, dead = 1107.0, 611 + 320 + 72
    assert abs((claimed + dead) / 239 - 8.8) < 0.1   # corrected endpoint ≈ 8.8×
```

---

### [R-004] [MODERATE] [PLAN:31-56 — A0 RESULT section] — Superseded interpretation carries no forward pointer

**Category:** DEVIATION (stale interpretation, unmarked)

**Description:**
The A0 section's bullet 1 (*"Track A GO… wait scales super-proportionally with the exchange-round rate"*) and bullet 3 (*"A1–A3 remain right, but the mechanism is fewer exposure events"*) record the pre-A1 reading. A0's **measurements** all stand; its **interpretation** was revised by A1, and nothing in the section says so. It precedes the work list, so a reader absorbs the dead model first.

**Trigger:** Reading the plan in order.

**Actual behavior:** The superseded causal reading carries the same authority as the still-valid measurements.

**Expected behavior:** One annotation at the section head.

**Suggested fix:**
```diff
 ## A0 RESULT (2026-07-24, jobs 52416603 + 52422891) — first measured milestone
+
+> **Interpretation revised by A1 (52472765):** every measurement below stands, but the causal
+> reading — "wait scales with round rate" — was a confound; A1 separated the variables and the
+> driver is SYNC rate. Bullet 1's "Track A GO" and bullet 3's "A1–A3 remain right" are superseded.
```

**Test case:**
```bash
grep -q "Interpretation revised by A1" document/kernel_dev/PLAN_performance_program_2026-07-23.md
```

---

### [R-005] [MODERATE] [PLAN:15-29 — Direction section] — "Communication first" no longer matches the plan's own live work list

**Category:** DEVIATION (direction statement vs. work-list inconsistency)

**Description:**
The direction says *"1. Communication first."*, justified in part by *"It is the only bitwise-safe track. Merging exchange rounds is identity-preserving by construction."* After A1: the identity-preserving merge is proven worthless; Track A's live items are A5 ("no new gain at np=256") and A6 (unsized); Track B holds the best-supported next action (B1: measured target, benched mechanism). The plan's own tables now lead with B. Justifications 1 and 3 survive — so this is a **scoped restatement within the settled direction**, not a direction change, and should say so explicitly to prevent re-litigation.

**Trigger:** Deciding the next action from the direction statement alone.

**Actual behavior:** Direction says comm-first; the evidence-ranked next action is B1 plus an A6 *sizing measurement*.

**Expected behavior:** Direction restated: comm remains the largest line item, but its only live lever is unsized; next actions are **B1 now, A6 sizing in parallel**; comm work resumes when A6 is sized.

**Suggested fix:**
```diff
-1. **Communication first.**
-2. **Kernels second**, predictor before faces, cheapest lever first.
+1. **B1 now** (the best-supported item program-wide: measured target, benched mechanism), with the
+   **A6 sizing measurement in parallel** — comm remains the largest line item (~40 % of wall), but
+   after A1 its only live lever is unsized, and committing to it unsized would repeat A1's mistake.
+2. **Remaining kernels next**, predictor before faces, cheapest lever first.
@@
-- **It is the only bitwise-safe track.** Merging exchange rounds is identity-preserving by
-  construction; kernel work needs tolerance gates.
+- ~~It is the only bitwise-safe track~~ — true, and A1 proved bitwise-safe work can still buy
+  nothing. Safety alone no longer orders the tracks.
```

**Test case:**
```bash
grep -q "B1 now" document/kernel_dev/PLAN_performance_program_2026-07-23.md
```

---

### [R-006] [MODERATE] [PLAN:91; RESULTS_a1; RESULTS_a0] — MPI_Waitall CALLS mislabeled as "rounds" (2× everywhere)

**Category:** BUG (mislabeled quantity, factor 2)

**Description:**
The A1 banner says *"Rounds 12,812 → 6,556"*; RESULTS_a1's table row reads *"exchange rounds/rank"*; RESULTS_a0 labels 67,746/42,686 the same way. All are **Waitall calls = 2× rounds** (`ExchangeFaceNbrData` posts two Waitalls). Check: 12,812/2/52 = 123.2 rounds/sync ≈ the designed 125; 6,556/2/52 = 63.0 ≈ 64. Every ratio is unaffected (the 2 cancels) so no conclusion changes — but a reader reconciling these counts against `BuildTickTable::n_collectives` will hit a clean 2× discrepancy, and this project already produced one wrong verdict from a misread Waitall row.

**Trigger:** Reconciling documented counts against the tick table or Caliper.

**Actual behavior:** Absolute counts labeled "rounds" are calls.

**Expected behavior:** Label "Waitall calls/rank (= 2× exchange rounds)" in the plan banner and both results docs, noting ratios are unaffected.

**Suggested fix:**
```diff
-> **A1 RESULT (job 52472765) — the merge is correct and buys NOTHING.** Rounds 12,812 → 6,556
+> **A1 RESULT (job 52472765) — the merge is correct and buys NOTHING.** Waitall calls/rank
+> (= 2× exchange rounds) 12,812 → 6,556
```

**Test case:**
```python
def test_R006_calls_are_2x_rounds():
    assert abs(12812 / 2 / 52 - 125) < 2 and abs(6556 / 2 / 52 - 64) < 1.5
```

---

### [R-007] [MODERATE] [PLAN:89 — A6 row] — A6 is "Unsized" but its bounds and sizing experiment are known and absent

**Category:** DEVIATION (omission that blocks the next decision)

**Description:**
The A6 row says only "Unsized." Both bounds are computable from job 52472765's OFF leg and bracket whether A6 is worth anything: **pessimistic** (wait floors at the Min rank, 638 s) → **1.07×**; **optimistic** (imbalance fully removed, wire-only remains) → **1.67×** — a 9× spread. The decisive, cheap sizing experiment is also known and unstated: per-rank arrival trace at tick boundaries, distinguishing **systematic** lateness (same ranks always late → rebalancing works) from **jitter** (different ranks each tick → rebalancing cannot help). Without these, A6 is a name, not a phase — and committing to it unsized repeats the exact mistake A1 just paid for.

**Trigger:** Attempting to schedule or budget A6.

**Actual behavior:** No bounds, no experiment, no gate.

**Expected behavior:** The row carries the [1.07, 1.67]× band and the sizing gate; the sizing run rides with the FLOP-counter leg (known-unknown #3) so one job answers both.

**Suggested fix:**
```diff
-| **A6** | **rank rebalancing (NEW, now the main Track-A item)** | the wait is imbalance at tick boundaries: skew 16.5 %, **Max/Avg 2.18**. Only rebalancing or fewer tick boundaries touch it. Unsized. |
+| **A6** | **rank rebalancing — SIZE BEFORE BUILDING** | bounds from 52472765: **1.07×** (wait → Min-rank floor) to **1.67×** (wait → wire-only), a 9× spread. Sizing gate: per-rank tick-boundary arrival trace — **systematic** lateness (same ranks late) → proceed; **jitter** → drop A6. Ride the FLOP-counter leg (one job, two unknowns). |
```

**Test case:**
```python
def test_R007_a6_bounds():
    wall, wait, wmin = 1898.15, 764.36, 638.18
    assert abs(wall / (wall - (wait - wmin)) - 1.07) < 0.01
    assert abs(wall / (wall - wait * 0.9985) - 1.67) < 0.01
```

---

### [R-008] [LOW] [PLAN:58-66 — "Where we are" table] — Missing the current measured baseline (λ=1)

**Category:** DEVIATION

**Description:** The measured-state table still anchors on 4356 s/sim-s (λ=0.63) and omits the current production state: **3123 s/sim-s** (52422891, 2 s), independently reproduced at **3164** (52472765 OFF leg, 0.6 s; +1.3 %). Future gates' denominator is the λ=1 state; it should be a row, not prose. (The +1.3 % run-to-run agreement is itself worth recording — it is the first repeatability figure this program has.)

**Suggested fix:**
```diff
 | **MFEM LTS** | **4356** | job 52344266 |
+| **MFEM LTS, λ=1 — current baseline** | **3123** (reproduced 3164, +1.3 %) | jobs 52422891 / 52472765 |
```

---

### [R-009] [LOW] [document/comm_dev/PLAN_lts_comm_reduction_2026-07-22.md] — Falsified companion plan carries no supersession banner

**Category:** QUALITY (stale governing document)

**Description:** The comm plan still presents merge→overlap with payoff *"2341 → ~150–400 s/sim-s (1.8–2.0×)"* as a live program. Its Phase 1 was A1 (measured: nothing); Phases 2–3 are devalued by the same result. Both kernel-plan predecessors carry `⛔ SUPERSEDED` banners; this document is the only falsified plan in the tree without one.

**Suggested fix:** Prepend a banner pointing at `RESULTS_a1_merge_2026-07-24.md`; note that the Phase-0 exchange inventory and the `NbrExchangerSplit` design remain valid as *reference material* while the payoff model is falsified.

---

### [R-010] [LOW] [PLAN:44-47, 89] — Skew figures quoted without their config/window tags

**Category:** QUALITY

**Description:** Three legitimate skew measurements now exist and differ: 24.6 % (λ=0.63, 2 s, 52422891), 19.0 % (λ=1, 2 s, 52422891), 16.5 % (λ=1, 0.6 s, 52472765 OFF). The plan quotes 24.1→24.6 % in the A0 section and 16.5 % in the A6 row without tags; a reader will read them as disagreement rather than different configurations. Tag each quote with (λ, window, job).

---

## Summary
- Critical issues: 3
- Moderate issues: 4
- Low issues: 3
- Plan compliance: PARTIAL — the work-list strikes and the A1 banner are correct and current, but the plan was not swept downstream of the falsified mechanism: a known-unknown still marked RESOLVED with the falsified answer (R-001), a directly self-contradicting section pair (R-002), and an appendix endpoint overstated ~1.9× (R-003).
- Verdict: **FAIL — one documentation pass required before the plan can govern the next actions.** No code is implicated: A1's implementation, its tests, and its disposition (keep, flag OFF) all survive this review untouched. Every fix is an edit to two markdown files, addressable in a single `/code-fix` round.

## What survives the audit (so the fix round does not overreach)
- The A1 banner and work-list strike-throughs (PLAN:83–100): correct and consistent with job 52472765.
- A0's **measurements** (1.370×, rupture stability, the LTS stage split, wire share 0.15 %) — only the causal *interpretation* needs the R-004 annotation.
- The corrected model in the A1 banner ("wait scales with how often ranks must meet") — this is the statement R-001/R-002 must be made consistent WITH, not revised.
- Track B's table and gates; the B3 guard section; the Prerequisites section (the FLOP leg now doubles with A6 sizing, per R-007).
- The "Neither track suffices alone" floors (7.8× / 5.2×) — recomputed against the post-A0 state, still correct.
- The plan's structural rule (measured body / quarantined projections). The appendix failed not because the rule is wrong but because the post-A1 sweep stopped at the work list.

## Unreviewed Areas
- `RESULTS_a0_wiggle_2026-07-24.md` beyond the rounds-label issue — a dated results document, not the governing plan; its "A1–A3 remain correct" consequence text shares R-004's staleness but is historically accurate to its date.
- The A1 design doc's projection section — a dated design artifact whose own status section already records the outcome.
- Track B's technical content (the B3 guard, bench transferability) — audited by the 2026-07-22 283-agent round; nothing in A1's result touches it.
