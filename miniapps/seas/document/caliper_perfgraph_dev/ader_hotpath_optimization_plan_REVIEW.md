# Code Review: ADER hot-path optimization plan — 2026-05-31

Adversarial review of an *optimization plan* (no implementation yet). Mandate
(from the request): (1) confirm the bottleneck is real, (2) verify the proposed
fixes cannot change the physics results, (3) add tooling to measure the speedup.
Every claim below was checked against the source in this worktree, not against
the plan's prose.

## Review Scope
- Plan: `document/caliper_perfgraph_dev/ader_hotpath_optimization_plan_2026-05-31.md`
- Source verified: `dynamic/wave_operator.{hpp,inl}`, `drivers/spatial_dyn_driver.cpp`,
  `Makefile`, `tests/unit/test_wave_operator_spatial_derivative.cpp`
- Domain context: `miniapps/seas/CLAUDE.md` (extreme-care files, no-local-mesh,
  Brent friction), memories (no-local-full-mesh, makefile-no-header-deps)
- Tools added this review (the third mandate): `tests/unit/bench_ader_hotpath.cpp`
  (+ Makefile `make bench`), `scripts/perfgraph_speedup.py`

## Bottleneck: CONFIRMED (against source AND by measurement)
- `ApplySpatialDerivative` (`wave_operator.inl:892`) recomputes
  `CalcShape`/`CalcPhysDShape`/`Tr->Weight()` per QP per call and re-derives the
  stiffness action by quadrature; all of it is geometry-only (constant in time).
  **Confirmed.**
- It is invoked from TWO independent CK recursions per macro-step —
  `ComputeADERSubStepStates` (`:1201`, via driver `:427`) and
  `ComputeADERTimeIntegrated` (`:1101`, via `AdvanceADER` `:5374`, driver `:473`)
  — both called with the **same `Q`, `dt_step`, `ader_order`** (driver lines
  427 & 473). The imposed-fault states set between them (`:467`) feed only
  `ComputeADERFaceFluxRHS`, **not** the recursion. **Confirmed: the two
  recursions are redundant and mergeable.**
- **Empirical confirmation** (`seas_bench_ader_hotpath 2 3 12`, p2/O3, 10368
  tets): `ApplySpatialDerivative` = 166 ms ×3-dir; `ComputeADERSubStepStates`
  = 334 ms ≈ 2×ASD (O3 ⇒ 2 CK levels); `AdvanceADER` = 480 ms; macro-step
  aggregate = 813 ms/step, of which the two recursions' ASD ≈ 668 ms (**~82%**,
  even higher than the perfgraph's 67.6% because the micro-bench has no
  fault/MPI overhead). The two recursions cost ≈ equal (334 ms each) — the
  redundancy Lever 2 targets is real and measured.

So the plan's diagnosis is sound. The findings below are guardrails that the
implementation MUST honor so the optimization does not (a) deadlock or
(b) silently change the physics — plus gaps in the plan's verification/cost model.

## Findings

### [R-001] [CRITICAL] [spatial_dyn_driver.cpp:427 / wave_operator.inl:1705] — `ComputeADERSubStepStates` must stay UNCONDITIONAL on every rank; its output feeds a per-substep MPI collective. Skipping it on "fault-less" ranks deadlocks the run.

**Category:** ASSUMPTION / DEADLOCK

**Description:**
`Q_per_node` (output of `ComputeADERSubStepStates`, driver `:427`) is consumed
by `EvaluateBulkAtFaultQPsCanonical` (driver `:437`), which at `np>1` performs a
per-substep **pairwise collective** `q_gf.ExchangeFaceNbrData` that **every rank
with face neighbours must participate in, even ranks with zero local fault QPs**
(explicit invariant R-1600, `wave_operator.inl:~1726`; also: a fault-less rank
still owns elements adjacent to a *neighbour's* fault face and must supply its
`Q_per_node` there). A natural-looking optimization — "most ranks own no fault
face at np=800, so skip `ComputeADERSubStepStates` there" — is therefore a
**deadlock trap**: the skipping ranks never enter `ExchangeFaceNbrData` while the
fault-owning ranks block in it.

The plan's Lever 2 wording ("compute `D(k)` once, form both outputs") could be
mis-implemented as a fault-less skip. It must instead **merge** the two
recursions, run unconditionally on all ranks, and only eliminate the *second*
(redundant) recursion.

**Trigger:** np>1 with any rank that has no local fault QPs (the common case at
np=800), if a fix makes `ComputeADERSubStepStates` (or its `Q_per_node`)
rank-conditional.

**Actual behavior (if mis-fixed):** hang at `ExchangeFaceNbrData` in
`EvaluateBulkAtFaultQPsCanonical`.

**Expected behavior:** `ComputeADERSubStepStates`/`Q_per_node` produced on all
ranks every step; only the redundant second recursion (`ComputeADERTimeIntegrated`
inside `AdvanceADER`) is removed.

**Suggested fix (Lever 2, collective-safe shape):**
```diff
  // driver, per macro-step — UNCHANGED call sequence, all ranks:
- wave.ComputeADERSubStepStates(Q, dt_step, ader_order, tau_nodes, Q_per_node);
+ wave.ComputeADERSubStepStatesAndIntegral(Q, dt_step, ader_order, tau_nodes,
+                                          Q_per_node, /*out*/ I_buf);  // one CK recursion
  ... EvaluateBulkAtFaultQPsCanonical(Q_per_node[o], ...)   // collective, unchanged
  ... friction Advance ...
- wave.AdvanceADER(Q, dt_step, ader_order, Q_new);
+ wave.AdvanceADER(Q, dt_step, ader_order, Q_new, /*I_precomputed=*/&I_buf); // skips rec #2
```
`ComputeADERSubStepStatesAndIntegral` runs the existing element-local recursion
ONCE and accumulates both `Q_per_node` (τ-node weights) and `I` (dt-integral
weights). No MPI inside it; the collective ordering (`ExchangeFaceNbrData` then
the `AdvanceADER` corrector collectives) is unchanged.

**Test case:**
```python
def test_R001_substepstates_unconditional_no_deadlock():
    # mpirun -np 4 on a TPV31-like fixture where >=1 rank owns 0 fault faces.
    # Assert: the merged-recursion build completes >=2 macro-steps without hang
    #         (timeout => FAIL), and a deliberately fault-less-skipping variant
    #         is REJECTED by a guard / documented as forbidden.
    run_mpi("seas_test_ader_substep_collective_np4", np=4, timeout_s=120)
```

---

### [R-002] [MODERATE] [wave_operator.inl:892 ApplySpatialDerivative] — Lever 1 (precompute `D_d^e = M_e⁻¹K_d^e`) is NOT bit-reproducible end-to-end; the plan's acceptance must use a physical tolerance, and Lever 2's "byte-identical" gate must baseline against a Lever-1-ON build.

**Category:** BUG (verification gap that can mask a real regression)

**Description:**
Today: `dQ_c = M⁻¹·(K_d·Q_c)` evaluated on the fly. Lever 1: `dQ_c = D_d·Q_c`
with `D_d = M⁻¹·K_d` formed by a matrix product at setup. These differ in
floating point (matrix-product re-association + the setup `Mult` rounding). The
plan correctly claims **machine-eps per call** (≤1e-12), not bit-identity — but
two consequences are unstated:
1. Over millions of steps of a **nonlinear** rupture, round-off diverges:
   late-time station traces from a Lever-1 build will NOT be bit-identical and
   may differ visibly while both remain physically valid. Acceptance must be a
   **physical** tolerance (rupture-front arrival, peak `V`, final slip within a
   stated %), never bit-identity, for any Lever-1 build.
2. The plan sequences Lever 1 → Lever 2 and calls Lever 2 "byte-identical." That
   holds only relative to the **Lever-1-ON** baseline. If the Lever-2 parity test
   compares against the *original* (pre-Lever-1) output it will falsely fail on
   Lever-1 round-off. The test gate must rebaseline after Lever 1.

**Trigger:** end-to-end run (or a multi-thousand-step integration test) of a
Lever-1 build; or a Lever-2 parity test that uses a pre-Lever-1 reference.

**Actual behavior:** plan implies traces "match"; an implementer may assert
bit-identity and either (a) fail a correct optimization or (b) pick a tolerance
that hides a real bug.

**Expected behavior:** explicit per-lever tolerances — Lever 1: unit test
≤1e-12 per call AND end-to-end physical metrics within stated %; Lever 2:
bit-identical vs the Lever-1-ON baseline.

**Suggested fix (plan §5 amendment):**
```diff
- 1. Lever 1 equivalence (numerical): ... matches the current quadrature-based
-    result to machine precision (≤ 1e-12 relative) on a small multi-element mesh
+ 1a. Lever 1 unit equivalence: |D_d·Q - M⁻¹(K_d·Q)| ≤ 1e-12·‖·‖ per call.
+ 1b. Lever 1 end-to-end: NOT bit-reproducible.  Acceptance = physical metrics
+     (rupture time, peak V, final slip) within <tol>% vs baseline — define <tol>.
+ 2.  Lever 2 parity: bit-identical vs the LEVER-1-ON baseline (not the original).
```

**Test case:**
```python
def test_R002_lever1_machine_eps_and_lever2_bitexact_vs_lever1():
    dq_fused   = apply_spatial_derivative_fused(Q)      # D_d·Q
    dq_onthefly= apply_spatial_derivative_quad(Q)       # M^-1 (K_d Q)
    assert rel_err(dq_fused, dq_onthefly) <= 1e-12      # 1a
    qnew_l2 = advance_merged(Q)        # Lever 2 on top of Lever 1
    qnew_l1 = advance_separate_lever1(Q)
    assert bitwise_equal(qnew_l2, qnew_l1)              # 2 (vs Lever-1 baseline)
```

---

### [R-003] [MODERATE] [wave_operator.inl:1267-1306 & 1146-1182] — Lever 2 merge must preserve the TWO distinct factorial-weight phases (R-1503); a naive loop-merge mis-weights `I` or `Q_per_node`.

**Category:** BUG (would silently change results)

**Description:**
The two recursions share `D(k)` but apply **different** Taylor weights, with the
factorial parameterised at different phases (the code explicitly warns against
"harmonising" them — `wave_operator.inl:1283`, `:1161`):
- `ComputeADERSubStepStates`: `fac=1.0`, `denom=k+1`, k=0 term `Q_per_node += D(0)`
  pre-added unscaled, then `Q_per_node[o] += (τ_o^k/k!)·D(k)`.
- `ComputeADERTimeIntegrated`: `fac=dt`, `denom=k+2`, k=0 term `I += dt·D(0)`,
  then `I += (dt^{k+1}/(k+1)!)·D(k)`.
A merged routine MUST keep BOTH accumulations with their own `fac`/`denom`/k=0
handling inside the single `k`-loop. Collapsing to one denominator (the obvious
"cleanup") changes the integral or the nodal states → wrong predictor → wrong
physics, on EVERY step.

**Trigger:** any Lever-2 implementation that unifies the factorial update.

**Actual/Expected:** see weights above; expected = both phases preserved.

**Suggested fix (merged kernel skeleton):**
```cpp
D_curr = Q;
for (o) Q_per_node[o].Add(1.0, D_curr);        // SubStep k=0 (unscaled)
I.Add(dt, D_curr);                              // Integral k=0 (dt-scaled)
real_t facI = dt; std::vector<real_t> facS(O, 1.0);
for (int k = 0; k < order-1; ++k) {
  D_next = 0; for (d) { ApplySpatialDerivative(d,D_curr,t); ApplyElementJacobian_(d,t,D_next,-1.0); }
  for (o) { facS[o] *= tau[o]/(k+1); Q_per_node[o].Add(facS[o], D_next); }   // denom k+1
  facI *= dt/(k+2); I.Add(facI, D_next);                                     // denom k+2
  mfem::Swap(D_curr, D_next);
}
```

**Test case:**
```python
def test_R003_merged_matches_both_separate_bitexact():
    qpn_m, I_m = compute_substepstates_and_integral(Q, dt, order, tau)
    qpn_s      = compute_substepstates(Q, dt, order, tau)        # current
    I_s        = compute_time_integrated(Q, dt, order)           # current
    for o: assert bitwise_equal(qpn_m[o], qpn_s[o])
    assert bitwise_equal(I_m, I_s)                               # orders {2,3,4}
```

---

### [R-004] [MODERATE] [POSSIBLE] [wave_operator setup, near :5552] — Lever 1 cache has no memory bound; high order / dense partition can OOM with no fail-loud fallback.

**Category:** EDGE_CASE / ASSUMPTION

**Description:**
`elem_deriv_op_` costs `3·ne_local·ndof²·8 B`. The plan tabulates p2 (≈18 MB/rank)
and p3 (≈72 MB/rank) but stops there. p4 tet (`ndof=35`) ≈ 29.4 kB/elem → at
15k elem/rank ≈ 440 MB/rank ≈ 22 GB/node at 50 ranks/node — and a denser
partition pushes further. There is no proposed runtime guard, so a high-order or
imbalanced run silently OOMs a Frontera allocation.

**Suggested fix:** at cache-build, compute the footprint and either fail loud
with the number, or fall back to the on-the-fly kernel above a budget:
```cpp
const size_t bytes = size_t(3)*ne_*size_t(ndof_per_el_)*ndof_per_el_*sizeof(real_t);
MFEM_VERIFY(bytes <= deriv_cache_budget_bytes_,
            "elem_deriv_op_ cache would use " << bytes/1e6 << " MB/rank > budget "
            << deriv_cache_budget_bytes_/1e6 << " MB; raise --deriv-cache-mb or "
            "build without the precomputed-derivative cache.");
```

**Test case:**
```python
def test_R004_cache_budget_guard():
    # construct at p4 on a mesh sized so the cache exceeds a tiny budget;
    # assert MFEM_VERIFY aborts with the byte count (not a bad_alloc).
    assert aborts_with_message(build_wave_op(order=4, budget_mb=1), "elem_deriv_op_ cache")
```

---

### [R-005] [LOW] — Perfgraph's actual mesh order (p2 vs "p3") is unconfirmed; ndof-dependent estimates hinge on it.

**Category:** QUALITY (doc clarity)

**Description:** the sbatch is `tpv31_p2_aderO3` (FE order 2) but the request
called it "p3"; the plan hedges "p≥2". The memory/speedup numbers scale with
`ndof²`, which differs 2.5× between p2 (ndof=10) and p3 (ndof=20).
**Suggested fix:** confirm `[mesh].order` of the profiled run and quote one
number; the new benchmark resolves it empirically:
`./seas_bench_ader_hotpath 2 3` vs `./seas_bench_ader_hotpath 3 3`.
**Test case:** n/a (doc).

---

### [R-006] [LOW] [wave_operator.inl:3539-3553] — Lever 4.1 (hoist the per-step consensus `Allreduce`) is verified SAFE, but the hoist must keep fail-loud semantics.

**Category:** DEVIATION (guardrail on a confirmed-good item)

**Description:** confirmed `has_bulk_bg_` is set once at setup
(`SetAbsorbingBackground`, `spatial_dyn_driver.cpp:1656`, before the time loop)
and never mutated during the run — so hoisting the `MPI_Allreduce(MIN)` consensus
check out of the per-step `ComputeADERFaceFluxRHS` to a one-time post-setup check
is correct. The fix must (a) still abort on ANY rank lacking the background,
(b) run the single check after the collective point where all ranks have called
`SetAbsorbingBackground` — not merely delete the per-step check.
**Suggested fix:** move the existing R-1505 Allreduce block verbatim into a
`VerifyBackgroundConsensus()` called once after setup; assert in
`ComputeADERFaceFluxRHS` (debug-only) that it ran.
**Test case:**
```python
def test_R006_missing_background_still_fails_loud_after_hoist():
    # np=2, one rank skips SetAbsorbingBackground -> hoisted check aborts BOTH ranks.
    assert mpi_run_aborts("seas_test_bg_consensus_np2", np=2, msg="SetAbsorbingBackground")
```

---

## Measurement tooling added (third mandate)
1. **`tests/unit/bench_ader_hotpath.cpp`** + Makefile (`make bench` →
   `seas_bench_ader_hotpath`). Local small-tet micro-benchmark of
   `ApplySpatialDerivative` / `ComputeADERSubStepStates` / `AdvanceADER`; prints
   ns/call and a macro-step aggregate. **Built + ran**; baseline (p2/O3, 10368
   tets): ASD×3 = 166 ms, SubStepStates = 334 ms, AdvanceADER = 480 ms, aggregate
   = 813 ms/step. Run baseline vs optimized with identical args; ratio = speedup.
   (Small mesh — honors the no-local-full-mesh rule; it is a unit-scale timer.)
2. **`scripts/perfgraph_speedup.py`** — diffs two Caliper `runtime-report` dumps
   (the format pasted from the `*.err`), summing duplicate region names (so the
   two `ApplySpatialDerivative` subtrees combine to 953.25 s) and printing
   per-region speedup + a TOTAL-step line. Validated on the supplied perfgraph.
   Use `--metric max` for a wall-time (critical-path) proxy.

## Summary
- Critical issues: 1 (R-001 — deadlock trap; the single most important guardrail)
- Moderate issues: 3 (R-002, R-003, R-004)
- Low issues: 2 (R-005, R-006)
- Plan compliance: N/A (plan-stage review; no implementation yet)
- **The bottleneck claim and the results-preserving claim are CONFIRMED** for
  Levers 1 & 2 **provided** R-001 (no rank-conditional skip), R-002 (physical
  tolerance + correct rebaseline), and R-003 (preserve both weight phases) are
  honored. Lever 4.1 is confirmed safe (R-006).
- Verdict: **PASS WITH FIXES** — amend the plan with R-001…R-004 before
  implementation; then proceed with Lever 1 (test-first), measuring with the two
  tools above.

## Unreviewed Areas
- Lever 1's batched-gemm layout variant (gather/scatter for the component-major
  `Q`) — correctness identical to the mat-vec form, but its own perf claim is
  unmeasured; the micro-benchmark will quantify it when implemented.
- Lever 3 (`ComputeVolumeRHS` geometry reuse) and Lever 4.4 (partition
  rebalancing) — lower priority; not deeply audited here.
- The bimaterial path's `ApplyElementJacobian_` interaction with a cached `D_d^e`
  was reasoned (orthogonal) but not exercised; the parity test in the plan's §5
  (run through `BimaterialWaveOperator`) remains required.
