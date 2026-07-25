# A1 results — per-tick exchange merge: CORRECT, but NO SPEEDUP

**Date:** 2026-07-24 · **Job 52472765** (COMPLETED, 1:17:40) · TPV104-200m, np=256, λ=1, 4 legs
**Verdict: A1 is bitwise correct at production scale and delivers ZERO speedup. Do not enable it.
The result invalidates the mechanism behind A2 and A3 as well.**

## What the run says

**Phase 1 — bitwise identity at np=256: PASS.** 6/6 fault-output files byte-identical, 0 differing,
`V_max_global = 1.00359e-16` identical on both legs. The design's correctness argument holds at
production scale, on the fault-interleave stepper, with the real 6-cluster layout and partition.

**Phase 2 — timing.** Both legs: same deck, same 52 syncs, one CLI flag apart.

| | OFF (per-correct) | ON (merged) | ratio | predicted |
|---|---:|---:|---:|---:|
| exchange rounds/rank | 12,812 | 6,556 | **1.954×** ✓ | 1.953× |
| exposed `MPI_Waitall` | 764.4 s | 787.4 s | **0.971×** ✗ | ~1.8× |
| wall (`step`) | 1898.1 s | 1920.9 s | **0.988×** ✗ | ~1.24× |
| wait share of wall | 40.2 % | 40.8 % | — | — |

**The merge did exactly what it was built to do — it halved the rounds — and the wait did not move.**
(The 1.2 % wall regression is at the edge of run-to-run noise; the robust statement is "no gain".)

## Why — and why A0 worked

Per-round wait went from **59.7 ms to 120.1 ms, a 2.01× increase**: each surviving round absorbs
precisely the wait the removed rounds used to carry. **The wait was never *in* the rounds.**

The distinguishing variable between the two experiments:

| | syncs | rounds/sync | rendezvous points in time | wait |
|---|---:|---:|---|---:|
| **A0** (λ 0.63→1) | 271 → **171** | 125 → 125 | **fewer** | **1.807× less** |
| **A1** (merge) | 52 → 52 | 125 → **64** | **unchanged** | 0.971× (flat) |

A0 lengthened `dt_base`, so the ranks had to meet *less often per simulated second*. A1 kept every
tick boundary and merged exchanges that were already happening back-to-back **at the same boundary,
with no compute between them** — so it removed messages but removed no rendezvous. Ranks that were
already stalled waiting for the slowest peer keep waiting exactly as long.

**Corrected model of the LTS wait:** it scales with **how often ranks must synchronise**, not with
how many messages they exchange when they do. That is consistent with everything measured so far —
A0's gain, the 0.15 % wire-time share, and GTS's ~2× higher round rate at ~0 wait — and it is the
first time the two hypotheses have been separated by an experiment that varies one and not the other.

## Consequences — A2 and A3 are devalued by the same mechanism

- **A2 (64 → 32, one wait phase per tick)** merges the two remaining rounds *within a tick*. That is
  the same operation A1 just showed does nothing. **Expected payoff ≈ 0.** Do not build it on the
  strength of the old model; it needs its own justification or it should be dropped.
- **A3 (split-post overlap)** hides wire time behind compute — but wire time is **0.15 %** of the
  wait (A0). There is nothing to hide. **Devalued.**
- **A4** was already de-prioritised for the same reason.
- The projected Track-A recovery of ~2041 s/sim-s (later 611 for A1) **has no measured support and
  should not be carried in any plan.** The only Track-A lever with a measured payoff is **A0**
  (banked, 1.370×).

**What remains, and it is now the main thing:** the wait is rank imbalance exposed at tick
boundaries. Skew is 16.5 % of the wait with **Max/Avg = 2.18** — the slowest rank waits over twice
the mean. Only two things touch that: **rebalancing** (make ranks arrive together) and **fewer tick
boundaries** (the A0 lever, already taken). Neither is in the current Track-A phase list.

## Disposition of the A1 code

**Keep, do not enable.** `--lts-merge-exchanges` stays opt-in and OFF. It is proven bitwise-identical
(local: both steppers, 2/3/4 clusters; production: np=256), it cuts MPI message volume 1.95× at no
correctness cost, and it costs ~1 % wall — so it is a genuine option if a future machine or scale is
message-rate-bound rather than skew-bound. Nothing about it is wrong; the *premise* it was built on
was wrong, and the run is what established that.

Cost of finding out: one 1:17 job, after ~1 day of implementation that also produced the first tests
`LtsFaultSyncStepper` has ever had.

## Honest note on the estimate

The plan projected 611 s/sim-s for A1 from A0's scaling law. That extrapolation stretched **one data
point** across a **different lever** (syncs → rounds-within-a-sync) — a caveat the design doc did
flag, and it is exactly where the projection broke. The lesson is not "projections are bad" but that
a projection which changes the *mechanism* being varied is not an extrapolation at all; it is a new
hypothesis, and it needed the experiment it just got.
