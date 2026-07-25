# GTS speed plan — MFEM p3, TPV104-200m

**Date:** 2026-07-25 · **Scope: GTS only.** LTS, comm merging and rebalancing are out of scope here;
GTS strips the LTS machinery away and leaves pure per-element kernel cost, which is where the
measured gap lives. **Companion:** `PLAN_performance_2026-07-25.md` (whole program).

## The one table this plan is built on

**MFEM p3 GTS = 180.25 µs-core/update; SeisSol o4 = 12.01; gap 15.0×** (jobs 52495071 / 52495072).
Stage shares are the measured GTS split (B0, job 52379131 Leg 2, same lever set):

| stage | share | µs-core | vs SeisSol's **entire** step |
|---|---:|---:|---:|
| **predictor (ADER-CK)** | 56.2 % | **101.4** | **8.4×** |
|  └ `ApplySpatialDerivative` | 27.1 % | 48.9 | |
| face interior | 13.0 % | 23.4 | 1.9× |
| face shared | 10.6 % | 19.0 | 1.6× |
| volume | 7.9 % | 14.2 | 1.2× |
| unattributed | 11.7 % | 21.1 | 1.8× |
| fault/friction | 0.65 % | 1.2 | 0.1× |

**The structural fact that orders everything below:** the predictor alone is **8.4× SeisSol's whole
step**. Zero every other stage and MFEM is still 8.4× behind; zero the predictor and keep the rest,
6.6×. There is no single-stage fix — but there is no ambiguity about ordering either. A 4× win on the
predictor moves 15.0× → 8.7×. The same 4× on faces moves it to 13.6×.

---

## G0 — the deciding measurement (one job, two questions, no code)

**Why first.** Everything downstream forks on two unknowns, and both are answered by one short GTS
run. Building B2 before this is how A1 failed: committing engineering to a mechanism nobody measured.

### G0a — is MFEM FLOP-**volume**-bound or FLOP-**rate**-bound?

*What:* `perf stat` counters on GTS legs at **p1 and p3** (both, not one), reporting FLOPs/update and
achieved GFLOP/s/core.

*Why it decides the plan:* SeisSol's achieved rate *rises* 2.62 → 4.54 → 7.84 across o2→o4 while its
per-step cost rises only 2.59× — so its FLOP volume grows ×7.75 and it absorbs that by vectorising
better. Two possible MFEM stories, with opposite fixes:
- **Rate-bound** (we do similar FLOPs, slowly) → traffic reduction is the right tool → **B1/B2/G3
  proceed as written**.
- **Volume-bound** (we do far more FLOPs) → tiling cannot help; the answer is M1/M3, *element-local
  dense blocks and per-element GEMM over precomputed operators* — **which is in no phase of any
  current plan** and would be a larger redesign.

*Gate:* report both numbers at both orders. **Blocking** — no ownership of the gap may be assigned to
B1/B2 until it returns.

### G0b — is the unconditional `SetCurvature` a real cost?

*What:* the same job, with an arm that skips `pmesh.SetCurvature(order)`.

*Why:* `spatial_dyn_driver.cpp:1733` applies it **unconditionally to a straight-sided tet mesh**, so
MFEM evaluates a 4/10/20-node H1 nodal transformation per quadrature point where SeisSol is always
affine. It is **order-dependent** (4→10→20 nodes), so it sits inside the measured order law, and on a
straight-sided mesh it may be removable outright.

*Interpretation rule, stated in advance:* a **null result with the levers on does NOT mean "geometry
is not in the gap."** Under `--deriv-cache` both volume paths are quadrature-free and `--face-cache`
removes interior-face transforms, leaving only fault/shared/boundary faces to pay it. Run **one arm
with levers off** to price the raw confound and one with levers on to price what remains.

*Gate:* p3 per-update compute falls, or the confound is priced and closed.

**Cost:** ~2 short GTS legs at p1 + p3 ≈ 40–80 core-h. Trivial against the 978 already spent.

---

## G1 — fuse the predictor's accumulation sweeps (start now, in parallel with G0)

*What:* the predictor's `O+1` full-length AXPYs and zeroing per expansion level are folded into the
loop that already touches the data.

*Why this one is safe to start before G0 returns:* it removes *measured* traffic — `Vector::Add`
18.9 % + `Vector::operator=` 9.3 % = **28 % of perf self-time in pure whole-vector movement** — and
that is a win under *either* G0a outcome, because it deletes work rather than re-scheduling it. It is
also the smallest change in the program.

*Gate:* ≥1.2× on the predictor stage; ≤1e-12 parity (round-off lever convention, opt-in flag,
defaults bit-identical). **Stop and reassess below 1.1×.**

*What it is worth if it lands at its gate:* predictor 101.4 → 84.5 µs-core, total 180 → 163 = 15.0× →
13.6×. Real but modest — G1 is a down payment on G2, not a headline.

---

## G2 — restructure the predictor (the only item that can halve the headline)

*What:* cache-blocked tiles (E ≈ 32–128 elements) with all expansion levels fused, so each element's
operators are read once per step instead of ~40 times; optionally per-element small GEMM.

*Why:* the microbenchmark measured 2.3–6.3× (tiled) and 2.8–8.3× (dgemm) on this exact structure
**under Rome contention** — where the current structure degrades 33 → 79 µs/element while the tiled
one stays flat at ~13.

*Gated on G0a.* If MFEM is volume-bound, this mechanism is wrong and the item is replaced by M1/M3.

*Note G1 and G2 do NOT add* — same stage, counted once. G1 is a cheaper route to part of G2's win.

*What it is worth at the bench central (4×):* predictor 101.4 → 25.4, total 180 → 104 = **15.0× →
8.7×.**

---

## G3 — faces and volume (only after the predictor is done)

Face interior 23.4 + face shared 19.0 + volume 14.2 = 56.6 µs-core combined. Worth doing, but note
**all three together are smaller than the predictor alone**, and each is only 1.2–1.9× SeisSol's
entire step. Same tiling mechanism; same G0a gate. B3 (face-cache → LTS corrector) is *not* in this
plan — it is an LTS item.

---

## G4 — the unattributed 11.7 %

21.1 µs-core is unaccounted for by the named regions. That is larger than the volume stage. Before
optimising it, **instrument it** — the sub-stage Caliper scopes inside the corrector (interior-face,
volume, friction) do not exist yet. Cheap, and it may reveal a target bigger than G3.

---

## Honest ceiling for GTS-only work

| state | µs-core | vs SeisSol |
|---|---:|---:|
| today | 180.25 | 15.0× |
| G1 at gate | 163 | 13.6× |
| G2 at bench central | 104 | 8.7× |
| G2 + G3 + G4 all at target | ~55–70 | **4.6–5.8×** |
| perfect predictor, rest untouched | 78.9 | 6.6× |

**Parity is not reachable by these steps.** Getting to 12 µs-core needs every stage on generated GEMM
chains over precomputed operators — SeisSol's M1/M2 architecture — not incremental tuning of a
quadrature-based path. This plan's realistic destination is **~5×**, and it says so up front.

## Sequencing

1. **G0** (one job, both questions) **and G1** (local code) start together — G0 is measurement, G1 is
   the change that pays under either answer.
2. **G2** only after G0a says rate-bound.
3. **G4 instrumentation** before **G3**, so G3 targets the right thing.

Nothing beyond G0/G1 is authorised by current evidence.
