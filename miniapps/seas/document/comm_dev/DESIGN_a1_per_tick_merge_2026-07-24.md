# A1 — per-tick exchange merge (125 → 64 rounds/sync): design + correctness proof

**Date:** 2026-07-24 · **Status: design complete, preconditions VERIFIED IN SOURCE, ready to implement.**
**Target:** 125 → 64 exchange rounds/sync at Nc=6 (comm plan Phase 1). A2 takes 64 → 32.
**Estimated payoff:** wait 1251 → ~641 s/sim-s ⇒ **~1.24× end-to-end** (3123 → 2512 s/sim-s,
13.1× → 10.5× SeisSol). Projected from A0's measured scaling law, linear case.

## Current structure

Per tick, the stepper loops over `tk.correct_clusters` (sorted ascending = fine→coarse). For **each**
correcting cluster it calls `AdvanceADERClusterBulk` → step 2c → `ComputeADERClusterSeamFaceFluxRHS`,
which does its **own** pack + `ExchangeFaceNbrData`:

| exchange | fired by | per sync (Nc=6) |
|---|---|---:|
| forecast (`wave_operator.inl:2706`) | every Correct(c), c < Nc−1 | 62 |
| batched I (`:2730`) | every Correct(c) | 63 |
| **total** | | **125** |

Round count verified against measurement: 171 syncs × 125 × 2 Waitall/round = 42,750 vs the
**42,686** calls/rank Caliper recorded in job 52422891 (0.15 %; the remainder is the truncated final sync).

## The change

Hoist both packs+exchanges out of the per-cluster call into a **per-tick pre-pass** that packs *all*
correcting clusters into one buffer and exchanges **once** per buffer:

```
per tick:
   ... all Predicts ...
   PrepareClusterSeamExchangeTick(correct_clusters, {I_[c]}, {dt_step[c]}, ...)   // 2 exchanges
   for c in correct_clusters:            // unchanged order, fine -> coarse
       AdvanceADERClusterBulk(c, ...)    // step 2c now READS the pre-exchanged buffers
```

⇒ **2 rounds/tick × 32 ticks = 64/sync**, independent of how many clusters correct.

New method (`wave_operator.hpp`, implemented in the `.inl`):

```cpp
void PrepareClusterSeamExchangeTick(
    const int *correct_clusters, int n_correct,
    const mfem::Vector *const *I_per_cluster,   // I_[c], one per correcting cluster
    const real_t *dt_step_per_cluster,
    real_t t_s, real_t dt_base, int tick, real_t T_actual,
    const real_t *dk_data, int dk_order, const int *provider_slot_of_elem,
    bool exchange_forecast) const;
```

It fills two members, `tick_fc_all_` and `tick_nbr_all_` (deep copies, since the second exchange
clobbers `FaceNbrData()` — this is why the existing code already deep-copies `fc_all`).
`ComputeADERClusterSeamFaceFluxRHS` gains a "buffers already prepared" path that skips its own
pack+exchange and reads those members. Behind an opt-in flag; default = today's per-correct path.

## Why this is bitwise identical — three preconditions, all VERIFIED IN SOURCE

**(1) The I buffer's extra blocks are never read.** The face loop classifies every shared face by
*both* sides' cluster ids (`wave_operator.inl:~2759`):

```
if      (c1 == cluster_c     && c_nbr == cluster_c)     mode = 0;   // reads nbr_all @ cluster_c
else if (c1 == cluster_c     && c_nbr == cluster_c + 1) mode = 1;   // coarse side from fc_all
else if (c1 == cluster_c + 1 && c_nbr == cluster_c)     mode = 2;   // reads nbr_all @ cluster_c
else                                                    continue;   // not incident -> skipped
```

`nbr_all` is read **only** at ghost elements whose cluster == `cluster_c` (modes 0 and 2). Mode 1's
coarse side comes from the *forecast* buffer, not `nbr_all`. Each element carries exactly one
`lts_cluster_id_`, so packing other correcting clusters writes **disjoint** slots that cluster c's
loop never touches.

**(2) The forecast buffer's extra blocks are disjoint.** Cluster `c` writes only
`seam_coarse_elems_[c+1]`. Different `c` ⇒ different index ⇒ disjoint element sets (one cluster id
per element). Each cluster keeps its own `[a,b]` sub-interval, applied to its own elements.

**(3) No correct mutates a hoisted pack's inputs.** This is the hazard that would have killed A1 —
hoisting reads *earlier* than today, so any write by an earlier-correcting cluster would change the
bytes. Verified it cannot happen:
- `I_cluster` is **`const Vector &`** in `AdvanceADERClusterBulk` (`:2194`) — a correct never
  mutates its own or any other cluster's `I`. `I_` is per-cluster storage (`lts_bulk_stepper.hpp:148`).
- `seam_coarse_dk_` has exactly **one** write site, `PrepareSeamCoarseForecast` (`:3033`), which is
  called from the **Predict** phase only (`lts_bulk_stepper.hpp:93`, `lts_fault_stepper.hpp:293`).

⇒ Every input to the merged pack is final at the end of the tick's predict phase. **Inserting the
pre-pass after all Predicts and before the first Correct reads exactly the bytes the per-correct
packs read today.**

**Ordering is untouched.** Only the pack+exchange (pure data movement) is hoisted. The per-cluster
face loops, their fine→coarse order, step 2c's scatter and step 2d's consume all stay where they are.

## Acceptance

**A1 needs BOTH a local gate and an Expanse run. They answer different questions.**

*Local (np ≤ 10, small fixture — never the production mesh):*
1. **Bitwise identity**, flag ON vs OFF, on a small multi-cluster fixture: identical output bytes.
   Bitwise identity is a *machine-independent* property, so a laptop settles it as well as a
   supercomputer — this is the primary correctness gate and A1 claims no tolerance. Pattern to follow:
   `tests/parallel/test_bimaterial_seam_fault_np2.cpp` (there is no LTS-specific parallel test yet;
   A1 must add one).
2. **Liveness** at np=2: the merged pre-pass is a collective and must fire on **every** rank at every
   tick, matched. A rank that skips it because its local `correct_clusters` looks empty hangs the job
   (the recorded R-1600 mode, `wave_operator.inl:3661-3672`). This reproduces at np=2 — no cluster
   needed. The tick table is rank-uniform, so assert that rather than assume it.
3. **Round count** — `GhostExchangeCount()` = 2/tick, and `BuildTickTable`'s `n_collectives` updated
   to match, else the driver's matched-collective `MFEM_VERIFY` aborts.

*Expanse (np=256, TPV104-200m) — REQUIRED, not optional:*
4. **Production bitwise identity** at the real scale and rank count. The local fixture cannot exercise
   the partition, the 6-cluster layout, or rank seams at np=256.
5. **The payoff.** ~1.24× is a *projection*; wait behaviour at 256 ranks sharing 2 nodes cannot be
   inferred from 10 local ranks. This project has been burned twice by exactly that transfer — the
   kernel bench **reversed its verdict** between laptop and Rome, and the stage split inverted between
   GTS-proxy and the real LTS leg. **No speed claim for A1 is admissible without this run.**
Reuse the A0 harness for (4)+(5) — same deck, same 256 ranks, flag ON vs OFF as the two legs, so the
comparison is a same-job A/B with everything else identical (the cleanliness A0 established).

## Risks

| risk | mitigation |
|---|---|
| the matched-collective invariant aborts (it hardcodes today's count) | update `n_collectives` in `BuildTickTable` for the merged mode; it is the same code that predicts the count |
| a rank with no correcting clusters skips the collective → hang | tick table is rank-uniform; assert `correct_clusters` identical across ranks in debug |
| deep-copy cost of two full halo buffers per tick | replaces 125 halo exchanges with 64 + 2 copies; the copies are local memcpy vs network round trips |
| structural skew floor (19–24 % of wait, flat through rupture) | A1 cannot remove it — if it binds, A1 lands nearer 1.15× than 1.24×. Not a correctness risk, a payoff risk. |


## Implementation status (2026-07-24) — A1 code-complete

Steps 1-5 landed: `PrepareClusterSeamExchangeTick` + prepared path (75cffba); merge-window
hooks + bulk stepper (d3bd683); merged `n_collectives` + `--lts-merge-exchanges` (8e98ae9);
bitwise ON-vs-OFF gate (1032e68); fault-interleave stepper wiring (this commit).

**Local gate PASSED — bitwise identity is a result, not an argument:**

| fixture | checks | collectives | saving |
|---|---:|---|---:|
| np=2, 2 clusters | 138/138 | 38 -> 30 | 1.27x |
| np=3, 3 clusters | 207/207 | 102 -> 62 | 1.65x |
| np=4, 4 clusters | 276/276 | 230 -> 126 | 1.83x |

trending to the predicted 125 -> 64 = 1.95x at production Nc=6.

### ~~⚠ COVERAGE GAP~~ — CLOSED 2026-07-24

`LtsFaultSyncStepper` has **no test anywhere in the repo** (`grep -l LtsFaultSyncStepper
tests/` is empty), and the production TPV104 deck runs **that** stepper, not the bulk one.
So A1's fault-path wiring is **enabled but locally unproven**. Its `BeforeCorrects` is
line-for-line the bulk stepper's, with the same forecast predicate copied from its own
`Correct`, and the merge logic it calls lives in the wave operator (which *is* covered) —
but identical code is not tested code.

**Consequence for the Expanse run:** acceptance item (4), production bitwise identity, is
**not optional and not merely a scale check** — it is the *only* evidence that will exist for
the fault path. Run the A/B as ON-vs-OFF and compare fault output bytes **before** reading any
timing. A speed number from an unverified path is worth nothing.


### Coverage gap CLOSED (2026-07-24)

`tests/unit/test_lts_a1_merge_fault_np2.cpp` now gates the fault-interleave stepper with the
same ON-vs-OFF `memcmp`. **Both LTS steppers are proven bitwise-identical**, with identical
results:

| fixture | bulk stepper | fault stepper | collectives |
|---|---:|---:|---|
| np=2, 2 clusters | 138/138 | **138/138** | 38 -> 30 (1.27x) |
| np=3, 3 clusters | 207/207 | **207/207** | 102 -> 62 (1.65x) |
| np=4, 4 clusters | 276/276 | **276/276** | 230 -> 126 (1.83x) |

The fixture deliberately has **zero fault faces**: A1 does not touch the fault half (D-2 keeps
fault faces rank-interior, so it fires no seam collectives), so a fault-free fault-stepper
exercises 100 % of A1's delta there while needing no friction physics. The mock iterator's
`Advance` aborts if ever reached, so the assumption fails loudly if the fixture drifts.

Consequence: the Expanse production-bitwise check is back to being a *scale* check rather than
the only evidence for the fault path. Still run it ON-vs-OFF before reading timing.
