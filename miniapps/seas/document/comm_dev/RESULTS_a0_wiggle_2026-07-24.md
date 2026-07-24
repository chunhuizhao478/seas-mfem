# A0 results — wiggle A/B (Track A go/no-go)

**Date:** 2026-07-24 · **Jobs:** 52416603 (0.6 s pre-rupture), **52422891 (2.0 s through-rupture, 4:20:02, COMPLETED)**
**Deck:** `jobs/comm_dev/a0_wiggle_tpv104_200m_expanse/` · TPV104-200m (2,464,689 tets), p3/ADER-O4,
2 nodes × 128 ranks, cfl 0.5, `--deriv-cache --shared-ck-recursion --face-cache`, ParaView OFF.
**Method:** two clustered-LTS legs in the same job, identical flags, differing in exactly one
non-comment config line (`lts_wiggle`). λ=0.63 is the shipped default (`"scan"`); λ=1.0 is `"off"`.
Both cluster to Nc=6, so rounds/sync = 125 is identical and only the sync *length* changes.

## Verdict

**Track A GO, on measurement.** Exposed wait scales with the exchange-round rate — super-proportionally.
**`lts_wiggle="off"` is a measured 1.370× end-to-end speedup for one config line, stable through
rupture.** Recommend making it the default for np≥256 runs (see the caveat on the objective function).

## 1. End-to-end (2 s window, through rupture)

| | λ=0.63 (default) | λ=1.0 (off) | |
|---|---:|---:|---|
| wall (`seas::spatial_dyn::step`) | 8554.3 s | **6245.2 s** | **1.370× faster** |
| s/sim-s | 4277 | **3123** | |
| behind SeisSol (239 as-run) | 17.9× | **13.1×** | |
| exposed `MPI_Waitall` | 4523.3 s (52.9 %) | 2502.7 s (40.1 %) | **1.807×** |
| exchange rounds/rank | 67,746 | 42,686 | **1.5871×** (predicted 1.5873) |
| sync intervals | 271 | 171 | |
| messages | 585 M | 229 M | |

The round-count prediction reproduces to **4 decimal places**, so the traffic model is exact.
Both windows agree: the 0.6 s leg gave 1.377×/R_wait 1.838, the 2 s leg 1.370×/R_wait 1.807.

## 2. Stability through rupture — **PASS**

This was the question the 2 s run existed to answer, because λ=1 runs *at* the CFL limit while
λ=0.63 ran 37 % below it.

| | λ=0.63 | λ=1.0 |
|---|---:|---:|
| final V_max | 17.83 m/s | 17.6871 m/s |
| completion | t = 2 s, 271 syncs | t = 2 s, 171 syncs |
| NaN / Inf / ABORT / divergence | none | none |

**V_max agrees to 0.80 %** through a fully developed rupture (V_max rose 0.098 → 17.3 → 17.8 m/s).
Both legs ran the same `dt_cfl = 3.66128e-4 s`; λ scales the LTS *base* step, not the CFL limit, so
λ=1 is at the scheme's intended margin, not beyond it. **No stability objection to λ=1.**

*Not established here:* agreement of the on-fault station waveforms against the SCEC TPV104
reference. V_max agreement is a strong smoke test, not a benchmark-physics acceptance.

## 3. Skew through rupture — **the retraction is resolved**

An earlier report of "skew 24 % < 40 %, comm-plan payoff holds" was measured at `tfinal=0.5 s`,
**before** nucleation (t₀ = 1.0 s), where fault work is near zero *by construction*. That number was
withdrawn as not a valid bound. This is the honest measurement:

| | 0.6 s (pre-rupture) | 2.0 s (through rupture) |
|---|---:|---:|
| λ=0.63 skew `(Avg−Min)/Avg` | 24.1 % | **24.6 %** |
| λ=0.63 `Max/Avg` | 1.73 | **1.72** |
| λ=1.0 skew | 20.0 % | **19.0 %** |
| λ=1.0 `Max/Avg` | 2.22 | **2.19** |

**Skew is flat through rupture** (24.1→24.6 %, 1.73→1.72). The imbalance is therefore **structural —
a property of the cluster/partition layout — not a rupture-driven fault-work effect.** This both
retires the withdrawn gate *and* weakens the hypothesis that fault-work imbalance is the main
residual: if it were, skew would grow when the fault becomes active. It does not.

## 4. The wait is **not** bandwidth — it is synchronization exposure

`MPI_Isend` + `MPI_Irecv` total **6.9 s of the 4523.3 s wait = 0.15 %**. Wire time is negligible.

Corroborating comparison: GTS pushes ~35,500 halo rounds/sim-s against LTS's ~16,900 — **twice as
many** — yet GTS's exposed wait is ~0.00002 % of wall. Same code, same mesh, same ranks; LTS sends
*fewer* messages and waits vastly more.

**Diagnosis: the LTS wait is ranks arriving at sync points at different times, not data on the wire.**
That is why removing 1.587× of rounds removed 1.807× of wait — what is saved is *exposure events*,
not bytes. Track A's merge levers (A1–A3) remain correct, but the mechanism is "fewer places to
stand around", not "fewer bytes". A4 (sparse payload) is accordingly de-prioritised: it attacks
0.15 % of the wait.

## 5. First measured LTS stage split — **the plan's GTS proxy was wrong**

Made possible by commit `884349a` (the LTS step region + corrector scope). Every Track-B per-stage
number in the plan was previously a GTS proxy borrowed from a different leg.

| stage (% of step) | λ=0.63 | λ=1.0 |
|---|---:|---:|
| predictor `…SubStepStatesAndIntegralCluster` | 18.65 % | 27.89 % |
| bulk corrector `AdvanceADERClusterBulk` | 80.13 % | 70.69 % |
|  └ seam `ComputeADERClusterSeamFaceFluxRHS` | 64.87 % | 47.94 % |
|    └ of which exposed wait | 52.86 % | 40.05 % |
|    └ of which compute | 12.01 % | 7.89 % |
| **true compute (step − wait)** | **47.1 %** (2015 s/sim-s) | **59.9 %** (1871 s/sim-s) |

**The predictor is 39.6–46.6 % of true compute, not the 56.2 % the plan assumed from the GTS leg.**
Track B's B1/B2 (predictor work) therefore address a smaller share than written, and the
seam-corrector's *compute* half (12.0 % / 7.9 % of step) is a newly visible target that appears in
no plan phase.

Coincidence worth noting: the λ=0.63 true-compute figure is **2015 s/sim-s**, exactly the value the
plan carried as its (derived) compute term. That is a genuine independent confirmation of the budget's
compute side — arrived at by direct measurement rather than by subtraction.

## 6. Consequences for the plan

1. **A0 is DONE and its assumption is retired.** "Whether exposed wait scales with round count" was
   known-unknown #1; it does, at 1.81× per 1.59× of rounds.
2. **`lts_wiggle="off"` should become the default for np≥256** — 1.370× for one line, stable through
   rupture. *But the principled fix is different:* `lts_clustering.cpp:compute_cost` minimises element
   updates with **no communication term**, which is why it picked a λ that costs 1.37×. Hardcoding
   λ=1 is right for this configuration and may be wrong at another rank count or mesh. Add a comm
   term to the objective (new item **A5**).
3. **Known-unknown #4 (comm skew) is resolved** and points *away* from fault imbalance: skew is flat
   through rupture, so it is structural. Known-unknown #5 (which rank is on the critical path) stays
   open and is now the more interesting question.
4. **Track B rescales.** The predictor is a smaller share of compute than assumed; re-derive B1/B2
   gates against this measured split, and consider the seam corrector's compute half.
5. **A4 de-prioritised** — wire time is 0.15 % of the wait.

## Reproduce

```bash
sbatch -t 8:00:00 --export=ALL,A0_TFINAL=2.0s \
  jobs/comm_dev/a0_wiggle_tpv104_200m_expanse/run_a0_wiggle_expanse.sbatch
```

Caliper reports: `runs/52422891/run_lambda{063_baseline,100_wiggleoff}.cali-region-report.txt`.

*Extractor note:* the job's auto-parsed `R_wait` was wrong on the first run (52416603) — a report
holds several `MPI_Waitall` rows and it took a negligible 11-call row instead of the real 67,746-call
one. Fixed in `69c5fd8` (selects the largest-Avg row); the 52422891 output above is from the fixed
version and its printed verdict is correct.
