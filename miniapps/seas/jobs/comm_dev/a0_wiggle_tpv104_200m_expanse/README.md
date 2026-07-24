# A0 — wiggle A/B (Track A go/no-go)

**Plan:** `document/kernel_dev/PLAN_performance_program_2026-07-23.md`, Track A phase **A0**.
**One-line experiment, highest value in the plan.**

## What this decides

Two clustered-LTS runs of SCEC TPV104-200m (2.46 M tets, p3/ADER-O4), in the **same job**, on
the **same 256 ranks**, with **identical flags**, differing in **exactly one config line**:

| leg | deck | `lts_wiggle` | λ | Nc |
|---|---|---|---:|---:|
| baseline | `../../lts_phase5/tpv104_200m_lts_speed_expanse/tpv104_200m_lts_rate2.toml` | *(unset → "scan")* | 0.63 | 6 |
| wiggle-off | `tpv104_200m_lts_wiggleoff.toml` | `"off"` | 1.0 | 6 |

`lts_nc_cap` defaults to 6 in both, and the λ=1 nc_cap merge is verified to fire, so **both cluster
to Nc=6 → 125 exchange rounds/sync is identical.** The only thing that changes is the sync *length*
(T_s ∝ λ), so exchange rounds **per sim-second** falls exactly 1/0.63 = **1.587×** and nothing else
about the run changes. That isolates the wiggle lever perfectly.

**A0 answers two questions:**
1. *(sanity)* rounds/sim-s falls 1.587× — near-certain by construction; confirmed via the `[lts]`
   banner (Nc) and sync count in each log.
2. **(the real test)** does exposed `MPI_Waitall` fall with the round rate? This is the assumption
   the **entire communication track** rests on, and **no run has ever varied the round rate** — the
   comm plan's payoff was fitted to a hypothetical. The job prints `R_wait = wait(0.63)/wait(1.0)`:
   - `≈1.587×` → wait scales with rounds → the A1–A3 merge/overlap levers are real → **Track A GO**.
   - `≈1.0×` → wait is latency/skew-bound → A1–A3 are worth far less → revisit Track A ordering.

   Either way we learn it in **one job** instead of after weeks of engineering.

λ=1 is also SeisSol's own configuration on this benchmark, so wiggle-off is not exotic.

## Scope — what A0 does NOT measure

A0 measures the **bulk seam-exchange wire scaling**, a steady-state property, so a short
**pre-rupture** window (`tfinal=0.6 s`; nucleation is at t₀=1.0 s) suffices and is cheap. A0 does
**not** measure fault/friction skew — that needs the rupture window (plan known-unknowns 4 & 5). The
earlier 24 % skew proxy was taken pre-nucleation and is **not** a bound on fault skew; A0 does not
claim to fix that.

## Cleanliness

- **ParaView OFF** for both legs — the B0 measurement was contaminated by ParaView + `perf`; A0
  removes that so `MPI_Waitall` is clean.
- **Caliper ON** (required — it is A0's instrument; the job aborts if the binary lacks libcaliper)
  and identical in both legs, so it cancels in the ratio.
- A0 is a **ratio** test: anything identical in both legs cancels.

## Prerequisites

- Binary built **with Caliper** and carrying commit `884349a` (the LTS `seas::spatial_dyn::step`
  region — without it `MPI_Waitall` reports as % of program, not % of step):
  `USE_CALIPER=YES bash <root>/build_expanse.sh`
- The 200 m TPV104 mesh staged on Expanse (`tpv104/mesh/tpv104_200m.msh`, ~126 MB, gitignored).

## Submit

```bash
cd <root>/miniapps/seas
sbatch jobs/comm_dev/a0_wiggle_tpv104_200m_expanse/run_a0_wiggle_expanse.sbatch
```

Cost: two ~0.6 s-window LTS legs (~45 min each) + `--lts-report` + init ≈ under 2 h; `-t 4:00:00`.

## Read the result

The job's stdout ends with an `=== [a0] RESULT ===` block: the `[lts] ENABLED` banner (λ, Nc) per
leg, the **raw** `MPI_Waitall` rows from each Caliper report (read these — the auto-ratio is
best-effort over a fixed-width format), and `R_wait` with a best-effort verdict. **Confirm `R_wait`
by hand from the raw rows before acting on it.** Full reports:
`runs/<jobid>/run_lambda063_baseline.cali-region-report.txt` and `…lambda100_wiggleoff…`.
