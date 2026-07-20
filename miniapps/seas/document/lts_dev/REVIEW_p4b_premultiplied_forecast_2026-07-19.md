# P4b — flux-PREMULTIPLIED (forecast) exchange: impl + review — 2026-07-19

Branch `safs-v4_0_0-alt-case1-mfem-speed` (local, unpushed).  This increment
replaces the raw-D(k) diff-1 ghost exchange of the batched 4b
(REVIEW_p4b_batched_exchange_2026-07-19.md) with the **premultiplied forecast**
exchange — the deeper D-7 payload — in its **byte-identical** form.

## What changed (the EDGE innovation, byte-identical form)
Batched-4b (and 4a) exchanged the coarse provider's **raw Taylor stack D(k)** at
each PREDICT of a cluster c>=1; the finer neighbour then integrated that ghost
D(k) over its closed-form sub-interval `[a,b]` at its own correct.  This increment
moves the integration to the side that OWNS the D(k):

- **Coarse side** integrates its OWN retained D(k) into the forecast block
  `IntegrateTaylor(a,b,D(k))` and exchanges the **already-projected forecast**.
- **Fine side** (diff-1 seam mode 1) reads the exchanged forecast **directly**
  from `FaceNbrData` via the byNODES vdof map — no local integration.

The exchange therefore moves from **PREDICT → CORRECT**: Predict is now
retain-only (`PrepareSeamCoarseForecast` does the extra element-local CK and
records the epoch, but does NOT exchange); the forecast is integrated + exchanged
once, in `ComputeADERClusterSeamFaceFluxRHS`'s new pre-pass, at the fine cluster's
Correct — one batched (all-NUM_STATE) collective per correcting cluster that has a
coarser neighbour (`c < num_clusters-1`).

### Why it is byte-identical to batched-4b
`IntegrateTaylor(a,b,D(k))` is deterministic; both of its inputs are bit-identical
on the two ranks:
- **D(k)**: a bit-copy — the coarse RETAINS the same raw stack the fine used to
  receive (dt-independent element-local CK); epoch-guarded against staleness.
- **[a,b]**: closed-form from the rank-uniform schedule (`t_s, dt_base, tick,
  T_actual` + the coarse cluster id `c+1`), so identical on both ranks and for
  every coarse-seam element at a given Correct(c) — computed ONCE in the pre-pass.

Hence the coarse-computed forecast equals, bit-for-bit, the forecast the fine used
to compute locally, and `fc_all[vdof(c,i)]` (new) == `forecast_blk[c*ndof+i]`
(old).  Mode 2 (coarse side of a diff-1 seam) still integrates its own D(k) with
the identical `[a,b]`, so both sides consume the identical coarse forecast ⇒
conservation exact.

## Files touched
- **dynamic/wave_operator.inl** — `ComputeADERClusterSeamFaceFluxRHS`: new
  `exchange_forecast` param + the forecast **pre-pass** (integrate local
  seam-coarse D(k) → pack byNODES → 1 batched exchange → DEEP-COPY `fc_all` before
  the I-exchange reuses the ghost GF).  Mode 1 reads `fc_all[vdof]`; mode 2
  unchanged (local integrate).  `PrepareSeamCoarseForecast`: retain-only + record
  `seam_coarse_dk_epoch_`.  `AdvanceADERClusterBulk`: forwards `exchange_forecast`.
  Deleted `ExchangeClusterProviderDkGhost` (the raw-D(k) primitive).
- **dynamic/wave_operator.hpp** — decls for the new param; removed `ghost_dk_` /
  `ghost_dk_epoch_` members + `GhostDkForTest` accessor + the
  `ExchangeClusterProviderDkGhost` decl; added `seam_coarse_dk_epoch_`.
- **dynamic/lts_stepper.hpp** — `BuildTickTable`: the diff-1 term is now
  `+1 per correcting cluster c < num_clusters-1` (was `+ader_order per predicting
  c>=1`).  `LtsGlobalMeta.exchange_bulk_provider_dk` docstring updated.
- **dynamic/lts_bulk_stepper.hpp / lts_fault_stepper.hpp** — Correct forwards
  `exchange_forecast = exchange_dk_ && (c < num_clusters-1)`.
- **tests/unit/test_lts_mpi_seam.cpp** — removed the raw-D(k) round-trip check
  (the primitive it exercised no longer exists; the forecast path is covered
  end-to-end by `test_lts_mpi_seam_diff1`).

## Byte-vs-4a / vs-batched-4b verification (GREEN — Δ identical)
| Gate | result | == committed 4b? |
|---|---|---|
| diff-0 seam np=2 | Δenergy=4.441e-16, Δcen=0.0, Δmax=0.0 (58/58) | YES (58 vs 60: −2 = removed round-trip CHECK ran on both ranks) |
| diff-0 seam np=4 | 4.441e-16 / 4.388e-16 / 0.0 (116/116) | YES (116 vs 120: −4 = removed CHECK ran on 4 ranks) |
| diff-1 Nc=3 np=3 | 1.614e-16 / 1.780e-16 / 0.0 (111/111) | YES (identical Δ) |
| fault interleave np=1-vs-np=2 (TPV104-1000m, Nc=6, 13 syncs) | worst finite max\|Δfield\| = 6.907e-26 (tol 1.2e-4); 9/9 SCEC stations byte-identical; V_max=1.21545e-14 both | round-off floor (== committed-4b order 1.666e-27) |
| np=10 fault interleave vs np=1 (no-hang; fault-FREE-rank path) | worst max\|Δfield\|=1.333e-26; **cycle-1=1.666e-27 (== committed-4b's reported value)**; 9/9 stations byte-identical; V_max=1.00358e-16 both; clean completion, no hang/abort/matched-collective violation | YES |
| regression | predictor 31 · reorder 14 · scheduler 31 · layout 102 · clustering 4043 · partition 19 · tpv102-ADER 4 | all match |

**All gates GREEN.**  Every seam Δ reproduces committed 4b bit-for-bit; the only
count changes are the removed raw-D(k) round-trip CHECK (−nranks per np).  The
fault np=1-vs-np{2,10} interleave differentials sit at the round-off floor
(≤7e-26), the np=10 run exercised the fault-free-rank path + matched-collective
contract cleanly, and cycle-1 at np=10 reproduces committed-4b's exact 1.666e-27.

## Adversarial review — self (byte-identity argument + gates)
- The forecast deep-copy (`fc_all`) is taken BEFORE the I-exchange reuses
  `q_gf_full` — no aliasing.  `nbr_all` (the I) is aliased once and read across the
  whole face loop, never re-exchanged inside it.
- The forecast exchange is rank-uniform (fires per Correct(c) from the shared tick
  table iff `exchange_dk_ && c<num_clusters-1`); a rank with no local seam-coarse
  elems still calls `ExchangeFaceNbrData` (matched collectives, P-007).  The tick
  table's diff-1 term was updated to match (correct-side +1, predict-side dropped).
- GAP-A3 staleness: the pre-pass asserts `seam_coarse_dk_epoch_[c+1] ==
  (tick/period)*period` — the retained D(k) is from the coarse cluster's current
  step.  Mode-1 asserts `have_forecast` (fail-loud if the stepper wiring omits the
  exchange for a Correct that has a mode-1 seam).

## Independent second-pass review (2026-07-19) — CONFIRMED
A second adversarial pass re-read the changeset from scratch, hunting for: `fc_all`
aliasing vs the I-exchange repack, matched-collective divergence for the forecast
term at scale, `[a,b]` inconsistency between pre-pass and mode-1/mode-2, retain-only
correctness of `PrepareSeamCoarseForecast`, a dangling raw-D(k) reference, and the
`n_collectives` term when Nc=1.  Every invariant above was re-verified independently
(deep-copy ordering at inl:2706 precedes the I-exchange repack at inl:2715+; the
forecast `ExchangeFaceNbrData` at inl:2703 is outside the pack gate ⇒ unconditional
under the rank-uniform `exchange_forecast`; pre-pass `[a,b]` at inl:2668-2674 is the
identical closed form as mode-2 at inl:2828-2834; `ghost_dk_`/
`ExchangeClusterProviderDkGhost`/`GhostDkForTest` fully removed — only a comment
mentions the old name; both steppers gate identically).  Two LOW perf items noted
(mode-2 re-integrates the same local D(k) the pre-pass already integrated; mode-2
recomputes the constant `[a,b]` per face) — NOT applied (byte-verified path, not
worth the churn/risk).  Re-ran the end-to-end diff-1 gate on the on-disk state:
np=1 33/33, np=3 111/111, `d(energy)=1.614e-16 d(cen)=1.780e-16 d(max)=0.000e+00`
(byte-identical).  **Verdict: PASS — no fixes required.**
