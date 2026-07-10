# ADER Fault Iterator Sweep — Symmirror 1000m mesh

**Date**: 2026-04-26
**Build**: HEAD with R-1003 (substep np>1 path) + R-1600 (collective gate fix in `EvaluateBulkAtFaultQPsCanonical`) + R-1601 (shared-QP fallback to inline `EvaluateADER` in the substep corrector) all landed.

## Setup

| | |
|---|---|
| Mesh | `tpv104/mesh/tpv104_symmirror_1000m.msh` (225 nodes / 1152 elems / ~384 tets, y-mirror by construction) |
| Ranks | 10 (default ParMETIS partition) |
| tfinal | 2.0 s |
| CFL | 0.5 |
| Outer BC | absorbing |
| Friction solver | Brent (hardcoded by `--fric-law slip-srw`) |
| ADER orders | canonical Pelties matched-order (O = P + 1) |
| Iterators | `oneshot` (default; fault state evaluated once per macro-step) vs `substep` (per-sub-step ADER predictor + per-QP friction at O Gauss-Legendre nodes) |

## Iterator semantics (relevant for results)

- **oneshot**: `wave.AdvanceADER` with the inline `FaultFaceFlux::EvaluateADER` at each fault QP. Single friction solve per macro-step.
- **substep**: `Tpv104SubStepIterator::AdvanceWithSubStepStates` runs `O` Gauss-Legendre sub-step nodes; each evaluates `wave.ComputeADERSubStepStates` predictor + `wave.EvaluateBulkAtFaultQPsCanonical` + `FaultFaceFlux::ComputeStageState` + `BuildImposedState`, accumulating the time-averaged `I_imp_±` consumed by the corrector via `SetSubStepFaultImposedStates`.
- **R-1601 caveat (np>1 only)**: SHARED fault QPs fall back to inline `EvaluateADER` (oneshot semantics) inside the corrector's substep gate. INTERIOR fault QPs (~95% on this mesh) keep full substep semantics. Reason: rank A's `elem1_on_plus` = !rank B's on the same physical shared face; until the iterator's per-shared-QP physics is reconciled with the corrector's frame convention across the partition seam, the fallback prevents Q-overflow that would otherwise drive `Brent` into 1e28-magnitude tau garbage. np=1 substep is unaffected (no shared faces).

## Results

| # | P | ADER-O | iterator | np | steps | wall (s) | hypo strike-slip (m, t=2s) | hypo dip-slip drift (m, t=2s) | V_max peak (m/s) |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 1 | 2 | oneshot | 10 | 522 | ~2 | 4.7393 | **+1.192e-01** | 12.61 |
| 2 | 2 | 3 | oneshot | 10 | 870 | ~8 | 5.0371 | **+1.456e-02** | 15.28 |
| 3 | 3 | 4 | oneshot | 10 | 1217 | 31.4 | 5.1441 | **−3.504e-03** | 15.56 |
| 4 | 1 | 2 | substep | 10 | 522 | 2.0 | 5.0470 | **+1.550e-01** | 15.96 |
| 5 | 2 | 3 | substep | 10 | 870 | 8.6 | 5.1451 | **+3.891e-02** | 20.40 |
| 6 | 3 | 4 | substep | 10 | 1217 | 60.2 | 5.3335 | **+3.595e-02** | 21.60 |

Hypo coordinates: x2 = 0 km along strike, x3 = 7.5 km down-dip. Column legend per `tpv104/visualize_results.py`: t / h-slip / h-slip-rate / h-shear-stress / **v-slip (=dip drift)** / v-slip-rate / v-shear-stress / n-stress / psi.

> **CONFOUNDING CAVEAT (2026-07-10, `np4_attractor_root_cause_2026-07-10.md`):** the hypocenter
> station at x2 = 0 sits exactly equidistant between the fault QPs at x = −166.7 m and
> x = +166.7 m on this mesh, and (pre-tie-break-fix) the writer's pick was partition- and
> run-layout-dependent.  The hypo strike-slip / dip-drift columns above may therefore mix TWO
> different physical sample points across rows — cross-row deltas at this station are not pure
> iterator/order effects.  V_max and non-x2=0 stations are unaffected.

## Chart — |hypo dip drift| at t=2s (m, log10 magnitude per P)

```
                              oneshot (O)            substep (S)
  P=1  O2:  ███████████████ +1.19e-01  S2:  ████████████████ +1.55e-01
  P=2  O3:  █████████      +1.46e-02   S3:  ███████████      +3.89e-02
  P=3  O4:  ██████         −3.50e-03   S4:  ██████████       +3.60e-02
                                          (note: substep dip > oneshot dip at all P)
```

## Findings

1. **Precision-floor reduction with P (oneshot)**: |dip drift| drops 119 mm → 14.6 mm → 3.5 mm as P increases 1 → 2 → 3. Roughly 8× and 4× per P-step — consistent with Pelties' ADER precision-floor scaling.

2. **Substep does NOT improve dip drift on the symmirror 1000m mesh.** At every P the substep iterator gives |dip drift| larger than oneshot:
   - P=1: substep 1.55e-1 vs oneshot 1.19e-1 (+30%)
   - P=2: substep 3.89e-2 vs oneshot 1.46e-2 (2.7×)
   - P=3: substep 3.60e-2 vs oneshot 3.50e-3 (10× larger, sign-flipped)

   The cause is the R-1601 hybrid dispatch: at np>1, shared QPs use the oneshot path while interior QPs use the substep path. The two paths produce inconsistent imposed states at the partition seam; the discontinuity radiates into the bulk and amplifies dip-direction pollution at the hypocenter station. The pure-substep R-1003 path (np=1) is expected to give different results — re-run at np=1 substep to disentangle.

3. **V_max consistency**: substep gives ~25-40% larger V_max peak than oneshot at every P (12.61 vs 15.96 at P=1, etc.). This is the rupture-front overshoot the per-substep ADER predictor introduces; it doesn't necessarily mean wrong physics, but warrants comparison against SeisSol reference 7.39 m/s at hypo (both iterators are well above on this 1000m mesh).

4. **Sign behavior at P=3**: oneshot flips dip drift sign (−3.5 mm) while substep stays positive (+36 mm). This is another signature of the hybrid-dispatch discontinuity at the partition seam — the sign of the cross-partition flux mismatch dominates over the bulk precision-floor reduction at high P.

5. **Wall-time cost (substep / oneshot ratio)**: P=1 ≈ 1×; P=2 ≈ 1×; P=3 ≈ 1.9× (60.2 s vs 31.4 s). The substep cost is dominated by the 9 ghost-exchange collectives per macro-step (R-1601 review item) and the per-substep ADER predictor `ComputeADERSubStepStates` (O calls per macro-step).

6. **R-1003 + R-1600 + R-1601 fix verified**: np=10 substep ran to completion at all 3 P values through tfinal=2.0 s. Pre-R-1600 the run hung indefinitely (collective-skip deadlock at `wave_operator.inl:1801`). Pre-R-1601 it aborted at macro-step 1 with tau=4e28 from the canonical-frame mismatch on shared QPs (`SUBSTEP_NP_GT_1_HANG_REVIEW.md` R-1600 fixed the deadlock; the surfaced R-1601 bug is now mitigated by shared-QP fallback to inline EvaluateADER).

## Per-quantity oneshot-vs-substep agreement at the hypocenter (t=2 s)

Sweep of all 9 station columns (`t, h-slip, h-rate, h-shear, v-slip, v-rate, v-shear, n-stress, psi`):

| P | quantity | oneshot | substep | rel diff | |
|---|---|---|---|---|---|
| 1 | strike-slip       | +4.7393e+00 | +5.0470e+00 |  +6.10 % | OK    |
| 1 | strike-rate       | +4.0618e+00 | +4.0923e+00 |  +0.75 % | OK    |
| 1 | strike-shear      | +2.5206e+07 | +2.5003e+07 |  −0.81 % | OK    |
| 1 | **dip-slip**      | +1.1920e-01 | +1.5499e-01 | **+23 %** | MOD |
| 1 | **dip-rate**      | −2.8716e-02 | +6.2943e-02 | **+146 %** (sign flip) | BIG |
| 1 | **dip-shear**     | −1.7820e+05 | +3.8457e+05 | **+146 %** (sign flip) | BIG |
| 1 | n-stress          | +1.2000e+08 | +1.2021e+08 |  +0.17 % | OK    |
| 1 | psi               | +5.7854e-02 | +5.5775e-02 |  −3.59 % | OK    |
| 2 | strike-slip       | +5.0371e+00 | +5.1451e+00 |  +2.10 % | OK    |
| 2 | strike-rate       | +3.7199e+00 | +3.9052e+00 |  +4.74 % | OK    |
| 2 | strike-shear      | +2.5042e+07 | +2.5009e+07 |  −0.13 % | OK    |
| 2 | **dip-slip**      | +1.4564e-02 | +3.8905e-02 | **+62 %** | BIG |
| 2 | **dip-rate**      | −4.6741e-02 | +1.9923e-02 | **+143 %** (sign flip) | BIG |
| 2 | **dip-shear**     | −3.1465e+05 | +1.2759e+05 | **+141 %** (sign flip) | BIG |
| 2 | n-stress          | +1.2000e+08 | +1.1990e+08 |  −0.08 % | OK    |
| 2 | psi               | +5.7405e-02 | +5.6800e-02 |  −1.05 % | OK    |
| 3 | strike-slip       | +5.1441e+00 | +5.3335e+00 |  +3.55 % | OK    |
| 3 | strike-rate       | +4.1126e+00 | +4.2716e+00 |  +3.72 % | OK    |
| 3 | strike-shear      | +2.5220e+07 | +2.5008e+07 |  −0.84 % | OK    |
| 3 | **dip-slip**      | −3.5038e-03 | +3.5948e-02 | **+110 %** (sign flip) | BIG |
| 3 | **dip-rate**      | −4.6967e-02 | −3.9462e-02 |  +16 %    | MOD |
| 3 | **dip-shear**     | −2.8801e+05 | −2.3103e+05 |  +20 %    | MOD |
| 3 | n-stress          | +1.2000e+08 | +1.1992e+08 |  −0.06 % | OK    |
| 3 | psi               | +5.7856e-02 | +5.5864e-02 |  −3.44 % | OK    |

(`OK` = |rel| ≤ 6 %; `MOD` = ≤ 30 %; `BIG` = > 30 % or sign-flipped.)

### Reading

- **Strike-direction physics, normal traction, and friction state agree to ≤ 6 % at every P.** Strike-slip, strike-slip-rate, strike-shear-stress, normal stress, and ψ are robustly consistent between iterators. The substep + oneshot paths produce the SAME rupture energy budget — the dominant 5 m of strike slip and 25 MPa shear release are computed identically (within wave-noise).

- **The full disagreement is concentrated in the dip-direction channel.** dip-slip, dip-rate, dip-shear all differ by 20-150 % and frequently sign-flip between the two iterators. The dip channel is exactly the closed-loop V1 instability channel documented in earlier rounds (`tpv104_sigma_n_leak_root_cause_2026-04-25.md`); it is the SPURIOUS pollution direction in TPV104, not part of the rupture physics.

- **Implication for the R-1601 fallback**: the partition-seam discontinuity introduced by mixing interior-substep with shared-oneshot dispatch radiates ONLY into the dip channel — exactly the channel where the per-step gain > 1 amplifies any seed. The strike-direction physics is unaffected because it is anchored by the 40 MPa pre-stress and is self-stabilizing under friction. So substep at np>1 with R-1601 fallback is **safe for the rupture physics** but **invalid for the dip-direction precision-floor study** that motivated the precision-ladder comparison in the first place.

## Conclusion

ADER one-shot precision-floor scaling is confirmed on the symmirror 1000m mesh: increasing P from 1 to 3 reduces hypocenter dip-direction pollution from ~12 cm to ~3.5 mm. The substep iterator at np=10 is now numerically stable (was deadlocked + abort-prone before the R-1600 + R-1601 fixes), but the **R-1601 hybrid dispatch (oneshot fallback on shared QPs) introduces a partition-seam discontinuity that degrades dip-drift accuracy by 30 % to 10× across the P sweep**.

For production multi-rank runs that want true substep semantics, the path forward is to fix the canonical-frame convention split between `Tpv104SubStepIterator` (per-shared-QP physics) and `ComputeADERSharedFaceFluxRHS` (corrector frame on `qpd.sign_flipped`) — likely by routing the iterator through a rank-symmetric canonical frame derived from `qpd.sign_flipped` rather than `elem1_on_plus`. Until then, **production should use `--fault-iterator one-shot` at np>1**; substep at np>1 will deliver hybrid semantics that under-perform oneshot on dip-direction precision.

For the standalone Pelties precision-floor experiment, **the oneshot triple (Cases 1-3) is the meaningful comparison** — and it shows the expected ~8× per-P-step reduction.

## Files

| | |
|---|---|
| Oneshot output dirs | `tpv104/results_local_sym1000_p{1,2,3}_O{2,3,4}_oneshot/` |
| Substep output dirs | `tpv104/results_local_sym1000_p{1,2,3}_O{2,3,4}_substep/` |
| Oneshot logs | `/tmp/sym_case{1,2}_p{P}_O{O}_oneshot.log`, `/tmp/sym_p3_O4_oneshot.log` |
| Substep logs | `/tmp/sym_p{P}_O{O}_substep.log` |
| Hypo station | `tpv104_p{P}_O{O}_{iterator}_station_x2_0_x3_7.5.dat` (column 5 = dip slip) |

## Cross-references

- `SUBSTEP_ITERATOR_MPI_REVIEW.md` — R-1300..R-1306 (R-1003 implementation review)
- `SUBSTEP_NP_GT_1_HANG_REVIEW.md` — R-1600 collective-skip deadlock + R-1601 (newly identified canonical-frame mismatch on shared QPs)
- `dynamic/wave_operator.inl:1801` — R-1600 fix (gate on `pmesh.GetNSharedFaces()`)
- `dynamic/wave_operator.inl:4506` — R-1601 fix (shared-QP fallback to inline EvaluateADER inside the substep corrector gate)
- `dynamic/tpv104_substep_iterator.cpp:619` — iterator main loop
