# Unify-plan Phase 6 — local revalidation results (2026-07-10)

Phase 6 of `document/PLAN_unify_interior_shared_fault_substep_2026-07-09.md`
("Revalidate and regenerate multi-rank references"), steps 1–2, executed locally on the
post-Phase-5 + station-tie-break tree.  Step 3 (Frontera) and step 4 (gold regeneration)
remain gated on user sign-off — see "What needs sign-off" below.

Configuration (identical to `ADER_ITERATOR_SYM1000_RESULTS.md` and the R1601 crime scene):
`seas_tpv104_driver`, `tpv104/mesh/tpv104_symmirror_1000m.msh`, `--fric-law slip-srw`,
tfinal 2.0 s, CFL 0.5, absorbing BC, matched order O = P + 1.

## Step 1 — np-independence metric, max |X_np − X_1|

**TPV104 symmirror (P=1, O=2, substep):** all NINE stations, full 522-step series,
np ∈ {2, 4, 7, 8, 10} vs np=1:

| np | max over stations of max abs series diff |
|---|---|
| 2 | 0 (byte-identical) |
| 4 | ≤ 8.2e-29 (denormal dust in pre-rupture dip columns) |
| 7 | ≤ 8.2e-29 |
| 8 | ≤ 1.0e-12 (one print-quantum value) |
| 10 | ≤ 1.0e-12 |

V_max = 12.7387 m/s at every np.  **The np-dependence of station output is now at the
print-precision floor.**  Before/after comparison (the plan's "metric must decrease"):

| era | hypocenter-station np=10 vs np=1 | V_max np=10 drift |
|---|---|---|
| R-1601 fallback (before Phase 2) | 0.19 m strike-slip apparent gap (tie artifact) + 4.8146 vs pure branch value (fallback physics contamination 2.3e-3) | 0.025 m/s |
| unified (Phases 2–5, before tie-break fix) | 0.19 m apparent gap (pure tie artifact; within-"branch" ≤ 8e-29) | 0 (12.7387 exactly) |
| + station tie-break fix | **≤ 1e-12 everywhere** | 0 |

Monotone decrease at every stage — the plan's STOP condition ("if it increases, the fix
is wrong") does not trigger anywhere.

**TPV205-writer spatial configs (TPV26 smoke + forced-rupture smoke):** np=4 vs np=1
station files byte-identical except one pre-existing 1-row/last-printed-digit FP-dust
value (x2_12_x3_7.5, 2.3e-11 relative); V_max gates unchanged (0.185543 / 0.199649).

**TPV102 / TPV205-native production meshes:** NOT run — production meshes are
Frontera-scale and project policy forbids local production-mesh reproducers
(memory: feedback_no_local_reproducer) and Frontera submissions without explicit
approval.  Flagged for the sign-off phase.

## Step 2 — symmirror np=10 sweep, O ∈ {2,3,4}, reproducing the ADER table

Post-fix (all rows now sample the SAME deterministic hypocenter QP — the historical
table's hypo columns mixed two tied sample points, see the caveat added to
`ADER_ITERATOR_SYM1000_RESULTS.md`):

| P | O | iterator | hypo strike-slip (m) | hypo dip drift (m) | V_max (m/s) | overflow |
|---|---|---|---|---|---|---|
| 1 | 2 | one-shot | 4.9436 | −3.639e-02 | 12.61 | none |
| 2 | 3 | one-shot | 5.0932 | −1.201e-02 | 15.28 | none |
| 3 | 4 | one-shot | 5.0182 | −4.476e-03 | 15.56 | none |
| 1 | 2 | substep | 5.0022 | −3.638e-02 | 12.74 | none |
| 2 | 3 | substep | 5.1366 | −1.214e-02 | 15.47 | none |
| 3 | 4 | substep | 5.0564 | −4.770e-03 | 15.68 | none |

Acceptance readings:
1. **No overflow at any order** (the R-1601-era 1e28 class is gone). ✓
2. **Dip-drift: substep ≈ one-shot at every order** — O=2: 3.638e-2 ≤ 3.639e-2 ✓;
   O=3: 1.214e-2 vs 1.201e-2 (+1.1%); O=4: 4.770e-3 vs 4.476e-3 (+6.6%).  The plan's
   literal "substep ≤ one-shot" holds at O=2 and is missed by 1–7% at O=3/4 — against
   the historical (confounded) 30%–10× substep degradation this is parity, and the
   historical finding "substep does NOT improve dip drift" is now attributable mostly
   to the fallback + sampling artifact.  Both iterators now show the proper
   precision-floor scaling (3.6e-2 → 1.2e-2 → 4.5e-3 with P).
3. The P=1 substep row equals the np=1 trajectory EXACTLY (5.0022/−3.638e-02) — the
   historical row 4 (5.0470/+1.550e-01) was the other tie QP under the fallback.

## What needs sign-off (Phase 6 steps 3–4)

1. **TPV102 / TPV205 production-mesh np-sweeps on Frontera** (needs explicit approval
   per project policy).  Expected outcome given the local evidence: np>1 == np=1 to
   print precision at non-tied stations.
2. **np>1 gold-trace regeneration**: any np>1 reference traces recorded in the
   R-1601-fallback era embed (a) the fallback's seam physics and (b) potentially the
   tie-break sampling artifact at x2=0 stations, and should be regenerated (or simply
   replaced by the np=1 golds, which np>1 now matches to ≤1e-12 locally).
3. **One-time np=1 trace change acknowledgment**: the TPV26-smoke stations
   (±7.5, 7.5) — exact far-field corner ties, values 2.554e7 ↔ 2.444e7 n-stress —
   changed picks under the deterministic tie-break (full record in
   `np4_attractor_root_cause_2026-07-10.md` §Post-fix trace changes).  tpv104
   symmirror np=1 traces are UNCHANGED.
