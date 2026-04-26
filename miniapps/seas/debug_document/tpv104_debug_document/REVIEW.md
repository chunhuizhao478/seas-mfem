# Code Review: TPV104 — Round 9: SeisSol reference comparison settles R-803/R-805 (2026-04-25)

## Bottom line up front

I read all 9 SeisSol benchmark trace files in `tpv104/benchmark_data/seisol/`. **The agent's "SeisSol claims slip_dip ≡ 0 to 8 sig figs" framing is wrong.** Actual SeisSol numbers:

- At the **hypocenter station** (0, 0, -7.5 km), within t ≤ 2 s: peak |v-slip| = **50 µm** at t = 1.97 s (`tpv104_seisol_x2_0_x3_7.5.txt`). Peak |v-slip-rate| = 2.86e-3 m/s.
- At **off-axis stations** (±9 km along strike, hypo depth), within t ≤ 12 s: peak |v-slip| = **7.7 cm**.
- At **corner stations** (±12 km along strike, 3 km depth), within t ≤ 12 s: peak |v-slip| = **41 cm**.
- All non-hypocenter stations have ZERO h-rate within t ≤ 2 s — the rupture front hasn't propagated to them yet by t = 2 s.

**MFEM P=2 hypocenter slip_dip = 58 µm at t = 2 s is essentially matching SeisSol's hypocenter slip_dip (50 µm peak).** The 174× reduction got us *to* SeisSol's numerical accuracy on this metric, not a tiny fraction of it. **The investigation can stop here on the hypocenter slip_dip metric.**

Two concerns remain (R-901, R-902 below): the V_max metric is comparing different things, and we have no MFEM data at off-axis stations to compare against SeisSol's *real* (non-zero) slip_dip there.

---

## Audit of the agent's interpretation

The R-805 / R-803 narrative — "MFEM P=2 has 60 µm at hypocenter, SeisSol claims 8 sig figs of zero, gap is 3 orders of magnitude" — was built on a misread. I verified the actual SeisSol traces:

```
$ awk 'NR>22 && NF>=6 && $1<=2.001 {abs5=($5<0)?-$5:$5; if(abs5>m){m=abs5; t=$1}}
       END{print "peak |v-slip| within t≤2s =", m, "m at t=", t, "s"}' \
   tpv104_seisol_x2_0_x3_7.5.txt
peak |v-slip| within t≤2s = 5.0438e-05 m at t = 1.9697 s
```

**SeisSol's peak |v-slip| at the hypocenter within the first 2 seconds is 50 µm.** The 8-sig-figs claim was about the FINAL-time slip_strike value, or a different comparison metric, not about slip_dip being identically zero. Re-checking the prior debug docs:

> `tpv104_bulk_asymmetry_followup_2026-04-25_pm.md` §13: "[12 s] hypocenter ... slip_dip │ −7 mm ← spurious │ ≈ 0"

The document said "≈ 0" as a baseline expectation, and the agent operationalized that as "exactly zero to FP precision." The actual SeisSol numbers are far less stringent.

**Side-by-side at the hypocenter, t ≤ 2 s:**

| metric | SeisSol (O5 ADER-DG) | MFEM P=1 | MFEM P=2 |
|---|---|---|---|
| peak \|v-slip\| (slip_dip) | 50 µm @ t=1.97 s | 1.0e-2 m = 10 mm @ t=2 s | 58 µm @ t=2 s |
| peak \|v-slip-rate\| (V_dip) | 2.9e-3 m/s @ t=1.96 s | 2.6e-2 m/s | 1.7e-2 m/s |
| peak \|h-slip-rate\| (V_strike) | 7.39 m/s @ t=1.07 s | (not reported per-station) | (not reported per-station) |

**MFEM P=2 slip_dip at hypocenter ≈ SeisSol slip_dip at hypocenter to within ~20%.** Both are O(50 µm). Not 4 orders apart; not 1 order apart; same order of magnitude.

MFEM P=2's V_dip rate (17 mm/s) is ~6× SeisSol's (2.9 mm/s) — the rupture transit is more abrupt in MFEM, but integrates to roughly the same total dip slip. That's a separate, lower-priority observation.

**MFEM P=1 was the outlier** — 10 mm slip_dip is 200× SeisSol's level. P=1 was simply under-resolved.

---

## What this resolves

| round-8 finding | round-9 status |
|---|---|
| R-803: V_max gap unexplained — convergence vs precision story | partial — see R-901 below |
| R-805: is 60 µm acceptable for SCEC verification? | **resolved.** SeisSol's own numerical accuracy at hypocenter is 50 µm. MFEM P=2's 58 µm is at SeisSol's level. Acceptable. |
| R-806: P=2 production-scale validation needed | unchanged. Still needed before shipping P=2 default. |
| R-801: H2 V1=0 lock retraction | unchanged. Still must come out of any candidate set; the actual SeisSol slip_dip is **not** zero, so an exact-zero clamp is provably wrong. |

The categorically-rejected H2 lock proposal was a fix for a problem that doesn't exist as stated. SeisSol's slip_dip is not zero, and any TPV104 implementation that forces V1 = 0 would be **less accurate** than the current MFEM P=2 result (which oscillates V1 around the SeisSol envelope, but at the right magnitude). H2 is not just rejected because it's masking — it's actually wrong physics.

---

## What this does NOT resolve

### R-901 — V_max metric mismatch

MFEM reports `V_max = max over all fault DOFs over all time` (driver line 1505+). SeisSol's per-station files give peak slip rate at SCEC's specific (x, z) sample points. These are different metrics.

Within t ≤ 2 s, SeisSol's peak |h-slip-rate| over all 9 sampled stations is 7.39 m/s (at the hypocenter — the only station the rupture has reached by t=2 s). MFEM P=1 reports V_max = 15.18 m/s, P=2 reports 17.71 m/s — both 2–2.4× SeisSol's hypocenter peak.

This could be:
- MFEM's V_max picks up an edge / corner DOF that SeisSol's stations don't sample. Most likely.
- MFEM has a real V_strike overshoot at the rupture front, sampled only by the V_max statistic.
- Different mesh resolution: SeisSol uses adaptive curved meshes at O5 polynomial order; MFEM is 500 m unstructured tets at P=1/P=2.

To audit: print MFEM's slip_rate AT the hypocenter QP (via the existing station writer) instead of relying on V_max. Compare against SeisSol's 7.39 m/s. If MFEM hypocenter peak h-rate is also ~7.4 m/s, V_max is just sampling an edge artifact and the rupture itself is fine. If MFEM hypocenter peak h-rate is ~15-17 m/s, the rupture is genuinely overshooting SeisSol.

### R-902 — Off-axis stations untested

SeisSol shows large, genuine, physically-real slip_dip at off-axis and corner stations:

| SeisSol station | peak |slip_dip| (full 12s run) | physical interpretation |
|---|---|---|
| (±9, 0, -7.5) — same depth as hypo, off-strike | 7.7 cm | rupture front lobe |
| (±12, 0, -3) — corner, top of fault | **41 cm** | strong dip-direction motion at the fault edge |
| (±12, 0, -12) — corner, bottom | 5.8 cm | edge effect at deep boundary |
| (0, 0, -3) — top center | 0.55 mm | small but non-zero |
| (0, 0, -12) — bottom center | 0.30 mm | small but non-zero |

**MFEM has not been compared at any off-axis station.** A "fix" that drives MFEM's hypocenter slip_dip to zero would be a regression at every other station where SeisSol shows non-zero values. This makes the H2 lock proposal even more wrong than I previously characterized — it would WORSEN agreement with SeisSol at the corner stations, not improve it.

For SCEC verification compliance, the hypocenter (0, 0, -7.5) is one station among many. **The full SCEC TPV104 cross-comparison evaluates all 9 stations.** If MFEM matches at hypocenter but fails at (12, 0, -3), the submission won't pass.

The only way to know is to compare. The MFEM driver writes per-station outputs (`tpv104_x2_*_x3_*.dat`) — those files exist in `plots_results_*` from prior runs. The comparison is a Python script away.

### R-903 — V_dip rate is 6× SeisSol's

SeisSol peak |v-slip-rate| at hypocenter = 2.86e-3 m/s. MFEM P=2 = 1.67e-2 m/s. Ratio: 5.8×.

The TIME-INTEGRAL of MFEM P=2's V_dip matches SeisSol's slip_dip (~50 µm). But the RATE is 6× higher and concentrated in a narrower time window. Means MFEM's rupture front transit at the hypocenter is more abrupt than SeisSol's.

This may or may not pass SCEC's transient comparison. Depends on which metric SCEC weights. Not a blocker for shipping P=2 default, but worth a head-to-head time-series plot (MFEM P=2 hypocenter trace overlaid on SeisSol) to assess shape similarity.

---

## Findings (round 9)

---

### [R-901] [MODERATE] [drivers/tpv104_driver.cpp::V_max statistic + benchmark comparison] — V_max comparison is metric-mismatched, audit per-station instead

**Category:** ASSUMPTION (the V_max numbers people are comparing aren't the same thing)

**Description:**
The agent's round-8 narrative ("P=1 V_max = 15.18, P=2 V_max = 17.71, 17% shift suggests rupture under-resolution") compared against SeisSol's hypocenter peak h-rate of 7.39 m/s — that's a 2× gap, and we don't know if the gap is real overshoot vs an MFEM-side edge-DOF artifact picked up by `max-over-all-fault-DOFs`.

The driver writes per-station `*.dat` files via `station_writer.WriteStep(t, dof_data)` (line 1523). The hypocenter station's per-step `slip_rate` is in those files. Comparing MFEM hypocenter slip_rate trace against SeisSol's `tpv104_seisol_x2_0_x3_7.5.txt` column 3 directly answers the metric-match question.

**Trigger:** any decision based on V_max comparison.

**Suggested fix (no code change — analysis):**
```bash
# After the P=2 run, locate the hypocenter station file:
ls miniapps/seas/tpv104/plots_results_*/tpv104_x2_0_x3_7.5.dat

# Find peak h-slip-rate (SCEC column 3 ↔ MFEM "V_strike" or |V| approximation):
awk 'NR>2 && $3>m {m=$3; t=$1} END {print "MFEM peak h-srate", m, "at t =", t}' tpv104_x2_0_x3_7.5.dat

# Compare to SeisSol:
awk 'NR>22 && $3>m {m=$3; t=$1} END {print "SeisSol peak h-srate", m, "at t =", t}' \
    benchmark_data/seisol/tpv104_seisol_x2_0_x3_7.5.txt
# (Already computed above: 7.39 m/s at t=1.07 s)
```

If MFEM hypocenter peak h-srate is ~7-8 m/s, P=1 / P=2 V_max statistics are just picking up an edge DOF and the rupture is converged. If MFEM hypocenter peak h-srate is 15-18 m/s, the rupture genuinely overshoots SeisSol and the P=2 ≠ P=1 V_max gap is a real convergence concern.

**Test case:** the awk above is the test.

---

### [R-902] [CRITICAL] [BLOCKING the "ship P=2 default" decision] — Off-axis stations not yet compared

**Category:** ASSUMPTION (untested)

**Description:**
P=2 hypocenter slip_dip ≈ SeisSol hypocenter slip_dip. That's one of nine SCEC-spec stations. The off-axis stations have orders-of-magnitude larger genuine slip_dip in SeisSol (up to 41 cm at (12, 0, -3)). MFEM's behavior at these stations under P=1 vs P=2 has not been audited in this debug-document chain.

If MFEM P=2 at (12, 0, -3) produces slip_dip very different from SeisSol's 41 cm, then the "P=2 default fixes TPV104" claim is premature. Could be:
- MFEM P=2 matches at all 9 stations → ship P=2 default, fix is real and global.
- MFEM P=2 matches only at the hypocenter; fails elsewhere → bug is more subtle than a noise-floor effect.
- MFEM P=2 over-shoots dip slip at corner stations → opposite-direction problem; need diagnostic.

The user has been running 12 s simulations earlier (job 7677822/7677831 per round-1 docs). Those runs presumably wrote per-station files. Compare MFEM peak |v-slip| at each of the 9 SCEC stations against the SeisSol numbers in this round's table.

**Trigger:** decision to ship P=2 as TPV104 default.

**Actual behavior:** unverified at off-axis stations.

**Expected behavior:** if MFEM P=2 matches SeisSol at all 9 stations to within 1 order of magnitude (or whatever SCEC tolerance specifies), then ship. If not, that's the real bug to investigate.

**Suggested fix (no code change — analysis):**
```bash
# For each of the 9 SCEC stations, peak |v-slip| comparison:
for station in -12_x3_12 -12_x3_3 -9_x3_7.5 0_x3_12 0_x3_3 0_x3_7.5 12_x3_12 12_x3_3 9_x3_7.5; do
  seisol_peak=$(awk 'NR>22 {abs5=($5<0)?-$5:$5; if(abs5>m){m=abs5}} END{print m}' \
                 benchmark_data/seisol/tpv104_seisol_x2_${station}.txt)
  mfem_peak=$(awk 'NR>2 {abs5=($5<0)?-$5:$5; if(abs5>m){m=abs5}} END{print m}' \
              plots_results_p2_run/tpv104_x2_${station}.dat)
  echo "$station: SeisSol=$seisol_peak  MFEM_P2=$mfem_peak  ratio=$(echo $mfem_peak / $seisol_peak | bc -l)"
done
```

Pass criterion: ratio in [0.5, 2.0] across all 9 stations. Fail: any ratio outside that range deserves investigation.

**Test case:** the bash above. The actual run that produced `plots_results_p2_run` would be a 12 s P=2 sbatch (probably needs Frontera unless the local 4×2×10 km repro takes 12 s of physical time).

---

### [R-903] [LOW] [v-slip-rate timing concern] — MFEM P=2 V_dip rate at hypocenter is 6× SeisSol

**Category:** EDGE_CASE

**Description:**
MFEM P=2 V_dip = 1.67e-2 m/s; SeisSol peak |v-slip-rate| at hypocenter = 2.86e-3 m/s. Ratio 5.8×. Total integrated dip slip matches (both ~50 µm), but MFEM's rupture-front transit at the hypocenter is concentrated in a narrower time window. Could cause a poor SCEC waveform-match score even though the integrated slip is correct.

Not a blocker for shipping; flag for follow-up after R-902 confirms global convergence.

**Suggested action:** overlay MFEM P=2's hypocenter `v-slip-rate` time-series against SeisSol's `tpv104_seisol_x2_0_x3_7.5.txt` column 6. Visual check of shape similarity in addition to peak-magnitude comparison.

---

### [R-904] [CRITICAL] [REJECT — fourth time] [H2 V1=0 lock and similar masking proposals] — Definitively wrong physics

**Category:** DEVIATION (provably wrong, not just masking)

**Description:**
Round 5 R-502, round 7 R-702, round 8 R-801 each rejected the V1=0 lock as masking. Round 9 strengthens this: SeisSol's reference data shows slip_dip is genuinely non-zero everywhere on the fault, especially at off-axis and corner stations (up to 41 cm). Forcing V1 = 0 in the MFEM solver is **provably wrong** — it would degrade agreement with SeisSol at every off-axis station even as it makes the hypocenter look perfect. This is not a TPV104-specific safeguard; it's a TPV104 verification regression.

The proposal must come out of every candidate set permanently. Any future round that lists it should be treated as a process error.

**Suggested fix:** if any production code was modified to implement H2, **revert immediately**. If the env-gate version was used for a one-off measurement, that's fine, but the result must not be cited as a "candidate fix that achieves zero slip_dip" — it achieves wrong slip_dip.

---

## Summary

- Critical issues: **2** (R-902 off-axis stations untested; R-904 H2 categorically wrong)
- Moderate issues: **1** (R-901 V_max metric audit)
- Low issues: **1** (R-903 V_dip rate timing)
- **Plan compliance:** R-805 (acceptable slip_dip) RESOLVED — 60 µm at hypocenter is at SeisSol's own level. R-803 (V_max convergence) PARTIALLY resolved — needs per-station comparison, R-901.
- **Verdict:** the np=1 hypocenter slip_dip investigation is **scientifically resolved.** P=2 produces SeisSol-level accuracy at the hypocenter. Round-6 sub-step iteration (R-602/603/604) is **NOT JUSTIFIED** by these data. The H2 lock proposal is **definitively wrong physics**, not just masking. Two follow-ups remain: R-901 V_max audit (5 minutes), R-902 off-axis station comparison (a longer run + 9 awks).

---

## Direct answers to the standing questions

**"Run E-3 (--order 3)?"** **No.** P=2 already at SeisSol's level on slip_dip at the hypocenter. E-3 cost is high (~5× P=2) and benefit is small (further reduction below an already-acceptable noise floor). Worth running ONLY if R-902 reveals systematic shortfall at off-axis stations AND the shortfall is plausibly precision-floor-bound rather than physics-bound.

**"Is sub-step iteration (round-6) justified?"** **No, by these data.** Sub-step iteration is a structural change targeting loop-gain regulation. The data shows:
- P=2 reduces slip_dip from 10 mm to 58 µm — at SeisSol's level.
- Sub-step iteration would reduce loop gain further but at significant cost (~250 LOC + full BP5 regression gate).
- ROI calculation: sub-step iteration's benefit is "marginal accuracy improvement past SeisSol's level"; cost is days of engineering. Defer.

**"Is 60 µm acceptable?"** **Yes — it's at SeisSol's own numerical accuracy at the hypocenter.** The agent's earlier framing ("SeisSol claims 8 sig figs of zero, gap is 3 orders of magnitude") was wrong. SeisSol's own slip_dip at hypocenter peaks at 50 µm within t≤2s. We're matching them.

**"Should we adopt SeisSol-style sub-step iteration?"** **Not for this bug.** The original motivation (R-501 confirmed loop-gain instability) is real, but P=2 default already brings the symptom below SeisSol's own level. Sub-step iteration becomes worth doing if/when a NEW concern surfaces that P=2 can't address — not from this investigation.

---

## Recommended next actions, in order

1. **Retract H2 V1=0 lock from any code path it touched** (R-904). Revert if present in production. Acknowledge in any follow-up document that the proposal is wrong physics, not just masking.

2. **Run R-901 V_max audit** (5 minutes, awk only). Compare MFEM hypocenter peak h-slip-rate to SeisSol's 7.39 m/s. If MFEM ≈ 7-8 m/s at the station, the V_max metric mismatch is benign and we can declare P=2 fully converged at hypocenter. If MFEM ≈ 15-17 m/s at the station, the rupture genuinely overshoots SeisSol — separate investigation.

3. **Run R-902 off-axis station comparison.** Re-run MFEM P=2 for tfinal=12 s (or grab the existing `plots_results_*` data if a 12s P=2 run was already done). Compare peak |v-slip| at all 9 SCEC stations against SeisSol's reference. Pass: ratios in [0.5, 2.0]. Fail: any station outside that range deserves drill-down.

4. **If R-901 + R-902 both pass:** ship `--order 2` as the TPV104 default after the round-6 R-605 regression suite (BP5 / TPV102 smoke at P=2). Document the resolution in a new debug doc; close this investigation.

5. **If R-902 fails at off-axis stations:** new investigation. Different mechanism than the np=1 hypocenter loop-gain bug.

6. **R-602/603/604 sub-step iteration adoption:** not justified by current data. Defer until and unless a future investigation shows P=2 alone can't reach SCEC compliance globally. Round-6 plan is recorded; it's not in the queue.

---

## Lessons (round 9, eight rounds in)

1. **Verify reference numbers before building a story on them.** The agent's "SeisSol claims slip_dip ≡ 0 to 8 sig figs" was wrong. Reading the actual benchmark traces took ~5 minutes and overturned three rounds of analysis built on that premise. Always read the reference data, not the claim about the reference data.
2. **A fix is wrong if it disagrees with the reference at unrelated metrics.** H2's "achieves slip_dip = 0 exact" looked like a strict cure; the reference data shows zero is the WRONG answer at 8 of 9 stations. Mask proposals don't just hide the bug — they introduce *new* disagreement that wasn't there.
3. **Per-station comparison dominates global statistics.** V_max-over-all-DOFs vs SeisSol's per-station peak is incomparable. The right comparison is hypocenter-vs-hypocenter, station-by-station. Should have been the first move; we got there in round 9.
4. **Investigation rounds shrink as the data narrows the question.** Round 1 had no data; round 5 had R-501 which narrowed to "loop gain"; round 8 had E-2 which narrowed to "precision floor regulator"; round 9 has the SeisSol comparison which CLOSES the slip_dip-at-hypocenter question. Each round halved the search space. The cost of NOT reading the reference data earlier was 4 rounds of structural-change speculation.

---

## Unreviewed Areas

- MFEM hypocenter peak h-slip-rate (per-station, NOT V_max) at P=1, P=2, on the dx=500 m repro — gates R-901.
- MFEM peak |v-slip| at all 9 SCEC stations under P=2 — gates R-902.
- The P=2 production-scale run (Frontera np=400 / dx=200 m / 12 s) — gates the actual ship decision (R-806 carry-forward).
- SCEC TPV104 official acceptance tolerance (the spec PDF). Worth a one-time look-up if the program goal is "pass SCEC verification" rather than "match SeisSol numerically."
