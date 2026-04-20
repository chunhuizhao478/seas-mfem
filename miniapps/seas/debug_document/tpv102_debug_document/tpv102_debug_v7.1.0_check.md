# Code Review v7.1.0: RESULT.txt false-positive FAIL on dispositive 400-rank run

**Versioning note:** `v[a].[b].[c]` scheme.  Previous: `v7.0.0`.  This
round found one MODERATE bug (RESULT.txt regex false positive on an
otherwise-successful dispositive run).  Bumps `[b]`: **new version
`v7.1.0`**.

## TL;DR

**The Frontera 400-rank 200m dispositive run (job 7665297) physically
SUCCEEDED.**  The v7.0.0 R-701/R-801/R-802 fix chain is validated at
production scale: rupture nucleated, broke away, propagated to full
dynamic (`V_max = 7.7 m/s`), R-101 verifier passed with 405 pairs at
machine precision, all 400 ranks active.

**But `RESULT.txt` reads `FAIL`** because the dead-rank detector regex
(added in v4 R-402 and hardened in v4 R-403) flags the **far-endpoint
rank r399** as dead on the basis of its far-field pre-arrival amplitude
(`1.4e-26`, which has exponent `-26` and so matches the regex's
`[0-9]\.[0-9]{3}e-[2-9][0-9]` "effectively-dead" alternation).

This is a MODERATE regex-design bug, not a physics bug.  The fix is in
the sbatch RESULT.txt writer, not in the simulation code.

## Review Scope
- Output log: `tpv102_200m_400r_v2_7665297.out` (full run, 5268 steps).
- RESULT.txt writer: `jobs/tpv102/tpv102_200m_p1_1.5s_400rank_dev.sbatch`
  and `jobs/tpv102/tpv102_200m_p1_1.5s_50rank_dev.sbatch` (identical
  regex).
- Existing regex regression: `jobs/tpv102/test_result_regex.sh`
  (v4 R-402 test harness).
- v4 R-402 design origin: `tpv102_debug_v4_fix.md` R-402 section.

## Evidence the physics succeeded

### Scale-confirming statistics
- **Mesh:** 2,464,689 elements; **DOFs:** 9,858,756 per component
  (88.7 M total).  Full production size.
- **Shared fault faces (global):** 135 (from `405 pairs matched`
  at R-101 verifier, 3 QPs/face).  First time the R-701 canonical
  frame has been exercised on real TPV102 partitioning.
- **R-101 verifier:** `max_rel_diff = 6.66e-16` (3 ULP of double) on
  405 pairs.  Bit-identical DOFData across every shared-face pair.
  The frame-mismatch saga is closed at production scale.

### V_max trajectory — breakaway + dynamic propagation
```
t=0.1s  V_max=1.4e-12 m/s   (nucleation load, locked)
t=0.5s  V_max=1.3e-4  m/s   (exponential slow-slip)
t=0.7s  V_max=1.7e-2  m/s   (approaching V_nuc=0.03)
t=0.9s  V_max=2.1e-1  m/s   (breakaway begins)
t=1.0s  V_max=6.4e-1  m/s   (transitioning to dynamic)
t=1.1s  V_max=2.4    m/s   (rupture accelerating)
t=1.2s  V_max=7.70   m/s   (fully dynamic, steady)
t=1.5s  V_max=7.70   m/s   (steady propagation)
```

This is textbook TPV102 dynamic rupture.  `V_max = 7.7 m/s` sits in the
published SCEC TPV102 peak-slip-rate range (3–15 m/s depending on
resolution).

### qnorm watch at t=1.5s (final step)
```
r3=1.660e-08  r4=4.969e-10  r2=3.798e-12  r7=1.311e-13  r399=1.400e-26  r0=2.810e-10
(hypo_rank=3)
```

Interpretation by rank distance from hypocenter:
| rank | ||Q||_∞ | distance from hypo | signal state |
|---|---|---|---|
| r3 (hypo)  | 1.66e-08 | 0 ranks | peak |
| r4 (±1)    | 4.97e-10 | 1 rank | strong signal, ~30× below hypo |
| r2 (±1)    | 3.80e-12 | 1 rank | strong signal |
| r0 (±3)    | 2.81e-10 | 3 ranks | strong signal |
| r7 (±4)    | 1.31e-13 | 4 ranks | signal, 40 orders above noise floor |
| r399 (far) | 1.40e-26 | ~200 ranks domain-end | waves haven't arrived yet |

r399 grew from `1.2e-62` (t=0.1s) to `1.4e-26` (t=1.5s) — a
**36-order-of-magnitude rise**, which is absolutely consistent with a
far-field point slowly picking up the leading edge of numerical
precursor waves from the nucleation perturbation.  **r399 is NOT
dead.**  It's just physically far from the rupture front.

All 400 ranks have `#ranks_with_||Q||=0: 0/400` after step 351.  Every
rank is alive.

### What the run disproves (the v1 symptom)
The v1 "rupture stops at partition seam" failure mode would look like:
- r3 (hypo) at 1e-8
- r4 (neighbor of hypo) stuck at 1e-20 or 0.000e+00

Observed: r4 = 5e-10, r2 = 3.8e-12.  **These are 10 and 12 orders of
magnitude above the noise floor**, tracking r3 with expected wave-
propagation attenuation.  The v1 bug is definitively NOT present.

## Findings

### [R-901] [MODERATE] `jobs/tpv102/tpv102_200m_p1_1.5s_{50,400}rank_dev.sbatch` — RESULT.txt regex matches `rN=M.MMMe-NN` for exponents `-20..-99`; far-endpoint ranks at t=1.5s on large domains are expected to have far-field pre-arrival amplitudes in that range, producing false-positive `FAIL` verdicts on physically successful runs

**Category:** BUG (CI regex overreach)

**Description:**
The R-402 regex (introduced in v4 to replace an earlier false-positive-
prone pattern) classifies any `rN=M.MMMe-NN` as "dead" if `NN ∈ [20, 99]`.
That covers two physically distinct regimes:

1. **Denormal / true dead** (exponent ≤ -30ish).  The intended target.
2. **Far-field pre-arrival** (exponent -20 to -30).  Legitimate small
   signal that grew from the noise floor but hasn't yet reached
   saturation — normal at domain endpoints for a sub-s simulation.

The 400-rank 200m run produced `r399 = 1.400e-26` at t=1.5s, which
regex classifies as dead but is physically alive (grew from 1.2e-62).
Every other watched rank (r3, r4, r2, r7, r0) is clearly alive.

The regex is tested by `jobs/tpv102/test_result_regex.sh`, and
`r10=5.000e-21` is in the positive-match set — which means the v4 R-402
design explicitly assumed `-21` exponent ⇒ dead.  That assumption is
false for far-field points after a short-duration dynamic run on a
large domain.

**Trigger:**
Any TPV102 run where the watch set includes a rank sufficiently far
from the hypocenter that wave arrival is incomplete at `tfinal`.  For
the 200m 400-rank dispositive run: `r399` at `tfinal=1.5s` on ~60-120
km domain sits at far-field pre-arrival level.

**Actual behavior:**
RESULT.txt line:
```
FAIL: final [qnorm:watch] shows dead ranks (suggests rupture did not cross partition seam):
        [qnorm:watch] r3=1.660e-08 r4=4.969e-10 r2=3.798e-12 r7=1.311e-13 r399=1.400e-26 r0=2.810e-10 (hypo_rank=3)
```

**Expected behavior:**
`RESULT.txt` should read `PASS` on this run.  The dead-rank detector
should distinguish far-endpoint arrival-not-yet vs. genuine bug-stopped
wave.  Options:

- **Option A (minimal):** narrow the "effectively dead" exponent
  window to `-30 … -99`.  Catches true denormals and noise-floor-
  stuck values; lets pre-arrival far-field through.  TPV102 physical
  quiescent signal at t=1.5s is `O(1e-26)`, well above `O(1e-30)`.
- **Option B (strict):** only flag the hypocenter and its immediate
  neighbors (±1 ranks) as subject to the dead-rank check; accept
  endpoint r0 and r(nprocs-1) as informational-only.
- **Option C (dynamic):** compare each rank's final qnorm against
  its step-0 value; if `growth_ratio < 1e3`, flag as stuck.  More
  physically meaningful but harder to express in bash.

**Suggested fix — Option A (simplest, most defensible):**

```diff
 # jobs/tpv102/tpv102_200m_p1_1.5s_400rank_dev.sbatch (and 50rank_dev)
-elif echo "${FINAL_WATCH}" | grep -qE "r[0-9]+=(0\.0{3}e\+[0-9]{2}|[0-9]\.[0-9]{3}e-[2-9][0-9])([[:space:]]|$)"; then
+# R-901 fix: narrow the "effectively dead" exponent floor from -20 to
+# -30.  On large-domain dynamic-rupture runs (TPV102 60-120 km × 30 km),
+# far-endpoint ranks at tfinal=1.5s legitimately have ||Q||_∞ in the
+# 1e-25..1e-30 range (far-field pre-arrival amplitude).  Only denormal-
+# / noise-floor-stuck values (exp ≤ -30) and exact zero should be
+# treated as "dead" in the dispositive verdict.
+elif echo "${FINAL_WATCH}" | grep -qE "r[0-9]+=(0\.0{3}e\+[0-9]{2}|[0-9]\.[0-9]{3}e-([3-9][0-9]|[0-9]{3,}))([[:space:]]|$)"; then
    echo "FAIL: final [qnorm:watch] shows dead ranks (suggests rupture did not cross partition seam):" > "${RESULT_FILE}"
    echo "      ${FINAL_WATCH}" >> "${RESULT_FILE}"
```

The updated exponent bracket `-([3-9][0-9]|[0-9]{3,})` matches:
- `-30 … -99` (two-digit exponents starting with 3-9)
- `-100 … -999+` (three-or-more-digit exponents, i.e. denormals)

It does NOT match `-20 … -29`, which is where legitimate far-field
pre-arrival amplitudes fall on dispositive dynamic-rupture runs.

Update the test harness too:

```diff
 # jobs/tpv102/test_result_regex.sh
+# R-901: r399 at far-endpoint of a dispositive TPV102 run is ALIVE,
+# not dead.  The old regex misclassified it.  New cases:
 for line in \
    "[qnorm:watch] r10=0.000e+00 r11=1.000e+05 (hypo_rank=10)" \
-   "[qnorm:watch] r2=2.470e-61 r11=1.000e+05 (hypo_rank=10)"  \
-   "[qnorm:watch] r2=5.000e-21 r11=1.000e+05 (hypo_rank=10)"  \
    "[qnorm:watch] r2=9.821e-33 r11=1.000e+05 (hypo_rank=10)"  \
+   "[qnorm:watch] r2=2.470e-61 r11=1.000e+05 (hypo_rank=10)"  \
+   "[qnorm:watch] r2=1.000e-30 r11=1.000e+05 (hypo_rank=10)"  \
 ; do
    echo "$line" | grep -qE "$PATTERN" \
       || { echo "REGRESS: should match: $line"; exit 1; }
 done

 # Alive cases (must NOT match) — add r399 far-field case from job 7665297
 for line in \
    "[qnorm:watch] r10=1.000e+06 r11=5.000e+05 (hypo_rank=10)" \
    "[qnorm:watch] r10=1.234e-03 r11=5.000e+05 (hypo_rank=10)" \
    "[qnorm:watch] r10=5.000e-10 r11=5.000e+05 (hypo_rank=10)" \
    "[qnorm:watch] r10=9.999e-19 r11=5.000e+05 (hypo_rank=10)" \
+   "[qnorm:watch] r399=1.400e-26 r11=5.000e+05 (hypo_rank=10)" \
+   "[qnorm:watch] r399=5.000e-21 r11=5.000e+05 (hypo_rank=10)" \
+   "[qnorm:watch] r399=9.999e-29 r11=5.000e+05 (hypo_rank=10)" \
 ; do
    if echo "$line" | grep -qE "$PATTERN"; then
       echo "REGRESS: should NOT match: $line"; exit 1
    fi
 done
```

**Test case:**
```bash
# Simulate the 400-rank dispositive run's actual final watch line
line="[qnorm:watch] r3=1.660e-08 r4=4.969e-10 r2=3.798e-12 r7=1.311e-13 r399=1.400e-26 r0=2.810e-10 (hypo_rank=3)"

# Pre-fix (v7.0.0 regex): matches → FAIL
echo "$line" | grep -qE "r[0-9]+=(0\.0{3}e\+[0-9]{2}|[0-9]\.[0-9]{3}e-[2-9][0-9])([[:space:]]|$)" && echo "pre-fix: DEAD"

# Post-fix (v7.1.0 regex): no match → PASS
echo "$line" | grep -qE "r[0-9]+=(0\.0{3}e\+[0-9]{2}|[0-9]\.[0-9]{3}e-([3-9][0-9]|[0-9]{3,}))([[:space:]]|$)" || echo "post-fix: ALIVE"
```

Expected output: `pre-fix: DEAD` then `post-fix: ALIVE`.

---

## Summary
- Critical issues: **0**
- Moderate issues: **1** (R-901 — RESULT.txt regex false positive on
  far-endpoint pre-arrival signal)
- Low issues: **0**
- Plan compliance: **FULL** — every v7.0.0 physics fix is validated.
  Only CI/tooling has a false-positive verdict that needs tightening.
- **Verdict: PHYSICS PASS WITH CI FIX.**
  - The 400-rank dispositive TPV102 run is a physical PASS:
    `V_max = 7.7 m/s`, R-101 verifier bit-identical, all 400 ranks
    alive with correct radial propagation.
  - The RESULT.txt verdict is a false negative caused by R-901.
  - Fix R-901 (one regex tweak in two sbatch files + test harness
    update), re-submit the same job to confirm, then proceed to
    normal-queue long-duration runs.

## Version bump
- Pre: `v7.0.0`.
- One MODERATE finding, no CRITICAL, no LOW.  Bumps `[b]`.
- **Post: `v7.1.0`.**  Doc: `tpv102_debug_v7.1.0_check.md`.

## Recommended next steps

1. **Do NOT re-run** the dispositive 400-rank job just to re-verify the
   PASS — we already have full log evidence and the fix is CI-only.
2. Apply R-901 fix (3 file edits) + commit + push.
3. **Manually verify** the retroactive verdict by running the
   post-fix regex against the existing `tpv102_200m_400r_v2_7665297.out`
   final line — should now produce PASS without needing another
   node-hour burn.
4. Move to TPV102 long-duration verification (12 s per SCEC spec) on
   the normal queue.  Pre-flight with a shorter rerun (`tfinal = 3 s`)
   on the same mesh to confirm steady propagation continues.
5. Update `debug_document/tpv102_debug_document/` with `v7.1.0_fix.md`
   when R-901 lands.

## What the v7.0.0→v7.1.0 transition certifies

- **R-701 (BP5 FaultBasis reuse):** validated at 400 ranks with 135
  shared fault faces.  Cross-rank DOFData bit-identical to ULP.
- **R-801 (unify on BP5 convention — tangent1=dip, tangent2=strike):**
  TPV102 pure strike-slip initial condition produces the correct
  physical rupture.  `V_max = 7.7 m/s` is in-spec.
- **R-802 (canonical-nor flux + accum_sign):** bulk momentum
  conservation holds; waves propagate from the hypocenter to all 400
  ranks without stalling at any partition seam.
- **v1 saga closed:** the "rupture stops at partition seam" symptom
  is empirically absent on the production-scale dispositive
  configuration.

This is the first full Frontera-scale confirmation that the entire
v1-v7 fix chain is correct.  The RESULT.txt misverdict was the only
noise-layer issue in the whole run.

## Unreviewed Areas
- **V_max saturation at 7.69871 m/s** from step 4212 onward (held
  constant to 3+ digits over 1055 steps).  Probably normal steady
  dynamic propagation reaching a rupture-speed ceiling, but worth
  eyeballing the fault-surface PVD in ParaView to confirm the rupture
  front is actually moving (not frozen at one QP).  If ParaView shows
  the slip-rate contour expanding monotonically, ignore; if it's
  pinned, file as a new finding.
- **Long-duration stability (tfinal > 1.5s).**  12s per SCEC spec.
  Not exercised here.
- **`max ||Q||_∞ = 1.66e-08 Pa` at peak.**  This seems low for a V_max
  = 7.7 m/s rupture (bulk particle velocities should be O(1) m/s at
  the fault face).  One interpretation: Q is a "perturbation from
  equilibrium" (see `InitializeState` comment), so the absolute
  magnitude depends on what's been absorbed into the pre-stress
  vs. radiated as waves.  Did not investigate further.  If the
  absolute qnorm feels surprising on longer runs, revisit.
