# TPV102 / TPV104 strike-slip deficit (RK45 + mixed-flux spatial driver) — findings & recommendations

Date: 2026-06-02
Author: Claude (code-debug)
Driver: `seas_spatial_dyn_driver` (all runs below)

Runs analysed (MFEM-vs-SCEC overlay PNGs + station `.dat`):
  * **Subject (RK45):**
    - `tpv102/plot_combined_7762877_7763025`  (p1 · RK45 · mixed-flux adjacent · scalar)
    - `tpv104/plot_combined_7762891_7763034`  (p1 · RK45 · mixed-flux adjacent · scalar)
    - `tpv205/gold/plot_combined_7762879_7763019` (p1 · RK45 · mixed-flux adjacent · LSW) — **matches benchmark**
    - station data: `~/Downloads/seas-mfem/tpv{102,104}_spatial/combined_*`
  * **Controls (ADER gold, same driver):**
    - `tpv104/gold/results_nonsymmetric_mesh_mixedflux_adjacent_O2_job7680545`  (p1 · **ADER-O2** · mixed-flux adjacent — *same mesh+flux as subject*)
    - `tpv104/gold/results_nonsymmetric_mesh_mixedflux_none_O2_job7680536`      (p1 · ADER-O2 · upwind)
    - `tpv104/gold/results_symmetric_mesh_O2_job7678882`  /  `..._O3_job7678716` (ADER-O2 / **O3**, symmetric mesh)
    - `tpv102/gold/results_mixed_flux_adjacent_p1_O2_dev_job7681589`            (p1 · **ADER-O2** · mixed-flux adjacent)
  * **Reference:** `tpv104/benchmark_data/{DRDG3D(O4),seisol(O5)}`, `tpv102/benchmark_data/{scec_drdg3d(O4),scec_pylith}`
  * Repro script: `tpv104_debug_document/compare_rk45_vs_ader_slip.py` (this dir)

Jobs: `jobs/tpv{102,104,205}_spatial/tpv*_p1_rk45_mixedflux_8N_400r_dev.sbatch`

> **No code or config was changed.** This is a diagnosis + recommendations report
> (per CLAUDE.md: report findings and ask before applying fixes). User elected
> "just the report" on 2026-06-02.

---

## 1. Bottom line

The smaller strike slip in the **TPV102 / TPV104 RK45 runs is numerical
under-resolution of the *stiff* rate-and-state rupture front in time — NOT a
friction-model, sign-convention, config, or σ_n bug.** Two superposed effects:

1. **(dominant) The RK45 path runs at a 3× larger macro-step than the validated
   ADER-O2 run and does not sub-cycle the friction within the step.** On the
   *identical* p1 mesh + adjacent flux, **ADER-O2 matches TPV102 to −1.5%** while
   **RK45 is −16%**; the RK45 deficit **grows with distance from the hypocenter**
   (−6 % center → −28 % corners) — the signature of a rupture front progressively
   over-damped step-by-step.
2. **(secondary, TPV104 only) genuine p1 spatial under-resolution** of TPV104's
   strong-rate-weakening front: ADER-O2 −3 % → **ADER-O3 −0.1 %**.

**ADER-O3 reproduces the benchmark to <1 % at every station**, which proves the
rate-state physics, configs, σ_n, sign conventions, and nucleation are all
correct. **TPV205 (LSW) matches at the identical RK45/CFL** because LSW friction
is *non-stiff* and tolerates the large `dt`; rate-state friction is stiff and
does not.

---

## 2. Symptom (quantified)

Strike slip (horizontal/`h-slip`) vs reference, interpolated at a common time
(TPV104 t≈9.87 s, the gold ADER-O2 end; TPV102 t=10 s):

### TPV104  (deficit vs mean of DRDG3D-O4 + SeisSol-O5)
| station | RK45 p1 mfadj | ADER-O2 p1 mfadj | ADER-O2 sym | **ADER-O3 sym** |
|---|---|---|---|---|
| x2_0_x3_7.5 (hypo) | −7.9 % | −3.1 % | −3.9 % | **−0.1 %** |
| x2_0_x3_3          | −11.2 % | −6.8 % | −8.8 % | **−0.1 %** |
| x2_-9_x3_7.5       | −16.6 % | −3.9 % | −5.1 % | **−0.7 %** |
| x2_-12_x3_3        | −19.5 % | −6.0 % | −7.3 % | **−0.7 %** |

### TPV102  (deficit vs mean of DRDG3D-O4 + PyLith)
| station | RK45 p1 mfadj | **ADER-O2 p1 mfadj** |
|---|---|---|
| flt_0_7.5 (hypo) | −6.4 % | **−1.7 %** |
| flt_0_12         | −15.6 % | **−2.2 %** |
| flt_±12_3        | −16.3 % | **−1.8 %** |
| flt_±12_12 (corner) | −28 % | **−0.5 %** |
| **MEAN (9 stns)** | **−16.0 %** | **−1.5 %** |

Rupture front also arrives **~0.3 s late** in ours vs benchmark (both ADER and
RK45) — consistent with a slightly over-damped/under-resolved front.

---

## 3. The three discriminating controls

All three TPV runs use the SAME driver, SAME numerics (p1 / mixed-flux adjacent /
scalar Riemann / `--deriv-cache`), SAME 200 m mesh resolution. Only the friction
law differs. That isolates the cause:

1. **Order ladder → physics is correct.** O3 ≈ ref > O2 > RK45-p1, monotone in
   accuracy. **ADER-O3 matches the benchmark to <1 %** → no model/config/sign bug;
   the p1/O2 gap is discretization error, not wrong physics.
2. **Time integrator, everything else fixed.** Same p1 + same adjacent flux + same
   nonsymmetric 200 m mesh: **ADER-O2 = −1.5 % (TPV102), RK45 = −16 %.** So p1
   spatial resolution is *not* the limiter for the RK45 gap; the RK45 path itself
   loses ~14 %. The deficit **grows with propagation distance** (uniform spatial
   error would not).
3. **Friction-law sensitivity.** TPV205 (LSW) matches at the identical RK45/CFL.
   LSW μ(δ) is bounded/linear (non-stiff); rate-state has the stiff direct effect
   `a·ln(V/V₀)` + state evolution (and TPV104's strong-rate-weakening collapse),
   which needs fine time resolution at the front.

---

## 4. Mechanism — the `dt` is 3× the validated ADER step

`WaveOperator::ComputeMaxDt` returns `cfl_mixed_flux_factor · cfl · h_min / c_p`
(`dynamic/wave_operator.inl:5880-5910`). The driver feeds it
`cfl · {RkCflFactor | CflSafetyFactor}` (`drivers/spatial_dyn_driver.cpp:1903-1906`):

| | base `cfl` | order factor (N=1) | adjacent flux factor | **effective CFL** |
|---|---|---|---|---|
| **RK45** (subject) | 0.25 | `RkCflFactor = 3/(2N+1)=1.0` | 0.6 (`cfl_rk_aware_`) | **0.15** |
| **ADER-O2** (gold) | 0.5 | `CflSafetyFactor = 1/(3(2N+1))=1/9` | 0.9 | **0.05** |

→ **`dt_RK / dt_ADER = 0.15 / 0.05 = 3.0`.** (`dt_max=0.05 s` in the configs never
binds; real `dt ~ ms`.) On top of the 3× larger step, the ADER path additionally
sub-steps the friction iterator *within* each (already 3× smaller) macro-step,
whereas the RK path solves friction only at its DP45 stage nodes. Net: the RK45
run resolves the stiff friction front at ≥3× coarser effective time resolution →
the slip-rate spike at the front is truncated → lower peak V → less accumulated
slip + slower rupture, accumulating with propagation distance.

These CFL numbers are exactly what the job headers and `ComputeMaxDt`/`RkCflFactor`
comments flag as **"a STARTING calibration … read the empirical stable dt before
committing normal-queue hours"** — the production validation was deferred to the
Frontera run, which is the run analysed here.

---

## 5. What is NOT the cause (ruled out)

- **σ_n drift / fault-locality partition** (see `tpv31_debug_document/
  tpv31_sigma_n_drift_findings_2026-06-02.md`). Measured on-fault σ_n drift is
  **0.005 % (TPV104: 120.000→120.006 MPa)** and **~0.2 % (TPV102→120.2 MPa)** —
  three orders too small to explain an 8–28 % slip deficit, and **present in the
  ADER runs too** (which match). It is a real, separate cleanliness issue but is
  **not** why slip is low here.
- **Spurious dip slip**: sub-mm (TPV104 −0.54 mm, TPV102 −1.5 mm) against 8–13 m
  of strike slip. Negligible; present in ADER too.
- **Config / parameter / sign / nucleation errors**: excluded by ADER-O3
  (TPV104) and ADER-O2 (TPV102) matching the benchmark with the same configs.

---

## 6. Code review — the RK45 stepper is correct

`AdvanceRKCoupled_Spatial` (`dynamic/rk_time_stepper.hpp:182-328`) was reviewed:
stage-local ψ written before `wave.Mult` (so the ψ-stateless fault Riemann sees
the staged value), correct b-row combine for ψ and `slip{1,2} += dt·Σ b_i V_i`,
and FSAL endpoint handling (DP45 skips the re-eval; classical RK4 re-evaluates at
`Q_new`). No accumulation / sign / phase-lag bug. The under-prediction is a CFL
*calibration* issue, not a coding error.

---

## 7. Recommendations (need a Frontera run to verify — no local full-mesh runs)

Ordered by cost:

1. **CFL sweep (decisive, cheap).** Re-run one TPV102 and one TPV104 dev job with
   `cfl ≈ 0.08` (≈ `0.25/3`, matching ADER's effective 0.05). Expected:
   - slip climbs to ADER-O2 levels ⇒ confirmed *pure-`dt`*; lower the production CFL.
   - slip stalls short ⇒ the RK path needs **within-step friction sub-cycling**
     (structural; keeps the large `dt` for throughput).
2. **TPV104 last few %**: needs higher spatial order (O3) or a finer mesh; this is
   genuine p1 limitation, independent of (1). The benchmark codes are O4/O5.
3. **Strategic**: for **Adjacent** flux (~5–10 % central faces) **ADER-O2 is
   already stable and more accurate** (`RkCflFactor`/`ComputeMaxDt` comments + the
   gold runs confirm). The RK integrator is only *required* for **AllContinuous**
   (~95 % central) flux. Consider reserving RK45 for AllContinuous and using
   ADER-O2 for adjacent-flux rate-state production.

---

## 8. Reproduce

```
python3 miniapps/seas/debug_document/tpv104_debug_document/compare_rk45_vs_ader_slip.py
```
Prints the per-station strike-slip deficit tables in §2 for TPV104 (and the TPV102
companion). Column note: **TPV104** `.dat` use the SCEC layout
(`t h-slip h-slip-rate …`, strike = col idx 1); **TPV102** `.dat` use the internal
layout (`time slip1(dip) slip2(strike) …`, strike = col idx **2**); benchmark
files always have strike = `h-slip` = col idx 1.
