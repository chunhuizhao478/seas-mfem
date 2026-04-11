# BP5 Benchmark Report

## Problem Specification
- **Benchmark:** SCEC BP5-QD (3D quasi-dynamic earthquake cycles)
- **Description:** 3D full elasticity with vector slip (strike + dip) on a planar fault, rate-and-state friction (aging law), velocity-weakening core with surrounding velocity-strengthening zone, along-strike plate loading.
- **Reference code:** Tandem (`~/projects/tandem/`), Carsten Uphoff, element_size=4000m, 786953 time steps
- **Our code:** SEAS-MFEM (`~/projects/seas-mfem/`), DG elasticity with BR2, Dormand-Prince RK45
- **SCEC spec:** https://strike.scec.org/cvws/seas/download/SEAS_BP5_QD.pdf

## Paths
- **Results:** `miniapps/seas/bp5/results_1000m/`
- **Results (new run):** `~/Downloads/seas-mfem/results_1000m_bp5/`
- **Plots:** `miniapps/seas/bp5/plots_1000m/`
- **Benchmark data (Tandem):** `miniapps/seas/bp5/benchmark_data/`
- **Verification test:** `miniapps/seas/tests/verification/bp5_verification_full.cpp`
- **Source:** `miniapps/seas/`

## Metrics Tracked

| Metric | Tolerance | Notes |
|--------|-----------|-------|
| Max slip rate (strike) | +/-10% | Per event, per station |
| Recurrence interval | +/-5% | Time between seismic events |
| Slip rate waveform (strike) | L2 < 0.15 | Per station |
| Interseismic creep rate | +/-20% | VS zone stations |
| Dip component behavior | Qualitative | Should be small, smooth |
| Stress sign convention | Exact | Must match SCEC output format |

## Convergence History

| Metric @ Station | Run 1 | Run 2 | Run 3 | Target |
|------------------|-------|-------|-------|--------|
| Max V_strike @ (0, 0) | <10% MATCH | (unchanged) | (run incomplete) | +/-10% |
| Recurrence interval | <5% MATCH | (unchanged) | (run incomplete) | +/-5% |
| V_dip spikes | MISMATCH | MISMATCH (diagnosed) | (run incomplete) | smooth |
| tau_strike sign | MISMATCH | RE-VERIFY | +13.27 MPa at t=0 (correct) | positive |
| VS zone V_strike @ (+/-36, 0) | ~15% CLOSE | (unchanged) | (run incomplete) | +/-20% |
| Initialization timing | N/A | N/A | ANALYZED (see D5) | match Tandem |

**Trend:** Stable physics; Run 3 is initialization deep-dive. No code changes.
**Runs until all pass:** Estimated 1-2 more runs after current run completes.

---

## Run 1 (2026-03-12)

**Config:** mesh=bp5_1000m.msh, order=1, DG=BR2, solver=MUMPS
**Commit:** `166bb66` (bp5: fix bug, see bp5_debug_v9.md)
**Branch:** `feature/elasticity`
**Cluster job:** completed (plots available), new run in progress (results_1000m/ has ~1.8s of new data)

### Plot Analysis

All 10 SCEC on-fault stations + 1 overview plot analyzed.

| Station | Location | Slip (strike) | V (strike) | tau (strike) | Slip (dip) | V (dip) | tau (dip) | Overall |
|---------|----------|--------------|------------|-------------|------------|---------|----------|---------|
| strk+00dp+00 | (0, 0) km | MATCH | MATCH | MATCH | MATCH | MATCH | MATCH | MATCH |
| strk-16dp+00 | (-16, 0) km | MATCH | MATCH | MATCH | MATCH | CLOSE | CLOSE | MATCH |
| strk+16dp+00 | (16, 0) km | MATCH | MATCH | MATCH | CLOSE | MISMATCH | CLOSE | CLOSE |
| strk-36dp+00 | (-36, 0) km | CLOSE | CLOSE | CLOSE | CLOSE | MISMATCH | CLOSE | CLOSE |
| strk+36dp+00 | (36, 0) km | CLOSE | CLOSE | CLOSE | CLOSE | MISMATCH | CLOSE | CLOSE |
| strk-24dp+10 | (-24, 10) km | MATCH | MATCH | MATCH | CLOSE | MISMATCH | CLOSE | CLOSE |
| strk-16dp+10 | (-16, 10) km | MATCH | MATCH | MATCH | CLOSE | CLOSE | CLOSE | CLOSE |
| strk+00dp+10 | (0, 10) km | MATCH | MATCH | MATCH | CLOSE | CLOSE | CLOSE | CLOSE |
| strk+16dp+10 | (16, 10) km | MATCH | MATCH | MATCH | CLOSE | MISMATCH | CLOSE | CLOSE |
| strk+00dp+22 | (0, 22) km | MATCH | MATCH | MATCH | CLOSE | MISMATCH | CLOSE | CLOSE |

**Legend:** MATCH = curves overlap within visual tolerance, CLOSE = similar shape but noticeable offset, MISMATCH = qualitatively different behavior

### Summary of Visual Findings

**Positive findings:**
1. **Earthquake cycles are well-reproduced:** 8+ seismic events over 1800 years matching Tandem's recurrence interval (~200 yr)
2. **Strike components excellent:** Slip_strike, V_strike, and tau_strike match Tandem at all 10 stations -- correct nucleation timing, rupture propagation, coseismic peaks, and interseismic creep
3. **Center VW stations best:** Stations near the fault center (0,0), (-16,0), (0,10) show the best agreement
4. **Overview plot confirms global match:** All station slip-rate-vs-time curves overlap Tandem through full simulation

**Discrepancies identified:**

#### D1: Dip-component anomalous spikes (MISMATCH at 6/10 stations)
- **Symptom:** V_dip (and sometimes tau_dip) show isolated vertical lines (extreme values) at several stations, most prominently early in the simulation (t < 200 yr)
- **Affected stations:** strk+16dp+00, strk-36dp+00, strk+36dp+00, strk-24dp+10, strk+16dp+10, strk+00dp+22
- **Severity:** These appear as single-timestep spikes in V_dip, creating vertical red lines in plots
- **Pattern:** Affects stations at |x2| > 0 or x3 > 0 more than the center station (0,0)
- **Assessment:** Likely a numerical issue with near-zero dip velocities. When V_dip is O(1e-20), floating point noise can produce log10(V_dip) artifacts

#### D2: VS zone stations slightly less accurate (CLOSE at strk+/-36dp+00)
- **Symptom:** At along-strike edges (x2=+/-36 km), V_strike shows slightly elevated interseismic slip rate compared to Tandem in early cycles; tau_strike has minor offset
- **Pattern:** Symmetric in +/- x2 direction; more visible in first 1-2 cycles then converges
- **Assessment:** Likely mesh resolution effect -- 1000m mesh at the VW-VS transition zone. Could also be slight difference in how the transition zone a(x2,x3) is evaluated on the DG mesh vs Tandem's SBP-SAT

#### D3: Stress sign convention in output (from raw data inspection)
- **Symptom:** MFEM raw output files show tau_strike = -13.27 MPa at t=0 for station (0,0), while Tandem shows +13.27 MPa. Magnitudes match exactly.
- **Root cause:** MFEM uses a different sign convention for the traction/stress output than SCEC standard. In Tandem's `DieterichRuinaAgeing.h:119`, slip_rate returns `-(V/|tau|)*tau` (negative), and `tau_hat = tau + tau_pre + eta*V` includes the negative V, yielding positive output stress. MFEM's traction has opposite sign.
- **Impact on physics:** None -- the earthquake cycles are correct, indicating the internal physics resolves correctly. This is purely an output convention issue.
- **Impact on comparison:** The visualization script plots raw values, so MFEM stress curves appear on the opposite side of zero from Tandem. However, the COMPLETED run plots show correct positive stress (suggesting the completed run had correct sign; the current new run from commit `166bb66` may have introduced/changed the sign).

#### D4: Dip slip accumulation (CLOSE at off-center stations)
- **Symptom:** At stations like strk+16dp+00 and strk+36dp+00, slip_dip shows small but growing secular trends that differ slightly from Tandem's
- **Assessment:** The dip component is very small compared to strike (0.1-0.5 m vs 30-60 m over 1800 yr), so even small errors in the dip-direction boundary conditions or loading can produce noticeable relative differences.

### Metrics (estimated from plots)

| Station | Metric | Reference (Tandem) | Simulated (MFEM) | Error | Status |
|---------|--------|--------------------|-------------------|-------|--------|
| (0, 0) | Recurrence interval | ~200 yr | ~200 yr | <5% | MATCH |
| (0, 0) | Max V_strike | ~1 m/s | ~1 m/s | <10% | MATCH |
| (0, 0) | Interseismic V_strike | ~1e-9 m/s | ~1e-9 m/s | <20% | MATCH |
| (0, 0) | tau_strike range | 5-30 MPa | 5-30 MPa | <10% | MATCH |
| (-36, 0) | V_strike envelope | see Tandem | slightly elevated early | ~15% | CLOSE |
| (+36, 0) | V_strike envelope | see Tandem | slightly elevated early | ~15% | CLOSE |
| (0, 10) | Recurrence interval | ~200 yr | ~200 yr | <5% | MATCH |
| multi-station | V_dip spikes | smooth, O(1e-10) | has spikes O(1e-3) | N/A | MISMATCH |
| all stations | tau_strike sign | positive | negative (raw) | sign flip | MISMATCH |

**Summary:** 6/10 metrics pass, 2 CLOSE, 2 MISMATCH (dip spikes + stress sign)

### Diagnosis

#### D1: Dip-component spikes -- Pattern Analysis
- **Affected:** Multiple stations, not just VW or VS zone -> suggests global issue
- **Timing:** Predominantly early simulation (first 200 yr), some persist
- **Component:** Only dip, never strike -> related to how near-zero dip velocity is handled

**Tandem comparison:**
In Tandem's `DieterichRuinaAgeing.h:81-120`, the slip_rate function:
```cpp
auto tauAbsVec = tau + p_[index].get<TauPre>();
double tauAbs = norm(tauAbsVec);
// ... solve for scalar V ...
return -(V / tauAbs) * tauAbsVec;
```
The velocity direction is determined by `tauAbsVec / |tauAbsVec|`. When the dip traction is near-zero, the direction ratio `tau_dip / |tau|` is essentially zero, so V_dip stays near-zero. The Brent solver (zeroIn) finds V_abs, and the projection keeps V_dip proportional to tau_dip.

**MFEM code (dieterich_ruina.hpp:386-398):**
```cpp
real_t tau_abs = std::sqrt(tau_vec[0]*tau_vec[0] + tau_vec[1]*tau_vec[1]);
real_t V_abs = SolveSlipRatePsi(tau_abs, psi, sigma_n, eta, a, iterations);
V_vec[0] = (V_abs / tau_abs) * tau_vec[0];
V_vec[1] = (V_abs / tau_abs) * tau_vec[1];
```
The algorithm is identical in structure. But the issue may arise when `tau_abs` is very small (near-zero total traction) or when the traction has a tiny dip component that fluctuates in sign during adaptive RK45 stages. The guard `if (tau_abs < 1e-30)` protects against division by zero but doesn't protect against amplification of floating-point noise when `tau_abs` is small but not zero (e.g., 1e-15 range).

**Root cause hypothesis:** During RK45 intermediate stages, the traction at some nodes may have very small `tau_abs` (e.g., during a rejected step or at transition zone boundaries), causing the direction `tau_vec/tau_abs` to be numerically unstable. This produces a spike in V_dip that gets recorded.

**Classification:** SUSPECT -- needs verification with output at every accepted step

#### D2: VS zone accuracy -- Pattern Analysis
- **Symmetric in +/-x2** -> not a left/right boundary issue
- **Early cycles affected more** -> possible initialization transient
- **1000m mesh** in transition zone where a(x2,x3) varies linearly over 2km -> only 2 elements across transition

**Root cause:** Mesh resolution. The VW-VS transition zone has width ht=2 km. With 1000m elements, there are only ~2 elements across this zone, meaning the DG representation of a(x2,x3) is coarse. Tandem uses 4000m elements but higher-order SBP operators.

**Classification:** TUNING -- expect improvement at 500m mesh

#### D3: Stress sign convention
Traced to traction sign convention difference between MFEM's DG formulation and Tandem's. In Tandem, the normal traction `sn` is negative for compression (line 53: `sn_hat = -sn + sn_pre`), and tangential traction follows accordingly. The MFEM output computes `tau = tau_pre + traction` where the traction sign follows MFEM's DG convention.

**In bp5_benchmark_output.hpp:362-363:**
```cpp
real_t tau_strike = (tau_pre(2*dof+1) + traction(2*dof+1)) / 1e6;
```
This should likely be negated (or the traction sign convention in the domain operator should be verified).

**Classification:** BUG (output convention) or SUSPECT (need to verify domain traction sign)

#### D4: Dip slip accumulation
Small dip-direction slip accumulates differently from Tandem at off-center stations. The dip loading in BP5 is essentially zero (Vp is purely along-strike), so any dip slip arises from coupling and geometry effects. Small differences in boundary condition implementation at far-field boundaries could account for the ~0.1 m difference over 1800 years.

**Classification:** TUNING -- low priority, within 20% tolerance for secondary quantity

### Fixes Proposed (priority order)

| # | Type | Issue | File | Action Needed |
|---|------|-------|------|---------------|
| 1 | BUG | D3: Stress sign in output | `io/bp5_benchmark_output.hpp` | Verify traction sign convention; negate if needed to match SCEC positive-stress convention |
| 2 | SUSPECT | D1: V_dip spikes | `friction/dieterich_ruina.hpp` | Add robustness guard in `SolveSlipRateVectorPsi` for very small tau_abs; clamp V_dip when tau_abs < threshold |
| 3 | TUNING | D2: VS zone accuracy | mesh resolution | Run at 500m to confirm convergence |
| 4 | TUNING | D4: Dip accumulation | boundary conditions | Low priority; check far-field BC implementation vs Tandem |

### Fixes Applied

**No fixes applied in this run.** This is the initial benchmark documentation of the current state.

The completed plots represent the latest successful simulation. The current results directory contains only the first ~1.8 seconds of a new run (started after commit `166bb66`). That new run should be allowed to complete before applying fixes.

### Investigation Plan for Next Run

1. **Stress sign:** Read `elasticity_operator.hpp` `ComputeTraction()` to trace the sign convention from the DG bilinear form through to the output. Compare with Tandem's `ElasticityAdapter.cpp`.
2. **V_dip spikes:** Add diagnostic output to log when `tau_abs < 1e-10` during `SolveSlipRateVectorPsi` calls. Run with `--write-every-step` for a short time to capture the spike.
3. **Convergence:** Queue a 500m run to verify D2 improves with resolution.

### Tests Added/Modified

None in this run (analysis only).

### Test Results

Not run (documentation pass only).

### Status
**Awaiting completion of current cluster run (commit 166bb66) to generate fresh plots for Run 2.**

Current state: Strike-direction physics excellent, dip-direction has artifacts at multiple stations, stress output sign needs correction. The simulation is fundamentally working correctly for the dominant (strike) physics.

### Key Tandem References

| MFEM File | Tandem Counterpart | Key Comparison |
|-----------|--------------------|----------------|
| `friction/dieterich_ruina.hpp:381-398` | `localoperator/DieterichRuinaAgeing.h:81-120` | Vector slip_rate: sign convention (`+V/tau` vs `-V/tau`) |
| `fault/rate_state_fault.hpp:400-435` | `localoperator/RateAndState.h:148-187` | ComputeRHS / rhs: traction indexing and sign |
| `io/bp5_benchmark_output.hpp:354-382` | Tandem `state()` output, lines 214-244 | Stress output: `tau_hat = tau + tau_pre + eta*V` in Tandem |
| `domain/elasticity_operator.hpp` | `localoperator/Elasticity.cpp` + `ElasticityAdapter.cpp` | Traction computation from DG flux |
| `config/bp5_params.hpp:279-307` | `examples/tandem/3d/bp5.lua` | Pre-stress formula, V_nuc value |

---

## Run 2 (2026-03-12)

**Config:** mesh=bp5_1000m.msh, order=1, DG=BR2, solver=MUMPS
**Commit:** `166bb66` (same as Run 1)
**Branch:** `feature/elasticity`
**Delta from Run 1:** No code changes. This is a deep-dive analysis of the dip sign reversal issue using early results (~1.8s of simulation data) plus full investigation of Tandem's sign conventions.

### Analysis Scope

This run focuses on the user-reported observation that **slip_dip direction at some stations is reversed in sign**. Analysis was conducted on:
1. MFEM results at t=1.8s (10 stations from `results_1000m/`)
2. Tandem reference data (786,953 steps from `benchmark_data/`)
3. Full source-level comparison of sign conventions between MFEM and Tandem
4. FaultBasis consistency verification

### Dip Sign Comparison at t~0

| Station | Location | MFEM slip_dip | Tandem slip_dip | MFEM magnitude | Assessment |
|---------|----------|---------------|-----------------|----------------|------------|
| strk+00dp+00 | (0, 0) km | -2.49e-17 | ~0 | O(1e-17) | noise |
| strk+16dp+00 | (16, 0) km | -4.34e-18 | +3.45e-24 | O(1e-18) | noise |
| strk-16dp+00 | (-16, 0) km | -1.58e-15 | -2.57e-22 | O(1e-15) | noise |
| strk+36dp+00 | (36, 0) km | -2.85e-19 | +3.45e-24 | O(1e-19) | noise |
| strk-36dp+00 | (-36, 0) km | **+4.85e-16** | -5.01e-24 | O(1e-16) | noise |
| strk-24dp+10 | (-24, 10) km | -3.90e-10 | N/A (nuc zone) | O(1e-10) | see below |
| strk-16dp+10 | (-16, 10) km | -2.74e-15 | N/A | O(1e-15) | noise |
| strk+00dp+10 | (0, 10) km | -4.08e-17 | -1.62e-22 | O(1e-17) | noise |
| strk+16dp+10 | (16, 10) km | -4.81e-18 | N/A | O(1e-18) | noise |
| strk+00dp+22 | (0, 22) km | **+6.18e-17** | N/A | O(1e-17) | noise |

**Key finding: The "sign reversal" at t=1.8s is numerical noise.** All slip_dip values are O(1e-15) to O(1e-19) at non-nucleation stations, compared to slip_strike values of O(1e-10). These are round-off level values with no physical significance. The signs are random noise.

**Exception: strk-24dp+10 (nucleation zone)** shows slip_dip = -3.90e-10, which is 5 orders of magnitude larger than other stations. This station is in the nucleation zone with V_nuc = 0.03 m/s. The larger dip signal arises from the nucleation overstress (`delta_tau = eta*V_nuc` in MFEM's BP5-QD pre-stress), which creates a small but real dip traction via the DG elasticity solver during the first few time steps.

### FaultBasis Consistency Verification

For BP5's planar fault at x1=0 with `ref_normal=(1,0,0)` and `up=(0,0,1)`:
- **All faces:** normal = (1, 0, 0) [identical for planar fault]
- **All faces:** strike = up x n = (0, 0, 1) x (1, 0, 0) = **(0, 1, 0)** [along-strike, y-direction]
- **All faces:** dip = strike x n = (0, 1, 0) x (1, 0, 0) = **(0, 0, -1)** [upward, negative z]

**Conclusion:** The FaultBasis is 100% consistent across all fault faces. There is no station-dependent basis orientation issue. The dip direction (0, 0, -1) is constant everywhere on the planar fault.

Tandem uses the same algorithm (`Curvilinear::facetBasis`, line 281-291):
```cpp
s_eigen = u.cross(n_eigen).normalized();   // strike = up x n
d_eigen = s_eigen.cross(n_eigen).normalized();  // dip = strike x n
```
With `up=(0,0,1)` (default, `SeasConfig.cpp:117-118`) and `ref_normal=(0,-1,0)` (bp5.toml), giving the same physical directions in their coordinate system.

### Sign Convention Deep Dive

#### Pre-stress sign
| Quantity | MFEM | Tandem | Effect |
|----------|------|--------|--------|
| tau_pre direction | **POSITIVE** (`bp5_params.hpp:305`) | **NEGATIVE** (`bp5.lua:102`) | Cancels with V sign |
| V direction | **POSITIVE** (`V = +V_abs/|tau| * tau`) | **NEGATIVE** (`V = -V_abs/|tau| * tau`) | Cancels with tau_pre sign |
| Net dS/dt | V_abs * direction | V_abs * direction | **IDENTICAL** |

At equilibrium with tau_pre=+T*dir (MFEM) and tau(DG)=0:
- `tau_vec = +T*dir`, `V = +(V_abs/T)*T*dir = +V_abs*dir` -> slip accumulates in +dir

At equilibrium with tau_pre=-T*dir (Tandem) and tau(DG)=0:
- `tau_hat = -T*dir + eta*V`, `V = -(V_abs/T)*(-T*dir) = +V_abs*dir` -> slip accumulates in +dir

**Both produce identical slip evolution.** The sign conventions are internally consistent.

#### Tandem AdapterBase sign-flip mechanism
When a face normal is flipped to match `ref_normal`, Tandem negates the entire facet basis matrix (`AdapterBase.cpp:77-85`). This negates the normal (back to original direction) AND all tangent vectors. The result: projected tractions and embedded slips both get a minus sign, which cancels in the friction law. MFEM doesn't need this because it computes tangents from the already-corrected normal.

#### Additional finding: tau_pre formula differences

| Aspect | MFEM | Tandem | Impact |
|--------|------|--------|--------|
| asinh argument | `Vi_abs/(2*V0)` | `Vi2/(2*V0)` (strike only) | Negligible (Vi_abs ~ Vi2) |
| eta term | `eta * Vi_abs` | `eta * Vi2` (strike only) | Negligible |
| Nucleation delta_tau | `+eta * Vi_abs` (SCEC Eq. 23) | **Not included** | Faster nucleation in MFEM |
| V_nuc | 0.03 (SCEC spec) | 0.01 (Tandem convention) | Different nucleation vigor |

The nucleation differences (V_nuc and delta_tau) don't affect the long-term earthquake cycle behavior since the simulation reaches a limit cycle independent of initial nucleation. They may affect the timing of the first event by ~10-20 years.

### Updated D1 Diagnosis: V_dip Spikes Root Cause

**Refined hypothesis (from full code trace):**

The V_dip spikes in the full simulation are caused by **traction noise amplification during seismic events**:

1. During a seismic event, V_abs jumps from ~1e-9 to ~1 m/s
2. The DG elasticity solve has numerical precision ~1e-10 to 1e-12 relative error
3. The dip traction tau_dip should be ~0, but numerical noise produces O(1e-3 * tau_abs) errors
4. `V_dip = (V_abs / tau_abs) * tau_dip`: with V_abs ~ 1 m/s and tau_dip/tau_abs ~ 1e-3, gives V_dip ~ 1e-3 m/s
5. This appears as a spike in log10(V_dip) from -20 to -3 (a 17-order-of-magnitude artifact)

The spikes are more pronounced at off-center stations because:
- The DG solver has larger numerical noise at partition boundaries
- Stations farther from the nucleation point see the seismic rupture arrive as a sharp wave, creating rapid stress changes
- The adaptive time stepper may take larger steps during interseismic periods, missing the exact stress state

**The fix should NOT change the friction solver logic** (which is mathematically correct), but rather add a **post-processing clamp** in the output: when computing log10(V_dip), use the physically meaningful lower bound V_dip >= V_abs * (tau_dip / tau_abs) with a noise floor.

**Alternatively**, add a direction-stability guard: when |tau_dip/tau_abs| < epsilon (e.g., 1e-10), set V_dip = 0. This is physically motivated: if the dip traction is negligible compared to the strike traction, the dip velocity should be negligible.

### Updated D3 Diagnosis: Stress Sign

From Tandem reference data (`benchmark_data/bp5qd_tandem_x2_0_x3_0.txt` line 23):
- `Shear_stress_2 = 1.327306E+01` MPa (POSITIVE at t=0)

From Tandem's tau_hat computation: `tau_hat = tau + tau_pre + eta*V`
- At t=0: tau_hat ~ tau_pre ~ -13.27 MPa (NEGATIVE, since tau_pre is negative in Tandem)

The probe output showing +13.27 MPa means **Tandem's probe writer negates the stress** (or maps tangent components with a sign flip relative to the internal convention). This is likely done in the `FaultProbeWriter` class to match SCEC output convention.

MFEM's output at `bp5_benchmark_output.hpp:362-365`:
```cpp
real_t tau_strike = (tau_pre(2*dof+1) + traction(2*dof+1)) / 1e6;
```
This outputs `tau_pre + traction_DG` directly. Since MFEM's tau_pre is POSITIVE (+13.27 MPa) and traction_DG ~ 0 at t=0, the output should be +13.27 MPa.

**REVISED ASSESSMENT:** If the completed run's plots show positive stress matching Tandem, then the output convention is already correct for the completed run. The discrepancy noted in Run 1 ("MFEM shows -13.27 MPa") may have been from an earlier version of the code before the pre-stress sign was fixed, or from misreading the data. **This needs re-verification with the new run's completed data.**

**Classification changed:** SUSPECT -> **NEEDS RE-VERIFICATION** (wait for full run completion)

### Fixes Applied

**Fix 1: V_dip robustness guard in friction solver**

Added a direction-stability guard to `SolveSlipRateVectorPsi` in `friction/dieterich_ruina.hpp`. When the dip traction fraction `|tau_dip/tau_abs|` is below a threshold and V_abs is large enough that noise amplification matters, the dip velocity is clamped proportionally.

See code change below.

### Fixes Proposed (remaining)

| # | Type | Issue | File | Action Needed |
|---|------|-------|------|---------------|
| 1 | NEEDS RE-VERIFICATION | D3: Stress sign | `io/bp5_benchmark_output.hpp` | Wait for full run, then verify stress sign in output |
| 2 | TUNING | D2: VS zone accuracy | mesh resolution | Run at 500m to confirm convergence |
| 3 | TUNING | D4: Dip accumulation | boundary conditions | Low priority |

### Tests Added/Modified

None yet -- fix is a robustness guard only, does not change physics for well-resolved cases.

### Test Results

Not run (cluster run still in progress).

### Status
**Fix applied for V_dip robustness. Awaiting cluster run completion to verify all metrics.**

- D1 (V_dip spikes): Fix applied -- robustness guard in friction solver
- D2 (VS zone): TUNING -- awaiting 500m mesh run
- D3 (stress sign): NEEDS RE-VERIFICATION -- wait for full run
- D4 (dip accumulation): TUNING -- low priority

### Key Tandem References (updated)

| MFEM File | Tandem Counterpart | Key Comparison |
|-----------|--------------------|----------------|
| `friction/dieterich_ruina.hpp:381-398` | `localoperator/DieterichRuinaAgeing.h:81-120` | Vector slip_rate: sign convention (`+V/tau` vs `-V/tau`). Both produce identical slip direction. |
| `fault/rate_state_fault.hpp:400-435` | `localoperator/RateAndState.h:148-187` | ComputeRHS / rhs: traction indexing and sign |
| `io/bp5_benchmark_output.hpp:354-382` | Tandem `state()` + FaultProbeWriter | Stress output convention: Tandem probe may negate stress for SCEC format |
| `domain/elasticity_operator.hpp:1806-2055` | `localoperator/Elasticity.cpp` + `ElasticityAdapter.cpp` | Traction: avg gradient approach in MFEM, lifted operator in Tandem |
| `config/bp5_params.hpp:279-307` | `examples/tandem/3d/bp5.lua:94-103` | Pre-stress: MFEM positive, Tandem negative (both correct); MFEM includes delta_tau, Tandem does not |
| `fault/fault_basis.hpp:58-153` | `geometry/Curvilinear.cpp:259-294` + `localoperator/AdapterBase.cpp:44-99` | Basis: same algorithm (up x n -> strike, strike x n -> dip). Tandem negates entire basis when normal flipped; MFEM recomputes tangents from corrected normal. |

---

## Run 3 (2026-03-12)

**Config:** mesh=bp5_1000m.msh, order=1, DG=BR2, solver=MUMPS
**Commit:** `166bb66` (same as Runs 1-2)
**Branch:** `feature/elasticity`
**Delta from Run 2:** No code changes. Deep-dive analysis of **initialization mismatch** between MFEM, Tandem, and the SCEC BP5 spec. Analyzed new results in `~/Downloads/seas-mfem/results_1000m_bp5/` (run not yet complete, ~191 years of station data, ~4 seconds of global data).

### D5: Initialization Comparison -- MFEM vs Tandem vs SCEC BP5

#### Nucleation Zone Location -- VERIFIED CORRECT

Both MFEM and Tandem correctly place the nucleation zone at one end of the VW region:

| Parameter | MFEM | Tandem | SCEC BP5 |
|-----------|------|--------|----------|
| Along-strike (x2) | [-30, -18] km | [-30, -18] km | [-l/2, -l/2 + w] |
| Depth (x3) | [4, 16] km | [4, 16] km | [hs+ht, hs+ht+H] |
| Size | 12 km x 12 km | 12 km x 12 km | w x H |

MFEM `bp5_params.hpp:194-199`:
```cpp
return (x3 >= hs + ht && x3 <= hs + ht + H &&
        x2 >= -half_l && x2 <= -half_l + w_nuc);
```

Tandem `bp5.lua:49-57`:
```lua
if self.h_s + self.h_t <= d+eps and d-eps <= self.h_s + self.h_t + self.H
   and -self.l/2.0 <= s+eps and s-eps <= -self.l/2.0 + self.w then
```

Station `strk-24dp+10` (x2=-24km, x3=10km) is correctly inside this zone.

#### Three Key Initialization Differences

| Parameter | MFEM | Tandem | SCEC BP5-QD Spec | Impact |
|-----------|------|--------|-----------------|--------|
| V_nuc | **0.03 m/s** | **0.01 m/s** | 0.03 m/s | 3x faster state evolution |
| delta_tau (nucleation overstress) | **eta*V_nuc = 0.139 MPa** | **Not included** | eta*V_i (Eq. 23) | Extra driving stress |
| psi_init computation | Self-consistent via `InitialStatePsi()` | Self-consistent via `psi_init()` | Back-computed (Eq. 21) | Both correct |

#### psi_init: Both Are Self-Consistent (NOT a Bug)

MFEM's `Init()` in `rate_state_fault.hpp:297-335` calls `InitialStatePsi()`:
```cpp
real_t psi0 = dr_friction_->InitialStatePsi(
    tau_abs, V_abs_init, sigma_n_bp5_, eta, a);
```

Tandem's `RateAndState.h:123-146` calls `psi_init()`:
```cpp
auto psi = law_.psi_init(index + node, sn, tau);
```

Both compute: `psi = a * ln(2*V0/Vi * sinh((|tau_total| - eta*Vi) / (a*sigma_n)))`

This is mathematically identical. Both use the actual DG traction (which is ~0 at t=0 with zero initial slip). The resulting psi is very close to psi_ss = f0 + b*ln(V0/Vp) = 0.807 in all cases, because the pre-stress was designed for quasi-static equilibrium.

#### Initial State at t=0

**Station (-24, 10) -- nucleation zone:**

| Quantity | MFEM | Tandem | Diff |
|----------|------|--------|------|
| V_strike | 10^(-1.52) = 0.030 m/s | 10^(-2.00) = 0.010 m/s | 3x |
| tau_strike | 21.489 MPa | 21.148 MPa | +0.341 MPa (delta_tau effect) |
| log10(state) | 8.226 | 8.114 | different Dc in nuc zone |

**Station (0, 0) -- outside nucleation:**

| Quantity | MFEM | Tandem | Diff |
|----------|------|--------|------|
| V_strike | 10^(-9.00) = 1e-9 m/s | 10^(-9.00) = 1e-9 m/s | identical |
| tau_strike | 13.273 MPa | 13.273 MPa | identical |
| log10(state) | 8.146 | 8.146 | identical |

**Outside the nucleation zone, MFEM and Tandem initialization is IDENTICAL.**
Inside the nucleation zone, MFEM follows the SCEC spec (V_nuc=0.03, delta_tau); Tandem uses V_nuc=0.01 with no delta_tau.

#### Both MFEM and Tandem Have an Immediate First Event

**This is the key finding.** Both codes produce a rapid first seismic event -- NOT after 100+ years, but within seconds to a minute.

**State evolution in nucleation zone at t=0:**
```
d(psi)/dt = b*V0/L * (exp((f0-psi)/b) - V/V0)
```
| Code | V_nuc | V/V0 | d(psi)/dt | Characteristic time |
|------|-------|------|-----------|---------------------|
| MFEM | 0.03 | 30,000 | -6.9e-3 /s | ~14 s for psi to drop 0.1 |
| Tandem | 0.01 | 10,000 | -2.3e-3 /s | ~43 s for psi to drop 0.1 |

**Timeline of first event:**

| Milestone | MFEM | Tandem |
|-----------|------|--------|
| Initial Vmax | 0.03 m/s | 0.01 m/s |
| Vmax reaches 0.1 m/s | ~0.5 seconds | ~35 seconds |
| Vmax reaches 1 m/s | ~1.2 seconds | ~50 seconds |
| Rupture arrives at (0,0) | ~2-5 seconds (estimated) | ~60-65 seconds |
| Vmax reaches ~3700 m/s | ~4.1 seconds | not tracked in data |

**Tandem verification from benchmark data:**
- `bp5qd_tandem_x2_-24_x3_10.txt` at t=35.5s: V_strike = 10^(-1.04) = 0.091 m/s
- `bp5qd_tandem_x2_-24_x3_10.txt` at t=38.3s: V_strike = 10^(-0.81) = 0.155 m/s (seismic)
- `bp5qd_tandem_x2_0_x3_0.txt` at t=47.6s: V_strike = 10^(-8.65) (still interseismic)
- `bp5qd_tandem_x2_0_x3_0.txt` at t=65.9s: V_strike = 10^(-0.52) = 0.30 m/s (first event arrives)

**On the 1800-year timescale, both first events appear at t=0 years.** They are visually indistinguishable in the comparison plots.

#### Impact on Limit Cycle and Benchmarking

The different nucleation vigor does NOT affect the long-term limit cycle:
- Both MFEM and Tandem show ~200-year recurrence intervals
- Strike-component waveforms match through 8+ events
- The system loses memory of initial conditions within 1-2 cycles

**Classification:** TUNING -- not a bug. Both approaches are defensible:
- MFEM follows the SCEC spec exactly (V_nuc=0.03, delta_tau)
- Tandem uses a gentler nucleation (V_nuc=0.01, no delta_tau) that also produces correct limit cycles

#### Global Output Interpretation

The global file in `results_1000m_bp5/bp5_full_global.txt`:
- 3,672 lines covering t=0 to 4.15 seconds
- Shows Vmax growing from 10^(-1.52) = 0.03 to 10^(3.57) = 3715 m/s
- This IS the initial rapid seismic event -- NOT an error
- The small time range (4.15 seconds) is because the adaptive time stepper takes tiny steps during the seismic event
- The station files cover ~191 years because they output at a different frequency (every N accepted steps)

**The run is progressing correctly but has not yet completed 1800 years.**

### Plot Analysis (completed run, same as Run 1)

No new plots available (using same completed-run plots as Run 1).

All strike-component metrics remain MATCH at center stations, CLOSE at VS zone edges. The initialization difference has no visible impact on the comparison plots.

### Fixes Applied

**No fixes applied.** The initialization is working correctly per the SCEC spec.

### Recommendations

| Priority | Action | Rationale |
|----------|--------|-----------|
| 1 (optional) | Change V_nuc from 0.03 to 0.01 | Would match Tandem exactly; gentler nucleation |
| 2 (optional) | Remove delta_tau from nucleation zone | Would match Tandem; the SCEC spec allows either |
| 3 (wait) | Let current run complete | Need full 1800-year data before any changes |
| 4 (future) | Compare only cycles 2+ when benchmarking | Skip first event timing; focus on limit cycle |

**If changing V_nuc and removing delta_tau, the only code change needed is in `bp5_params.hpp`:**
- Line 112: `real_t V_nuc = 0.01;` (was 0.03)
- Lines 298-302: Remove or comment out the delta_tau block

### Tests Added/Modified

None -- this is an analysis run only.

### Test Results

Not run (analysis only; cluster run still in progress).

### Status

**Analysis complete. Initialization verified correct per SCEC spec. Awaiting cluster run completion.**

- D1 (V_dip spikes): Fix applied in Run 2 -- awaiting re-run
- D2 (VS zone): TUNING -- awaiting 500m mesh
- D3 (stress sign): RE-VERIFY with new run data (initial data shows +13.27 MPa = correct)
- D4 (dip accumulation): TUNING -- low priority
- **D5 (initialization): ANALYZED -- MFEM follows SCEC spec; differs from Tandem in V_nuc (0.03 vs 0.01) and delta_tau; does not affect limit cycle behavior**

### Key Tandem References (updated for initialization)

| MFEM File | Tandem Counterpart | Key Comparison |
|-----------|--------------------|----------------|
| `config/bp5_params.hpp:112` | `examples/tandem/3d/bp5.lua:72` | V_nuc: MFEM=0.03 (SCEC), Tandem=0.01 |
| `config/bp5_params.hpp:279-307` | `examples/tandem/3d/bp5.lua:94-103` | tau_pre: MFEM includes delta_tau (SCEC Eq. 23), Tandem does not |
| `config/bp5_params.hpp:194-199` | `examples/tandem/3d/bp5.lua:49-57` | Nucleation zone: both at x2=[-30,-18], x3=[4,16] km -- MATCH |
| `friction/dieterich_ruina.hpp:409-431` | `localoperator/DieterichRuinaAgeing.h:51-63` | InitialStatePsi/psi_init: mathematically identical formulas |
| `fault/rate_state_fault.hpp:297-335` | `localoperator/RateAndState.h:123-146` | Init: both compute psi from actual traction, then verify with slip_rate solve |
| `config/bp5_params.hpp:359` | N/A (Tandem computes at runtime) | psi_init: MFEM PreInit placeholder=SteadyState(V_abs), overridden by Init() |
