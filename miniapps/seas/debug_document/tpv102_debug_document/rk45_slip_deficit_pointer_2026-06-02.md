# TPV102 RK45 strike-slip deficit — see the unified report

Date: 2026-06-02

The TPV102 (and TPV104) RK45 / mixed-flux spatial-driver strike-slip deficit is
diagnosed in the unified report:

  **`../tpv104_debug_document/tpv102_tpv104_rk45_slip_deficit_findings_2026-06-02.md`**

TPV102 one-liner: at the *identical* p1 + mixed-flux-adjacent discretization,
**ADER-O2 matches the SCEC benchmark to −1.5 %** but the **RK45 run is −16 %**
(growing −6 % at the hypocenter → −28 % at the corners). Root cause: the RK45
macro-step is **3× the validated ADER-O2 step** (effective CFL 0.15 vs 0.05) and
the RK path does not sub-cycle the stiff rate-state friction → the rupture front
is progressively over-damped. **Not** a friction/config/σ_n/sign bug (ADER-O2
with the same config matches). No code changed; recommendations + a CFL-sweep
experiment are in the unified report §7.
