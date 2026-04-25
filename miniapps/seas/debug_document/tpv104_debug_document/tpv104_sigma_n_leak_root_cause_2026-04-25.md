# TPV104 σ_n perturbation root-cause finding — stress-channel leak

Date: 2026-04-25
Job: Frontera 7677329 (`tpv104_normaltrace_dev.sbatch`)
Probe: C-1n NORMAL in `dynamic/fault_face_flux.cpp:271-303`
Probed QP: rank 3, DOF 68, fault coord (-26.6 m, 0, -7550 m)
  — closest QP to nominal hypocenter (0, 0, -7500 m) on a 200 m mesh; offset is intra-element, expected.

## Bottom line

The ±0.5 MPa σ_n_trial perturbation observed at the hypocenter during rupture transit is **dominated by the stress-average channel of the trial-traction formula, not the velocity-jump channel**. Specifically, the bulk wave operator's compressional stress field (SXX, SYY, SZZ in the fault-local frame) acquires non-zero values at fault-adjacent QPs that should remain ~0 for a pure strike-slip rupture in the TPV104 spec.

## Trial-traction formula recap

```
sigma_n_trial = etaP * [(Q_minus[VX] − Q_plus[VX])              ← v_jump_term
                       + (Q_plus[SXX]/Zp + Q_minus[SXX]/Zp_neig)] ← stress_term
```

(VX, SXX in this comment are fault-LOCAL components: VX = fault-normal velocity, SXX = fault-normal stress in the per-QP local basis.)

## Frontera run 7677329 — raw probe output

| Quantity | Value |
|---|---|
| Total `[C-1n NORMAL]` lines | 7024 (= 2.0 s ÷ ~2.85e-4 s, expected) |
| Tagging line | `[diag] rank 3 tagging hypo DOF 68 at (-26.6, 0.0, -7550.0)` (single line, MPI_MINLOC selected exactly one rank, expected) |
| Peak \|σ_n_trial\| | **5.279 × 10⁵ Pa** ≈ 0.53 MPa |
| Peak v_jump_term | +7.156 × 10⁴ Pa |
| Peak stress_term | **−5.994 × 10⁵ Pa** |
| Peak occurred at | line 3766 of 7024 → t ≈ 1.07 s (mid-rupture transit) |

### Channel decomposition at peak |σ_n_trial|

| Channel | Magnitude | Fraction of total |
|---|---|---|
| v_jump_term | 7.16 × 10⁴ Pa | 14% |
| **stress_term** | **5.99 × 10⁵ Pa (opposite sign)** | **113%** |

The two channels partially cancel; the net σ_n_trial = stress_term + v_jump_term ≈ −5.99e5 + 0.72e5 ≈ −5.28e5 Pa. The stress channel is the dominant carrier of the leak; v_jump is small but nonzero.

### Late-time samples (last 3 macro steps)

```
sigma_n_trial=-1.6717e+05  v_jump=+2.0220e+03  stress=-1.6919e+05  | Q_p[VX]=-8.02e-3  Q_m[VX]=-7.77e-3  Q_p[SXX]=-1.60e+05  Q_m[SXX]=-1.78e+05
sigma_n_trial=-1.6714e+05  v_jump=+2.0476e+03  stress=-1.6919e+05  | Q_p[VX]=-8.17e-3  Q_m[VX]=-7.92e-3  Q_p[SXX]=-1.60e+05  Q_m[SXX]=-1.78e+05
sigma_n_trial=-1.6711e+05  v_jump=+2.0543e+03  stress=-1.6916e+05  | Q_p[VX]=-8.29e-3  Q_m[VX]=-8.03e-3  Q_p[SXX]=-1.60e+05  Q_m[SXX]=-1.78e+05
```

Two key observations from the late-time data:

1. **Q_*[SXX] are not converging to zero post-rupture.** Both sides hold ≈ −1.7 × 10⁵ Pa steadily, indicating a persistent compressional residue in the fault-adjacent bulk, not a transient.
2. **Q_plus[SXX] and Q_minus[SXX] are asymmetric** by ~1.8 × 10⁴ Pa (-1.60e5 vs -1.78e5). This asymmetry, weighted by etaP/Zp, is exactly what stress_term accumulates as the σ_n leak.
3. **Q_plus[VX] vs Q_minus[VX] are almost equal** (−8.02e-3 vs −7.77e-3, diff ~2.5e-4). The fault opening rate is small relative to the slip rate (~7-8 m/s in this regime, so V_normal/V_strike ratio is O(10⁻⁴) — physically tiny).

### Early-time samples (first 3 macro steps)

```
step 1:  sigma_n_trial=+0  v_jump=0  stress=0   | Q_p[VX]=0  Q_m[VX]=0  Q_p[SXX]=0  Q_m[SXX]=0
step 2:  sigma_n_trial=−1.79e-29  v_jump=−1.62e-29  stress=−1.72e-30  | Q_p[VX]=+1.68e-36  Q_m[VX]=−3.37e-37  Q_p[SXX]=−3.43e-30  Q_m[SXX]=−2.38e-44
step 3:  sigma_n_trial=+1.60e-13  v_jump=+1.24e-13  stress=+3.64e-14  | Q_p[VX]=−7.91e-21  Q_m[VX]=+7.56e-21  Q_p[SXX]=+4.75e-14  Q_m[SXX]=+2.53e-14
```

Step 1: identically zero (initial Q = 0 + nucleation hasn't started).
Step 2-3: Q_*[SXX] are populated by O(10⁻³⁰) → O(10⁻¹⁴) **before** Q_*[VX] grow comparably. So compressional stress radiates into the fault-adjacent bulk *first*, not as a consequence of fault opening.

## What this rules out

| Hypothesis | Status | Reason |
|---|---|---|
| Basis-rotation bug (per-QP sign-flip leaking strike→normal) | **Ruled out as primary** | Such a bug would primarily corrupt the velocity channel; v_jump is only 14% of the swing. The earlier per-QP audit ([tpv104_review_round?_2026-04-24.md], 52/48 sign-flip bimodal) verified can_n=(0,−1,0) on velocity rotation. The stress-rank-2 rotation was *not* audited and is now the prime suspect. |
| Mesh non-mirror-symmetry | **Ruled out as primary** | Mesh check (0/215143 tets have mirror partner) would cause ±-symmetric leak in both channels, not a stress-channel-only leak with sign that follows the rupture front. |
| Nucleation pre-stress contamination | **Ruled out** | `data.sigma_n_nuc` is added to `sigma_n_fric` *after* `Evaluate` returns the trial; the probe sees the trial value, which has no nucleation contribution. The fact that Q_*[SXX] is non-zero while in the loading phase (steps 2-3, before t=0.5 s nucleation onset) confirms this. |
| ψ cadence (macro-step vs per-sub-step ψ update) | **Ruled out as primary** | Earlier dt-halving experiment (job 7677079) made the σ_n perturbation *worse*, not better. ψ cadence affects the friction solve, not the trial traction. |

## What this points to: stress-rotation in the fault frame

The trial-traction probe shows compressional stress (Q_*[SXX] in fault-local frame) acquiring values of order 10⁵ Pa at fault-adjacent QPs, with a persistent ±-asymmetry between the two sides. For a pure strike-slip rupture on a y=0 plane with normal=(0,−1,0), the bulk should radiate predominantly shear waves (SXY in fault-local), not compressional. The presence of SXX — and, more diagnostically, the asymmetry between sides — implies the stress field that the wave operator carries is being seeded with a fault-normal compressional component on every macro step.

The most likely upstream source: the **stress rotation in the fault face flux's imposed state**.

The Riemann path computes a fault-frame imposed σ_n (from `data.sigma_n_fric`) and rotates it back to global to construct the imposed Q* state. For a pure strike-slip fault on y=0 with normal=(0,−1,0), this rotation should produce only σ_yy in the global frame. **However, if the per-QP basis rotation matrix has any off-diagonal contamination — even at the level of floating-point sign-flip bimodality — the imposed σ_n will leak into σ_xx (and σ_zz), and that σ_xx propagates into the bulk via the wave operator's RHS injection at the fault face**.

Cross-check: the velocity rotation was audited and shown to be clean (can_n=(0,−1,0) for all 100 of 100 QPs after sign normalization). But:

- The velocity rotation is rank-1 (vector): a per-QP sign-flip in (can_n, can_t1, can_t2) self-cancels because v is bilinear in basis.
- The stress rotation is rank-2 (tensor): σ_global = R · σ_local · R^T. A per-QP sign-flip in *only one* axis of (can_n, can_t1, can_t2) does *not* self-cancel in a rank-2 rotation; it produces an off-diagonal sign error.

This asymmetry between rank-1 and rank-2 rotations explains why velocity diagnostics (V_strike, V_dip) look clean while stress diagnostics (σ_n) leak.

## Why this matches the observed station data

- **σ_n perturbation magnitude**: 0.5 MPa peak in stations matches 0.53 MPa peak |σ_n_trial| in probe.
- **σ_n perturbation timing**: the 1.07 s peak in the probe matches the rupture-front passage time at the hypocenter.
- **V_dip leakage** (3 cm/s peak) and **slip_dip drift** (3.5 mm by t=2s): these would be downstream consequences of σ_n perturbation perturbing the friction equilibrium, which couples Theta = sqrt(tau1² + tau2²) and produces small dip-component imbalance.
- **V_strike clean (0.35% post-fw-fix)**: the velocity rotation is correct, so strike-slip propagates normally.

## Probe overhead / cost

- 7024 probe lines × ~250 chars = ~1.8 MB stderr, well within Frontera limits.
- Single rank emits, no MPI, no deadlock at any rank count.
- Probe was strict zero outside `#ifdef SEAS_DIAG_FAULT_FLUX` (verified: no-flag binary contains 0 diag strings).

## Build-system fixes that enabled this finding

(For posterity — these blocked the previous probe attempt 7677305.)

1. `Makefile:2613, 2617`: added `$(SEAS_EXTRA_CPPFLAGS)` to TPV104_DRIVER_OBJ and TPV104_SUBSTEP_ITERATOR_OBJ rules. Without this, the driver TU was silently built without `-DSEAS_DIAG_FAULT_FLUX` even when other TUs picked it up, leading to undefined `g_seas_my_rank` symbol at link time.
2. `jobs/tpv104/tpv104_normaltrace_dev.sbatch:138`: changed `MFEM_CPPFLAGS=-DSEAS_DIAG_FAULT_FLUX` → `SEAS_EXTRA_CPPFLAGS=-DSEAS_DIAG_FAULT_FLUX`. The Makefile reads `SEAS_EXTRA_CPPFLAGS` (Makefile:284); MFEM_CPPFLAGS is silently ignored.
3. `drivers/tpv104_driver.cpp:661-703`: added the hypocenter-tagging block (mirrored from tpv102_driver.cpp:540-588). Without this, no DOF had `diag_print = true` and the probe could not fire even if compiled in. MPI_Allreduce(MPI_DOUBLE_INT, MPI_MINLOC) selects the rank whose closest fault QP wins; only that rank sets `diag_print = true` on exactly one DOF.

Commit: `8aac5fa` on `feature/elasticity-inertia`.

## Recommended next step

Probe target: **stress rotation in the fault face flux imposed state**.

Concrete action:
1. Add a probe inside `dynamic/godunov_flux.cpp` (or wherever `BuildImposedState` constructs Q_imposed for the fault Riemann path), printing at the hypocenter QP:
   - The fault-frame σ_n (input to rotation)
   - The global-frame σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz (output of rotation)
   - The per-QP rotation matrix R (3×3)
2. Run the same `tpv104_normaltrace_dev.sbatch` shape (no parameter changes), tfinal=2.0 s, 400 ranks.
3. Expected signature if hypothesis is right:
   - σ_yy carries the imposed σ_n cleanly
   - σ_xx is non-zero by the same magnitude as Q_*[SXX] swing (~10⁵ Pa)
   - The rotation matrix R has off-diagonal entries near machine zero **except** for some QPs where one entry is O(1) — those are the sign-flipped QPs.

If the hypothesis is confirmed, the fix is to enforce a consistent stress rotation that does not depend on per-QP basis sign — likely by replacing the on-the-fly rotation with a precomputed R once at init (after sign-normalization), and reusing it across all macro steps and all QPs of each face.

## Files / artifacts referenced

- Probe source: `dynamic/fault_face_flux.cpp:256-305` (commit 8aac5fa).
- Probe driver block: `drivers/tpv104_driver.cpp:661-703` (commit 8aac5fa).
- Sbatch: `jobs/tpv104/tpv104_normaltrace_dev.sbatch` (commit 8aac5fa).
- Raw probe data on Frontera: `tpv104_normaltrace_7677329.err` in submit dir, ~1.8 MB.
- Earlier per-QP velocity-rotation audit: see `tpv104_review_round[3-7]_2026-04-24.md` (sign_flipped bimodality).
- Mesh-symmetry check: `tpv104_mesh_asymmetry_finding_2026-04-24.md`.
- ψ cadence review: `tpv104_psi_cadence_review_2026-04-24.md`.
- f_w fix that resolved V_strike over-shoot but left σ_n leak: commit `4d05df2`.
