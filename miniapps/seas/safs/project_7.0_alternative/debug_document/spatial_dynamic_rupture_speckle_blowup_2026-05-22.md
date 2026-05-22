# Spatial dynamic rupture — on-patch "speckle" blow-up (job 7744292, 2026-05-22)

Second SAFS dynamic-rupture blow-up, **distinct** from the 2026-05-20
north-tip sliver event (`spatial_dynamic_rupture_postevent_debug_2026-05-20.md`).
Run: `spatial_dyn_smoke` job 7744292, 400 ranks, config
`config/spatial_friction_slip_weakening_safs_projected_stress.toml`
(psi=37° / az351 H&Z stress, 500 m `z0embed` sliver-free mesh, μ_s=0.65 /
μ_d=0.30 / D_c=1.0, gradual_overstress δτ=6.5 MPa).

---

## Symptom

`slip_rate_strike` on the fault shows a **smooth Gaussian nucleation patch
(correct)** surrounded by **salt-and-pepper** — isolated, grid-scale DOFs
with huge slip rate, NOT a coherent propagating front (ParaView, t≈0.84 s).
V_max trace:

| step | t (s) | V_max (m/s) | regime |
|---|---|---|---|
| 1000–2300 | 0.35–0.81 | ~3 | healthy nucleation patch slipping |
| 2400 | 0.84 | 975 | blow-up onset |
| 2400–12900 | 0.84–4.51 | 975 → 5.7e145 | exponential growth ~3.2%/step |
| 13000–13200 | 4.55–4.62 | 5.4e136 → 1.5e89 → **3.17** | **recovery** |

Key features distinguishing this from the 2026-05-20 event:
1. **Early** (t≈0.8 s, during the nucleation regime) — not post-event (t≈3 s).
2. **On / next to the nucleation patch** — not at the far north tip.
3. **Recovers** to the physical ~3 m/s — the 05-20 event went to NaN.
4. Mesh is the new sliver-free `500m z0embed`, so it is **not** the slivers.

The recovery (finite the whole way, 1e145 → 3, never NaN) rules out a global
CFL / scheme instability (those never recover). It is a **localized,
self-limiting** instability: once the speckled DOFs slip past D_c = 1.0 m,
μ pins to μ_d, the slip-weakening gradient that feeds the growth is
exhausted, the corrupted stress radiates away, and the genuine rupture front
(~3 m/s, elsewhere) becomes the global max again.

---

## Evidence — the dynamic-rupture CODE is correct (TPV205 comparison)

Traced both drivers. The fault machinery is the **same code path** for the
proven planar TPV205 driver and the SAFS spatial driver:

`Tpv205SubStepIterator::AdvanceWithSubStepStates` → `StepOneQP_` →
`FaultFaceFlux::ComputeTrialTraction` / `SolveLSW_TPV205` /
`BuildImposedState`, with shared-fault QPs through
`FaultFaceFlux::EvaluateADER_LSW`, corrector via `WaveOperator::AdvanceADER`.

- `SolveLSW_TPV205` (tpv205_friction.hpp:93) is the standard radiation-damped
  balance `V = max(0, (|τ|−μσ_n)/η_s)`, V ∥ τ — stable, no sign error.
- Slip accumulation `d.slip{1,2} += V*dt_sub` (tpv205_substep_iterator.cpp:122)
  is correct and single-owner (the iterator; `EvaluateADER_LSW` deliberately
  does NOT re-accumulate — fault_face_flux.cpp:798).
- CFL is reduced by the identical `cfl/(3·(2N+1))` factor as TPV205
  (spatial_dyn_driver.cpp:1218; the job-7743554 lesson). dt is NOT the issue.
- Interior flux is plain homogeneous Godunov (`--no-sidecar-material`,
  λ=μ=32 GPa); the Phase-R bi-material path is inactive.
- Initial stress is smooth in ParaView before rupture → the per-DOF stress
  projection (`FaultGeometry::ComputeSAFSParams`) and the gradual_overstress
  resolver are NOT injecting salt-and-pepper into the IC.

**No bug found in the dynamic-rupture code.** TPV205 is stable on the same
code because of three protections SAFS removed:

| Protection | TPV205 (stable) | SAFS run 7744292 (blows up) |
|---|---|---|
| Cohesive-zone resolution | well-resolved | `L_nuc/h_min` min = **4.79** (criterion ≥ 10) |
| Strength barrier outside rupture area | `mu_s_barrier`=1e4 → `V≡0` | none on the active fault (only z∈[−20,−15] km) |
| Background criticality | uniform, subcritical | `\|τ_pre\|/(μ_s·σ_n_eff)` outside = **0.986** |

---

## Root cause

A **two-factor** numerical instability (same class as the 05-20 doc, but
moved early and onto the patch by the new psi37 calibration):

1. **SEED** — the cohesive zone is under-resolved (`L_nuc/h_min` min 4.79 ≪ 10;
   ~2.5 elems in the patch). An under-resolved slip-weakening front emits
   grid-scale (salt-and-pepper) slip-rate oscillations of a few %.
2. **AMPLIFIER** — the background is 1.4% below static yield (outside ratio
   0.986) with no barrier to clamp it. The few-% grid noise tips near-critical
   background DOFs over yield; slip-weakening positive feedback then grows
   each one exponentially (~3.2%/step) until D_c saturation halts it.

The `[derived]` well-posedness gate (`spatial_print_derived.cpp`) gates only
on `outside_max ≥ 1.0` (0.986 squeaked by) and `L_nuc/h_min` is **printed,
never gated** — so the run launched straight into the known blow-up regime.

---

## Fix (applied) — lower background criticality, keep spontaneous rupture

Per run directive (chose: lower criticality; keep μ_d=0.30 strong stress drop;
ensure spontaneous rupture). Edited
`config/spatial_friction_slip_weakening_safs_projected_stress.toml`:

| param | old | new | why |
|---|---|---|---|
| `mu_s_default` | 0.65 | **0.85** | outside ratio 0.641/0.85 ≈ **0.75** (25% below yield) → removes the amplifier |
| `mu_d_default` | 0.30 | 0.30 | keep strong +14.4 MPa dynamic stress drop (propagation) |
| `delta_tau_strike_pa` | 6.5e6 | **20.0e6** | reach the higher static yield AND build a supercritical core > L_nuc |
| `d_c_default` | 1.0 | 1.0 | unchanged |

**Spontaneous-rupture verification** (patch σ_n_eff≈49.3 MPa, |τ_pre|≈29.2 MPa,
μ_bulk=32 GPa, Gaussian `F(r)=exp(−(r/radius)²)`, radius=3000 m):

- L_nuc = 32e9·1.0/((0.85−0.30)·49.3e6) = **1.18 km** (was 1.86 km; smaller → easier to exceed).
- yield gap at center = μ_s·σ_n_eff − |τ_pre| = 41.9 − 29.2 = 12.7 MPa.
- supercritical core radius r_c = 3000·√(−ln(12.7/20)) = **2.02 km**.
- margin **r_c/L_nuc ≈ 1.71** (> 1 required) → ruptures spontaneously, does not arrest at T_nuc.
- center overshoot = 29.2 + 20 − 41.9 = **+7.3 MPa** ≥ 0 → `[derived]` TRIGGER gate passes.
- in-patch dynamic ratio = 29.2/(0.30·49.3) = **1.97 > 1** → STRESS-DROP gate passes.
- δτ=20 MPa is below the 25 MPa "gross overdrive" ceiling.

All four `[derived]` gates green by hand; the gate aborts at startup
(`abort_on_failure=true`) if any prediction is off, so a bad recalibration
fails fast rather than wasting a multi-hour run.

---

## Next-run checklist (`--print-derived`)

Expect (within mesh/projection rounding):
- `max OUTSIDE ratio` ≈ **0.75** (no longer the 0.986 WARNING).
- `L_nuc / h_min` min ≈ **3.0** (max ≈ 9.7) — *worse* than before (see risk).
- `nucleation overshoot` ≈ **+7 MPa (sufficient)**.
- `in-patch min |τ_pre|/(μ_d·σ_n_eff)` ≈ **1.97 (> 1)**.
- `in-patch max |τ_pre|/(μ_s·σ_n_eff)` ≈ **0.70 (< 1, locked)**.

## Known remaining risk + next lever
This fix removes the **amplifier** (criticality), not the **seed**
(resolution) — raising μ_s shrank L_nuc, so `L_nuc/h_min` is now ~3.0
(worse). The bet is that a 25%-subcritical background can no longer grow the
grid noise. **If residual front speckle reappears** (watch `slip_rate_strike`
near the front in ParaView, and V_max for any super-linear growth), the next
lever is the **250 m z0embed mesh** (`[mesh].path`, ~2× finer →
`L_nuc/h_min` ~6–19), which attacks the seed directly. Pairing both is the
robust combination.

## Status
- CONFIRMED (code read): fault machinery identical to proven TPV205; LSW
  solve, slip accumulation, CFL all correct. Root cause = under-resolution
  (seed) + near-critical unbarriered background (amplifier).
- APPLIED (config): μ_s 0.65→0.85, δτ_strike 6.5→20 MPa; spontaneous-rupture
  margin r_c/L_nuc ≈ 1.71 verified by hand.
- PENDING (re-run): Frontera `--print-derived` to confirm the gate numbers,
  then a full run to confirm the speckle is gone and the rupture propagates.
