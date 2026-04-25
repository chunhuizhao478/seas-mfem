# TPV104 V_strike over-shoot — divergence localised to t=0.8s, NOT a wave-operator logic bug

**Date:** 2026-04-24
**Scope:** Localising the V_strike over-shoot (1.7-2.7×) at the propagating-front stations.

## Audit conclusion: wave-operator interior-face flux is logically correct

Read both ADER (`wave_operator.inl:2117-2226`, `:2902-2942`) and Mult
(`wave_operator.inl:1255-1334`, `:1417-1429`) paths plus the
`GodunovFlux::Interior` kernel (`godunov_flux.cpp:350-374`) +
eigendecomposition setup (`godunov_flux.cpp:24-119`).

* **Non-fault interior face**: standard upwind DG.  `F_h = A_n^+ Q_self
  + A_n^- Q_nbr`; `rhs[Elem1] -= F · shape1`, `rhs[Elem2] += F · shape2`.
* **Fault interior face**: per-side imposed-state flux from the Riemann
  solver.  `F_h_plus = A_{can_n}·Q_imp_plus`, `F_h_minus =
  A_{can_n}·Q_imp_minus` via `Interior(can_n, Q, Q)` identity.
  Routing `rhs[+side] -= F_h_plus`, `rhs[-side] += F_h_minus` is
  consistent for both `elem1_on_plus = true` and `false` cases.
* **Shared (MPI-seam) fault face**: same per-side pattern with
  `assemble_sign = elem1_on_plus ? -1 : +1`.
* **`GodunovFlux::Interior` kernel**: textbook `T·A_n^±·T^{-1}·Q`;
  eigendecomposition built from analytic eigenvectors of the elastic
  wave Jacobian; `cp/cs/0` eigenvalue diagonal split into Λ_plus and
  Λ_minus correctly.

**No logic bug found in the wave operator's flux assembly.**

## Time-resolved divergence at the hypocenter (x2_0_x3_7.5)

| t (s) | MFEM V | SeisSol V | ratio | MFEM τ | SeisSol τ | MFEM ψ | SeisSol ψ |
|-------|--------|-----------|-------|--------|-----------|--------|-----------|
| 0.40 | 1.90e-7 | 1.86e-7 | 1.02 | 65.64 | 65.67 | 0.564 | 0.564 |
| 0.50 | 4.64e-5 | 4.57e-5 | 1.01 | 72.24 | 72.26 | 0.564 | 0.564 |
| 0.60 | 2.83e-3 | 2.81e-3 | 1.01 | 77.17 | 77.21 | 0.564 | 0.564 |
| 0.70 | 4.71e-2 | 4.72e-2 | 1.00 | 80.49 | 80.74 | 0.563 | 0.563 |
| 0.80 | **2.72e-1** | **2.58e-1** | **1.06** | **81.54** | **82.84** | 0.554 | 0.555 |
| 0.90 | 2.45e+0 | 1.52e+0 | 1.62 | 70.34 | 82.37 | 0.437 | 0.485 |
| 1.00 | 9.45e+0 | 5.98e+0 | 1.58 | 25.42 | 74.46 | 0.051 | 0.232 |
| 1.10 | 9.05e+0 | 7.32e+0 | 1.24 | 13.62 | 62.75 | -0.046 | 0.078 |
| 1.20 | 7.97e+0 | 6.51e+0 | 1.22 | 12.66 | 55.80 | -0.053 | 0.053 |

* **Through t = 0.7 s**, MFEM and SeisSol agree to 0.3% on τ and 1% on V
  across SEVEN ORDERS OF MAGNITUDE in V (1e-7 → 5e-2 m/s).
* **At t = 0.8 s, the divergence opens**: MFEM's τ drops 1.30 MPa below
  SeisSol's; this 1.6% τ difference is amplified by the friction
  non-linearity (V ≈ V_w_in = 0.1 m/s, dV/dτ steep) into a 5.4% V
  difference.
* **By t = 0.9 s** the over-shoot is 62%; by t = 1.0 s it's 58% — the
  positive-feedback ψ-relaxation cascade is in full effect.

## What this localises

| eliminated by | from |
|---|---|
| **Cadence (R7-001(b))** | dt-halving made things worse; t<0.7 already matches at 1% |
| **Friction-input arithmetic** | τ matches to 0.3% through t=0.7 (FaultFaceFlux is fed correctly) |
| **D2/D3 unit tests** | ComputeTrialTraction is symmetric and correct |
| **D1 rotation pipeline** | wave_operator.inl per-QP rotation is consistent |
| **Wave-operator flux assembly** | interior-face logic is textbook DG, both paths |

What's NOT eliminated:

1. **Numerical dispersion** at MFEM p=1, h=200m vs SeisSol p=4, h≈200m
   (TPV5 mesh).  At p=1 the elastic wave operator has ~6 elements per
   wavelength at the rupture-front frequency (5 Hz), borderline.  At p=4
   the same mesh has equivalent-resolution closer to ~24 modes per
   wavelength.  **Could account for a ~1-2% τ deficit at t=0.8 from
   numerical dispersion alone.**
2. **Mesh asymmetry** discretization residue (already established as
   the σ_n leak source).
3. **Nucleation-accumulator timing**: per-macro-step cumulative
   (MFEM) vs per-sub-step cumulative (SeisSol).  Both telescope to the
   same total but the rate of approach within any macro-step differs
   by O(dt_macro).  At t=0.8, dG/dt ≈ 0.4 → over a 100 µs macro-step
   that's ΔG ≈ 4e-5 → Δτ_2_nuc ≈ 45·4e-5 = 1.8 kPa, **two orders of
   magnitude smaller** than the observed 1.30 MPa gap.  Nucleation
   accumulator timing CANNOT account for it.
4. **Imposed-state velocity-jump amplitude** — the back-flux from the
   fault to bulk Q via `T·A_n·T^{-1}·Q_imp`.  Need to instrument
   directly.

Most likely: combination of (1) and (2).  These are the discretization
residues that would NOT show up at low V (because the slip rate is too
small for them to feed back) but appear sharply once V crosses
~V_w_in = 0.1 m/s (the friction non-linearity).

## Fix candidates

| candidate | mechanism | est. cost | symmetry guarantee |
|---|---|---|---|
| **Increase polynomial order p** to 2 or 3 | reduces O(h^{p+1}) discretization residue ~5×-50× | low (CLI flag); compute cost ~5×-25× | nope; just smaller residue |
| **Run the symmetric-mesh experiment** | rule out (or confirm) (2) — see `tpv104_mesh_asymmetry_finding_2026-04-24.md` | medium (1-2 day mesh work) | exact |
| **Instrument τ_2_trial directly** | dump the bulk-radiation contribution to fault τ at every macro-step at the hypocenter; compare against analytical estimate | small (10-line printf, no [C2] now) | diagnostic |
| **Run with deeper free surface or larger lateral domain** | rule out boundary-reflection contamination | small (CLI / mesh) | n/a |

Recommended order:

1. **τ_2_trial instrumentation** at the hypocenter (cheapest; tells us
   how much the bulk wave operator is delivering).  If MFEM's τ_2_trial
   at t=0.8 is 1.3 MPa less than SeisSol's published τ at t=0.8 minus
   the analytic nucleation contribution at t=0.8 = 82.84 - 43.16 = 39.68
   MPa, the over-shoot is bulk-side; if MFEM's τ_2_trial = 39.68 MPa
   too, the discrepancy is in the friction-side bookkeeping (the
   reported τ includes nucleation + something else).
2. **Increase polynomial order to p=2** as a sanity check.  If the
   over-shoot drops by ~5×, dispersion is the dominant cause.  If
   unchanged, look elsewhere.
3. **Symmetric-mesh experiment** to isolate mesh asymmetry's
   contribution.

## Tooling

* `tpv104/scripts/tpv104_mesh_mirror_check.py` — mesh symmetry verifier.
* `tests/unit/test_tpv104_sigma_n_invariance.cpp` — D2/D3 invariance.
* `wave_operator.inl` D1 instrumentation block (gated by
  `-DSEAS_DIAG_TPV104_FAULT_BASIS`) — per-QP basis dump.
* This file: V_strike divergence point analysis.

## Not-bugs ruled out in this round

- σ_n perturbation (1-4 MPa peaks) explained by mesh-asymmetry
  discretization residue at p=1.  Within SCEC tolerance (10% on σ_n);
  not the V over-shoot driver.
- Per-QP `sign_flipped` mismatch — wave_operator.inl rotation correctly
  normalises `can_n = (0,-1,0)` for ALL QPs.
- FaultFaceFlux::ComputeTrialTraction — algebraically symmetric in σ_n.
- macro-step ψ cadence — dt-halving made things worse; pre-divergence
  agreement at 1% rules cadence out as the source of the divergence.
- Wave-operator flux assembly logic — textbook DG, no error.
