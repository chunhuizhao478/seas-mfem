# Building Mixed (Central) Flux + Runge–Kutta for Dynamic Rupture

**Date:** 2026-05-28
**Scope:** SEAS-MFEM dynamic-rupture miniapp (`miniapps/seas/dynamic/`, driver
`drivers/spatial_dyn_driver.cpp`).
**Companion / motivation:** `mixed_flux_dev/drdg3d_mixed_flux_comparison_2026-04-28.{md,pdf}`
(the root-cause analysis showing why mixed/central flux blows up under ADER at p=1).
**Reference implementation:** DRDG3D (`/Users/chunhuizhao/projects/drdg3d`), which runs
central flux + RK54 stably at CFL ≈ 0.3.

---

## 0. Executive summary

The central-flux instability at p=1 is **not** fixable inside ADER: at spatial order
p, the Cauchy–Kovalevskaya predictor is bounded to time order `O ≤ p+1`
(`wave_operator.hpp:421-424`, recursion `D(k+1) = −Σ_d A_d ∂_{x_d} D(k)`, element-local),
so p=1 is locked to ADER-O2, whose stability function `R(z)=1+z+z²/2` has
`|R(iy)|² = 1+y⁴/4 > 1` for every imaginary mode the central flux creates. The only
way to keep central flux at p=1 is an integrator whose stability region **contains a
segment of the imaginary axis** — i.e. an explicit Runge–Kutta of order ≥ 3.

**The good news: most of the machinery already exists.** This is an assembly job, not a
from-scratch build:

| Capability | Status | Location |
|---|---|---|
| Method-of-lines RHS `dQ/dt = M⁻¹(−Face+Vol)` | **exists** | `WaveOperator::Mult`, `wave_operator.inl:885` |
| Mixed-flux dispatch (central on `central_flux_face_set_`) **in the Mult path** | **exists** | `ComputeFaceFluxRHS:2939`, `ComputeSharedFaceFluxRHS:3505` |
| Instantaneous fault friction in the Mult path | **exists** | `fault_flux_->Evaluate`, `wave_operator.inl:2715` |
| `Evaluate` is ψ-stateless (RK-stage-safe by design) | **exists** | `fault_face_flux.cpp:328`, guard `:432` |
| Instantaneous dψ/dt for SAFS slip-law-SRW | **exists** | `SlipLawSRWPsi::Rate_SRW`, `slip_law_srw_psi.hpp:260` |
| Classical RK4 idiom on `wave.Mult` | **exists** (unit test) | `DoRK4Step`, `test_rk4_conservation.cpp:67` |
| **Coupled RK4-on-(Q,ψ) stepper** | **existed, removed** | commit `8461c67` (TPV102 v9.2.0), replaced by ADER in `9d79554` |
| Coupled RK stepper wired into a production driver | **missing** | — to add to `spatial_dyn_driver.cpp` |
| Nucleation evaluated at RK stage times (absolute, not telescoped) | **missing** | — adapt `spatial_nucleation.cpp` |
| CFL factors calibrated for RK + central flux | **missing** | `ComputeMaxDt:5638` factors are "interim placeholders" tuned for ADER |

**Recommended approach:** a **hand-written, fixed-step coupled RK stepper** that
integrates the augmented state `(Q, ψ, slip)` together — porting the proven v9.2.0
`8461c67` RK4 loop, swapping the aging-law ψ-rate for `Rate_SRW`, evaluating nucleation
at stage times, and calling `wave.SetMixedFluxMode(MixedFluxMode::Adjacent)` before the
loop. A generic `mfem::ODESolver` (e.g. `DormandPrinceRK45`) cannot be used directly
because the fault state `(ψ, slip)` lives in `DOFData`, not in the MFEM state `Vector`.

---

## 1. Background: why RK, and why it must be order ≥ 3

(Full derivation in the companion analysis doc; condensed here.)

1. **Central flux ⇒ imaginary eigenvalues.** `GodunovFlux::Central` is
   `F* = ½·(A⁺+A⁻)·(Q_L+Q_R)` (`godunov_flux.cpp:382-423`) — zero dissipation. Its
   header records the identity `Interior − Central = +½·|A_n|·(Q_self−Q_nbr)`
   (`godunov_flux.hpp:64-68`): central is upwind with the entire dissipation term
   deleted. The fault-adjacent modes it touches become purely imaginary, `λ = iω`.

2. **ADER-O2 amplifies them unconditionally.** `R(z)=1+z+z²/2 ⇒ |R(iy)|²=1+y⁴/4>1`
   for all `y>0`. No Δt stabilizes a purely imaginary mode.

3. **At p=1 you cannot raise the ADER order.** The predictor's k-th Taylor term needs
   the k-th element-local spatial derivative (`wave_operator.hpp:421-424`); a degree-1
   field has none beyond first order, so `D(2)=D(3)=0` identically and
   `--ader-order 3/4` collapse bit-for-bit to O2. `O ≤ p+1`: p=1→O2, p=2→O3, p=3→O4
   (`order` is `MFEM_VERIFY`'d to `{2,3,4}`, `wave_operator.inl:5398`).

4. **RK decouples time order from p.** A method-of-lines RK is order-q in time
   regardless of spatial p. Required: its stability region must contain a segment of
   the imaginary axis.
   - Forward Euler / RK2: `R(z)=1+z+z²/2` (RK2 ≡ ADER-O2) — **fails** (no imaginary
     coverage).
   - RK3: `|R(iy)|≤1` for `y≤√3≈1.73` (conditional).
   - **RK4 (classical): `y ≤ 2√2 ≈ 2.83`.**
   - **RK54 (Kennedy–Carpenter 5-stage low-storage, DRDG3D's choice): `y ≲ 2.8–3.4`.**
   - Dormand–Prince RK45 (adaptive 7-stage 4(5)): similar RK4-class imaginary coverage.

   ⇒ Use RK order ≥ 3; RK4 or RK54 are the natural picks.

---

## 2. What already exists (inventory with code)

### 2.1 `WaveOperator::Mult` — the RK entry point is complete

`wave_operator.inl:884-967` (it `override`s `TimeDependentOperator::Mult`,
`wave_operator.hpp:155`):

```cpp
void WaveOperator<MeshType>::Mult(const Vector &Q, Vector &dQdt) const
{
   dQdt.SetSize(height);  dQdt = 0.0;
   ComputeVolumeRHS(Q, dQdt);            // strong-form volume term
   ComputeFaceFluxRHS(Q, dQdt);          // interior + boundary + FAULT faces
   if constexpr (IsParallelMesh<MeshType>::value)
      ComputeSharedFaceFluxRHS(Q, dQdt); // MPI seam faces (incl. shared fault)
   if (pml_layer_) ApplyPMLDamping(Q, dQdt);
   ApplyMassInverse(dQdt);               // dQ/dt = M⁻¹(...)
}
```

This is exactly what an RK stage calls. It already assembles everything an RK stage
needs.

### 2.2 The Mult path **already does mixed-flux dispatch**

`ComputeFaceFluxRHS`, interior non-fault faces, `wave_operator.inl:2939-2947`:

```cpp
if (mf_on && central_flux_face_set_.count(f) > 0)
   flux_.Central(nor, Q_self, Q_nbr, F_h);   // zero-dissipation, on adjacent faces
else
   flux_.Interior(nor, Q_self, Q_nbr, F_h);  // upwind elsewhere
```

Same predicate on the parallel seam (`ComputeSharedFaceFluxRHS:3505-3512`). `mf_on` is
`mf_on_`, set by `SetMixedFluxMode(...)`, which also (re)builds `central_flux_face_set_`
via `BuildCentralFluxFaceSet_` (`wave_operator.inl:1504`). **So enabling mixed flux on
the RK path is a single call** — `wave.SetMixedFluxMode(MixedFluxMode::Adjacent)` — with
no change to the flux code. (The identical predicate exists on the ADER path at `:4486`
/ `:5254`; the only difference there is it operates on the time-integrated state `I`.)

### 2.3 The Mult path **already handles the fault face** (instantaneously)

`ComputeFaceFluxRHS` detects the fault attribute on interior faces
(`wave_operator.inl:2597`) and runs the per-side friction Riemann solve
(`:2715-2718`):

```cpp
fault_flux_->Evaluate(fdata, Q_plus_local, Q_minus_local,
                      Q_imp_plus, Q_imp_minus);   // instantaneous solve at Q
```

each side then fluxes its imposed state (`flux_.Interior(can_n, Q_imp_*, ...)`,
`:2801`). The shared-fault counterpart calls the same `Evaluate` at `:3355`. **This is
the per-stage friction kernel RK needs** — it solves for slip-rate `V` and the imposed
traction at the *instantaneous* stage state, which is precisely the method-of-lines
semantics.

### 2.4 `Evaluate` is ψ-stateless **by design, for RK staging**

`fault_face_flux.cpp:328`, with the revealing invariant comment at `:333-339` and the
debug guard at `:432`:

```cpp
// v9.2.0 F01+F02 invariant (REVIEW R-V92-H07): the driver's coupled
// RK4 on (Q, psi) ... reads `data.psi` at entry to this call as the
// stage-local psi and relies on it being UNCHANGED on return. A future
// refactor that adds a psi write here would silently invalidate the
// RK4 stage arithmetic (the driver would double-integrate psi).
...
MFEM_ASSERT(data.psi == psi_at_entry, "FaultFaceFlux::Evaluate mutated data.psi ...");
```

`Evaluate` → `WriteBackState` (`:310-323`) writes `data.{slip_rate,V1,V2,tau1_corr,
tau2_corr,sigma_n_corr}` but **never** `data.psi`, `data.slip1`, or `data.slip2`. The
state variables are integrated *by the stepper*. This is the seam that makes a coupled
RK on `(Q, ψ, slip)` clean.

### 2.5 SAFS state-evolution rate exists

`SlipLawSRWPsi` (`slip_law_srw_psi.hpp`) implements both:
- analytic update `UpdateStateAnalyticSlipLawSRW(...)` (`:92`) — used by the ADER substep
  iterator; and
- **instantaneous rate** `Rate(V,psi,L)` (`:185`) and the SRW-aware
  `Rate_SRW(V,psi,L,V_w,a,b,V0,f0,muW)` (`:260`) — `dψ/dt` for RK staging.

⚠️ Production-mode guard (R-001, header `:16-21`): the bare `Rate(V,psi,L)` aborts in
production; the coupled RK loop must call `Rate_SRW(...)` with the per-QP `V_w`
(directly analogous to the template's `aging_law.Rate(...)`).

### 2.6 The proven template: coupled RK4-on-(Q,ψ), commit `8461c67`

TPV102 ran a hand-written classical RK4 that integrated `(Q, ψ)` together before being
replaced by ADER. Recover it with `git show 8461c67:./drivers/tpv102_driver.cpp`
(loop ~lines 783–965). Annotated skeleton:

```cpp
Vector k1, k2, k3, k4, Q_tmp;                  // bulk stage rates
std::vector<real_t> psi_n, sr_k1..sr_k4, psi_k1..psi_k4;
AgingLawPsi aging_law(...);                    // SAFS: use SlipLawSRWPsi + Rate_SRW

for (int step = 0; step < nsteps; step++) {
   for (i) psi_n[i] = dof_data[i].psi;         // save ψ at step start

   // STAGE 1 @ t,  ψ = ψ_n
   ApplyNucleation(dof_data, ..., t);          // nucleation at STAGE TIME (absolute)
   wave.Mult(Q, k1);                           // does mixed flux + friction solve;
                                               //   writes V_i = dof_data[i].V1/V2/slip_rate
   for (i) {
      sr_k1[i] = dof_data[i].slip_rate;
      psi_k1[i] = aging_law.Rate(sr_k1[i], psi_n[i], Dc);   // dψ/dt at stage 1
      dof_data[i].psi = psi_n[i] + 0.5*dt*psi_k1[i];        // stage-2 ψ input
   }
   // STAGE 2 @ t+dt/2,  ψ = ψ_n + (dt/2)ψ_k1
   ApplyNucleation(dof_data, ..., t + dt/2);
   add(Q, dt/2, k1, Q_tmp);  wave.Mult(Q_tmp, k2);
   for (i) { sr_k2..; psi_k2[i]=Rate(sr_k2,dof_data[i].psi,Dc);
             dof_data[i].psi = psi_n[i] + 0.5*dt*psi_k2[i]; }
   // STAGE 3 @ t+dt/2 (same time → no new nucleation call)
   add(Q, dt/2, k2, Q_tmp);  wave.Mult(Q_tmp, k3);
   for (i) { sr_k3..; psi_k3[i]=Rate(sr_k3,dof_data[i].psi,Dc);
             dof_data[i].psi = psi_n[i] + dt*psi_k3[i]; }
   // STAGE 4 @ t+dt
   ApplyNucleation(dof_data, ..., t + dt);
   add(Q, dt, k3, Q_tmp);  wave.Mult(Q_tmp, k4);
   for (i) { sr_k4..; psi_k4[i]=Rate(sr_k4,dof_data[i].psi,Dc); }

   // COMBINE — Butcher (1,2,2,1)/6 weights, applied to BOTH Q and ψ:
   for (n) Q[n]        += dt/6*(k1[n]+2*k2[n]+2*k3[n]+k4[n]);
   for (i) dof_data[i].psi = psi_n[i] + dt/6*(psi_k1+2*psi_k2+2*psi_k3+psi_k4);
   //   slip += V_avg*dt with the same (1,2,2,1)/6 weights  ⇔ RK4 integral of V
   t += dt;
}
```

The key structural facts this template establishes (all directly reusable):
- **ψ is a first-class RK state**, integrated with the same Butcher weights as `Q`
  (O(dt⁴) coupling), not operator-split.
- **`wave.Mult` writes the per-stage `V` into `DOFData`**, which the loop reads to form
  `dψ/dt` and to accumulate slip.
- **Nucleation is applied at each stage time** (`t`, `t+dt/2`, `t+dt`), as an *absolute*
  value — not the telescoped per-substep increment the ADER path uses.
- **slip accumulation** uses the RK4-weighted `V_avg` (since `dslip/dt = V`).

---

## 3. What does NOT carry over from ADER (and why that's fine)

The companion analysis worried that the friction coupling is "welded to the ADER
predictor." It is — but that machinery is a **parallel path the RK route never touches.**
The ADER-only pieces are:

| ADER-only machinery | File | RK replacement |
|---|---|---|
| Time-averaged friction `Q̄ = I/dt` then `I_imp = Q_imp·dt` | `EvaluateADER`, `fault_face_flux.cpp:651` | none — use instantaneous `Evaluate` (already in `Mult`) |
| Cauchy–Kovalevskaya substep states `Q(τ)=Σ τ^k/k! D(k)` | `ComputeADERSubStepStates:1302` | RK stage states `Q + Σ a_{ij}k_j` |
| Side-channel `I_imp` pointers (`substep_I_imp_*_flat_`) | `wave_operator.hpp:836` | none — fault flux computed live inside each `Mult` |
| Gauss–Lobatto substep nodes/weights `deltaT_`, `time_weights_` | `tpv104_substep_iterator.cpp` | RK Butcher tableau |
| Substep iterators (`Tpv104SubStepIterator`, `friction_iterator.hpp`) | `dynamic/*_substep_iterator.cpp` | the coupled RK loop (§2.6) |
| Telescoped per-substep nucleation increment | `spatial_nucleation.cpp:135` | absolute `SmoothStep(t_stage)` at stage time |

None of the six "ADER coupling assumptions" apply because **RK uses the instantaneous
`Evaluate`, not `EvaluateADER`.** The instantaneous path is older, simpler, and already
on the critical path of `Mult`. We are not fighting the ADER coupling — we are using the
*other* fault path that already exists.

> One inherited caveat: the ADER shared-fault path (np>1) deliberately bypasses its own
> side-channel because of a frame-reconciliation issue (R-1600/R-1601,
> `wave_operator.inl:4863-4920`). The **Mult** shared-fault path uses `Evaluate`
> directly (`:3355`), so it is *not* subject to that ADER-specific workaround — but
> multi-rank fault correctness on the RK path must still be tested explicitly (§6).

---

## 4. Design: coupled explicit RK on the augmented state `(Q, ψ, slip)`

### 4.1 State partition

| Component | Lives in | Integrated by |
|---|---|---|
| Bulk field `Q` (σ, v) | MFEM `Vector` | `wave.Mult` + RK weights |
| State variable `ψ` (per fault QP) | `DOFData[i].psi` | `Rate_SRW` + RK weights |
| Slip `slip1/slip2` (per fault QP) | `DOFData[i].slip1/2` | `V_avg·dt` (RK weights) |

Because `ψ`/`slip` are **not** in the `Vector`, a stock `mfem::ODESolver` cannot drive
the whole system — hence the hand-written coupled loop. (Alternative, §5.0 option B:
pack `ψ`/`slip` into the `Vector` to enable a stock adaptive solver; larger change.)

### 4.2 Per-stage recipe (the contract every stage obeys)

For RK stage `i` at stage time `t + cᵢ·dt` with stage bulk state `Q⁽ⁱ⁾`:

1. **Nucleation:** write the *absolute* overstress at `t + cᵢ·dt` into
   `DOFData.tau1_nuc/tau2_nuc` (not the telescoped increment).
2. **Stage ψ:** ensure `DOFData[i].psi` holds the stage value `ψ⁽ⁱ⁾` (set at the end of
   the previous stage block).
3. **`wave.Mult(Q⁽ⁱ⁾, kᵢ)`** — internally: mixed-flux on `central_flux_face_set_`,
   instantaneous friction `Evaluate` at `(Q⁽ⁱ⁾, ψ⁽ⁱ⁾)` → writes `V⁽ⁱ⁾`, imposed
   traction; `kᵢ = dQ/dt`.
4. **Capture rates:** `sr_kᵢ = DOFData.slip_rate`; `psi_kᵢ = Rate_SRW(V⁽ⁱ⁾, ψ⁽ⁱ⁾, …)`.
5. **Advance ψ to next stage input** per the Butcher tableau.

Final combine: `Q`, `ψ`, and `slip` each updated with the same tableau weights.

### 4.3 Mixed flux

One line before the time loop:

```cpp
wave.SetMixedFluxMode(MixedFluxMode::Adjacent);   // central on fault-adjacent faces
```

`Adjacent` (Zhang Fig 4b) is the right mode (central only near the fault); `AllContinuous`
(central everywhere) is far more aggressive and needs CFL ≈ 0.4× (§5.4). Match the
existing `--mixed-flux` CLI string handling already present for the ADER path.

---

## 5. Phased implementation plan

### Phase 0 — Choose the RK scheme  *(decision required)*

| Scheme | Stages | Imag-axis | Pros | Cons |
|---|---|---|---|---|
| **Classical RK4** (recommended first cut) | 4 | y≤2.83 | template exists (`8461c67`, `DoRK4Step`); simplest | fixed step; 4 RHS/step |
| **RK54** (Kennedy–Carpenter, low-storage) | 5 | y≲2.8–3.4 | matches DRDG3D exactly; low memory | new coefficients to add |
| **Dormand–Prince RK45** (adaptive) | 7 | RK4-class | adaptive dt; solver exists (`seas_driver.cpp:440`) | needs ψ/slip in the `Vector` (option B) or a custom embedded-error loop |

**Recommendation:** start with **classical RK4** (reuse `8461c67`), validate, then switch
coefficients to **RK54** to track DRDG3D. Defer adaptive DP45 unless step-size control is
needed — for a CFL-bound hyperbolic system the adaptivity is largely wasted.

### Phase 1 — Bulk-only coupled stepper (no mixed flux yet)
Port the `8461c67` RK4 loop into `spatial_dyn_driver.cpp` as an alternative to
`AdvanceADERWithSubStep_Spatial` (the time loop is at `:2214-2250`; the ADER advance at
`:2244`). Add a CLI selector (e.g. `--time-integrator {ader,rk4,rk54}`, default `ader`).
Keep `mixed_flux=none`. **Acceptance:** on a linear (frictionless) problem, RK4 and
ADER-O2 agree to O(dt²) (the `DoRK4Step` conservation test, `test_rk4_conservation.cpp`,
is the template).

### Phase 2 — Wire SAFS fault state into the RK stages
Replace the template's `aging_law.Rate(...)` with `SlipLawSRWPsi::Rate_SRW(V, ψ, Dc,
V_w[i], a, b, V0, f0, muW)` (`slip_law_srw_psi.hpp:260`); source per-QP `V_w`, `a`, `b`
exactly as the substep iterator does (`tpv104_substep_iterator.cpp:864`). Accumulate slip
with RK4-weighted `V_avg`. **Acceptance:** with `mixed_flux=none`, the RK4 run reproduces
the ADER-O2 upwind SAFS run (job 7751974 baseline: V_max ≈ 10.5 m/s) within integration-
order tolerance.

### Phase 3 — Nucleation at stage times
The current `gradual_overstress` resolver is a *telescoped per-substep increment*
(`ApplyGradualOverstressIncrement`, `spatial_nucleation.cpp:135`;
`SmoothStepIncrement = SmoothStep(t) − SmoothStep(t−dt)`, `:31`). Under RK this would
**double-apply** across stages that revisit overlapping sub-intervals. Add an *absolute*
form, `ApplyGradualOverstressAbsolute(dof_data, params, T_nuc, t_stage)`, that writes
`tau1_nuc = SmoothStep(t_stage)·F(r)·Δτ_dip` (and strike) directly — mirroring TPV102's
`ApplyNucleation(..., t_stage)` in the template. Call it at each stage time.
**Acceptance:** summed nucleation forcing over `[0, T_nuc]` matches the ADER telescoped
total to round-off.

### Phase 4 — Enable mixed flux + recalibrate CFL
Add `wave.SetMixedFluxMode(MixedFluxMode::Adjacent)`. Recalibrate `dt` (§5.4).
**Acceptance:** the central-flux SAFS run (the configuration that ran away under ADER,
job 7751975) is now **stable** at the RK CFL, and shows the mixed-flux accuracy benefit
the user measured at p=1.

### Phase 5 — Verification & regression (§6).

### 5.0 Note — Option B (stock adaptive solver)
If adaptive RK45 is wanted: extend the evolved `Vector` to append `[ψ, slip1, slip2]` per
fault QP, have `Mult` read/write those entries (instead of `DOFData`), and hand the whole
thing to `mfem::DormandPrinceRK45`. This is the only way to reuse MFEM's adaptive control,
but it touches the fault-flux data plumbing broadly and is **not** recommended for the
first cut.

### 5.4 CFL / stability calibration

Current `dt` (`ComputeMaxDt`, `wave_operator.inl:5638-5723`):

```
dt = cfl_mixed_flux_factor · cfl · h / c_p          c_p = √((λ+2μ)/ρ),  h = inscribed diam
cfl_mixed_flux_factor =  None 1.0 | Adjacent 0.9 | AllContinuous 0.4
```

plus a driver-side order de-rating `cfl / (3·(2N+1))`, **N = mesh order p**
(`spatial_dyn_driver.cpp:1434`).

⚠️ Both knobs are ADER-calibrated. The `0.9/0.4` factors are flagged in-source as
"interim placeholders pending a multi-step stability calibration" (`:5648-5655`), and the
`1/(3(2N+1))` de-rating is a TPV205-derived **ADER** safety conversion. For RK:
- Replace the `1/(3(2N+1))` de-rating with the RK stability bound: keep
  `max|λ|·dt < y_max` with `y_max ≈ 2.83` (RK4) or `≈ 2.8` (RK54, conservative).
- DRDG3D's empirical mixed-flux number is **CFL ≈ 0.3** — use as the starting target for
  `Adjacent`.
- Recalibrate `cfl_mixed_flux_factor` against the RK stability region (the ratio of
  central faces still matters; `Adjacent` ≈ near-upwind, `AllContinuous` needs the
  tighter ≈0.3–0.4).

**ψ explicit-RK stability caveat** (from the `8461c67` review, R-V92-H06): explicit RK on
the state law is only *conditionally* stable, with local timescale `~ Dc / (b·V0·exp(...))`.
For a dynamic event the wave CFL dt is far below that timescale (TPV102: dt/τ ~ 3e-3), so
it is safe — but keep the `ψ` tripwire from the template, and for slip-law-SRW verify the
margin on the SAFS envelope (low-`ψ`, high slip-rate). If a future cycle simulation drives
`ψ` into the stiff regime, the state update must go implicit/analytic.

---

## 6. Verification & regression

1. **Linear conservation (Phase 1):** RK4 vs ADER-O2 on frictionless bulk —
   `test_rk4_conservation.cpp` pattern; energy drift bounded.
2. **Upwind equivalence (Phase 2):** RK + `mixed_flux=none` reproduces the ADER-O2 upwind
   SAFS baseline (job 7751974) within time-integration tolerance.
3. **Stability with central flux (Phase 4):** the previously-divergent central-flux
   config (job 7751975) runs to `t_final` with bounded `V_max`.
4. **Mixed-flux benefit (Phase 4):** confirm the p=1 accuracy improvement the user
   measured, now stable.
5. **MPI fault correctness:** multi-rank run with a fault crossing a rank seam — confirm
   the Mult shared-fault path (`Evaluate` at `:3355`) is conservation-consistent across
   the seam (the ADER R-1600 frame issue does **not** automatically transfer, but must be
   re-checked on the RK path).
6. **Nucleation total (Phase 3):** integrated overstress matches the ADER telescoped sum.
7. **Byte-exact contract:** none of these touch the TPV*/BP5 ADER paths; the existing
   regression contract (`make test`) must stay green — the RK path is additive behind a
   CLI selector.

---

## 7. Open questions / risks

- **Slip-law-SRW explicit-RK stability margin** on the full SAFS envelope (not just the
  TPV102 numbers) — quantify before production. Mitigation: keep the analytic
  `UpdateStateAnalyticSlipLawSRW` as an operator-split fallback for stiff QPs.
- **RK54 coefficients** — confirm the exact Kennedy–Carpenter low-storage tableau to
  match DRDG3D bit-for-bit (DRDG3D `mod_para.F90` RK54 default), if exact cross-comparison
  is desired.
- **`Adjacent` vs `AllContinuous` for the measured p=1 benefit** — if the benefit comes
  from reduced dissipation broadly (not just near the fault), `AllContinuous` may be
  wanted, at the tighter CFL. Test both.
- **Output cadence** — the substep `slip_rate_substep_max` diagnostic
  (`spatial_dyn_driver.cpp:2230`) is ADER-substep-specific; define the RK-stage analogue
  (per-stage max, or final-stage value).
- **Per-stage MPI cost** — RK does N `Mult` calls/step (each with its face-neighbour
  exchange) vs ADER's predictor-once; profile to confirm the extra exchanges are
  acceptable.

---

## 8. References

- **Root-cause analysis:** `mixed_flux_dev/drdg3d_mixed_flux_comparison_2026-04-28.{md,pdf}`.
- **Reference code:** DRDG3D `/Users/chunhuizhao/projects/drdg3d` (central flux + RK54 @ CFL≈0.3).
- **Proven template:** `git show 8461c67:./drivers/tpv102_driver.cpp` (coupled RK4-on-(Q,ψ)).
- **Key source:**
  - `dynamic/wave_operator.inl` — `Mult:885`, mixed-flux dispatch `:2939/:3505`,
    fault `Evaluate` `:2715/:3355`, `ComputeMaxDt:5638`, `BuildCentralFluxFaceSet_:1504`.
  - `dynamic/fault_face_flux.cpp` — `Evaluate:328`, `WriteBackState:310`, ψ-guard `:432`.
  - `dynamic/godunov_flux.cpp` — `Central:382`, `Interior:350`.
  - `friction/slip_law_srw_psi.hpp` — `Rate:185`, `Rate_SRW:260`,
    `UpdateStateAnalyticSlipLawSRW:92`.
  - `dynamic/spatial_nucleation.cpp` — `ApplyGradualOverstressIncrement:135`,
    `SmoothStepIncrement:31`.
  - `drivers/spatial_dyn_driver.cpp` — time loop `:2214`, ADER advance `:429/:2244`,
    CFL `:1422`.
  - `tests/unit/test_rk4_conservation.cpp` — `DoRK4Step:67`.
```
