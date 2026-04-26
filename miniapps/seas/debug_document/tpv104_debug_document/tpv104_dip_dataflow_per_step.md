# TPV104 Dip-Direction Per-Step Data Flow (2026-04-25)

Purpose: trace every read and write of dip-direction quantities in one TPV104 macro-step. Used by Round-3 R-201 Direction 1 to localize where the np=1 baseline `slip_dip = −1.0078e-02 m` is born — so we can put falsifiable probes at each stage.

## What "dip" means in this codebase

| name | direction (TPV104, vertical y=0 fault, ref_normal=(0,−1,0), up=(0,0,1)) | fault-local Q channel |
|---|---|---|
| `tangent1 = can_t1 = dip` | `(0, 0, −1)` (down-dip) | `VY` (velocity), `SXY` (shear stress) |
| `tangent2 = can_t2 = strike` | `(+1, 0, 0)` (along strike) | `VZ` (velocity), `SXZ` (shear stress) |
| `normal = can_n` | `(0, −1, 0)` | `VX`, `SXX` |

So in the **fault-local frame** that `FaultFaceFlux::Evaluate` operates on:
- `Q[VY]` = velocity in the **dip** direction
- `Q[SXY]` = shear stress in the **dip** direction
- `tau1_*` (in `DOFData` / `EvalStageState`) = **dip** shear traction
- `V1` = **dip** slip rate
- `slip1` = **dip** accumulated slip

The naming `tau1`/`V1`/`slip1` is BP5's `t1=dip, t2=strike`. The fault-local slot `VY/SXY` is the rotated Q's "first tangent" component, which corresponds to dip after `Tinv_can · Q_global`.

For pure strike-slip (TPV104 expected): `V1 ≡ slip1 ≡ tau1_corr ≡ 0`. Any non-zero value is the pollution we are tracking.

---

## Macro-step skeleton (`drivers/tpv104_driver.cpp:1397+`)

For each `step in [0, nsteps)`:

```cpp
// (M.0) Save psi at step start
for (i in num_fault_total)   psi_n[i] = dof_data[i].psi;

// (M.1) Nucleation increment (driver-side, per-DOF, no MPI)
if (!disable_nucleation && num_fault_total > 0) {
    ApplyNucleationIncremental_TPV104(dof_data, fault_coords,
                                      t + dt_step, dt_step);
}

// (M.2) One ADER predictor-corrector step.  Reads and writes Q;
//       internally invokes the fault Riemann at every fault QP and
//       writes back to dof_data[i].V1/V2/slip_rate/tau1_corr/tau2_corr.
wave.AdvanceADER(Q, dt_step, ader_order, Q_new);
Q.Swap(Q_new);
t += dt_step;

// (M.3) Per-DOF state-variable update (analytic SRW)
for (i in num_fault_total) {
    dof_data[i].psi = UpdateStateAnalyticSlipLawSRW(psi_n[i],
                                                    dof_data[i].slip_rate,
                                                    dof_data[i].Dc, dt_step,
                                                    V_w[i], dof_data[i].a,
                                                    b, V0, f0, f_w);
    // (M.4) Per-DOF slip integration
    dof_data[i].slip1 += dof_data[i].V1 * dt_step;     // <-- DIP slip update
    dof_data[i].slip2 += dof_data[i].V2 * dt_step;
}

// (M.5) V_max reduction (no dip-specific math)
// (M.6) Station / ParaView output (read-only consumer of dof_data.slip1)
```

The dip-channel is read or written at exactly five sites: M.1 (no — TPV104 sets only `tau2_nuc`), M.2 (the bulk story — see §A below), M.3 (no — psi update consumes `slip_rate=|V|` only), M.4 (the slip1 integrator), M.6 (read-only).

So per macro-step, the dip-channel quantities `tau1_corr` / `V1` / `slip1` are written exclusively inside `wave.AdvanceADER`. Localizing the asymmetry reduces to localizing where inside `AdvanceADER` `V1 ≠ −V1_mirror` first appears.

---

## §A — `wave.AdvanceADER` internal flow (`wave_operator.inl:3415+`)

```
AdvanceADER(Q, dt, order, Q_new):
   I = ComputeADERTimeIntegrated(Q, dt, order)                    // (A.1) predictor
   rhs = 0
   ComputeADERVolumeUpdate(I, rhs)                                 // (A.2) bulk volume
   ComputeADERFaceFluxRHS(I, dt, rhs)                              // (A.3) interior faces
   if parallel:
       ComputeADERSharedFaceFluxRHS(I, dt, rhs)                    // (A.4) shared faces
   if pml_layer_:
       apply PML damping                                            // (A.5) (n/a TPV104)
   ApplyMassInverse(rhs)                                           // (A.6)
   Q_new = Q + rhs                                                 // (A.7)
```

### Dip-channel touches per stage

| stage | reads dip-related Q channels (VY, SXY) | writes dip outputs | MPI | notes |
|---|---|---|---|---|
| A.1 ComputeADERTimeIntegrated | yes — `D_curr = Q`, then `D_next = -A_d · ∂_d D_curr` recursion (3 directions × order−1 sweeps) | writes `I[:]` (time-integrated state); `I[VY]`, `I[SXY]` carry dip | NONE — fully element-local | predictor; no mass-inverse coupling between elements |
| A.2 ComputeADERVolumeUpdate | reads `I[VY]`, `I[SXY]`, ... at QPs | adds `(K_d · I)` to `rhs[*]` | NONE | element-local volume integral via `ComputeVolumeRHS` |
| A.3 ComputeADERFaceFluxRHS (non-fault else-branch) | reads `I_self`, `I_nbr` from local elements; calls `flux_.Interior(nor, I_self, I_nbr, F_h)` | adds `−w · shape · F_h` to `rhs` on both elements | NONE | the standard Godunov flux `F_h` mixes all 9 channels via the Jacobian |
| A.3 ComputeADERFaceFluxRHS (fault branch) | reads `I_self`, `I_nbr` at fault QPs; rotates to canonical frame; calls `fault_flux_->EvaluateADER` | per-side imposed state → `F_h_plus / F_h_minus` → `rhs[..dof_offset1..]` and `rhs[..dof_offset2..]`; ALSO writes `dof_data[i].V1`, `tau1_corr`, `slip_rate` (via `EvaluateADER → Evaluate → WriteBackState`) | NONE for interior fault faces | the only stage that writes `dof_data` per QP |
| A.4 ComputeADERSharedFaceFluxRHS | per-component `q_gf.ExchangeFaceNbrData()` × NUM_STATE; reads `nbr_data[c]` at shared-face QPs | writes `rhs[..dof_offset1..]`; for shared-fault faces also writes `dof_data` (via `EvaluateADER`) | YES — ghost exchange | the only MPI stage in the time step |
| A.5 PML | n/a | n/a | NONE | TPV104 doesn't use PML |
| A.6 ApplyMassInverse | reads `rhs[..]` | writes `rhs[..] := M⁻¹ · rhs` | NONE | element-local block-diagonal solve |
| A.7 `Q_new = Q + rhs` | reads Q, rhs | writes Q_new | NONE | trivial |

After A.7, `Q_new` becomes Q at the next macro-step. The dip-channel components of Q are `Q[c=VY,SXY,...]`. They are written ONCE per macro-step here.

---

## §B — Fault Riemann internals at one fault QP (`fault_face_flux.cpp`, called from A.3 / A.4)

This is where TPV104's per-QP dip update happens. Sequence at QP `q` of fault face `f` on rank R:

### B.0 — Frame reconstruction (caller, `wave_operator.inl:1264–1273`)

```cpp
// Per QP at fault face:
//   qpd = fault_basis_->GetBasis(fb_idx).qp_data[q]
//   — qpd.normal, tangent1, tangent2 already include the FaultBasis Step-5
//     sign_flipped negation (set at constructor time).
real_t can_n[3], can_t1[3], can_t2[3];
for (d in 0..3):
   can_n[d]  = qpd.sign_flipped ? -qpd.normal[d]   : qpd.normal[d];
   can_t1[d] = qpd.sign_flipped ? -qpd.tangent1[d] : qpd.tangent1[d];
   can_t2[d] = qpd.sign_flipped ? -qpd.tangent2[d] : qpd.tangent2[d];
```

Result for TPV104 (any QP, both +y and −y mirror partner faces):
```
can_n  = (0, −1, 0)   (= ref_normal)
can_t1 = (0, 0, −1)   (= dip, down)
can_t2 = (+1, 0, 0)   (= strike)
```

(The `qpd.sign_flipped` negation undoes Step 5 of `FaultBasis::ComputeOrientedFrame`, recovering the canonical frame that's bit-identical on both sides of any face.)

### B.1 — Side label

```cpp
// post-R-101: per-face, geometric, bit-equal to the old !sign_flipped
//             flag on this axis-aligned mesh
const bool elem1_on_plus = interior_fault_elem1_on_plus_[fb_idx];
```

### B.2 — Build rotation, rotate Q to fault-local

```cpp
GodunovFlux::BuildRotation(can_n, can_t1, can_t2, T_can);
GodunovFlux::BuildRotationInverse(can_n, can_t1, can_t2, Tinv_can);
// I_self_can[c] = sum_k Tinv_can(c, k) * I_self[k]   (and same for I_nbr)
// Pick:
//   I_plus_local  = elem1_on_plus ? I_self_can : I_nbr_can
//   I_minus_local = elem1_on_plus ? I_nbr_can  : I_self_can
```

After this, `I_plus_local[VY] = velocity in dip direction on + side` and `I_plus_local[SXY] = shear traction in dip direction on + side`.

### B.3 — Time-average for ADER

```cpp
// fault_face_flux.cpp:567–573
Q_avg_plus[c]  = I_plus_local[c]  / dt;
Q_avg_minus[c] = I_minus_local[c] / dt;
```

### B.4 — `Evaluate` pipeline (the actual physics, `fault_face_flux.cpp::Evaluate` → 4 stages)

#### B.4.a `ComputeStageState → ComputeTrialTraction` (line 40–66)

Reads `Q_avg_plus[VY], Q_avg_minus[VY], Q_avg_plus[SXY], Q_avg_minus[SXY]`:

```cpp
// Eq. (7b)  fault-local "tangent1" channel = dip
s.tau1_trial = data.eta_s * (Q_minus[VY] - Q_plus[VY]
                             + Q_plus[SXY]  / Zs_plus
                             + Q_minus[SXY] / Zs_minus);
```

This is the **first place** in the per-step flow where `tau1_trial` is nonzero. For a true y-mirror Q field with `Q_plus = mirror(Q_minus)`, the dip-channel components on the two sides should satisfy `Q_minus[VY] = +Q_plus[VY]` (mirror) and `Q_minus[SXY] = +Q_plus[SXY]` (mirror), giving:
```
s.tau1_trial = eta_s · (Q_plus[VY] − Q_plus[VY] + Q_plus[SXY]·(1/Zs_p + 1/Zs_m))
             = eta_s · 2·Q_plus[SXY]/Zs        (if Zp=Zs equal both sides)
```

NOTE: this is NOT zero even for a y-mirror state — `tau1_trial = 2·η_s·SXY_plus/Zs` — but it is the SAME for both members of any y-mirror pair of fault faces. So if the fault face flux is y-mirror invariant, both partners produce equal `tau1_trial` → equal `V1` → equal contribution. The dip pollution requires a y-asymmetry in `Q[VY]` or `Q[SXY]` *between* +y and −y mirror faces.

For TPV104 the IC has `Q[SXY] ≡ 0` (only `tau1_0=0`, `tau2_0=29.38 MPa` are nonzero, and the static prestress lives in DOFData, not bulk Q in fluctuation mode). So at t=0 step 1: `Q ≡ 0` everywhere, `s.tau1_trial = 0`, `V1 = 0`, `slip1 += 0`. Slip_dip starts at 0. 

#### B.4.b `CompleteFromTrial` (line 124–139)

```cpp
s.tau1_total = data.tau1_0 + data.tau1_nuc + s.tau1_trial;
s.tau2_total = data.tau2_0 + data.tau2_nuc + s.tau2_trial;
s.Theta = sqrt(tau1_total² + tau2_total²);
```

For TPV104: `tau1_0 = 0` (pure strike-slip; CLAUDE.md "Fault-local tangent frame"), `tau1_nuc = 0` (TPV104 nucleation only writes `tau2_nuc`, see `tpv104_nucleation.hpp:134`). So `tau1_total = s.tau1_trial`. Cleanly inherits the trial value.

#### B.4.c `CompleteFromTheta → solver_.Solve` (Brent on `|V|`) (line 141–152)

Reads `Theta = sqrt(tau1_total² + tau2_total²)`, `psi`, `sigma_n_total`, `eta_s`, `a`. Returns `V_abs = |V|`. **No dip-specific dispatch** — Brent operates on the magnitude. The output `V_abs` couples back to dip via:

#### B.4.d `CompleteFromVabs` (line 155–184)

```cpp
// Eq. (9): slip-rate decomposition
strength = |sigma_n_total| · a · asinh(V_abs · C);
s.V1 = V_abs · tau1_total / (strength + eta_s · V_abs);     // <-- DIP slip rate
s.V2 = V_abs · tau2_total / (strength + eta_s · V_abs);

// Eq. (10): corrected traction
s.tau1_corr = s.tau1_trial - eta_s · s.V1;                 // <-- DIP corrected traction
s.tau2_corr = s.tau2_trial - eta_s · s.V2;
```

`V1` is computed as `V_abs · (tau1_total / Theta)` modulo the friction-strength denominator. For a y-mirror invariant input (`tau1_total` same on +y and −y mirror partners), `V1` is the same on both partners. **The dip pollution can only arise here if `tau1_total` differs between mirror partners — i.e., if `Q[VY]` or `Q[SXY]` from the bulk side differs in a non-y-mirror way.**

#### B.4.e `BuildImposedState` (line 186–223)

```cpp
// Eq. (11)/(12) imposed states.  Both sides see s.tau1_corr / s.tau2_corr.
Q_imp_minus[VY] = Q_minus[VY] - (1/Zs_m) · (s.tau1_corr - Q_minus[SXY]);
Q_imp_plus[VY]  = Q_plus[VY]  + (1/Zs_p) · (s.tau1_corr - Q_plus[SXY]);
Q_imp_minus[SXY] = s.tau1_corr;     Q_imp_plus[SXY] = s.tau1_corr;
```

Same `s.tau1_corr` on both sides. Asymmetric Q-input → symmetric tau1_corr output → asymmetric Q_imp[VY] (because the velocity adjustment depends on side-specific Q[VY], Q[SXY] state).

#### B.4.f `WriteBackState` (line 225–238)

```cpp
data.slip_rate    = s.V_abs;
data.V1           = s.V1;                                      // <-- DOFData write
data.V2           = s.V2;
data.tau1_corr    = data.tau1_0 + data.tau1_nuc + s.tau1_corr;  // stored as TOTAL
data.tau2_corr    = data.tau2_0 + data.tau2_nuc + s.tau2_corr;
data.sigma_n_corr = data.sigma_n0 + data.sigma_n_nuc + s.sigma_n_corr;
```

This is the only `dof_data.V1` write per step. The next M.4 step reads it as `dof_data[i].slip1 += dof_data[i].V1 * dt_step`.

### B.5 — Rotate imposed state back to global, accumulate flux

```cpp
// I_imp = Q_imp · dt   (post-EvaluateADER scaling)
// I_imp_plus_g[c] = sum_k T_can(c, k) · I_imp_plus[k]
// I_imp_minus_g[c] = sum_k T_can(c, k) · I_imp_minus[k]
flux_.Interior(can_n, I_imp_plus_g,  I_imp_plus_g,  F_h_plus);
flux_.Interior(can_n, I_imp_minus_g, I_imp_minus_g, F_h_minus);

// Per-side DG assembly (post-R-103-retraction: this is verified bit-symmetric
//                       in IEEE-754; both branches arithmetically identical)
if (elem1_on_plus) {
    rhs[Elem1] -= w · shape1 · F_h_plus;
    rhs[Elem2] += w · shape2 · F_h_minus;
} else {
    rhs[Elem1] += w · shape1 · F_h_minus;
    rhs[Elem2] -= w · shape2 · F_h_plus;
}
```

`F_h_plus` and `F_h_minus` are 9-vectors in the global frame; their dip-channel content is in `F_h[VY]` and `F_h[SXY]`.

---

## §C — Where the dip channel can be polluted

Compositing the trace, the dip-direction asymmetry can ENTER the per-step flow at exactly these read sites:

| site | reads | y-mirror sensitivity |
|---|---|---|
| C.1 IC: `Q[*]` at t=0 | written by `InitializeState_TPV104` | should be `≡ 0` for fluctuation mode; IC asymmetry would persist |
| C.2 `data.tau1_0`, `data.tau1_nuc`, `data.sigma_n0` per fault QP | written by `InitializeFaultDOFs_TPV104` and `ApplyNucleationIncremental_TPV104` | `tau1_0 = 0`, `tau1_nuc = 0` per spec; if either is non-zero on +y vs −y mirror QP, polluted |
| C.3 `data.psi` per fault QP | written by `M.3 UpdateStateAnalyticSlipLawSRW` (consumes `slip_rate = |V|`) | psi propagates the `|V|` magnitude only; mirror-symmetric for y-mirror |
| C.4 `Q_avg_plus[VY], Q_avg_minus[VY], Q_avg_plus[SXY], Q_avg_minus[SXY]` at fault QP | bulk Q at the QP, rotated by `Tinv_can` | the FIRST load-bearing dip read per step; the source of `s.tau1_trial` |
| C.5 `s.tau1_total = data.tau1_0 + data.tau1_nuc + s.tau1_trial` | (C.2) and (C.4) | inherits asymmetry from C.4 (since C.2 should be 0) |
| C.6 `s.V1 = V_abs · tau1_total / (...)` at fault QP | (C.5) and `V_abs` from Brent | Brent's `V_abs` is mirror-symmetric (depends on `Theta`, `psi`, `sigma_n`, all mirror-symmetric); `V1` inherits sign and magnitude of `tau1_total` |
| C.7 `data.V1` per fault QP, written by `WriteBackState` | (C.6) | the ONLY per-step dip writeback to DOFData |
| C.8 `slip1 += V1 · dt` | (C.7) | accumulator in M.4 |

**Conclusion of the trace.** The dip-direction pollution at the hypocenter slip_dip output is the time-integral of `V1` at the hypocenter QP. `V1` is determined ENTIRELY by `tau1_total = tau1_0 + tau1_nuc + tau1_trial` (and by `V_abs`, which is mirror-symmetric). Per the TPV104 spec, `tau1_0 = tau1_nuc = 0`. So `V1 ≠ 0` requires `tau1_trial ≠ 0`. And `tau1_trial` (B.4.a) is determined entirely by `Q_avg_plus[VY], Q_avg_minus[VY], Q_avg_plus[SXY], Q_avg_minus[SXY]` at the fault QP.

**The intrinsic baseline asymmetry must be born in one of:**
1. Bulk Q channels `VY` and/or `SXY` becoming y-asymmetric *between* mirror-paired fault QPs *during one ADER step* (i.e., A.1 predictor + A.2 volume + A.3 non-fault face flux + A.4 shared-face flux + A.6 mass inverse fail to commute with y-mirror at FP precision).
2. `data.tau1_0` or `data.tau1_nuc` actually being non-zero at some QPs (config bug).
3. The IC `Q` actually being non-zero in the dip channels (config bug).

The next probe (Direction 1 in REVIEW.md) should cleanly separate these. Specific probe points:

| probe | code site | env-gate | falsification target |
|---|---|---|---|
| **P-1** | end of `InitializeFaultDOFs_TPV104` | `SEAS_TPV104_AUDIT_INIT=1` | dump `(qp_idx, x, y, z, tau1_0, tau1_nuc, sigma_n0, psi)`; assert `tau1_0 == 0` and `tau1_nuc == 0` at every QP — rules out C.2 |
| **P-2** | end of `InitializeState_TPV104` | `SEAS_TPV104_AUDIT_IC=1` | for every DOF i, assert `Q[VY, i]·Q[VY, mirror(i)] ≥ 0` and `|Q[VY, i] − Q[VY, mirror(i)]| < 1e−14` — rules out C.1 |
| **P-3** | start of `Evaluate` (after rotation, before `ComputeTrialTraction`) | `SEAS_TPV104_AUDIT_INPUT=1` | for every QP `q` and its y-mirror partner `q'`, dump `(Q_avg_plus[VY], Q_avg_minus[VY], Q_avg_plus[SXY], Q_avg_minus[SXY])` and check mirror — first hit identifies whether asymmetry exists at QP read time |
| **P-4** | start of step N=10 (post-rupture-front), at one fault QP at hypocenter and its mirror | as above | by step ~10 the wave has propagated 0.5–1 cell radius; if P-2 passes (IC clean) but P-3 at step 10 fails, the asymmetry is accumulating in bulk Q during AdvanceADER |

If P-3 fails at step 1 with `Q ≡ 0` IC: there is something in `AdvanceADER` that produces non-zero asymmetric `Q[VY]` or `Q[SXY]` from a zero IC — that is a clear bug in A.1 / A.2 / A.3 / A.6. A `assert(rhs == 0)` probe inside each stage of A.1–A.6 with zero input identifies the responsible kernel.

If P-3 passes at step 1 but the slip_dip eventually grows: the asymmetry comes from accumulated bulk-side noise as the rupture propagates. Then the next probe target is the `flux_.Interior` non-fault Godunov flux's y-mirror commutativity at a non-fault face adjacent to the fault.

If P-3 fails at step 1 because the IC is not bit-mirror: P-2 should have caught it; revise `InitializeState_TPV104`.

---

## Quick-reference: which file, which line

| stage | file | line range |
|---|---|---|
| M.1 nucleation | `dynamic/tpv104_nucleation.hpp` | 104–138 |
| M.2 AdvanceADER | `dynamic/wave_operator.inl` | 3415+ |
| M.3 psi update | `drivers/tpv104_driver.cpp` | 1488–1500 |
| M.4 slip integration | `drivers/tpv104_driver.cpp` | 1501–1502 |
| A.1 ComputeADERTimeIntegrated | `dynamic/wave_operator.inl` | 910+ |
| A.1.1 ApplySpatialDerivative | `dynamic/wave_operator.inl` | 756–855 |
| A.2 ComputeADERVolumeUpdate | `dynamic/wave_operator.inl` | 972+ |
| A.3 ComputeADERFaceFluxRHS (interior fault) | `dynamic/wave_operator.inl` | 2304–2470 |
| A.4 ComputeADERSharedFaceFluxRHS | `dynamic/wave_operator.inl` | 2470+ |
| A.6 ApplyMassInverse | `dynamic/wave_operator.inl` | (search "ApplyMassInverse") |
| B.4.a ComputeTrialTraction | `dynamic/fault_face_flux.cpp` | 40–66 |
| B.4.b CompleteFromTrial | `dynamic/fault_face_flux.cpp` | 124–139 |
| B.4.c CompleteFromTheta | `dynamic/fault_face_flux.cpp` | 141–153 |
| B.4.d CompleteFromVabs | `dynamic/fault_face_flux.cpp` | 155–184 |
| B.4.e BuildImposedState | `dynamic/fault_face_flux.cpp` | 186–223 |
| B.4.f WriteBackState | `dynamic/fault_face_flux.cpp` | 225–238 |
| B.5 EvaluateADER scaling | `dynamic/fault_face_flux.cpp` | 558–606 |

---

## What the trace tells us about the round-2 retracted fixes

- **R-101 (per-face elem1_on_plus)** would change which side gets which flux at B.5. Empirically a no-op on the symmirror mesh because per-QP `sign_flipped` is constant on planar fault faces. No-op confirmed.
- **R-103 (sign-symmetric DG assembly)** would change the FP-rounding sequence at B.5's per-side accumulation. Bit-equivalent under IEEE-754 (`a -= b` and `a += -b` are identical when b is finite). No-op.
- **R-104 (mesh reorient mirror)** would change the per-tet vertex tuples — affects `CalcOrtho` normals consumed at A.3/A.4 face flux loops. Empirically negligible because, on this mesh, both the per-QP `sign_flipped` flag and the per-face `interior_fault_elem1_on_plus_` flag agree. No detectable downstream effect.
- **R-102 (stable shared-face iteration order)** would change the FP-summation order at A.4 only. Cannot affect the np=1 baseline (no shared faces at np=1). May still matter for the (2,4) → (4,2) gap at np=8.

This trace **does not tell us where the np=1 baseline pollution is born**. It tells us where to put the probes (P-1..P-4 above) so that the next experimental round localizes it.
