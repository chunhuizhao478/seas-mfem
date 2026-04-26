# TPV104 C-1s STRESS-ROT probe results — Frontera job 7677547

Date: 2026-04-25
Probe site: `dynamic/wave_operator.inl:2304-2362` (ADER interior-fault path)
Probed QP: rank 3, DOF 68, fault coord (-26.6 m, 0, -7550 m) — closest to nominal hypocenter (0, 0, -7500) on tpv104_200m mesh.
Probe count: 7024 ADER-INT BASIS lines, 0 ADER-SHR BASIS lines (production reaches the interior-fault path; no rank-boundary shared-fault dispatch at this QP).

This document records ONLY the verbatim probe outputs and direct numerical comparisons that follow from them. It does NOT extrapolate to claims about whether MFEM matches SeisSol or whether σ_n behavior is "physical" or a "bug" — those questions need separate evidence (a SeisSol reference trace at this hypocenter QP, and a numerical check of etaP/Zp arithmetic) that has not yet been collected.

## 1. BASIS line — rotation matrix at the hypocenter QP

```
[C-1s ADER-INT BASIS] rank=3
   can_n=(-2.9440e-31, -1.0000e+00, -6.1232e-17)
   can_t1=(+1.8027e-47, +6.1232e-17, -1.0000e+00)
   can_t2=(+1.0000e+00, -2.9440e-31, +0.0000e+00)
   sign_flipped=1   elem1_on_plus=0   dt=2.8477e-04
```

Numerical observations:
- can_n is `(0, -1, 0)` to within 1e-17 (machine precision).
- can_t1 is `(0, 0, -1)` to within 1e-17.
- can_t2 is `(1, 0, 0)` exactly.
- All off-diagonal entries are at 1e-17 to 1e-31, which are floating-point round-off, not physical.
- This QP is on the sign_flipped=1 half of the bimodal distribution found in earlier audits (52/48 split).
- elem1_on_plus=0 means Elem1 is on the canonical-minus side at this rank.
- dt = 2.8477e-04 s (the macro-step size used by EvaluateADER).

## 2. Step 1 (early — first macro step, before nucleation)

```
GLOB-IN+   SXX=+0  SYY=+0  SZZ=+0  SXY=+0  SYZ=+0  SXZ=+0  VX=+0  VY=+0  VZ=+0
GLOB-IN-   SXX=+0  SYY=+0  SZZ=+0  SXY=+0  SYZ=+0  SXZ=+0  VX=+0  VY=+0  VZ=+0
LOC-IN+    SXX=+0  SYY=+0  SZZ=+0  SXY=+0  SYZ=+0  SXZ=+0  VX=+0  VY=+0  VZ=+0
LOC-IN-    SXX=+0  SYY=+0  SZZ=+0  SXY=+0  SYZ=+0  SXZ=+0  VX=+0  VY=+0  VZ=+0
LOC-IMP+   SXX=+0  SYY=+0  SZZ=+0  SXY=+0  SYZ=+0  SXZ=-4.6244e-10  VX=+0  VY=+0  VZ=-5.0000e-17
LOC-IMP-   SXX=+0  SYY=+0  SZZ=+0  SXY=+0  SYZ=+0  SXZ=-4.6244e-10  VX=+0  VY=+0  VZ=+5.0000e-17
GLOB-OUT+  SXX=+2.7229e-40  SYY=-2.7229e-40  SZZ=+0
           SXY=+4.6244e-10  SYZ=-8.3365e-57  SXZ=+2.8317e-26
           VX=-5.0000e-17  VY=+1.4720e-47  VZ=+0
GLOB-OUT-  SXX=+2.7229e-40  SYY=-2.7229e-40  SZZ=+0
           SXY=+4.6244e-10  SYZ=-8.3365e-57  SXZ=+2.8317e-26
           VX=+5.0000e-17  VY=-1.4720e-47  VZ=+0
```

Direct numerical observations:
- All bulk inputs (GLOB-IN, LOC-IN) are exactly zero. Q starts at rest.
- LOC-IMP (output of `EvaluateADER` in fault-local frame) has a non-zero SXZ = -4.62e-10 Pa and VZ = ±5.0e-17 m/s on both sides. The remaining 7 components are exactly zero.
- GLOB-OUT (LOC-IMP rotated back to global) has a non-zero SXY = +4.62e-10 Pa. The other off-diagonal entries are at 1e-26 to 1e-57 (effectively zero).
- The LOC-IMP SXZ → GLOB-OUT SXY mapping is consistent with the rotation: with can_t2 = (+1, 0, 0), `σ_xy_global` is read from `σ_xz_local` (SXZ in our index naming where local x=normal, z=t2; σ_n_t2_local has y_global×x_global content after rotation back).
- The +4.62e-10 Pa is six orders of magnitude smaller than typical TPV104 stresses (~1e7 Pa) but is non-zero.

What this rules in / out at step 1:
- The wave operator does not introduce non-zero σ_yy_global at step 1 from a zero input. The fault path's `EvaluateADER` produces a tiny SXZ_local = SXY_global of order 1e-10 Pa, of unidentified origin (likely friction-corrected traction at zero V; not yet investigated).
- σ_yy_global is at 1e-40 (round-off), not at 1e-10. So the imposed-state at step 1 does NOT inject σ_yy.

## 3. Step 3766 (peak |σ_n_trial| per C-1n NORMAL probe)

C-1n NORMAL at peak (job 7677329, same QP, prior run): peak |σ_n_trial| = 5.279e+05 Pa, stress_term = -5.994e+05 Pa, v_jump_term = +7.156e+04 Pa.

C-1s ADER-INT at step 3766 (job 7677547):

```
GLOB-IN+   SXX=-2.9088e+06  SYY=-5.9973e+05  SZZ=-3.1242e+05
           SXY=+5.4903e+07  SYZ=+6.4927e+04  SXZ=-1.9532e+06
           VX=-3.8126e+00   VY=+2.3801e-02   VZ=+1.2869e-02

GLOB-IN-   SXX=+3.0322e+06  SYY=-5.9912e+05  SZZ=-7.4568e+05
           SXY=+5.5070e+07  SYZ=+1.0027e+04  SXZ=+2.2795e+06
           VX=+3.8217e+00   VY=+1.4867e-02   VZ=-4.1228e-03

LOC-IN+    SXX=-5.9973e+05  SYY=-3.1242e+05  SZZ=-2.9088e+06
           SXY=+6.4927e+04  SYZ=+1.9532e+06  SXZ=-5.4903e+07
           VX=-2.3801e-02   VY=-1.2869e-02   VZ=-3.8126e+00

LOC-IN-    SXX=-5.9912e+05  SYY=-7.4568e+05  SZZ=+3.0322e+06
           SXY=+1.0027e+04  SYZ=-2.2795e+06  SXZ=-5.5070e+07
           VX=-1.4867e-02   VY=+4.1228e-03   VZ=+3.8217e+00

GLOB-OUT+  SXX=-2.9088e+06  SYY=-5.2786e+05  SZZ=-3.1242e+05
           SXY=+5.4382e+07  SYZ=+5.4387e+04  SXZ=-1.9532e+06
           VX=-3.7562e+00   VY=+1.9315e-02   VZ=+1.4009e-02

GLOB-OUT-  SXX=+3.0322e+06  SYY=-5.2786e+05  SZZ=-7.4568e+05
           SXY=+5.4382e+07  SYZ=+5.4387e+04  SXZ=+2.2795e+06
           VX=+3.7473e+00   VY=+1.9315e-02   VZ=+6.7345e-04
```

### 3a. Rotation verification (LOC = Tinv · GLOB)

For can_n=(0,-1,0), can_t1=(0,0,-1), can_t2=(+1,0,0), the local-frame stress is a permutation of global-frame stress (with sign rules from the squared/cross products in the rank-2 rotation):

| Local | = | Global | Plus side check | Minus side check |
|---|---|---|---|---|
| LOC SXX | = | GLOB σ_yy | -5.9973e5 = -5.9973e5 ✓ | -5.9912e5 = -5.9912e5 ✓ |
| LOC SYY | = | GLOB σ_zz | -3.1242e5 = -3.1242e5 ✓ | -7.4568e5 = -7.4568e5 ✓ |
| LOC SZZ | = | GLOB σ_xx | -2.9088e6 = -2.9088e6 ✓ | +3.0322e6 = +3.0322e6 ✓ |
| LOC SXY | = | GLOB σ_yz | +6.4927e4 = +6.4927e4 ✓ | +1.0027e4 = +1.0027e4 ✓ |
| LOC SXZ | = | -GLOB σ_xy | -5.4903e7 = -(+5.4903e7) ✓ | -5.5070e7 = -(+5.5070e7) ✓ |
| LOC SYZ | = | -GLOB σ_xz | +1.9532e6 = -(-1.9532e6) ✓ | -2.2795e6 = -(+2.2795e6) ✓ |
| LOC VX | = | -GLOB v_y | -2.3801e-2 = -(+2.3801e-2) ✓ | -1.4867e-2 = -(+1.4867e-2) ✓ |
| LOC VY | = | -GLOB v_z | -1.2869e-2 = -(+1.2869e-2) ✓ | +4.1228e-3 = -(-4.1228e-3) ✓ |
| LOC VZ | = | GLOB v_x | -3.8126 = -3.8126 ✓ | +3.8217 = +3.8217 ✓ |

All 18 entries (9 components × 2 sides) match the expected permutation/sign at peak step 3766. The rotation Tinv produces fault-local components that are exact permutations of global components, to all printed digits (4-5 significant figures).

### 3b. Plus/minus symmetry per component (peak)

| Component (global frame) | Plus | Minus | (Plus − Minus) | rel.diff |
|---|---|---|---|---|
| σ_xx | -2.9088e6 | +3.0322e6 | -5.94e6 | sign-anti-symmetric |
| **σ_yy** | **-5.9973e5** | **-5.9912e5** | **-6.1e2** | **0.1%** |
| σ_zz | -3.1242e5 | -7.4568e5 | +4.33e5 | not symmetric |
| σ_xy | +5.4903e7 | +5.5070e7 | -1.67e5 | 0.3% |
| σ_yz | +6.4927e4 | +1.0027e4 | +5.49e4 | not symmetric |
| σ_xz | -1.9532e6 | +2.2795e6 | -4.23e6 | sign-anti-symmetric |
| v_x | -3.8126 | +3.8217 | -7.63 | sign-anti-symmetric |
| v_y | +2.3801e-2 | +1.4867e-2 | +8.93e-3 | 60% |
| v_z | +1.2869e-2 | -4.1228e-3 | +1.69e-2 | sign-flipped |

Direct observation: σ_yy_global is symmetric across the fault to 0.1% at peak (Plus = -5.9973e5, Minus = -5.9912e5). Other components show sign-flips (σ_xx, σ_xz, v_x) or larger asymmetries (σ_zz, σ_yz, v_y, v_z) at the same time step.

### 3c. GLOB-IN → GLOB-OUT change (what the fault path imposes back)

| Component | GLOB-IN+ | GLOB-OUT+ | Δ (out − in) |
|---|---|---|---|
| σ_xx | -2.9088e6 | -2.9088e6 | 0 |
| **σ_yy** | **-5.9973e5** | **-5.2786e5** | **+7.187e4** |
| σ_zz | -3.1242e5 | -3.1242e5 | 0 |
| σ_xy | +5.4903e7 | +5.4382e7 | -5.21e5 |
| σ_yz | +6.4927e4 | +5.4387e4 | -1.054e4 |
| σ_xz | -1.9532e6 | -1.9532e6 | 0 |
| v_x | -3.8126 | -3.7562 | +5.64e-2 |
| v_y | +2.3801e-2 | +1.9315e-2 | -4.49e-3 |
| v_z | +1.2869e-2 | +1.4009e-2 | +1.14e-3 |

| Component | GLOB-IN- | GLOB-OUT- | Δ (out − in) |
|---|---|---|---|
| σ_xx | +3.0322e6 | +3.0322e6 | 0 |
| **σ_yy** | **-5.9912e5** | **-5.2786e5** | **+7.126e4** |
| σ_zz | -7.4568e5 | -7.4568e5 | 0 |
| σ_xy | +5.5070e7 | +5.4382e7 | -6.88e5 |
| σ_yz | +1.0027e4 | +5.4387e4 | +4.44e4 |
| σ_xz | +2.2795e6 | +2.2795e6 | 0 |
| v_x | +3.8217 | +3.7473 | -7.44e-2 |
| v_y | +1.4867e-2 | +1.9315e-2 | +4.45e-3 |
| v_z | -4.1228e-3 | +6.7345e-4 | +4.80e-3 |

Direct observations on the imposed change:
- σ_xx (= LOC σ_zz), σ_zz (= LOC σ_yy), σ_xz (= LOC σ_yz): Δ = 0 on both sides — these three components are passed through unchanged from input to imposed state. Cross-check with `BuildImposedState` source (`fault_face_flux.cpp:163-189`): only the three traction components (LOC SXX, LOC SXY, LOC SXZ = σ_n, τ_1, τ_2) and the three velocities are overwritten; the three "passed-through" stress components (LOC SYY, LOC SZZ, LOC SYZ) are kept verbatim. The data confirms: GLOB σ_xx (=LOC SZZ), GLOB σ_zz (=LOC SYY), GLOB σ_xz (=LOC SYZ) are unchanged.
- σ_yy (= LOC SXX = σ_n_local): Δ = +7.187e4 Pa on plus side, +7.126e4 Pa on minus side. Both sides receive the same change to within 8.5%. The imposed σ_yy is identical on the two sides (-5.2786e5 = -5.2786e5). Before the imposition, plus and minus differ by 6.1e2 Pa; after, they are exactly equal.
- σ_xy (= LOC SXZ, the τ_2 traction): Δ = -5.21e5 Pa on plus, -6.88e5 Pa on minus. Imposed σ_xy is identical on both sides (+5.4382e7 = +5.4382e7).
- σ_yz (= LOC SXY, the τ_1 traction): Δ = -1.054e4 plus, +4.44e4 minus. Imposed σ_yz is identical on both sides (+5.4387e4 = +5.4387e4).

This is consistent with `BuildImposedState` overwriting the three traction components on both sides with the friction-corrected values (`s.sigma_n_corr`, `s.tau1_corr`, `s.tau2_corr`).

## 4. Last sample (post-rupture)

```
GLOB-IN+   SXX=-9.2233e+05  SYY=-1.6010e+05  SZZ=-4.8167e+05
           SXY=+5.9866e+07  SYZ=+3.1145e+04  SXZ=-8.3998e+05
GLOB-IN-   SXX=+8.8855e+05  SYY=-1.7822e+05  SZZ=-6.5022e+04
           SXY=+5.9876e+07  SYZ=+3.1247e+04  SXZ=+8.5472e+05

LOC-IN+    SXX=-1.6010e+05  ...  (= GLOB σ_yy ✓)
LOC-IN-    SXX=-1.7822e+05  ...  (= GLOB σ_yy ✓)

GLOB-OUT+  SXX=-9.2233e+05  SYY=-1.6711e+05  SZZ=-4.8167e+05
GLOB-OUT-  SXX=+8.8855e+05  SYY=-1.6711e+05  SZZ=-6.5022e+04
```

Late-time σ_yy_global comparison:
| | Plus | Minus | (Plus − Minus) | rel.diff |
|---|---|---|---|---|
| GLOB-IN σ_yy | -1.6010e5 | -1.7822e5 | +1.812e3 | 10.6% |
| GLOB-OUT σ_yy | -1.6711e5 | -1.6711e5 | 0 | 0% (equal by imposition) |

The plus/minus σ_yy asymmetry at the GLOB-IN side has grown from 0.1% (at peak) to 10.6% (post-rupture). The imposed σ_yy at GLOB-OUT remains identical on both sides (as expected from `BuildImposedState`).

## 5. What can be concluded directly from this probe

Verified by data (numerical-comparison observations, not extrapolations):

1. The rotation matrix at the hypocenter QP is a clean axis permutation `(x_local→n=−y_global, y_local→t1=−z_global, z_local→t2=+x_global)` to within 1e-17. There is no measurable off-diagonal contamination.
2. LOC-IN components are exact permutations/sign-changes of GLOB-IN components, verified at peak step 3766 across all 18 component values to 4-5 sig fig. The Tinv rotation does not introduce or remove any signal.
3. At step 1 (Q ≡ 0 input), LOC-IMP has SXZ = -4.62e-10 Pa and VZ = ±5.0e-17 m/s; σ_yy_imp_global is at 1e-40 (round-off). The fault path produces tiny non-zero values at the start, magnitudes ~6 orders below typical run values, of unidentified origin.
4. At peak step 3766, σ_yy_global is symmetric between plus and minus sides (0.1% asymmetry, -5.9973e5 vs -5.9912e5 Pa). σ_xx, σ_xz, v_x are sign-anti-symmetric. σ_zz, σ_yz, v_y, v_z have larger asymmetries.
5. `BuildImposedState` overwrites three components (σ_yy, σ_xy, σ_yz in global frame, corresponding to σ_n, τ_2, τ_1 in local frame at this QP) with values that are identical on both sides. The other six components (σ_xx, σ_zz, σ_xz, v_x, v_y, v_z) are modified per the per-side velocity-update formula, NOT enforced equal across sides.
6. The σ_yy_global asymmetry at GLOB-IN is 0.1% at peak and grows to 10.6% by the end of the run. The σ_yy at GLOB-OUT is enforced equal on both sides at every step.

## 6. What cannot be concluded from this probe alone

The probe captures the per-step pipeline at one QP. It does NOT, by itself, answer:

- Whether σ_yy ≈ -0.6 MPa at the hypocenter is what SeisSol's TPV104 produces at the same time. (Need: SeisSol reference σ_n trace at hypocenter station; not collected.)
- Whether the C-1n NORMAL stress_term -5.994e5 Pa numerically equals etaP·(Q+[SXX]/Zp + Q-[SXX]/Zp_neig) using the actual TPV104 etaP and Zp values. (Need: arithmetic check using `TPV104Params::rho`, `TPV104Params::cp`, `TPV104Params::cs` to derive etaP, Zp; not done.)
- Whether the bulk σ_yy at this QP comes from radiation (waves arriving from elsewhere on the fault) or from imposed-state injection (this QP's own previous step pushing σ_yy back into the bulk). (Need: probe at non-fault bulk DOFs adjacent to the hypocenter, or comparison of GLOB-OUT_t and GLOB-IN_{t+1}; not done.)
- Whether per-QP `sign_flipped=1` vs `sign_flipped=0` produces a measurable difference in σ_yy at neighboring fault QPs. (Need: tagging a sign_flipped=0 QP as a probe DOF for comparison; only one QP currently tagged.)

## 7. Files

- Probe source: `dynamic/wave_operator.inl:2304-2362` (commit 72d4110).
- Raw probe data on Frontera: `tpv104_normaltrace_7677547.err` in submit dir.
- Hypocenter station data on Frontera: `tpv104/results_normaltrace_job7677547/stations/tpv104_normaltrace_station_x2_0_x3_7.5.dat`.
- Earlier C-1n NORMAL finding: `tpv104_sigma_n_leak_root_cause_2026-04-25.md`.
