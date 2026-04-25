# TPV104 σ_n leak — follow-up (afternoon 2026-04-25)

Picks up where `tpv104_sigma_n_leak_root_cause_2026-04-25.md` left off, after that doc concluded the σ_n leak was driven by stress-channel asymmetry in bulk Q. This doc records:

1. Two diagnostics built and run today: the env-gated σ_n freeze, and the C-2 bulk-asymmetry probe set (C-2A BULK-DOF, C-2B NONFAULT-FACE, C-2C BULK-DELTA).
2. Two real bugs found in those diagnostics themselves while interpreting Frontera results: (a) freeze gate doubled σ_n in fluctuation-Q mode; (b) C-2 probes were placed on dead code paths and then on a DOF mapping that picks the apex DOF.
3. Two more nuanced findings revealed by the second Frontera run (job 7677719): the σ_xx asymmetry is mostly physical (anti-symmetric for strike-slip), and the bug-relevant channel is σ_yy (fault-normal).
4. The mapping question's wider context — how SeisSol and our own BP5 quasi-dynamic path determine plus/minus side and per-side DOFs, and what we should adopt for the future BP5↔TPV104 merge.

All commits cited are on branch `feature/elasticity-inertia`.

---

## 1. The σ_n freeze diagnostic — what it is, and the doubling bug

### Goal

The previous doc identified that `Q_global[SXX]` (in fault-local frame, after rotation) on the two sides of the fault becomes asymmetric during rupture, producing a 0.5 MPa σ_n perturbation feeding the friction equation. Plan: pin σ_n at the SCEC TPV104 background of 120 MPa to test how much of the residual V_dip / slip_dip leakage is downstream of the σ_n channel through friction physics, vs independently driven.

### Initial implementation (commit `fada7e0`)

Env-var gate `SEAS_TPV104_FREEZE_SIGMA_N` added to `dynamic/fault_face_flux.cpp::ComputeStageState`, the chokepoint reached by both `Evaluate` (RK4 stages) and `EvaluateADER` (ADER corrector) on the live TPV104 dispatch path. Sets `s.sigma_n_trial = override` immediately after `ComputeTrialTraction`.

Local verification: 25/25 unit tests pass with env unset (byte-identical baseline); `seas_test_tpv102_total_locked_fault` gate-2 worst |sigma_n_corr| pinned exactly to 1.200e+08 Pa with `=120e6` (was 2.083e+08 baseline drift).

Companion sbatch `jobs/tpv104/tpv104_freeze_sigma_n_120MPa.sbatch` mirrored the fwfix run shape (200 m, 8N×400r, 2.0 s, dev queue), with the env var exported before `ibrun`.

### Frontera result that exposed the doubling bug — job 7677661

The freeze run produced `sigma_n_corr = 240 MPa` instead of 120 MPa. ParaView free-surface output at the hypocenter station shows σ_n stuck at 240 MPa for the entire run, while SeisSol reference is at 120 MPa.

### Root cause

`CompleteFromTrial` (`fault_face_flux.cpp:97`) does:

```cpp
s.sigma_n_total = data.sigma_n0 + data.sigma_n_nuc + s.sigma_n_trial;
```

For TPV104 (fluctuation-Q mode, confirmed by `drivers/tpv104_driver.cpp:650, 710, 726`), `data.sigma_n0 = 120 MPa` is in DOFData (not in bulk Q). Setting `s.sigma_n_trial = 120e6` then yields `120 + 0 + 120 = 240 MPa`.

The local `seas_test_tpv102_total_locked_fault` test passed because it calls `ZeroDOFDataPreStressTotal(...)` (its line 176) — TOTAL mode. So the local verification only exercised the case where `data.sigma_n0 = 0`.

### Fix (commit `fa97e12`)

Target `sigma_n_total` directly, not `sigma_n_trial`. One-line change:

```cpp
s.sigma_n_trial = override - data.sigma_n0 - data.sigma_n_nuc
```

- TPV104 fluctuation: override=120e6, sigma_n0=120e6 → s.sigma_n_trial = 0 → sigma_n_total = 120 MPa, Q_imp[SXX] = 0.
- TPV102 total: override=120e6, sigma_n0=0 → s.sigma_n_trial = 120e6 → sigma_n_total = 120 MPa, Q_imp[SXX] = 120e6.

Added regression test `tests/unit/test_tpv104_freeze_sigma_n_gate.cpp` with 6 cases covering both modes (fluctuation at rest, total at rest, fluctuation with bulk perturbation, env unset, =50e6 override). 6/6 pass; the at-rest fluctuation case fails on the pre-fix code (240 MPa instead of 120) — the regression test catches the exact bug we hit.

### Status

Freeze diagnostic is correct and tested. It's a single tool, not a fix; intercepts at LOC-IN (sets sigma_n_trial = 0 for fluctuation), so the friction Riemann sees a frozen σ_n_total = 120 MPa regardless of bulk Q perturbations. Did NOT re-submit the freeze run after the fix was pushed — superseded by the bulk-asymmetry probe (Section 2) which is the more diagnostic-rich tool.

---

## 2. C-2 bulk-asymmetry probe set

### Design (commit `813a2eb`)

Three new probes under `SEAS_DIAG_FAULT_FLUX`, gated single-rank by the existing hypocenter MPI_MINLOC selection. Goal: localize whether the σ_n-channel asymmetry is born in the bulk wave operator (interior non-fault face Godunov flux on the non-mirror mesh) or self-perpetuated by the fault-Riemann back-injection.

| Probe | Location | What it dumps | Volume estimate (2 s, 7000 macro-steps) |
|---|---|---|---|
| **C-2A BULK-DOF** | `wave_operator.inl::Mult` end (later moved to `AdvanceADER` end) | per-Mult: Q± and k± (pre-mass-inverse RHS) at the hypocenter face DOFs of E_plus / E_minus tets | ~7 MB |
| **C-2B NONFAULT-FACE** | `wave_operator.inl::ComputeFaceFluxRHS` interior non-fault branch (later also `ComputeADERFaceFluxRHS`) | per-Mult on each of 6 non-fault faces: Q_self, Q_nbr, F_h, nor, w | ~50 MB |
| **C-2C BULK-DELTA** | `tpv104_driver.cpp` after `Q.Swap(Q_new)` | per-macro-step: Q^{n+1} − Q^n at diag DOFs, diff between sides | ~2 MB |

Driver-side wiring extends the existing hypocenter `MPI_MINLOC` block to also identify the fault face index, the two adjacent tets (via canonical +/- from `FaultBasis::sign_flipped`), the closest local DOF index per tet to the hypocenter QP, and the 6 non-fault interior faces of those tets. All pushed to WaveOperator via three new public setters (`SetDiagBulkElems`, `SetDiagBulkFaceDofs`, `SetDiagNonFaultFaces`).

Local verification: 31/31 unit tests pass with diag flag unset (production byte-identical); diag-on binary has all 11 expected probe format strings present.

Companion sbatch `jobs/tpv104/tpv104_bulk_asymmetry_probe.sbatch` clones `normaltrace_dev` shape, explicitly unsets `SEAS_TPV104_FREEZE_SIGMA_N` (we want the natural asymmetry, not a frozen one), captures each probe stream into a separate `.log` file, runs an awk summary, and produces ParaView output (free-surface at 0.25 s, bulk at 0.5 s).

### First Frontera run — job 7677698 — only C-2C fired

After build, log file showed:

| Tag | Lines |
|---|---|
| C-1n NORMAL | 5500 |
| C-1s ADER-INT | 49500 |
| **C-2A BULK-DOF** | **0** |
| **C-2B NONFAULT-FACE** | **0** |
| C-2C BULK-DELTA | 5500 |
| diag-c2 | 1 |

C-2C fired (it's in the driver loop). C-2A and C-2B fired zero times. Investigation revealed:

### Bug — C-2A and C-2B placed on dead code paths

The TPV104 driver uses `wave.AdvanceADER` (one-shot ADER predictor-corrector at `tpv104_driver.cpp:1084`), which calls `ComputeADERTimeIntegrated + ComputeADERVolumeUpdate + ComputeADERFaceFluxRHS` — never `WaveOperator::Mult` and never `ComputeFaceFluxRHS`. My original placement of C-2A in `Mult` end-of-function and C-2B in `ComputeFaceFluxRHS` interior branch was on dead code for this driver.

### Fix (commit `7285f36`)

Duplicated both probe blocks at the equivalent ADER call sites:
- C-2A: end of `AdvanceADER` body, before `ApplyMassInverse(rhs)`. Reads bulk Q (input) and time-integrated corrector rhs at the diag DOFs.
- C-2B: `ComputeADERFaceFluxRHS` runtime interior non-fault else-branch, right after `flux_.Interior(nor, I_self, I_nbr, F_h)`. Dumps `I_self / I_nbr` (time-integrated states, dt-scaled) and `F_h` (Godunov flux), only on `q==0` (first QP per face per AdvanceADER call) — matches the per-Mult cadence the user chose.

Existing copies in `Mult` and `ComputeFaceFluxRHS` are kept so a future RK4-based driver also gets coverage.

---

## 3. Second Frontera run — job 7677719 — all probes fire

Verified counts:

| Probe | Lines | Per macro-step |
|---|---|---|
| C-1n NORMAL | 6,681 | 1 |
| C-1s ADER-INT (9 tags) | 60,129 | 9 |
| C-2A BULK-DOF (4 lines) | 26,720 | 4 |
| C-2B NONFAULT-FACE (4 lines × 6 faces) | 160,344 | 24 |
| C-2C BULK-DELTA | 6,680 | 1 |

Total stderr 45 MB. All probes fire as designed.

### Channel-by-channel verdict at peak step (n=3721, t=1.06 s)

```
|Q+_SXX − Q−_SXX|_state           = 1.065e+07 Pa  (10.65 MPa, peak)
|Q+_SYY − Q−_SYY|_state at peak   = 4.62e+05 Pa  (peak at t=0.977 s)
|Q+_SZZ − Q−_SZZ|_state           = 1.83e+06 Pa  (peak at t=1.16 s)
d_diff_SXX/step at peak           = -333 Pa/step  (saturated; no longer growing)
```

### Per-step driver asymmetry — C-2A k± SXX

```
k+_SXX = -1.67e+08    k-_SXX = +6.98e+08    diff = -8.64e+08
```

Opposite signs, large magnitudes — the two sides are being driven in opposite directions by the time-integrated rhs.

### Fault-path contribution (b) — C-1s ADER-INT GLOB-IN/GLOB-OUT

```
GLOB-IN+ [SXX]  = -2.875e+06 Pa     GLOB-IN-  = +3.087e+06 Pa
GLOB-OUT+ [SXX] = -2.875e+06 Pa     GLOB-OUT- = +3.087e+06 Pa     (BIT-EXACT MATCH)
b = GLOB-OUT+ − GLOB-OUT- = (input back through, no change in global SXX)
```

The fault back-rotation leaves global SXX EXACTLY UNCHANGED at every step. Reason: BuildImposedState modifies σ_local at the SXX/SXY/SXZ slots (fault-normal-related components in fault-local frame), but for TPV104's vertical y=0 fault with `can_n=(0,−1,0)`, σ_global_xx is reconstructed from the LOCAL tangential channels (SYY/SZZ in fault-local) which BuildImposedState passes through unchanged from input. b ≈ 0 on the global SXX channel — the fault face is conservative for this component.

### Bulk-path contribution (a) — C-2B per-face F_h_SXX at peak step

| E+ side (elem 5611) | F_h_SXX (Pa) | Mirror | E− side (elem 409) | F_h_SXX (Pa) | Δ |
|---|---|---|---|---|---|
| face 10813 | +6.77e+07 | ↔ | face 1339 | +6.44e+07 | +0.33e+07 |
| face 11470 | -4.17e+07 | ↔ | face 1341 | -4.36e+07 | +0.19e+07 |
| **face 12535** | **+3.80e+07** | ↔ | **face 1342** | **−3.87e+07** | **+7.67e+07** ← |

Two pairs are mirror-symmetric (~3% F_h difference, machine-noise level); one pair is anti-symmetric (face 12535 ↔ face 1342, signs opposite, magnitudes equal). The single asymmetric pair contributes ~94% of the bulk asymmetry at peak.

### Mass-balance verdict (initial reading)

```
a / (a + b)  ≈  100%   ALL of the bulk SXX asymmetry is born in the bulk wave op.
b / (a + b)  ≈    0%   The fault Riemann is exactly conservative on global SXX.
```

This **rules out** the previous doc's prime suspect ("rank-2 stress rotation in the fault face flux's imposed state"). The friction-side rotation is clean.

---

## 4. Two follow-up corrections to the verdict above

### 4a. The Q_global[SXX] asymmetry is mostly physical, not a leak

Closer inspection of I_self/I_nbr at the asymmetric pair (face 1342 vs face 12535) revealed:

- Both face normals point in **+x direction** (not mirror normals — same sign!).
- I_self[VX] has opposite signs on the two sides: face 1342 (E−) VX = +1.01e−03, face 12535 (E+) VX = −8.28e−04. **This is the natural strike-slip rupture velocity pattern** (one side slides in +x, the other in −x).
- The Godunov flux for the SXX channel scales as `−(λ+2μ)·v_x`. Opposite VX → opposite F_h_SXX. So the F_h_SXX sign flip is a *consequence* of the strike-slip rupture, not a numerical bug.

For a strike-slip rupture, σ_global_xx is **physically anti-symmetric in y** by mirror reflection (one side under tension, the other under compression as the rupture propagates). The 10.65 MPa Q+_SXX − Q−_SXX asymmetry is the natural anti-symmetric strike-slip σ_xx response — not a numerical artifact.

**The actual bug-relevant channel is global σ_yy** (= fault-normal for TPV104's vertical y=0 fault with `can_n=(0,−1,0)`), which by mirror physics SHOULD be symmetric across the fault. Earlier C-2C showed |Q+_SYY − Q−_SYY|_max = 0.46 MPa at t=0.977 s — this is the actual leak.

Per-face C-2B SYY at the σ_yy peak step (n=3431):

| E− face | F_h_SYY (Pa) | E+ face | F_h_SYY (Pa) | Mirror? |
|---|---|---|---|---|
| 1339 | +6.70e+06 | 10813 | +9.23e+06 | same sign ✓ (38% mag mismatch) |
| 1341 | −8.88e+06 | 11470 | −8.39e+06 | same sign ✓ (6% — near-mirror) |
| **1342** | **−7.88e+06** | **12535** | **+5.04e+06** | **OPPOSITE SIGN ✗** |

Same offending pair as in σ_xx. The asymmetry geometry is consistent across channels — single mesh-topology defect, not per-channel.

I_self / I_nbr at the asymmetric pair (peak σ_yy step):

| | Face 1342 (e1=409) | Face 12535 (e1=5611) | Mirror expectation |
|---|---|---|---|
| I_self[SYY] | −83.8 Pa·s | −58.8 Pa·s | same sign ✓ but 30% mag mismatch |
| **I_nbr[SYY]** | **−13.6 Pa·s** | **+454.6 Pa·s (TENSION!)** | should be same sign ✗ |

The smoking gun: I_nbr at face 12535 (= elem 5960's state) has σ_yy in tension; I_nbr at face 1342 (= elem 550's state) has slight compression. **Outer neighbor tets 5960 and 550 are at qualitatively different stress states** even though they should be mirror partners.

### 4b. The diag-DOF mapping itself may be contaminated

Reading the very early-step C-2C data (BEFORE any rupture activity):

```
step 1 (t=2.85e-4):  Q+_SXX = +1.15e-41   Q-_SXX = +9.11e-42   ratio 1.26
step 2 (t=5.69e-4):  Q+_SXX = -3.81e-13   Q-_SXX = -1.59e-13   ratio 2.40
step 2 SYY:          Q+_SYY = +6.4e-14    Q-_SYY = +4.6e-15    ratio 14× !
step 3 SXX:          Q+_SXX = -1.08e-12   Q-_SXX = -4.59e-13   ratio 2.36
```

For a homogeneous linear-elastodynamic wave equation with `Q = 0` IC and uniform loading, two **mirror-image DOFs** must produce **bit-identical Q[c]**. A consistent 2.4× ratio at step 2 is **not roundoff** — it indicates the two diag DOFs are NOT at mirror positions.

Most likely cause: in `tpv104_driver.cpp` the `closest_dof_idx` lambda finds the closest tet DOF by full Euclidean distance. For a P=1 GaussLobatto L2 tet with 4 nodes (3 face vertices at y=0, 1 apex at |y|≈141 m for 200 m mesh), the apex DOF could in principle win the contest if `hpos` is near a fault-face edge — putting C-2A on E+ at a face vertex (y=0) and on E− at the apex (y≈−141 m). The "asymmetry" then includes a position offset.

### Fix (commit `2fb28a3`)

`closest_face_dof_idx` restricts the search to on-face DOFs (`|dy| < 1 mm` tolerance), breaking ties on (x, z) distance only. Falls back to closest-by-full-distance if no on-face DOF is found, with a non-zero `y_off` reported.

Also adds `[diag-c2-pos]` print emitting the actual physical coordinates of both diag DOFs:

```
[diag-c2-pos] rank=R DOF+ at (X1, 0, Z1) y_off=0   DOF- at (X2, 0, Z2) y_off=0   dx=DX dz=DZ
```

Pass conditions (next run, before drawing any conclusion about bulk asymmetry):
- `y_off=0` on BOTH sides → both DOFs are on the fault face.
- `dx ≈ 0` and `dz ≈ 0` (< 1 m on 200 m mesh) → both at same physical position modulo fault-normal mirror.
- `face_dof+` and `face_dof−` may differ — MFEM tet-local vertex ordering differs per tet (expected).

Production-path byte-identical (closest_face_dof_idx is only called inside `#ifdef SEAS_DIAG_FAULT_FLUX`); 31/31 unit tests pass with diag unset.

---

## 5. ParaView output added (commit `9abb046`)

The bulk-asymmetry sbatch now emits `${RESULT_DIR}/ParaView/` (free surface velocity at 0.25 s — 9 frames over [0, 2 s]) and `${RESULT_DIR}/ParaView_bulk/` (full-volume velocity + sigma_yy + sigma_xy + sigma_xz + mpi_rank at 0.5 s — 5 frames). Lets us visualize the bulk SXX/SYY asymmetry directly on the volumetric mesh, complementing the single-DOF C-2A/C-2C probes.

---

## 6. How SeisSol determines plus/minus (looked at production code)

`/Users/chunhuizhao/projects/SeisSol/src/Geometry/MeshReader.cpp::extractFaultInformation` (lines 102–235), called once at startup:

```cpp
// For each fault face:
//   1. Compute face normal from fixed-ordering of 3 face vertices.
//   2. Read user-supplied refPoint from parameters file.
//   3. isPlus = dot(refPoint - vertex[0], normal) * dot(apex - vertex[0], normal) > 0;
//      (Point method — refPoint side = +side; element is +side iff its apex
//       is on the same side as refPoint.)
//      Or isPlus = dot(refPoint, normal) > 0;  (Normal method — direction)
//   4. f.element = +side element; f.neighborElement = -side element;
//      f.side = local face index in +side tet;
//      f.neighborSide = local face index in -side tet.
```

**Per-face, global, NOT per-QP.** The `Fault` struct stores `element` (always +side), `neighborElement` (always −side), `side` and `neighborSide` (local face indices), all computed at mesh-read time using a USER-SUPPLIED reference point.

For the per-side DOF correspondence at face QPs, SeisSol uses purely topological tables (`FACE2NODES`, `NEIGHBORFACENODE2LOCAL`, `sideOrientations`) computed once at mesh-read by matching the global vertex IDs of the shared face. **No closest-by-Euclidean-distance search anywhere.**

### Robustness comparison

| | Ours (per-QP `sign_flipped`) | SeisSol (per-face refPoint + topological tables) |
|---|---|---|
| Plus/minus assignment | Per-QP, derived from CalcOrtho-vs-ref_normal dot | Per-face, global, user-supplied refPoint |
| Per-side DOF mapping | Closest-by-Euclidean-distance (the diag bug) | Topological table via shared global vertex IDs |
| Apex-DOF risk | Real (fixed by y_tol filter today) | None — table never returns an apex DOF |
| Per-QP sign-flip bimodality | Real (52/48 split observed in earlier audit) | None — assignment is per-face |
| FP drift sensitivity | CalcOrtho ULP near grazing faces | Insensitive (boolean dot product, large margin) |

SeisSol's approach is meaningfully more robust for the mapping problem, and the robustness comes from three design choices we don't currently use:
1. Plus/minus per-face, not per-QP.
2. DOF correspondence by topology, not Euclidean distance.
3. A reference point baked into mesh-read, applied uniformly.

---

## 7. How BP5 (our quasi-dynamic SEAS path) handles this

`fault/fault_nodes.hpp` and `domain/elasticity_operator_traction.inl`:

BP5 stores **per-FACE fault DOFs** (P1 nodes ON the face). Slip is one vector per face DOF — no `+`/`−` copy, no per-side identification. The `+`/`−` asymmetry is handled IMPLICITLY by MFEM's standard DG framework via the IP integrator:

```cpp
// elasticity_operator_traction.inl:557–681 (fault-face traction extraction):
fes_->GetElementVDofs(FTr->Elem1No, vdofs1);     // ALL DOFs of Elem1
fes_->GetElementVDofs(FTr->Elem2No, vdofs2);     // ALL DOFs of Elem2
displacement.GetSubVector(vdofs1, u1_all);
displacement.GetSubVector(vdofs2, u2_all);

trac_integ.ComputeTractionAtQuadPoints(*fe1, *fe2, *FTr,
                                       u1_all, u2_all, delta_u_quad_t,
                                       T_quad_new, ...);
//                  ↑     ↑      ↑    ↑       ↑
//                  Elem1 FE    Elem2 FE    FaceElementTransformations
```

The IP integrator walks face QPs internally; at each QP uses `Loc1.Transform(ip, ip1)` to map face QP → Elem1's reference coords, evaluates `fe1->CalcShape(ip1, shape1)`, dot with `u1_all` → u_at_QP_from_side1. Same for side 2. Computes `[u]` jump, `{σ}` average, traction — all via standard DG framework. **No closest-DOF search anywhere.**

### What this means for our diag probe

| | BP5 (production) | TPV104 production fault Riemann | TPV104 C-2A/C-2C diag probe (pre-fix) |
|---|---|---|---|
| Per-side data extraction | At face QPs via shape1·u1_all, shape2·u2_all | At face QPs via shape1·Q_data, shape2·Q_data ✓ | At a SINGLE TET DOF via closest-by-distance ✗ |
| Robustness to non-mirror mesh | ✓ Standard FE framework | ✓ Same framework | ✗ Closest-distance can pick apex |
| Robustness to per-tet vertex re-ordering | ✓ Loc1/Loc2 handles it | ✓ Same | ✗ |

The TPV104 production fault Riemann path (C-1s GLOB-IN, lines 1014–1022 and 1136–1144 of `wave_operator.inl`) already follows the BP5/SeisSol pattern. **The C-2A/C-2C single-DOF reads are the fragile addition** I introduced for the bulk-side probe; the y_tol fix is OK for catching the apex case, but the more robust replacement is to make C-2A/C-2C also operate at face QPs.

---

## 8. Future BP5 ↔ TPV104 merge architecture

The user asked whether we could call a BP5 face-QP-read function directly without modifying BP5. Found:

- **No single BP5 function** exposes "read bulk Q on each side at a face QP." The pattern is INLINED in `domain/elasticity_operator_traction.inl:557–681` and many other places, plus inside the IP integrator (which is in MFEM library).
- The pattern uses **pure MFEM API** (`FaceElementTransformations::Loc1/Loc2`, `FiniteElement::CalcShape`). So calling MFEM directly from `dynamic/` is equivalent to "calling BP5's pattern" — same primitives.

### Proposed shared helper (to be added later if needed)

```cpp
// common/face_qp_bulk_read.hpp  (NEW — shared, header-only, no BP5 edit)
template <typename FESpaceType>
void ReadBulkQAtFaceQP(FESpaceType &fes, int face_idx,
                       const IntegrationPoint &face_qp_ip,
                       const real_t *Q_data, int ndof_total, int nstate,
                       real_t *Q_self_out, real_t *Q_nbr_out);
```

Implements the standard MFEM/BP5/SeisSol pattern (Loc1/Loc2 + CalcShape + dot) in one place. Both code paths can call it identically. Single point of maintenance for future generalizations (p > 1, precomputed shape-table caching).

Plan when we're ready to add it (only if the y_tol fix isn't enough):
1. Add `common/face_qp_bulk_read.hpp` (new file, no BP5 edit).
2. Replace C-2A and C-2C single-DOF reads with calls to `ReadBulkQAtFaceQP(...)` at the diag fault face's QPs.
3. Refactor C-2B to call the helper too (cosmetic — already at face QPs).
4. Production fault Riemann (`wave_operator.inl:1014–1022`, ADER analogues): leave inline for now, audited.
5. BP5 inline reads: leave alone (no-touch rule). Future BP5 refactor can adopt the helper as opportunistic cleanup with full regression coverage.

---

## 9. Open question, awaiting next Frontera run

`tpv104_bulk_asymmetry_probe.sbatch` was updated with the y_tol filter and ParaView output. Pending re-submission and analysis.

**First diagnostic to read in the next .err file**:

```
[diag-c2-pos] rank=R  DOF+ at (X1, 0, Z1) y_off=0   DOF- at (X2, 0, Z2) y_off=0   dx=DX dz=DZ
```

Branches:

- **`y_off=0` on both sides AND `dx ≈ 0`, `dz ≈ 0`**: diag DOFs are at mirror positions. The σ_yy asymmetry observed in 7677719 (0.46 MPa) is a real bulk leak. Localized source: outer-neighbor tet pair (5960, 550). Next steps: (a) inspect ParaView_bulk for visual confirmation of the mesh non-mirror geometry; (b) consider regenerating the mesh with explicit y↔−y symmetry; (c) optionally implement the BP5-style face-QP read in C-2A/C-2C and re-run for a fully-robust measurement.
- **`y_off ≠ 0` on either side**: y_tol filter still missed; the previous σ_yy result was a probe artifact. Switch C-2A/C-2C to face-QP reads (proposed `ReadBulkQAtFaceQP` helper) before any conclusion.
- **`y_off=0`, `dx, dz ≈ 0`, but σ_yy asymmetry shrunk dramatically**: the previous result was the apex-DOF artifact. Real bulk asymmetry is much smaller; freeze approach can be revisited if needed but probably no longer the priority.

---

## 10. Commit timeline (today, after morning C-1n / C-1s probe work)

| Commit | Title |
|---|---|
| `fada7e0` | TPV104: env-gated σ_n freeze diagnostic + 120 MPa sbatch |
| `fa97e12` | TPV104: fix σ_n freeze gate doubling bug in fluctuation-Q mode |
| `813a2eb` | TPV104: bulk-side C-2 probes (BULK-DOF / NONFAULT-FACE / BULK-DELTA) |
| `7285f36` | TPV104: move C-2A/C-2B probes to ADER paths (Mult/ComputeFaceFluxRHS unused) |
| `9abb046` | TPV104 bulk-asym sbatch: enable ParaView free-surface (0.25s) + bulk (0.5s) |
| `2fb28a3` | TPV104 C-2 probe: restrict diag DOFs to fault-face nodes (apex bug fix) |

All on `feature/elasticity-inertia`. Production-path byte-identical at every commit; 31/31 unit tests pass with diag flag unset throughout.

---

## 11. Files referenced

Source (this branch):
- `dynamic/fault_face_flux.cpp::ComputeStageState` — env-gated σ_n freeze.
- `dynamic/wave_operator.inl::AdvanceADER` end — C-2A BULK-DOF (ADER path).
- `dynamic/wave_operator.inl::ComputeADERFaceFluxRHS` non-fault else-branch — C-2B NONFAULT-FACE (ADER path).
- `dynamic/wave_operator.inl::Mult` end — C-2A duplicate (RK4 path; dead for TPV104).
- `dynamic/wave_operator.inl::ComputeFaceFluxRHS` non-fault else-branch — C-2B duplicate (RK4 path; dead for TPV104).
- `dynamic/wave_operator.hpp` — `SetDiagBulkElems`/`SetDiagBulkFaceDofs`/`SetDiagNonFaultFaces` setters + storage.
- `drivers/tpv104_driver.cpp` — extended hypocenter MPI_MINLOC block, `closest_face_dof_idx` with y_tol, `[diag-c2-pos]` print, C-2C BULK-DELTA per-step.
- `tests/unit/test_tpv104_freeze_sigma_n_gate.cpp` — 6-case regression test for the freeze gate doubling bug.
- `jobs/tpv104/tpv104_freeze_sigma_n_120MPa.sbatch` — freeze diagnostic run.
- `jobs/tpv104/tpv104_bulk_asymmetry_probe.sbatch` — bulk-asymmetry probe run with ParaView output.

Frontera artifacts:
- Job 7677547 (C-1n / C-1s probe) — established 5.97 MPa GLOB-IN SXX asymmetry.
- Job 7677661 (σ_n freeze) — exposed the doubling bug (σ_n stuck at 240 MPa).
- Job 7677698 (first C-2 run) — exposed the dead-code-path bug (C-2A/C-2B silent).
- Job 7677719 (second C-2 run, all probes fire) — established 10.65 MPa SXX, 0.46 MPa SYY, 1.83 MPa SZZ DOF-level asymmetries; localized to face pair (12535, 1342) and outer neighbor tets (5960, 550). σ_xx component recognized as physical anti-symmetric strike-slip response; σ_yy is the actual bug-relevant leak.

External references:
- `/Users/chunhuizhao/projects/SeisSol/src/Geometry/MeshReader.cpp::extractFaultInformation` — refPoint-based +/− assignment.
- `tpv104_sigma_n_leak_root_cause_2026-04-25.md` — predecessor doc (morning), established the σ_n channel as stress-average dominated.
- `tpv104_mesh_asymmetry_finding_2026-04-24.md` — earlier mesh check (0/215143 mirror partners).
