# Mixed-Flux + Heterogeneous (Bi-material) Riemann Solver — Code Exploration

*Worktree:* `.claude/worktrees/mixed-flux-hetero-riemann` (branch `feat/mixed-flux-hetero-riemann`, off `system/spatial_dyn_driver` @ 04674fc)
*Reference:* `/Users/chunhuizhao/projects/drdg3d` — Zhang, Liu & Chen (2023), *A Mixed-Flux-Based Nodal DG Method for 3D Dynamic Rupture* (JGR Solid Earth, 128, e2022JB025817).
*Scope:* document how "mixed flux" and the "heterogeneous (bi-material) Riemann solver" are realized in drdg3d and in this SEAS-MFEM tree, and pin down why they are currently mutually exclusive.

---

## Overview

Two numerical ingredients are at play, and drdg3d uses them **together**:

1. **Heterogeneous (bi-material) Riemann solver** — the interior-face numerical flux uses
   *per-side* impedances `(Z⁻, Z⁺)` instead of a single impedance. This is the exact
   linearised elastic Riemann solver across a material contrast; it reduces to the standard
   upwind/Godunov flux when the two sides match.
2. **Mixed flux** — a per-face switch between the **upwind (Godunov)** flux on most interior
   faces and a **central (non-dissipative) averaged flux** on faces *adjacent to the fault*.
   The central flux removes upwind dissipation in the near-fault region (Zhang 2023 §3.2,
   Fig. 5), which is what gives the method its accuracy on rupture fronts.

In this SEAS-MFEM tree **both ingredients already exist** but are wired on **separate,
mutually-exclusive code paths**: mixed flux lives on the scalar `WaveOperator` (single
impedance), and the bi-material Riemann solver lives on `BimaterialWaveOperator` (per-element
materials, no mixed flux). The combination — which is exactly what drdg3d does — is currently
forbidden by a structural guard.

---

## Architecture

### Reference: drdg3d (Fortran, nodal DG)

| File | Role |
|------|------|
| `src/mod_eqns.F90` | Physical flux Jacobians `Flux1/2/3`; the 1-D characteristic solver `riemannSolver_continuous` (lines 187–209). |
| `src/mod_wave.F90` | `RHS` (volume + LIFT) and `get_flux` (lines 347–1040) — the numerical-flux kernel for every face, including the mixed-flux switch and the fault. |
| `src/mod_fault.F90` | Fault arrays init / IO. |
| mesh fields `fluxtype(Nfaces,nelem)`, param `flux_method` | Per-face mixed-flux selector + global method id (read in `mod_read.F90:105`, default **2**). |

### This tree: `miniapps/seas/dynamic/`

| File | Role |
|------|------|
| `godunov_flux.{hpp,cpp}` | **Homogeneous** Godunov flux. `Interior()` = upwind `A⁺Q_self + A⁻Q_nbr`; `Central()` = `½·A·(Q_self+Q_nbr)` (the mixed-flux averaged flux). Single impedance `Zp_,Zs_` set at construction. |
| `godunov_flux_bimaterial.{hpp,cpp}` | **Heterogeneous** Riemann solver (SeisSol/Pelties 2012 formulation). Builds per-face 9×9 projectors from **per-side** `(λ,μ,ρ)`; runtime apply is `F_h_self = fluxLocal·Q_self + fluxNeighbor·Q_nbr`. Reduces byte-exactly to `GodunovFlux::Interior` in the homogeneous limit (test R.1.T-1). |
| `godunov_flux_pool.{hpp,cpp}` | De-duplicated table of per-element `GodunovFlux` objects (one per distinct material triple). |
| `heterogeneous_material.{hpp,cpp}` | `MaterialField` — Constant / GridFunction / Coefficient modes; supplies `(λ,μ,ρ)` per element or per QP. |
| `wave_operator.{hpp,inl}` | Scalar base operator. Owns the mixed-flux machinery: `MixedFluxMode` enum, `SetMixedFluxMode`, `BuildCentralFluxFaceSet_`, and the per-face dispatch `InteriorFaceFlux_`/`SharedInteriorFaceFlux_`. |
| `bimaterial_wave_operator.{hpp,inl}` | Subclass of `WaveOperator`. Overrides the material-dependent hooks to apply the per-face bi-material matrices. **`SetMixedFluxMode` is overridden to ABORT on any non-None mode.** |
| `fault_face_flux.{hpp,cpp}` | Fault interface flux + friction. **Already bi-material**: per-DOF `Zp_±,Zs_±`, harmonic means `eta_p,eta_s`, trial traction from both sides. |
| `precomputed_face_fluxes.{hpp,cpp}` | Optional precomputed homogeneous per-face flux matrices (`use_precomputed_face_fluxes_`); mutually exclusive with mixed flux (R-1203). |

Driver: `drivers/spatial_dyn_driver.cpp` constructs **either** `WaveOperator` (config
`[numerics].interior_flux="scalar"`) **or** `BimaterialWaveOperator`
(`interior_flux="matrix"`) — lines 1099–1138 — and then calls
`wave.SetMixedFluxMode(ParseMixedFlux(cfg.numerics.mixed_flux))` at line **1245**.

---

## Call Graph (interior-face numerical flux)

```
spatial_dyn_driver::main
  ├─ wave_ptr = WaveOperator | BimaterialWaveOperator        (1099–1138, config switch)
  ├─ wave.SetMixedFluxMode(ParseMixedFlux(cfg…))             (1245)  ← ABORTS on matrix+non-None
  └─ time loop → wave.Mult → …Apply face flux
       WaveOperator::InteriorFaceFlux_ (inl:1193)
         if (mf_on_ && central_flux_face_set_.count(face))
              flux_.Central(nor,Q_self,Q_nbr,F_h)            (½ A (Q_self+Q_nbr))
         else flux_.Interior(nor,Q_self,Q_nbr,F_h)           (upwind, single Z)
       BimaterialWaveOperator::InteriorFaceFlux_ (inl:471)
         BimaterialFlux::ApplyPerFaceFlux(mat_e1_local, mat_e1_nbr, Q_self, Q_nbr, F_h_e1)
         BimaterialFlux::ApplyPerFaceFlux(mat_e2_local, mat_e2_nbr, Q_self, Q_nbr, F_h_e2)
         (no central branch — mixed flux is structurally disabled here)
```

---

## Detailed Walkthrough

### drdg3d `get_flux` (mod_wave.F90:347–1040) — the canonical combination

Per face quadrature node the kernel pulls the **neighbour** material
(`rho_out,cp_out,cs_out` from `mesh%rho(neigh)` / MPI neighbour, lines 438–467), forms
per-side impedances `Zp=ρc_p, Zs=ρc_s` and `Zp_out,Zs_out` (lines 416–417, 484–485), then:

- **Mixed-flux central branch** (lines 587–618): if `flux_method==2 .and. fluxtype==1`
  (or `flux_method==1` for any interior face) it sets
  `fstar = ½·(ΔT, n·ΔV)` — the **averaged physical-flux jump**, *impedance-independent* —
  and jumps past the Riemann solve. This is the central flux used on fault-adjacent faces.
- **Bi-material Godunov branch** (lines 675–688): otherwise it forms the upwind state with
  **both** impedances:

  ```
  dF = FR - FL
  alpha(1) = (dF(1) + Zp_out·dF(4)) / (Zp + Zp_out)      ! P
  alpha(2) = (dF(2) + Zs_out·dF(9)) / (Zs + Zs_out)      ! S (m)
  alpha(3) = (dF(3) + Zs_out·dF(8)) / (Zs + Zs_out)      ! S (l)
  fstar(1)=alpha1·Zp; fstar(4)=alpha1; fstar(8)=alpha3; fstar(9)=alpha2; …
  ```

  These are the same `(Z⁻,Z⁺)` weights as `riemannSolver_continuous`
  (`eta=z_p·z_m/(z_p+z_m)`, mod_eqns.F90:200). Equal impedances → standard upwind.
- **Fault branch** (lines 697+): `eta = zs_in·zs_out/(zs_in+zs_out)` (harmonic-mean shear
  impedance) couples slip rate to the trial "locked" traction
  `Tau_lock = fL + fstar + Tau0`, then the friction return map.

**Takeaway:** in drdg3d the central (mixed-flux) branch and the bi-material upwind branch are
two arms of *one* per-face `if`. They are not mutually exclusive — `fluxtype` selects central
on fault-adjacent faces and bi-material Godunov everywhere else, and the fault uses the
harmonic-mean impedance. This is the target behaviour for SEAS-MFEM.

### SEAS-MFEM `GodunovFlux::Central` (godunov_flux.cpp:382–423)

Computes `F_h = ½·(Ax_plus_ + Ax_minus_)·(Q_self_rot + Q_nbr_rot)` in the face frame, i.e.
`½·A_self·(Q_self+Q_nbr)`. Because the scalar operator is **homogeneous**, `A_self == A_nbr`,
so this equals the averaged physical flux `½(A_self Q_self + A_nbr Q_nbr)` — drdg3d's central.
**This equivalence breaks under heterogeneity** (see Gotchas).

### SEAS-MFEM `BimaterialFlux` (godunov_flux_bimaterial.hpp)

`BuildPerFaceFluxMatricesGlobal(nor, flux_self, flux_nbr, fluxLocal, fluxNeighbor)` precomputes
two global-frame 9×9 matrices per face/side from each side's `GodunovFlux`. There is **no
central counterpart** — only the upwind (Godunov) projector is built. The rotation machinery
(`T·A_facelocal·T⁻¹`) is exactly what a central builder would reuse.

### The mutual-exclusion guard

`BimaterialWaveOperator::SetMixedFluxMode` (bimaterial_wave_operator.inl:567–576):

```cpp
MFEM_VERIFY(m == MixedFluxMode::None,
  "BimaterialWaveOperator::SetMixedFluxMode: mixed flux is incompatible with the "
  "heterogeneous (bimaterial) interior_flux=\"matrix\" path …");
```

Plus the scalar-path cross-check `mixed flux ⟂ use_precomputed_face_fluxes_`
(wave_operator.inl:1613–1621, R-1203), and the design note that R-003 made the
matrix×mixed-flux exclusion *structural* (wave_operator.inl:1623–1626). So today
`interior_flux="matrix"` + `mixed_flux≠none` **aborts at driver line 1245**.

### `BuildCentralFluxFaceSet_` (wave_operator.inl:1650+)

Already mode-aware and **operator-agnostic** (it lives on the base `WaveOperator`): `Adjacent`
inserts interior/shared faces that touch a fault element but are not the fault face itself;
`AllContinuous` inserts every interior face. It uses rank-symmetric fault-face sets and
excludes non-fault BC faces (R-1404). **This builder is reusable as-is by the bi-material
path** — nothing in it depends on the homogeneity of the flux.

---

## Conventions

- **State vector** `Q` is 9-component: `(ρVx,ρVy,ρVz, Exx,Eyy,Ezz,Eyz,Exz,Exy)` (strain-velocity),
  matching drdg3d's `U`.
- **Sign / frame:** fault-local frame is the BP5/Tandem convention `tangent1=dip, tangent2=strike`
  (see `miniapps/seas/CLAUDE.md`). Bi-material `matR` follows SeisSol's sign convention verbatim
  (header of `godunov_flux_bimaterial.hpp`) — *not* MFEM's `GodunovFlux::R`.
- **Parallelism:** MPI; mixed-flux mode must be set with the **same mode in the same order on
  every rank** or the post-walk `Allgatherv` deadlocks (wave_operator.inl:1605). Bi-material
  precompute is per-rank with a 4 GiB soft cap warning.
- **Config:** TOML `[numerics] interior_flux = "scalar"|"matrix"`,
  `mixed_flux = "none"|"adjacent"|"all_continuous"` (CLI overrides at driver:563/679).

---

## Gotchas

- **Central flux is NOT impedance-weighted, but its bi-material form is subtle.** drdg3d's
  central branch averages the *physical* traction/velocity jumps: `F* = ½(A_L Q_L + A_R Q_R)`
  (single-valued, conservative). MFEM's `GodunovFlux::Central` computes `½·A_self·(Q_self+Q_nbr)`
  — which uses `A_self` for *both* terms. These agree **only when `A_L==A_R`**. A naive reuse of
  `Central()` per side under heterogeneity gives a *non-single-valued* flux (different on each
  side) and breaks conservation. The correct bi-material central is `F* = ½(A_self Q_self +
  A_nbr Q_nbr)`, deposited identically to both sides — per-face matrices `½A_self` (mult Q_self)
  and `½A_nbr` (mult Q_nbr). This is the single most important design decision for the merge.
- **The fault interface is already bi-material** (`fault_face_flux`: per-side `Zp_±/Zs_±`,
  harmonic `eta`). The merge only needs to fix the *interior* fault-adjacent faces, not the fault
  flux itself.
- **CFL.** `BimaterialWaveOperator::ComputeMaxDt` (inl:515) currently assumes mixed flux is off
  (factor = 1). The scalar path applies a mixed-flux CFL factor; enabling mixed flux on the
  bi-material path must apply the same factor in the heterogeneous CFL walk.
- **`use_precomputed_face_fluxes_` is a third, separate path** (homogeneous, precomputed) also
  mutually exclusive with mixed flux (R-1203). The merge concerns the `matrix` path, not this one.
- **Homogeneous byte-exactness is a hard contract.** `test_bimaterial_wave_operator_parity`
  (22/22) asserts the bi-material operator matches the scalar operator to ~1e-9 on homogeneous
  material; any central-flux addition must preserve `mixed_flux=none` parity exactly and reduce
  to `GodunovFlux::Central` in the homogeneous limit.

---

## Open Questions

- **Heterogeneity location for SAFS:** is the material contrast primarily *bulk* (CVM-H velocity
  model, smooth) with the fault embedded in near-uniform rock, or a genuine *across-fault*
  contrast (different blocks)? This determines whether fault-adjacent central faces actually
  straddle an impedance jump in production runs (affects test priority, not correctness).
- **Does Adjacent mode need the bi-material central on shared (cross-rank) fault-adjacent faces?**
  `SharedInteriorFaceFlux_` has only side-0 matrices populated; a shared fault-adjacent face needs
  the neighbour's `A_nbr` from the face-nbr material exchange — verify that path supplies it.
- **CFL factor value** for mixed flux on the bi-material path — confirm the scalar operator's
  factor and whether it is mode-dependent (Adjacent vs AllContinuous).
