# drdg3d vs SEAS-MFEM: bimaterial central / mixed flux comparison

Date: 2026-06-06
Question: drdg3d has a similar bimaterial path with central flux — what is missing
in our code that makes TPV31 (matrix path) leak sigma_n at the material interfaces
while their central flux "works"?

drdg3d source: `/Users/chunhuizhao/projects/drdg3d/src/mod_wave.F90` (subroutine
`get_flux`, lines ~347-1000).  All math ASCII (no LaTeX) per project convention.

---

## TL;DR (the answer)

There is **no magic bimaterial central flux in drdg3d**.  Their central flux is the
*same non-dissipative average ours is*.  The bimaterial impedance contrast in drdg3d
is handled entirely by the **upwind (impedance-weighted Riemann)** flux — which we
**also already have** (our Pelties bimaterial Godunov).

What is missing in our code is the **restriction**: drdg3d's central (non-dissipative)
flux is only ever applied on faces *within homogeneous material*; the material
contrast is handled by the impedance-weighted upwind.  Our
`BimaterialWaveOperator::BuildPerFaceCentralFluxMatrices_` applies the central flux to
**every** fault-adjacent corridor face with **no material-jump guard**
(`dynamic/bimaterial_wave_operator.inl:515-690`), so corridor faces that straddle the
TPV31 material steps at 2400/5000/10000 m get a non-dissipative, non-impedance-weighted
flux across an impedance contrast → the sigma_n leak localized to exactly those depths.

Fix = keep central inside homogeneous material; use the impedance-weighted bimaterial
upwind (or a theta-blend toward it) on corridor faces that cross a material contrast.

---

## Architecture (drdg3d `get_flux`)

drdg3d computes ALL interior/fault face fluxes in one routine, `get_flux`
(`mod_wave.F90:347`).  Two global switches + one per-face tag select the scheme:

- `flux_method` (config `get_params("flux_method",...,2)`, **default 2**):
  - `== 1`  -> PURE CENTRAL at every interior (BC_IN) face (`mod_wave.F90:604-618`).
  - `== 2`  -> MIXED: central only where `fluxtype==1`, impedance-weighted upwind
              everywhere else (`mod_wave.F90:587-602` + fall-through `675-688`).
- `mesh%fluxtype(is,ie)` : per-(face,element) tag, value `1` marks the central
  corridor.  **Set by the external mesh generator and READ from the mesh NetCDF
  file** (`mod_mesh.F90:452-456`, `part_mesh.F90:709-712`) — NOT computed at runtime.
- `mesh%bctype(is,ie)`  : BC_IN (interior), BC_FAULT (fault), boundaries.

So drdg3d default == our `--mixed-flux adjacent`: central in a fault-adjacent corridor,
upwind in the bulk, Riemann+friction on the fault face.

### The three flux forms

1. CENTRAL (corridor), `mod_wave.F90:591-600`:

       fstar(1:3) = dT            ! traction jump  T_R - T_L (each side's physical T)
       fstar(4..9)= n .* dV       ! velocity jump  V_R - V_L
       fstar      = 0.5 * fstar   ! = 0.5*(F_R - F_L), the DG central contribution

   Non-dissipative.  Uses each side's *physical* traction (material-correct) but does
   NOT impedance-weight.  **Same character as our `BimaterialFlux::*Central*`**
   (`dynamic/godunov_flux_bimaterial.cpp:312-381`,
   F* = 0.5*(A_self*Q_self + A_nbr*Q_nbr)).

2. UPWIND / impedance-weighted Riemann (bulk + base of fault), `mod_wave.F90:675-688`:

       dF       = FR - FL                                   ! rotated traction/vel jumps
       alpha(1) = (dF(1) + Zp_out*dF(4)) / (Zp + Zp_out)    ! P  characteristic
       alpha(2) = (dF(2) + Zs_out*dF(9)) / (Zs + Zs_out)    ! S1 characteristic
       alpha(3) = (dF(3) + Zs_out*dF(8)) / (Zs + Zs_out)    ! S2 characteristic
       fstar(1) = alpha(1)*Zp ; fstar(4) = alpha(1)         ! traction ; velocity ...

   This is THE bimaterial treatment: the (Zs+Zs_out) denominator and Zs_out numerator
   are the impedance weighting.  **Equivalent to our Pelties bimaterial Godunov**
   (`godunov_flux_bimaterial.cpp:24-33`):
   `Q*[v] = (2 Z_L/(Z_L+Z_R)) v`, `Q*[sigma] = -(2 Z_L Z_R/(Z_L+Z_R)) v`.

3. FAULT Riemann + friction (BC_FAULT), `mod_wave.F90:697-976`:
   `eta = Zs_in*Zs_out/(Zs_in+Zs_out)` (harmonic-mean shear impedance);
   `Tau_lock = fL + fstar + Tau0` with `fstar` from form (2); LSW/RS solve
   `Vel=(|Tau_lock|-|Tau_str|)/eta`.  Mirrors our `FaultFaceFlux`.

### Homogeneous identity (why central is "upwind minus dissipation")

With Zs==Zs_out, form (2) gives `fstar(2) = 0.5*dT + 0.5*Zs*dV`, while central form (1)
gives `fstar(2) = 0.5*dT`.  So `upwind - central = 0.5*Z*dV` — exactly the dissipation
term, matching our `godunov_flux.hpp` identity `Interior - Central = 0.5*|A|*[[Q]]`.
Both codes' mixed flux drops the SAME `0.5*Z*[[v]]` dissipation in the corridor.

---

## Call graph (drdg3d critical path)

    rk_step (mod_rk.F90)
      -> wave_deriv / compute_rhs (mod_wave.F90:~40-345)
           -> get_flux(mesh,u,ie,qi,fluxs)        (mod_wave.F90:347)
                -> extract_traction_velocity (uL, each side's rho,cp,cs)   :540,578
                -> branch on flux_method / bctype / fluxtype:
                     central  0.5*(dT,dV)          :587-602   (corridor)
                     upwind   alpha impedance Riem :675-688   (bulk)
                     fault    Riemann + friction   :697-976
                -> rotate_u(invTv,invTs,fstar)      :979      (local -> global)
           -> ru += LIFT * fluxs                    :331

---

## Side-by-side

| aspect                      | drdg3d                              | SEAS-MFEM (matrix path)                         |
|-----------------------------|-------------------------------------|------------------------------------------------|
| mixed-flux structure        | central corridor + upwind bulk      | same (`--mixed-flux adjacent`)                 |
| central flux                | 0.5*(dT,dV), non-dissipative        | 0.5*(A_self Q_self + A_nbr Q_nbr), non-dissip. |
| bimaterial upwind           | alpha impedance Riemann (675-688)   | Pelties Godunov (godunov_flux_bimaterial.cpp)  |
| corridor selection          | `fluxtype` from **mesh generator**  | `central_flux_face_set_` from fault adjacency  |
| **central across a material jump** | **never** (corridor is homogeneous; shipped examples tpv5/102/104/105/stepover are homogeneous) | **yes, unguarded** (bimaterial_wave_operator.inl:515-690) |
| material contrast treatment | impedance-weighted upwind only      | impedance-weighted upwind EXCEPT we override it with central in the corridor |

---

## What is missing in our code

1. **A material-jump guard on the central corridor.** `BuildPerFaceCentralFluxMatrices_`
   loops over `central_flux_face_set_` and calls `BuildPerFaceCentralMatricesGlobal`
   for every face, never checking whether the two elements share a material
   (`bimaterial_wave_operator.inl:550-571`).  drdg3d's corridor never crosses a
   contrast, so it never hits this case.

2. **Consequence.** At the TPV31 depth steps (2400/5000/10000 m) the fault-adjacent
   corridor crosses an impedance contrast.  There we apply a flux that is BOTH
   non-dissipative AND not impedance-weighted.  On the asymmetric production mesh this
   seeds a `[[v_n]]` that LSW friction ratchets into the sigma_n collapse — observed
   only at those three depths (mid-layer stations flat).  See
   `tpv31_debug_document/` sigma-n-by-depth table.

3. **It is not a drdg3d feature we lack.** We already have the correct bimaterial
   tool (Pelties impedance Riemann == drdg3d alpha).  We are simply *not using it* on
   the contrast-crossing corridor faces; we override it with central.

---

## Recommended fix (consistent with drdg3d)

On corridor faces, choose the flux by whether the face crosses a material contrast:

- (a) HARD: face with `mu_self != mu_nbr` (or Z_self != Z_nbr) -> use the bimaterial
      upwind (Pelties) instead of central.  Equivalent to drdg3d (homogeneous corridor
      + upwind contrasts).  Smallest change; at those faces it reverts to the already-
      validated upwind.
- (b) BLEND: `F = central + theta*(bimaterial_upwind - central)` on corridor faces
      (theta in (0,1]), restoring theta of the impedance-weighted dissipation while
      keeping the low-dispersion benefit elsewhere.  Matches PLAN_upwind_unstructured
      Phase 2 (theta-blend), generalized to the bimaterial upwind.

Default must be OFF / byte-exact (theta=1 == pure upwind on those faces == today's
None there; central elsewhere unchanged).  Local guard test: bimaterial central across
an impedance contrast is NOT energy-dissipative while the upwind is.  Production gate:
Frontera A/B, TPV31 matrix on the normal mesh, central-corridor-with-guard vs without.

---

## Gotchas / non-obvious

- drdg3d `fluxtype` is a *mesh-generation artifact* (NetCDF field), not runtime logic —
  the central-corridor placement decision is made by their mesher, which is where any
  "don't mark contrast faces as central" rule would live (not visible in the Fortran).
- drdg3d central `fstar` is assembled in the GLOBAL frame (dTx,dVx) then passed through
  `rotate_u(invTv,invTs,...)` at :979 like the rotated-frame upwind `fstar`; verify the
  frame bookkeeping if porting their exact kernel (not relevant to the leak diagnosis).
- Both codes' fault FACE always uses the upwind/Riemann + friction (never central);
  central is corridor-only.

## Open questions

- Does drdg3d's mesher ever mark a contrast-crossing face as `fluxtype==1`?  Needs the
  mesh-gen tool (not in `drdg3d/src`).  Working assumption: no (homogeneous corridors).
- Final isolation of central-vs-mesh-asymmetry for our leak still wants the one Frontera
  A/B: TPV31 matrix, `--mixed-flux none`, normal mesh.
