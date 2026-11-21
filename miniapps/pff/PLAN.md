# pff port plan: farms_cdms phase-field models to MFEM

Goal: reproduce two MOOSE/FARMS model families of farms branch `cdms` in MFEM,
faster than MOOSE, with every phase gated against an analytic result or an
existing MOOSE run.

| Model | MOOSE inputs (farms@cdms) | Physics |
| --- | --- | --- |
| spe_hftc | `pulsepower/spe_hftc/code_{2d,3d}/elasticity_absorbing.i` + `fracture.i` | Dry acrylic cube with a borehole, 20 pulses of 144 J. Elastodynamics with Kelvin–Voigt damping, AT1 phase field, spectral split, Lysmer dashpots, time-switched quench |
| cmame_revision | `pulsepower/cmame_revision/{2d,3d}_{solid,hydromech}/…/elasticity*.i` + `fracture*.i` | Granite disc (2D) or pierced cylinder (3D), 5 J pulses, 1 MPa confinement. Solid: as above, without quench. Hydromech: u–p Biot with damage-dependent α, φ, M and permeability, plus an enhanced history energy H |

The longer write-up (formulation, MOOSE cost profile, MFEM mapping, risks) is
the Claude Doc "farms_cdms phase field: MOOSE formulation and MFEM port plan".
This file is the working plan.

## Discretisation (decided in the write-up, revisit at the gates)

- **u**: H1 P1 (vector), lumped mass, central difference (Newmark β = 0,
  γ = 1/2), as in `seas` `nitsche_wave_operator.hpp`. Dashpots enter as
  (M + Δt/2 C_b) with a lumped, diagonal C_b. CFL from a power iteration on
  M⁻¹K. Hydromech uses the undrained P-wave speed √((K + α²M + 4G/3)/ρ).
- **d**: H1 P1, a linear SPD operator per step (diffusion 2 Gc l/c0, and
  reaction g''·D for AT1 with quadratic g), solved as a bound-constrained
  problem d_lb ≤ d ≤ 1. The parity path is PETSc `vinewtonrsls`; a native
  active-set solver with HypreBoomerAMG comes later. d_lb = d_old, and is 0.9
  on seeded boxes at the first step.
- **p** (hydromech): H1 P1, backward Euler, the same mass equation as
  `ElkPorousFlowFullySaturatedMassTimeDerivative`:
  (p − p_old)/(Δt M) + α(ε_v − ε_v,old)/Δt − ∇·(k/μ ∇p) = 0, with ε_v taken
  from u at n+1. SPD; PCG + BoomerAMG.
- **Coupling**: single pass per step: u(n+1) from d(n) and p(n), then p(n+1),
  then H and d(n+1). A fixed-point option exists for the parity study.
- **Quench (spe_hftc)**: ζ_quench = 1e-5 s stiffness damping limits the
  explicit Δt about 600×. The options are `quench = "mass"` (C = α_M M) or
  `"implicit"` (Newmark on the ring-down only), chosen at gate 4 on 2D
  evidence. `"stiffness"` stays available for parity at a small Δt.

## Phases

Order: R first (quasi-static benchmarks), then the CPU phases C1–C5, each
followed by its GPU phase G1–G5 on the same source (`--device cuda`, Vista
GH200). The CPU result is the reference of each GPU gate.

0. **Pointwise kernels** (done): `constitutive/`, `phasefield/`, `loading/`,
   `config/`, `multiphysics/poroelastic/poro_damage.hpp`,
   `io/gmsh_physical_names.hpp`, `tests/unit/`. Gate: `make test`. The split
   satisfies σ = ∂ψ/∂ε to 1e-9 on 1000 random strains; the pulse peak matches
   206/215/223 MPa; the poro laws match the Elk formulas.
   **C0** (done): layout restructure; the same 98 checks pass.
1. **R: RACCOON quasi-static benchmarks** (`raccoon/`), RACCOON 9fd3a6f:
    - R1 `tutorials/mode1_brittle_fracture`: AT2, no split, 30×15 QUAD4 with
      two initial refinements in [0.4,1]×[0,0.05]. Gates: 888 elements and
      994 nodes as in the gold (passes, same coordinates); the 2 gold steps
      (Fy = 5.69042, 11.38077) and fields within relative 1e-4 (passes:
      `tests/verification/raccoon_mode1_gold.sh`); the full Fy–u curve
      against MOOSE M-R1 within 1e-4 (passes: 1.6e-5, peak equal to 1e-11,
      `tests/verification/raccoon_mode1_moose.sh`).
    - R2 `tutorials/mode2_brittle_fracture`: spectral split (Newton,
      analytic tangent checked by finite differences), slit mesh of two
      stitched halves. Gates: 2826 elements and 3026 nodes; Fx = 0.923682,
      1.847363 at the 2 gold steps; the full curve against MOOSE.
    - R3 three-point bending: deferred (RACCOON's case is finite-strain J2
      elastoplastic), revisited after R1 and R2.
2. **C1/G1 explicit elastodynamics, no damage**: force kernel, lumped mass,
   dashpots, central difference, CFL. Gates: energy conserved to 0.1 % with
   ζ = 0; the cavity hoop stress (spe_hftc PLAN A.8); dashpot reflection
   under 2 %; v_max(t) of the 2D section against MOOSE with d = 0. G1: one
   step matches CPU to 1e-12, no state-sized copies inside a step.
3. **C2/G2 damage**: native active-set solver against PETSc VI to 1e-8;
   spe_hftc 2D one pulse to 10 µs (η = 1: d > 0.5 out to 1.76 mm and
   77 J/m; η = 0.3: d_max = 0.18); cmame 2D solid first pulse. G2: d matches
   the CPU solver to 1e-8.
4. **C3/G3 pulse trains and the quench**: 20-pulse spe_hftc 2D η sweep and
   100 pulses of the cmame 2D solid benchmark; extents per pulse within one
   element, the energy budget closes.
5. **C4/G4 hydromech** (`multiphysics/poroelastic/`): Terzaghi, Mandel, the
   undrained column speed; cmame 2D hydromech pp(t) and damage pattern; the
   static initial state uses the phase-R quasi-static solver.
6. **C5/G5 3D and timing**: spe_hftc cube and cmame 3D em_5J first pulses
   against MOOSE; a wall-time table against Vista/LS6 MOOSE and between CPU
   and GPU.
7. **Optional**: AMR, Weibull/Gaussian Gc fields, the other permeability and
   porosity variants.

## Parity details to keep (not to "fix")

- RACCOON transfers `psie_active` to the fracture app as a CONSTANT
  MONOMIAL, so the damage source sees the element average of ψ_active, not
  its quadrature-point values. Phase R reproduces this.
- MOOSE runs the RACCOON inputs with `automatic_scaling = true`; pff
  reproduces it in the fixed-point test (`[solver] automatic_scaling`). With
  it the mode I gold file (one pass at steps 1 and 2) matches to 1e-10, and
  the pass counts of M-R1 match on 143 of 175 steps; the rest point to a
  constant factor of about 3/4 in MOOSE's norm, still open (raccoon/README.md).

- The driving-energy conventions differ. spe_hftc: ψ = αGc/(c0 l) + g ψ⁺ with
  no history. cmame hydromech: ψ = αGc/(c0 l) + ½ g H with
  H = max_t⟨2ψ⁺ + c_p p²⟩₊. For a dry history both reduce to ψ⁺.
- λ = K − 2G/3 in 2D too (MOOSE `LIBMESH_DIM` = 3): the 2D cases are plane strain.
- Seeded damage: MOOSE writes 0.9 on the box block, then 0 on the other
  blocks, so nodes on the box boundary end at 0. Reproduce node-interior seeding.
- The MOOSE inputs' "#137 MPa" comment next to ft is stale. The fparse value
  for E = 50 GPa, Gc = 100 and l = 0.2 mm is 96.8 MPa; w_c = Gc/ft uses 96.8 MPa.
- cmame 3D hydromech uses HHT α = 0.1. The explicit scheme has no HHT
  dissipation, so compare at matched resolution and record the difference.
