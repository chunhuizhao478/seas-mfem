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

Each phase lists its directory and its gate.

1. **Pointwise kernels** (done in this commit): `materials/`, `loading/`,
   `config/`, `io/gmsh_physical_names.hpp`, `tests/unit/`. Gate: `make test`.
   The split satisfies σ = ∂ψ/∂ε to 1e-9 on 1000 random strains. The pulse
   peak matches the write-up's 206/215/223 MPa. The rates match the input
   comments. The poro laws match the Elk formulas.
2. **Explicit elastodynamics, no damage**: `integrators/split_elasticity_integrator.hpp`,
   `integrators/dashpot_integrator.hpp`, `domain/explicit_elasticity_operator.hpp`,
   `solver/cfl.hpp`, `drivers/pff_driver.cpp`. Gates:
    - With ζ = 0 and free faces, energy conserves to 0.1 %.
    - The early-time hoop stress matches the pressurised cavity (spe_hftc PLAN A.8).
    - Dashpot reflection is under 2 %.
    - v_max(t) of the 2D section matches MOOSE with d frozen at 0.
3. **Damage**: `domain/damage_operator.hpp`, `solver/bound_constrained_solver.hpp`.
   Gates:
    - The RACCOON quasi-static gold cases on the branch (`cdm_examples/pf_raccoon_examples`) reproduce.
    - spe_hftc 2D, one pulse to 10 µs: η = 1 gives d > 0.5 out to 1.76 mm and 77 J/m; η = 0.3 gives d_max = 0.18.
    - cmame 2D solid benchmark: the first pulse matches the MOOSE extent and energies.
4. **Pulse trains and the quench**: the 20-pulse spe_hftc 2D η sweep (Vista
   results) and 100 pulses of the cmame 2D solid benchmark. Gate: extents per
   pulse within one element, and the energy budget closes.
5. **Hydromech**: `domain/pressure_operator.hpp`, `materials/poro_damage.hpp`
   (done), the static initial state (one implicit solve with confinement and
   p0, as `static_solve.i`). Gates:
    - Terzaghi and Mandel (analytic).
    - The undrained wave speed in a 1D column.
    - cmame 2D hydromech benchmark: pp(t) at the hole and the damage pattern
      against MOOSE; the fluid energy terms of the MOOSE CSV.
6. **3D and performance**: per-rank mesh files or parallel refinement,
   HMIS AMG (MFEM's default), skipping the damage solve when ψ⁺ ≤ ψ⁺ of the
   last solve everywhere, and timing against the Vista/LS6 numbers. Gates:
    - spe_hftc cube and cmame 3D em_5J, first pulses against MOOSE.
    - A wall-time table.
7. **Optional**: GPU kernels (partial assembly of the split force), AMR,
   Weibull/Gaussian Gc fields, the other permeability and porosity variants of
   `cmame_revision/2d_hydromech/permeability_formula/`.

## Parity details to keep (not to "fix")

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
