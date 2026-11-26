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

- **u**: H1 P1 (vector), consistent mass (lumped optional), central difference (Newmark β = 0,
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
  **Decided 2026-11-22 (Andy):** `quench = "mass"` with α_M = 2e5 /s after
  t_free (decay time 1/α_M = 10 µs, the 1/ζ_quench design target), and
  ζ_pulse = 0. The MOOSE runs use ζ_pulse = 1e-8 s, which removes about 12 %
  of the input energy during the pulse in the MHz band (pff kv_work 147 J/m of
  1190 at 10 µs, η = 1; MOOSE's budget does not record it) and lowers the
  damage (η = 1, 10 µs: crack energy 76.6 → 97.5 J/m, d ≥ 0.5 extent 1.88 →
  3.24 mm). Parity for ζ_pulse = 0 is M-C2c (MOOSE M-C2a with the CLI override
  zeta_pulse=0); the production configs no longer reproduce the published
  MOOSE runs step for step.
  **Mass matrix (2026-11-22):** at ζ_pulse = 0 the lumped mass grows the
  crack faster than MOOSE's consistent mass after 6 µs (η = 1, 10 µs: crack
  energy 97.5 vs 82.5 J/m in M-C2c, 79.7 in M-C2d at dt ≤ 3e-8). pff with
  `-cm` (consistent mass, PCG per step, +23 % per step in 2D) gives 80.35 J/m,
  d ≥ 0.5 extent 2.41 mm (M-C2d 2.17, M-C2c 2.45). Not the time step (pff
  cfl 0.2: 96.1 J/m; MOOSE dt 1e-7 → 3e-8 moves it 3 %). With ζ_pulse = 1e-8
  the high frequencies are damped and the two masses agree. The result is
  sensitive to the mass discretization, i.e. the mesh (h = 0.3 mm, l = 0.6 mm)
  is not converged for undamped crack growth.
  With the CFL from the consistent mass (power iteration on M_c^{-1} K,
  dt_crit 5.84e-8 vs 9.53e-8 lumped) pff -cm gives at 10 µs crack 79.57 J/m,
  d_max 0.988, extents 2.170 / 1.358 mm vs M-C2d 79.66, 0.983, 2.170 /
  1.378 mm; 2x the lumped wall time. MOOSE single pass (M-C2e, 76.2 J/m)
  rules out the coupling. **Decided (Andy, 2026-11-22): `[time] mass =
  "consistent"` (the default) for production; the C1/C2 verification cases keep
  `mass = "lumped"` as validated.**

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
      against MOOSE M-R1 (passes: 4.4e-12 with equal pass counts, gate 1e-9,
      `tests/verification/raccoon_mode1_moose.sh`).
    - R2 `tutorials/mode2_brittle_fracture`: spectral split (Newton,
      analytic tangent checked by finite differences), slit mesh of two
      stitched halves. Gates: 2826 elements and 3026 nodes; Fx = 0.923682,
      1.847363 at the 2 gold steps; the full curve against MOOSE.
      Status (2026-11-21): passes. Mesh 2826/3026 with the 28 slit pairs;
      gold steps to 5e-12; full curve against M-R2 to 6e-12 with the same
      pass count on all 1000 steps (548 capped at 20 in both). Needs the
      MOOSE Newton norm (scaled, NodalBC rows) and preset = false.
    - R3 three-point bending: deferred (RACCOON's case is finite-strain J2
      elastoplastic), revisited after R1 and R2.
2. **C1/G1 explicit elastodynamics, no damage**: force kernel, lumped mass,
   dashpots, central difference, CFL. Gates: energy conserved to 0.1 % with
   ζ = 0; the cavity hoop stress (spe_hftc PLAN A.8); dashpot reflection
   under 2 %; v_max(t) of the 2D section against MOOSE with d = 0. G1: one
   step matches CPU to 1e-12, no state-sized copies inside a step.
   Status (2026-11-21): code done (`pff_driver -fd`). Energy 0.048 % over
   10 021 steps (free faces, ζ = 0, 7 pulses); cavity: explicit vs static FE
   0.35 %, static FE vs Lamé 3.5 % at h = 0.3 mm and 1.74 % at 0.15 mm (P1
   stresses, first order); dashpot reflection 0.24 % (P and S). Against
   M-C1: KE and SE within 1 % after 1 us, v_max 13 % higher in the first
   microsecond and 0.4 % at 25 us (spe_hftc/gold/moose_c1_section_d0_*). The kernel force equals SplitElasticityIntegrator to 1e-10
   for the spectral split with damage (unit test).
3. **C2/G2 damage**: native active-set solver against PETSc VI to 1e-8;
   spe_hftc 2D one pulse to 10 µs (η = 1: d > 0.5 out to 1.76 mm and
   77 J/m; η = 0.3: d_max = 0.18); cmame 2D solid first pulse. G2: d matches
   the CPU solver to 1e-8.
   Status (2026-11-21): code done (`pff_driver` without `-fd`): one damage
   solve per explicit step on ψ⁺(u^{n+1}) (element values), bounds
   max(d_old, seed) ≤ d ≤ 1 (farms VariableConstantIrreversibleBounds is
   initial_damage_aux at step 1 and d_old after, so MOOSE damage is
   irreversible), fracture.i tolerances with its automatic scaling; a step
   whose start point already meets the absolute tolerance is skipped, as
   PETSc stops at iteration 0. Energy budget with the crack energy closes to
   0.13 J of 1190 J/m (η = 1, 10 µs); halving dt changes the crack energy by
   0.1 %. The gate numbers above come from the 2026-10-21 MOOSE runs
   (farms local_runs/README.md, 3cea9bdc: free faces, ζ = 0, mesh made then,
   8 ranks). pff at ζ = 0: η = 1 d > 0.5 to 3.24 mm, 97.5 J/m; η = 0.3
   d_max = 0.25, 6.8 J/m, i.e. more damage than that run. pff at the current
   ζ_pulse = 1e-8: η = 1 1.88 mm, 76.6 J/m; η = 0.3 d_max = 0.166, 5.3 J/m.
   Against M-C2a/b (current inputs, ζ = 1e-8, the pff mesh, Mac 2026-11-21):
   η = 1 crack energy 76.57 vs 76.58 J/m, d ≥ 0.9 to 0.609 mm in both, the
   d ≥ 0.5 front on the same node through 9.5 us (then one node further in
   pff: 1.88 vs 1.65 mm), d_max 0.949 vs 0.951; η = 0.3 crack 5.336 vs
   5.340 J/m, d_max 0.1656 vs 0.1657 (spe_hftc/gold/moose_c2_section_*).
   These replace the gate numbers above.
4. **C3/G3 pulse trains and the quench**: 20-pulse spe_hftc 2D η sweep and
   100 pulses of the cmame 2D solid benchmark; extents per pulse within one
   element, the energy budget closes.
   Status (2026-11-26): done. Andy 2026-11-24: bit-for-bit agreement with
   MOOSE is not required; the gate is pff reliability (convergence in dt,
   energy closure), MOOSE is a sanity check. All runs: ζ_pulse = 0, mass
   quench α = 2e5 /s after t_free = 25 µs (spe_hftc; cmame has no quench, as
   in MOOSE), consistent mass, cfl 0.4, single pass. energy_error is the
   budget residual, quoted against W_ext.
   - spe_hftc 2D section, 20 pulses (extents are d ≥ 0.5 from the wall):
     | η | pulse 1 crack (J/m) | pulse 20 crack (J/m) | extent (mm) | energy_error |
     |---|---|---|---|---|
     | 0.3 | 5.54 | 17.4 | 3.23 (first d ≥ 0.5 at pulse 19) | 0.005 % |
     | 0.45 | (d_max 0.43, no crack) | 887.9 | equal at cfl 0.4 and 0.2 (gold README) | −38.8 J/m; cfl 0.2: 895.6, −10.5 J/m |
     | 1.0 | 80.28 | 2660 | 32.28 = corner, reached at pulse 6 | 0.21 % |
     η 0.45: cfl 0.4 vs 0.2 differ 0.9 % in crack energy with equal extents
     and equal crack counts every pulse; the error falls at second order.
     MOOSE (Vista 1162848 + 1164181, recover-corrected) gives 703.8 J/m with
     a different crack pattern (6/3 cracks at 8/15 mm against 4/5), a
     bifurcation and not a pff time-step error
     (spe_hftc/gold/moose_c3_section_eta0p45_zetapulse0_*). η 1.0 pulse 1
     against M-C2d: 80.28 vs 79.66 J/m, 2.17 mm in both. No MOOSE 20-pulse
     runs exist at η 0.3 or 1.0 with ζ_pulse = 0. Runs: Mac, 2 ranks, commit
     85ac281, 3055 s (η 0.3) and 6015 s (η 1.0); configs =
     section_homog_eta0p45_absorb.toml with only [pulse] efficiency, [case]
     name and the mesh path changed.
   - cmame 2D solid benchmark, 100 pulses (cloud, 4 ranks, commit 85ac281,
     8074 s): pulse 1 crack 1.882 J/m and front at 5.94 mm (M-C2f 1.883,
     5.94); the d ≥ 0.5 front reaches the outer edge (12.5 mm) at pulse 5;
     crack energy less the 0.2645 J/m seed 5.40 / 17.6 / 51.3 J/m at pulses
     10 / 50 / 100; energy_error < 0.003 % through pulse 40 and 0.31 % at
     pulse 100. No MOOSE output beyond pulse 1 is reachable (Frontera).
   - The energy residual grows with the amount of damage: η 1.0 0.04 % at
     pulse 2, 0.14 % at pulse 5 and 0.21 % at pulse 20; η 0.45 0.27 % at
     cfl 0.4 and 0.07 % at cfl 0.2 (it falls with the step); cmame 0.3 % only
     after pulse 50, when the disc is fragmented.
5. **C4/G4 hydromech** (`multiphysics/poroelastic/`): Terzaghi, Mandel, the
   undrained column speed; cmame 2D hydromech pp(t) and damage pattern; the
   static initial state uses the phase-R quasi-static solver.
6. **C5/G5 3D and timing**: spe_hftc cube and cmame 3D em_5J first pulses
   against MOOSE; a wall-time table against Vista/LS6 MOOSE and between CPU
   and GPU.
   Status (2026-11-26), spe_hftc cube only: cube_d7p94_mid η 0.45, 20 pulses
   on Vista (job 1168452, commit 85ac281, 4 gg nodes, 576 ranks, 14.6 M
   tets): 2 h 37 min, 73,860 steps, 0.126 s/step; crack 4.917 J, d ≥ 0.5
   extent 11.01 mm, energy_error 0.016 %. Large meshes are partitioned once
   (`pff_partition`, `pff_driver -pm`). The MOOSE run of the same case takes
   days on 16 nodes (owned by the 3D SPE thread). cmame 3D and GPU: open.
7. **Optional**: AMR, Weibull/Gaussian Gc fields, the other permeability and
   porosity variants.

## Parity details to keep (not to "fix")

- RACCOON transfers `psie_active` to the fracture app as a CONSTANT
  MONOMIAL, so the damage source sees the element average of ψ_active, not
  its quadrature-point values. Phase R reproduces this.
- MOOSE runs the RACCOON inputs with `automatic_scaling = true`; pff
  reproduces it in the fixed-point test (`[solver] automatic_scaling`),
  with the factor from the constrained (hanging-node) Jacobian diagonal as
  libMesh assembles it. The mode I gold file matches to 1e-10 and M-R1 to
  4.4e-12 with the same pass count on every step.

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
