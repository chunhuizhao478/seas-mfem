# integrators (phase 2–5)

- `split_elasticity_integrator.hpp`: `NonlinearFormIntegrator::AssembleElementVector`.
  ∫ σ(ε(u), d) : ∇w with `materials/strain_split.hpp`; g(d) at the quadrature point; −α(d) p I in hydromech.
- `dashpot_integrator.hpp`: boundary matrix ρ[c_p n⊗n + c_s(I − n⊗n)], row-sum lumped (`FarmsNonReflectDashpotBC`).
- pressure traction: `VectorBoundaryFluxLFIntegrator` with the `loading/pulse_load.hpp` coefficient (−p n).
- damage: diffusion 2 Gc l/c0 and the linearised source of `PhaseFieldModel::Source`.
- hydromech: Darcy k(d, ε)/μ (`materials/poro_damage.hpp`), storage 1/M, Biot coupling α ε̇_v.
