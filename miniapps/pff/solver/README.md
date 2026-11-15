# solver (phase 2–5)

- `cfl.hpp`: λ_max(M⁻¹K) by power iteration (undrained moduli for hydromech).
- `explicit_stepper.hpp`: central difference; landing on pulse onsets k·T_w; ζ(t) windows; quench modes.
- `bound_constrained_solver.hpp`: PETSc `vinewtonrsls` (parity) and a native primal–dual active set
  with HypreBoomerAMG (HMIS coarsening, MFEM's default).
- `staggered_step.hpp`: u → p → H → d per step; optional fixed-point passes.
