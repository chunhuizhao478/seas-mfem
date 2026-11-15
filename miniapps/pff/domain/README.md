# domain (phase 2–5)

- `explicit_elasticity_operator.hpp`: lumped M, lumped dashpot C_b, F_int(u, d, p), essential DOFs
  (u_z on the end faces, the corner pin). It follows `seas` `nitsche_wave_operator.hpp`'s
  (M + Δt/2 C) A = F − K U − C U* step.
- `damage_operator.hpp`: SPD operator and RHS per step, bounds d_lb ≤ d ≤ 1.
- `pressure_operator.hpp`: backward-Euler u–p mass equation (hydromech).
