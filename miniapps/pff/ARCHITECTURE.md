# pff architecture

The layout mirrors `miniapps/seas` (see its ARCHITECTURE.md on the active
branches): config, then pointwise physics, operators, the time stepper,
I/O, drivers, and tests split into unit, parallel and verification. The
cases and cluster jobs live in the repo, so a run is reproducible from a
commit. Namespace: `mfem::pff`.

| Module | Directory | Status | Contents |
| --- | --- | --- | --- |
| Common | `common/` | started | `pff_types.hpp`: enums for decomposition, crack model, driving energy, permeability, porosity |
| Config | `config/` | done (phase 1) | `pff_config.hpp` (CaseConfig), `pff_config_parser.hpp` (strict TOML, toml11 v3) |
| Materials | `materials/` | done (phase 1) | `phase_field_functions.hpp` (g, α, c0, source, crack energy), `strain_split.hpp` (none/spectral/voldev σ, ψ±), `poro_damage.hpp` (α(d), φ(d), M, H, permeability models) |
| Loading | `loading/` | done (phase 1) | `pulse_load.hpp` (ElkPulseLoadExperiment), rates from (t_r, t_d) |
| Integrators | `integrators/` | phase 2–5 | split-elasticity internal force (element loop, no tangent), lumped dashpot, pressure traction, damage reaction/source, Darcy and Biot coupling |
| Domain | `domain/` | phase 2–5 | explicit elasticity operator (M, C_b, F_int(u, d, p)), damage operator (SPD + bounds), pressure operator (backward Euler) |
| Solver | `solver/` | phase 2–5 | CFL estimate, central-difference stepper with the pulse-window schedule, bound-constrained damage solver (PETSc VI and native active set), staggered step |
| I/O | `io/` | started | `gmsh_physical_names.hpp` (names to attributes); next: energy budget, fracture metrics, ParaView output, checkpoint/restart |
| Drivers | `drivers/` | phase 2 | `pff_driver.cpp -c cases/<model>/<case>.toml` |
| Tests | `tests/unit/` | 5 files, 98 checks | `make test` |
|  | `tests/parallel/`, `tests/verification/` | phase 2+ | MPI tests; analytic and MOOSE-comparison benchmarks |
| Cases | `cases/spe_hftc/`, `cases/cmame_revision/` | 6 files | one TOML per MOOSE case |
| Jobs | `jobs/<cluster>/` | phase 4+ | sbatch scripts with the commit guard of the repo CLAUDE.md |
| Archive | `archive/qs_prototype/` | frozen | earlier quasi-static AT2 prototype |
| Extern | `extern/toml11` | submodule | same commit as `miniapps/seas/extern/toml11` |

## Per-step data flow (planned, phase 2–5)

```
t_n --> u^{n+1} = central difference with F_int(u^n, d^n) - alpha(d^n) p^n, pulse p(t), dashpots, quench
    --> p^{n+1}  backward Euler with eps_v(u^{n+1}), k(d^n, eps), M(d^n)        [hydromech]
    --> psi+ (u^{n+1});  H = max(H, <2 psi+ + c_p p^2>_+)                        [hydromech]
    --> d^{n+1}  bound-constrained SPD solve, d_old <= d <= 1 (skipped when no point exceeds its last driving energy)
    --> energies, metrics, output
```

All pointwise functions are `MFEM_HOST_DEVICE`, so the element kernels can
move to partial assembly on GPUs without changing the physics code.
