# pff architecture

pff is a standalone miniapp. It shares no code with `miniapps/seas` and
includes nothing from it; where a seas pattern is useful (directory roles,
benchmark folders, the GPU rules) pff follows the pattern with its own code.
Namespace: `mfem::pff`. The restructure proposal and the CPU/GPU plans are the
Claude Doc "pff restructure proposal"; PLAN.md is the working plan.

## Directory layout

| Directory | Status | Contents |
| --- | --- | --- |
| `common/` | started | `pff_types.hpp`: enums for decomposition, crack model, driving energy, permeability, porosity |
| `config/` | phase 1 done | `pff_config.hpp` (CaseConfig, dynamic cases), `qs_config.hpp` (QsCaseConfig, phase R), their parsers, `toml_reader.hpp` (strict TOML helpers, toml11 v3) |
| `constitutive/` | phase 1 done | `strain_split.hpp`: none / spectral / voldev σ and ψ± (MOOSE branches). The material interface is pff's own, in MFEM's style: integrators and coefficients over quadrature data, with the `MFEM_HOST_DEVICE` kernels inside |
| `phasefield/` | phase 1 done | `phase_field_functions.hpp`: g(d), α(d), c0, source, crack energy; `initial_damage.hpp`: node-interior seeding of initial_damage boxes |
| `loading/` | phase 1 done | `pulse_load.hpp` (ElkPulseLoadExperiment), rates from (t_r, t_d) |
| `domain/` | C1 done | `qs_mesh.hpp`, `qs_damage_operator.hpp` (phase R); `explicit_elasticity_operator.hpp`: lumped mass, lumped Lysmer dashpot blocks, pulse load vector, power iteration for ω_max |
| `integrator/` | started | split-elasticity residual and tangent, damage diffusion / reaction / source, dashpot and pressure traction later |
| `solver/` | C1 done | `bound_constrained_qp.hpp` (native active set), `central_difference.hpp` (half-step, dashpot at vⁿ, Kelvin–Voigt as the backward difference of the stored internal force), `cfl.hpp` |
| `kernels/` | C1 done | `split_force_kernels.hpp`: split internal force as one `mfem::forall` over E-vectors (P1 simplices) |
| `multiphysics/poroelastic/` | started | `poro_damage.hpp`: α(d), φ(d), M, c_p, H, permeability models. Pressure operator and Biot coupling later (C4) |
| `io/` | started | `gmsh_physical_names.hpp`; next: reaction force CSV, ParaView, energies, checkpoint |
| `drivers/` | started | `qs_driver.cpp` (quasi-static staggered, phase R); `pff_driver.cpp` (explicit dynamics, C1: `-fd` holds d = 0) |
| `raccoon/` | phase R | RACCOON quasi-static benchmarks: `mode1_brittle_fracture/`, `mode2_brittle_fracture/` |
| `spe_hftc/`, `cmame_revision/` | configs | `configs/*.toml` per MOOSE case; `mesh/`, `gold/`, `benchmark_data/`, `benchmark_document/` as they are added |
| `tests/unit/` | 7 files, 164 checks | `make test` |
| `tests/parallel/`, `tests/verification/` | started | `raccoon_mode1_gold.sh`, `raccoon_mode1_moose.sh`, `dashpot_reflection.cpp`, `cavity_hoop.cpp` (README lists the C1 gates) |
| `jobs/<cluster>/` | later | sbatch scripts with the commit guard of the repo CLAUDE.md |
| `archive/qs_prototype/` | frozen | earlier quasi-static AT2 prototype (does not match the benchmarks) |
| `extern/toml11/` | submodule | toml11 v3.8 (commit d4eb5f3) |

Each benchmark folder (`raccoon/<case>/`, `spe_hftc/`, `cmame_revision/`)
holds `configs/`, and as needed `mesh/`, `gold/<run>_jobid<N>_git<hash>/`,
`benchmark_data/` and `benchmark_document/`.

## Per-step data flow

Quasi-static (phase R, and the static initial state of cmame):

```
load step n --> repeat until the relative change < tol (at most max_passes):
                   u  = argmin E(u; d)            linear (no split) or Newton (spectral)
                   psi_active(u) at quadrature points (or element averages, as MOOSE transfers)
                   d  = bound-constrained solve, d_old <= d <= 1
            --> reaction force, output
```

Dynamic (phases C1–C5):

```
t_n --> u^{n+1} = central difference with F_int(u^n, d^n) - alpha(d^n) p^n, pulse p(t), dashpots, quench
    --> p^{n+1}  backward Euler with eps_v(u^{n+1}), k(d^n, eps), M(d^n)        [hydromech]
    --> psi+ (u^{n+1});  H = max(H, <2 psi+ + c_p p^2>_+)                        [hydromech]
    --> d^{n+1}  bound-constrained SPD solve, d_old <= d <= 1
    --> energies, metrics, output
```

## CPU and GPU

One source, `--device cpu|cuda`. Per-element and per-quadrature-point loops
are `mfem::forall` over `Read()`/`Write()` views with no virtual dispatch,
STL or `DenseMatrix` calls inside. A path without a device kernel calls
`MFEM_ABORT`; it never falls back to the host silently.
