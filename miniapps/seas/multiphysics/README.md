# multiphysics

Physics that couples another field to the elastic bulk of the Nitsche
dynamic-rupture solver, one subfolder per physics. Each is an option that is
off by default; with it off every existing run is bit-identical.

| Folder | Physics | Option | Plan |
|--------|---------|--------|------|
| `poroelastic/` | Biot poroelasticity, the relative fluid displacement `w` in the same DG space as `u` | `[multiphysics.poroelastic] enabled = true`, or `--poroelastic 1` | `document/poroelastic_dev/poro-plan-rev3.pdf`, equations in `poro-derivation-rev1.2.pdf` |

Inside a physics folder the files follow the plan's step list: the
configuration first (`poro_config.hpp`), then the material sample, the
operators, the stepper and the fault coupling, and `tests/` with one test
file per phase (`make test-poro-p0` and so on).

The papers and reference data for the poroelastic benchmarks stay in
`miniapps/seas/poroelastic/benchmark_document/` and `benchmark_data/`.
