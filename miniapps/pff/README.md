# pff: dynamic phase-field fracture miniapp

MFEM port of the pulse-power phase-field fracture models of FARMS (MOOSE app,
`chunhuizhao478/farms`, branch `cdms`):

- **spe_hftc**: acrylic cube under a 20-pulse train (dry, AT1, spectral split,
  absorbing faces).
- **cmame_revision**: granite disc/cylinder under 5 J pulses, solid and u–p
  hydromechanical variants.

Status and next steps: [PLAN.md](PLAN.md). Module layout: [ARCHITECTURE.md](ARCHITECTURE.md).
Constraints for anyone editing the code: [CLAUDE.md](CLAUDE.md).

## Build and test

```bash
git submodule update --init miniapps/pff/extern/toml11
cd miniapps/pff
make test                                   # MFEM built in-tree (../../)
make MFEM_BUILD_DIR=/path/to/build test     # or an out-of-tree MFEM build
```

`make test` builds and runs the unit tests in `tests/unit/` (pointwise
kernels and case files, no mesh needed).

## Cases

`<model>/configs/*.toml`, one file per MOOSE case. Each file's header names the
MOOSE input it reproduces, and the values are copied from it. See
[spe_hftc/README.md](spe_hftc/README.md) and [cmame_revision/README.md](cmame_revision/README.md).

The quasi-static AT2 prototype that used to live here is in
`archive/qs_prototype/`, kept for reference. Its `CONVERGENCE_DEBUG.md` lists
defects the new code does not inherit.
