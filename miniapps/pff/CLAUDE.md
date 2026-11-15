# pff miniapp: rules for editing

The repository CLAUDE.md applies (GitHub first, commit guard in every job
script, debugging behaviour). In addition, for pff:

## Source of truth

- The MOOSE/FARMS inputs and C++ on farms branch `cdms` define the models.
  Local copy: `/Users/chunhuizhao/projects/farms_cdms/pulsepower/{spe_hftc,cmame_revision}`.
  Every case TOML names the input it reproduces, and every ported function
  names the MOOSE object in its header comment.
- Port formulas exactly, including clamps, floors and defaults (for example
  g_min = 1e-6 in the Biot coefficient, and the 1e-20 floor of 1/M). A
  difference you believe is a MOOSE bug gets reported and discussed, not
  silently fixed. See "Parity details" in PLAN.md.
- Phase R (quasi-static benchmarks) reproduces RACCOON at the commit pinned
  by the farms@cdms `raccoon` submodule (hugary1995/raccoon 9fd3a6f,
  `tutorials/`). Same rule: inputs ported exactly, the input cited.

## Standalone

- pff shares no code with `miniapps/seas` and includes nothing from it.
  Follow a seas pattern by writing it in pff.
- Hydromech code lives in `multiphysics/poroelastic/`.

## Numerics that must not regress

- Damage irreversibility comes from bounds (d_old ≤ d ≤ 1), not from
  projection after an unconstrained solve (the archived prototype's defect).
- The damage diffusion coefficient is 2 Gc l / c0 (ADPFFDiffusion).
- λ = K − 2G/3 in 2D as well (plane strain, MOOSE LIBMESH_DIM = 3).
- The stiffness-proportional quench cannot be used with explicit steps at
  production Δt (about 600× smaller stable step). Pick the quench mode
  explicitly in the case file.

## Tests

- `make test` must pass before every push. New pointwise physics needs a
  unit test against a value taken from a MOOSE input, a MOOSE output or a
  closed form, with the source cited in the test.
- Verification against MOOSE runs cites the run (cluster, job id, commit).
