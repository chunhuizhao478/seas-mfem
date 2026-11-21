# raccoon: quasi-static benchmarks (phase R)

The RACCOON tutorials at hugary1995/raccoon 9fd3a6f (the `raccoon` submodule
of farms@cdms; identical copies in farms `cdm_examples/pf_raccoon_examples`),
run with `qs_driver`:

```bash
make qs_driver
./qs_driver -c raccoon/mode1_brittle_fracture/configs/mode1.toml      # full curve, out/mode1/
tests/verification/raccoon_mode1_gold.sh                              # gold check (numpy, netCDF4)
tests/verification/raccoon_mode1_moose.sh                             # full curve against MOOSE M-R1
```

| Case | RACCOON input | Status |
| --- | --- | --- |
| `mode1_brittle_fracture/` | `tutorials/mode1_brittle_fracture` | R1 done: gold check and full curve against MOOSE pass |
| `mode2_brittle_fracture/` | `tutorials/mode2_brittle_fracture` | R2: next (slit mesh, spectral tangent) |
| three-point bending | `tutorials/three_point_bending` | deferred (finite-strain J2 elastoplastic) |

Each case folder: `configs/` (one TOML per run, header cites the input),
`benchmark_data/raccoon_9fd3a6f/` (RACCOON's gold exodus, LGPL-2.1, copied
unchanged), and later `gold/<run>_jobid<N>_git<hash>/` for full MOOSE runs.

## mode1_brittle_fracture (R1)

- Mesh: the 30×15 QUAD4 mesh with two initial refinement passes reproduces the
  gold mesh exactly (888 elements, 994 nodes, same coordinates).
- Gold file (2 steps): every nodal value of disp_x, disp_y and d matches to
  about 1e-10 relative, and Fy = 5.690424519 and 11.38077285 (gold
  5.69042452, 11.38077285), with the input's fixed-point settings
  (`mode1_gold.toml`) and with one pass per step (`mode1_gold_single_pass.toml`).
  The fixed-point test uses MOOSE's automatic scaling
  (`[solver] automatic_scaling = true`, as in elasticity.i): one factor per
  displacement component, 1/max|J_ii| of the unconstrained volumetric
  Jacobian at the first solve (here 1/((λ+3μ)·4/3) = 2.0635e-6), applied to
  the kernel rows, not to the Dirichlet rows. With it the test passes after
  one pass at steps 1 and 2, as in MOOSE, which is why the gold file looks
  like a single pass. Without the scaling pff took 3 passes there.
- Full run (`mode1.toml`, 175 steps to t = 3.5e-3, np = 1): about 45 s. Peak
  Fy = 691.00227 at u = 2.86e-3; the load drops to about 14 by u = 3.14e-3.
- Against MOOSE (M-R1, `gold/moose_mode1_mac_np4_farms46ad2eb-dirty/`):
  the peak agrees to 1e-11 and the largest relative difference over the curve
  is 1.6e-5. The pass counts per step agree on 143 of 175 steps (MOOSE caps
  25 steps at 20 passes, pff the same 25). On the other 32, all at the
  boundary between two counts, pff takes one pass more. The pff/MOOSE pass
  counts are consistent only with MOOSE's fixed-point norms being 0.749–0.753
  times pff's on every step, so a constant factor of about 3/4 is left in the
  scaling or the norm. Open: MOOSE's printed scaling factors (verbose run
  requested) decide it.
