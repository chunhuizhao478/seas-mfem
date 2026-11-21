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
- Gold file (2 steps): with one staggered pass per step
  (`mode1_gold_single_pass.toml`) every nodal value of disp_x, disp_y and d
  matches to about 1e-10 relative, and Fy = 5.690424519 and 11.38077285
  (gold 5.69042452, 11.38077285). With the input's fixed-point iterations
  (`mode1_gold.toml`, up to 20 passes to 1e-8) the differences are at most
  2e-3 relative on disp_x near zero and 8e-5 on d, inside the exodiff
  tolerance of the RACCOON test (rel 1e-4 with abs_zero 1e-9: no value fails).
  So the gold file holds one pass per step. The MOOSE run M-R1 explains it:
  MOOSE does iterate, but its fixed-point test (on the automatically scaled
  residual) is met after one pass at steps 1 and 2; pff, without the scaling,
  takes 3 passes there.
- Full run (`mode1.toml`, 175 steps to t = 3.5e-3, np = 1): 42 s. Peak
  Fy = 691.0 at u = 2.86e-3; the load drops to about 14 by u = 3.2e-3.
  52 steps after the onset reach the 20-pass limit and are accepted, as the
  input allows (`accept_on_max_fixed_point_iteration = true`).
- Against MOOSE (M-R1, `gold/moose_mode1_mac_np4_farms46ad2eb-dirty/`):
  the peak Fy = 691.00227 at t = 2.86e-3 agrees to 1e-11, and the largest
  relative difference over the curve is 2.0e-5 (step 2). Where both codes hit
  the 20-pass cap (MOOSE steps 132–156) the values agree to round-off; the
  other differences come only from the number of passes. MOOSE caps 25 steps,
  pff 55, because pff's fixed-point test uses the unscaled residual
  (automatic_scaling is not reproduced yet).
