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
| `mode2_brittle_fracture/` | `tutorials/mode2_brittle_fracture` | R2 done: mesh, gold steps and full curve against MOOSE pass |
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
  displacement component, 1/max|J_ii| of the volumetric Jacobian at the
  first solve, taken from the constrained (hanging-node) system as libMesh
  assembles it: 1.547619e-6 here, equal to the factor MOOSE prints in the
  verbose run M-R1b (the unconstrained element diagonal gives 2.0635e-6,
  4/3 of it). It is applied to the kernel rows, not to the Dirichlet rows.
  With it the test passes after one pass at steps 1 and 2, as in MOOSE,
  which is why the gold file looks
  like a single pass. Without the scaling pff took 3 passes there.
- Full run (`mode1.toml`, 175 steps to t = 3.5e-3, np = 1): about 45 s. Peak
  Fy = 691.00227 at u = 2.86e-3; the load drops to about 14 by u = 3.14e-3.
- Against MOOSE (M-R1, `gold/moose_mode1_mac_np4_farms46ad2eb-dirty/`):
  the curve agrees to 4.4e-12 relative (gate 1e-9), the peak 691.00227 at
  u = 2.86e-3 in both, and the fixed-point pass count is the same on every
  step (1, 1, 2, ..., the cap of 20 on steps 132-156, then 8, 5, 5, 5, ...).
  The earlier 1.6e-5 and the 32 steps with one pass more came from the
  unconstrained scaling diagonal (a factor 3/4 in the norm), now fixed.

## mode2_brittle_fracture (R2)

- Mesh (`[mesh] type = "slit"`, domain/qs_mesh.hpp MakeSlitMesh): two 30×15
  QUAD4 halves sharing the y = 0 vertices with x ≥ 0.5
  (StitchedMeshGenerator), two OrientedBoxMarker refinements: 2826 elements,
  3026 nodes and 28 duplicated node pairs on the slit (x ≤ 0.4917), as in
  M-R2 and RACCOON's gold file.
- Spectral split with the consistent tangent (SplitTangentApply, checked
  against finite differences pointwise and per element).
- Newton as MOOSE runs it: scaled residual norm with the NodalBC rows u − g,
  and `preset = false` on xdisp, so the first Newton step also moves the
  boundary. With the preset approach pff stopped one Newton iterate earlier
  and step 2 differed by 2e-8.
- Gold steps (`mode2_gold.toml`): Fx = 0.923681826558, 1.84736335602 against
  0.9236818265627, 1.8473633560272; max |u_y| and max d agree to the printed
  digits. A nodal comparison needs the gold exodus (1.84 MB, on the Mac,
  sha256 a76cd020…), not yet in the repo.
- Full run (`mode2.toml`, 1000 steps, np 4, 14.5 min): against MOOSE M-R2
  (`gold/moose_mode2_mac_np2_farms46ad2eb/`) the curve agrees to 6e-12
  relative on every step, peak 423.02086 at t = 0.0097 in both, and the
  fixed-point pass count is the same on all 1000 steps, including the 548
  steps (453–1000) that both codes accept at the 20-pass cap.
  `tests/verification/raccoon_mode2_moose.sh` (gate 1e-9).
