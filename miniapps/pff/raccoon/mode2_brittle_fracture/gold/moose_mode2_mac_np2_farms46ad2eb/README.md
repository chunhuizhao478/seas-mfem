# M-R2: MOOSE run of RACCOON mode2_brittle_fracture (full curve)

- Run: 2026-11-21, Andy's Mac (arm64), `mpiexec -n 2`, 5284 s wall
  (MOOSE: 5276 s, 167 MB). Inputs byte-copied (sha256 checked) from
  farms_cdms `raccoon/tutorials/mode2_brittle_fracture/{elasticity,fracture}.i`
  (raccoon 9fd3a6f) into `~/pff_moose_refs/M-R2_mode2/`; the inputs build
  their own mesh.
- Executable: farms_cdms `farms-opt`, farms 46ad2eb, moose 94b6528.
- Command (in the run folder):
  `PATH=$CONDA/envs/moose/bin:$PATH SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX15.4.sdk mpiexec -n 2 farms-opt -i elasticity.i Outputs/csv=true`.
- Mesh: 900 elements / 976 nodes before, 2826 / 3026 after the 2 initial
  adaptivity steps (the R2 mesh gate).
- `elasticity_out.csv`: time, Fx (1000 steps, dt = 2e-5, to t = 0.02).
  Steps 1-2: Fx = 0.92368182656251, 1.8473633560272 (the gold values).
  Peak Fx = 423.02086 at t = 0.0097 (step 485).
- Fixed-point passes per step (cap 20, `accept_on_max_fixed_point_iteration
  = true`): 1-11:1 12-88:2 89-172:3 173-237:4 238-285:5 286-321:6 322-348:7
  349-368:8 369-384:9 385-397:10 398-408:11 409-416:12 417-423:13 424-429:14
  430-435:15 436-439:16 440-443:17 444-446:18 447-449:19 450-1000:20.
  Steps 450-452 converge (CONVERGED_ABS) on the 20th pass; steps 453-1000
  (548 steps, including the peak and the whole post-peak branch) end at the
  cap with REACH_MAX_ITS and are accepted unconverged. A pff comparison past
  step 452 therefore compares capped iterations; it is exact only if pff
  reproduces each pass.
