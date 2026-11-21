# M-R1: MOOSE run of RACCOON mode1_brittle_fracture (full curve)

- Run: 2026-11-21, Andy's Mac (arm64), `mpiexec -n 4`, 124.5 s wall
  (MOOSE: 112.54 s, 156 MB). Inputs byte-copied from
  farms_cdms `raccoon/tutorials/mode1_brittle_fracture/{elasticity,fracture}.i`
  (raccoon 9fd3a6f) into `~/pff_moose_refs/M-R1_mode1/`.
- Executable: farms_cdms `farms-opt` (mtime 2026-06-05 22:55:49), farms_cdms
  HEAD 46ad2eb0dab1bc27f22908421c857ab0ddd9c1c6 with a dirty tree (364 files
  changed, almost all deleted cdm_examples outputs), built against
  ~/projects/moose 94b6528fe26f481eb2d84c3a729a3d3c05f01ec7 (2026-04-01).
- Command (in the run folder):
  `PATH=$CONDA/envs/moose/bin:$PATH SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX15.4.sdk mpiexec -n 4 farms-opt -i elasticity.i Outputs/csv=true`.
  SDKROOT is needed only for the ADFParser JIT link (the conda ld64 cannot
  link against the default MacOSX26.5 SDK); it does not change the numerics.
- `elasticity_out.csv`: time, Fy (175 steps, dt = 2e-5).
- Fixed-point passes ("Beginning fixed point iteration" lines): steps 1–3
  take 1, 1, 2; steps 132–156 hit the cap of 20 (133–156 report
  REACH_MAX_ITS and are accepted); steps 157–160 take 8, 5, 5, 5; the other
  151 steps report CONVERGED_ABS. Peak Fy = 691.00227 at t = 2.86e-3 (step 143).
