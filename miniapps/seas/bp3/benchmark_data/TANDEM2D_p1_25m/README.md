# Tandem 2-D, p = 1, BP3-QD dip 60 thrust on the 25 m mesh (reference for the seas-mfem p1 gold)

Run by Claude for Chunhui Zhao on his Mac (2-D in-plane design, gate G3, decision D4 of 2026-10-26: "Tandem 2D
p=1 and p=2 on the Mac, compared with seas-mfem at matching orders"). It is not a SCEC submission: the SCEC
platform's Tandem entry for BP3 is p = 8 on 250 m elements.

- Code: Tandem (github.com/TEAR-ERC/tandem) a26753d645af96af68a869d402a8e751a4aba3e0 (v1.0-201-ga26753d); app/
  and src/ unmodified (one job script in the checkout was modified and not built). Built by run/build_tandem.sh
  (run/build_tandem.out): CMake Release, DOMAIN_DIMENSION 2, POLYNOMIAL_DEGREE 1, ARCH noarch, with PETSc 3.21.6
  (MUMPS, ParMETIS), clang 22 and Open MPI 5.0.10 from the conda env tandem-petsc321.
- Input: run/bp3.toml (sha256 7423dcffa299570ee8dc185e0d2ad42d32733f95965575ef84885e387448915a) with Tandem's
  example run/bp3.lua (scenario bp3_d60_reverse; sha256 49a50b86...f162), mode QD, ref_normal [1, 0], final time
  4.7304e10 s, on miniapps/seas/bp3/mesh/bp3_d60_25m.msh (sha256 fd08af92...6d3c, the mesh of the seas-mfem gold).
  The toml adds the 12 SCEC fault probes (same points as the seas-mfem stations) and 8 free-surface domain probes;
  fn-00 and fn+00 sit at x1 = -1 m and +1 m.
- Run: `mpirun -np 4 tandem bp3.toml --mode QDGreen --petsc -options_file lu_mumps.cfg -options_file rk45.cfg`
  (run/run_info.txt, run/restart_qdgreen.sh; RK 5(4) Dormand-Prince, atol 1e-7, rtol 1e-50, infinity norm; LU by
  MUMPS), OMP_NUM_THREADS=1, 2026-10-26 22:51:00 to 2026-10-27 02:32:47 PDT, wall 13,307 s; 742,981 steps and
  34,165 rejected; 16 events (first at 176.7347 yr). Two earlier starts were stopped and are not used: one with a
  7 MiB stack, and one in the default QD mode (about 0.5 s per step).
- Output sampling: Tandem's adaptive rule (app/tandem/AdaptiveOutputStrategy.h): a row when the fault-wide max
  slip rate has changed by 1 % since the last row (rtol 0.01, atol 1e-50, the defaults) or after t_max = 9460800 s
  (0.3 yr, set in bp3.toml). Every probe file has the same 67,699 rows, t = 0 to 4.730110667668105e10 s
  (1499.92 yr). The rows are dense (ms) where max V changes fast and leave gaps of up to about 2 s on the
  mid-rupture plateau.
- Files: fltst_dp*.csv ("Time", "state", "slip0", "traction0", "slip-rate0", "normal-stress"; SI with stresses in
  MPa, traction0 negative for this thrust) and srfst_fn*.csv ("Time", "u0", "u1"; u1 up). probes.sha256 and
  run.sha256 hold their sha256; each file was moved off the Mac gzipped and checked against the sha256 of the Mac's
  copy.
- Compared with seas-mfem in bp3/gold/qd-2d_p1_d60_thrust_25m_mac_git9503fe2 (README there).
