# TPV102 gold: Nitsche dynamic rupture, Kelvin-Voigt damping, p1, 200 m (Frontera job 7995743)

- Code: seas-mfem branch nitsche-dynamic-rupture, commit 2e0c04a (built on Frontera as b0c0228 + kv_node.patch,
  byte-identical to 2e0c04a; verified, see damping_tests/github_vs_frontera). Driver seas_nitsche_dyn_driver.
- Deck and job script: jobs/ (= configs/ at 972f139): Prakash-Clifton v* 1 m/s, L 0.01 m; [numerics] kv_eta 0.1,
  kv_scale "node", kv_length "inscribed" (beta_e = 0.1 (1/3) (6V/A)/(p c_s)); cfl 0.5.
- Mesh: tpv102_200m.msh (2,466,377 tets, sha256 eea25a27813dafd60a5cb5bc3bf1c8de4674e9b07858c3f9b7b8d159fa5be219; not in the repo).
- Run: 8 nodes / 448 ranks (sbatch -N 8 -n 448; script says 4/224), 2026-10-15 14, 0-12 s, 56-16 steps, dt 09 PDT s, central difference ("cd"), --station-copies first.
- results/: 9 SCEC station files (t slip1 slip2 V1 V2 tau1 tau2 sigma_n log10_theta; 1 = dip, 2 = strike; SI) + driver.log.
  Shear includes the damping traction (ClosureGamma output fix).
- plots/: per-station PNGs run vs DRDG3D / PyLith, overview of h-slip-rate; comparison PDFs *_undamped_damped_vs_drdg3d.pdf (undamped Nitsche PC reference 7994409 / 7994410 vs this run vs DRDG3D / PyLith, with a summary table) and *_damped_vs_benchmark.pdf (this run vs DRDG3D / PyLith).
- vs DRDG3D (median of 9 stations): arrival 47530, final slip 0.000252475, peak h-slip-rate -0.141 s;
  noise RMS sigma_n / tau_h / tau_v -0.12 %:-31.8 %:0.051 / 0.040 / 0.047 (DRDG3D 0.016 / 0.035 / 0.050) MPa.
- Requested by Chunhui Zhao (damping-tests thread, 2026-10-16).
