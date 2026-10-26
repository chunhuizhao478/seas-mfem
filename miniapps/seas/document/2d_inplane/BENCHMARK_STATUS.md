# 2D in-plane benchmarks: status

Branch `2d-inplane` (design report rev 3, `/mnt/project-files/2d-inplane/2d-inplane-design-rev3.pdf`,
section 6).  One row per gate; a benchmark is completed when Chunhui Zhao says so and its gold is in
the repo.

| Gate | Benchmark | Status | Gold | Reference data |
|---|---|---|---|---|
| G1 | TPV205-2D (slip weakening, plane strain) | **completed** 2026-10-26 | `tpv205/gold/nitsche-2d-kv0p1_p2_cd_pc_100m_cloud_gitb0ce3cb` (p2, 100 m, b0ce3cb6) | `tpv205/benchmark_data/PYLITH2D_200m` |
| G2 | TPV105-2D (rate and state, thermal pressurization) | **completed** 2026-10-26 | `tpv105/gold/nitsche-2d-kv0p1_p2_cd_pc_tp_100m_cloud_gitb0ce3cb` (p2, 100 m, b0ce3cb6) | `tpv105/benchmark_data/MDSBI2D_50m`, `tpv105/benchmark_data/SPECFEM2D_50m` |
| G3 | BP3-QD (dipping fault, quasi-dynamic cycles) | phase E: the 2D quasi-dynamic solve and the decks, inputs and outputs in the branch (`config/fd_cycle/bp3_qd_*.toml`, `cycle/scec_bp3_export.py`); the 25 m runs not started | - | SEAS platform BP3-QD (on the Mac, not in the repo) |
| G4 | BP3-FD (fully dynamic cycles) | phase F, after G3 | - | - |

Chunhui Zhao, 2026-10-26: "I think we can gold tpv205-2d, tpv105-2d results p=2, and mark as
completed".

Summary of the two completed gates (details and every number in the gold READMEs):

- TPV205-2D, p2 100 m vs PyLith 200 m: final slip within 1.4 % from -12 to +4.5 km and +3.1 / +3.3 % at
  +7.5 / +12 km; arrivals within one PyLith sample (0.06 s) when read at its sample times; peak slip
  rates above PyLith's (not converged at this resolution, ours rise with refinement).  Open: +7.5 km
  (slip +3 %, peak +22 % or more at every resolution compared); PyLith is the only reference code.
- TPV105-2D, p2 100 m vs MDSBI 50 m, with SPECFEM2D 50 m as the second reference: from -15 to +15 km
  arrivals within 0.015 s and final slip within 2.2 %; peak slip rate within -3.6 / +5.3 %, peak
  temperature within 2 %, pore pressure within 1 MPa; final slip between the two references from -18 to
  +9 km and below both at +15 and +18 km (-2.2 and -6.1 % against MDSBI).
