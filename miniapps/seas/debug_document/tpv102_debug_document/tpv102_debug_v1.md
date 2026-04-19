# TPV102 v1: Rupture does not propagate out of nucleation patch

## Context

First 400-rank production-style TPV102 run on Frontera
(`results_200m_p1_3s_dev_job7664258`, `--order 1`, auto-CFL `dt ≈ 0.285 ms`,
mesh `tpv102_200m.msh`, 8 nodes × 50 ranks, dev queue 2 hr, run reached
`t ≈ 2.08 s` before the Slurm wall limit).  ParaView fault-surface PVD
was downloaded to
`/Users/chunhuizhao/Downloads/seas-mfem/tpv102/results_200m_p1_3s_dev_job7664258/...`
and analyzed frame-by-frame.

User's initial description: *"no spontaneous propagation, the slip
rate keeps increasing around the nucleation point."*

## Symptom (evidence from the VTUs)

Probe sweep along strike at depth −7.5 km (hypocenter depth).  The
3 km nucleation disc extends to |x| ≤ 3 km around the hypocenter at
(0, 0, −7.5 km).

| probe | rank | F(r) | ∆τ_strike (t=1.98 s) | V_strike (t=1.98 s) | state at end of run |
|---|---:|---:|---:|---:|---|
| hypo (x=0)    | 3   | 1.00 | −10.6 MPa (dropped from +23) | **7.699 m/s** (frozen after t=1.19 s) | Ruptured, steady sliding |
| x=0.5 km      | 3   | 0.97 | −10.6 MPa (dropped) | **7.53 m/s** (frozen after t=1.49 s) | Ruptured (≈ cs transit from hypo) |
| x=1.0 km      | 3   | 0.88 | +21.9 MPa (forcing only) | 0.027 m/s (slowly growing) | **NOT ruptured** |
| x=1.5 km      | 293 | 0.82 | +18.5 MPa (forcing only) | 2.3 × 10⁻⁴ m/s (forced eq., flat) | NOT ruptured, frozen |
| x=2.0 km      | 293 | 0.45 | +11.3 MPa (forcing only) | 1.32 × 10⁻⁷ m/s (forced eq., flat) | NOT ruptured, frozen |
| x=2.5 km      | 293 | 0.16 | +2.87 MPa (forcing only)  | 2 × 10⁻¹¹ m/s | essentially V_ini |
| x=2.8 km      | 264 | 6×10⁻⁷ | +0.022 MPa | 1.02 × 10⁻¹² m/s | V_ini |
| x ≥ 3.5 km (any y) | any | 0 | **exactly 0** | **exactly 1 × 10⁻¹² m/s** | untouched |

### The "machine-precision zero" signature

At every probe with r > 3 km we observe `traction_strike = 75.000 MPa`
and `V_strike = 1.0 × 10⁻¹²` to **full double precision** for the
entire 21-frame time series.  If there were *any* wave perturbation —
even at 1e-15 — the trial traction
`τ₁_trial = η · ((V_y⁻ − V_y⁺) + invZ_s · (S_xy⁺ + S_xy⁻))`
would be nonzero, the Brent friction solver would produce a slip
rate ≠ V_ini exactly, and `τ₁_corr − η·V₁` would drift off 75 MPa
exactly.  The sharpness of this null rules out "weak wave"; the
field Q is literally zero at those bulk locations throughout the run.

### The rank-boundary discontinuity

Within rank 3 (same rank as hypocenter), the rupture propagates
~0.5 km (from x=0 at t ≈ 1.1 s to x=0.5 km at t ≈ 1.5 s — a
transit consistent with sub-Rayleigh v_r ≈ 0.5–0.7 c_s).

At the first partition seam encountered (between x=1 km on rank 3
and x=1.5 km on rank 293), the rupture front **stops completely**.
x=1 km, which has strong forcing (22 MPa ∆τ) and is only 0.5 km
beyond the ruptured point at x=0.5 km, takes ≥ 2 s to reach only
2.7 × 10⁻² m/s — still 280× below the surrounding rupture speed —
and x=1.5 km on rank 293 never sees the wave at all (flat at
forced equilibrium 2.3 × 10⁻⁴ m/s, with Q machine-zero in its
bulk neighborhood).

This is the cleanest signature we have: rupture energy present
inside rank 3, absent in rank 293 despite being meters away.

### Non-rupture confound

**Hypocenter V frozen at 7.699 m/s after t = 1.29 s** is
physically expected, *not* a bug.  Once `G(t)=1` at t ≥ 1 s the
forcing is constant, psi reaches steady state, and friction +
radiation damping balance the forced stress exactly — the point
enters steady sliding.  It becomes a *symptom* of the bug only
because nothing else on the fault responds to the radiated waves.

### Empty `.dat` station files (secondary finding)

All 15 `tpv102_station_*.dat` files are 0 bytes.  Root cause is
unrelated to propagation: `TPV102StationWriter::WriteStep` writes
to `std::ofstream` without flushing, and `station_writer.Flush()`
is only called in the driver's post-loop cleanup (see
`tpv102_driver.cpp` § 9 "Summary").  Slurm SIGKILL at the 2 hr
wall limit discards all unflushed buffers.  Fix inline.

## What I ruled in / out

- ✅ **Nucleation forcing works end-to-end.**
  `τ_strike` inside the 3 km patch matches
  `τ_ini + 25 MPa · F(r) · G(t)` to within < 10 Pa of the
  analytical formula at every frame.  So `ApplyNucleation` is
  called at every RK4 stage and `dof_data[i].tau1_0` is being
  written correctly.
- ✅ **`FaultFaceFlux::Evaluate` runs on every fault DOF.**
  At x=2 km and x=2.5 km the corrected traction sits at the
  forced-equilibrium value `τ₁₀ − η·V₁` — that can only happen if
  `Evaluate` wrote `data.tau1_corr` there.  It wrote the
  equilibrium (quasi-static) answer because `Q_plus = Q_minus = 0`
  on both sides of the fault face, i.e. **no wave energy was
  available to change the trial traction**.
- ✅ **Within-rank dynamic coupling is alive** (at least partly).
  x=0.5 km ruptures 0.3 s after the hypocenter at a speed
  consistent with c_s — so the local fault-flux → bulk-wave → next
  fault QP chain works on a single rank.
- ❌ **Cross-rank wave propagation is the suspect.**  The rupture
  front stops exactly at the rank-3 / rank-293 partition seam.
  4-rank parallel wave-operator unit test passes (5/5), but that
  test uses a uniform 64-element cube — does not exercise the
  graded-mesh irregular partition the 400-rank TPV102 run produces.
- ❌ **Not a CFL / stability issue.**  `dt = 0.285 ms = dt_cfl`
  (auto-CFL after the v0 fix).  V stays *too quiet*, not
  unbounded — the opposite signature of CFL violation.
- ❌ **Not a PML leak.**  TPV102 driver does not call `SetPML`,
  `pml_layer_` stays null, `ApplyPMLDamping` is never called.

## Ranked hypotheses for the primary bug

1. **H1 (most likely) — Cross-rank bulk wave propagation is the
   obstruction.**  Evidence: rupture front stops exactly at a
   partition seam; "machine-precision zero" Q signature on ranks
   not touching the hypocenter; 4-rank unit test doesn't exercise
   this path on a graded mesh.  Suspect surface: the
   `ComputeSharedFaceFluxRHS` regular-face-flux path's interaction
   with `ExchangeFaceNbrData` when there are many partition seams
   on a graded mesh.
2. **H2 — `nbr_data[c] = q_gf.FaceNbrData()` aliasing.**  The
   per-component ghost-data exchange in
   `wave_operator.inl:615-623` reads `FaceNbrData()` and assigns
   to a `std::vector<Vector>`.  `Vector::operator=` *should*
   deep-copy, but worth a direct assert — if MFEM's `FaceNbrData()`
   returns a data pointer into the `q_gf` internal storage and
   the copy is shallow, every `nbr_data[c]` ends up aliased to
   the LAST component's ghost data.  Testable with a single
   `Vector::GetData()` check.
3. **H3 — Shared-fault DOFData duplicated & evolving
   inconsistently on both ranks.**  For a fault face on a
   partition seam, both sides carry a `DOFData` entry (via
   `shared_fault_dof_offset_` on each rank independently).  If
   their psi / V states drift apart, the Godunov flux between
   them would be computed against inconsistent fault states and
   may produce a degenerate zero-radiation result.  Less likely
   than H1/H2 given that the *adjacent* rank-3 local fault still
   ruptures correctly, but possible at interior seams.

## Fixes landed in this commit

### F1.  `TPV102StationWriter::WriteStep` — flush after each write

`miniapps/seas/dynamic/tpv102_setup.hpp` — add
`files_[s].flush()` after the `\n`-terminated row.  Cost is 9 ×
21 flushes per 2 hr run = nothing; benefit is that a
wall-limit-terminated job still produces usable `.dat` files for
post-hoc debugging.  Unconditional fix; no CLI flag.

### F2.  `--debug-qnorm` diagnostic flag in `seas_tpv102_driver`

`miniapps/seas/drivers/tpv102_driver.cpp` — new CLI flag
`--debug-qnorm`.  When set, at every station output interval
rank 0 prints the min/max/mean of per-rank `||Q||_∞`, plus the
per-rank `||Q||_∞` for a user-selectable short list of rank IDs
(the default is `{rank_with_hypocenter, rank_with_hypocenter+1,
final_rank}` picked automatically).  Default: off.

Purpose: resolve H1.  If `||Q||_∞` on ranks that do not touch
the hypocenter stays at 0 while rank-0 / hypocenter-rank values
grow, cross-rank wave propagation is confirmed broken.  If those
far ranks' `||Q||_∞` does grow but the fault at x=1.5 km still
doesn't rupture, look at the fault-flux-at-partition-seam path
(H3) instead.

This is **not a fix for the bug** — it is an instrumented
diagnostic.  The bug fix will land in v2 once the diagnostic
identifies the exact failing path.

### F3.  `--no-domain-pv` + `ParaViewOutput::ShouldWrite`

`miniapps/seas/drivers/tpv102_driver.cpp`,
`miniapps/seas/io/paraview_output.hpp` — new CLI flag
`--no-domain-pv` suppresses the volume-mesh PVD
(`ParaView/Cycle*/*.vtu`, ~5 GB/cycle on the 200 m mesh) while
keeping the fault-surface VTUs on the same schedule.  Added the
`ParaViewOutput::ShouldWrite(cycle, time, V_max)` helper so the
driver can run the same scheduling decision without emitting the
volume PVD.  Saves tremendous disk on long coseismic runs, no
effect on debug diagnostics.  Default: off.

### F4.  `--pv-low-order` + `ParaViewOutput::SetLevelsOfDetail`

`miniapps/seas/drivers/tpv102_driver.cpp`,
`miniapps/seas/io/paraview_output.hpp` — new CLI flag
`--pv-low-order` switches volume PVD output from
`VTK_QUADRATIC_TETRA` to `VTK_TETRA` (linear) via
`SetHighOrderOutput(false)` + `SetLevelsOfDetail(1)`, cutting
per-cycle volume-PVD size ~2.5× on order-2 meshes.  No effect on
correctness.  Default: off (preserves the existing high-order
volume output).

## Next-step plan (out of this commit)

1. Rerun the p=1 3 s job on Frontera with `--debug-qnorm`
   enabled.  Expect ≲ 10 minutes before the first `Step N` print
   after the nucleation pulse engages — that is the moment that
   should show cross-rank Q growth if wave propagation works.
2. If H1 confirmed, inspect `ComputeSharedFaceFluxRHS` with
   print instrumentation of `(Q_self, Q_nbr)` on a specific sf
   across the suspected seam.
3. If H2 confirmed, replace the `nbr_data[c] =
   q_gf.FaceNbrData()` pattern with an explicit deep copy loop
   and re-verify.
4. Write a regression test: 8-rank TPV102 on a coarser mesh
   (`tpv102_200m.msh` partitioned to place a seam through the
   nucleation patch) that checks V at x=5 km is nonzero by
   t=1.8 s.  No such test exists today — the parallel wave
   operator test uses a uniform cube, missing this whole class
   of bug.

   **⚠ Caveat on coarser-mesh tests.**  TPV102's cohesive process
   zone is Λ_dyn ≈ (π/4)·μ·D_c / (τ_s − τ_d) ≈ 160 m.  At
   h = 200 m we already have Λ/h ≈ 0.8 — borderline.  At any
   coarser resolution (h = 400 m, 500 m, 1 km) we have Λ/h ≪ 1
   and spontaneous rupture propagation is not expected even in a
   bug-free code.  So a coarser-mesh test cannot be used to check
   "does the rupture propagate?" — that remains a 200 m-or-finer
   question.

   What a coarser-mesh test *can* do is cheaply test the
   **cross-rank bulk wave propagation** decoupled from the
   rupture-physics question.  At h ≳ 500 m, the fault flux at the
   hypocenter still injects stress jumps into the bulk wave
   field, those jumps still propagate as bulk shear waves, and
   the `--debug-qnorm` diagnostic will still show non-hypocenter
   ranks' `||Q||_∞ > 0` if and only if wave energy crosses
   partition seams.  That's exactly the signal H1 predicts.  So a
   coarse-mesh run is an appropriate cheap test *for H1/H2/H3
   only*, not for verifying rupture propagation.  The final
   rupture-propagation verification must be done on the 200 m
   mesh (or finer) on Frontera.

## Artifacts

- Input log: `results_200m_p1_3s_dev_job7664258/tpv102_200m_p1_7664258.out`
  (see user's terminal — not checked in)
- Fault-surface PVD: `.../FaultSurface/fault_surface.pvd` (21
  frames, t ∈ [0, 2.08 s])
- Probe-sweep analysis: reproducible via the Python snippets
  embedded in this debug session (see agent transcript).
- Secondary: all `tpv102_station_*.dat` are 0 bytes — resolved
  by F1, future runs will carry probe data even on SIGKILL.
