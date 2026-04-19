# Code Review: tpv102_debug_v1 and related changes

## Review Scope
- Plan/debug doc: `miniapps/seas/debug_document/tpv102_debug_document/tpv102_debug_v1.md`
- Files reviewed:
  - `miniapps/seas/drivers/tpv102_driver.cpp` (recently modified — `--debug-qnorm`, `--no-domain-pv`, `--dt`, `--pv-low-order`, PV integration)
  - `miniapps/seas/dynamic/tpv102_setup.hpp` (F1 flush fix in `TPV102StationWriter::WriteStep`, surrounding `TPV102SurfaceStationWriter`)
  - `miniapps/seas/io/paraview_output.hpp` (orphan-face `-1` guards, `SetLevelsOfDetail`, `ShouldWrite`)
  - `miniapps/seas/dynamic/wave_operator.inl` (shared-face fault path — H1/H2/H3 suspects)
  - `miniapps/seas/dynamic/fault_face_flux.cpp` (trial-traction sign convention)
- Domain context: `miniapps/seas/CLAUDE.md` (numerical-invariant list), SCEC TPV101/102 spec

## Findings

### [R-001] [CRITICAL] `dynamic/wave_operator.inl:ComputeSharedFaceFluxRHS` + `dynamic/fault_face_flux.cpp:Evaluate` — Shared fault DOFData evolves inconsistently on the two ranks sharing a fault face (debug-doc H1/H3 confirmed by code inspection)

**Category:** BUG

**Description:**
A fault face that crosses a rank partition is represented as a shared face by MFEM. `ParMesh::GetSharedFaceTransformations(sf)` returns `Elem1No = local element, Elem2No = ghost` on **every** rank. So both ranks that share the face call `fault_flux_->Evaluate(fdata_local, Q_plus=Q_self, Q_minus=Q_nbr, …)` with:
- Rank A (A-side local): `Q_plus = Q_A`, `Q_minus = Q_B`
- Rank B (B-side local): `Q_plus = Q_B`, `Q_minus = Q_A`

`ComputeTrialTraction` (Eq. 7b) is NOT invariant under `(+,−)` swap:

```
τ₁_trial(Q+,Q−) = η_s·(V−,t1 − V+,t1) + η_s·(τ+/Z+ + τ−/Z−)
```

The velocity-jump term flips sign, the stress term does not. So `tau1_trial_A ≠ tau1_trial_B` in general (they differ by `2η_s·(V_A,t1 − V_B,t1)`). Cascading through `Evaluate`:

- `tau1_total = tau1_0 + tau1_trial` → different on the two ranks
- `V1 = V_abs · tau1_total / (strength + η_s·V_abs)` → different magnitude AND potentially opposite sign
- `data.tau1_corr`, `data.V1`, `data.slip1` → drift apart step by step

Each rank stores an **independent** `DOFData` entry for the same physical fault QP (`shared_fault_dof_offset_` is populated independently on each rank in `wave_operator.hpp:138-146`). Each rank updates its own copy with its own `+` / `−` interpretation. No reconciliation step exists.

This directly matches hypothesis **H3** in `tpv102_debug_v1.md` and explains the **H1** symptom ("rupture front stops exactly at the rank-3 / rank-293 seam"): at a partition seam, the two halves of the fault QP compute different corrected tractions, feed different `Q_imp` into `flux_.Interior`, and the net momentum/energy transferred across the seam is not the true physical value — so the bulk shear wave that *should* advance the rupture into rank 293 is either attenuated or cancelled.

It also explains the station-output-empty secondary finding only **partially** — the .dat files are empty primarily because of the flush-on-SIGKILL issue (F1), but the data that *would* have been written at fault QPs straddling seams would also be inconsistent.

**Trigger:**
Any run with a ParMesh partition that places a seam through the fault surface (i.e., every realistic production run). 4-rank uniform-cube unit tests miss this because their graded-mesh partition is too regular to cross the fault.

**Actual behavior:**
Per-rank `DOFData.tau1_corr`, `DOFData.V1`, `DOFData.slip1`, `DOFData.psi` for a shared fault QP diverge across ranks as the simulation advances. `F_h` accumulated into the volume RHS on the two sides is not Godunov-consistent, so wave energy does not propagate correctly across the partition seam.

**Expected behavior:**
The two ranks must agree on the corrected fault state and the imposed `Q_imp` side states. Per-QP state must be single-sourced.

**Suggested fix (architectural — several acceptable variants; the simplest):**
Give a single rank canonical ownership of each shared fault QP. Adopt the "lower rank ID owns" convention (easily available from MPI). Only the owning rank calls `Evaluate`; the non-owner receives `{V1, V2, tau1_corr, tau2_corr, sigma_n_corr, psi}` via an MPI exchange and uses those to construct `Q_imp_plus/minus` locally. Concretely, add an exchange step after all local+shared fault-flux computations in `WaveOperator::Mult`:

```diff
 ComputeFaceFluxRHS(Q, dQdt);

 if constexpr (IsParallelMesh<MeshType>::value)
 {
+   // Canonical owner (lower rank ID) computes; non-owner receives.
+   // Skip the Evaluate call on non-owner shared fault faces and instead
+   // copy {psi,V1,V2,tau*_corr,sigma_n_corr} from owner's DOFData over MPI.
    ComputeSharedFaceFluxRHS(Q, dQdt);
+   SyncSharedFaultDOFData(*fault_dof_data_);
 }
```

Implementation detail: add a `shared_fault_owner_rank_` array (size `GetNSharedFaces()`) that stores the peer rank for each shared fault face; non-owners `continue` past the friction solve and `Evaluate` inside `ComputeSharedFaceFluxRHS`, then post-pass synchronises the 6 scalars with `MPI_Sendrecv`.

Below-the-fix alternative (simpler and possibly enough): in `ComputeSharedFaceFluxRHS`, on the non-owner side, call `Evaluate` with `(Q_plus=Q_nbr, Q_minus=Q_self)` so both ranks use the same global `(+,−)` convention. Requires agreeing on which rank is the "+" side (again, lowest rank ID). Then the DOFData updates match bit-for-bit. This is a single swap at `line 734`:

```diff
-fault_flux_->Evaluate(fdata, Q_plus_local, Q_minus_local,
-                      Q_imp_plus, Q_imp_minus);
+const bool i_am_plus_owner = my_rank < peer_rank;
+if (i_am_plus_owner) {
+   fault_flux_->Evaluate(fdata, Q_plus_local, Q_minus_local,
+                         Q_imp_plus, Q_imp_minus);
+} else {
+   // Peer has the canonical "+" side. Swap so we use the same convention.
+   fault_flux_->Evaluate(fdata, Q_minus_local, Q_plus_local,
+                         Q_imp_minus, Q_imp_plus);
+}
```

Either path must be accompanied by a documented invariant: *"For shared fault QPs, `DOFData.tau1_corr / V1 / slip1` are the same on both ranks."*

**Test case:**
```python
def test_R001_shared_fault_dofdata_consistency():
    """
    Run 8-rank TPV102 on coarse mesh (h=500m is fine since we are not
    checking rupture physics; see v1 doc caveat).  After one RK4 step,
    MPI_Gather every shared fault QP's DOFData.{V1, tau1_corr} from both
    ranks that own it; assert they are bitwise equal.
    """
    run_driver(mesh="tpv102_500m.msh", nranks=8, dt=5e-4, nsteps=1,
               extra_flags=["--dump-shared-fault-dofdata"])
    owner, peer = load_shared_fault_dofdata_dump()
    # For every shared fault QP the two ranks' DOFData must agree
    assert_array_equal(owner["V1"],          peer["V1"])
    assert_array_equal(owner["tau1_corr"],   peer["tau1_corr"])
    assert_array_equal(owner["sigma_n_corr"], peer["sigma_n_corr"])
    assert_array_equal(owner["psi"],         peer["psi"])
```

---

### [R-002] [MODERATE] `drivers/tpv102_driver.cpp:main` — `--debug-qnorm` output deviates silently from what the debug doc advertises

**Category:** DEVIATION

**Description:**
`tpv102_debug_v1.md` (F2) promises two things from the `--debug-qnorm` flag:
1. min/max/mean of per-rank `||Q||_∞`
2. per-rank `||Q||_∞` for a user-selectable rank list — default `{rank_with_hypocenter, rank_with_hypocenter+1, final_rank}` picked automatically

The implementation (`tpv102_driver.cpp:811-815`) only prints item (1) plus a count of ranks with `||Q||_∞ == 0`. It does NOT print per-rank values for the hypocenter neighborhood and does NOT auto-detect the hypocenter rank. The promised diagnostic is the *primary* test for H1 (is the hypocenter neighbor rank's Q growing while far ranks' Q stays zero?), and the aggregate min/max cannot answer that — the per-rank breakout is the entire point.

**Trigger:**
Run with `--debug-qnorm` enabled. Examine the `[qnorm]` output lines.

**Actual behavior:**
```
  [qnorm] min=0 max=1.2e-6 mean=3e-9 #ranks_with_||Q||=0: 395/400
```

Cannot tell which 5 ranks have nonzero Q — in particular, cannot distinguish "hypocenter + 4 neighbors ruptured, nothing else" from "5 random ranks".

**Expected behavior:**
Also print per-rank `||Q||_∞` for a short list of ranks, defaulting to
hypocenter-owning rank + a couple of neighbors + the highest rank (the "farthest-away" ones the diagnostic cares about).

**Suggested fix:**
Add a per-rank breakout that identifies, at setup time, the rank that owns the hypocenter QP, then prints `{hypo, hypo+1, hypo+Δ, last}`:

```diff
 if (rank == 0)
 {
    real_t qmin = qn_all[0], qmax = qn_all[0], qsum = 0.0;
    int nzero = 0;
    for (int r = 0; r < nprocs; r++)
    {
       qmin = std::min(qmin, qn_all[r]);
       qmax = std::max(qmax, qn_all[r]);
       qsum += qn_all[r];
       if (qn_all[r] == 0.0) { nzero++; }
    }
    std::cout << "  [qnorm] min=" << qmin
              << " max=" << qmax
              << " mean=" << (qsum / nprocs)
              << " #ranks_with_||Q||=0: " << nzero
              << "/" << nprocs << "\n";
+   // Per-rank breakout for H1 diagnosis: hypocenter owner (if detected) +
+   // neighbors + final rank.  hypo_rank is computed once up front via
+   // MPI_Allreduce on which rank contains the hypocenter fault QP.
+   std::vector<int> watch = {hypo_rank,
+                             std::min(hypo_rank + 1, nprocs - 1),
+                             std::min(hypo_rank + 4, nprocs - 1),
+                             nprocs - 1};
+   std::cout << "  [qnorm:watch]";
+   for (int r : watch) {
+      std::cout << " r" << r << "=" << qn_all[r];
+   }
+   std::cout << "\n";
 }
```

To detect `hypo_rank` at init, do one MPI_Allreduce(MAX) of `{distance_to_origin == 0 ? my_rank : INT_MIN}` after fault-DOF setup.

**Test case:**
```python
def test_R002_debug_qnorm_prints_hypocenter_rank():
    out = run_driver(extra=["--debug-qnorm", "--tfinal", "0.1",
                            "--dt", "0.01"])
    assert "[qnorm] min=" in out
    assert "[qnorm:watch]" in out        # per-rank breakout line
    # must include the hypocenter-owning rank explicitly:
    assert re.search(r"\[qnorm:watch\].*\br\d+=.*\br\d+=", out)
```

---

### [R-003] [MODERATE] `dynamic/tpv102_setup.hpp:TPV102SurfaceStationWriter::WriteStep` — Missed F1 fix; surface-station files still lost on SIGKILL

**Category:** BUG (missed fix)

**Description:**
F1 in `tpv102_debug_v1.md` added `files_[s].flush();` inside `TPV102StationWriter::WriteStep` (line 299) so that Slurm SIGKILL at the wall limit doesn't discard buffered probe data. The identical pattern in `TPV102SurfaceStationWriter::WriteStep` (lines 441-442) was **not** updated, even though it has exactly the same failure mode. The 6 SCEC surface stations (velocity time series) will be 0 B on SIGKILL — this is a silent data-loss hazard for the very next production run.

**Trigger:**
Any wall-limited Frontera run (i.e., every run). The bug does not surface locally because local runs exit cleanly.

**Actual behavior:**
```cpp
files_[s] << std::scientific << std::setprecision(10)
          << t << " " << vx << " " << vy << " " << vz << "\n";
// NO FLUSH
```

**Expected behavior:**
Same pattern as the F1 fix.

**Suggested fix:**
```diff
 files_[s] << std::scientific << std::setprecision(10)
           << t << " " << vx << " " << vy << " " << vz << "\n";
+files_[s].flush();  // F1 parity: keep probe data on Slurm SIGKILL
```
At line 442 of `tpv102_setup.hpp`.

**Test case:**
```python
def test_R003_surface_station_survives_sigkill():
    """Launch driver, SIGKILL it mid-run, verify surface station .dat files
       contain non-zero data."""
    proc = start_driver(flags=["--tfinal", "10.0", "--dt", "0.01"])
    time.sleep(30)  # let several output cycles complete
    proc.send_signal(signal.SIGKILL)
    proc.wait()
    for name in ["surf_0_9", "surf_0_n9", "surf_12_6",
                 "surf_12_n6", "surf_n12_6", "surf_n12_n6"]:
        p = Path(f"out/tpv_station_{name}.dat")
        assert p.stat().st_size > 0, f"{p} is empty — flush missing"
```

---

### [R-004] [MODERATE] `drivers/tpv102_driver.cpp:main` — `--no-domain-pv` is undocumented scope creep relative to the v1 commit

**Category:** DEVIATION

**Description:**
`tpv102_debug_v1.md` § "Fixes landed in this commit" lists exactly two changes: F1 (flush) and F2 (`--debug-qnorm`). The driver actually introduces a third CLI flag `--no-domain-pv` (lines 101-103, 113, 122 of `tpv102_driver.cpp`) together with a new `ParaViewOutput::ShouldWrite` helper (lines 649-672 of `paraview_output.hpp`). Neither is mentioned in the debug document. These are nontrivial behaviour changes (they touch the fault-surface VTU scheduling path) and should either be (a) documented in the v1 doc or (b) split into a separate commit so the "debug v1" record is actually auditable.

Scope creep is flagged MODERATE here because it affects the auditability of the debug record, not runtime correctness — but it would bite whoever tries to reproduce the v1 experiment from the doc alone.

**Trigger:**
Review of the v1 doc vs the commit diff.

**Actual behavior:**
Commit contains three new CLI flags (`--debug-qnorm`, `--no-domain-pv`, `--pv-low-order`). The doc mentions only `--debug-qnorm`.

**Expected behavior:**
Doc either mentions all three flags under "Fixes/changes landed" with rationale, or the commit is split into two.

**Suggested fix:**
Amend `tpv102_debug_v1.md` "Fixes landed in this commit" section to add:

```markdown
### F3.  `--no-domain-pv` + `ParaViewOutput::ShouldWrite`

Skip writing the volume PVD (ParaView/Cycle*/*.vtu, ~5 GB/cycle on the
200 m mesh) while keeping the fault-surface VTUs on the same schedule.
Saves tremendous disk on long coseismic runs, has no effect on debug
diagnostics.

### F4.  `--pv-low-order` + `ParaViewOutput::SetLevelsOfDetail`

Write VTK_TETRA (linear) instead of VTK_QUADRATIC_TETRA for the volume
PVD, ~2.5x cycle-size reduction on order-2 meshes.  No effect on
correctness.
```

**Test case:** N/A — documentation-only finding. Tracked via doc review checklist, not unit test.

---

### [R-005] [MODERATE] [POSSIBLE] `dynamic/wave_operator.inl:ComputeSharedFaceFluxRHS` — `nbr_data[c] = q_gf.FaceNbrData()` may silently alias (H2 unverified)

**Category:** ASSUMPTION

**Description:**
`paraview_output.hpp:622` stores the result of `q_gf.FaceNbrData()` into `std::vector<Vector>`. MFEM's `Vector::operator=` should deep-copy, but:
- `ParGridFunction::face_nbr_data` is a `mutable Vector` internal member, returned by reference from `FaceNbrData()`.
- The next iteration (`c+1`) reloads `q_gf` and calls `ExchangeFaceNbrData()`, which overwrites the *same* `face_nbr_data` buffer.
- If at any point MFEM's `Vector::operator=` takes the "view / MakeRef" branch (not "own and deep-copy"), every `nbr_data[c]` ends up pointing at the same buffer containing the LAST component's ghost data — a silent 9-way aliasing that would cause cross-component corruption consistent with the "Q literally zero far from hypocenter" v1 symptom.

This is the exact H2 hypothesis from the debug doc. The diagnostic is trivial (1 line) and should be added unconditionally.

**Trigger:**
Any parallel run. Triggering is a property of MFEM's `Vector::operator=` implementation; since we cannot trivially verify its branch behaviour from here, add a runtime assert.

**Actual behavior:**
Code assumes deep copy without verification.

**Expected behavior:**
Either assert deep-copy at startup (one-shot, negligible cost), OR replace with an explicit copy that cannot be miscompiled into aliasing.

**Suggested fix:**
Replace the assignment with an explicit deep copy that takes ownership:

```diff
-for (int c = 0; c < NUM_STATE; c++)
-{
-   for (int i = 0; i < ndof_total_; i++)
-   {
-      q_gf[i] = Q_data[c * ndof_total_ + i];
-   }
-   q_gf.ExchangeFaceNbrData();
-   nbr_data[c] = q_gf.FaceNbrData();
-}
+for (int c = 0; c < NUM_STATE; c++)
+{
+   for (int i = 0; i < ndof_total_; i++)
+   {
+      q_gf[i] = Q_data[c * ndof_total_ + i];
+   }
+   q_gf.ExchangeFaceNbrData();
+   const Vector &src = q_gf.FaceNbrData();
+   nbr_data[c].SetSize(src.Size());   // force own-storage
+   std::memcpy(nbr_data[c].GetData(), src.GetData(),
+              src.Size() * sizeof(real_t));
+}
+// Runtime assert: the 9 buffers must not alias the q_gf internal storage.
+for (int c = 0; c < NUM_STATE; c++)
+{
+   MFEM_ASSERT(nbr_data[c].GetData() != q_gf.FaceNbrData().GetData(),
+               "nbr_data["<<c<<"] aliases q_gf.FaceNbrData() — H2 hit.");
+}
```

**Test case:**
```python
def test_R005_nbr_data_not_aliased():
    """Run parallel driver with --dump-nbr-data-ptrs; assert the 9
       components point at distinct buffers."""
    out = run_driver(nranks=4, flags=["--dump-nbr-data-ptrs"])
    ptrs = parse_nbr_data_ptrs(out)
    assert len(set(ptrs)) == 9, f"nbr_data aliasing: {ptrs}"
```

---

### [R-006] [LOW] `tpv102_debug_v1.md` — Stale line-number reference

**Category:** QUALITY

**Description:**
The debug doc says `station_writer.Flush() is only called in the driver's post-loop cleanup (tpv102_driver.cpp:799)`. After the same commit added `--debug-qnorm` (which inserted ~40 lines before the cleanup block), the actual line is `851-852`. Not a correctness issue, but line references in debug docs rot the moment any neighbouring code changes — better to point at the symbol, not the line.

**Trigger:**
Reading the debug doc and trying to navigate to the referenced line.

**Actual behavior:**
`tpv102_driver.cpp:799` is inside the `--debug-qnorm` MPI_Gather block, not the Flush call.

**Expected behavior:**
Reference the symbol (e.g., "`station_writer.Flush()` in the post-loop cleanup") rather than a specific line number.

**Suggested fix:**
```diff
-`station_writer.Flush()` is only called in the driver's post-loop cleanup
-(`tpv102_driver.cpp:799`).
+`station_writer.Flush()` is only called in the driver's post-loop cleanup
+(`tpv102_driver.cpp` § 9 "Summary").
```

**Test case:** N/A — documentation-only.

---

### [R-007] [LOW] `drivers/tpv102_driver.cpp:paraview_write` — Fault-field packing loop always runs, even on cycles that do not write

**Category:** QUALITY

**Description:**
`paraview_write` (lines 575-628) packs `pv_local_slip` / `pv_local_traction` / … from `dof_data[i]` on every call (every RK4 step), regardless of whether `ShouldWrite`/`Save` will actually emit a snapshot. For a typical setup (`--paraview-every 100`), ~99% of the packing loop iterations are thrown away. On a 227k-QP fault this is ~1.8 M pointless stores per step. Not a bug — but on Frontera scale, per-step wasted work adds up.

**Trigger:**
Any run with `paraview_step_interval > 1` or `paraview_dt > dt`.

**Actual behavior:**
```cpp
auto paraview_write = [&](int step_num, real_t time, real_t V_max) {
    if (!pv_out) { return; }
    // ... always copies velocity + packs fault arrays ...
    bool wrote = pv_out->Save(step_num, time, V_max);
    if (wrote) { pv_out->WriteFaultSurfaceVTU(...); }
};
```

**Expected behavior:**
Check `ShouldWrite` first; only then pack.

**Suggested fix:**
```diff
 auto paraview_write = [&](int step_num, real_t time, real_t V_max)
 {
    if (!pv_out) { return; }
+   // Skip packing on cycles that won't write anyway.
+   if (!pv_out->ShouldWrite(step_num, time, V_max)) { return; }
    // ... existing body ...
 };
```

Note: this requires `ShouldWrite` to become idempotent (not advance `last_write_time_` when called from here). Easier alternative: add a `PeekShouldWrite` that only reads.

**Test case:** N/A — performance only, not correctness.

---

## Summary
- Critical issues: 1 (R-001 — shared-fault DOFData divergence; root cause of v1 symptom)
- Moderate issues: 4 (R-002 debug-qnorm deviation, R-003 surface-station missed flush, R-004 undocumented scope creep, R-005 possible H2 aliasing)
- Low issues: 2 (R-006 stale doc line, R-007 wasted packing)
- Plan compliance: **PARTIAL** — F1 applied to the wrong/incomplete scope (missed `TPV102SurfaceStationWriter`); F2 implemented but degraded from doc spec; also two extra features landed (`--no-domain-pv`, `--pv-low-order`) that the doc doesn't document.
- Verdict: **FAIL — must fix before proceeding.** R-001 is the actual rupture-propagation bug the v1 doc is chasing. R-003 blocks the next wall-limited Frontera run from producing valid surface-station data.

## Unreviewed Areas
- Behaviour of MFEM's `Vector::operator=` when the RHS is a reference to `ParGridFunction::face_nbr_data` — flagged as R-005 [POSSIBLE] because direct verification requires stepping through MFEM internals.
- Actual orientation of `CalcOrtho(ftr->Face->Jacobian())` on shared faces (does MFEM guarantee same geometric direction on both ranks?). R-001's reasoning assumes the same direction; the bug is worse if MFEM flips it for one rank, better if MFEM somehow arranges consistent `(+,−)` — but either way `Evaluate` is called with swapped Q arguments so R-001 stands. A ground-truth test on a 2-element shared-fault mesh would settle it.
- Performance of per-step `MPI_Gather` in `--debug-qnorm` at 400 ranks (negligible in practice; not investigated further).
