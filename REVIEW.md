# Code Review: split-bulk-solutions refactor (commit 7c16e62)

**Date:** 2026-05-12 (supersedes 2026-05-10 review)

## Review Scope
- Plan: `miniapps/seas/document/io_dev/PLAN_split_bulk_solutions_2026-05-12.md`
- Files reviewed:
  - `miniapps/seas/drivers/tpv102_driver.cpp`
  - `miniapps/seas/drivers/tpv104_driver.cpp`
  - `miniapps/seas/drivers/tpv205_driver.cpp`
  - `miniapps/seas/tests/verification/bp5_verification_full.cpp`
  - `miniapps/seas/io/fault_vtkhdf_writer.hpp`
  - `miniapps/seas/io/paraview_output.hpp` (read-only — not modified by the commit but referenced by the plan)
  - `miniapps/seas/tests/unit/test_fault_surface_vtkhdf{,_mpi,_zfp}.cpp`
  - `miniapps/seas/document/codebase_explorer/paraview_io_walkthrough_2026-05-11.md`
- Domain context: `miniapps/seas/CLAUDE.md`, the plan doc, `dynamic/wave_state.hpp` (Q-vector ordering)

## Findings

### [R-001] [MODERATE] [tpv102/tpv104/tpv205/bp5 drivers] — `kinematics.vtkhdf` carries 12 unwanted L2-p0 fault projections; plan + commit message + walkthrough doc all claim it carries only velocity + mpi_rank

**Category:** DEVIATION

**Description:**
The plan, commit message, and walkthrough doc all state that `kinematics.vtkhdf` holds only velocity + mpi_rank (TPV*) or displacement + mpi_rank (BP5). The actual code still calls `pv_out->InitFaultOutputBP5(...)` after registering velocity + mpi_rank. `InitFaultOutputBP5` registers **12 additional fields** with `pv_dc_` (paraview_output.hpp:784–795):

```
slip_dip, slip_strike, slip_rate_dip, slip_rate_strike,
traction_dip, traction_strike, state_variable, normal_stress,
param_a, param_Dc, fault_x2, fault_x3
```

All four drivers do this (TPV102 line 1934, TPV104 ~line 1880, TPV205 ~line 1957, BP5 line 1723). So every `kinematics.vtkhdf` contains **14 fields**, not 2. The 12 extra fields duplicate data already in `fault.vtkhdf`.

The commit message even contradicts itself:
> "kinematics.vtkhdf holds velocity + mpi_rank only. No fault L2-p0 projections (still allocated for now, but a follow-up commit will drop them...)"

"No fault L2-p0 projections" and "still allocated" cannot both be true. The fields aren't merely allocated — they're allocated AND registered, so they're written to the .vtkhdf on every save.

**Trigger:**
Run any TPV* or BP5 driver, dump the resulting `kinematics.vtkhdf` with `h5ls -r`, count the fields under `/VTKHDF/PointData` (or `CellData`).

**Actual behavior:**
`kinematics.vtkhdf` contains 14 fields per timestep. The 12 L2-p0 projections are duplicates of data already in `fault.vtkhdf`.

**Expected behavior (per the plan / commit message / walkthrough doc):**
`kinematics.vtkhdf` contains exactly velocity + mpi_rank (TPV) or displacement + mpi_rank (BP5).

**Suggested fix:**
Gate the registrations in `InitFaultOutputBP5` behind a new opt-in toggle. Default OFF means the projections allocate (the accessors `GetFaultSlipDip()` etc. keep working in unit tests that read them directly) but are NOT registered with `pv_dc_`, so they don't appear in `kinematics.vtkhdf`.

In `miniapps/seas/io/paraview_output.hpp` add a member around line 1646 (next to the other fault-output member state):

```diff
@@ paraview_output.hpp:~1648 @@
    std::unique_ptr<GF> fault_param_a_;
+
+   /// PLAN_split_bulk_solutions_2026-05-12: opt-in toggle for
+   /// registering the 12 L2-p0 fault projections with the primary
+   /// (kinematics) collection.  Default OFF — the projections
+   /// duplicate data in fault.vtkhdf and used to be the largest
+   /// contributor to volume.vtkhdf file size.  Enable via
+   /// `SetRegisterFaultProjectionsInVolumePV(true)` if you really
+   /// want them rendered alongside velocity / displacement.
+   bool register_fault_projections_in_volume_pv_ = false;
+
+public:
+   void SetRegisterFaultProjectionsInVolumePV(bool enable)
+   { register_fault_projections_in_volume_pv_ = enable; }
+   bool GetRegisterFaultProjectionsInVolumePV() const
+   { return register_fault_projections_in_volume_pv_; }
+private:
```

Then wrap the 12 registrations at paraview_output.hpp:784–795:

```diff
@@ paraview_output.hpp:783 @@
-      // Register fields (Phase 6.1: forwards through abstract `pv_dc_`).
-      pv_dc_->RegisterField("slip_dip",         fault_slip_dip_.get());
-      pv_dc_->RegisterField("slip_strike",      fault_slip_strike_.get());
-      pv_dc_->RegisterField("slip_rate_dip",    fault_slip_rate_dip_.get());
-      pv_dc_->RegisterField("slip_rate_strike", fault_slip_rate_strike_.get());
-      pv_dc_->RegisterField("traction_dip",     fault_trac_dip_.get());
-      pv_dc_->RegisterField("traction_strike",  fault_trac_strike_.get());
-      pv_dc_->RegisterField("state_variable",   fault_state_.get());
-      pv_dc_->RegisterField("normal_stress",    fault_normal_stress_.get());
-      pv_dc_->RegisterField("param_a",          fault_param_a_.get());
-      pv_dc_->RegisterField("param_Dc",         fault_param_Dc_.get());
-      pv_dc_->RegisterField("fault_x2",         fault_coord_x2_.get());
-      pv_dc_->RegisterField("fault_x3",         fault_coord_x3_.get());
+      // PLAN_split_bulk_solutions_2026-05-12: register the 12 L2-p0
+      // projection fields with the kinematics collection ONLY when
+      // the opt-in toggle is set.  Default OFF — these duplicate the
+      // data in fault.vtkhdf and used to be the largest contributor
+      // to volume.vtkhdf file size.  The GFs are still allocated so
+      // accessors / direct reads keep working.
+      if (register_fault_projections_in_volume_pv_)
+      {
+         pv_dc_->RegisterField("slip_dip",         fault_slip_dip_.get());
+         pv_dc_->RegisterField("slip_strike",      fault_slip_strike_.get());
+         pv_dc_->RegisterField("slip_rate_dip",    fault_slip_rate_dip_.get());
+         pv_dc_->RegisterField("slip_rate_strike", fault_slip_rate_strike_.get());
+         pv_dc_->RegisterField("traction_dip",     fault_trac_dip_.get());
+         pv_dc_->RegisterField("traction_strike",  fault_trac_strike_.get());
+         pv_dc_->RegisterField("state_variable",   fault_state_.get());
+         pv_dc_->RegisterField("normal_stress",    fault_normal_stress_.get());
+         pv_dc_->RegisterField("param_a",          fault_param_a_.get());
+         pv_dc_->RegisterField("param_Dc",         fault_param_Dc_.get());
+         pv_dc_->RegisterField("fault_x2",         fault_coord_x2_.get());
+         pv_dc_->RegisterField("fault_x3",         fault_coord_x3_.get());
+      }
```

No driver-side change required: by default the fields stop appearing in `kinematics.vtkhdf`. The toggle exists for anyone who wants to opt back in (e.g. unit tests that inspect the volume mesh's fault projections).

**Test case:**
```cpp
// miniapps/seas/tests/unit/test_kinematics_field_set.cpp
TEST(KinematicsContents, OnlyVelocityAndMpiRankByDefault) {
   // Stand up a tiny ParaViewOutput on a serial mesh; call
   // RegisterDomainField("velocity") + RegisterDomainField("mpi_rank")
   // + InitFaultOutputBP5(...); Save once; open the resulting .vtkhdf
   // and verify the field-name set.
   std::set<std::string> registered =
       PvOutFieldNames(pv);   // helper that walks pv->pv_dc_->GetFieldMap()
   EXPECT_EQ(registered, std::set<std::string>{"velocity", "mpi_rank"});

   pv.SetRegisterFaultProjectionsInVolumePV(true);
   pv.InitFaultOutputBP5(/*reinit*/...);  // or expose a re-register API
   registered = PvOutFieldNames(pv);
   EXPECT_GT(registered.size(), 13);  // 2 + 12
}
```

---

### [R-002] [MODERATE] [paraview_io_walkthrough_2026-05-11.md] — Walkthrough doc claims duplication is "RESOLVED" but L2-p0 projections still duplicate fault.vtkhdf

**Category:** DEVIATION

**Description:**
The "Gotchas" section reads:

> **(RESOLVED 2026-05-12)** ... kinematics.vtkhdf (renamed from volume.vtkhdf) holds velocity + mpi_rank (no stresses, no fault projections).

This is false. The 12 L2-p0 fault projections (R-001) are still registered. The doc misleads any future reader into thinking the cleanup is complete.

**Trigger:**
A reader of the walkthrough doc relies on the RESOLVED note to conclude that no duplication remains between `kinematics.vtkhdf` and `fault.vtkhdf`.

**Actual behavior:**
The 12 L2-p0 fields in `kinematics.vtkhdf` duplicate the per-DOF fault data already in `fault.vtkhdf`.

**Expected behavior:**
Doc accurately describes what's in each file.

**Suggested fix:**
If R-001 fix lands, the doc is correct as-is (the "no fault projections" claim becomes true). If R-001 is deferred, retract the RESOLVED claim:

```diff
- - **(RESOLVED 2026-05-12)** Earlier versions of the codebase wrote `velocity`
+ - **(PARTIALLY RESOLVED 2026-05-12; full fix tracked as R-001)** Earlier
+   versions of the codebase wrote `velocity`
```

and adjust the "no fault projections" wording.

Since R-001 has a concrete one-file fix above, the recommendation is to apply R-001 and leave the doc claim accurate.

**Test case:** N/A — documentation finding; verified by reading the doc and comparing against `paraview_output.hpp:784–795`.

---

### [R-003] [LOW] [tpv102/tpv104/tpv205] — Stale doc-comment block references the old "wave_bulk" field list

**Category:** QUALITY

**Description:**
Each TPV driver has a comment block summarizing what `pv_bulk_out` contains. The comments still describe the pre-refactor field list:

- `tpv102_driver.cpp:1857`: `// - pv_bulk_out (output_dir/ParaView_bulk): velocity + sigma_yy + sigma_xy + sigma_xz`
- `tpv104_driver.cpp:1805`: same
- `tpv205_driver.cpp:1884`: same

After the refactor, `pv_bulk_out` is the stress collection (file = `stress.vtkhdf`) holding the 6-component symmetric tensor.

**Trigger:**
A reader of the driver source.

**Actual behavior:**
Comments mislead about what's actually written.

**Expected behavior:**
Comments match the actual `RegisterDomainField` calls a few hundred lines below.

**Suggested fix:**
For each TPV driver, replace the stale line. tpv102 example:

```diff
@@ tpv102_driver.cpp:1857 @@
-   //       - pv_bulk_out (output_dir/ParaView_bulk): velocity + sigma_yy +
-   //         sigma_xy + sigma_xz at coarser cadence (--paraview-bulk-dt).
+   //       - pv_bulk_out (output_dir/ParaView_bulk/stress.vtkhdf):
+   //         full symmetric stress tensor (sigma_xx/yy/zz/xy/xz/yz)
+   //         at coarser cadence (--paraview-bulk-dt).  velocity and
+   //         mpi_rank are NOT in this file (they live in pv_out
+   //         only after PLAN_split_bulk_solutions_2026-05-12).
```

Mirror to tpv104 and tpv205.

**Test case:** N/A — comment correctness; visual inspection.

---

### [R-004] [LOW] [fault_vtkhdf_writer.hpp:107, :234] — Stale doc-comments still say `<prefix>/fault_surface.vtkhdf`

**Category:** QUALITY

**Description:**
Two docstring lines in the fault writer header still reference the old filename:

- Line 107: `///                <prefix>/fault_surface.vtkhdf).  The chunk filter`
- Line 234: `///                `<prefix>/fault_surface.vtkhdf`.`

The actual collection name was renamed to "fault" at line 189; the on-disk path is `<prefix>/fault.vtkhdf`.

**Trigger:**
Reading the header docstrings.

**Actual behavior:**
Docstrings disagree with the code.

**Expected behavior:**
Docstrings reflect the actual filename.

**Suggested fix:**
```diff
@@ fault_vtkhdf_writer.hpp:107 @@
-   ///                <prefix>/fault_surface.vtkhdf).  The chunk filter
+   ///                <prefix>/fault.vtkhdf).  The chunk filter
```
```diff
@@ fault_vtkhdf_writer.hpp:234 @@
-///                `<prefix>/fault_surface.vtkhdf`.
+///                `<prefix>/fault.vtkhdf`.
```

**Test case:** N/A — doc-comment correctness.

---

### [R-005] [LOW] [paraview_output.hpp:376, :1343] — Stale comments still reference "fault_surface"

**Category:** QUALITY

**Description:**
Two doc-comments in `paraview_output.hpp`:

- Line 376: mentions "fault_surface" in passing in the constructor docstring.
- Line 1343: comment in `WriteFaultSurfaceVTU` saying `<prefix>/fault_surface.vtkhdf`.

The file was renamed; these comments are stale.

**Suggested fix:**
```diff
@@ paraview_output.hpp:376 @@
- /// "fault_surface").
+ /// "fault").
```
```diff
@@ paraview_output.hpp:1343 @@
-      //    and emits a single `<prefix>/fault_surface.vtkhdf` rather
+      //    and emits a single `<prefix>/fault.vtkhdf` rather
```

**Test case:** N/A — doc-comment correctness.

---

### [R-006] [LOW] [tpv102/tpv104/tpv205] — Banner text wording change ("Bulk collection" → "Stress collection") may break out-of-tree log parsers

**Category:** QUALITY

**Description:**
The runtime banner text changed from `"  Bulk collection: ON ..."` to `"  Stress collection: ON ..."`. Any out-of-tree log-parsing tool that greps `"Bulk collection: ON"` to verify the secondary writer is active will silently miss the new banner.

**Trigger:**
Out-of-tree post-run script that greps the SLURM stdout for the old wording.

**Actual behavior:**
The banner is now "Stress collection: ON ...".

**Expected behavior:**
Not a correctness bug — just a wording change. Documented here so any downstream user knows.

**Suggested fix:**
None required in this repo. Out-of-tree users update their greps.

**Test case:** N/A.

---

### [R-007] [POSSIBLE LOW] [test_fault_surface_vtkhdf_mpi.cpp] — MPI test was rebuilt but not executed under mpirun in the review

**Category:** ASSUMPTION

**Description:**
The serial test `seas_test_fault_surface_vtkhdf` was run and passed 45/45. The MPI variant (`seas_test_fault_surface_vtkhdf_mpi`) was rebuilt cleanly but `mpirun` was not available in the review shell. The rename from "fault_surface" → "fault" should work identically in serial and MPI (the path string is on a single rank), but unconfirmed.

**Trigger:**
Run on any environment with mpirun available.

**Actual behavior:**
Unknown.

**Expected behavior:**
Test should pass collectively.

**Suggested fix:**
The fix agent should attempt `make test-fault-surface-vtkhdf-mpi` (or invoke the binary via the seas Makefile's MPI runner). Treat as a sanity check, not a blocker.

**Test case:**
Already exists as `test_fault_surface_vtkhdf_mpi.cpp`; just needs to run.

---

## Summary
- Critical issues: **0**
- Moderate issues: **2** (R-001, R-002 — both about kinematics-vs-fault data duplication still being present despite the doc/commit claiming it's resolved)
- Low issues: **5** (R-003 through R-007)
- Plan compliance: **PARTIAL** — kinematics file contents don't match the plan; L2-p0 fault projections still registered.
- Verdict: **PASS WITH FIXES** — no critical bugs in the build or runtime behaviour; the file split + rename work and binaries build clean. Spec deviation R-001 should be resolved before users start consuming the new files (otherwise the kinematics file is still ~50–70% duplicate data and the walkthrough doc is misleading).

## Unreviewed Areas
- `miniapps/seas/document/io_dev/PLAN_split_bulk_solutions_2026-05-12.md` — the plan doc itself. Read-only.
- MPI behaviour of the renamed fault writer (R-007).
- The rest of `paraview_output.hpp` — only the InitFaultOutputBP5 registration block and surrounding state were inspected.
- The `seas_driver.cpp` TOML wiring — explicitly deferred per R-313 and the plan doc.
