# Code Review: 2026-05-17 round 2 — paraview-compaction merge audit (post R-001..R-009 fix sweep)

## Review Scope

- **Plan (implicit):** Verify (a) `feature/paraview-compaction` merged correctly into
  `feature/safs-quasi-dynamic`; (b) prior fix sweep (commit `1a805ef`, addressing
  R-001..R-009 from `REVIEW.md` round 1) did not regress anything and did not leave
  silent gaps; (c) no critical bugs remain on the merged surface area.
- **Files freshly re-audited (three-pass adversarial):**
  - `miniapps/seas/io/paraview_output.hpp` (2293 lines — full re-read of the post-merge
    state, including R-001..R-007 round-1/2/3 fixes baked into header)
  - `miniapps/seas/io/fault_vtkhdf_writer.hpp` (308 lines — Init invariants, ZFP plugin
    probe, geometry-static checks)
  - `miniapps/seas/io/fault_vtu_binary.hpp` (lines 1..500 + raw-appended emit;
    cross-rank field-name hash check)
  - `miniapps/seas/fault/fault_geometry.hpp` (the R-001 split-counter + t1-fallback
    accounting at lines 820..1010)
  - `miniapps/seas/fault/rate_state_fault.hpp` (SetSAFSMode guard + SAFS routing at
    lines 70..230, 370..520, 1130..1190)
  - `miniapps/seas/fault/fault_geometry_safs.inl` (ComputeSAFSParams; tau_pre_ overwrite
    contract)
  - `miniapps/seas/config/seas_config_bridge.hpp` (ApplySAFSMode; construct→Compute→Set
    ordering)
  - `miniapps/seas/domain/domain_operator.hpp` (base default for
    `ComputeTractionDiagnostics` — round-1 R-003 fix)
  - `miniapps/seas/solver/seas_operator.hpp` (Mult dispatch; face_tracer call sites at
    lines 339..380, 536..544)
  - `miniapps/seas/trace/face_trace_logger.hpp` (RecordMult body at lines 345..490 —
    indexing pattern audit, post-merge)
  - `miniapps/seas/drivers/tpv102_driver.cpp` (paraview_write lambda at lines
    2195..2298 — CommitSchedule 1-arg/2-arg dispatch)
  - `miniapps/seas/drivers/tpv104_driver.cpp` (same, lines 2230..2330)
  - `miniapps/seas/drivers/tpv205_driver.cpp` (same, lines 2180..2290)
  - `miniapps/seas/tests/unit/test_safs_mode_wiring.cpp` (T_66_1..T_66_5 — round-1 R-001
    positive control)
  - `miniapps/seas/tests/unit/test_R003_antiplane_diagnostics_sizing.cpp` (round-1
    R-003 base-default sizing test)
  - MFEM-core diffs in `mesh/vtkhdf.cpp,hpp` + `fem/datacollection.cpp,hpp`
    (HDFCompression selector, ZFP filter dispatch)
- **Domain context consulted:** `miniapps/seas/CLAUDE.md` (sign conventions, ParaView
  Phase 6 plan, SAFS no-hardcoding), root `CLAUDE.md` ("REPORT and ASK before falling
  back"), `feedback_complete_sign_sites.md`, `feedback_no_hardcoded_numbers.md`,
  `feedback_debug_docs_as_plans.md`.
- **Prior review consumed:** `REVIEW.md` round 1 (R-001..R-009) and the fix commit
  `1a805ef seas: address REVIEW.md R-001..R-009 findings`.  Findings R-001..R-009 are
  treated as resolved; this round audits the resulting code as a *fresh* artifact.

---

## Findings

### [R-201] [MODERATE] [drivers/tpv102_driver.cpp:2288, tpv104_driver.cpp:2316, tpv205_driver.cpp:2278] — 1-arg `CommitSchedule(time)` in fault-only PV branch silently breaks cadence with `hysteresis_factor > 1`

**Category:** BUG (latent; manifests only when CLI overrides the schedule)

**Description:**
The `paraview_write` lambda in all three TPV drivers branches on
`pv_out->GetVolumeSaveEnabled()`:

```cpp
if (pv_out->GetVolumeSaveEnabled())
{
   pv_out->UpdateFaultFieldsBP5(...);
   pv_out->ForceSave(step_num, time);
   pv_out->CommitSchedule(time, V_max);     // 2-arg form — correct
}
else
{
   pv_out->CommitSchedule(time);            // 1-arg form — uses stale last_v_max_
}
```

The 1-arg shim `CommitSchedule(time)` delegates to `CommitSchedule(time, last_v_max_)`.
Inside the fault-only path (`--no-volume-pv` and no `--volume-pv-dt`), `Save` /
`ShouldWrite` / `ForceSave` are never called by the caller, so `last_v_max_` is never
refreshed and stays at its previous value (defaults to 0.0 at construction).  As soon
as the driver flips into hysteresis territory by setting
`pv_out->GetSchedule().hysteresis_factor > 1.0`, `NextRegime(stale_V, current_regime_)`
no longer agrees with `NextRegime(actual_V, current_regime_)` — the regime state
machine drifts away from the actual `V_max` trajectory.

The header explicitly warns about this contract at `io/paraview_output.hpp:1707-1709`:

> Callers that want hysteresis to advance through a PeekShouldWrite-gated path MUST
> call the two-arg overload instead.

The drivers gate I/O on `PeekShouldWrite` (see `tpv104_driver.cpp:2234`,
`tpv102_driver.cpp:2198`, `tpv205_driver.cpp:2186`) — exactly the situation the header
forbids 1-arg for.  All four drivers currently use the default `hysteresis_factor = 1.0`
so the bug is dormant, but `AdaptiveSchedule::hysteresis_factor` is publicly settable
via `GetSchedule()` and the CLI does not block changing it; any future flag (or a TOML
config that ships) would trip silently.

The header text plus the existing 2-arg form one line above prove the author knew the
contract.  The 1-arg branch is a leftover from before `PeekShouldWrite` existed.

**Trigger:**
1. Set `pv_out->GetSchedule().hysteresis_factor = 2.0` (or any > 1.0).
2. Run with `--no-volume-pv`.
3. After a coseismic burst that drops back through the nucleation band, the fault-only
   schedule sticks at the post-burst regime that was computed against a stale `V_max`
   (0.0 by default → interseismic) rather than the true `V_max`.

**Actual behavior:**
Stale-`V_max`-driven `NextRegime` in fault-only mode; regime state lags the actual
slip-rate trajectory.

**Expected behavior:**
Both branches should use `CommitSchedule(time, V_max)`.  The 1-arg shim is appropriate
for legacy paths that already called `Save`/`ShouldWrite` immediately before — not for
this `PeekShouldWrite`-gated pattern.

**Suggested fix:**
```diff
--- a/miniapps/seas/drivers/tpv102_driver.cpp
+++ b/miniapps/seas/drivers/tpv102_driver.cpp
@@ -2286,5 +2286,7 @@
       else
       {
-         pv_out->CommitSchedule(time);
+         // R-201: PeekShouldWrite-gated path MUST use the 2-arg form so
+         // last_v_max_ / current_regime_ track the live trajectory.
+         pv_out->CommitSchedule(time, V_max);
       }
```

Identical diff applies to `tpv104_driver.cpp:2316` and `tpv205_driver.cpp:2278`.

**Test case:**
```cpp
// tests/unit/test_R201_fault_only_commit_uses_v_max.cpp
void test_R201_fault_only_path_records_live_v_max()
{
   auto mesh = MakeSmallBP5Mesh();
   seas::ParaViewOutput<mfem::Mesh> pv("out", *mesh, 1);
   pv.GetSchedule().hysteresis_factor = 2.0;   // enable hysteresis
   pv.GetSchedule().Validate();
   pv.SetTotalRunTime(1.0);
   pv.SetVolumeSaveEnabled(false);

   // Simulate paraview_write fault-only branch with the BUGGY 1-arg form:
   pv.CommitSchedule(0.5, 5e-2);   // V_max coseismic
   // Driver bug: 1-arg form ignores the V we just observed.
   pv.CommitSchedule(0.6);
   TEST_ASSERT_NEAR(pv.GetLastVMax(), 5e-2, 0.0,
                    "FIXED: last_v_max_ MUST equal the V passed to PeekShouldWrite, "
                    "not the previously-committed value.  If this is 0.0 the driver "
                    "fault-only branch silently downgraded the regime.");
}
```

---

### [R-202] [MODERATE] [trace/face_trace_logger.hpp:409-416] — `RecordMult` unconditionally uses `2*dof` indexing, defeating round-1 R-003 fix for Antiplane + face_tracer

**Category:** BUG (latent — composition is currently unreached, but R-003 specifically
called this out as a future combination that the base default should make safe)

**Description:**
Round-1 R-003 added a base-class default for
`DomainOperator<MeshType>::ComputeTractionDiagnostics` that sizes
`traction_stress` / `traction_correction` / `jump_residual` to match `traction.Size()`
with zeros, so downstream `FaceTraceLogger::RecordMult` can iterate without OOB on the
Antiplane + face_tracer composition.  The fix correctly sizes the vectors, but
`face_trace_logger.hpp:409-416` still indexes those vectors with BP5 2-component
strides:

```cpp
real_t tau_d = staged_traction_(2 * dof);
real_t tau_s = staged_traction_(2 * dof + 1);
real_t ts_d = staged_stress_(2 * dof);
real_t ts_s = staged_stress_(2 * dof + 1);
real_t tc_d = staged_corr_(2 * dof);
real_t tc_s = staged_corr_(2 * dof + 1);
real_t jr_d = staged_jump_res_(2 * dof);
real_t jr_s = staged_jump_res_(2 * dof + 1);
```

For Antiplane (`NumSlipComponents == 1`), `traction_.Size() == num_owned_dofs`.  After
R-003, the base default sets the three diagnostic vectors to that same size (zeros).
The tracer then reads `(2*dof)` up to `2*(num_dofs-1)` — twice the legal index range —
which is an OOB read on the second half of every owned DOF iteration.  Even with R-003
in place, instantiating `SEASQuasiDynamicOperator<Mesh, AntiplaneDomainOperator<>, ...>`
with an active `FaceTraceLogger` still crashes (or silently reads garbage) on the
first `RecordMult`.

Why this slipped past R-003: R-003's test (`test_R003_antiplane_diagnostics_sizing.cpp`)
verifies the *base default's* output sizing in isolation; it does not actually call
`face_tracer_->RecordMult` on Antiplane data.  No code currently composes Antiplane +
tracer, so this is latent.  But R-003 was explicitly motivated by "downstream
`FaceTraceLogger::RecordMult` reading past the end of an empty Vector", so the obvious
follow-up — make `RecordMult` honour the actual stride — was never done.

**Trigger:**
Construct `SEASQuasiDynamicOperator<Mesh, AntiplaneDomainOperator<>, RateStateFaultOperator<Mesh, 1>>`,
call `SetFaceTracer(...)` with an active tracer, call `Mult()`.

**Actual behavior:**
`staged_traction_(2*dof)` reads past the end of a 1-component-sized vector; on every
iteration for `dof >= n/2`.  Crash or wrong data.

**Expected behavior:**
Tracer should adapt to the actual slip-component count, or document that it is
BP5-only (and add a compile-time guard that prevents Antiplane instantiation).

**Suggested fix:**
The minimum-invasive fix is to make the tracer template-parameterise over slip
components (mirroring `RateStateFaultOperator`).  Until that lands, gate `RecordMult`
on a runtime stride check and bail with a clear error rather than crashing:

```diff
--- a/miniapps/seas/trace/face_trace_logger.hpp
+++ b/miniapps/seas/trace/face_trace_logger.hpp
@@ -345,6 +345,16 @@
    void RecordMult(
       const Vector &traction,
+      // R-202: Tracer assumes BP5 2-component layout. Bail loudly rather
+      // than OOB-read if a future caller wires Antiplane (1 comp) here.
+      // Note round-1 R-003 sized the diagnostics to traction.Size() for
+      // exactly this composition; the safe path is to refuse the call
+      // until the tracer is templated over SlipComponents.
       ...
    {
       if (!IsActive()) { return; }
+      const int n_owned = static_cast<int>(traced_dofs_.size());
+      MFEM_VERIFY(n_owned == 0 || traction.Size() >= 2 * n_owned,
+                  "FaceTraceLogger::RecordMult: traction.Size()="
+                  << traction.Size() << " < 2 * n_owned=" << 2*n_owned
+                  << ".  Tracer requires BP5 2-component layout; Antiplane "
+                  "(1 comp) instantiation must template-specialise the tracer "
+                  "or skip tracer wiring.");
```

The long-term fix is a `template <typename MeshType, int SlipComponents = 2>` on
`FaceTraceLogger` so the indexing math is generated correctly for each path.

**Test case:**
```cpp
// tests/unit/test_R202_antiplane_tracer_refuses.cpp
void test_R202_antiplane_tracer_indexing_refused()
{
   auto mesh = BP2MeshGenerator::Create({...});
   AntiplaneDomainOperator<Mesh> domain(*mesh, /*order=*/1, ...);
   Vector traction(domain.GetNumFaultDOFs());   // size = 1 * n_dofs
   traction = 0.0;
   FaceTraceLogger<Mesh> tracer(cfg, /*rank=*/0);
   // Register a single owned DOF so n_owned > 0 and the guard triggers.
   tracer.RegisterDof(...);
   bool aborted = false;
   try { tracer.RecordMult(traction, traction, traction, traction,
                           Vector(), Vector(), Vector(), traction, 0.0); }
   catch (...) { aborted = true; }
   TEST_ASSERT(aborted,
               "RecordMult must refuse 1-component traction or template-specialise");
}
```

---

### [R-203] [MODERATE] [fault/rate_state_fault.hpp:214-223] — SAFS silent-corruption guard only blocks zero-normal fallbacks; t1-degenerate slots have undefined sign that flips SAFS pre-stress

**Category:** ASSUMPTION (sign-of-shear may flip silently per DOF)

**Description:**
Round-1 R-001 added a guard inside `SetSAFSMode` that aborts when
`FaultGeometry::NumZeroNormalFallbacks() > 0`.  The companion class doc at
`fault_geometry.hpp:330-358` separates fallbacks into two kinds:

- **zero-normal** — basis slot is zeroed; SAFS would silently project to zero.
- **t1-degenerate** — basis slot is valid *but the FaultBasis sign convention is
  undefined for that DOF*.

The R-001 guard rejects only the first kind.  The class doc says t1-degenerate "may
proceed with a runtime warning" — but the actual runtime warning at
`fault_geometry.hpp:1000-1007` only flags **zero-normal** as a silent-corruption
hazard:

```text
... The zero-normal slots have zeroed basis and would silently project
to zero pre-stress and zero sigma_n in SAFS mode.  SetSAFSMode refuses
to enable when zero-normal fallbacks > 0.
```

It does not warn that t1-degenerate DOFs may have *sign-flipped* pre-stress
projections.  The repo CLAUDE.md is explicit about this hazard:

> Slip rate direction: PARALLEL to traction tau (not antiparallel).  Antiparallel sign
> creates positive feedback -> unbounded growth (debug v8).

A flipped basis t2 at a single DOF projects the sidecar stress tensor with reversed
sign in the strike (or dip) tangent direction, putting that DOF into the antiparallel
regime — the exact failure mode CLAUDE.md flags as causing unbounded growth.

This is consistent with the round-1 R-101 / R-003 commentary at
`fault_geometry.hpp:961-976`:

> The cross-product form `t2 = n × t1` *loses* the FaultBasis sign-flip convention
> (Tandem: `t2_face = -strike_can` when the CalcOrtho normal is anti-aligned with
> `ref_normal_`...).  Dropping the sign-flip would place the SAFS pre-stress in a
> different frame from the elastic traction...

The author understood the sign-flip risk while writing the t2 recovery; the SAFS
guard didn't follow through.

**Trigger:**
A SAFS run on a mesh whose `GetFaultDOFBasis` produces any DOF with `t1` magnitude
below 1e-12 after Gram-Schmidt against the dip-projected up-vector — common on
curvilinear faults near the strike/dip-aligned corners that the SAFS plan flags as
production targets.

**Actual behavior:**
SetSAFSMode succeeds; ComputeRHS at the t1-fallback DOF reads a sidecar projection
whose tangential components have sign opposite to what FaultBasis would produce
elsewhere; slip rate at that DOF acquires antiparallel feedback → simulation either
diverges or quietly produces wrong slip direction.

**Expected behavior:**
Either (a) extend the guard to also reject `NumT1Fallbacks() > 0` (the conservative
fix that mirrors R-001's "refuse rather than corrupt silently" intent), OR (b)
re-derive `t2` on the fallback path so it preserves the FaultBasis sign convention,
OR (c) at minimum upgrade the runtime warning to explicitly flag the sign-flip risk
for t1-fallback DOFs.

The most defensible fix is (a) — the warning text below already pushes the user
toward fixing the upstream operator rather than tolerating any fallbacks.

**Suggested fix:**
```diff
--- a/miniapps/seas/fault/rate_state_fault.hpp
+++ b/miniapps/seas/fault/rate_state_fault.hpp
@@ -213,12 +213,18 @@
             // allowed; user is warned via FaultGeometry's stderr message.
             MFEM_VERIFY(geom_ == nullptr ||
-                        geom_->NumZeroNormalFallbacks() == 0,
+                        (geom_->NumZeroNormalFallbacks() == 0 &&
+                         geom_->NumT1Fallbacks() == 0),
                         "SetSAFSMode: cannot enable SAFS mode — "
                         "FaultGeometry has "
                         << geom_->NumZeroNormalFallbacks()
-                        << " fault DOF(s) with a zeroed basis (zero-length "
-                        "normal during ComputePerDOFCoordsAndBasis_).  "
-                        "Their sidecar pre-stress and sigma_n would silently "
-                        "project to zero.  Fix the operator's "
-                        "GetFaultDOFBasis to populate all owned fault DOFs.");
+                        << " zero-normal + " << geom_->NumT1Fallbacks()
+                        << " t1-degenerate fault DOF(s).  zero-normal DOFs "
+                        "would silently project sidecar pre-stress and sigma_n "
+                        "to zero; t1-degenerate DOFs have a sign-flip-"
+                        "undefined basis that can place the projected stress "
+                        "ANTIPARALLEL to the elastic traction at that DOF, "
+                        "creating the positive-feedback regime CLAUDE.md "
+                        "documents (debug v8).  Fix the operator's "
+                        "GetFaultDOFBasis (zero-normal) and the upstream "
+                        "mesh / up-vector choice (t1-degenerate) before "
+                        "enabling SAFS.");
```

And update the runtime warning text in `fault_geometry.hpp:1000-1007`:

```diff
-         mfem::err << "FaultGeometry WARNING: "
-                   << num_zero_normal_fallbacks_ << " zero-normal + "
-                   << num_t1_fallbacks_ << " t1-degenerate (of "
-                   << num_owned << " owned fault DOFs) hit the basis "
-                   << "fallback path. The zero-normal slots have zeroed "
-                   << "basis and would silently project to zero pre-stress "
-                   << "and zero sigma_n in SAFS mode. SetSAFSMode refuses "
-                   << "to enable when zero-normal fallbacks > 0.\n";
+         mfem::err << "FaultGeometry WARNING: "
+                   << num_zero_normal_fallbacks_ << " zero-normal + "
+                   << num_t1_fallbacks_ << " t1-degenerate (of "
+                   << num_owned << " owned fault DOFs) hit the basis "
+                   "fallback path.  Zero-normal slots have a zeroed basis "
+                   "and would silently project to zero pre-stress + sigma_n "
+                   "in SAFS mode.  t1-degenerate slots have a sign-flip-"
+                   "undefined basis that can place the projected sidecar "
+                   "stress ANTIPARALLEL to the elastic traction, creating "
+                   "the positive-feedback failure documented in CLAUDE.md "
+                   "(debug v8).  SetSAFSMode refuses to enable in either "
+                   "case (see rate_state_fault.hpp:214).\n";
```

**Test case:**
Extend `test_safs_mode_wiring.cpp::T_66_5` with a negative-control fork test that
constructs a fixture with one t1-degenerate DOF and verifies `SetSAFSMode(true, ...)`
aborts:

```cpp
static void T_66_6_safs_mode_guard_refuses_t1_fallback(const Fixture &fix)
{
   // Construct a synthetic FaultGeometry with one t1-degenerate DOF
   // (use the same mechanism the production code triggers — fault face
   // exactly parallel to the up_ref vector so the projected dip vanishes).
   auto degen = MakeFaultGeometryWithT1Fallback(/*degen_idx=*/0);
   TEST_ASSERT(degen.NumT1Fallbacks() == 1,
               "Fixture has one t1-degenerate DOF (positive control)");
   RateStateFaultOperator<Mesh, 2> op(&degen, friction_.get(), evolution_.get(),
                                       params_);
   Vector tau_pd(2 * degen.NumFaultDOFs()), sn_pd(degen.NumFaultDOFs());
   tau_pd = 0.0; sn_pd = 50e6;
   // After R-203 fix: must abort because t1-fallback may flip stress sign.
   bool aborted = false;
   try { op.SetSAFSMode(true, &tau_pd, &sn_pd); }
   catch (...) { aborted = true; }  // MFEM_VERIFY → abort → caught in fork test
   TEST_ASSERT(aborted, "SetSAFSMode refuses t1-degenerate fixture (R-203)");
}
```

---

### [R-204] [LOW] [io/paraview_output.hpp:1557, 2247, 2253-2268] — `fault_pvd_entries_` grows unbounded; `WriteFaultPVD` rewrites the whole file every cycle → O(N²) cumulative I/O

**Category:** QUALITY (performance / disk write amplification at scale)

**Description:**
Each successful fault-surface write on rank 0 pushes one `PVDEntry { time, file }` into
the `fault_pvd_entries_` `std::vector` (lines 1557, 2247) and then re-emits the entire
PVD index by truncating and rewriting the file (lines 2253-2268).  With the BP5
production cap of `max_total_snapshots = 5000` (per
`bp5_debug_document/regime_budget_2026-05-16.md`), the cumulative PVD I/O over a run is:

  sum_{k=1..5000} k * ~80 bytes = ~1.0 GB of rewrites of a file that ends at ~400 KB.

For short runs (~hundred writes) this is invisible.  For TPV/SAFS production runs that
push the cap to 5000+, it's measurable on shared scratch.  More importantly, the
truncate+rewrite is not atomic — a job killed mid-write leaves a truncated PVD that
ParaView refuses to open, and the entries vector is the only authoritative copy.

There is no checkpoint of `fault_pvd_entries_` either (the V2 checkpoint at
`io/petsc_ts_checkpoint.hpp` does not capture this state); a restarted run starts the
list from empty and the new PVD only references post-restart entries.

**Trigger:**
Production run with `--paraview-max-snapshots 5000` and `--paraview-fault-hdf5`
disabled (so the VTU path's PVD is in play), or any TPV run with `paraview_max_snapshots`
set high.

**Actual behavior:**
~1 GB cumulative rank-0 PVD writes; possible truncation on kill mid-write; PVD
"forgets" pre-restart entries.

**Expected behavior:**
Append-only PVD writer that opens for append once and emits one `<DataSet>` line per
call, or at minimum (a) check-pointed `fault_pvd_entries_`, and (b) rename-after-write
for atomicity.

**Suggested fix:**
Atomic write-and-rename:

```diff
--- a/miniapps/seas/io/paraview_output.hpp
+++ b/miniapps/seas/io/paraview_output.hpp
@@ -2253,15 +2253,21 @@
    void WriteFaultPVD(const std::string &fault_dir)
    {
-      std::string pvd_name = fault_dir + "/fault_surface.pvd";
-      std::ofstream pvd(pvd_name, std::ios::trunc);
+      // R-204: write to .partial and atomic-rename so a job killed
+      // mid-write does NOT leave a truncated PVD that ParaView rejects.
+      // (The O(N^2) cumulative rewrite is unchanged — the followup is a
+      // streaming append, but the atomic rename costs nothing and fixes
+      // the corruption-on-kill hazard immediately.)
+      const std::string pvd_name = fault_dir + "/fault_surface.pvd";
+      const std::string tmp_name = pvd_name + ".partial";
+      std::ofstream pvd(tmp_name, std::ios::trunc);
       pvd << std::setprecision(17);
       pvd << "<?xml version=\"1.0\"?>\n";
       pvd << "<VTKFile type=\"Collection\" version=\"0.1\">\n";
       pvd << "<Collection>\n";
       for (const auto &e : fault_pvd_entries_)
       {
          pvd << "<DataSet timestep=\"" << e.time
              << "\" file=\"" << e.file << "\"/>\n";
       }
       pvd << "</Collection>\n</VTKFile>\n";
       pvd.close();
+      ::rename(tmp_name.c_str(), pvd_name.c_str());  // atomic on POSIX
    }
```

The atomic rename via `tmp + rename` is the minimum acceptable fix and is much smaller
than refactoring to a streaming append; it fixes the kill-mid-write hazard at zero
algorithmic cost.  The full N→N+1 streaming append is a follow-up.

**Test case:**
Inspection / runtime measurement only — there is no functional regression to assert.
Note as known scaling limit in `KNOWN_DISABLED_TESTS.md` once a guard is added.

---

### [R-205] [LOW] [fault/fault_geometry.hpp:840-845] — `ComputePerDOFCoordsAndBasis_` early-return path skips fallback-counter reset

**Category:** ASSUMPTION (maintainability landmine; not a current bug)

**Description:**
The early-return branch at line 840-845 handles the antiplane-fallback case (operator
did not populate `full_coords` / `full_basis`):

```cpp
if (full_coords.Size() == 0 || full_basis.Height() == 0)
{
   dof_coords_3d_.SetSize(0);
   dof_basis_.SetSize(0, 0);
   return;
}
```

It does **not** reset `num_dof_basis_fallbacks_`, `num_zero_normal_fallbacks_`, or
`num_t1_fallbacks_`.  These are class-member initialised to 0
(`fault_geometry.hpp:548-550`), so on the first invocation the counters happen to be
0 by accident — the SAFS guard at `rate_state_fault.hpp:214` passes.  No bug *today*.

But the explicit reset block at lines 883-885 inside the main loop is documentation of
intent: "callers should observe a freshly-reset counter".  If a future refactor
re-uses this function (or calls it twice, e.g. after `domain_op->PopulateFaultBasis()`
is added), the early-return would silently carry over a stale nonzero counter, and the
SAFS guard would incorrectly abort or pass.  Move the reset to the function preamble.

**Trigger:**
Future refactor that calls `ComputePerDOFCoordsAndBasis_` more than once, with one
invocation taking the early-return path between two main-loop invocations.

**Actual behavior:**
Counters leak across calls; SAFS guard sees stale value.

**Expected behavior:**
Counters reset on every call regardless of which branch taken.

**Suggested fix:**
```diff
--- a/miniapps/seas/fault/fault_geometry.hpp
+++ b/miniapps/seas/fault/fault_geometry.hpp
@@ -830,6 +830,12 @@
    void ComputePerDOFCoordsAndBasis_(DomainOperator<MeshType> &domain_op)
    {
+      // R-205: reset counters BEFORE the early-return branch so the
+      // antiplane-fallback path is observably "no fallbacks" rather
+      // than carrying over whatever a previous invocation left.
+      num_dof_basis_fallbacks_   = 0;
+      num_zero_normal_fallbacks_ = 0;
+      num_t1_fallbacks_          = 0;
       // Gather full local (interior + shared) coords + basis from the operator
       Vector full_coords;
       DenseMatrix full_basis;
@@ -881,9 +887,5 @@
       // R-006 Gram-Schmidt re-orthonormalisation, mirroring Phase 2
       // `basis_to_node` (project_to_fault_stress.py).
-      num_dof_basis_fallbacks_   = 0;
-      num_zero_normal_fallbacks_ = 0;
-      num_t1_fallbacks_          = 0;
       for (int i = 0; i < num_owned; i++)
```

**Test case:**
Not directly testable without the future refactor that re-uses the function.
Annotate the preamble reset with `// R-205 — keep at function entry` to prevent
regression.

---

### [R-206] [LOW] [fault/rate_state_fault.hpp:156 + fault/fault_geometry_safs.inl:58] — `tau_pre_` cached at operator construction; if `ComputeSAFSParams` runs before construction, BP5-mode reads silently use sidecar values

**Category:** ASSUMPTION (brittle construction-order contract; not enforced)

**Description:**
The BP5 constructor copies the FaultGeometry's `tau_pre_` into its own member at
construction time:

```cpp
// rate_state_fault.hpp:155-157
Dc_values_     = geom_->GetDcValues();
tau_pre_       = geom_->GetTauPre();      // ← cached
V_init_values_ = geom_->GetVInit();
```

`ComputeSAFSParams` later *overwrites* `geom_->tau_pre_` with the sidecar projection
(`fault_geometry_safs.inl:55-61`):

```cpp
FieldProjector::ProjectFaultPreStress<MeshType>(
   field, *this, sigma_n_per_dof_, tau_pre_, ...);   // ← overwrites geom_->tau_pre_
```

This works correctly only if the canonical order
`operator-ctor → ComputeSAFSParams → SetSAFSMode` is followed: operator's cached
`tau_pre_` keeps the BP5 analytic; geom's `tau_pre_` becomes sidecar; SAFS mode reads
the latter via the pointer.

If a user reverses the order — `ComputeSAFSParams → operator-ctor → SetSAFSMode` —
the operator's cached `tau_pre_` is the *sidecar projection*, and any toggle back to
BP5 mode (`SetSAFSMode(false)`) silently reads sidecar pre-stress while pretending to
be in BP5 mode.  No assertion catches this.  `ApplySAFSMode` happens to call them in
the right order, so production paths are fine, but tests that try
`SetSAFSMode(false)` after `ComputeSAFSParams` (e.g. for bit-exact comparison) would
not get BP5 behaviour.

**Trigger:**
Direct test code that calls `geom.ComputeSAFSParams(...)` BEFORE constructing the
`RateStateFaultOperator`, then constructs the operator and toggles SAFS mode off for
comparison.

**Actual behavior:**
BP5-mode (safs_mode_=false) reads cached `tau_pre_` which is actually the sidecar
projection.

**Expected behavior:**
Either (a) document the contract loudly in the constructor (in addition to the
existing doc at line 173-175 which only mentions pointer lifetime, not construction
order), OR (b) add a `MFEM_VERIFY(!geom_->HasSAFSParams())` to the constructor to
abort if the geom has been pre-projected, OR (c) re-fetch `tau_pre_` inside
`SetSAFSMode(false)` to undo the cache when toggling.

Recommended: option (b) — explicit guard.

**Suggested fix:**
```diff
--- a/miniapps/seas/fault/rate_state_fault.hpp
+++ b/miniapps/seas/fault/rate_state_fault.hpp
@@ -137,6 +137,14 @@
       MFEM_ASSERT(geom_ != nullptr && geom_->IsBP5(),
                   "BP5 constructor requires BP5 FaultGeometry");
+      // R-206: this constructor caches geom_->tau_pre_ for the BP5
+      // (non-SAFS) hot path.  ComputeSAFSParams *overwrites* that
+      // member with the sidecar projection, so constructing the
+      // operator AFTER ComputeSAFSParams would silently cache the
+      // sidecar values and make SetSAFSMode(false) read them — the
+      // opposite of what the user expects.  Abort early to enforce
+      // the canonical order: construct → ComputeSAFSParams → SetSAFSMode.
+      MFEM_VERIFY(!geom_->HasSAFSParams(),
+                  "RateStateFaultOperator: BP5 ctor must run BEFORE "
+                  "FaultGeometry::ComputeSAFSParams; otherwise BP5-mode "
+                  "tau_pre_ caches sidecar values silently.");
 
       if (num_nodes_ > 0)
```

(`HasSAFSParams()` already exists per `fault_geometry.hpp:365`.)

**Test case:**
```cpp
void test_R206_construct_after_compute_safs_aborts()
{
   auto fix = BuildBP5Fixture();
   fix.fault_geom->ComputeSAFSParams(...);   // overwrite tau_pre_ first
   bool aborted = false;
   try { auto op = std::make_unique<RateStateFaultOperator<Mesh, 2>>(
                     fix.fault_geom.get(), fix.friction.get(),
                     fix.evolution.get(), fix.params); }
   catch (...) { aborted = true; }
   TEST_ASSERT(aborted,
               "BP5 ctor must abort after ComputeSAFSParams (R-206)");
}
```

---

### [R-207] [LOW] [io/fault_vtu_binary.hpp:431-441] — Cross-rank `std::hash<std::string>` comparison assumes a process-deterministic hash policy

**Category:** ASSUMPTION (works on libstdc++/libc++ today; implementation-defined)

**Description:**
The gather's field-name uniformity check broadcasts a `std::hash<std::string>` value
from rank 0 and compares each rank's local hash:

```cpp
std::string concat;
for (const auto &n : local.field_names) { concat += n; concat += '|'; }
const uint64_t local_hash = std::hash<std::string>{}(concat);
uint64_t root_hash = local_hash;
MPI_Bcast(&root_hash, 1, MPI_UNSIGNED_LONG_LONG, 0, comm);
MFEM_VERIFY(local_hash == root_hash, ...);
```

`std::hash<std::string>` is intentionally not specified by the standard to be stable
across processes or releases.  libstdc++ and libc++ today use deterministic seeds
within a process and don't randomise across processes of the same binary, so this
works in practice.  But a future libstdc++ that picks up ASLR-driven seeding (rumoured
several times) or a different STL (libstdc++ vs libc++ across heterogeneous ranks —
unlikely in MPI jobs but possible in container setups) would break this silently.

The `MFEM_VERIFY` would fire correctly *if* the hashes ever diverge, so the failure
mode is loud abort, not silent corruption — that's why this is LOW.

**Suggested fix:**
Use a deterministic hash that's not implementation-defined:

```diff
--- a/miniapps/seas/io/fault_vtu_binary.hpp
+++ b/miniapps/seas/io/fault_vtu_binary.hpp
@@ -432,7 +432,13 @@
       std::string concat;
       for (const auto &n : local.field_names) { concat += n; concat += '|'; }
-      const uint64_t local_hash = std::hash<std::string>{}(concat);
+      // R-207: don't rely on std::hash<std::string> being deterministic
+      // across processes (implementation-defined).  Use FNV-1a 64-bit —
+      // 5 lines, deterministic, no dependencies.
+      uint64_t local_hash = 14695981039346656037ULL;
+      for (unsigned char c : concat)
+      {
+         local_hash ^= c;
+         local_hash *= 1099511628211ULL;
+      }
       uint64_t root_hash = local_hash;
       MPI_Bcast(&root_hash, 1, MPI_UNSIGNED_LONG_LONG, 0, comm);
```

---

## Summary

- **Critical issues:** 0 (no immediate correctness bugs that block production)
- **Moderate issues:** 3
  - R-201 fault-only `CommitSchedule` fragility under hysteresis
  - R-202 face_trace_logger Antiplane-OOB defeating R-003 fix
  - R-203 SAFS guard too narrow on t1-fallback sign-flip
- **Low issues:** 4
  - R-204 PVD O(N²) rewrites + non-atomic truncate
  - R-205 `ComputePerDOFCoordsAndBasis_` early-return doesn't reset counters
  - R-206 `tau_pre_` caching brittle vs `ComputeSAFSParams` order
  - R-207 `std::hash<std::string>` determinism assumption
- **Plan compliance:** PARTIAL
  - Merge correctness: FULL — no conflicts left, MFEM rebuilt with HDF5 + PETSc, fix
    commit `1a805ef` cleanly applied to all 9 round-1 findings.
  - Paraview capability: FULL — VTKHDF + ZFP + cap-aware schedule all live;
    `test_fault_surface_vtkhdf*`, `test_paraview_schedule_cap`, and
    `test_R003_antiplane_diagnostics_sizing` cover the new surface.
  - "No effect on SAFS development": **PARTIAL** — R-201 is fault-PV-only;
    R-203 is the SAFS-specific gap (t1-fallback) the round-1 R-001 fix left open;
    R-206 is SAFS construction-order brittleness.  These should be addressed before
    SAFS production runs go live.
  - No critical bugs: CONFIRMED for the immediate production path.  R-202 is
    technically a latent crash but the composition is unreached today.
- **Verdict:** PASS WITH FIXES — R-201 and R-203 should be applied before the next
  SAFS-mode sweep; R-202 should be applied before any test or driver couples
  Antiplane with `FaceTraceLogger`; lows are housekeeping.

## Unreviewed Areas

- **`miniapps/seas/io/data_field_3d.{hpp,cpp}` and `stress_field_3d.{hpp,cpp}`** — the
  sidecar readers that `ApplySAFSMode` calls into.  Re-read only at the interface level
  (signature shape, error paths from constructor); their internal HDF5 dataset reads
  were exercised by `test_data_field_3d`, `test_stress_field_3d`,
  `test_field_projector`, `test_project_fault_prestress` (all passed on the prior
  sweep) and not line-by-line audited this round.
- **`miniapps/seas/safs/` content** — fault-geometry meshes + CFM data files; no
  source code, only artifacts.
- **Round-1 test sweep `b1d5ehu8w`** — assumed to have completed clean since the fix
  commit landed; not re-run as part of this round-2 audit.
- **Drivers `tpv102_driver.cpp`, `tpv205_driver.cpp`** — the `paraview_write`
  lambda is structurally identical to `tpv104_driver.cpp`; R-201 diff applies to all
  three.  Cross-driver verification deferred to /code-fix.
- **`miniapps/seas/document/io_dev/PLAN_paraview_compaction_2026-04-28.md`** — the
  1942-line plan was scope-reviewed (chapter headings + Phase 6 sections only) for
  cross-check against implementation; line-by-line plan-vs-code audit not redone.
- **Performance impact of R-204 (~1 GB cumulative PVD writes)** — flagged but not
  measured.  Quantify on a 5000-snapshot BP5 run; pull `du -sb output_dir/FaultSurface/`
  before/after to confirm magnitude.
