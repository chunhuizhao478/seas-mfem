# Code Review: SAFS v3_4_1 THERMAL CASE2 port — Phases 2/3 (2026-07-09)

## Review Scope

- **Plan:** `miniapps/seas/safs/project_7.0_preferred/document/PLAN_thermal_case2_mixedflux_port_2026-07-08.md`
- **Files reviewed:**
  - `miniapps/seas/spatial/code/spatial_friction.hpp` / `.cpp`
  - `miniapps/seas/drivers/spatial_dyn_driver.cpp`
  - `miniapps/seas/Makefile`
  - `miniapps/seas/tests/unit/test_spatial_friction_sidecar.cpp` (new)
  - `miniapps/seas/tests/unit/test_spatial_friction_config.cpp`
  - `miniapps/seas/tests/unit/test_slip_law_srw_psi.cpp`
  - `miniapps/seas/safs/project_7.0_preferred/friction/code/build_pref_friction_sidecar.py` (new)
  - `miniapps/seas/safs/project_7.0_preferred/velocity/code/build_pref_velocity_sidecar.py`
  - `miniapps/seas/safs/project_7.0_preferred/config/*.toml` (2 new)
  - `miniapps/seas/jobs/safs/safs_expanse/spatial_dyn_ratestate_v3_4_1_pref_thermal_case2_*.sbatch` (2 new)
- **Domain context:** `miniapps/seas/CLAUDE.md` (sign conventions, "files requiring extreme care",
  byte-exact regression contract), repo-root `CLAUDE.md` (report errors, do not silently simplify),
  `safs/project_7.0_alternative/document/spatial_friction_config_schema.md`,
  `io/data_field_3d.hpp` (schema-v1 / `OOBPolicy` contract).
- **Known deviation declared by the implementer:** D-A — `RateStateSidecarFields` members are
  `std::function` rather than the plan's `std::unique_ptr<DataField3D>`, to keep HDF5 out of ~15 link
  targets. **Reviewed and accepted:** the justification is verifiable (`seas_test_spatial_friction_config`
  links `SPATIAL_FRICTION_OBJ` and neither `DATA_FIELD_3D_OBJ` nor `HDF5_LIBS`), every plan semantic is
  preserved (including the `if (fields->b)` bool contract), and the idiom matches two existing codebase
  precedents. Not a finding.

---

## Findings

### [R-001] CRITICAL [spatial_dyn_driver.cpp:PrintFrictionSidecarSummary] — Acceptance gate reports clamped edge values as if they were in-hull samples

**Category:** BUG

**Description:**
The `--print-derived` friction summary samples the sidecar at the nucleation centre and prints
`a` and `V_w` there. This is the *online counterpart of the offline G4 hypocentre-anchor gate* — the whole
point is to confirm that MFEM reads `a = 0.015`, `V_w = 0.05` at the hypocentre.

But the evaluators are built with `OOBPolicy::Clamp` whenever `far_field_clamp = true` (the default, and
what the production config sets). Under `Clamp`, a query **outside the sidecar hull** silently returns the
nearest-edge value instead of aborting. `PrintFrictionSidecarSummary` has no way to tell the two apart and
prints the clamped value with no annotation.

Because the CTM CASE2 field's global minima happen to be exactly the anchors it is checking
(`min a = 0.015`, `min V_w = 0.05`), an out-of-hull point clamps to **precisely the values the gate
expects**. The gate cannot fail. It is not a check; it is a constant.

**Trigger:**
Any run where the nucleation centre (or `[hypocenter]`) lies outside the friction sidecar's bounding box,
with `far_field_clamp = true`.

Observed for real during implementation, at np=4:

```
friction hull : x[303000,696000] y[3612000,3901500] z[-21000,200]
tpv104 nucleation centre (0, 0, -7500)  -> 300 km OUTSIDE the hull

[print-derived] friction sidecar (a, V_w)
  sidecar at [nucleation.gradual_overstress_compact_circular] centre (0, 0, -7500):
      a   = 0.015        <- the exact CASE2 anchor
      V_w = 0.05 m/s     <- the exact CASE2 anchor
```

Both printed numbers are the sidecar's global minima, reached by clamping. A reader — or an automated gate
scraping this output — concludes the hypocentre anchor passed.

**Actual behavior:**
Prints clamped edge values indistinguishably from genuine in-hull samples. The Phase-3 acceptance gate is
a false positive for any out-of-hull sample point.

Secondary defect from the same root cause: with `far_field_clamp = false` (`OOBPolicy::Abort`) and an
out-of-hull sample point, `fields.a(hx,hy,hz)` calls `MFEM_ABORT` and **kills the production run from
inside a diagnostic**. The existing `have_point` guard only covers the `(0,0,0)` default; it does not
cover a genuinely-specified but out-of-hull centre.

**Expected behavior:**
The summary must know the sidecar's bounding box, state explicitly whether the sample point is inside it,
and never abort the run from a diagnostic. An out-of-hull sample must be labelled as such (or skipped), so
it can never be mistaken for a passing anchor check.

**Suggested fix:**
Carry the hull and the policy alongside the evaluators. `DataField3D::BBox()` already exists and costs
nothing to copy.

```diff
--- a/miniapps/seas/spatial/code/spatial_friction.hpp
+++ b/miniapps/seas/spatial/code/spatial_friction.hpp
 struct RateStateSidecarFields
 {
    using FieldFn = std::function<real_t(real_t, real_t, real_t)>;
 
    FieldFn a;     ///< REQUIRED when the struct is supplied
    FieldFn V_w;   ///< REQUIRED when the struct is supplied
    FieldFn b;     ///< optional; empty ⇒ keep the already-resolved b
    FieldFn Dc;    ///< optional; empty ⇒ keep the already-resolved Dc
+
+   /// Inclusive data hull {xmin,xmax,ymin,ymax,zmin,zmax} of the underlying
+   /// sidecar, copied from DataField3D::BBox() at load time.  Diagnostics MUST
+   /// consult this before sampling: under OOBPolicy::Clamp an out-of-hull query
+   /// silently returns an edge value, and for the CTM CASE2 field those edge
+   /// values are exactly the hypocentre anchors a gate is trying to verify.
+   std::array<real_t, 6> bbox = {{0, 0, 0, 0, 0, 0}};
+   /// true ⇒ the evaluators clamp out-of-hull queries; false ⇒ they abort.
+   bool clamps_out_of_hull = false;
+
+   /// True iff (x,y,z) lies inside `bbox` (inclusive).
+   bool InHull(real_t x, real_t y, real_t z) const
+   {
+      return x >= bbox[0] && x <= bbox[1]
+          && y >= bbox[2] && y <= bbox[3]
+          && z >= bbox[4] && z <= bbox[5];
+   }
 };
```

```diff
--- a/miniapps/seas/drivers/spatial_dyn_driver.cpp
+++ b/miniapps/seas/drivers/spatial_dyn_driver.cpp
    auto open = [&spec, oob](const std::string &field)
    {
       auto reader = std::make_shared<seas::DataField3D>(spec.path, field, oob);
       return RateStateSidecarFields::FieldFn(
                 [reader](real_t x, real_t y, real_t z)
                 { return reader->Evaluate(x, y, z); });
    };
 
    auto fields = std::make_shared<RateStateSidecarFields>();
-   fields->a   = open(spec.a_field);
+   // Load `a` through a named reader so its bbox can be recorded; all fields in
+   // a schema-v1 sidecar share one grid (StressField3D asserts this), so one
+   // bbox describes them all.
+   {
+      auto a_reader = std::make_shared<seas::DataField3D>(spec.path, spec.a_field, oob);
+      fields->bbox = a_reader->BBox();
+      fields->clamps_out_of_hull = spec.far_field_clamp;
+      fields->a = RateStateSidecarFields::FieldFn(
+                     [a_reader](real_t x, real_t y, real_t z)
+                     { return a_reader->Evaluate(x, y, z); });
+   }
    fields->V_w = open(spec.V_w_field);
```

```diff
-   std::cout << "  sidecar at " << origin << " (" << hx << ", " << hy << ", "
-             << hz << "):\n"
-             << "      a   = " << fields.a(hx, hy, hz) << "\n"
-             << "      V_w = " << fields.V_w(hx, hy, hz) << " m/s\n";
+   if (!fields.InHull(hx, hy, hz))
+   {
+      // NEVER evaluate here: under Clamp we would print an edge value that is
+      // indistinguishable from a real sample (and, for the CTM CASE2 field,
+      // equals the very anchor a gate checks); under Abort we would kill the
+      // run from inside a diagnostic.
+      std::cout << "  sidecar at " << origin << " (" << hx << ", " << hy << ", "
+                << hz << "): OUTSIDE the sidecar hull "
+                << "x[" << fields.bbox[0] << "," << fields.bbox[1] << "] "
+                << "y[" << fields.bbox[2] << "," << fields.bbox[3] << "] "
+                << "z[" << fields.bbox[4] << "," << fields.bbox[5] << "]\n"
+                << "      NOT SAMPLED — an out-of-hull read would "
+                << (fields.clamps_out_of_hull ? "silently clamp to an edge value."
+                                              : "abort the run.")
+                << "\n" << std::flush;
+      return;
+   }
+
+   std::cout << "  sidecar at " << origin << " (" << hx << ", " << hy << ", "
+             << hz << ")  [in hull]:\n"
+             << "      a   = " << fields.a(hx, hy, hz) << "\n"
+             << "      V_w = " << fields.V_w(hx, hy, hz) << " m/s\n";
```

Add `#include <array>` is already present in `spatial_friction.hpp`.

**Test case:**
```cpp
// tests/unit/test_spatial_friction_sidecar.cpp
static void S_8_bbox_and_in_hull_predicate()
{
   // Fixture hull is x[0,2] y[-1,1] z[-2,0].
   const std::string path = make_friction_sidecar();
   auto f = make_fields(path, OOBPolicy::Clamp);

   TEST_ASSERT(f->InHull(1.0, 0.0, -1.0),  "centre of the hull is in-hull");
   TEST_ASSERT(f->InHull(0.0, -1.0, -2.0), "min corner is in-hull (inclusive)");
   TEST_ASSERT(f->InHull(2.0,  1.0,  0.0), "max corner is in-hull (inclusive)");
   TEST_ASSERT(!f->InHull(10.0, 10.0, 10.0), "far point is out of hull");
   TEST_ASSERT(!f->InHull(2.001, 0.0, -1.0), "just past x_max is out of hull");

   // The regression this guards: a clamped out-of-hull read returns the edge
   // value, which for this fixture is NOT distinguishable from a real sample.
   const real_t clamped = f->a(10.0, 10.0, 10.0);
   const real_t at_corner = f->a(kX.back(), kY.back(), kZ.back());
   TEST_ASSERT(clamped == at_corner,
               "clamped out-of-hull read == edge value (why InHull must gate it)");

   TEST_ASSERT(f->clamps_out_of_hull, "clamp flag recorded from the spec");
   auto g = make_fields(path, OOBPolicy::Abort);
   // NOTE: make_fields must thread far_field_clamp through; see fix.
   ::unlink(path.c_str());
}
```

---

### [R-002] MODERATE [spatial_friction.hpp:FrictionSidecarSpec] — `far_field_clamp` defaults to `true`, silently hiding a sidecar that does not cover the fault

**Category:** ASSUMPTION

**Description:**
`FrictionSidecarSpec::far_field_clamp` defaults to `true` (⇒ `OOBPolicy::Clamp`). Two things make this the
wrong default:

1. **It disagrees with its own sibling.** `VelocitySpec::far_field_clamp` defaults to **`false`**
   (`spatial_friction.hpp`), and its doc-comment says the strict default "keeps the strict
   interpolation-only containment contract". `DataField3D`'s own default is `OOBPolicy::Abort`. The friction
   block is the only schema-v1 consumer that opts into silent clamping by default.
2. **Friction is evaluated only at fault DOFs.** Unlike the velocity sidecar (queried across the whole mesh,
   including the far-field absorbing box that genuinely pokes outside the CVM hull), the friction resolver
   only samples the fault. A correct friction sidecar therefore never needs clamping. A sidecar that *does*
   trigger clamping is, by construction, one that fails to cover the fault — precisely the condition that
   should abort loudly.

With the default as-is, shipping a friction sidecar whose grid misses part of the fault produces a run that
silently edge-holds `a` and `V_w` on the uncovered facets. That is silent wrong physics.

**Trigger:** any config that omits `far_field_clamp` and supplies a sidecar that does not cover the fault.

**Actual behavior:** uncovered fault DOFs receive nearest-edge `a` / `V_w`; the run proceeds.

**Expected behavior:** abort at setup, naming the field and the offending coordinate (which
`DataField3D::Evaluate`'s abort message already does).

**Suggested fix:**
```diff
--- a/miniapps/seas/spatial/code/spatial_friction.hpp
+++ b/miniapps/seas/spatial/code/spatial_friction.hpp
-   /// true ⇒ `OOBPolicy::Clamp` (ASAGI nearest-edge hold) for queries outside
-   /// the sidecar hull; false ⇒ `OOBPolicy::Abort`.  The CTM hull covers the
-   /// SAFS fault but not the far-field absorbing box, so SAFS needs `true`.
-   bool        far_field_clamp = true;
+   /// true ⇒ `OOBPolicy::Clamp` (ASAGI nearest-edge hold) for queries outside
+   /// the sidecar hull; false ⇒ `OOBPolicy::Abort`.
+   ///
+   /// Defaults to FALSE, matching `VelocitySpec::far_field_clamp` and
+   /// `DataField3D`'s own `OOBPolicy::Abort`.  Friction is resolved ONLY at
+   /// fault DOFs, so a sidecar that covers the fault never needs clamping; one
+   /// that triggers clamping does not cover the fault, and that must abort
+   /// rather than silently edge-hold `a` / `V_w`.  SAFS sets it explicitly.
+   bool        far_field_clamp = false;
```

Both production TOMLs already set `far_field_clamp = true` explicitly, so this changes no shipped run.
Update the schema doc's default column accordingly.

**Test case:**
```cpp
// tests/unit/test_spatial_friction_config.cpp
static void T_76_sidecar_far_field_clamp_defaults_false()
{
   const std::string toml = MinimalLSWHeader(1, "rate_state") + MinimalRSBlock()
      + "[friction.rate_state.sidecar]\npath = \"f.h5\"\n";
   const auto& sc = ParseSpatialFrictionConfigString(toml).rate_state->sidecar;
   TEST_ASSERT(!sc.far_field_clamp,
               "far_field_clamp defaults to false (strict, like [velocity])");
}
```

---

### [R-003] MODERATE [config/*.toml:`[material]`] — `sidecar_path` is never opened; the comment says it resolves the material

**Category:** ASSUMPTION

**Description:**
Both production TOMLs carry:

```toml
[material]
kind         = "sidecar_hdf5"             # required: interior_flux="matrix" forbids kind="constant"
sidecar_path = ".../multiscale_statewise_cvm_v3_4_1/velocity_safs.h5"
```

with the comment *"Resolves to the same h5 as `[velocity].override_path`."* That is false.
`cfg.material.sidecar_path` appears in exactly two places in the whole tree — it is parsed
(`spatial_friction.cpp:1736`) and checked non-empty (`:1793`). **Nothing ever opens it.** The heterogeneous
material comes solely from the `[velocity]` bundle; `[material].kind` is a gate that only has to be
non-`constant` for `interior_flux = "matrix"`.

The comment invites a future maintainer to repoint `[material].sidecar_path` at a different CVM and expect
the run to change. It will not. Silent no-op config.

**Trigger:** editing `[material].sidecar_path` and expecting it to take effect.

**Actual behavior:** the value is validated non-empty and discarded.

**Expected behavior:** the comment states plainly that the key is an inert gate, and points at
`[velocity].override_path` as the load-bearing one.

**Suggested fix:**
```diff
 [material]
-# Required: interior_flux="matrix" forbids [material].kind="constant".
-# Resolves to the same h5 as [velocity].override_path.
+# `kind` is a GATE ONLY: interior_flux="matrix" forbids [material].kind="constant".
+# `sidecar_path` is parsed and checked non-empty, then DISCARDED — nothing opens it
+# (grep cfg.material.sidecar_path: parsed at spatial_friction.cpp:1736, verified at
+# :1793, never read).  The heterogeneous material is loaded ENTIRELY from the
+# [velocity] block above.  Change [velocity].override_path to change the material;
+# changing the path below has NO effect.
 kind         = "sidecar_hdf5"
 sidecar_path = "safs/project_7.0_preferred/velocity/results/multiscale_statewise_cvm_v3_4_1/velocity_safs.h5"
```

**Test case:**
```python
def test_R003_material_sidecar_path_is_inert():
    # Grep-level assertion: the only uses of cfg.material.sidecar_path are the
    # parse and the non-empty check.  If a third use ever appears, this test
    # fails and the TOML comment must be corrected.
    import subprocess
    hits = subprocess.run(
        ["grep", "-rn", "material.sidecar_path", "miniapps/seas/"],
        capture_output=True, text=True).stdout.strip().splitlines()
    code_hits = [h for h in hits if h.endswith(".cpp") or ".cpp:" in h]
    assert len(code_hits) == 2, f"material.sidecar_path gained a consumer: {code_hits}"
```

---

### [R-004] MODERATE [build_pref_friction_sidecar.py] — `sidecar` import can silently bind to the wrong module

**Category:** EDGE_CASE

**Description:**
```python
REPO = Path(__file__).resolve().parents[4]          # .../miniapps/seas
ALT = REPO / "safs" / "project_7.0_alternative"
sys.path.insert(0, str(ALT / "velocity" / "code"))
from sidecar import write_sidecar
```
`parents[4]` is a positional assumption about directory depth. If the script is moved one level (or copied
into a scratch dir, which is exactly how these builders get used), `ALT` silently points somewhere that does
not exist, `sys.path.insert` succeeds anyway, and `from sidecar import write_sidecar` then either raises a
bare `ModuleNotFoundError` with no hint about the real cause, or — worse — imports some *other* `sidecar.py`
that happens to be importable, producing a sidecar with the wrong schema.

The velocity builder has the identical construct.

**Trigger:** running the script from a copied/moved location, or any future directory reshuffle.

**Actual behavior:** cryptic `ModuleNotFoundError`, or a wrong-module import.

**Expected behavior:** fail immediately with a message naming the path it expected.

**Suggested fix:**
```diff
 REPO = Path(__file__).resolve().parents[4]          # .../miniapps/seas
 ALT = REPO / "safs" / "project_7.0_alternative"
-sys.path.insert(0, str(ALT / "velocity" / "code"))
-from sidecar import write_sidecar  # noqa: E402
+_SIDECAR_DIR = ALT / "velocity" / "code"
+if not (_SIDECAR_DIR / "sidecar.py").is_file():
+    raise SystemExit(
+        "build_pref_friction_sidecar: cannot find the shared schema-v1 writer at\n"
+        f"    {_SIDECAR_DIR / 'sidecar.py'}\n"
+        "This script locates it relative to its own path (parents[4] == miniapps/seas).\n"
+        "Run it from its home in the repo, or fix REPO above."
+    )
+sys.path.insert(0, str(_SIDECAR_DIR))
+from sidecar import write_sidecar  # noqa: E402
```
Apply the same guard to `build_pref_velocity_sidecar.py`.

**Test case:**
```python
def test_R004_import_guard_fires_when_moved(tmp_path):
    import shutil, subprocess, sys
    src = "miniapps/seas/safs/project_7.0_preferred/friction/code/build_pref_friction_sidecar.py"
    dst = tmp_path / "build_pref_friction_sidecar.py"
    shutil.copy(src, dst)
    r = subprocess.run([sys.executable, str(dst), "--in", "x.nc", "--out", "y.h5"],
                       capture_output=True, text=True)
    assert "cannot find the shared schema-v1 writer" in (r.stderr + r.stdout), \
        "moved script must fail with a path-naming message, not ModuleNotFoundError"
```

---

### [R-005] LOW [spatial_friction.cpp:parse_rate_state] — `MFEM_ASSERT` belt-and-braces compiles away in release

**Category:** QUALITY

**Description:**
```cpp
MFEM_ASSERT(!(out.sidecar.enabled && out.depth_profile.enabled),
            "[friction.rate_state] sidecar + depth_profile survived the "
            "table-presence guard");
```
`MFEM_ASSERT` is a no-op when `MFEM_DEBUG` is off, which is the production build. The comment presents it as
a safety net ("asserts the two stay in sync"), but in the build that matters it checks nothing. It is also
now fully redundant with the `MFEM_VERIFY` added to `resolve_rs_impl`, which *is* always live and covers the
hand-built-config case the assert was aimed at.

**Suggested fix:** delete it, and let the resolver's `MFEM_VERIFY` be the single real guard.
```diff
-   // Belt-and-braces: the table-presence guard above already rejected the
-   // ambiguous pair before any file was opened.  This asserts the two stay in
-   // sync if either parser ever learns to set `enabled` some other way.
-   MFEM_ASSERT(!(out.sidecar.enabled && out.depth_profile.enabled),
-               "[friction.rate_state] sidecar + depth_profile survived the "
-               "table-presence guard");
+   // (The table-presence guard above rejected the ambiguous pair before any
+   // file was opened.  ResolveRateState re-checks it with a live MFEM_VERIFY
+   // for callers that build a RateStateBlock by hand, bypassing the parser.)
```

---

### [R-006] LOW [build_pref_friction_sidecar.py:trilinear] — docstring claims it mirrors `DataField3D`, but it clamps where `DataField3D` aborts

**Category:** QUALITY

**Description:**
The comment says *"mirrors io/data_field_3d.cpp so the G4 gate below tests what MFEM will actually read"*.
`locate()` clamps at both ends (`p <= ax[0] -> (0, 0.0)`, `p >= ax[-1] -> (len-2, 1.0)`), i.e. it reproduces
`OOBPolicy::Clamp`, not the schema-v1 default `OOBPolicy::Abort`. It is harmless today because
`check_hypocentre` verifies containment first, but the comment is wrong and the ordering dependency is
undocumented.

**Suggested fix:** state the clamp explicitly and the precondition.
```diff
-# trilinear sampling — mirrors io/data_field_3d.cpp so the G4 gate below tests
-# what MFEM will actually read, not what numpy thinks the nearest node is.
+# trilinear sampling — mirrors io/data_field_3d.cpp's INTERPOLATION so the G4
+# gate tests what MFEM will read, not what numpy thinks the nearest node is.
+# NOTE: `locate()` CLAMPS at both ends (i.e. OOBPolicy::Clamp), whereas
+# DataField3D defaults to OOBPolicy::Abort.  Callers MUST verify containment
+# first — `check_hypocentre` does, before it calls this.
```

---

### [R-007] LOW [spatial_dyn_driver.cpp:PrintFrictionSidecarSummary] — `switch (nuc.kind)` has no `default:`

**Category:** EDGE_CASE

**Description:**
The switch covers all three current `NucleationKind` values, and `have_point = true` is set *before* it. If a
fourth kind is added, the code compiles (with a `-Wswitch` warning at most), `have_point` stays `true`, and
`hx/hy/hz` silently keep the `[hypocenter]` values — which may be the `(0,0,0)` default — while `origin`
still reads `"[hypocenter]"`. The result is a diagnostic sampling the origin of the UTM grid.

With R-001's `InHull` gate in place the consequence degrades to a confusing "OUTSIDE the hull" line rather
than a wrong anchor, so this is LOW once R-001 lands.

**Suggested fix:**
```diff
          case seas::spatial::NucleationKind::InstantaneousOverstressCircular:
             ...
             break;
+         default:
+            MFEM_ABORT("PrintFrictionSidecarSummary: unhandled NucleationKind "
+                       << static_cast<int>(nuc.kind)
+                       << "; add its centre to the switch.");
       }
```

---

### [R-008] LOW [sbatch: per-rank memory pre-flight] — `du -m` is a proxy that breaks if sidecars are ever compressed

**Category:** ASSUMPTION

**Description:**
The guard estimates resident sidecar memory from on-disk size:
```bash
_sidecar_mb=$(du -m "${VEL_H5}" "${STR_H5}" "${FRI_H5}" | awk '{s+=$1} END {print s}')
```
This is exact **only because** `sidecar.py::write_sidecar` and `csm_stress_nc_to_mfem_hdf5.py` write
uncompressed `float64`. If anyone adds `compression="gzip"` to the writers (an obvious future optimisation
for a 1.3 GB scp), disk size collapses while `DataField3D` still materialises the full `float64` array, and
the guard would wave through a job that OOMs.

**Suggested fix:** state the invariant where it can be seen, and make the failure mode explicit.
```diff
 # v3_4_1 raises this from ~780 MB/rank ... to ~1330 MB/rank
+#
+# NOTE: this uses on-disk size as a proxy for resident size.  That is EXACT only
+# while the sidecar writers emit UNCOMPRESSED float64 (they do: sidecar.py and
+# csm_stress_nc_to_mfem_hdf5.py pass no `compression=`).  If a writer ever gains
+# gzip/szip, disk size collapses but DataField3D still materialises the full
+# float64 array — this guard would then under-count.  Compute from the dataset
+# shapes instead if that day comes.
 _sidecar_mb=$(du -m "${VEL_H5}" "${STR_H5}" "${FRI_H5}" | awk '{s+=$1} END {print s}')
```

---

## Summary

- Critical issues: **1** (R-001)
- Moderate issues: **3** (R-002, R-003, R-004)
- Low issues: **4** (R-005, R-006, R-007, R-008)
- Plan compliance: **FULL** — every Phase 2/3/4/5/6/8 requirement is implemented. Phase 7 is deferred by
  decision D6. Deviation D-A is declared, justified and verifiable. Two additions beyond the plan
  (`SPATIAL_HEADERS` in the Makefile; the sbatch memory pre-flight) are direct consequences of defects the
  implementation itself exposed (plan G13 and G12) and are in scope.
- Verdict: **PASS WITH FIXES** — R-001 must be fixed before the Phase-3 acceptance gate can be trusted. It
  does not affect the *physics* of a production run (the resolver's per-DOF path is correct and independently
  verified: `dt_cfl` matched prediction to 3 s.f., `σ_n` at the hypocentre facet reproduced the deck's
  169.83 MPa, and 24 test binaries pass). It affects whether the gate that certifies the run means anything.

---

## Resolution (fix pass, 2026-07-09)

| ID | Severity | Status | Where |
|----|----------|--------|-------|
| R-001 | CRITICAL | **FIXED** | `spatial_friction.hpp` (`bbox`, `clamps_out_of_hull`, `InHull()`), `spatial_dyn_driver.cpp` (record bbox at load; gate the point sample) |
| R-002 | MODERATE | **FIXED** | `spatial_friction.hpp` + `.cpp` — `far_field_clamp` defaults to `false` |
| R-003 | MODERATE | **FIXED** | both production TOMLs — `[material]` comment now states the key is inert |
| R-004 | MODERATE | **FIXED** | both Python builders — path-naming import guard |
| R-005 | LOW | **FIXED** | `spatial_friction.cpp` — dead `MFEM_ASSERT` removed |
| R-006 | LOW | **FIXED** | `build_pref_friction_sidecar.py` — `trilinear` docstring corrected |
| R-007 | LOW | **FIXED (variant)** | `spatial_dyn_driver.cpp` — see note below |
| R-008 | LOW | **FIXED** | both sbatch — `du -m` invariant documented |

**R-007 deviation from the suggested fix.** The review proposed adding
`default: MFEM_ABORT(...)` to the `switch (nuc.kind)`. That would *silence* `-Wswitch`, trading a
compile-time error for a runtime one — strictly worse. Instead, `have_point` is now cleared before the
switch and set to `true` inside each `case`. The switch stays enum-complete (so `-Wswitch` still flags a
newly-added kind at compile time) **and** a kind that somehow slips through degrades to the
"point sample skipped" branch rather than silently sampling `(0,0,0)`.

**Verification of R-001 on its exact reproducer** (np=4, nucleation centre 300 km outside the hull):

```
before:  sidecar at [...] centre (0, 0, -7500):
             a   = 0.015          <- clamped edge value, read as a PASS
             V_w = 0.05 m/s

after :  sidecar at [...] centre (0, 0, -7500): OUTSIDE the sidecar hull
             x[303000, 696000] y[3.612e+06, 3.9015e+06] z[-21000, 200]
             NOT SAMPLED — an out-of-hull read would silently clamp to an edge value (a false PASS).
```

and the positive branch, with the centre moved inside the hull:

```
after :  sidecar at [...] centre (609063, 3.70953e+06, -10000)  [in hull]:
             a   = 0.015
             V_w = 0.05 m/s
```

The production hypocentre `(609062.8722, 3709528.1324, -10000)` is inside the hull
`x[303000,696000] y[3612000,3901500] z[-21000,200]`, so the production gate samples rather than skips —
and now genuinely *can* fail.

**Test suite after fixes:** 24/24 affected binaries pass from a forced-rebuild of the header dependency
chain. `test_spatial_friction_sidecar` 40/40 (was 25, +S-8), `test_spatial_friction_config` 314/314
(was 312, +T-76), `test_spatial_friction_resolver` 127/127 (unchanged — byte-exact legacy path),
`test_slip_law_srw_psi` 27/27.

---

## Unreviewed Areas

- `miniapps/seas/tests/unit/test_slip_law_srw_psi.cpp` additions were checked for correctness of the
  `f_w = 0` claim (they byte-match an independent inlined reference already in that file) but not
  re-derived from the SCEC TPV104 spec.
- The two `.sbatch` scripts were syntax-checked (`bash -n`) and their arm-detection and memory arithmetic
  exercised against both real configs, but **never submitted**. Scheduler-side behaviour
  (`SLURM_MEM_PER_NODE` / `SLURM_NTASKS_PER_NODE` population under `--mem=200000M` on Expanse) is
  unverified.
- Runtime behaviour beyond `--dry-run` setup is unverified: no time-stepping was executed on the production
  mesh. The nucleation-amplitude finding (75 MPa → predicted no nucleation) is the driver's own diagnostic,
  not a review finding, and is recorded in the plan's Gate 4 as a decision for the user.
