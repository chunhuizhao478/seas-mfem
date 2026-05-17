# Code Review: All implemented phases (1–5 Python + Phase 6 Tranche 1 C++)

## Review Scope
- Plan: `PLAN_onfaultstress.md` (Phases 1–5: lines 318–1536; Phase 6 §1, §2, §3: lines 1587–1702)
- Phase 1–4 source (Python, prior cycles): `code_preprocess/project_to_fault_stress.py`, `test_project_to_fault_stress.py`
- Phase 5 source (Python, prior cycle): `code_preprocess/build_stress_safs.py`, `test_build_stress_safs.py`
- Phase 6 Tranche 1 source (C++, this cycle): `io/stress_field_3d.{hpp,cpp}`,
  `io/stress_field_coefficient.hpp`, `io/field_coefficient.{hpp,cpp}` (modified),
  `tests/unit/test_stress_field_3d.cpp`, `Makefile` (modified),
  `PHASE6_DEVIATIONS.md`
- Domain context: `miniapps/seas/CLAUDE.md`, R-501/R-502 sign-flip rule,
  schema-v1 canonical component order (xx, yy, zz, xy, yz, xz), CLAUDE.md
  "Files Requiring Extreme Care" list
- Live test runs verified before review:
  - Python: **170 / 170 passed** (`test_project_to_fault_stress.py` + `test_build_stress_safs.py`)
  - C++ Phase 6 Tranche 1: **46 / 46 passed** (`seas_test_stress_field_3d`)
- Prior review status (Phase 1–4 round R-701..R-706 and Phase 5
  round R-801..R-806): **all resolved**, regression-guard tests added.

This pass is a fresh adversarial sweep across all implemented phases
with primary focus on the new Phase 6 Tranche 1 deliverables.

---

## Findings

### [R-901] MODERATE [tests/unit/test_stress_field_3d.cpp + Phase 6 §2] — `StressFieldCoefficient` has zero unit-test coverage

**Category:** DEVIATION / BUG (test gap on plan-required surface)

**Description:**
The plan calls out `StressFieldCoefficient` as a deliverable in
Phase 6 §2 (lines 1647–1673) with explicit contract requirements:

- vdim = 6 (set in the base-class ctor)
- `Eval(v, T, ip)` returns components in schema-v1 canonical order
  (xx, yy, zz, xy, yz, xz) — **NOT** Voigt (xx, yy, zz, yz, xz, xy)
- `scale` and `offset` applied uniformly to all six components

`test_stress_field_3d.cpp` covers `StressField3D` extensively (10
test groups, 46 assertions) but **does not include a single test
for `StressFieldCoefficient`**. The Eval path's component ordering
is one of the highest-value invariants in the entire phase — a
silent swap of (e.g.) `sigma_yz` with `sigma_xy` in the Eval body
would compile and pass every existing test, then poison the
downstream static-equilibration pipeline.

The implementer's completion report notes the deviation breakdown
but does NOT acknowledge missing coverage on `StressFieldCoefficient`.
PHASE6_DEVIATIONS.md likewise omits this gap.

**Trigger:**
Inspect the test catalog. `grep -c StressFieldCoefficient
tests/unit/test_stress_field_3d.cpp` → 0.

**Actual behavior:**
Phase 6 §2 ships untested.

**Expected behavior:**
At least one test that:
1. Constructs `StressFieldCoefficient` over a synthetic sidecar.
2. Builds a minimal mesh + IntegrationPoint + ElementTransformation.
3. Calls `Eval(v, T, ip)` and asserts vdim == 6 and that each of
   the six entries matches the corresponding schema-v1 component
   at the transformed coordinate.
4. Re-runs with non-trivial `scale` and `offset` and asserts they
   apply uniformly.

**Suggested fix:**
Add a `T_6_11_stress_field_coefficient_eval` test block to
`test_stress_field_3d.cpp`. A minimal sketch (the test uses a
unit-cube serial Mesh + an H1 FiniteElementSpace to obtain a
real `ElementTransformation`):

```cpp
#include "../../io/stress_field_coefficient.hpp"   // top of file

static void T_6_11_stress_field_coefficient_eval()
{
   std::cout << "\n[T-6-11] StressFieldCoefficient::Eval order + scale/offset\n";

   const auto spec = make_default_spec();   // bbox [0,2] x [0,2] x [-2,0]
   const std::string path = make_tmp_path("scfeval");
   write_stress_sidecar(path, spec);

   StressField3D field(path);
   StressFieldCoefficient coef(field);

   // A tiny serial cube that lies inside the synthetic sidecar bbox.
   // Vertices 0..7 of a unit cube in [0.5, 1.5]^2 x [-1.5, -0.5].
   mfem::Mesh mesh(1, 8, 1);
   const double verts[8][3] = {
      {0.5, 0.5, -1.5}, {1.5, 0.5, -1.5},
      {1.5, 1.5, -1.5}, {0.5, 1.5, -1.5},
      {0.5, 0.5, -0.5}, {1.5, 0.5, -0.5},
      {1.5, 1.5, -0.5}, {0.5, 1.5, -0.5},
   };
   for (int v = 0; v < 8; ++v) { mesh.AddVertex(verts[v]); }
   const int hex_idx[8] = {0,1,2,3,4,5,6,7};
   mesh.AddHex(hex_idx);
   mesh.FinalizeHexMesh(1);

   mfem::H1_FECollection fec(1, 3);
   mfem::FiniteElementSpace fes(&mesh, &fec);
   mfem::ElementTransformation& T = *mesh.GetElementTransformation(0);
   mfem::IntegrationPoint ip; ip.Set3(0.5, 0.5, 0.5);   // element center

   mfem::Vector x(3); T.Transform(ip, x);
   mfem::Vector v(6); coef.Eval(v, T, ip);

   const double xc = x[0], yc = x[1], zc = x[2];
   const double xx_ref = field.Field(0).Evaluate(xc, yc, zc);
   const double yy_ref = field.Field(1).Evaluate(xc, yc, zc);
   const double zz_ref = field.Field(2).Evaluate(xc, yc, zc);
   const double xy_ref = field.Field(3).Evaluate(xc, yc, zc);
   const double yz_ref = field.Field(4).Evaluate(xc, yc, zc);
   const double xz_ref = field.Field(5).Evaluate(xc, yc, zc);

   TEST_ASSERT(coef.GetVDim() == 6,
               "StressFieldCoefficient::GetVDim() == 6");
   TEST_NEAR(v(0), xx_ref, 1e-9, "v[0] == sigma_xx (schema-v1 canonical)");
   TEST_NEAR(v(1), yy_ref, 1e-9, "v[1] == sigma_yy");
   TEST_NEAR(v(2), zz_ref, 1e-9, "v[2] == sigma_zz");
   TEST_NEAR(v(3), xy_ref, 1e-9, "v[3] == sigma_xy (NOT Voigt!)");
   TEST_NEAR(v(4), yz_ref, 1e-9, "v[4] == sigma_yz");
   TEST_NEAR(v(5), xz_ref, 1e-9, "v[5] == sigma_xz");

   // Scale / offset apply uniformly.
   StressFieldCoefficient coef2(field, /*scale=*/2.0, /*offset=*/100.0);
   mfem::Vector v2(6); coef2.Eval(v2, T, ip);
   TEST_NEAR(v2(0), 2.0 * xx_ref + 100.0, 1e-9, "scale/offset on sigma_xx");
   TEST_NEAR(v2(3), 2.0 * xy_ref + 100.0, 1e-9, "scale/offset on sigma_xy");

   ::unlink(path.c_str());
}
```

Wire it into `main()` next to the other T_6_* calls.

**Test case:** the test above is itself the test case.

---

### [R-902] MODERATE [tests/unit/test_stress_field_3d.cpp + Phase 6 §"Acceptance Criteria"] — Phase 6 acceptance criteria #1 and #2 (interop with real Phase 5 output) not exercised

**Category:** DEVIATION / BUG (acceptance-criteria gap)

**Description:**
Plan §"Acceptance Criteria" (lines 1948–1956) requires two
integration assertions that the implementer's synthetic-sidecar
unit tests do NOT exercise:

1. (line 1949–1952) `StressField3D("stress_safs.h5")` constructs
   successfully against the **Phase-5 output** and exposes a
   `(3, 3)` symmetric tensor at a sample interior point that
   matches the H&Z `demo_safod` σ⁰ values to 1e-3 Pa.
2. (line 1953–1956) `FieldProjector::ProjectStress` returns six
   `ParGridFunction`s whose pointwise values at a sample of cell
   centroids match the sidecar values to within the trilinear
   interpolation noise (≤ 1 Pa for a constant field).

Acceptance #3–#6 cover deferred subscopes (correctly out of scope
for Tranche 1). Acceptance #1 and #2 cover **only the implemented
tranche** and could and should be tested now. Without them, the
Tranche-1 implementation is unit-test green but acceptance-blind
to a Phase-5/Phase-6 contract mismatch (e.g., the Phase 5 writer
emitting fields in a different name order than `StressField3D`
expects).

This is the same class of regression that the Phase 5 fix-cycle
R-802 surfaced (output-dir drift); only a cross-phase test
catches it.

**Trigger:**
The implementer's synthetic sidecar was crafted to match the
contract; a Phase 5 regression that violated the contract would
not be caught by the existing tests.

**Suggested fix:**
Add a Python-level integration test that:
1. Runs `build_stress_safs.build_stress_safs(...)` to produce a
   real `stress_safs.h5`.
2. Invokes a tiny C++ smoke binary (or, more pragmatically, a
   pytest that shells out to `seas_test_stress_field_3d_against_phase5`)
   that opens the produced sidecar via `StressField3D`, evaluates
   at the bbox center, and prints the (3, 3) tensor. The pytest
   parses the printout and compares to the analytic
   `bulk_stress_tensor_field(z) * 1e6` value.

Alternatively, add a pure-C++ test fixture that:
1. Invokes `python build_stress_safs.py ...` via `system()` (the
   test environment already loads `mfem-dev` + `pythonenv`).
2. Opens the produced HDF5 via `StressField3D`.
3. Compares `Evaluate(...)` at a sample point to the SAFOD H&Z
   reference values.

A minimal pytest-only sketch (assumes a small C++ helper binary
`seas_print_stress_at` is built; if not, the assertion can be
done in Python using `h5py` directly to confirm the field-name
contract, which is the largest plan-Phase-6 sensitivity):

```python
def test_R902_phase5_field_names_match_phase6_reader():
    """R-902: every field name expected by Phase 6's StressField3D
    must be present in Phase 5's stress_safs.h5 output."""
    import h5py
    import tempfile
    from pathlib import Path
    from build_stress_safs import build_stress_safs
    # Build a tiny sidecar via Phase 5.
    with tempfile.TemporaryDirectory() as td:
        mesh = _write_synth_bulk_vtu(Path(td))   # reuse the helper
        out = Path(td) / "stress_safs.h5"
        build_stress_safs(
            mesh_paths=[mesh], out_path=out,
            SHmax=113.0, Shmin=49.0, Sv=45.0, SHmax_az_deg=23.0,
        )
        # Phase 6 reader expects exactly these names; verify they
        # are present.
        with h5py.File(out, "r") as h5:
            fields = set(h5["fields"].keys())
        expected = {"sigma_xx", "sigma_yy", "sigma_zz",
                    "sigma_xy", "sigma_yz", "sigma_xz"}
        assert fields == expected, (
            f"R-902: Phase 5 field set {fields} does not match "
            f"Phase 6 StressField3D's expected set {expected}"
        )
```

This is the minimum-viable contract test that catches the
field-name and the canonical-order drift.

**Test case:** the test above is the test case.

---

### [R-903] LOW [io/stress_field_3d.hpp + Phase 6 §1 line 1644–1645] — `BBox()` returns sigma_xx's bbox; plan literal says "intersection of all six component bboxes"

**Category:** DEVIATION (plan-literal)

**Description:**
Plan §1 line 1644–1645:
> `BBox()` returns the intersection of all six component bboxes —
> they are identical by construction (same grid), but the
> intersection is the schema-conformant value.

Source `stress_field_3d.cpp:56`:
```cpp
bbox_ = sigma_xx_.BBox();
```

The source uses sigma_xx's bbox verbatim. The
`AssertConsistentGrid_()` check asserts equality across all six,
so functionally the stored value matches the intersection — but
only because the assertion holds. If a future maintainer relaxes
the assertion (e.g., to support meshes with per-component padding
in some hypothetical schema-v2), `BBox()` would silently return
sigma_xx's bbox rather than the actual intersection.

Documentation/robustness issue. Not a correctness bug under the
current schema-v1 contract.

**Suggested fix:**
Compute the actual intersection rather than relying on the
assertion:

```diff
 StressField3D::StressField3D(const std::string& sidecar_path,
                              OOBPolicy oob)
    : sigma_xx_(sidecar_path, "sigma_xx", oob)
    , sigma_yy_(sidecar_path, "sigma_yy", oob)
    , sigma_zz_(sidecar_path, "sigma_zz", oob)
    , sigma_xy_(sidecar_path, "sigma_xy", oob)
    , sigma_yz_(sidecar_path, "sigma_yz", oob)
    , sigma_xz_(sidecar_path, "sigma_xz", oob)
 {
-   // Pin the shared bbox (every component carries the identical
-   // grid by construction — Phase 5 writes one (Nx, Ny, Nz) per
-   // field). The assertion below catches a hand-modified sidecar.
-   bbox_ = sigma_xx_.BBox();
+   // Pin the bbox as the strict intersection of all six component
+   // bboxes (plan §1 line 1644-1645).  Under the schema-v1 contract
+   // the six are bit-exact identical (single /grid/* per file), so
+   // this reduces to the shared value; computing the intersection
+   // explicitly is defensive against future schema relaxation.
+   const DataField3D* comps[6] = {
+      &sigma_xx_, &sigma_yy_, &sigma_zz_,
+      &sigma_xy_, &sigma_yz_, &sigma_xz_
+   };
+   bbox_ = comps[0]->BBox();
+   for (int c = 1; c < 6; ++c)
+   {
+      const auto& bc = comps[c]->BBox();
+      // axis-aligned intersection: (max of mins, min of maxes)
+      bbox_[0] = std::max(bbox_[0], bc[0]);
+      bbox_[1] = std::min(bbox_[1], bc[1]);
+      bbox_[2] = std::max(bbox_[2], bc[2]);
+      bbox_[3] = std::min(bbox_[3], bc[3]);
+      bbox_[4] = std::max(bbox_[4], bc[4]);
+      bbox_[5] = std::min(bbox_[5], bc[5]);
+   }

    AssertConsistentGrid_();
 }
```

(`<algorithm>` needs to be included for `std::min` / `std::max`.)

**Test case:** N/A under the current schema-v1 contract (the
existing T-6-3 already verifies `BBox()` matches the
single-grid value). The fix is internally consistent under
both v1 and any future relaxed schema.

---

### [R-904] LOW [io/field_coefficient.hpp:184–193 + io/field_coefficient.cpp:366–410] — `ProjectStress` has no `scale` / `offset` parameters; `Project` / `ProjectVelocity` accept them

**Category:** QUALITY (missing flexibility, API inconsistency)

**Description:**
Existing siblings `Project(field, target_fes, scale=1.0,
offset=0.0)` and `ProjectVelocity` flow `scale` / `offset` to the
underlying `FieldCoefficient`. The new `ProjectStress` does NOT
accept them and uses defaults internally:

```cpp
Result rxx = Project(field.Field(0), target_fes);
```

This is fine for production (the sidecar already carries Pa), but
prevents the natural use case "project the sidecar Pa values into
GPa for visualisation" or applying a uniform shift to inspect
deviatoric components. A future caller would have to manually
post-process the six ParGridFunctions.

The plan §3 signature (line 1689–1692) does not mention scale /
offset, so this is plan-compliant. But callers familiar with
`ProjectVelocity` will be surprised. LOW.

**Suggested fix:**
Add `scale` and `offset` parameters with defaults that match the
existing `Project` defaults:

```diff
-   static StressFields ProjectStress(
-      const std::string&           sidecar_path,
-      mfem::ParFiniteElementSpace& target_fes,
-      InterpMode interp            = InterpMode::Trilinear);
+   static StressFields ProjectStress(
+      const std::string&           sidecar_path,
+      mfem::ParFiniteElementSpace& target_fes,
+      InterpMode interp            = InterpMode::Trilinear,
+      real_t scale                 = 1.0,
+      real_t offset                = 0.0);
```

Forward to inner `Project(...)` calls:

```diff
-   Result rxx = Project(field.Field(0), target_fes);
-   Result ryy = Project(field.Field(1), target_fes);
-   Result rzz = Project(field.Field(2), target_fes);
-   Result rxy = Project(field.Field(3), target_fes);
-   Result ryz = Project(field.Field(4), target_fes);
-   Result rxz = Project(field.Field(5), target_fes);
+   Result rxx = Project(field.Field(0), target_fes, scale, offset);
+   Result ryy = Project(field.Field(1), target_fes, scale, offset);
+   Result rzz = Project(field.Field(2), target_fes, scale, offset);
+   Result rxy = Project(field.Field(3), target_fes, scale, offset);
+   Result ryz = Project(field.Field(4), target_fes, scale, offset);
+   Result rxz = Project(field.Field(5), target_fes, scale, offset);
```

**Test case:** N/A (additive API change; no behavior change at
default values).

---

### [R-905] LOW [io/field_coefficient.cpp:366–410] — `ProjectStress` re-runs the mesh-bbox containment check six times

**Category:** QUALITY (performance, redundancy)

**Description:**
Each of the six inner `Project(field.Field(c), target_fes)`
calls performs `ComputeMeshBBoxParallel(*pmesh, ...)` followed
by `field.ContainsBBox(...)`. The parallel bbox computation costs
two MPI_Allreduce calls per invocation. Six redundant passes →
twelve unnecessary all-reduces per `ProjectStress` call. The
shared-grid invariant in `StressField3D::AssertConsistentGrid_()`
guarantees the check returns the same answer for all six
components, so five out of the six all-reduces are wasted work.

On a large parallel run this is measurable but not a correctness
issue.

**Suggested fix:**
Either (a) hoist the bbox check to the top of `ProjectStress`
and bypass the per-component checks, or (b) accept the cost as
"defensive consistency with sibling ProjectVelocity". Option (a)
diverges from the `ProjectVelocity` style; option (b) is
defensible. If (a) is chosen, factor the per-field assemble out
of `Project` and call it via a tighter loop. Recommended: defer
to a future perf-only refactor of `FieldProjector`. **Flag as
LOW.**

**Test case:** N/A (performance only).

---

### [R-906] LOW [io/stress_field_3d.hpp:11–16 docstring + cross-phase doc consistency] — Docstring claims sign flip happens in Phase 5 (`build_stress_safs.py:evaluate_stress_field_on_grid`); the flip actually happens in Phase 3

**Category:** QUALITY (documentation accuracy)

**Description:**
`stress_field_3d.hpp` lines 10–16:
> Sign convention: this reader is a strict pass-through.  The bulk
> path's single-source sign flip from H&Z continuum-mechanics
> (compression negative) to SEAS-internal (compression positive)
> happens once in the Python preprocessor (Phase 3
> `bulk_stress_tensor_field`) at the source; Phase 5 (sidecar
> writer) does the MPa->Pa unit conversion only; this C++ reader
> returns the on-disk values unchanged (R-501 / R-502 contract).

Good — this docstring is correct.

But the **plan §1 line 1604–1606** (and the implementer's
PHASE6_DEVIATIONS.md echoing the plan) reads:
> The full-tensor sign flip from geomechanics happens in the
> Python preprocessor at the schema-v1 sidecar write boundary
> (`code_preprocess/build_stress_safs.py:evaluate_stress_field_on_grid`);
> this reader is a pure pass-through.

The plan is **wrong**: `evaluate_stress_field_on_grid` (Phase 5
§4 line 1424–1428 of the plan) is a pure MPa → Pa unit conversion
with no sign flip. The actual sign flip lives in Phase 3
`bulk_stress_tensor_field` (Phase 3 §2 line 856–905 of the plan,
verified in `project_to_fault_stress.py:1182–1289`).

The source `stress_field_3d.hpp` docstring got this right
despite the plan-literal saying otherwise — which is good for
the source but means the plan and the deviation doc disagree
with the implemented header. Future maintainers reading the
plan first will be confused.

**Suggested fix:**
Update plan §1 line 1604–1606 to point to the correct phase:

```diff
       /// Evaluate the symmetric stress tensor at (x, y, z) in canonical
       /// CRS. Returns a (3, 3) dense matrix in Pa, **compression
       /// POSITIVE (SEAS internal convention)**.  The full-tensor sign
-      /// flip from geomechanics happens in the Python preprocessor at
-      /// the schema-v1 sidecar write boundary
-      /// (`code_preprocess/build_stress_safs.py:evaluate_stress_field_on_grid`);
+      /// flip from geomechanics happens in the Python preprocessor at
+      /// the source-site (Phase 3
+      /// `code_preprocess/project_to_fault_stress.py:bulk_stress_tensor_field`);
+      /// Phase 5's `build_stress_safs.py:evaluate_stress_field_on_grid`
+      /// is a pure MPa → Pa unit conversion;
       /// this reader is a pure pass-through.
```

PHASE6_DEVIATIONS.md uses a similar mid-flip phrasing; recommend
the same correction.

**Test case:** N/A (documentation correctness).

---

## Summary
- Critical issues:   **0**
- Moderate issues:   **2**
  - R-901 — no unit test for `StressFieldCoefficient` (silent
    Voigt vs schema-v1 ordering bug would not be caught)
  - R-902 — acceptance criteria #1 and #2 (cross-phase Python ↔ C++
    contract) not exercised
- Low issues:        **4**
  - R-903 — `BBox()` deviates from plan-literal "intersection"
  - R-904 — `ProjectStress` missing `scale`/`offset` parameters
    relative to `Project`/`ProjectVelocity`
  - R-905 — six redundant mesh-bbox containment all-reduces per
    `ProjectStress` call
  - R-906 — plan / deviation doc misattribute the sign-flip site
    (source docstring is correct)
- Plan compliance:   **PARTIAL** — every Phase 6 §1, §2, §3
  requirement implemented; the gap is in TEST coverage (R-901,
  R-902), not in production code.
- Verdict:           **PASS WITH FIXES** — no CRITICAL or
  blocking-MODERATE bugs. The two MODERATE findings are test-coverage
  gaps that should land before Phase 7 (which consumes the Tranche-1
  outputs and would surface any contract drift downstream).

## Unreviewed Areas
- **Phase 6 deferred subscopes (Phase 6.A, §4, §5, §6, §7):** out
  of this cycle per `PHASE6_DEVIATIONS.md`. These touch the BP5
  hot path and will be reviewed in their own implement → review →
  fix cycle (after they are implemented).
- **Real-data integration smoke (Phase 5 acceptance #1 against
  the SAFS 500/1000/2000 m bulk meshes):** test exists at the
  Python level (`TestRealDataSmoke` in `test_build_stress_safs.py`)
  and passes locally; not re-derived here.
- **MFEM library build status:** the SAFS-tree MFEM is not built;
  the new C++ code was verified to compile and link against
  `/Users/chunhuizhao/projects/seas-mfem/libmfem.a`. Once the
  in-tree library is built, `make test-stress-field-3d` will work
  out of the box. Not a code issue; flagged for environment
  awareness only.
- **Phase 1–4 / Phase 5 source:** previously reviewed (R-701..R-706
  and R-801..R-806). Re-running the existing test suites in this
  cycle confirmed **170/170 still passing** — no regression.
- **Sidecar-attr propagation (Phase 5 → Phase 6):** verified at
  the unit-test level (synthetic sidecars work) but not against
  the real Phase 5 binary — R-902 above is the missing
  integration test.
