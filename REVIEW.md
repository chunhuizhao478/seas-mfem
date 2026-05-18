# Code Review: Phases 6 (Tranches 2–3), 7, 8 — round 2 verification (post-R-001 fix)

## Review Scope
- Plan: `miniapps/seas/safs/project_7.0_alternative/code_preprocess/data_projection_onfaultstress/PLAN_onfaultstress.md`
- Prior review: `REVIEW.md` round 1 (R-001 through R-009, six fixes applied per `REVIEW_fix.md`)
- Files reviewed (fresh, full three-pass audit):
  - `miniapps/seas/fault/fault_geometry.hpp` (post-R-001 Gram-Schmidt fix)
  - `miniapps/seas/fault/fault_geometry_safs.inl`
  - `miniapps/seas/fault/rate_state_fault.hpp` (post-R-003 debug asserts)
  - `miniapps/seas/io/field_coefficient.{hpp,cpp}` (post-R-005 error msg fix)
  - `miniapps/seas/domain/domain_operator.hpp` (post-R-006 fallback removed)
  - `miniapps/seas/domain/elasticity_operator.hpp`
  - `miniapps/seas/domain/elasticity_operator_traction.inl`
  - `miniapps/seas/config/seas_config.hpp`
  - `miniapps/seas/config/seas_config_parser.hpp`
  - `miniapps/seas/config/seas_config_bridge.hpp`
  - `miniapps/seas/drivers/project_stress_to_mesh.cpp`
  - `miniapps/seas/safs/.../verify_onfault_stress.py` (post-R-009 warning fix)
  - All 7 C++ unit tests (Phase 6 Tranche 1 + 6.A/§4/§5/§6/§7 + Phase 7)
  - Phase 8 pytest (`test_verify_onfault_stress.py`)
- Test re-run on 2026-05-12:
  - C++: **116 / 116** (Phase 6 Tranche 1 = 56; Phase 6.A = 18; §4 = 12; §5 = 13; §6 = 5; §7 = 9; Phase 7 = 3)
  - Python: **182 / 182** (Phase 1-4 = 129; Phase 5 = 44; Phase 8 = 9)
- Domain context: `miniapps/seas/CLAUDE.md`, FaultBasis Tandem convention
  (`fault/fault_basis.hpp`), R-501/R-502 sign-flip rule,
  `bp5_params.hpp:tau0_vec` BP5 frame convention,
  `feedback_complete_sign_sites.md`, `PHASES_6_7_8_COMPLETION.md`,
  `REVIEW_fix.md`.

---

## Findings

### [R-101] LOW [fault/fault_geometry.hpp:789-791] — Stale docstring still describes the buggy `t2 = n × t1` Gram-Schmidt step

**Category:** QUALITY (documentation drift)

**Description:**
The class-level docstring for `ComputePerDOFCoordsAndBasis_` at lines
788-791 still reads:

```cpp
/// Mirrors Phase 2's `basis_to_node` convention:
///   n_i  ← n_i / |n_i|
///   t1_i ← (t1_i − (t1_i·n_i) n_i),  then  t1_i ← t1_i / |t1_i|
///   t2_i ← n_i × t1_i  (already unit, |t2_i| ≤ 1 by orthogonality)
```

But after the R-001 fix, the implementation does Gram-Schmidt
projection of the **input** `t2` against `(n, t1)`, NOT the cross
product `n × t1` shown here. A future maintainer reading this docstring
would believe the buggy convention is the intended one and might
"correct" the implementation back to the buggy form.

**Trigger:** Future maintainer reads the docstring and changes the
code to match it.

**Actual behavior:** Docstring contradicts implementation.

**Expected behavior:** Docstring describes the post-fix Gram-Schmidt-of-input-t2 algorithm.

**Suggested fix:**
```diff
    /// @brief Phase 6.A: populate per-DOF global 3-D coords and (n, t1, t2)
    /// basis from the domain operator, restrict to owned DOFs, and
    /// re-orthonormalise via Gram-Schmidt.
    ///
-   /// Mirrors Phase 2's `basis_to_node` convention:
-   ///   n_i  ← n_i / |n_i|
-   ///   t1_i ← (t1_i − (t1_i·n_i) n_i),  then  t1_i ← t1_i / |t1_i|
-   ///   t2_i ← n_i × t1_i  (already unit, |t2_i| ≤ 1 by orthogonality)
+   /// Gram-Schmidt re-orthonormalises ALL THREE input vectors so the
+   /// resulting (n, t1, t2) matches the elasticity operator's per-face
+   /// FaultBasis sign convention element-wise (R-001 contract: the
+   /// SAFS pre-stress must live in the same frame as the elastic
+   /// traction produced by FaultBasis::ProjectTraction):
+   ///
+   ///   n_i  ← n_i / |n_i|
+   ///   t1_i ← (t1_i − (t1_i·n_i) n_i),  then  t1_i ← t1_i / |t1_i|
+   ///   t2_i ← (t2_i − (t2_i·n_i) n_i − (t2_i·t1_i) t1_i),
+   ///                                          then  t2_i ← t2_i / |t2_i|
+   ///
+   /// Note: the cross-product form `t2 = n × t1` *loses* the
+   /// FaultBasis sign-flip when `sign_flipped == true` (see
+   /// fault_basis.hpp:463-471); projecting the input t2 preserves it.
    ///
    /// Degenerate-case fallback (plan §1768-1773): if the projected t1
```

**Test case:** N/A (documentation only).

---

### [R-102] LOW [safs/.../verify_onfault_stress.py:_cell_centroids, _cell_data_array] — Mixed-cell-type meshes can produce shape mismatch in `_verify_field`

**Category:** EDGE_CASE / POSSIBLE BUG

**Description:**
`_cell_centroids(m)` iterates **every** cell block in the meshio mesh
(line 388-393), but `_cell_data_array(m, name)` returns
`np.concatenate(m.cell_data.get(name))` (line 401), which only contains
arrays for cell blocks where the named field is defined. If the bulk
or fault VTU contains multiple cell types (e.g., tetrahedra for the
volume + triangles for the embedded fault surface), and `sigma_xx_Pa`
is only defined on tetrahedra, then:

- `_cell_centroids(m)` returns `N_tetra + N_triangle` centroids.
- `_cell_data_array(m, 'sigma_xx_Pa')` returns `N_tetra` values.

`_verify_field` would then raise `ValueError: shape mismatch on
'sigma_xx_Pa' — observed (N_tetra,) vs predicted (N_tetra +
N_triangle, 3)`, which would crash the verifier instead of producing a
useful diagnostic.

The current pytest fixture uses tetrahedra-only meshes, so this is not
triggered in CI — but the verifier is supposed to be resilient to "VTU
files written with either point-data or cell-data alone" (plan §2209)
which implies it should handle mixed-block meshes gracefully.

**Trigger:** A bulk VTU containing both tetrahedra and boundary
triangles where the σ_*_Pa cell-data is only populated on tetrahedra.

**Actual behavior:** `_verify_field` raises ValueError with a confusing
shape message.

**Expected behavior:** The function should restrict `_cell_centroids`
to the same cell blocks where the named field is present, OR raise an
informative warning and skip the field.

**Suggested fix:**
Refactor `_cell_centroids` and `_cell_data_array` to optionally filter
on a per-field basis:

```diff
 def _cell_centroids(m) -> np.ndarray:
     """Compute centroids (cell-mean of node coords) for all cell blocks."""
     pts = m.points
     out = []
     for cb in m.cells:
         out.append(np.mean(pts[cb.data], axis=1))
     if not out:
         return np.zeros((0, 3))
     return np.concatenate(out, axis=0)


+def _cell_centroids_with_field(m, name: str) -> np.ndarray:
+    """Compute centroids of cell blocks that carry the named field."""
+    pts = m.points
+    field_arrs = m.cell_data.get(name)
+    if field_arrs is None:
+        return np.zeros((0, 3))
+    out = []
+    for cb, arr in zip(m.cells, field_arrs):
+        if arr is None or arr.size == 0:
+            continue
+        out.append(np.mean(pts[cb.data], axis=1))
+    if not out:
+        return np.zeros((0, 3))
+    return np.concatenate(out, axis=0)
+
+
 def _cell_data_array(m, name: str) -> Optional[np.ndarray]:
```

Then in `verify_bulk_cell_data` / `verify_fault_cell_data`, replace the
single `centroids = _cell_centroids(m)` with a per-field
`_cell_centroids_with_field(m, name)` call inside the loop.

**Test case:**
```python
def test_R102_mixed_cell_blocks(tmp_path, safod_params):
    """Mixed tetra + triangle bulk VTU with sigma_*_Pa only on tetra
    must not crash with a ValueError shape mismatch."""
    import meshio
    pts = np.array([
        [0,0,0],[1,0,0],[0,1,0],[0,0,-1],   # one tetra
        [0,0,0],[1,0,0],[0,0,-1],            # one boundary triangle
    ])
    cells = [
        ("tetra",    np.array([[0,1,2,3]])),
        ("triangle", np.array([[4,5,6]])),
    ]
    # Only define the cell-data on the tetra block.
    cell_data = {"sigma_xx_Pa": [np.array([1.0e7]), None]}
    out = tmp_path / "mixed.vtu"
    meshio.write_points_cells(str(out), pts, cells, cell_data=cell_data)
    # Should not raise.
    results = vof.verify_bulk_cell_data(out, safod_params)
    # Either the field is verified on the tetra (1 cell) or skipped
    # with a warning; in neither case should an exception propagate.
```

---

### [R-103] LOW [fault/fault_geometry.hpp:879-925] — Degenerate-t1 fallback preserves the original input t2 instead of re-deriving it, may abort on degenerate input

**Category:** EDGE_CASE (defensive code path)

**Description:**
The fallback branch at lines 879-894 (entered when the input t1 is so
close to zero that Gram-Schmidt-ing it against n produces |t1_proj| <
1e-12) re-derives t1 from the reference up-vector via
`t1 = up − (up·n) n`. After this fallback `t1` is fresh, but the code
then Gram-Schmidt-projects the **original** input `t2` against the new
(n, t1_fallback) at lines 914-925.

If the input t1 was degenerate, the input `t2` from FaultBasis is
almost certainly degenerate too (FaultBasis derives `t1 = dip = strike
× n` and `t2 = strike = up × n`; both go bad when `up × n` is near
zero). The Gram-Schmidt projection of a near-zero t2 yields a
near-zero t2_proj, which then trips the MFEM_VERIFY at line 921.

In practice this fallback is unreachable from a clean FaultBasis
because `FaultBasis::ComputeOrientedFrame` already aborts on `s_len >
1e-12` at fault_basis.hpp:413 — but the fallback is still in the code
as a defensive path, and as written it cannot succeed when the input
is the kind of degenerate input it was meant to handle.

**Trigger:** Hypothetical mesh where the FaultBasis abort threshold is
relaxed and the operator pipes a face with `up · n ≈ ±1` (horizontal
face) through to FaultGeometry.

**Actual behavior:** MFEM_VERIFY abort at line 921.

**Expected behavior:** In the fallback case, re-derive t2 from
`n × t1_fallback` (sign convention is moot because the original
FaultBasis sign convention was undefined for the degenerate face).

**Suggested fix:**
Track whether the fallback fired and use `n × t1` in that case:

```diff
          if (t1_len < 1e-12)
          {
             // ...
             num_dof_basis_fallbacks_++;
+            // Fallback path: input t2 was associated with a degenerate
+            // FaultBasis t1, so derive t2 from the freshly-built
+            // (n, t1_fallback) frame.  Sign convention is moot — the
+            // original FaultBasis sign was undefined for this face.
+            t2[0] = n[1]*t1[2] - n[2]*t1[1];
+            t2[1] = n[2]*t1[0] - n[0]*t1[2];
+            t2[2] = n[0]*t1[1] - n[1]*t1[0];
          }
          const real_t inv_t1 = 1.0 / t1_len;
          t1[0] *= inv_t1; t1[1] *= inv_t1; t1[2] *= inv_t1;
@@
-         // R-001 (post-fix): Gram-Schmidt re-orthonormalise the input
-         // t2 against (n, t1) ...
+         // R-001 (post-fix): Gram-Schmidt re-orthonormalise the input
+         // t2 against (n, t1).  In the t1-fallback branch above, t2
+         // was already re-derived as n × t1_fallback (the original
+         // input t2 is unreliable when its companion t1 was degenerate).
+         // ...
          const real_t t2_dot_n  = t2[0]*n[0]  + t2[1]*n[1]  + t2[2]*n[2];
          const real_t t2_dot_t1 = t2[0]*t1[0] + t2[1]*t1[1] + t2[2]*t1[2];
```

**Test case:** N/A (defensive path; not reachable from a passing
FaultBasis::ComputeOrientedFrame in current code).

---

### [R-104] LOW [domain/elasticity_operator_traction.inl:GetFaultDOFCoords3D / GetFaultDOFBasis] — Substantial code duplication between interior-face and shared-face loops in two new methods

**Category:** QUALITY (duplication)

**Description:**
Both `GetFaultDOFCoords3D` (lines 144-206) and `GetFaultDOFBasis` (lines
208-310) have nearly identical interior-face and shared-face loop
bodies. For each method, the difference between the two loops is only
which face-transformation accessor is called
(`GetInteriorFaceTransformations` vs `GetSharedFaceTransformations`).
The per-DOF coordinate / basis computation is identical.

This is a maintenance risk: a future change to the per-DOF logic must
be applied twice in each method (four sites total), and any divergence
silently breaks parallel-only paths.

**Suggested fix:**
Extract a helper that takes a face-transformation pointer and the
output offset:

```cpp
template <typename MeshType>
void ElasticityDomainOperator<MeshType>::WriteFaultDOFCoords3DToFace_(
   FaceElementTransformations *FTr, int face_dof_offset, int nbf,
   const IntegrationRule &nir, Vector &dof_coords_3d) const
{
   for (int kk = 0; kk < nbf; kk++)
   {
      const IntegrationPoint &nip = nir.IntPoint(kk);
      FTr->SetAllIntPoints(&nip);
      Vector coords(3);
      FTr->Elem1->Transform(FTr->GetElement1IntPoint(), coords);
      const int dof_idx = face_dof_offset + kk;
      dof_coords_3d(3 * dof_idx + 0) = coords(0);
      dof_coords_3d(3 * dof_idx + 1) = coords(1);
      dof_coords_3d(3 * dof_idx + 2) = coords(2);
   }
}
```

Then the interior and shared loops each call this helper once instead
of duplicating its body. Similar refactor for `GetFaultDOFBasis`.

**Test case:** N/A (style; the existing tests cover both branches via
the parallel-mesh test fixture).

---

### [R-105] POSSIBLE LOW [drivers/project_stress_to_mesh.cpp:main] — `ProjectStress` may write partial output before aborting on bbox containment failure

**Category:** POSSIBLE EDGE_CASE

**Description:**
`FieldProjector::ProjectStress` performs the mesh-bbox containment
check inside each of the six per-component `Project` calls. If the
mesh extends outside the sidecar's bbox, only the FIRST component's
`Project` aborts — by the time it does, no `.pvd` / `.vtu` files have
been written (the writers run after the rotation, line 149-166 of the
driver). So there's no partial-output risk in the driver as written.

However, the `--out` directory may already have been touched by
`ParaViewDataCollection::SetPrefixPath` (mkdir-on-construct in some MFEM
versions). If the abort fires after the directory is created, the user
sees a leftover empty output directory — minor UX issue.

**Trigger:** Run the driver with a mesh whose bbox extends outside the
sidecar bbox.

**Actual behavior:** Driver aborts cleanly via MFEM_ABORT, but may
leave an empty output directory behind.

**Expected behavior:** Either pre-flight the bbox check before any
filesystem mutation, or document the leftover-directory behaviour.

**Suggested fix:**
Hoist the bbox check before constructing the ParaViewDataCollection by
calling `StressField3D` directly and invoking `ContainsBBox` against
the mesh's global bbox:

```diff
    mfem::H1_FECollection fec(order, pmesh.Dimension());
    mfem::ParFiniteElementSpace fes(&pmesh, &fec);

+   // Pre-flight: open the sidecar and check bbox containment BEFORE
+   // touching the filesystem.
+   {
+      mfem::seas::StressField3D probe(sidecar_path);
+      probe.SetInterpMode(interp_mode);
+      real_t mxmin, mxmax, mymin, mymax, mzmin, mzmax;
+      // ... compute mesh bbox via Allreduce ...
+      if (!probe.ContainsBBox(...)) {
+         if (rank == 0)
+            std::cerr << "ERROR: mesh bbox outside sidecar bbox.\n";
+         MPI_Finalize();
+         return 1;
+      }
+   }
+
    const long long global_ne = pmesh.GetGlobalNE();
```

**Test case:** N/A (UX-level diagnostic; existing tests pass).

---

## Verification of prior round-1 fixes

Each round-1 finding was re-audited against current source. Status:

| Finding | Status |
|---------|--------|
| R-001 (CRITICAL, sign-flip drop) | **FIXED** — Gram-Schmidt of input t2 preserves FaultBasis sign convention. T_6A_6 + T_64_8 are now non-tautological and pass. |
| R-002 (MODERATE, tautological T_64_2) | **FIXED** — uses FaultBasis as independent reference. |
| R-003 (MODERATE, lifetime contract) | **FIXED** — MFEM_ASSERT debug checks; BP5 bit-exact T_66_2 still passes. |
| R-004 (MODERATE, BP5 bit-exact) | **DEFERRED** per user instruction (requires full BP5 driver run). |
| R-005 (LOW, BP2 error msg) | **FIXED** — `MFEM_ABORT` with clear message at top of ProjectFaultPreStress. |
| R-006 (LOW, GetFaultDOFCoords3D default) | **FIXED** — default returns empty vector. |
| R-007 (LOW, ComputeSAFSParams irreversible) | **DEFERRED** per user instruction. |
| R-008 (LOW, Phase 7 round-trip values) | **DEFERRED** per user instruction. |
| R-009 (LOW, Phase 8 silent skip) | **FIXED** — `print(... file=sys.stderr)` warning. |

---

## Summary
- Critical issues:   **0**
- Moderate issues:   **0**
- Low issues:        **5** (R-101 through R-105 — all introduced /
  uncovered by the round-1 fix audit; none affect correctness)
- Plan compliance:   **FULL** — every Phase 6 (§4–§7), Phase 7, and
  Phase 8 deliverable is implemented; the previously-flagged
  sign-convention bug R-001 is closed.
- Test suite:        **PASS — 116 / 116 C++ + 182 / 182 Python**.
- Verdict:           **PASS WITH FIXES** — every CRITICAL/MODERATE
  finding from round 1 is resolved or explicitly deferred. The new
  round-2 findings are LOW-severity documentation, defensive-code, or
  style issues; none block merging. R-101 (stale docstring) and
  R-102 (mixed-cell-block resilience) should be addressed before the
  next release; R-103, R-104, R-105 may be deferred to a future
  refactor cycle.

## Unreviewed Areas
- **Full BP5 simulation regression** (plan §1972-1973). Still pending
  on a built SAFS-tree libmfem.a; tracked under deferred R-004.
- **MPI-parallel paths** (ParMesh / `MFEM_USE_MPI`). Compiles but not
  exercised at runtime in this environment.
- **Real SAFS sidecar (`stress_safs.h5`) projection on the 500 / 1000 /
  2000 m meshes** — verifier and projector exercised only against
  synthetic constant-value sidecars.
- **Linker / Makefile compatibility once SAFS-tree libmfem.a is built**
  — verified compilation but not end-to-end `make test` flow.
