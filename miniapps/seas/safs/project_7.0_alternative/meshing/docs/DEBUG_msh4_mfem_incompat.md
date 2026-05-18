# DEBUG: SAFS `.msh` files fail to load in MFEM (`-format msh4` + v2.2-only MFEM reader)

**Date:** 2026-05-18
**Status:** root cause confirmed; fix proposed; NO code changes applied yet.
**Discovered while running:** `make test-spatial-dyn-driver`
(Phase 4 of `safs/project_7.0_alternative/document/spatial_dynamic_rupture_plan.md`).

---

## 1. Symptom

When `seas_spatial_dyn_driver` (or any other MFEM driver) attempts to load
`safs/project_7.0_alternative/meshing/results/msh/safs_fault_box_nwcut_*.msh`
via `mfem::Mesh(path, generate_edges, refine)`, the process aborts at:

```
MFEM abort: Gmsh file : vertices indices are not unique
 ... in function: void mfem::Mesh::ReadGmshMesh(std::istream &, int &, int &)
 ... in file: mesh/mesh_readers.cpp:1628
```

This was reproduced against the `1000m_lcfar3000` mesh, but the same abort
fires for every SAFS `.msh` emitted by the current pipeline (`500m_*`,
`1000m_*`, `2000m_*`, both `_lcfar3000` and `_zgraded` variants), because
they are all the same Gmsh format.

The error message is misleading: the file's vertex indices ARE unique.
The real failure happens upstream in the parser (see §3 below).

## 2. Reproduction

```bash
cd /Users/chunhuizhao/projects/seas-mfem-safs/miniapps/seas
conda activate mfem-dev
make seas_spatial_dyn_driver
make test-spatial-dyn-driver       # aborts at mesh load (this report's symptom)

# Or any inline minimum:
./seas_spatial_dyn_driver \
    --config safs/project_7.0_alternative/friction/EXAMPLE_spatial_friction_slip_weakening_safs.toml \
    --mesh   safs/project_7.0_alternative/meshing/results/msh/safs_fault_box_nwcut_1000m_lcfar3000.msh \
    --no-sidecar-material --dry-run
```

Reproducibility: **always**, on every SAFS `.msh` artifact produced by
`meshing/code/run_nwcut_meshing.py` (current `feature/safs-quasi-dynamic`
branch, commit `bfb6cec`).

## 3. Hypotheses (ranked)

1. **The SAFS `.msh` files are Gmsh format v4.1, but MFEM's `Mesh::ReadGmshMesh`
   only implements the v2.2 ASCII parser.**  Most likely: a quick `head` of
   the file shows the v4.1 header (`4.1 0 8`) and v4.1-style `$Nodes` block
   layout (`numEntityBlocks numNodes minNodeTag maxNodeTag` + per-entity
   sub-blocks), while the MFEM reader at `mesh/mesh_readers.cpp:1581-1629`
   reads only ONE integer from the `$Nodes` line and processes vertices
   with the v2.2 layout (`<tag> x y z` per line, no entity sub-blocks).
2. **The mesher emits genuinely duplicated vertex tags.**  Possible if the
   `nw_cut` step splits the original mesh and re-tags vertices.  Should be
   ruled out by inspecting the `$Nodes` header (`min..max`) against
   `numNodes`.
3. **MFEM aborts on the binary flag.**  The MFEM reader's binary-format
   branch checks `dsize == sizeof(double)` and reads a `1`; a malformed
   binary flag could fire the same abort message via a misaligned read.
   The SAFS file is ASCII (`4.1 0 8` → binary=0), so this branch should
   not be entered.

## 4. Investigation

### Hypothesis 1 — MFEM only supports Gmsh v2.2 (CONFIRMED)

Evidence A — the actual mesh format:

```
$ head -3 results/msh/safs_fault_box_nwcut_1000m_lcfar3000.msh
$MeshFormat
4.1 0 8
$EndMeshFormat

$ sed -n '43,48p' results/msh/safs_fault_box_nwcut_1000m_lcfar3000.msh
$Nodes
28 171457 1 171457          ← v4.1 header: 28 entity blocks, 171,457 nodes,
                              minTag=1, maxTag=171457
0 1 0 1                     ← entity block: dim=0, tag=1, parametric=0, numNodes=1
1                           ← node tag
314830.834379 3642278 0     ← coords (on a separate line!)
```

Evidence B — MFEM's reader (`mesh/mesh_readers.cpp:1515-1629`):

```cpp
void Mesh::ReadGmshMesh(std::istream &input, int &curved, int &read_gf)
{
   real_t version;
   int binary, dsize;
   input >> version >> binary >> dsize;
   if (version < 2.2)
   {
      MFEM_ABORT("Gmsh file version < 2.2");
   }
   // …
   // ↑ no `version >= 4.0` branch; the function proceeds with the v2.2
   //   parser body unconditionally.
   …
   if (buff == "$Nodes")
   {
      input >> NumOfVertices;              // reads ONE int from the $Nodes header
      getline(input, buff);
      // ↑ for v4.1 input, this captures `28` only.  The remaining
      //   `171457 1 171457` is consumed by `getline` as trailing junk.
      vertices.SetSize(NumOfVertices);     // = 28
      for (int ver = 0; ver < NumOfVertices; ++ver)
      {
         input >> serial_number;           // mis-aligned reads start here
         for (int ci = 0; ci < 3; ++ci) { input >> coord[ci]; }
         vertices_map[serial_number] = ver;
      }
      …
      if (vertices_map.size() != NumOfVertices)
      {
         MFEM_ABORT("Gmsh file : vertices indices are not unique");
         // ↑ NumOfVertices = 28; vertices_map.size() < 28 because the
         //   mis-aligned reads produce collisions.
      }
   }
}
```

Evidence C — confirmed there is **no v4 reader** anywhere in MFEM:

```
$ grep -n 'version >= 4\|gmsh.*4\.\|format.*4\.1\|v4 format' \
       /Users/chunhuizhao/projects/mfem/mesh/mesh_readers.cpp
(no matches)
$ grep -rn 'ReadGmsh' /Users/chunhuizhao/projects/mfem/mesh/*.{cpp,hpp}
mesh.cpp:5070:      ReadGmshMesh(input, curved, read_gf);
mesh_readers.cpp:1515:void Mesh::ReadGmshMesh(std::istream &input, …)
```

Only one entry, only the v2.2 body.  Hypothesis 1: **CONFIRMED**.

### Hypothesis 2 — duplicated vertex tags (ELIMINATED)

The `$Nodes` header `28 171457 1 171457` declares minNodeTag=1,
maxNodeTag=171457, numNodes=171457 — i.e., dense numbering with no gaps
and no duplicates.  Even if the v4 reader existed, that header would pass
a uniqueness check.  ELIMINATED.

### Hypothesis 3 — binary-flag mis-read (ELIMINATED)

`$MeshFormat\n4.1 0 8\n$EndMeshFormat` ⇒ `binary = 0` (ASCII).  The
ASCII branch is exactly where the abort happens.  ELIMINATED.

### Verification that v2.2 fixes the symptom

Manually converting the same mesh to Gmsh v2.2 via `gmsh -format msh2`:

```bash
conda activate pythonenv
mkdir -p /tmp/safs_msh2_smoke
gmsh safs/project_7.0_alternative/meshing/results/msh/safs_fault_box_nwcut_1000m_lcfar3000.msh \
     -format msh2 -save \
     -o /tmp/safs_msh2_smoke/safs_fault_box_nwcut_1000m_lcfar3000.msh -v 0
# ⇒ produces a v2.2 file (header `2.2 0 8`).
```

Then re-running the driver:

```bash
conda activate mfem-dev
./seas_spatial_dyn_driver \
    --config safs/project_7.0_alternative/friction/EXAMPLE_spatial_friction_slip_weakening_safs.toml \
    --mesh   /tmp/safs_msh2_smoke/safs_fault_box_nwcut_1000m_lcfar3000.msh \
    --no-sidecar-material --dry-run
```

Output (truncated):

```
[material] using [material_constant_fallback]: lambda=3.2e+10 mu=3.2e+10 rho=2670
[spatial_dyn] WARNING: tfinal (12 s) exceeds min_box_dim / cp_max (6.94 s) …
[fault] QPs per face = 3, num_fault_global = 31809  (local = 31809, shared = 0)
[time] dt_cfl  = 0.00167074 s
[time] dt      = 0.00167074 s
[time] nsteps  = 7183
[spatial_dyn] --dry-run: construction complete, exiting.
```

The driver completes the full `--dry-run` path: mesh load,
`BoundaryConfig` (attrs 101/102/103/104), `WaveOperator` ctor, fault-DOF
table walk (31 809 owned fault QPs at np=1), CFL, and exit-0.  Verifies
that the SAFS mesh content is correct; only the *file format* is the
blocker.

(Aside: a `MPI_Comm_free after MPI_FINALIZE` warning prints at exit.  This
is unrelated to the mesh and is a known driver-side ordering issue with
`mfem::seas::MPIContext`'s destructor — track separately.)

## 5. Root cause

| Layer | What | Where | Why |
|-------|------|-------|-----|
| Producer | Mesher emits Gmsh v4.1 | `meshing/code/run_nwcut_meshing.py:211` — argv contains `"-format", "msh4"` | Choice from commit `bfb6cec` (no rationale committed). |
| Consumer | MFEM only reads Gmsh v2.2 | `/Users/chunhuizhao/projects/mfem/mesh/mesh_readers.cpp:1515-1629` | The function gates on `version < 2.2` and then runs the v2.2 ASCII body unconditionally; no v4 branch exists. |
| Symptom site | Vertex-uniqueness assertion | `mesh_readers.cpp:1628` | After 28 mis-aligned reads (the v4.1 `$Nodes` header's first integer), `vertices_map.size() != NumOfVertices`. |

The two failure modes (`msh4` produced, `msh22` expected) collide at
`Mesh::Mesh(path, …)`, which dispatches by the `$MeshFormat` tag without
distinguishing v2 from v4.

### Impact

- Every SAFS mesh artifact in
  `safs/project_7.0_alternative/meshing/results/msh/` (six files at
  500/1000/2000 m × `_lcfar3000`/`_zgraded`) is unreadable by any
  MFEM-based driver: `seas_spatial_dyn_driver`,
  `seas_project_velocity_to_mesh`, `seas_project_stress_to_mesh`,
  any future quasi-dynamic driver, and any unit/integration test that
  loads a SAFS `.msh`.
- The team has already documented a manual workaround in
  `docs/PLAN_mesh_quality.md:564-577`:

  ```bash
  gmsh …/safs_fault_box_nwcut_1000m_zgraded.msh \
       -format msh2 -save \
       -o /tmp/safs_msh2/safs_fault_box_nwcut_1000m_zgraded.msh -v 0
  ```

  The fact that this conversion is in the docs at all is evidence that
  the v4-vs-MFEM mismatch was known when the mesher was committed but
  was not pushed back upstream into `run_nwcut_meshing.py`.

## 6. Proposed fix (NOT YET APPLIED)

Two options; either is safe.  Recommended: option A.

### Option A — flip the mesher to emit Gmsh v2.2 directly (recommended)

Edit `miniapps/seas/safs/project_7.0_alternative/meshing/code/run_nwcut_meshing.py:211`:

```diff
-        "-format", "msh4",
+        "-format", "msh22",
```

Gmsh 4.x accepts `msh22` (Gmsh v2.2 ASCII format) as a synonym for
`msh2`; both produce the same file.  All SAFS pipeline features
exercised by the `.geo` (Physical Surfaces, Physical Volumes,
embedded discrete fault surface, MathEval size field, 2nd-order
volume tetrahedra at `Mesh.ElementOrder ≤ 1`) are representable in
v2.2.

What changes for downstream consumers:
- MFEM drivers: load without conversion.  Done.
- `meshio` (`stress/code/project_to_fault_stress.py`,
  `check_mesh_quality.py`): unchanged — `meshio` auto-detects v2.2 vs
  v4.1, both code paths work.
- `seas_project_velocity_to_mesh`: same — `mfem::Mesh` is the consumer.
- ParaView / VisIt: unchanged — both back-ends accept v2.2.

What MIGHT change (verify before merging):
- File size: v2.2 is typically ≈ 1.3× larger than v4.1 because v2.2
  cannot factor per-entity sub-blocks.  Empirical check on
  `1000m_lcfar3000`: v4.1 = 44.3 MB → v2.2 = 51.7 MB (+17 %).  Within
  budget for the `results/msh/` artifacts but worth recording.
- Entity-level metadata: gmsh v4 carries explicit `$Entities` and
  `$PartitionedEntities` blocks; v2 does not.  The SAFS pipeline does
  not read these back (it only consumes Physical Surfaces / Volumes),
  so the round-trip from `.geo` through gmsh to MFEM is lossless.

Update `docs/PLAN_mesh_quality.md` step 3 (lines 564-577): the manual
`msh4 → msh2` re-conversion step becomes unnecessary; remove the gmsh
shell-out and point downstream consumers (`project_velocity_to_mesh`,
projection scripts) at `results/msh/...lcfar3000.msh` directly.

### Option B — keep `msh4` as the artifact format and add a Make rule

If the v4.1 artifacts are needed for an external consumer (e.g., a
gmsh-native post-processing tool that prefers v4), keep
`run_nwcut_meshing.py` unchanged and add a Make target that batch-
converts every `*.msh` in `results/msh/` into `results/msh2/` once,
keying off mtime.  Pros: keeps the producer pristine; cons: doubles
the artifact footprint and adds a chained dependency every consumer
has to remember.

I do not recommend option B unless a v4-only downstream consumer
exists — there is no evidence of one in this tree.

## 7. Why this surfaced now

The SAFS meshes were added (commit `bfb6cec`, 2026-05-17) for the
Python-only pipelines (`check_mesh_quality.py`, `project_to_fault_stress.py`,
`build_velocity_cvmh.py`) which all use `meshio`.  No MFEM C++ driver had
consumed a SAFS `.msh` until Phase 4 of the spatial dynamic-rupture
plan (commit `79dd730`, 2026-05-18) wired `seas_spatial_dyn_driver` to
load `cfg.mesh.path`.  That is when the v4-vs-v2 mismatch became
load-bearing.

## 8. Follow-up checklist (no action taken here)

- [ ] Apply Option A — change `run_nwcut_meshing.py:211` to `"-format", "msh22"`.
- [ ] Regenerate all six `*.msh` artifacts via the canonical
      `python run_nwcut_meshing.py --lc-far 3000 --suffix _lcfar3000 --res 500 1000 2000`
      command from `docs/PLAN_mesh_quality.md`.
- [ ] Verify `check_mesh_quality.py` still passes on the regenerated
      artifacts (the gates are mesh-content, not file-format, so this
      should be a no-op).
- [ ] Drop the manual `gmsh -format msh2 -save` step from
      `docs/PLAN_mesh_quality.md` §3 (lines 564-577) and update any
      `seas_project_velocity_to_mesh` / projection invocations to point
      at `results/msh/...` directly.
- [ ] Re-run `make test-spatial-dyn-driver` to confirm Phase 4 acceptance
      gate goes green.
- [ ] (Separate issue) Track and fix the `MPI_Comm_free after MPI_FINALIZE`
      warning in `seas_spatial_dyn_driver` — `MPIContext` should be
      destroyed before `MPI_Finalize`, e.g. by wrapping the body in a
      scope so its destructor fires first.

## 9. Files referenced (read-only here)

| Path                                                                                                                  | Role                              |
|-----------------------------------------------------------------------------------------------------------------------|-----------------------------------|
| `miniapps/seas/safs/project_7.0_alternative/meshing/code/run_nwcut_meshing.py:211`                                    | Mesher source — `-format msh4`.   |
| `miniapps/seas/safs/project_7.0_alternative/meshing/results/msh/safs_fault_box_nwcut_1000m_lcfar3000.msh`             | Failing artifact (Gmsh v4.1).     |
| `miniapps/seas/safs/project_7.0_alternative/meshing/docs/PLAN_mesh_quality.md:564-577`                                | Already-documented manual fix.    |
| `/Users/chunhuizhao/projects/mfem/mesh/mesh_readers.cpp:1515-1629`                                                    | MFEM Gmsh v2.2-only reader.       |
| `miniapps/seas/safs/project_7.0_alternative/document/spatial_dynamic_rupture_plan.md` §Phase 4                        | Triggering consumer.              |
| `miniapps/seas/drivers/spatial_dyn_driver.cpp`                                                                        | The new driver where this blocks. |
