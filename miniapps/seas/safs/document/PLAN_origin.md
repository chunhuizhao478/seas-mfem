# Implementation Plan: Origin Convention and CFM→Gmsh Mapping for SAFS Mesh

## Overview
This document fixes the choice of model origin (0,0,0), the rigid mapping from raw
CFM-UTM coordinates into the local Gmsh model frame, and the recovery path from local
back to UTM. It supersedes Phase 2 step 1 of `PLAN.md` (which proposed an unstable
centroid-based origin). The decision affects every script in the SAFS mesh pipeline and
must be agreed before any code is written.

**Decision in one sentence:** the model origin is a **fixed UTM grid point**
`(E0, N0, 0) = (500 000, 3 765 000, 0)` in UTM Zone 11N / NAD83, with z = 0 = sea level;
the model frame is local-Cartesian metres aligned with UTM E/N (no rotation, no scaling),
so the only transformation is a rigid translation of `(−E0, −N0, 0)`.

## Background: Why Not the Vertex-Cloud Centroid?

The first draft (`PLAN.md` Phase 2 step 1) proposed translating by the centroid of the
combined fault vertex cloud. That has three problems:

1. **The origin moves when the input set changes.** Adding the Pinto Mountain fault
   shifts the centroid by ~10 km. Re-running the cleanup with a different `--res`
   (denser triangulation in the SBMT-SAF region than in others) changes vertex weighting
   and shifts the centroid by O(100 m). Either way, the same UTM point maps to a
   different local coordinate, so meshes generated under different selections cannot be
   overlaid or compared without re-running everything.
2. **Centroid coordinates are non-round.** With the current 8 faults the centroid lands
   near `(493 691, 3 765 638)` — values that no human will recognize, recall, or sanity
   check.
3. **Recovery is asymmetric.** With a centroid-based origin, recovering UTM requires
   reading `transform.json`. With a fixed grid origin, the inverse is trivial mental
   arithmetic: `+500 km E, +3 765 km N`.

A fixed UTM grid origin solves all three.

## Decision: The Origin and Why These Numbers

```
Origin in UTM Zone 11N (NAD83):
    E0 =   500 000 m   (UTM Zone 11N false easting → central meridian, 117° W)
    N0 = 3 765 000 m   (5 km-rounded value near the fault-system geographic centroid)
    Z0 =         0 m   (sea level)
```

### Justification per coordinate

- **`E0 = 500 000 m`.** UTM Zone 11N's false easting *is* 500 000 m — by construction of
  the projection, `E = 500 000` is the central meridian (longitude 117° W). The central
  meridian is the line of zero scale distortion in UTM, and the SAFS faults straddle it
  (Mojave SAF at 365–428 km E sits west of the meridian; Coachella Mission Creek at
  541–571 km E sits east). Snapping `E0` to the central meridian therefore (a) places
  the model at the spot of minimum projection distortion, (b) keeps post-translation X
  coordinates as nearly symmetric as the data allow, and (c) is a "geographically named"
  reference any future user can recognize without reading our metadata.

  After translation, the per-fault X extents become:
  ```
  MJVS-SAF (Mojave):           X ∈ [−134 946, −71 974] m
  SBMT-SAF (San Bern):         X ∈ [ −72 032, +38 743] m
  SBMT Mill Creek:             X ∈ [ −30 212, +28 501] m
  SBMT Mission Creek:          X ∈ [  −7 019, +43 146] m
  MULT-Banning:                X ∈ [ −10 597, +56 569] m
  COAV Mission Creek:          X ∈ [ +41 475, +70 749] m
  ETRA Pinto Mtn:              X ∈ [ +24 384, +106 749] m
  MULT Southern SAF + Banning: X ∈ [ +57 182, +122 329] m
  ```
  Max |x| = 134.9 km. Bulk box default with `buf_x = 30 km` reaches |x| = 165 km.

- **`N0 = 3 765 000 m`.** The selected faults span Y ∈ [3 692 278, 3 838 998] m, midpoint
  3 765 638 m. We round to the nearest 5 km grid line (3 765 000 m). After translation:
  ```
  Y ∈ [−72 722, +73 998] m   (max |y| ≈ 74 km)
  ```
  Symmetry to within ~1 km. Snapping to 5 km gives a memorable value.

- **`Z0 = 0` (no translation).** CFM's `ZPOSITIVE Elevation` and the SEAS-MFEM project
  convention (per `miniapps/seas/CLAUDE.md`: "Z=0 at surface, Z<0 is depth") agree.
  Translating Z would break the project-wide sign convention used by `fault_basis.hpp`
  (`dip direction = (0, 0, +1) downward`). We keep Z = elevation in metres, then clamp
  any z > 0 to z = 0 in the cleanup step (Phase 2 of `PLAN.md`) to flatten residual
  topography to the free-surface plane.

### Sign convention (write this down once, refer to it forever)

```
v_local = v_utm − origin_utm                  (apply origin to UTM ⇒ local)
v_utm   = v_local + origin_utm                (recovery)

with origin_utm = (E0, N0, Z0) = (500 000, 3 765 000, 0)
```

X grows east, Y grows north, Z grows upward (elevation). No axis flips, no rotation. The
local frame is right-handed, identical to UTM Zone 11N's right-handed (E, N, Up) frame.

### Numerical conditioning check

OpenCASCADE has an internal default tolerance of about 1e−7 in working units. With raw
UTM coordinates (~3.8 × 10^6 m), relative tolerance becomes ~2.6 × 10^−14, dangerously
close to double-precision noise. After translation, max |coordinate| is ~135 km =
1.35 × 10^5 m, so relative tolerance becomes ~7.4 × 10^−13 — three orders of magnitude
of headroom restored. Translation is therefore *required*, not just convenient.

## Constraints

- **Rigid translation only.** No rotation, no scaling, no projection change. UTM Zone 11N
  is preserved; the local frame is UTM with the origin shifted.
- **Fixed origin constants.** `E0`, `N0`, `Z0` are **literals**, not derived from the
  current input data. They live as named constants in exactly one place
  (`mesh/safs_origin.py`, see Phase 1) and every other script imports from there.
- **Z is not translated.** Free surface remains z = 0 (sea level). This preserves the
  existing project sign convention for `fault_basis.hpp` and `bp5_params.hpp`.
- **`transform.json` is the wire format.** Any tool that needs to convert between local
  and UTM reads `transform.json`, never reaches into `safs_origin.py` directly. The JSON
  is the data-side contract; the Python module is the code-side source of truth that
  generates the JSON.
- **Per-feedback memory (`feedback_no_hardcoded_numbers.md`).** The origin literals
  `(500000, 3765000, 0)` are *constants of the geographic reference frame*, not numerical
  parameters of the simulation. They live in one named-constant location with a comment
  explaining the choice, and they may be overridden via CLI for sensitivity studies. This
  is the correct pattern: not buried magic numbers, but a single declaratively-named
  reference.
- **Override mechanism.** Every script that consumes the origin accepts `--origin-easting
  <E0>`, `--origin-northing <N0>`, `--origin-elevation <Z0>` flags whose defaults are the
  fixed values above. Sensitivity studies (e.g., does origin choice affect mesh quality?)
  pass non-default flags; production runs use defaults.

## Phase 1: Origin Module and Metadata Schema

### Goal
After this phase a single Python module defines the origin constants and provides the
forward/inverse mapping helpers, plus a JSON-schema-validated `transform.json` is written
during Phase 2 of `PLAN.md` and read by every downstream tool.

### Files to Create
- `miniapps/seas/safs/mesh/safs_origin.py` — the source-of-truth Python module.
- `miniapps/seas/safs/mesh/transform.schema.json` — JSON Schema for `transform.json`.

### Files to Modify
- `miniapps/seas/safs/PLAN.md` — replace Phase 2 step 1 ("Translation to local frame")
  with a one-line reference to this document; remove all mention of vertex-cloud centroid.

### Detailed Requirements

1. **`safs_origin.py` interface.**
   ```python
   # Module-level constants. Do NOT hardcode these values anywhere else.
   UTM_ZONE = "11N"
   DATUM    = "NAD83"
   E0       = 500_000.0    # UTM easting of model origin, metres
   N0       = 3_765_000.0  # UTM northing of model origin, metres
   Z0       = 0.0          # elevation of model origin, metres (sea level)

   ORIGIN_UTM: tuple[float, float, float] = (E0, N0, Z0)

   def utm_to_local(v_utm: np.ndarray,
                    origin: tuple[float, float, float] = ORIGIN_UTM) -> np.ndarray:
       """v_utm shape (N, 3) → v_local shape (N, 3). v_local = v_utm − origin."""

   def local_to_utm(v_local: np.ndarray,
                    origin: tuple[float, float, float] = ORIGIN_UTM) -> np.ndarray:
       """Inverse of utm_to_local."""

   def write_transform_json(path: str | Path,
                            origin: tuple[float, float, float] = ORIGIN_UTM,
                            extra: dict | None = None) -> None:
       """Write transform.json conforming to transform.schema.json."""

   def read_transform_json(path: str | Path) -> dict:
       """Load and validate transform.json against the schema. Returns parsed dict."""
   ```

   Implementation notes:
   - Use `numpy` only; no third-party CRS libs (no `pyproj`, no `geopandas`). The origin
     decision is documented in metadata, but the math is plain subtraction.
   - `utm_to_local` and `local_to_utm` accept either `(N, 3)` arrays or `(3,)` single
     points; broadcast accordingly.
   - `write_transform_json` validates against the schema before writing (using the
     standard library `jsonschema` if available, else a minimal hand-coded check that
     verifies the required keys are present with correct types).

2. **`transform.schema.json` (JSON Schema Draft 7).**
   ```json
   {
     "$schema": "http://json-schema.org/draft-07/schema#",
     "title": "SAFS mesh local-frame transform",
     "type": "object",
     "required": ["version", "utm_zone", "datum",
                  "origin_utm_m", "convention", "comment"],
     "properties": {
       "version":   { "const": 1 },
       "utm_zone":  { "type": "string", "pattern": "^[0-9]{1,2}[NS]$" },
       "datum":     { "type": "string" },
       "origin_utm_m": {
         "type": "array",
         "items": { "type": "number" },
         "minItems": 3, "maxItems": 3,
         "description": "[E0, N0, Z0] in metres. v_local = v_utm − origin_utm_m."
       },
       "convention": {
         "type": "object",
         "required": ["x_axis", "y_axis", "z_axis", "z_sign"],
         "properties": {
           "x_axis":  { "const": "UTM_easting" },
           "y_axis":  { "const": "UTM_northing" },
           "z_axis":  { "const": "elevation" },
           "z_sign":  { "const": "positive_up" }
         }
       },
       "comment":   { "type": "string" },
       "generated_by": { "type": "string" },
       "git_commit":   { "type": "string" }
     }
   }
   ```

3. **`transform.json` (concrete example, written during Phase 2).**
   ```json
   {
     "version": 1,
     "utm_zone": "11N",
     "datum": "NAD83",
     "origin_utm_m": [500000.0, 3765000.0, 0.0],
     "convention": {
       "x_axis": "UTM_easting",
       "y_axis": "UTM_northing",
       "z_axis": "elevation",
       "z_sign": "positive_up"
     },
     "comment": "v_local = v_utm − origin_utm_m. Right-handed (E,N,Up) metres. Free surface at z=0.",
     "generated_by": "ts_to_stl.py",
     "git_commit": "<filled at runtime>"
   }
   ```

### Interfaces
- Public Python: `utm_to_local`, `local_to_utm`, `write_transform_json`,
  `read_transform_json`, plus the `ORIGIN_UTM` constant.
- File contract: `mesh/transform.json` validated by `mesh/transform.schema.json`.

### Edge Cases to Handle
- User passes `--origin-easting 0` (no translation, raw UTM): the helpers must work
  identically; `transform.json` records `origin_utm_m: [0, 0, 0]`. The numerical
  conditioning warning above no longer applies — emit a warning if any coord exceeds
  10^6 m.
- `transform.json` exists but predates schema v1 (no `version` field): refuse to load,
  print a migration message rather than silently using a stale convention.
- `read_transform_json` called when the file does not exist: raise `FileNotFoundError`
  with a message pointing the user to run `ts_to_stl.py` first.

### Acceptance Criteria
- [ ] `python -c "from safs_origin import utm_to_local, ORIGIN_UTM;
      import numpy as np;
      v = np.array([[500000.0, 3765000.0, 0.0]]);
      print(utm_to_local(v))"` prints `[[0. 0. 0.]]`.
- [ ] Round-trip identity: for `v` random in `R^3`, `local_to_utm(utm_to_local(v)) == v`
      bit-exact (subtraction is its own inverse for IEEE-754 floats *when no overflow*,
      which holds at our magnitudes).
- [ ] `write_transform_json(...)` produces a file that passes `read_transform_json(...)`
      schema validation.
- [ ] No script in the repo other than `safs_origin.py` contains the literal `500000`
      or `3765000` (grep gate; values must come from the module).

### Dependencies
- Depends on: nothing.
- Required by: `PLAN.md` Phases 2, 3, 5, 6.

---

## Phase 2: Apply the Mapping in the TSurf → STL Pipeline

### Goal
After this phase, `mesh/ts_to_stl.py` (the script defined in `PLAN.md` Phase 2) writes
STLs whose vertex coordinates are in the local frame defined above, and writes
`mesh/transform.json` documenting the choice.

### Files to Create
None.

### Files to Modify
- `mesh/ts_to_stl.py` — replace the centroid-based translation with a call to
  `safs_origin.utm_to_local(...)`.

### Detailed Requirements

1. **Replace centroid logic.** Inside `ts_to_stl.py`, immediately after parsing each
   `.ts` file's `V_utm` array, call:
   ```python
   from safs_origin import utm_to_local, write_transform_json, ORIGIN_UTM
   V_local = utm_to_local(V_utm)   # uses default origin
   ```
   No centroid is computed. No per-fault translation. The same `ORIGIN_UTM` is used for
   all 8 faults.

2. **Write `transform.json` once per run.** At the end of `ts_to_stl.py` (after all 8
   STLs are written), call:
   ```python
   write_transform_json(
       out_dir / "transform.json",
       origin=cli_origin,                       # default is ORIGIN_UTM, overridable
       extra={"generated_by": "ts_to_stl.py",
              "git_commit": _current_git_commit(),
              "input_resolution_m": args.res,
              "input_files": [str(p) for p in input_paths]}
   )
   ```
   The `extra` dict is merged into the JSON before schema validation.

3. **CLI flags.** `ts_to_stl.py` gains:
   ```
   --origin-easting   FLOAT   default 500000.0   (UTM E0)
   --origin-northing  FLOAT   default 3765000.0  (UTM N0)
   --origin-elevation FLOAT   default 0.0        (Z0)
   ```
   When any of these is passed, the override propagates into both `utm_to_local()` and
   `write_transform_json()`.

4. **Free-surface clamp ordering.** The flow is:
   ```
   parse → utm_to_local → clamp z>0 to 0 → collapse 1mm duplicates → drop zero-area
   ```
   The clamp happens in the local frame (z is unchanged by the translation since Z0=0
   by default; if user overrides `Z0`, the clamp threshold becomes `−Z0` in the local
   frame — adjust accordingly: clamp to `−Z0` not `0`).

### Interfaces
None new — extends an interface defined in `PLAN.md`.

### Edge Cases to Handle
- User passes `--origin-elevation 100` (a non-zero Z0): vertices originally at z = 50 m
  elevation are now at z = −50 m local; clamp threshold becomes `0 − Z0 = −100 m`.
  Implementation: `V_local[V_local[:,2] > -Z0, 2] = -Z0`. Document this in the script's
  `--help`. (For default `Z0 = 0`, clamp is the familiar `> 0 → 0`.)
- A vertex is exactly on the Z0 plane (z = 0 in UTM): no clamp applied (strict `>`).
- The git commit hash subprocess call fails (running outside a git checkout): record
  `"git_commit": "unknown"` and continue.

### Acceptance Criteria
- [ ] Running `python ts_to_stl.py --res 2000` (defaults) produces STLs whose union
      bounding box is `X ∈ [−134 946, +122 329], Y ∈ [−72 722, +73 998],
      Z ∈ [−18 014, 0]` (Z upper clamped to 0). Tolerance: ±1 m on each bound.
- [ ] `transform.json` validates against `transform.schema.json`.
- [ ] Re-running `ts_to_stl.py` with `--origin-easting 0 --origin-northing 0` produces
      STLs whose union bounding box matches the *original UTM* bounds (i.e., the
      translation is the only difference between the two outputs).
- [ ] `local_to_utm(V_local) == V_utm` for every vertex of every STL, to within 1e−6 m
      (round-trip test in `validate_msh.py`).

### Dependencies
- Depends on: Phase 1 (this document); `PLAN.md` Phase 1.
- Required by: `PLAN.md` Phase 3.

---

## Phase 3: Thread the Origin Through Gmsh

### Goal
After this phase the bulk box and physical-group bounding-box queries in `safs.geo` are
expressed in the local frame, with no UTM literals anywhere in the `.geo` source.

### Files to Create
None.

### Files to Modify
- `mesh/safs.geo` — bulk box uses local-frame extents (already the plan in `PLAN.md`
  Phase 3, but make explicit).
- `mesh/generate_safs_mesh.py` — reads `transform.json`, computes local-frame bounding
  boxes, passes them to Gmsh via `-setnumber`.

### Detailed Requirements

1. **`generate_safs_mesh.py` flow.**
   ```python
   tx = read_transform_json("mesh/transform.json")
   E0, N0, Z0 = tx["origin_utm_m"]                # local origin in UTM
   per_fault_bbox = read_bbox_json("mesh/bbox.json")  # already in local frame, written by ts_to_stl
   x_min = min(b["xmin"] for b in per_fault_bbox.values())
   ...
   subprocess.run(["gmsh", "-3", "mesh/safs.geo",
                   "-setnumber", "x_min", str(x_min),
                   "-setnumber", "x_max", str(x_max),
                   "-setnumber", "y_min", str(y_min),
                   "-setnumber", "y_max", str(y_max),
                   "-setnumber", "buf_x", str(args.buf_xy),
                   "-setnumber", "buf_y", str(args.buf_xy),
                   "-setnumber", "z_top", "0.0",
                   "-setnumber", "z_bot", str(args.z_bot),
                   "-setnumber", "res_f", str(args.res),
                   "-o", args.output])
   ```
   `bbox.json` is a per-fault bounding-box file written by `ts_to_stl.py` *in the local
   frame* (this is a small addition to `PLAN.md` Phase 2, listed there only as
   `output/cleanup_log_<res>.csv` — replace or supplement with `bbox.json`).

2. **`safs.geo` constants block (no UTM literals).**
   ```
   DefineConstant[ x_min  = {-150000, Name "Fault bbox xmin (m, local)"} ];
   DefineConstant[ x_max  = { 150000, Name "Fault bbox xmax (m, local)"} ];
   DefineConstant[ y_min  = { -80000, Name "Fault bbox ymin (m, local)"} ];
   DefineConstant[ y_max  = {  80000, Name "Fault bbox ymax (m, local)"} ];
   DefineConstant[ buf_x  = {  30000, Name "Buffer along ±X (m)"} ];
   DefineConstant[ buf_y  = {  30000, Name "Buffer along ±Y (m)"} ];
   DefineConstant[ z_top  = {      0, Name "Free surface elevation (m, local)"} ];
   DefineConstant[ z_bot  = { -40000, Name "Domain bottom (m, local)"} ];
   ```
   Bulk box:
   ```
   Box(1) = { x_min - buf_x, y_min - buf_y, z_bot,
              (x_max - x_min) + 2*buf_x,
              (y_max - y_min) + 2*buf_y,
              z_top - z_bot };
   ```
   The defaults are deliberately wide enough to mesh even without `-setnumber` overrides
   (so a developer can `gmsh safs.geo` interactively without running the driver), but
   the driver always sets explicit values from `bbox.json`.

3. **No UTM in the `.geo` source, ever.** This is enforced by a grep gate in CI:
   ```
   ! grep -E "(500000|3765000|36[0-9]{4}\.|37[0-9]{4}\.)" miniapps/seas/safs/mesh/safs.geo
   ```
   If a UTM literal slips into `safs.geo`, the test fails. (The pattern catches both the
   origin literals and any 6-digit UTM-magnitude number.)

4. **Gmsh's coordinate system is the local frame.** When the user opens
   `output/safs_2000m.msh` in Gmsh GUI, vertex 1 reads as `(x_local, y_local, z_local)`,
   not UTM. To re-georeference to UTM, downstream tools (paraview scripts, Python
   post-processors) load `transform.json` and call `local_to_utm`.

### Interfaces
- `generate_safs_mesh.py --origin-easting / --origin-northing / --origin-elevation` flags
  exist for symmetry with `ts_to_stl.py`, but are **not actually used by Gmsh** — the
  driver only forwards the local-frame bounding box. The flags exist solely so that if
  the user passes mismatched origins between `ts_to_stl.py` and `generate_safs_mesh.py`,
  the driver detects the mismatch (against `transform.json`) and aborts with a clear
  error.

### Edge Cases to Handle
- `transform.json` and CLI flags disagree → driver aborts with `"transform.json says
  origin = (E0, N0, Z0) but CLI flag says (E0', N0', Z0'); refusing to mesh against
  inconsistent metadata."`.
- `bbox.json` per-fault values exceed the `safs.geo` `DefineConstant` defaults → no
  problem, the driver overrides via `-setnumber`. But document that running
  `gmsh safs.geo` *without* the driver will use the conservative defaults and may not
  capture all faults — `ts_to_stl.py` warns if any fault's local bbox exceeds the
  defaults.

### Acceptance Criteria
- [ ] `grep -E "500000|3765000" miniapps/seas/safs/mesh/safs.geo` returns nothing.
- [ ] `python generate_safs_mesh.py --res 2000` and `gmsh safs.geo` (interactive,
      defaults) both produce a valid mesh that contains all 8 fault surfaces.
- [ ] Driver detects and rejects mismatched origin between flags and `transform.json`.

### Dependencies
- Depends on: Phase 2, `PLAN.md` Phase 3.
- Required by: `PLAN.md` Phase 5.

---

## Phase 4: Validation, Inverse Mapping, and Visualization Aids

### Goal
After this phase, `validate_msh.py` (defined in `PLAN.md` Phase 6) verifies that the
produced mesh is *re-georeferenceable* — i.e., applying `local_to_utm` to its vertices
recovers UTM coordinates that lie inside the union of the input `.ts` bounding boxes.

### Files to Create
- `mesh/visualize_traces_on_map.py` — optional visualization script that plots fault
  surface traces (z = 0 polylines) on a Southern California base map, using
  `transform.json` to convert back to UTM and then to lat/lon for the plot library.

### Files to Modify
- `mesh/validate_msh.py` — add a section that round-trips coordinates through
  `transform.json`.

### Detailed Requirements

1. **Round-trip recovery test in `validate_msh.py`.**
   ```python
   from safs_origin import read_transform_json, local_to_utm
   tx = read_transform_json("mesh/transform.json")
   verts_local = mesh.points  # (N, 3)
   verts_utm = local_to_utm(verts_local, origin=tuple(tx["origin_utm_m"]))
   assert verts_utm[:,0].min() >= MIN_E_OBSERVED - 100   # tolerance for buffer
   assert verts_utm[:,0].max() <= MAX_E_OBSERVED + buf_x + 100
   ...
   ```
   The reference bounds `MIN_E_OBSERVED` etc. come from the per-fault bbox audit done in
   `PLAN.md` Phase 1.

2. **Optional visualization.** `visualize_traces_on_map.py` extracts the z=0 trace of
   each fault (intersection of fault physical surface with the free-surface physical
   surface), converts to UTM via `local_to_utm`, then to lat/lon via `pyproj`
   (Zone 11N → EPSG:4326), and plots on top of a `cartopy` map of Southern California.
   The output PNG should visually match the red lines in the original SCEC CFM viewer
   image. This is a sanity check for the human, not an automated test.

3. **Document recovery in `safs/README.md`.** Add a one-paragraph "How to re-georeference
   the mesh" section that walks through:
   ```python
   from safs_origin import read_transform_json, local_to_utm
   import meshio
   m = meshio.read("output/safs_2000m.msh")
   tx = read_transform_json("mesh/transform.json")
   m.points = local_to_utm(m.points, origin=tuple(tx["origin_utm_m"]))
   meshio.write("output/safs_2000m_utm.msh", m)
   ```

### Interfaces
None new.

### Edge Cases to Handle
- Recovered UTM coordinates lie outside the per-fault bbox by more than the buffer
  width: indicates a bug in the translation. Validation script fails with the offending
  vertex index and coordinates.
- `transform.json` missing when `validate_msh.py` runs: degrade to "local-frame only"
  validation (do not run round-trip), warn loudly. This catches the case where the user
  copied the mesh but forgot the metadata.

### Acceptance Criteria
- [ ] On a freshly generated mesh, `validate_msh.py` round-trip recovery succeeds for
      every vertex.
- [ ] `visualize_traces_on_map.py` produces a PNG whose fault traces are visually
      consistent with the original SCEC CFM map (manual check, no automated tolerance).

### Dependencies
- Depends on: `PLAN.md` Phase 5; Phase 1, 2, 3 of this document.
- Required by: nothing.

---

## Testing Strategy

| Test | Where | What it proves |
| ---- | ----- | -------------- |
| `safs_origin.py` round-trip identity | unit test (`pytest`) | `local_to_utm(utm_to_local(v)) == v` exactly for representative magnitudes |
| `safs_origin.py` no-overflow at UTM magnitudes | unit test | `(v − origin) + origin == v` to ≤ 1e−6 m for v ∈ [3.5e5, 6.5e5] × [3.7e6, 3.85e6] |
| `transform.json` schema validation | unit test | invalid JSON is rejected with a useful message |
| `ts_to_stl.py` STL bbox in local frame | integration test | matches the expected post-translation bounds listed in this plan |
| `validate_msh.py` round-trip on `.msh` | integration test | mesh vertices round-trip to UTM within 1 m |
| Grep gate: no UTM literals in `.geo` | CI script | enforces single-source-of-truth |
| Visual map check | manual | traces match SCEC CFM viewer |

## Risk Assessment

1. **Origin choice locks downstream consumers.** Once meshes are generated and shared,
   changing `E0` or `N0` invalidates every cached mesh, every saved ParaView state, every
   plotting script. We mitigate by (a) versioning `transform.json` (`"version": 1`), (b)
   making the override mechanism work but logged, and (c) discouraging override in
   `safs/README.md`. If the team ever decides to change the origin, the change is a
   coordinated migration, not a quiet patch.

2. **Project's BP5 mesh uses a different convention** (X = along-strike, Y = fault-normal).
   SAFS uses (X = UTM E, Y = UTM N). Code in `fault/`, `domain/`, and `solver/` that
   currently assumes BP5's mesh convention will need review when SAFS meshes are first
   consumed (but that is solver work, out of scope here). The risk specific to this plan:
   a developer might assume the SAFS mesh uses BP5's convention and produce wrong fault
   normals. Mitigation: (a) `safs/README.md` explicitly states the convention, (b)
   `validate_msh.py` prints the convention in its summary output.

3. **UTM Zone 11N assumption.** All CFM data inspected so far falls inside Zone 11N. If
   a future SAFS-area fault outside Zone 11N is added (none of the eight selected do
   this), the UTM translation no longer makes geographic sense — distances would be
   distorted. Mitigation: `audit_ts_quality.py` already records bounding boxes; an
   easting < 200 km or > 800 km would indicate the file is in a different zone; flag
   such files at audit time.

4. **Origin literal proliferation.** A future contributor might inline `500000` somewhere
   (e.g., a one-off conversion script). Mitigation: the grep gate (acceptance criterion
   in Phase 1) makes this a CI failure, forcing them to import from `safs_origin.py`.
