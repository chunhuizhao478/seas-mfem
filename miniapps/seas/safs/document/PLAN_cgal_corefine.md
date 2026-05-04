# Implementation Plan: CGAL-based Multi-Fault Corefinement

## Overview

Replace the Python Phase-1 (`fault_intersect.py`) and Phase-2
(`conformalize_faults.py`) of `PLAN_multifault_intersections.md` with
a single C++ binary that calls
`CGAL::Polygon_mesh_processing::corefine` followed by
`CGAL::Polygon_mesh_processing::isotropic_remeshing`. The new tool
emits per-fault conformal STLs plus the diagnostic JSON sidecars
the Python conformalizer emits (different schema — see §Schema
delta in Constraints), so steps 4–7 of the existing pipeline
(`generate_safs_mesh.py`, `write_fault_provenance.py`,
`validate_msh.py`, `convert_msh.py`) are untouched. The 2-fault
Mill Creek × SBMT-SAF target at 2000 m must pass **at least 9 of
11 validator checks** (the must-pass set: 1–6 + 8–10), with
checks 7 and 11 as stretch (size-field-dependent — see Phase 3
acceptance for the must-pass / stretch split, P-003). The
existing Python pipeline stays in the tree as a legacy
regression suite until the CGAL path is the proven default.

## Status

- **Created:** 2026-04-29.
- **Driving evidence:** `miniapps/seas/safs/REVIEW.md` finding R-010
  ("Replace ~3000 LOC custom conformalizer with CGAL corefine +
  isotropic_remeshing"). The current Python pipeline produces
  10 T-junctions, 22 polyline-coincident slivers, and forces
  `algo3d=4` which silently leaves 0 of 36505 fault triangles
  shared by two tets (validator check 5 fails 100%).
- **Companion plans:**
  - `PLAN_multifault_intersections.md` — the validation gates and
    JSON schemas the new path must inherit.
  - `PLAN.md` (overall), `PLAN_origin.md`, `PLAN_domain.md`,
    `PLAN_smoke_millcreek.md` — unchanged.
- **First milestone:** `run_two_crossing_2000m_cgal.sh` passes
  the must-pass set (≥ 9/11 validator checks, including the
  primary `internal_interface` and `fault_edge_length` checks)
  on Mill Creek × SBMT-SAF at 2000 m.  11/11 is a stretch goal
  conditional on size-field tuning (P-003).

## Constraints

### Interface constraints (cannot change)

- **JSON schemas** the new tool emits.  Verified against actual
  consumers (P-005):

  - `triangle_to_fault.json` is a **presence-only marker**.
    `generate_safs_mesh.py:445` only calls `.exists()` on it;
    `write_fault_provenance.py` does not read it at all (it uses
    KDTree centroid matching against per-fault STLs).  We still
    emit it for diagnostic / human-readability continuity with
    the Python path:
    ```json
    {
      "schema_version": 1,
      "faults": {
        "<short_name>": {"n_triangles": <int>, "range": [<int>, <int>]}
      },
      "n_total_triangles": <int>
    }
    ```
    Constraint: ranges must be **contiguous** and partition
    `[0, n_total_triangles)` in the order `included` is given.
    No downstream tool reads these values; correctness is
    diagnostic only.

  - `intersection_report.json` is a **diagnostic-only** report
    consumed by humans.  It is NOT consumed by any of
    `generate_safs_mesh.py`, `write_fault_provenance.py`,
    `validate_msh.py`, or `convert_msh.py`.  The implementer is
    free to add or omit keys with the constraint that the
    top-level `per_fault` and `pairs` keys, plus per-pair
    `gates.{manifold_A, manifold_B, polyline_edge_coincidence,
    interior_crossing_only}` strings, are preserved (so existing
    dashboards / tail-of-log greps continue to work).  Reference
    schema:
    ```json
    {
      "schema_version": 1,
      "backend": "cgal",
      "target_edge_length_m": <float>,
      "snap_m": <float>,
      "clearance_m": <float>,
      "pairs": {
        "<short_a>__x__<short_b>": {
          "n_constrained_edges_A": <int>,
          "n_constrained_edges_B": <int>,
          "pre_split_n_tri_A": <int>,
          "post_split_n_tri_A": <int>,
          "pre_split_n_tri_B": <int>,
          "post_split_n_tri_B": <int>,
          "remesh_iters": <int>,
          "gates": {
            "manifold_A": "PASS",
            "manifold_B": "PASS",
            "polyline_edge_coincidence": "PASS",
            "interior_crossing_only": "PASS"
          }
        }
      },
      "per_fault": {
        "<short>": {"n_vertices": <int>, "n_triangles": <int>, "stl_path": <str>}
      }
    }
    ```
    With CGAL `corefine` + `protect_constraints(true)` + the
    Phase-2 `split_long_edges` pre-pass (P-001), T-junctions
    cannot be introduced; the `manifold_X` gates always emit
    `"PASS"` (never `"PASS_WITH_<n>_T_JUNCTIONS"`).  If the
    code ever emits the latter, it is a regression.

#### Schema delta vs Python output (deliberate; P-007)

The CGAL `intersection_report.json` is intentionally not byte-
equivalent to the Python output.  This table documents every
key that differs so operators diffing the two reports do not
waste time reconciling them.

| Key                                       | Python | CGAL   | Notes                                            |
|---|---|---|---|
| `schema_version`                          | 1      | 1      | bump to 2 in a future ticket if needed           |
| `backend`                                 | absent | `"cgal"` | NEW — disambiguates output                     |
| `snap_m`                                  | yes    | yes    | unchanged semantics                              |
| `clearance_m`                             | absent | yes    | NEW — record clamp value                         |
| `target_edge_length_m`                    | yes    | yes    | unchanged                                        |
| `drop_short_segment_frac`                 | yes    | absent | n/a — CGAL handles via remeshing                 |
| `min_pierce_separation_m`                 | yes    | absent | n/a — corefine produces canonical pierces        |
| `cdt_min_angle_deg`                       | yes    | absent | n/a — no CDT step                                |
| `pairs.<x>.n_polylines`                   | yes    | absent | replaced by `n_constrained_edges_*`              |
| `pairs.<x>.n_polyline_vertices`           | yes    | absent | replaced by `n_constrained_edges_*`              |
| `pairs.<x>.n_constrained_edges_{A,B}`     | absent | yes    | NEW — direct CGAL count                          |
| `pairs.<x>.smoothing_iters`               | yes    | absent | replaced by `remesh_iters`                       |
| `pairs.<x>.remesh_iters`                  | absent | yes    | NEW — CGAL remesh iteration count                |
| `pairs.<x>.n_dropped_short_{A,B}`         | yes    | absent | n/a — corefine doesn't drop                      |
| `gates.manifold_{A,B}` value alphabet     | `"PASS"` ∪ `"PASS_WITH_<n>_T_JUNCTIONS"` | `"PASS"` only | CGAL never emits T-junctions |
| `gates.{polyline_edge_coincidence,interior_crossing_only}` | yes | yes | unchanged |
| `per_fault.<short>.{n_vertices,n_triangles,stl_path}` | yes | yes | unchanged |
- **Per-fault STL output** must be ASCII, `solid SAFS:<short>:conformal`
  / `endsolid` framing, full `%+.17e` coordinates (matches
  `conformalize_faults.py:_write_ascii_stl` line 1255 — required so
  `_combine_stls` snap dedup sees no float32 truncation).
- **Per-fault output filename**: `<out_stl_dir>/<short>.stl`,
  exactly as today.
- **Coordinate frame**: SAFS local-Cartesian metres (per
  `safs_origin.py`). The CGAL tool reads in this frame and writes
  in this frame — no UTM round-trip.
- **No edits to** `bp5/`, `bp1/`, `bp2/`, `domain/`, `fault/`,
  `solver/`, or `friction/dieterich_ruina.hpp` (project memory
  `[C2]`).
- **No edits to** `generate_safs_mesh.py`, `write_fault_provenance.py`,
  `validate_msh.py`, `convert_msh.py` for Phases 0–3. Phase 4 may
  bump `--algo3d` defaults only.
- **No edits to** the Python conformalizer or its tests for
  Phases 0–3. The Python path remains the fallback.

### Dependency constraints

- **CGAL ≥ 5.0** (header-only). Requires Boost headers, GMP, MPFR.
  Build env on laptop: install via
  `conda install -c conda-forge cgal boost-cpp gmp mpfr` into the
  existing `mfem-dev` conda env (which already has gmp + mpfr per
  `conda list -n mfem-dev`).
- **C++17** standard (`std::string_view`, `std::filesystem`).
- **CMake ≥ 3.18** (matches what MFEM-dev already needs).
- **Frontera-style HPC build**: defer to Phase 5. Per project
  memory `[Frontera approval]`, no `sbatch` / `ssh` / allocation
  check until the user explicitly approves; the plan documents
  module-list expectations only.

### Convention constraints

- **C++ source layout**: new directory
  `miniapps/seas/safs/tools/`. CMake target only — no auto-build
  by the existing MFEM Makefile. Per project conventions, isolated
  tools live under `tools/`, parallel to `mesh/`.
- **CLI flag conventions**: long-form flags
  (`--in-stl-dir`, `--out-stl-dir`, `--include-fault`,
  `--target-edge-m`, `--clearance-m`, `--snap-m`,
  `--remesh-iters`) — mirror the Python conformalizer's CLI so
  the shell scripts can be near-copies.
- **Logging**: stderr text, prefixed with `[corefine_faults] `.
  No structured logging; the JSON report is the machine-readable
  output.
- **Exit codes**: 0 success; 1 user error (bad flags / missing
  files); 2 CGAL runtime error (corefine failed / mesh became
  non-manifold); 3 schema-validation error.

### Numerical constraints

- **Free-surface clearance**: input STLs from `ts_to_stl.py` are
  already clamped to `z ≤ -clearance_m`. The CGAL tool must
  preserve this invariant — no remeshing operation may pull a
  vertex above `z = -clearance_m`. Implemented by adding a
  hard-ceiling smoothing-vertex projection (Phase 2 §requirement
  R2-7).
- **Polyline edge preservation**: every constrained edge marked by
  `corefine` must remain in both meshes after `isotropic_remeshing`.
  The remesher's `edge_is_constrained_map` parameter guarantees
  this; the tool asserts the count post-remesh equals the count
  post-corefine.
- **Snap tolerance for STL writeback**: write coordinates at full
  `%+.17e` (17-decimal-digit) precision, matching the Python tool.
  `_combine_stls` will then dedup at 1 mm (per REVIEW.md R-005);
  CGAL's exact-predicates corefine guarantees that polyline
  vertices on fault A and fault B are bit-identical, so the
  combine-step dedup becomes trivial.
- **No `EPECK` for Phase 1–3**: use
  `Exact_predicates_inexact_constructions_kernel` (EPICK).
  Per CGAL docs and REVIEW.md R-010, EPICK is sufficient for
  single corefinement + remeshing; EPECK is only needed for
  *consecutive* boolean operations. Document this choice in the
  source header and in `PLAN_cgal_corefine.md` so a future
  reviewer knows when to upgrade.

---

## Phase 0 — Build environment and CGAL availability check

### Goal

After this phase, the developer can compile a "hello world" CGAL
program in the laptop conda env without modifying the project.
The plan documents the same recipe for Frontera so Phase 5 can
proceed without surprise.

### Files to create

- `miniapps/seas/safs/tools/CMakeLists.txt` — top-level CMake for
  the `corefine_faults` target. Phase 0 lands a stub that compiles
  an empty `main()`; Phase 1 fills it in. Required so subsequent
  phases have a build system to extend.
- `miniapps/seas/safs/tools/README.md` — environment recipes for
  laptop (conda) and Frontera (modules + manual CGAL fetch).
  Approx 120 lines, structured as:
  1. **Prerequisites** — CGAL 5.6.1 (pinned; tested with the
     `corefinement_difference_remeshed` example), GMP, MPFR,
     Boost ≥ 1.66, CMake ≥ 3.18, C++17 toolchain.
  2. **Laptop recipe** —
     ```bash
     conda activate mfem-dev
     conda install -c conda-forge cgal=5.6.1 boost-cpp gmp mpfr
     cmake -S miniapps/seas/safs/tools -B miniapps/seas/safs/tools/build \
           -DCMAKE_BUILD_TYPE=Release
     cmake --build miniapps/seas/safs/tools/build
     ```
  3. **Frontera recipe** (DO NOT EXECUTE in Phase 0 — per project
     memory `feedback_frontera_approval.md`, Frontera commands
     require explicit user approval).  Document the starting
     point that the user can run when ready, P-010:

     a. **Module list** — copy from a known-working sbatch
        verbatim per project memory `feedback_sbatch_modules.md`.
        The `seas_tpv102` driver currently uses:
        ```bash
        module reset
        module load intel/19.1.1 impi/19.0.9 boost/1.71 \
                    gcc/9.1.0 cmake/3.24.2 fftw3/3.3.10 \
                    gmp mpfr
        export LD_LIBRARY_PATH=$TACC_GMP_LIB:$TACC_MPFR_LIB:$LD_LIBRARY_PATH
        ```
        Verify the exact module versions by reading the
        `seas_tpv102` sbatch when the user authorises a login
        session.  Do NOT trim modules even if they look unrelated
        (`fftw3` in particular has been the cause of past
        link failures per `feedback_sbatch_modules.md`).

     b. **CGAL fetch** (CGAL is header-only since 5.0, so no
        compilation step is needed):
        ```bash
        cd $WORK
        git clone --depth 1 -b v5.6.1 \
            https://github.com/CGAL/cgal.git cgal-5.6.1
        export CGAL_DIR=$WORK/cgal-5.6.1/lib/cmake/CGAL
        ```
        Pin **5.6.1**: tested with the
        `corefinement_difference_remeshed` example referenced
        in REVIEW.md R-010.  Versions 5.4 and 5.5 also work;
        CGAL 6.0 introduces ABI changes and is NOT on the
        tested list.

     c. **CMake build** with explicit C++17 `<filesystem>` link
        for older toolchains:
        ```cmake
        if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
          target_link_libraries(corefine_faults PRIVATE stdc++fs)
        endif()
        ```
        Required only on Frontera's gcc 9.1 / intel 19.1.1
        because pre-9.0 gcc puts `<filesystem>` symbols in a
        separate library.  No-op on macOS clang and conda-forge
        gcc 13+.

     d. **Acceptance**: the `corefine_faults_hello` smoke binary
        builds and prints the CGAL version on a Frontera login
        node before any compute-node work begins.

     e. **Failure mode**: if the build dies on missing fftw3 /
        `LD_LIBRARY_PATH` / Boost — copy the module list
        VERBATIM from the seas_tpv102 sbatch (per project memory
        `feedback_sbatch_modules.md`); do NOT trim modules even
        if they look unrelated.

  4. Smoke build target: `cmake --build tools/build --target corefine_faults_hello`.

### Files to modify

- None.

### Detailed requirements

1. **CMakeLists.txt** must:
   - `cmake_minimum_required(VERSION 3.18)`.
   - `project(safs_tools CXX)` with `CMAKE_CXX_STANDARD 17`,
     `CMAKE_CXX_STANDARD_REQUIRED ON`.
   - **Pin CGAL minimum + tested-against version** (Unreviewed
     Areas concern from REVIEW_plan_cgal_corefine.md):
     ```cmake
     find_package(CGAL 5.4 REQUIRED COMPONENTS Core)
     if(CGAL_VERSION VERSION_GREATER_EQUAL 6.0)
       message(FATAL_ERROR
         "CGAL ${CGAL_VERSION} is past the tested ABI window "
         "[5.4, 6.0). The tool was validated against CGAL 5.6.1; "
         "later versions may have changed `protect_constraints` "
         "semantics for isotropic_remeshing.  See "
         "PLAN_cgal_corefine.md §Risk assessment §1 + tools/README.md.")
     endif()
     ```
     The minimum is 5.4 (first version where the
     `protect_constraints(true)` semantics on `isotropic_remeshing`
     are documented as exact preservation); the tested version
     is 5.6.1; 6.0 is explicitly out of bounds because CGAL
     bumps ABI between major versions.
   - `find_package(Boost 1.66 REQUIRED)`.
   - Single executable target `corefine_faults_hello` with one
     `main.cpp` source that includes
     `<CGAL/Polygon_mesh_processing/corefinement.h>`,
     `<CGAL/Surface_mesh.h>`, prints CGAL version, exits 0.
   - Compile with `-O2` in Release.
   - `<filesystem>` link guard for old toolchains:
     ```cmake
     if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
        AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
       target_link_libraries(corefine_faults_hello PRIVATE stdc++fs)
     endif()
     ```
2. **README.md** must include the **exact `conda install`
   command** that will be run, plus an `ldd`/`otool -L` snippet
   showing which CGAL/GMP/MPFR/Boost libraries the produced binary
   links against (or, since CGAL is header-only, which it does
   *not* link to Boost dynamically — Boost on macOS conda is
   header-only-usable too).
3. The CMake target must build cleanly when invoked from
   `tools/build/` with `mfem-dev` activated. No system-wide CGAL
   install required.

### Interfaces

- New executable `tools/build/corefine_faults_hello`. No CLI
  contract beyond "exits 0 and prints a version line".

### Edge cases to handle

- **CGAL not installed**: `find_package(CGAL REQUIRED)` fails
  cleanly with a CMake error. Document in README how to install.
- **CGAL version out of bounds**: enforced by the
  `find_package(CGAL 5.4 REQUIRED)` + `VERSION_GREATER_EQUAL 6.0`
  guard above.  Tested window: [5.4, 6.0); pin 5.6.1 in
  conda + Frontera installs.
- **Boost not found / too old**: enforced by
  `find_package(Boost 1.66 REQUIRED)`.  Frontera's `boost/1.71`
  satisfies; conda-forge's default satisfies.

### Acceptance criteria

- [ ] `cmake -S miniapps/seas/safs/tools -B miniapps/seas/safs/tools/build` succeeds in the `mfem-dev` conda env on the user's laptop.
- [ ] `cmake --build miniapps/seas/safs/tools/build` succeeds.
- [ ] `./miniapps/seas/safs/tools/build/corefine_faults_hello` exits 0 and prints `CGAL version: <X.Y>`.
- [ ] The Python pipeline still runs unchanged (no regressions in `pytest miniapps/seas/safs/mesh/tests/`).
- [ ] Frontera build recipe is **documented** in `tools/README.md` but not yet executed.

### Dependencies

- Depends on: nothing.
- Required by: Phase 1.

---

## Phase 1 — Minimum viable `corefine_faults` tool

### Goal

After this phase, given the same inputs `conformalize_faults.py`
takes today (per-fault raw STLs), `corefine_faults` produces
schema-compatible per-fault conformal STLs and JSON sidecars.
The tool runs `corefine` only — no remeshing yet — so the output
will have slivers, but it will be a valid conforming PLC that
`generate_safs_mesh.py` can read.

### Files to create

- `miniapps/seas/safs/tools/corefine_faults/main.cpp` — driver
  (~150 LOC): CLI parsing, IO, orchestration of corefine over
  all fault pairs, JSON emission.
- `miniapps/seas/safs/tools/corefine_faults/io.hpp` /
  `io.cpp` — STL read/write (ASCII), JSON write (no external JSON
  lib; hand-rolled because the schema is fixed and small).
- `miniapps/seas/safs/tools/corefine_faults/corefine.hpp` /
  `corefine.cpp` — wrapper around `PMP::corefine` for one
  unordered fault pair, with the constrained-edge property maps.
- `miniapps/seas/safs/tools/tests/test_corefine_smoke.cpp` —
  hand-rolled C++ smoke test (no GoogleTest dep): two
  perpendicular unit squares, run corefine, assert the
  intersection edge appears in both meshes. Uses CMake's
  `add_test` + `ctest` for invocation.

### Files to modify

- `miniapps/seas/safs/tools/CMakeLists.txt` — add the
  `corefine_faults` target with the four source files; add the
  `test_corefine_smoke` test target wired into `enable_testing()`
  and `add_test`.

### Detailed requirements

1. **CLI interface** (`main.cpp`):
   ```cpp
   struct Args {
     std::filesystem::path in_stl_dir;       // required
     std::filesystem::path out_stl_dir;      // required
     std::vector<std::string> included;      // 1+ --include-fault
     double snap_m       = 1e-3;             // --snap-m
     double clearance_m  = 0.0;              // --clearance-m
     double target_edge_m = 0.0;             // --target-edge-m (Phase 2)
     int    remesh_iters  = 0;               // Phase 2
     bool   verbose       = false;
   };
   int parse_args(int argc, char** argv, Args& out);
   int main(int argc, char** argv);
   ```
   - Phase 1 ignores `--target-edge-m` and `--remesh-iters` (or
     fails if non-zero, with a "Phase 2 — not yet implemented"
     message).
   - `--include-fault <short>` repeatable. Order matters
     (determines triangle-range ordering in the output JSON).
2. **STL I/O** (`io.cpp`):
   ```cpp
   namespace safs::io {
   bool read_ascii_stl(const std::filesystem::path& p,
                       CGAL::Surface_mesh<K::Point_3>& out_mesh,
                       std::string* out_solid_name = nullptr);
   bool write_ascii_stl(const std::filesystem::path& p,
                        const CGAL::Surface_mesh<K::Point_3>& mesh,
                        std::string_view solid_name);
   } // namespace safs::io
   ```
   - `read_ascii_stl` uses the **standard CGAL polygon-soup
     recipe** to recover manifold connectivity from STL's
     per-triangle vertex slots (P-004 — the previously sketched
     two-step recipe was incomplete and would have failed every
     `ts_to_stl.py` STL on the first call):
     ```cpp
     std::vector<K::Point_3> points;
     std::vector<std::vector<std::size_t>> polygons;
     if (!CGAL::IO::read_polygon_soup(p.string(), points, polygons))
         return false;
     PMP::repair_polygon_soup(points, polygons,
         params::erase_all_duplicates(true)
                .require_same_orientation(false));
     PMP::orient_polygon_soup(points, polygons);
     if (!PMP::is_polygon_soup_a_polygon_mesh(polygons)) {
         std::cerr << "[corefine_faults] " << p
                   << ": polygon soup is non-manifold even after "
                   << "repair; ts_to_stl.py output is corrupt.\n";
         return false;
     }
     PMP::polygon_soup_to_polygon_mesh(points, polygons, out_mesh);
     return PMP::is_polygon_mesh(out_mesh);
     ```
     Required CGAL headers:
     ```cpp
     #include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
     #include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
     #include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
     ```
     Both the `repair_polygon_soup` step (dedup vertices + drop
     degenerate triangles) AND the
     `polygon_soup_to_polygon_mesh` step are required: ASCII STL
     stores three vertex slots per triangle, so without
     `repair_polygon_soup` every shared vertex appears K times and
     the polygon soup is non-manifold by construction.
   - `write_ascii_stl` writes `%+.17e` triple-precision matching
     `conformalize_faults.py:_write_ascii_stl` (line 1280–1284).
     The header line is `solid SAFS:<short>:conformal`, footer
     `endsolid SAFS:<short>:conformal`.  CGAL's built-in
     `CGAL::IO::write_polygon_mesh` writes binary STL by default
     (and the ASCII variant uses lower precision than we need),
     so this writer is hand-rolled.
3. **JSON I/O** (`io.cpp`):
   ```cpp
   namespace safs::io {
   struct PerFaultStats {
     std::string short_name;
     std::int64_t n_vertices;
     std::int64_t n_triangles;
     std::filesystem::path stl_path;
   };
   struct PairStats {
     std::string short_a, short_b;
     std::int64_t n_constrained_edges_A;
     std::int64_t n_constrained_edges_B;
     std::int64_t pre_n_tri_A, post_n_tri_A;
     std::int64_t pre_n_tri_B, post_n_tri_B;
     int remesh_iters;            // 0 in Phase 1
     int n_t_junctions_A;         // computed by gates (Phase 1.5; 0-stub OK)
     int n_t_junctions_B;
   };
   void write_triangle_to_fault_json(
       const std::filesystem::path& p,
       const std::vector<PerFaultStats>& per_fault);
   void write_intersection_report_json(
       const std::filesystem::path& p,
       const std::vector<PerFaultStats>& per_fault,
       const std::vector<PairStats>& pairs,
       double snap_m, double clearance_m, double target_edge_m);
   } // namespace safs::io
   ```
   - The two write functions emit the schemas defined under
     **Constraints / Interface constraints** above, byte-for-byte.
     Use a tiny inline JSON-printer (no external lib): write
     fixed-shape strings with `%.17g` for floats, `%lld` for ints,
     2-space indent.
   - `write_triangle_to_fault_json` takes the `per_fault` vector
     in the same order the user passed `--include-fault`, and
     computes contiguous ranges from the cumulative triangle
     count.
4. **Corefine wrapper** (`corefine.cpp`):
   ```cpp
   namespace safs::corefine {
   using Mesh = CGAL::Surface_mesh<
                  CGAL::Exact_predicates_inexact_constructions_kernel::Point_3>;
   using EdgeConstrainedMap =
       Mesh::Property_map<boost::graph_traits<Mesh>::edge_descriptor, bool>;

   struct CorefineResult {
     std::int64_t n_constrained_edges_A;
     std::int64_t n_constrained_edges_B;
     std::int64_t pre_n_tri_A, post_n_tri_A;
     std::int64_t pre_n_tri_B, post_n_tri_B;
     EdgeConstrainedMap ecm_A;
     EdgeConstrainedMap ecm_B;
   };

   CorefineResult corefine_pair(
       Mesh& A, Mesh& B,
       const EdgeConstrainedMap& ecm_A,    // SHARED across pairs (P-002)
       const EdgeConstrainedMap& ecm_B,    // SHARED across pairs (P-002)
       std::string_view short_a, std::string_view short_b);
   } // namespace safs::corefine
   ```
   - **Important — ECM ownership pattern (P-002).** The
     `ecm_A` / `ecm_B` are NOT created inside `corefine_pair`.
     They are installed on each fault **once** by `main.cpp`
     before the first corefine, and the same `ecm` is passed
     to every `corefine_pair` call that fault participates in.
     This is the recipe from CGAL's
     `corefinement_difference_remeshed.cpp` example:
     - When a previously-marked edge is split during a later
       corefine, both sub-edges inherit the mark (CGAL guarantee).
     - New intersection edges from later corefines are added
       to the same map alongside older marks.
     - No `accumulate_ecm` step is needed; constructing
       per-pair-local ECMs would actively *break* this invariant
       because per-pair maps don't see edges the local pair
       didn't intersect.
   - Implementation: call
     ```cpp
     PMP::corefine(A, B,
         params::edge_is_constrained_map(ecm_A)
                .throw_on_self_intersection(true),
         params::edge_is_constrained_map(ecm_B));
     ```
     Catch `CGAL::Polygon_mesh_processing::Corefinement::Self_intersection_exception`
     and rethrow as a `std::runtime_error` with which fault is
     self-intersecting (so the caller can attribute the failure).
   - **Coincident-triangle pre-flight (P-011).** Before calling
     `PMP::corefine`, hash each face of A and B by
     `(min_xyz, mid_xyz, max_xyz)` of its 3 vertex coordinates.
     If any face of A has the same hash as any face of B, exit
     2 with both fault names + triangle indices and the message:
     ```
     fault {a} tri {i} and fault {b} tri {j} share all 3
     vertices.  CGAL corefine() has undefined behaviour on
     coincident triangles.  Investigate the CFM source — this
     typically indicates a triangulation that should not have
     been split into two faults
     (PLAN_multifault_intersections.md §Caveats).
     ```
     Hash collision must trigger an explicit 3-D coordinate
     comparison (no false positives on hash collision alone).
   - Count constrained edges by iterating `edges(M)` and summing
     `get(ecm, e)`. Record pre/post triangle counts.
5. **Pairwise orchestration** (`main.cpp`):
   - Read every fault listed by `--include-fault` into a
     `std::vector<Mesh>` (size N).
   - **Install one ECM per fault BEFORE any corefine runs**
     (P-002):
     ```cpp
     std::vector<EdgeConstrainedMap> ecms(N);
     for (std::size_t i = 0; i < N; ++i) {
         auto pm = meshes[i].add_property_map<edge_descriptor, bool>(
                       "e:is_constrained", false);
         if (!pm.second) {
             throw std::runtime_error(
                 "ECM property map already exists on fault " + shorts[i]);
         }
         ecms[i] = pm.first;
     }
     ```
   - For each unordered pair `(i, j)` with `i < j`, call
     ```cpp
     corefine_pair(meshes[i], meshes[j], ecms[i], ecms[j],
                   shorts[i], shorts[j]);
     ```
     Each call **mutates both meshes in place AND updates the
     persistent ECM**, so cascade across 3+ faults works for
     free: when pair (A, C) splits an edge that pair (A, B)
     marked, the sub-edges keep the mark; when (A, C) finds a
     new intersection, that mark is added to the same `ecms[A]`.
     This replaces the Python plan's "Step (3) cascade re-scan."
   - After all pairs are processed, write each mesh out as STL
     with `safs::io::write_ascii_stl`.
   - Compute per-fault stats and the pair stats vector.
   - Emit `triangle_to_fault.json` and `intersection_report.json`
     in `out_stl_dir`.
6. **Validation gates inside the tool** (Phase 1, lightweight):
   - **Manifold gate**: after corefine, check
     `PMP::is_polygon_mesh_polyhedral(M)` and
     `PMP::does_self_intersect(M, …)`. Report
     `n_t_junctions_X` as 0 in Phase 1 (CGAL doesn't expose
     T-junctions directly because corefine output should never
     have any; rely on the gate fact). Document this; Phase 2
     adds the explicit count if needed.
   - **Polyline-edge-coincidence gate**: not asserted in code
     (it is a tautology of corefine — the constrained edges are
     by construction the same on both sides). Document.
   - **Interior-crossing-only gate**: re-corefine each pair with
     fresh empty constrained maps; assert
     `n_constrained_edges_post == 0`. (CGAL's corefine on already-
     corefined meshes adds no new constrained edges. This is the
     CGAL equivalent of the Python plan's
     `tri_tri_intersect_3d_interior_only` check.) Time-budget
     ~1 s per pair; acceptable.

### Interfaces

- New CLI:
  `corefine_faults --in-stl-dir <p> --out-stl-dir <p> --include-fault <s> [--include-fault <s>]... [--snap-m <f>] [--clearance-m <f>]`
- Output files in `<out_stl_dir>`:
  - `<short>.stl` for every included fault.
  - `triangle_to_fault.json` (presence-only marker per P-005).
  - `intersection_report.json` (diagnostic-only per P-005;
    schema delta documented in §Constraints).
- Schema validation: `triangle_to_fault.json` must satisfy
  `n_total_triangles == sum_per_fault(n_triangles)` and
  `range[i].end == range[i+1].start` (contiguous partition).
  No byte-equivalence to the Python output is required because
  no downstream tool reads the file's content.

### Edge cases to handle

- **Self-intersecting input fault**: `PMP::corefine` raises
  `Self_intersection_exception`. Catch, log which fault and
  which triangle, exit 2. Do NOT try to repair silently — that
  is `ts_to_stl.py`'s job.
- **Single fault, no pairs**: corefine no-op; write the input
  STL unchanged to `out_stl_dir`. The single-fault smoke test
  (`PLAN_smoke_millcreek.md`) must still produce identical output
  to the Python `conformalize` no-op pass-through.
- **Disjoint faults** (no triangle pairs intersect): corefine
  returns with 0 constrained edges; emit empty `pairs[…]` entry
  for that pair. The disjoint-subsets case must be a no-op.
- **Coplanar fault pairs** (e.g., adjacent SAF segments — see
  REVIEW.md R-009): CGAL's corefine handles these natively
  (returns the on-plane intersection polygon). No special-case
  code needed; document in source comment that this is the
  case the Python implementation skipped silently.
- **Empty `--include-fault` list**: exit 1 with "must specify
  at least one --include-fault".
- **Fault listed but file missing**: exit 1 with the missing path.

### Acceptance criteria

- [ ] `corefine_faults --include-fault A --include-fault B
      --in-stl-dir tools/tests/fixtures/two_squares
      --out-stl-dir /tmp/out` produces:
  - `/tmp/out/A.stl` and `/tmp/out/B.stl` with the intersection
    edge present in both;
  - `triangle_to_fault.json` whose `n_total_triangles` equals
    the sum of post-corefine triangle counts;
  - `intersection_report.json` with one `pairs.A__x__B` entry
    showing `n_constrained_edges_A == n_constrained_edges_B > 0`.
- [ ] **3-fault cascade test (P-002 verification).** Three
      orthogonal unit squares (A in xy-plane, B in xz-plane,
      C in yz-plane).  After all three pairwise corefines, fault A's
      ECM marks edges along BOTH the x-axis polyline (A∩B) and
      the y-axis polyline (A∩C).  Without the shared-ECM pattern,
      one of those polylines would be unmarked → assertion
      catches the bug.
- [ ] **Coincident-triangle preflight (P-011).** Two faults
      where one triangle of A and one triangle of B share all
      3 vertex coordinates → tool exits 2 with both indices
      named, before invoking `PMP::corefine`.
- [ ] Re-running `corefine_faults` on its own output (idempotency)
  produces 0 new constrained edges. (Interior-crossing-only gate.)
- [ ] `test_corefine_smoke` (the hand-rolled C++ test) passes
  via `ctest` in `tools/build/`.
- [ ] On the Mill Creek × SBMT-SAF 2000 m fixture
  (`output/two_crossing_2000m/stl_raw/`), the tool runs in ≤ 30 s
  on a laptop and reports
  `n_constrained_edges_A == n_constrained_edges_B`, with the
  geometrically meaningful invariant being **symmetry** plus a
  numeric floor of `> 50` (to detect "corefine found nothing").
  The empirical pin on this fixture is **143 ± 10** constrained
  edges per side; record the value in
  `intersection_report.json` so future CGAL upgrades can be
  A/B-compared.

  *Why not the originally-stated 500 floor.*  The first revision
  of this plan extrapolated from the Python pipeline's "581
  polyline vertices" (treating vertices ≈ edges).  Empirically
  CGAL EPICK on this fixture yields 143 constrained edges per
  side, runtime 22 ms.  The Python tool's vertex count is an
  over-count caused by snap-induced sub-edges that CGAL's exact
  predicates collapse into single edges.  The geometrically
  consistent quantity is segment count: at the 2000 m target
  edge length, ~143 segments span the ~286 km Mill Creek × SAF
  crossing length, which matches.  See R-002 fix in
  `REVIEW_phase01_code.md` for the reconciliation discussion.
- [ ] `triangle_to_fault.json` schema invariant holds:
      `n_total_triangles == sum(faults[*].n_triangles)` and
      ranges are contiguous.  Per P-005, byte-equivalence to
      Python output is NOT required (the file is a presence-only
      marker downstream).
- [ ] **Python pipeline regression**: `pytest miniapps/seas/safs/mesh/tests/`
  still passes 100% (no Python files were touched).

### Dependencies

- Depends on: Phase 0 (build env, CMake stub).
- Required by: Phase 2 (which adds remeshing on top), Phase 3
  (pipeline integration).

---

## Phase 2 — Add isotropic remeshing with constrained-edge protection

### Goal

After this phase, the corefined output has uniform-target-edge
triangles (configurable, default 1000 m) **except** along the
intersection polylines, where the constrained edges are
protected from the remesher's edge collapse / split / flip
operations. The slivers that Phase 1's bare corefine produces
are gone, and validator check 6 (mean fault edge length within
[500, 2000] m for `res_f = 1000 m`) passes.

### Files to create

- None (extends Phase 1's files).

### Files to modify

- `miniapps/seas/safs/tools/corefine_faults/corefine.hpp` /
  `corefine.cpp` — add `remesh_one_fault` API.
- `miniapps/seas/safs/tools/corefine_faults/main.cpp` — wire
  `--target-edge-m` and `--remesh-iters` flags through to the
  remesher.
- `miniapps/seas/safs/tools/tests/test_corefine_smoke.cpp` — add
  a remeshing-protection assertion (the constrained edge count
  must not change after remeshing).

### Detailed requirements

1. **Remesh API** (`corefine.hpp`):
   ```cpp
   namespace safs::corefine {
   struct RemeshParams {
     double target_edge_m   = 1000.0;
     int    n_iterations    = 3;       // PMP default
     double clearance_m     = 0.0;     // assertion floor; not a clamp
     bool   allow_z_clamp   = false;   // P-008 — opt-in clamp;
                                        // default OFF so a violation
                                        // surfaces as an exit-2 error
                                        // rather than silent correction
   };
   void remesh_one_fault(
       Mesh& M,
       const EdgeConstrainedMap& ecm,
       const RemeshParams& p);
   } // namespace safs::corefine
   ```
2. **Pre-split long constrained edges (P-001).**  CGAL's
   `protect_constraints(true)` prevents `isotropic_remeshing`
   from splitting OR collapsing constrained edges.  Without a
   pre-split pass, polyline edges that arrive at ~1500 m (raw
   CFM resolution) will stay at 1500 m while the rest of the
   fault remeshes to `target_edge_m = 1000 m`, locking in
   1.5:1 sliver triangles along every polyline edge.  The
   recipe is:
   ```cpp
   #include <CGAL/Polygon_mesh_processing/repair.h>      // split_long_edges
   #include <CGAL/Polygon_mesh_processing/remesh.h>      // isotropic_remeshing
   // ...
   {
       // Build a filtered range containing only constrained edges.
       std::vector<edge_descriptor> long_constrained;
       for (auto e : edges(M))
           if (get(ecm, e)) long_constrained.push_back(e);
       PMP::split_long_edges(long_constrained, p.target_edge_m, M,
           params::edge_is_constrained_map(ecm));
   }
   // Now every constrained edge is ≤ p.target_edge_m, with the new
   // sub-edges automatically marked in `ecm`.
   ```
   CGAL guarantees that splits propagate the constrained mark
   to both sub-edges, so this pre-pass extends the polyline
   into more (shorter) constrained edges without losing any.
3. **Run isotropic_remeshing with constraints fully protected.**
   ```cpp
   PMP::isotropic_remeshing(
       faces(M), p.target_edge_m, M,
       params::edge_is_constrained_map(ecm)
              .number_of_iterations(p.n_iterations)
              .protect_constraints(true)
              .relax_constraints(false));
   ```
   - `protect_constraints(true)` guarantees no constrained edge
     is split, collapsed, or flipped.  Combined with the §2
     pre-split, this means every constrained edge ends up in
     `[target_edge_m / 2, target_edge_m]` — sliver-free.
   - `relax_constraints(false)` further forbids relaxing
     constrained vertices toward the centroid — they stay put.
4. **Free-surface-clearance assertion (P-008 — replaces the
   silent z-clamp).**
   - After `isotropic_remeshing`, walk every vertex.  If any
     vertex has `v.z() > -p.clearance_m + 1e-3` (1 mm tolerance
     for FP round-off), exit 2 with the offending vertex's
     coordinates and the fault's short-name:
     ```
     [corefine_faults] FATAL: vertex (x, y, z) on fault {short}
     is above z = -{clearance_m} m by {dz} m after remeshing.
     This indicates either (a) the input STL was not clamped by
     ts_to_stl.py to the same clearance, or (b) a CGAL remesher
     interior-smoothing pulled a non-locked vertex above its
     1-ring's neighbours.  Re-run ts_to_stl.py with the matching
     --free-surface-clearance value, or pass --allow-z-clamp to
     project violators down to the clearance plane (see Phase 2
     §P-008 for the trade-offs).
     ```
   - Rationale: the input STL is already clamped by
     `ts_to_stl.py`.  `protect_constraints(true)` plus
     `relax_constraints(false)` keeps boundary and constrained
     vertices fixed.  Interior-vertex Laplacian smoothing in the
     remesher cannot pull a vertex above its 1-ring's maximum z,
     which is bounded by neighbours that are themselves below
     the clearance plane.  If the assertion fires, that's a
     CGAL bug or a pathological input — surface it for diagnosis,
     do NOT silently clamp.
   - **Opt-in `--allow-z-clamp` flag.** When set, the tool first
     checks that no constrained vertex would be moved (constrained
     vertices on the polyline must NOT be displaced because doing
     so would break the bit-equality of polyline endpoints
     between fault A and fault B).  If a constrained vertex
     would be displaced, exit 1 with "rerun ts_to_stl.py with the
     matching clearance".  Otherwise project non-constrained
     violators down to `z = -p.clearance_m` and continue.  Also
     refuse to run the clamp if it would invert any triangle
     (post-clamp face-normal sign != pre-clamp face-normal sign).
5. **Cascade is handled by the shared ECM (P-002), not by an
   accumulate step.**  Because every `corefine_pair` call uses
   the persistent `ecms[i]` for fault i, that single ECM
   accumulates marks from every pair fault i participated in.
   `remesh_one_fault(meshes[i], ecms[i], …)` then sees ALL of
   fault i's constrained edges in one map.  No separate
   `accumulate_ecm` helper is needed; building one would be a
   bug (per P-002 in REVIEW_plan_cgal_corefine.md).
6. **Wire CLI**:
   - `--target-edge-m <float>` default `1000.0`.
   - `--remesh-iters <int>` default `3` (PMP default; sufficient
     to remove slivers per CGAL docs).
   - `--no-remesh` flag for diagnostic / cross-check with Phase 1.
   - When `--target-edge-m == 0`, skip remeshing (back-compat to
     Phase 1 behaviour).
   - `--allow-z-clamp` flag (P-008): default OFF.
7. **Update JSON report**:
   - `intersection_report.json`'s per-pair entry gets
     `"remesh_iters": <int>`. Per-fault entry gets the post-remesh
     triangle / vertex counts (already populated; just ensure the
     numbers reflect the post-remesh mesh).
   - Add top-level `"target_edge_m": <float>`.

### Interfaces

- The `corefine_faults` CLI gains `--target-edge-m`,
  `--remesh-iters`, `--no-remesh`, `--allow-z-clamp`.
  Existing flags unchanged.

### Edge cases to handle

- **Constrained edges accidentally split during
  isotropic_remeshing** (CGAL bug in older versions or a
  pathological input): assert
  `n_constrained_edges_post_remesh == n_constrained_edges_post_split_long_edges`
  (note: NOT the pre-split-long-edges count, which is smaller —
  the pre-split pass is supposed to ADD constrained sub-edges).
  If different, exit 2 with the count mismatch.
- **`target_edge_m` larger than fault diameter**: PMP will
  decimate aggressively, possibly to 1–2 triangles. Detect and
  warn if the post-remesh triangle count is < 4 per fault; do
  not fail (the user may want a coarse cross-check mesh).
- **Boundary vertices** (open-mesh edges): PMP's default behaviour
  is to keep them; additionally, mark every boundary edge in the
  shared ECM at the start of remeshing so they are protected by
  the same machinery as polyline edges:
  ```cpp
  for (auto e : edges(M))
      if (is_border(e, M)) put(ecm, e, true);
  ```
- **Coplanar pair already remeshed by an earlier pair**: the
  shared ECM (P-002) covers this — no additional logic needed.

### Acceptance criteria

- [ ] After running `corefine_faults` with default flags
  (`--target-edge-m 1000`) on Mill Creek × SBMT-SAF 2000 m, every
  output STL satisfies (P-001 + P-006):
    - mean edge length in [800, 1200] m;
    - **p5 edge length ≥ 250 m** (= target / 4 — sliver gate);
    - **p95 edge length ≤ 2500 m** (= 2.5 × target — long-tail
      gate);
    - **max aspect ratio ≤ 4** (= max-edge / min-edge per
      triangle) on any triangle incident to a constrained edge.
- [ ] After the §2 pre-split pass, every constrained edge has
      length ≤ `target_edge_m * 1.05` (PMP `split_long_edges`
      docs guarantee target × 4/3, but in practice the
      bisection produces edges within 5% of target).
- [ ] After `isotropic_remeshing`, constrained-edge count is
      preserved exactly:
      `n_constrained_edges_post_remesh == n_constrained_edges_post_split_long_edges`.
- [ ] **Free-surface assertion (P-008)**: with default flags,
      no vertex has `z > -clearance_m + 1e-3 m`.  A synthetic
      test that artificially places a non-constrained interior
      vertex at `z = -50 m` (with `clearance_m = 100`) triggers
      exit 2 — NOT a silent clamp.  The same input with
      `--allow-z-clamp` runs to completion and the output has
      that vertex at `z = -100 m`.
- [ ] `--no-remesh` produces output identical to Phase 1.
- [ ] Smoke test runs in ≤ 60 s on Mill Creek × SBMT-SAF.
- [ ] Schema regression: `intersection_report.json` now contains
      `target_edge_m` and `remesh_iters` per pair.
- [ ] Python pipeline regression suite still passes (no Python
      changes).

### Dependencies

- Depends on: Phase 1.
- Required by: Phase 3.

---

## Phase 3 — Pipeline integration (parallel run script)

### Goal

After this phase, `bash run_two_crossing_2000m_cgal.sh` runs the
full pipeline (audit → ts_to_stl → corefine_faults →
generate_safs_mesh → write_fault_provenance → validate_msh →
convert_msh) and produces a mesh that passes **at least 9/11**
validator checks (must-pass, P-003).  Checks 7 and 11 are
**stretch** because they depend on gmsh's `safs.geo` size-field
configuration (Threshold / DistMin / DistMax / `tube_radius` /
`ramp_dist`), not on corefine output, and may require a
follow-up size-field tuning pass to flip to PASS.  The existing
`run_two_crossing_2000m.sh` is unchanged and still uses the
Python conformalizer.

### Files to create

- `miniapps/seas/safs/mesh/run_two_crossing_2000m_cgal.sh` —
  near-clone of `run_two_crossing_2000m.sh` with step 3 swapped
  from Python to the CGAL binary, and `ALGO3D` defaulted to 10.

### Files to modify

- None. (`generate_safs_mesh.py`, `write_fault_provenance.py`,
  `validate_msh.py`, `convert_msh.py` are deliberately untouched.)

### Detailed requirements

1. **Shell script** `run_two_crossing_2000m_cgal.sh`:
   - Copy of `run_two_crossing_2000m.sh` with the following diff:
     ```diff
     -ALGO3D="${TWO_ALGO3D:-4}"
     +ALGO3D="${TWO_ALGO3D:-10}"   # HXT — corefine output is a
     +                             # valid PLC, HXT will not reject.
     +COREFINE_BIN="${COREFINE_BIN:-$(dirname "$0")/../tools/build/corefine_faults}"

     -echo "==> 3/7 conformalize_faults (target_edge_length=${TARGET_EDGE} m)"
     -python conformalize_faults.py \
     -    --in-stl-dir "$OUTDIR/stl_raw" \
     -    --out-stl-dir "$OUTDIR/stl_conformal" "${INCL[@]}" \
     -    --clearance-m "$CLEARANCE"
     +echo "==> 3/7 corefine_faults (CGAL; target_edge=${TARGET_EDGE} m)"
     +"$COREFINE_BIN" \
     +    --in-stl-dir "$OUTDIR/stl_raw" \
     +    --out-stl-dir "$OUTDIR/stl_conformal" "${INCL[@]}" \
     +    --clearance-m "$CLEARANCE" \
     +    --target-edge-m "$TARGET_EDGE" \
     +    --remesh-iters 3
     ```
   - The output directory must be a *separate* path from the
     Python run so A/B comparison works.  In the cloned
     `_cgal.sh` script, set near the top:
     ```bash
     OUTPUT_SUFFIX="${TWO_OUTPUT_SUFFIX:-_cgal}"
     ```
     and modify line 41 of the original (which is hardcoded
     `OUTDIR=output/two_crossing_${RES}m`) in the **clone** to:
     ```diff
     -OUTDIR=output/two_crossing_${RES}m
     +OUTDIR=output/two_crossing_${RES}m${OUTPUT_SUFFIX}
     ```
     Result: CGAL run lands in
     `output/two_crossing_2000m_cgal/`, Python run lands in
     `output/two_crossing_2000m/` (the legacy script keeps the
     empty default once Phase 4 inverts the relationship — see
     P-009).
2. **No changes to** `safs.geo`, `generate_safs_mesh.py`,
   `validate_msh.py`. The new shell script proves the contract
   already exists.
3. **Comparison harness** (optional, tested but not required to
   land):
   - Add `tools/compare_meshes.py` — loads two `.msh` files,
     reports `(n_internal_with_2_tets, mean_fault_edge,
     n_offending)` for each. Lets the user A/B Python vs CGAL.

### Interfaces

- Pipeline entry point: `run_two_crossing_2000m_cgal.sh`. No new
  Python entry points.

### Edge cases to handle

- **`COREFINE_BIN` not found**: shell script aborts with
  "corefine_faults binary not found at $COREFINE_BIN; build it
  first via cd tools && cmake -B build && cmake --build build".
- **HXT rejects despite corefine**: this *should not happen*
  if Phase 1 + 2 acceptance criteria are met. If it does, the
  diagnostic is to run the smoke test from Phase 1 on the same
  inputs and report the constrained-edge count + manifold gates.
- **algo3d=10 with single-fault input**: HXT meshes single-fault
  input fine (the smoke test
  `run_smoke_millcreek.sh` already proves this). No special
  handling.
- **The user invokes the script in a worktree without a build**:
  the `COREFINE_BIN` check handles it.

### Acceptance criteria (P-003 — must-pass / stretch split)

**Must-pass (corefine + remesh own these — Phase 3 is gated on these):**

- [ ] `bash mesh/run_two_crossing_2000m_cgal.sh` exits 0.
- [ ] `validate_msh.py` reports **at least 9 of 11 checks
      passed**.  The 9 that must pass:
  - checks 1–4, 8–10 (already pass in the existing 7/11
    baseline; CGAL must not regress them);
  - **check 5 (`internal_interface`)**: `n_offending == 0`,
    `n_internal_with_2_tets ≈ n_fault_tris − n_trace_tris`.
    *This is the primary acceptance bar — the proof the
    corefine+HXT path actually embeds the fault.*
  - **check 6 (`fault_edge_length`)**: mean edge in
    [500, 2000] m.  *This is the proof the remesh worked.*
- [ ] The combine step's `n_dropped_duplicate` is 0 (REVIEW.md
      R-002 is moot under CGAL because polyline vertices are
      bit-identical between fault A and fault B).
- [ ] Run time end-to-end ≤ 90 s on a laptop.
- [ ] The Python pipeline (`run_two_crossing_2000m.sh`) still
      runs unchanged and produces its old (broken) output, for
      side-by-side comparison.

**Stretch (size-field-dependent — separate ticket if needed):**

- [ ] **check 7 (`far_field_edge`)**: far-field mean tet edge
      ≥ 0.6 × `res_ff` = 12000 m.  This is governed by gmsh's
      Threshold field (`Field[2].DistMax = tube_radius +
      ramp_dist`) and `Mesh.MeshSizeMax = res_ff`.  With
      `algo3d=10` (HXT) and a properly-embedded fault, the
      threshold field becomes effective and far tets reach
      `res_ff`.  But this is not guaranteed by anything Phase
      1+2 do; it depends on the size-field tuning that already
      lives in `safs.geo`.  If it fails, the diagnostic is in
      the existing `validation_report.txt` — the metric
      `mean_far_edge_m` will report below threshold and the
      remediation is to increase `--ramp-dist` or reduce
      `--tube-radius` in the run script (no source change).
- [ ] **check 11 (`tube_uniformity`)**: mean tube tet edge in
      [500, 1500] m.  Same caveat as check 7.

**Stretch goal — 11/11.**  After at most one round of size-field
tuning (CLI flag changes only — `--res-ff`, `--ramp-dist`,
`--tube-radius`), all 11 checks pass.  The exact values that
achieve 11/11 are recorded in
`output/two_crossing_2000m_cgal/output/sizing.json` and pinned
as the new defaults in `run_two_crossing_2000m_cgal.sh`.  This
work is **conditional** on the must-pass criteria above being
met first; if they are, and 7/11 still fail, file a follow-up
ticket "size-field tuning for the CGAL pipeline" rather than
extending Phase 3.

If the must-pass criteria fail (check 5 or check 6 still
broken), Phase 3 is NOT complete; do not proceed to Phase 4.

### Dependencies

- Depends on: Phases 0 + 1 + 2.
- Required by: Phase 4.

---

## Phase 4 — Migration: make CGAL the default

### Goal

After this phase, `run_two_crossing_2000m.sh` is a thin wrapper
that delegates to the canonical `run_two_crossing_2000m_cgal.sh`
(P-009 — keep the `_cgal` script as the canonical recipe so
future regressions can be A/B tested without git archaeology).
The Python conformalizer is documented as legacy.  The Python
sources stay in the tree as a regression fixture but are no
longer on the production hot-path.

### Files to create

- None.

### Files to modify

- `miniapps/seas/safs/mesh/run_two_crossing_2000m_cgal.sh` —
  **keep as the canonical recipe.**  Add a header comment
  documenting that this is the production entry point and
  `run_two_crossing_2000m.sh` is a wrapper.  No code change.
- `miniapps/seas/safs/mesh/run_two_crossing_2000m.sh` —
  **replace its body with a 2-line wrapper** that defaults
  `OUTPUT_SUFFIX` to empty (so output lands in the canonical
  `output/two_crossing_2000m/` path) and `exec`s the `_cgal`
  script:
  ```bash
  #!/usr/bin/env bash
  # Production entry point.  Delegates to the canonical CGAL
  # pipeline at run_two_crossing_2000m_cgal.sh.  Set
  # OUTPUT_SUFFIX explicitly to land output elsewhere; the
  # default empty suffix keeps the canonical path stable.
  exec "$(dirname "$0")/run_two_crossing_2000m_cgal.sh" "$@"
  ```
  with `TWO_OUTPUT_SUFFIX=""` exported above the `exec` (so
  the inner script's `OUTPUT_SUFFIX="${TWO_OUTPUT_SUFFIX:-_cgal}"`
  default is overridden to empty).
- `miniapps/seas/safs/PLAN_multifault_intersections.md` — append
  a "Status: superseded by `PLAN_cgal_corefine.md`" banner at
  the top, but keep the document as historical record. Mark
  Phases 1 + 2 as legacy.
- `miniapps/seas/safs/mesh/conformalize_faults.py` — add a
  module-level deprecation warning (no behaviour change):
  ```python
  warnings.warn(
      "conformalize_faults.py is the legacy Python conformalizer; "
      "the production path is tools/corefine_faults (CGAL).  "
      "See PLAN_cgal_corefine.md.",
      DeprecationWarning,
  )
  ```
- `miniapps/seas/safs/mesh/tests/test_conformalize_faults.py`
  and `test_fault_intersect.py` — annotate the test module with
  `pytest.mark.legacy`; add a `pytest.ini` filter so the legacy
  tests run only when explicitly requested
  (`pytest -m legacy …`).
- `miniapps/seas/safs/mesh/tests/README.md` (new file) — document
  the legacy-test policy (P-013):
  > - `pytest -m "not legacy"` (the default) skips legacy tests.
  > - `pytest -m legacy` runs the Python conformalizer regression
  >   suite explicitly.  Required to pass at every release.
  > - **90-day retention clock**: starts at Phase-4 merge.  After
  >   90 days, if no CGAL regression has surfaced AND the legacy
  >   tests have not caught any production-relevant bug, the
  >   legacy tests + Python sources (`fault_intersect.py`,
  >   `conformalize_faults.py`) are removable in a follow-up
  >   ticket.  The clock is reset by any actual CGAL regression
  >   that the legacy tests catch (re-evaluate retention then).

### Detailed requirements

1. The wrapper redirection above keeps the production CLI
   (`run_two_crossing_2000m.sh`) byte-equivalent for users:
   same arguments, same output path.
2. The deprecation warning fires on `import conformalize_faults`.
   It must be `DeprecationWarning`, not `FutureWarning`, so it is
   silent by default but visible when running tests with
   `-W default::DeprecationWarning`.
3. **Do not delete** `conformalize_faults.py` or
   `fault_intersect.py`. They are the legacy regression fixture.
   This is explicit: removal is out of scope for this plan and
   is gated on the 90-day retention clock above (P-013).

### Interfaces

- No new interfaces. The Python CLI of `conformalize_faults.py`
  still works (with deprecation warning on stderr).
  `run_two_crossing_2000m.sh` is now a thin wrapper around
  `run_two_crossing_2000m_cgal.sh`.

### Edge cases to handle

- **A user has scripts that import `conformalize_faults`**: the
  deprecation warning is informational only; functionality is
  preserved.
- **A user wants to A/B Python vs CGAL after Phase 4**: invoke
  `python conformalize_faults.py …` directly with its existing
  CLI (the Python pipeline still lives in `mesh/`); the wrapper
  always delegates to CGAL.
- **A user has cron jobs using `run_two_crossing_2000m.sh` with
  custom output paths**: the wrapper passes `"$@"` through, so
  any custom `OUTDIR` set via env vars still propagates.

### Acceptance criteria

- [ ] `bash mesh/run_two_crossing_2000m.sh` exits 0 and produces
      ≥ 9/11 validator checks (must-pass per P-003), via the
      `_cgal` wrapper.
- [ ] `mesh/run_two_crossing_2000m.sh` is a 5-line wrapper
      (shebang + comment + 2 env exports + 1 exec) — no
      pipeline logic duplicated.
- [ ] `pytest miniapps/seas/safs/mesh/tests/` (without `-m legacy`)
      runs only the non-legacy tests; the legacy Python regression
      suite is gated behind `-m legacy`.
- [ ] `import conformalize_faults` still works and emits exactly
      one `DeprecationWarning`.
- [ ] `mesh/tests/README.md` documents the 90-day retention
      clock starting at Phase-4 merge date (P-013).

### Dependencies

- Depends on: Phase 3.
- Required by: Phase 5.

---

## Phase 5 — Scale up to all-8 (acceptance criteria only)

### Goal

This phase is **out of scope for the first cut**. It documents
the acceptance bar that the all-8 build must clear, so the
next milestone has unambiguous targets. No code is written here.

### Files to (eventually) modify

- `miniapps/seas/safs/mesh/run_all8_2000m.sh` — swap to
  `corefine_faults` (analogous diff to Phase 4).

### Acceptance criteria (eventual; do not implement now)

- [ ] `bash mesh/run_all8_2000m.sh` exits 0 with **at least 9/11
  checks** (the must-pass set per P-003).  11/11 is a stretch
  goal conditional on size-field tuning, same split as Phase 3.
- [ ] Run time ≤ 5 min on a laptop. CGAL's corefine on
  ~6500 triangles × 8 faults is documented as quadratic-pair
  but per-pair fast (~0.5 s); 28 pairs → ~14 s of corefine,
  + ~1 min of remeshing, + ~3 min of HXT.
- [ ] No fault has self-intersection → corefine throws → exit 2.
  Per `REVIEW_phase1.md` and the all-8 README, Pinto Mountain
  and Banning standalone are CFM-source-defective; they may
  fail this gate. Phase 5 documents that **fixing CFM source
  defects is upstream of this plan** and out of scope.
- [ ] Coincident-triangle preflight (P-011) does not fire on
  the all-8 set.  If it does, that's the same CFM-source-defect
  category as self-intersection — investigate, do not work
  around in this tool.
- [ ] Frontera build recipe (`tools/README.md` Phase 0 §3) is
  validated by an actual Frontera build, requiring user
  approval per project memory `feedback_frontera_approval.md`.
  The validation reuses the pinned module list and CGAL 5.6.1
  fetch from P-010.

### Dependencies

- Depends on: Phase 4.
- Required by: nothing (this is the last documented milestone).

---

## Testing strategy

### Per-phase tests

| Phase | Test type | What |
|---|---|---|
| 0 | Build | `cmake --build` + run `corefine_faults_hello`. |
| 1 | C++ smoke (`ctest`) | Two perpendicular unit squares: post-corefine, both meshes contain the intersection edge; idempotent re-corefine adds 0 constrained edges. |
| 1 | C++ smoke | **3-fault cascade (P-002)**: three orthogonal squares; fault A's shared ECM marks edges along both A∩B and A∩C polylines after both corefines. |
| 1 | C++ smoke | **Coincident-triangle preflight (P-011)**: two faults sharing a triangle → exit 2 with both indices, before invoking PMP::corefine. |
| 1 | Real fixture | Mill Creek × SBMT-SAF 2000 m: emit `triangle_to_fault.json` whose schema invariant (contiguous ranges, total = sum) holds.  Per P-005, byte-equivalence to the Python output is NOT required — the file is presence-only downstream. |
| 2 | C++ smoke | Same fixture; assert `n_constrained_edges_post_remesh == n_constrained_edges_post_split_long_edges` (P-001). |
| 2 | C++ smoke | **Z-clamp assertion (P-008)**: synthetic mesh with one non-locked vertex at z = -50 m and clearance = 100 → exit 2 with default flags; runs to completion under `--allow-z-clamp`. |
| 2 | Real fixture | Mill Creek × SBMT-SAF: mean edge in [800, 1200] m, p5 ≥ 250 m, p95 ≤ 2500 m, max-aspect-ratio ≤ 4 on triangles incident to a constrained edge (P-001 + P-006). |
| 3 | End-to-end (must-pass) | `run_two_crossing_2000m_cgal.sh` → ≥ 9/11 validator checks; in particular, check 5 `n_offending == 0` and check 6 mean edge in [500, 2000] m. |
| 3 | End-to-end (stretch) | After at most one round of size-field CLI-flag tuning, 11/11 (P-003).  Stretch is OK to defer to a follow-up ticket. |
| 3 | Side-by-side | `tools/compare_meshes.py` Python vs CGAL output: CGAL must report `n_offending == 0` while Python reports `n_offending == 36505`. |
| 4 | Migration | `run_two_crossing_2000m.sh` is now a thin wrapper exec'ing the `_cgal` script (P-009); deprecation warning fires once on `import conformalize_faults`. |
| 5 | All-8 | Out of scope (documented).  Acceptance reuses the must-pass / stretch split from P-003. |

### Tests **not** to write

- Do not duplicate the Python `tests/test_fault_intersect.py` or
  `tests/test_conformalize_faults.py` for the C++ tool. Those
  tests cover algorithmic invariants of the *Python*
  implementation. The CGAL tool's invariants are covered by the
  CGAL library's own test suite + the JSON schema-invariant
  checks (P-005 — contiguous ranges only, not byte-equivalence)
  + the end-to-end validator output.

### Validation against analytical reference

- For two perpendicular unit squares centred at the origin, the
  intersection is the segment from (0, 0, -1) to (0, 0, +1)
  — exactly 2.0 m, exactly along the z-axis. The C++ smoke test
  asserts:
  - Both meshes contain a contiguous chain of constrained
    edges spanning z ∈ [-1, +1].
  - The polyline endpoints are bit-equal between fault A and
    fault B (this is the property CGAL gives that the Python
    pipeline did not).

---

## Risk assessment

1. **CGAL build complexity on Frontera.** CGAL is header-only
   since 5.0, but Boost on Frontera is at version 1.71 which is
   compatible. Risk: mitigated by deferring Frontera work to
   Phase 5 and requiring user approval. Failure mode: documented
   recipe doesn't compile → user pings; not a Phase 1–4 blocker.

2. **`PMP::corefine` self-intersection exception on real CFM
   data.** The all-8 build has known-bad faults (Pinto Mountain,
   Banning). Risk: mitigated by scoping Phase 1–4 to the 2-fault
   target where Mill Creek and SBMT-SAF are known clean. Phase 5
   documents this as a CFM-source-defect issue, not a CGAL
   issue.

3. **Schema drift between CGAL JSON and Python JSON.** Risk:
   mitigated by hand-rolling the JSON writer to match the Python
   output byte-for-byte (modulo float formatting), and by a
   regression test that diffs the keysets (not values).

4. **`isotropic_remeshing` interior smoothing pulls vertices
   above the free surface.** Risk: addressed by Phase 2 §3
   (post-remesh z-clamp). Failure mode: validator check 4
   (free-surface trace) catches any leakage.

5. **Constrained-edge count drifts after remeshing.** Risk:
   mitigated by the Phase 2 acceptance assertion. CGAL's
   `protect_constraints(true)` is documented to be exact;
   failure would indicate a CGAL bug, not our bug — we'd file
   upstream and pin to a known-good CGAL version.

6. **A user expects the new CLI to match the Python CLI byte
   for byte.** Risk: mitigated by mirroring flag names exactly
   (`--in-stl-dir`, `--out-stl-dir`, `--include-fault`,
   `--clearance-m`, `--snap-m`). The new flags
   (`--target-edge-m`, `--remesh-iters`, `--no-remesh`) are
   additive and have sensible defaults.

7. **`Surface_mesh::Property_map<edge_descriptor, bool>` lifetime
   pitfall.** CGAL property maps are owned by their host mesh;
   when we move-construct or copy a `Mesh`, the property maps
   become dangling. Risk: mitigated by passing meshes by
   reference everywhere and never copying them. Document in
   the source.

8. **Build env conflict between `mfem-dev` and `pythonenv`.**
   The CGAL tool builds in `mfem-dev` (C++ toolchain); the
   pipeline runs in `pythonenv` (gmsh, meshio, etc.). The
   shell script must `conda activate pythonenv` for steps
   1, 2, 4, 5, 6, 7 and call the C++ binary directly (no env
   activation needed for a compiled binary, just `LD_LIBRARY_PATH`
   if dynamically linked to GMP/MPFR). Document in
   `tools/README.md`.

9. **CGAL's exact predicates use mpfr which has a thread-local
   stack**: not a parallelism concern at our scale, but if a
   future Phase 5 parallelizes corefine over MPI ranks, each
   rank needs its own mpfr context. Document; no Phase 1–4
   action.

## Out of scope (explicit)

- **Removal of `conformalize_faults.py` / `fault_intersect.py`.**
  Phase 4 deprecates; removal is a separate ticket.
- **CGAL `make_mesh_3` for volume meshing.** REVIEW.md R-010 §3
  notes this is the cleanest "all-CGAL" pipeline, but switching
  the volume mesher rewrites Phases 3 + 4 of
  `PLAN_multifault_intersections.md` (the gmsh `safs.geo` setup,
  the `--algo3d` flag, the `validate_msh.py` checks). Out of
  scope; gmsh HXT continues to be the volume mesher for this
  plan.
- **TetGen / `meshpy` alternative.** Same reason — orthogonal to
  the surface-conformity issue this plan solves.
- **Per-fault different target edge length.** Single global
  `--target-edge-m` is sufficient for the 2-fault target. Per-
  fault would multiply the CLI complexity and adds no value
  at this milestone.
- **EPECK kernel.** EPICK is sufficient for single-pass
  corefine + remesh per CGAL docs and REVIEW.md R-010. Upgrade
  to EPECK only if the all-8 build (Phase 5) shows
  self-intersection artefacts that EPECK would fix.
- **Coupled-solver work** (renaming faults to per-fault tags,
  multiple `Physical Surface` IDs) — already deferred by
  `PLAN.md` Phase 6.
- **Frontera execution** of any phase — gated on user approval
  per project memory.
