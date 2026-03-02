---

# MFEM SEAS Miniapp: System Update Plan

**Prepared:** 2026-02-28
**Scope:** Preparatory refactoring of the existing 2D antiplane SEAS miniapp — input file infrastructure, fault physical group detection, interface cleanup, unit test design, and build/run reference. This work is a prerequisite for the BP5 full elasticity extension; see `fullelasticity_implementation_plan_02282026.md` for the BP5 implementation plan.

---

## Part A: Current Code Structure Assessment

### A.1 What Is Well-Structured and Reusable

**Class hierarchy is already forward-thinking.** The base class `DomainOperator<MeshType>` in `/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/domain/domain_operator.hpp` was written with explicit attention to future generalization. Two virtual methods already signal the intent:

```cpp
virtual int NumComponents() const = 0;  // 1 (antiplane) or 3 (full elasticity)
virtual int Dimension() const = 0;      // 2 (2D) or 3 (3D)
```

The interface methods `Solve()`, `ComputeTraction()`, `GetNumFaultDOFs()`, `GetFaultDepths()` are all expressed in terms of `Vector` (flat arrays), not scalar-specific types, so the same interface works for multi-component problems without modification.

**Template mesh duality is complete.** The `<MeshType>` template parameter and the type-alias machinery in `common/seas_types.hpp` (`FESpaceForMesh`, `GridFunctionForMesh`, etc.) fully supports serial (`Mesh`) and parallel (`ParMesh`) execution. This dual-mode design needs no modification for BP5.

**The SEAS time-integration operator (`seas_operator.hpp`) is domain-agnostic.** `SEASQuasiDynamicOperator<MeshType>` holds a pointer to `DomainOpType` and a pointer to `RateStateFaultOperator<MeshType>`. The only coupling is through the `Vector`-based interfaces `Solve(t, slip, u)` and `ComputeTraction(u, slip, traction)`. A 3D elasticity domain operator that satisfies the same interface will slot in without modifying the SEAS time integrator.

**The friction layer is completely physics-agnostic.** `FrictionLaw` (abstract base), `DieterichRuinaFriction`, and `StateEvolution`/`AgingLaw` operate on scalar quantities: `V` (slip speed norm), `theta`, `a`, `eta`, `sigma_n`. For BP5, the friction law itself does not change; only the `V` fed into it changes from a scalar to the Euclidean norm of a 2-component slip rate vector. The existing classes are therefore fully reusable.

**The BR2 integrator infrastructure exists.** `dg_br2_integrator.hpp` implements `BR2InteriorFaceIntegrator` and `BR2BoundaryFaceIntegrator` for the scalar Laplace problem. The lifting-operator logic (mass matrix inverse, face-integral tensors, L_q computation) is geometrically general; only the stress-computation step needs to be generalized from scalar gradient to Voigt-form stress tensor.

**I/O and output framework is largely reusable.** `BenchmarkOutput`, `ProbeOutput`, `ParaviewOutput` use `Vector`-based interfaces and time-series file formats. The BP5 output format requires additional fields (`slip_2`, `slip_3`, `shear_stress_2`, `shear_stress_3`, off-fault displacements and velocities), but the file-writing infrastructure (header format, adaptive output frequency logic, probe interpolation) can be directly extended.

**`FaultGeometry` handles parallelism cleanly.** The gather/scatter and deduplication logic in `fault_geometry.hpp` is already parameterized by `num_fault_dofs_`. For BP5 the fault is 2D (a surface), so `fault_geometry.hpp` will need to store 2D coordinates `(x2, x3)` per DOF instead of just `depths`. The MPI gather/deduplication logic is directly applicable.

### A.2 What Is Hardcoded for 2D/Scalar and Needs Generalization

**`AntiplaneDomainOperator` is entirely scalar.** The FE space is created with a scalar `H1_FECollection`/`DG_FECollection` (1 DOF per node), the stiffness matrix uses `DiffusionIntegrator` (Laplacian), and `ComputeTraction` computes a scalar traction `tau = mu * {{d u / d x_normal}}`. All of these must be replaced for 3D elasticity.

**Fault DOF identification uses `x ≈ 0` geometry.** In `antiplane_operator.hpp`, `IsFaultFace()` checks `std::abs(x) < tol`. In 3D, the fault is the plane `x1 = 0` and interior faces on that plane need the same test. The coordinate check generalizes directly to 3D (same `x1` test), but the DOF layout changes: each fault node has 3 displacement DOFs rather than 1.

**`RateStateFaultOperator` assumes scalar slip.** The state layout `[slip_0, theta_0, slip_1, theta_1, ...]` stores one `slip` per node (scalar). For BP5 there are two tangential slip components: `slip_2` (along-strike) and `slip_3` (along-dip). The state layout must expand to `[s2_0, s3_0, psi_0, s2_1, s3_1, psi_1, ...]` with `StatePerNode = 3`. The friction law call `SolveSlipRate(tau, theta, sigma_n, eta, a)` takes scalar `tau`; for BP5 it must take the Euclidean norm `||tau||` and decompose the result into component-wise rates `Vi = V * tau_hat / ||tau_hat||` (as in Tandem's `DieterichRuinaAgeing::slip_rate`).

**`BP2Params` encodes 1D fault geometry.** The `a_of_z(z)` function depends only on depth. BP5's `a(x2, x3)` depends on both along-strike position `x2` and depth `x3` (see Eq. 14 in the spec). A new `BP5Params` struct will need a 2D parameter function.

**`FaultGeometry` stores only depth.** The member `depths_` (Vector of z-coordinates) must become a 2D coordinate set for BP5.

**`BenchmarkOutput` is wired to BP2 fields.** Column names, probe depths (1D), and field extraction are BP2-specific. BP5 output has 10 on-fault stations at 2D coordinates `(x2, x3)`, plus 9 off-fault stations at 3D positions `(x1, x2, x3)`, plus source parameter and earthquake catalog files.

**The mesh generator (`BP2MeshGenerator`) is 2D only.** It uses `Mesh::MakeCartesian2D`. BP5 requires a 3D Gmsh `.geo`/`.msh` file with the fault plane embedded as an interior surface.

**`SEASQuasiDynamicOperator` hardcodes `DomainOpType = AntiplaneDomainOperator<MeshType>`.** The line `using DomainOpType = AntiplaneDomainOperator<MeshType>` in `seas_operator.hpp` locks the operator type. This can be generalized by adding a second template parameter or using the abstract base class pointer.

**The BR2 integrator penalty parameter is hardcoded as `D+1 = 3` (2D).** For 3D it must be `D+1 = 4`.

**`FaultNodes` needs generalization for 3D.** `FaultNodes<MeshType>` in `fault/fault_nodes.hpp` maps fault interior face indices to local DOF pairs and precomputes 2×2 face mass-matrix inverses for L2 projection of tractions. In 3D (BP5), fault faces are 2D surfaces (triangles or quads) rather than 1D segments, so the face mass matrices become larger (number of face DOFs × number of face DOFs), and the L2 projection from face integrals to DOFs is more complex. The DOF-to-face mapping must also account for vector-valued displacement (3 DOFs per node) rather than scalar (1 DOF per node).

### A.3 Suggested File/Class Rearrangements

The existing structure is good. The recommended additions, not replacements, are:

```
miniapps/seas/
├── config/
│   ├── bp1_params.hpp       (existing)
│   ├── bp2_params.hpp       (existing)
│   └── bp5_params.hpp       (NEW: BP5Params with 2D a(x2,x3), tau_pre vector)
├── domain/
│   ├── domain_operator.hpp  (existing base — minor extensions needed)
│   ├── antiplane_operator.hpp   (existing)
│   ├── antiplane_bdrload_operator.hpp  (existing)
│   ├── bp2_mesh.hpp         (existing 2D mesh)
│   └── elasticity3d_operator.hpp   (NEW: 3D full elasticity domain operator)
├── fault/
│   ├── fault_geometry.hpp   (extend to store 2D coords + generalize to N-component slip)
│   ├── rate_state_fault.hpp (extend StatePerNode from 2 to 3 for BP5 via template param or subclass)
│   └── fault_geometry3d.hpp (NEW or extend: 2D fault coordinate set, 2D probe matching)
├── friction/
│   ├── friction_law.hpp     (extend: add vector-valued SolveSlipRateVector)
│   ├── dieterich_ruina.hpp  (extend: add vector slip rate solver following Tandem)
│   └── state_evolution.hpp  (unchanged)
├── integrator/
│   ├── dg_br2_integrator.hpp   (extend penalty for 3D; add elasticity variant)
│   └── dg_elasticity_br2_integrator.hpp  (NEW: BR2 for linear elasticity)
├── solver/
│   ├── seas_operator.hpp        (generalize DomainOpType template)
│   └── seas_operator3d.hpp      (optional: BP5-specific typed alias)
├── io/
│   ├── benchmark_output.hpp     (extend for BP5 2D probe coords)
│   └── bp5_benchmark_output.hpp (NEW: BP5-specific output format)
└── bp5/
    ├── benchmark_document/  (existing spec PDF)
    ├── bp5.geo              (NEW: Gmsh geometry — modeled on Tandem bp5.geo)
    └── bp5_verification.cpp (NEW: BP5 driver program)
```

The key insight is that `domain_operator.hpp` (the abstract base) needs only one addition to support the off-fault displacement output required by BP5:

```cpp
// Addition to DomainOperator interface
virtual void GetOffFaultDisplacement(const std::vector<Vector> &points,
                                     Vector &displacements) const {}
```

---


## Part B: Structured Input File Design (MOOSE-Style Driver)

### B.1 Motivation: Current Driver Complexity

The current driver files (`bp1_verification_full.cpp`, `bp1_bdrload.cpp`) embed every configuration decision directly in `main()`. A user who wants to change the mesh resolution, friction parameters, or output probe locations must edit C++ source code and recompile. The two files are also nearly identical — 90% of their ~900 lines is shared boilerplate, with the only real differences being which domain operator and fault adapter are constructed.

Concrete problems identified:

| Concern | Where it lives today | Should be |
|---|---|---|
| Mesh file path, scale factor | `--mesh`, `--mesh-scale` argv | `[mesh]` block |
| DG method (BR2 vs IP), poly order | hardcoded `order=1, DGMethod::BR2` | `[domain]` block |
| Vp, Wf, mu, rho | `MakeBP1Params()` compiled in | `[material]` + `[loading]` blocks |
| Friction a(z), b, Dc, V0, f0 | `MakeBP1Params()` compiled in | `[friction]` block |
| Far-field BC type (Natural vs Dirichlet) | choice of executable (`seas_bp1_full` vs `seas_bp1_bdrload`) | `[boundary_conditions]` block |
| RK45 tolerances, dt limits | hardcoded `AbsTol=1e-7`, `DtMax=0.5yr` | `[time_integration]` block |
| Output dir, prefix, probe depths | argv + `GetBP1ProbeDepths()` compiled in | `[output]` block |
| Earthquake detection thresholds | hardcoded `1e-3`, `1e-6` m/s | `[output]` or `[run]` block |
| Checkpoint interval | `--checkpoint-interval` argv | `[run]` block |

The result is that switching from BP1 to BP2 to BP5 requires writing an entirely new ~900-line driver file. This scales poorly.

> **Implementation timing note:** The full TOML infrastructure (parser integration, `SEASConfig` hierarchy, validation, factories, `seas_driver.cpp`, regression testing against old drivers) is a substantial standalone effort. It is **not on the critical path to BP5**. The BP5 elasticity operator, 3D fault detection, and multi-component slip are the hard problems. Consider deferring the full input file system until after BP5 is working, when the full shape of `SEASConfig` is known from experience with 3+ benchmarks. Designing the config struct now risks rework when BP5 reveals parameters not yet anticipated (e.g., 3D mesh partitioning options, multi-component output fields, fault orientation vectors).
>
> **Lighter-weight interim approach:** To reduce the 95% code duplication between `bp1_verification_full.cpp` and `bp1_bdrload.cpp` immediately, extract the shared ~800 lines into a `bp1_common.hpp` with template/callback hooks for the differing parts (domain operator type, BC setup). This achieves the deduplication goal without the TOML dependency.
>
> If the TOML system is built first, prefer `toml++` over `toml11` — it is a single header, TOML v1.0 compliant, C++17, and more actively maintained.

### B.2 Reference: Tandem's .toml + .lua Input Design

Tandem solves this problem cleanly. Its `bp1.toml` has seven lines:

```toml
final_time = 94608000000
mesh_file  = "bp1.msh"
mode       = "QD"
type       = "poisson"
lib        = "bp1.lua"
scenario   = "bp1"
ref_normal = [-1, 0]
boundary_linear = true
```

The `type` field selects the physics (`"poisson"` = antiplane, `"elasticity"` = full 3D). The `lib`/`scenario` fields point to a Lua file that provides all parameter functions (`a(x,y)`, `mu(x,y)`, `boundary(x,y,t)`) as closures — no recompilation needed to change parameters.

This gives Tandem a single binary that runs BP1, BP2, BP3, BP5 by swapping the input file. The MFEM implementation should adopt the same philosophy.

### B.3 Proposed SEAS Input File Format

Use **TOML** (same as Tandem, already in the ecosystem) with a consistent block structure that works identically for antiplane 2D, antiplane 3D, and full elasticity 3D. Each block maps directly to one C++ struct and one operator in the code.

```toml
# ---- SEAS simulation input file ----
# Works for: antiplane_2d, antiplane_3d, elasticity_3d

[mesh]
file  = "bp1/mesh/bp1_ss_100m.msh"
scale = 1000.0        # km -> m

[domain]
physics = "antiplane"          # "antiplane" | "plane_strain" | "elasticity_3d"
dg_method = "BR2"              # "BR2" | "IP"
order = 1

[material]
mu     = 32.038e9              # Pa
lambda = 32.038e9              # Pa (ignored for antiplane)
rho    = 2670.0                # kg/m^3 (for radiation damping)

[loading]
Vp = 1.0e-9                   # m/s  plate velocity
Wf = 40.0e3                   # m    rate-state fault depth

[fault]
physical_tag = 3               # Gmsh Physical Curve/Surface tag for fault
components   = 1               # 1 for antiplane, 2 for strike+dip (BP5)

[boundary_conditions]
type = "farfield_dirichlet"    # "farfield_dirichlet" | "natural"
# farfield_dirichlet: u(x=±Lx) = ±Vp/2*t  (Tandem-style, removes null space)
# natural: all-Neumann outer BCs (legacy, may fail for structured quad meshes)

[friction]
law   = "dieterich_ruina_aging"
V0    = 1.0e-6
f0    = 0.6
b     = 0.015
Dc    = 8.0e-3                 # m

# Spatially varying a(z): piecewise linear vs depth
[[friction.a_segments]]
depth_km = [0.0,  15.0]
a        = [0.010, 0.010]

[[friction.a_segments]]
depth_km = [15.0, 18.0]
a        = [0.010, 0.025]      # linear ramp

[[friction.a_segments]]
depth_km = [18.0, 9999.0]
a        = [0.025, 0.025]

# NOTE: The a_segments piecewise-linear format works for BP1/BP2 where a(z) depends
# only on depth. BP5's a(x2, x3) depends on two spatial coordinates with a smooth
# function (SCEC spec Eq. 14), which cannot be expressed as depth segments.
# Options for BP5:
#   (a) a_function = "bp5" — dispatches to a compiled function in bp5_params.hpp
#   (b) Lua callback (like Tandem) — adds a dependency but maximum flexibility
#   (c) Separate [friction.a_2d] block with analytic function parameters
# For now, a_segments is BP1/BP2-only. BP5 will need a different mechanism.

[time_integration]
solver    = "dormand_prince_rk45"
abs_tol   = 1.0e-7
rel_tol   = 0.0                # pure absolute tolerance (Tandem convention)
dt_min    = 1.0e-6             # s
dt_max    = 1.577e7            # s = 0.5 yr
dt_init   = 1.0e3              # s
t_final   = 9.461e10           # s = 3000 yr

[run]
max_steps           = 10000000
checkpoint_interval = 5000
earthquake_V_on     = 1.0e-3   # m/s  threshold to declare earthquake
earthquake_V_off    = 1.0e-6   # m/s  threshold to end earthquake

[output]
dir          = "bp1/results_ss_100m"
prefix       = "bp1_ss_100m"
ref_dir      = "bp1/benchmark_data"
probe_depths_km = [0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5,
                   20.0, 25.0, 30.0, 35.0]
```

### B.4 Single Universal Driver Binary

With the above input format, a single `seas_driver` binary replaces all of `seas_bp1_full`, `seas_bp1_bdrload`, `seas_bp2`, `seas_bp5`, etc.:

```
ibrun ./seas_driver --input bp1/input/bp1_ss_100m.toml
```

The `main()` function reduces to ~50 lines:

```cpp
int main(int argc, char* argv[]) {
    MPIContext mpi(&argc, &argv);
    SEASConfig cfg = SEASConfig::Load(argv[1]);   // parse TOML

    auto mesh   = MeshLoader::Load(cfg.mesh, mpi);
    auto domain = DomainFactory::Create(cfg, *mesh); // dispatches on cfg.domain.physics
    auto fault  = FaultFactory::Create(cfg, *domain, mpi);
    auto seas   = SEASFactory::Create(cfg, *domain, *fault, mpi);
    auto output = OutputFactory::Create(cfg, *fault, mpi);

    SEASTimeLoop(cfg, seas, output, mpi);         // generic RK45 loop
    return 0;
}
```

**`DomainFactory::Create`** dispatches:
- `"antiplane"` + `"natural"`         → `AntiplaneDomainOperator`
- `"antiplane"` + `"farfield_dirichlet"` → `AntiplaneBdrLoadOperator`
- `"elasticity_3d"` + `"farfield_dirichlet"` → `ElasticityDomainOperator` (BP5)

This eliminates the code duplication between `bp1_verification_full.cpp` and `bp1_bdrload.cpp` and makes adding BP5 as simple as implementing `ElasticityDomainOperator` — the driver, time loop, output, and fault operator are shared.

### B.5 Implementation Steps for the Input File System

**Step B1 (before Phase 1 in the BP5 plan):** Define `SEASConfig` struct hierarchy in `config/seas_config.hpp`:
- `MeshConfig`, `DomainConfig`, `MaterialConfig`, `LoadingConfig`, `FaultConfig`, `BCConfig`, `FrictionConfig`, `TimeConfig`, `RunConfig`, `OutputConfig`
- Use a lightweight TOML parser (e.g., `toml++` header-only, or `toml11`)

**Step B2:** Implement `SEASConfig::Load(filename)` which reads the TOML and populates the structs. Provide validation (missing required fields → clear error messages).

**Step B3:** Implement factory functions `DomainFactory::Create`, `FaultFactory::Create`, `OutputFactory::Create`. These replace the inline construction logic currently in `main()`.

**Step B4:** Write a single `seas_driver.cpp` using the factory pattern. Keep the old `bp1_verification_full.cpp` and `bp1_bdrload.cpp` building for regression testing until the new driver is validated.

**Step B5:** Convert all `bp1_*.sbatch` and `bp1_ss_*.sbatch` jobs to use `seas_driver --input <file>.toml`.

---

## Part C: Fault Detection via Gmsh Physical Groups

### C.1 Why Coordinate-Based Detection Must Be Replaced

The current `IsFaultFace()` in `antiplane_operator.hpp` uses:

```cpp
bool IsFaultFace(int face) const {
    real_t x, z;
    GetFaceCenter(face, x, z);
    const real_t tol = 1e-10 * std::max(Wf_, 1.0);
    return std::abs(x) < tol;  // fault at x = 0
}
```

This fails for:
1. **Dipping faults** (BP3: 60° dip — fault plane not aligned with any coordinate axis)
2. **Curved fault traces** (BP6, real earthquake geometries)
3. **Multiple fault segments** (conjugate faults, fault systems)
4. **Offset faults** (fault not passing through origin)
5. **3D fault surfaces** (BP5: fault is a plane in 3D, not a line)

It also has a correctness issue: the tolerance `1e-10 * Wf_` is in meters (after km→m scaling), making it ~4×10⁻⁶ m, which may miss fault faces if element centers aren't exactly at x=0 in floating-point arithmetic.

### C.2 Tandem's Approach: Physical Group Tag as BC

Tandem maps Gmsh `Physical Curve/Surface` tags directly to the `BC` enum at mesh load time (in `GlobalSimplexMeshBuilder.cpp`):

```
Physical Curve(1)  → BC::None
Physical Curve(3)  → BC::Natural  (free surface, bottom)
Physical Curve(5)  → BC::Fault    (rate-state fault)
Physical Curve(7+) → BC::Dirichlet (far-field loading)
```

The `BC::Fault` tag (integer 3) is written in the Gmsh `.geo` file by the user. The code never checks coordinates — it only checks `face_bc == BC::Fault`. Adding a new fault geometry requires only editing the `.geo` file, not the C++ source.

The data structure is a flat array `std::vector<BC> boundaryConditions` indexed by local face number, populated once during mesh loading. `SetupFaultInfo()` iterates over all faces and appends those with `BC::Fault` to the fault face list.

### C.3 Proposed MFEM Implementation

**Step 1: Add fault tag to Gmsh `.geo` files.**

In `bp1.geo` and `bp1_selfsimilar.geo`, add a `Physical Curve` for the fault:

```geo
// Fault attribute tag (must match SEAS_FAULT_TAG in C++ code)
Physical Curve(3) = {7, 8, 9, 10};    // fault line segments (rate-state zone + below-Wf)
```

This is a one-line addition to each `.geo` file. The fault lines are the interior edges at x=0.

For BP5 in 3D (`bp5.geo`):
```geo
Physical Surface(3) = {fault_surface};  // 2D fault plane at x1=0
```

**Step 2: Create `FaultBoundaryData` — face-to-tag map.**

Add a new class `FaultBoundaryData` (or extend `BP2MeshGenerator::LoadGmshMesh`) that, after loading the Gmsh mesh, builds a mapping from face index to fault attribute:

```cpp
// In bp2_mesh.hpp or a new fault_mesh.hpp
class FaultBoundaryData {
public:
    // After loading the Gmsh mesh, scan boundary elements for the fault tag.
    // Interior faces adjacent to a boundary element with fault_tag are fault faces.
    // NOTE: Return mfem::Array<int> (not std::vector<int>) for consistency with
    // the existing fault_interior_faces_ and fault_shared_faces_ members.
    static Array<int> FindFaultInteriorFaces(
        const Mesh& mesh, int fault_tag);
    static Array<int> FindFaultSharedFaces(
        const ParMesh& pmesh, int fault_tag);
};
```

**Step 3: Replace `IsFaultFace` with tag-based lookup in `SetupFaultInfo`.**

Current code (to be replaced):
```cpp
// OLD: coordinate check
for (int f = 0; f < num_faces; f++) {
    if (GetInteriorFaceTransformations(f) && IsFaultFace(f))
        fault_interior_faces_.Append(f);
}
```

New code (tag-based):
```cpp
// NEW: physical group tag check
fault_interior_faces_ = FaultBoundaryData::FindFaultInteriorFaces(mesh_, fault_tag_);
fault_shared_faces_   = FaultBoundaryData::FindFaultSharedFaces(pmesh_, fault_tag_);
```

The `fault_tag_` is passed at construction time (from TOML config if available, or as a constructor argument defaulting to `FAULT_TAG = 5`).

**Step 4: How MFEM accesses Gmsh interior boundary elements.**

When MFEM reads a Gmsh `.msh` v2.2 file, `Physical Curve` elements that lie on interior edges (not on the domain boundary) are loaded as boundary elements with `BdrAttribute = tag`. They appear in `mesh.GetNBE()` alongside the true outer boundary elements.

> **Gmsh format version note:** The above behavior is confirmed for Gmsh `.msh` v2.2 format. Behavior may differ for v4.1 format. The `.geo` files should explicitly specify `Mesh.MshFileVersion = 2.2;` to ensure consistent behavior, or the implementation should be tested with both formats before relying on v4.1.

To find the interior face corresponding to a boundary element with a given tag:

```cpp
// Scan boundary elements for the fault tag
for (int be = 0; be < mesh.GetNBE(); be++) {
    if (mesh.GetBdrAttribute(be) != fault_tag) continue;
    // Get the face index this boundary element corresponds to
    int face, ori;
    mesh.GetBdrElementFace(be, &face, &ori);
    // Verify it is an interior face (has two adjacent elements)
    auto* FTr = mesh.GetInteriorFaceTransformations(face);
    if (FTr != nullptr)
        fault_interior_faces_.Append(face);
}
```

This is the exact approach Tandem uses, adapted to MFEM's mesh API.

### C.4 Required Changes to `.geo` Files

The two existing `.geo` files need the fault tagged explicitly. Note that in the current full-domain `.geo` files, the fault line segments (the edges at x=0) are already defined but not assigned to a Physical Curve. They just exist as geometric lines shared between the two Physical Surfaces. We need to add:

**`bp1.geo`** — add after the existing Physical groups:
```geo
// Fault interface (interior, rate-state zone: z=0 to z=-Wf)
Physical Curve(3) = {7, 8, 9, 10};    // Lines at x=0, z = [0, -d4]

// Note: Line 11 (grading zone, z=[-d4,-d4-d5]) and Line 12 (below, z=[-d4-d5,-D])
// are the below-Wf region — include them in fault tag so the code can apply Vp there
// OR keep a separate tag for below-Wf if the code distinguishes rate-state vs locked
Physical Curve(3) = {7, 8, 9, 10, 11, 12};   // full fault including below-Wf
```

> **Recommended approach: single fault tag, programmatic RS/locked distinction.** Tag the entire fault interface (rate-state zone + below-Wf) with one Physical Curve tag. Determine rate-state vs. locked programmatically from `a(z) < b` (velocity-weakening condition), as currently done in `RateStateFaultOperator::ComputeRHS()`. Do not use separate Gmsh tags for RS vs. locked — this couples mesh generation to physics parameters and breaks when `a(z)` changes. For BP5, the RS zone boundary is defined by a smooth `a(x2, x3)` function, not a simple depth cutoff, making tag-based distinction even less practical.

**`bp1_selfsimilar.geo`** — add:
```geo
// Fault interface (Lines 17 and 18: x=0, z=0 to z=-Zf, then z=-Zf to z=-D)
Physical Curve(3) = {17, 18};          // fault at x=0
```

### C.5 Fault Tag Convention

Following Tandem's BC enum mapping, adopt this convention in the MFEM code:

```cpp
// In a new file: domain/seas_boundary_tags.hpp
struct SEASBoundaryTags {
    static constexpr int FARFIELD_LEFT  = 1;   // existing, unchanged
    static constexpr int FARFIELD_RIGHT = 2;   // existing, unchanged
    static constexpr int FREE_SURFACE   = 3;   // existing (top, z=0)  ← CONFLICT
    static constexpr int BOTTOM         = 4;   // existing, unchanged
    static constexpr int FAULT          = 5;   // NEW: fault physical group
};
```

**Note on the FREE_SURFACE conflict:** The current code uses tag 3 for `FREE_SURFACE` (top boundary, z=0). Tandem uses tag 3 for `BC::Natural` (both free surface and bottom) and tag 5 for `BC::Fault`. To avoid ambiguity when the free-surface tag (3) and the Tandem fault tag (5) would be confused, adopt tag `5` for `FAULT` and keep tags 1–4 for the outer boundaries.

**Interim solution (before TOML input system):** Hard-code `FAULT = 5` as a compile-time constant in `seas_boundary_tags.hpp` and pass it as a constructor argument defaulting to 5. This decouples the fault tag from the TOML input system. If the TOML system is built later, `[fault] physical_tag = 5` makes the tag user-configurable.

The `.geo` files must then use `Physical Curve(5)` (or whatever tag the user specifies) for the fault:

```geo
Physical Curve(5) = {7, 8, 9, 10, 11, 12};   // FAULT (matches Tandem's Curve 5)
```

### C.6 Impact on Existing Code

The coordinate-based `IsFaultFace` can be kept as a **fallback** for backward compatibility when no fault physical tag is specified. The `AntiplaneDomainOperator` constructor gains an optional `int fault_tag = -1` parameter:

- `fault_tag == -1`: use legacy coordinate check (current behavior, no `.geo` change needed)
- `fault_tag >= 1`: use physical group tag (new behavior, requires updated `.geo` files)

This means existing sbatch jobs and meshes continue working without any changes, while new BP5 work uses the tag-based approach from the start.

### C.7 Generalizing to Complex Fault Geometries

With the physical-group approach, complex fault geometries require only Gmsh file changes:

| Geometry | `.geo` change | C++ change |
|---|---|---|
| Dipping fault (BP3 60°) | Tag the dipped fault surface with `Physical Curve(5)` | None |
| Two fault segments | Tag both with `Physical Curve(5)` | None |
| Curved fault trace | Tag with `Physical Curve(5)`, mesh normally | None |
| Multiple faults (e.g. main + splay) | Tag with separate tags, e.g. `Physical Curve(5)` and `Physical Curve(6)` | Pass `fault_tags = {5, 6}` to constructor |
| 3D fault surface (BP5) | `Physical Surface(5)` on the fault plane | Switch from `Physical Curve` to `Physical Surface` in the finder |

The only C++ generalization needed beyond the tag check is allowing `fault_tags` to be a list rather than a single integer, which enables multi-fault simulations.

---

## Part D: Revised Phase 1: Preparatory Refactoring (Updated)

Phase 1 from the BP5 implementation plan is now expanded to include fault detection, interface cleanup, and driver deduplication. These are the critical-path prerequisites for BP5.

> **Revised ordering rationale:** The original ordering (1a: TOML → 1b: fault tags → 1c: interface) placed the TOML input system first. The revised ordering below puts fault detection and interface cleanup first, since these directly unblock BP5 implementation. The TOML input system is deferred to a later phase (see Part B timing note), because:
> - Fault tag detection is a small, self-contained change with immediate value — it unblocks BP3/BP5 mesh design.
> - Interface cleanup is also small and directly enables the BP5 elasticity operator.
> - The TOML infrastructure is large, not on the critical path, and its design benefits from knowing what BP5 actually needs.

### Phase 1a: Fault Physical Group Detection (PRIORITY — do first)
1. Add `Physical Curve(5)` tags to `bp1.geo` and `bp1_selfsimilar.geo` for the fault edges.
2. Define `SEASBoundaryTags` in a new `domain/seas_boundary_tags.hpp` with `FAULT = 5`.
3. Implement `FaultBoundaryData::FindFaultInteriorFaces()` and `FindFaultSharedFaces()` in `domain/seas_boundary_tags.hpp`.
4. Update `AntiplaneDomainOperator::SetupFaultInfo()` to accept `fault_tag` and use tag-based lookup when `fault_tag >= 1`.
5. Keep `fault_tag = -1` default for backward compatibility (triggers legacy coordinate check).
6. Run full BP1 verification to confirm the tag-based fault detection produces identical results to the coordinate-based approach on the existing meshes.

### Phase 1b: `DomainOperator` Interface Cleanup
1. Add `NumSlipComponents()` virtual method (returns 1 for antiplane, 2 for BP5).
2. Add `GetFaultCoords2D()` virtual method for 2D fault coordinate sets.
3. Generalize `SEASQuasiDynamicOperator` to use `DomainOperator<MeshType>*` (abstract base) instead of hardcoded `AntiplaneDomainOperator<MeshType>`.
4. Write `test_domain_operator_interface.cpp` to verify the interface.

### Phase 1c: Driver Deduplication
1. Extract shared code from `bp1_verification_full.cpp` and `bp1_bdrload.cpp` (~800 shared lines) into `bp1_common.hpp` with template/callback hooks for the differing parts (domain operator type, BC setup).
2. Refactor both drivers to use the shared code.
3. Verify identical output from both refactored drivers.

### Phase 1d: Input File Infrastructure (DEFERRED — do after BP5 works)
1. Add `toml++` as a header-only dependency in the build system.
2. Define `SEASConfig` hierarchy in `config/seas_config.hpp`.
3. Implement `SEASConfig::Load(filename)` with validation.
4. Implement `DomainFactory::Create`, `FaultFactory::Create`, `OutputFactory::Create`.
5. Write `seas_driver.cpp` using the factories. Build alongside existing drivers initially.
6. Write `.toml` input files for all existing BP1 runs (100m, 50m, 25m, bdrload variants, self-similar variants).
7. Validate that `seas_driver` produces identical output to `seas_bp1_full` and `seas_bp1_bdrload`.

The Phase 1a→1b→1c ordering ensures: fault detection comes first (so BP5 `.geo` design is clean from the start), interface cleanup comes second (so the elasticity operator can be implemented against a stable, finalized interface), and driver deduplication comes third (reducing maintenance burden before adding more benchmarks). The TOML input system (1d) is deferred until BP5 reveals the full parameter space.

---

## Part E: Unit Test Design

Every new infrastructure piece introduced in Phases 1a, 1b, and 1c must have unit tests that verify correctness independently of the full SEAS simulation. This section specifies the test files, test cases, inputs, and pass criteria for each component.

The project uses a lightweight custom test framework (no external library). The pattern established in `tests/unit/test_fault_operator.cpp` and `tests/unit/test_antiplane.cpp` is:

```cpp
static int num_tests  = 0;
static int num_passed = 0;
static int num_failed = 0;

#define TEST_ASSERT(condition, message) \
   do { num_tests++; \
        if (!(condition)) { std::cerr << "FAILED: " << message \
                             << " (line " << __LINE__ << ")\n"; num_failed++; } \
        else { std::cout << "  PASSED: " << message << "\n"; num_passed++; } \
   } while (0)

#define TEST_NEAR(value, expected, tol, message) \
   do { num_tests++; \
        real_t _v = (value), _e = (expected), _t = (tol); \
        if (std::abs(_v - _e) > _t) { \
            std::cerr << "FAILED: " << message << " got=" << _v \
                      << " expected=" << _e << " tol=" << _t \
                      << " (line " << __LINE__ << ")\n"; num_failed++; } \
        else { std::cout << "  PASSED: " << message << "\n"; num_passed++; } \
   } while (0)
```

All test files are placed in `miniapps/seas/tests/unit/` and added to the `Makefile` alongside existing tests.

---

### E.1 `test_seas_config.cpp` — TOML Input File Parsing

> **Dependency note:** This test depends on the TOML input system (Phase 1d). If Phase 1d is deferred (see Part D), this test is also deferred. Prioritize E.3 (fault detection) and E.4 (domain interface) which are on the critical path.

**Purpose**: Verify that `SEASConfig::Load()` correctly parses valid `.toml` files, reports errors on invalid files, and populates all struct fields.

**File location**: `miniapps/seas/tests/unit/test_seas_config.cpp`

**Dependencies**: `config/seas_config.hpp`, `toml++` header (no MFEM required for these tests).

#### Test E.1.1 — Load minimal valid antiplane config

```
Input file (written inline as a string, not from disk):
  [mesh]
  file = "bp1/mesh/bp1_ss_100m.msh"
  [domain]
  physics = "antiplane"
  loading = "farfield"
  [material]
  mu = 32.04e9
  [fault]
  Wf = 15000.0
  [friction]
  type = "dieterich_ruina"
  [time_integration]
  dt_min = 0.01
  dt_max = 3.15e7
  t_final = 3.15e9
  [output]
  dir = "results_test"
  interval = 100

Pass criteria:
  TEST_ASSERT(cfg.mesh.file == "bp1/mesh/bp1_ss_100m.msh")
  TEST_ASSERT(cfg.domain.physics == PhysicsType::Antiplane)
  TEST_ASSERT(cfg.domain.loading == LoadingType::FarField)
  TEST_NEAR(cfg.material.mu, 32.04e9, 1e3)
  TEST_NEAR(cfg.fault.Wf, 15000.0, 1e-6)
  TEST_ASSERT(cfg.friction.type == FrictionType::DieterichRuina)
  TEST_NEAR(cfg.time.dt_min, 0.01, 1e-10)
  TEST_NEAR(cfg.time.dt_max, 3.15e7, 1.0)
  TEST_NEAR(cfg.time.t_final, 3.15e9, 1.0)
  TEST_ASSERT(cfg.output.dir == "results_test")
  TEST_ASSERT(cfg.output.interval == 100)
```

#### Test E.1.2 — Missing required field reports error

```
Input: same as E.1.1 but without [material] block.

Pass criteria:
  TEST_ASSERT(SEASConfig::Load(str).has_error() == true)
  (or: Load() throws std::runtime_error with a message mentioning "material.mu")
```

#### Test E.1.3 — Invalid physics type reports error

```
Input: [domain] block with physics = "antiplane_typo"

Pass criteria:
  TEST_ASSERT(Load(str).has_error() == true)
  (error message contains "unknown physics type")
```

#### Test E.1.4 — Default values populated when optional fields absent

```
Input: valid config without [fault] fault_tag field.

Pass criteria:
  TEST_ASSERT(cfg.fault.fault_tag == -1)   // -1 means: use coordinate fallback
```

#### Test E.1.5 — Friction segment list parsed correctly

```
Input: [friction] block with segment list:
  [[friction.segments]]
  depth_min = 0.0
  depth_max = 15000.0
  a = 0.010
  [[friction.segments]]
  depth_min = 15000.0
  depth_max = 20000.0
  a = 0.025

Pass criteria:
  TEST_ASSERT(cfg.friction.segments.size() == 2)
  TEST_NEAR(cfg.friction.segments[0].a, 0.010, 1e-6)
  TEST_NEAR(cfg.friction.segments[1].depth_min, 15000.0, 1e-3)
```

#### Test E.1.6 — Round-trip: save then reload produces identical struct

```
Setup:
  SEASConfig cfg1 = SEASConfig::Load(valid_toml_string);
  cfg1.Save("test_roundtrip.toml");
  SEASConfig cfg2 = SEASConfig::Load("test_roundtrip.toml");

Pass criteria:
  TEST_ASSERT(cfg1.mesh.file == cfg2.mesh.file)
  TEST_NEAR(cfg1.material.mu, cfg2.material.mu, 1e-3)
  TEST_ASSERT(cfg1.domain.physics == cfg2.domain.physics)
  TEST_ASSERT(cfg1.friction.segments.size() == cfg2.friction.segments.size())
  // clean up: remove test_roundtrip.toml
```

---

### E.2 `test_domain_factory.cpp` — Factory Dispatch Correctness

> **Dependency note:** This test depends on `DomainFactory` (Phase 1d). Deferred along with Phase 1d.

**Purpose**: Verify that `DomainFactory::Create()` instantiates the correct operator subclass for each (physics, loading) combination, and that the returned object reports correct metadata.

**File location**: `miniapps/seas/tests/unit/test_domain_factory.cpp`

**Dependencies**: `config/seas_config.hpp`, `domain/antiplane_operator.hpp`, `domain/antiplane_bdrload_operator.hpp`, `domain/domain_factory.hpp`, MFEM (for mesh).

**Test mesh**: A small programmatic 2D mesh is built inline (no file I/O) using `mfem::Mesh::MakeCartesian2D` with appropriate boundary attributes.

#### Test E.2.1 — antiplane + natural → AntiplaneDomainOperator

```cpp
SEASConfig::DomainConfig dom_cfg;
dom_cfg.physics = PhysicsType::Antiplane;
dom_cfg.loading = LoadingType::Natural;

auto op = DomainFactory::Create(dom_cfg, mesh, params);

TEST_ASSERT(dynamic_cast<AntiplaneDomainOperator*>(op.get()) != nullptr)
TEST_ASSERT(op->NumSlipComponents() == 1)
TEST_ASSERT(op->Dimension() == 2)
```

#### Test E.2.2 — antiplane + farfield → AntiplaneBdrLoadOperator

```cpp
dom_cfg.loading = LoadingType::FarField;
auto op = DomainFactory::Create(dom_cfg, mesh, params);
TEST_ASSERT(dynamic_cast<AntiplaneBdrLoadOperator*>(op.get()) != nullptr)
TEST_ASSERT(op->NumSlipComponents() == 1)
```

#### Test E.2.3 — Unknown physics type throws

```cpp
dom_cfg.physics = PhysicsType::Unknown;
TEST_ASSERT(throws_exception([&]{ DomainFactory::Create(dom_cfg, mesh, params); }))
```

#### Test E.2.4 — Factory respects fault_tag field

```cpp
dom_cfg.physics = PhysicsType::Antiplane;
dom_cfg.loading = LoadingType::Natural;
dom_cfg.fault_tag = 3;   // pass tag through config
auto op = DomainFactory::Create(dom_cfg, mesh, params);
// Access the operator's internal fault_tag_ (via accessor or friend)
TEST_ASSERT(op->FaultTag() == 3)
```

---

### E.3 `test_fault_detection.cpp` — Tag-Based vs Coordinate-Based Fault Identification

**Purpose**: Verify that `FaultBoundaryData::FindFaultFaces()` (tag-based) and the legacy `IsFaultFace()` (coordinate-based) return **identical sets** of interior face indices on meshes where both methods are applicable.

**File location**: `miniapps/seas/tests/unit/test_fault_detection.cpp`

**Dependencies**: `domain/seas_boundary_tags.hpp`, `domain/antiplane_operator.hpp`, MFEM.

#### Test setup — programmatic minimal mesh

To avoid file I/O, the test builds a small mesh programmatically:

```cpp
// 2D Cartesian mesh: x in [-0.5, 0.5], z in [-1.0, 0.0], 4x4 elements
// The interior faces at x=0 are the "fault"
mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(4, 4, mfem::Element::QUADRILATERAL,
                                               false, 1.0, 1.0);
// Shift so x in [-0.5, 0.5]
{
    for (int i = 0; i < mesh.GetNV(); i++) {
        double* v = mesh.GetVertex(i);
        v[0] -= 0.5;
    }
}
// Assign boundary attribute 5 to all boundary faces on x = 0.0
// (In the real mesh, x=0 is an interior fault interface, not a boundary.
//  We test using a tiny 2-column mesh where the left half and right half
//  share interior faces exactly at x=0.)
```

For the fault detection test, a better mesh is two side-by-side rectangles sharing a common edge, where the shared edge is tagged as `Physical Curve 5` (fault tag):

```cpp
// Build mesh with explicit interior face tags using Mesh::AddBdrElement
// OR: use a minimal .msh file embedded as a string constant in the test
// The .msh string describes a 2x1 grid with a tagged interior edge at x=0
```

#### Test E.3.1 — Tag-based finder returns only interior faces

```
Setup: mesh with fault_tag=5 assigned to interior faces at x=0.

auto faces = FaultBoundaryData::FindFaultFaces(mesh, 5);

TEST_ASSERT(faces.size() > 0)
for each face index f in faces:
  TEST_ASSERT(mesh.FaceIsInterior(f))
```

#### Test E.3.2 — Coordinate-based fallback returns same set

```
Setup: same mesh.

auto faces_tag   = FaultBoundaryData::FindFaultFaces(mesh, 5);
auto faces_coord = FaultBoundaryData::FindFaultFacesByCoord(mesh, 0.0, 1e-8);

TEST_ASSERT(faces_tag.size() == faces_coord.size())
// Check sets are equal (sort both, then compare)
for each i:
  TEST_ASSERT(faces_tag[i] == faces_coord[i])
```

#### Test E.3.3 — No fault faces returned when tag absent

```
Setup: mesh where no interior face has attribute 5.

auto faces = FaultBoundaryData::FindFaultFaces(mesh, 5);
TEST_ASSERT(faces.size() == 0)
```

#### Test E.3.4 — fault_tag = -1 triggers coordinate fallback in operator

```
Setup: AntiplaneDomainOperator constructed with fault_tag = -1.

// After SetupFaultInfo(), the operator must still find the correct fault faces
// using the legacy IsFaultFace() coordinate check.
auto op = AntiplaneDomainOperator(mesh, params, /*fault_tag=*/-1);
op.SetupFaultInfo();
TEST_ASSERT(op.NumFaultFaces() > 0)
```

#### Test E.3.5 — tag-based and coordinate-based produce same NumFaultFaces on bp1_100m mesh

```
Setup: load bp1/mesh/bp1_ss_100m.msh (if available at test time, else SKIP).

auto op_coord = AntiplaneDomainOperator(mesh, params, /*fault_tag=*/-1);
auto op_tag   = AntiplaneDomainOperator(mesh, params, /*fault_tag=*/5);

TEST_ASSERT(op_coord.NumFaultFaces() == op_tag.NumFaultFaces())
```

This regression test guards against any future refactoring changing the fault face count.

---

### E.4 `test_domain_operator_interface.cpp` — `DomainOperator` Base Interface

**Purpose**: Verify the virtual interface functions (`NumSlipComponents()`, `Dimension()`, `GetFaultCoords2D()`) return physically correct values for the antiplane operator.

**File location**: `miniapps/seas/tests/unit/test_domain_operator_interface.cpp`

**Dependencies**: `domain/antiplane_operator.hpp`, `domain/domain_operator.hpp` (base class), MFEM.

#### Test E.4.1 — NumSlipComponents for antiplane = 1

```cpp
auto op = make_antiplane_operator(small_mesh, params);
TEST_ASSERT(op->NumSlipComponents() == 1)
```

#### Test E.4.2 — Dimension for antiplane = 2

```cpp
TEST_ASSERT(op->Dimension() == 2)
```

#### Test E.4.3 — GetFaultCoords2D returns sorted depth array

```
Setup: antiplane operator on a mesh with N fault faces.

auto coords = op->GetFaultCoords2D();

TEST_ASSERT(coords.Size() == op->NumFaultFaces())
// Check monotonically non-decreasing (sorted by depth)
for i in 1..coords.Size()-1:
  TEST_ASSERT(coords[i] >= coords[i-1])
// Check range: all depths in [0, Wf + delta]
for i in 0..coords.Size()-1:
  TEST_ASSERT(coords[i] >= 0.0)
  TEST_ASSERT(coords[i] <= params.Wf + 1.0)
```

#### Test E.4.4 — Polymorphic dispatch works via base pointer

```cpp
std::unique_ptr<DomainOperator> op = DomainFactory::Create(cfg, mesh, params);

// Call via base pointer — must not crash or return garbage
TEST_ASSERT(op->NumSlipComponents() >= 1)
TEST_ASSERT(op->Dimension() >= 2)
```

---

### E.5 `test_toml_config_integration.cpp` — End-to-End Config → Operator Construction

**Purpose**: Verify that reading a `.toml` file and constructing an operator via the factory chain produces an operator that passes all interface sanity checks. This integrates E.1, E.2, and E.4.

**File location**: `miniapps/seas/tests/unit/test_toml_config_integration.cpp`

#### Test E.5.1 — Full pipeline: toml string → SEASConfig → DomainOperator

```
1. Write a minimal valid .toml string for antiplane/farfield (no disk file needed).
2. SEASConfig cfg = SEASConfig::LoadFromString(toml_string);
3. mfem::Mesh mesh = build_small_test_mesh();
4. auto op = DomainFactory::Create(cfg.domain, mesh, cfg.material);

TEST_ASSERT(op != nullptr)
TEST_ASSERT(op->NumSlipComponents() == 1)
TEST_ASSERT(op->Dimension() == 2)
```

#### Test E.5.2 — Invalid config does not construct operator

```
1. Write an invalid .toml (missing mu field).
2. TEST_ASSERT(SEASConfig::LoadFromString(bad_toml).has_error())
3. Confirm DomainFactory::Create is never called (exception before factory).
```

---

### E.6 `test_fault_detection_3d.cpp` — 3D Fault Surface Detection (Phase 3 prerequisite)

**Purpose**: Verify that the 3D generalization of `FindFaultFaces()` correctly identifies fault surface faces when a 3D mesh has a `Physical Surface` fault tag.

> **Implementation timing:** This test should be written when Phase 3 (3D elasticity operator) begins, not during Phase 1. Writing a test for code that doesn't exist yet creates a broken test target. During Phase 1, create at most a stub file with a comment indicating the planned tests, or gate the test with `#ifdef SEAS_HAS_ELASTICITY3D`. Add the test to the Makefile only when the 3D code path is implemented.

**File location**: `miniapps/seas/tests/unit/test_fault_detection_3d.cpp`

#### Test E.6.1 — 3D mesh with tagged internal face returns non-zero face count

```
Setup: two 3D brick elements sharing a common face; the shared face is tagged
       with boundary attribute 5.

auto faces = FaultBoundaryData::FindFaultFaces3D(mesh, 5);
TEST_ASSERT(faces.size() > 0)
for each face f in faces:
  TEST_ASSERT(mesh.FaceIsInterior(f))
```

#### Test E.6.2 — NumSlipComponents for full elasticity 3D = 2

```
// After Phase 3, ElasticityDomainOperator exists:
auto op = ElasticityDomainOperator(mesh_3d, params_3d);
TEST_ASSERT(op->NumSlipComponents() == 2)
TEST_ASSERT(op->Dimension() == 3)
```

---

### E.7 Test Build and Makefile Integration

All new tests are compiled and run via `make test`:

```makefile
# In miniapps/seas/tests/unit/Makefile (or appended to the existing one):

UNIT_TESTS += test_seas_config
UNIT_TESTS += test_domain_factory
UNIT_TESTS += test_fault_detection
UNIT_TESTS += test_domain_operator_interface
UNIT_TESTS += test_toml_config_integration
UNIT_TESTS += test_fault_detection_3d   # added when Phase 3 starts

test_seas_config: test_seas_config.cpp ../../config/seas_config.hpp
	$(CXX) $(CXXFLAGS) -I../../.. $< -o $@ $(LDFLAGS)
	./$@

test_domain_factory: test_domain_factory.cpp ../../domain/domain_factory.hpp
	$(CXX) $(CXXFLAGS) $(MFEM_FLAGS) $< -o $@ $(MFEM_LIBS) $(LDFLAGS)
	./$@

# ... similar rules for each test
```

Each test binary exits with code 0 on all-pass and code 1 on any failure, enabling CI integration.

---

### E.8 Test Coverage Summary by Phase

| Phase | New infrastructure | Test file | Key assertions | Priority |
|---|---|---|---|---|
| 1a | `FindFaultFaces()` (tag) | `test_fault_detection.cpp` | all faces interior, same count as coord-based, empty on absent tag | **HIGH** — write first |
| 1a | `IsFaultFace()` fallback | `test_fault_detection.cpp` | fault_tag=-1 triggers coord, identical result to tag-based | **HIGH** |
| 1b | `DomainOperator` interface | `test_domain_operator_interface.cpp` | NumSlipComponents==1, Dimension==2, coords sorted and in range | **HIGH** |
| 1d | `SEASConfig::Load()` | `test_seas_config.cpp` | field values, error on bad input, defaults, round-trip | Deferred with 1d |
| 1d | `DomainFactory::Create()` | `test_domain_factory.cpp` | correct subclass, NumSlipComponents, fault_tag wired | Deferred with 1d |
| 1d | Integration | `test_toml_config_integration.cpp` | full pipeline from string→config→operator→interface check | Deferred with 1d |
| 3 | 3D fault surface + elasticity | `test_fault_detection_3d.cpp` | 3D face detection, NumSlipComponents==2 | Write when Phase 3 starts |

**Regression guard**: tests E.3.5 and E.4.3 explicitly compare against the current bp1_100m mesh results, so any accidental change to fault face counting or coordinate ranges will be caught immediately.

**Design rule**: every function added in Phase 1–3 has at least one positive test (valid input → correct output) and one negative test (invalid input → correct error). No new infrastructure is merged until its test file passes with `num_failed == 0`.

---

## Part F: Build and Run Reference

This section is a self-contained quick-reference for building the code, generating meshes, running simulations, and running tests — both locally and on TACC. It requires no prior knowledge of the build system.

---

### F.1 Prerequisites

| Dependency | Version used | Purpose |
|---|---|---|
| MFEM | 4.7+ (built with MPI) | FEM framework — must be built before SEAS |
| HYPRE | 2.31.0 | AMG preconditioner (MFEM dependency) |
| MUMPS | 5.3 | Direct sparse solver |
| ParMETIS | any | Parallel mesh partitioning |
| Gmsh | 4.x | Mesh generation from `.geo` files |
| MPI | OpenMPI or MVAPICH | Parallel execution |

**Check MFEM was built with MPI** (required for all parallel executables):
```bash
grep "MFEM_USE_MPI" ../../config/config.mk
# Should show: MFEM_USE_MPI = YES
```

---

### F.2 Building

All commands are run from the SEAS miniapp directory:
```bash
cd /path/to/seas-mfem/miniapps/seas
```

#### Build a specific executable

| Command | What it builds |
|---|---|
| `make seas_bp1_bdrload` | BP1 far-field Dirichlet loading (main BP1 driver) |
| `make seas_bp1_full` | BP1 all-Natural BCs (legacy, do not use for self-similar meshes) |
| `make seas_bp2_full` | BP2 full multi-cycle simulation |
| `make seas_bp2_verify` | BP2 first-cycle verification |
| `make seas_pseas` | General parallel SEAS driver |

#### Build everything

```bash
make all          # builds all sequential + parallel targets
make -j4 all      # parallel build (use -jN for N cores)
make clean        # remove all .o files and executables
make clean-build  # same as clean
```

#### Build on TACC (inside a job script)

The sbatch scripts rebuild automatically before running:
```bash
make seas_bp1_bdrload    # inside the job, after module loads
```

If you want to build interactively on a login node (no MPI execution):
```bash
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}
make seas_bp1_bdrload
```

---

### F.3 Mesh Generation with Gmsh

Meshes are `.msh` files generated from `.geo` geometry files. The `-setnumber` flag sets a parameter in the `.geo` file (mesh resolution `h` in km).

#### BP1 unstructured triangle mesh

```bash
# h = mesh size in km on the fault (100m = 0.1 km)
gmsh -2 bp1/mesh/bp1.geo -o bp1/mesh/bp1_100m.msh -setnumber hf 0.1
gmsh -2 bp1/mesh/bp1.geo -o bp1/mesh/bp1_50m.msh  -setnumber hf 0.05
gmsh -2 bp1/mesh/bp1.geo -o bp1/mesh/bp1_25m.msh  -setnumber hf 0.025
```

#### BP1 self-similar structured quad mesh

```bash
# h = uniform mesh size in km throughout the domain
gmsh -2 bp1/mesh/bp1_selfsimilar.geo -o bp1/mesh/bp1_ss_100m.msh -setnumber h 0.1
gmsh -2 bp1/mesh/bp1_selfsimilar.geo -o bp1/mesh/bp1_ss_50m.msh  -setnumber h 0.05
gmsh -2 bp1/mesh/bp1_selfsimilar.geo -o bp1/mesh/bp1_ss_25m.msh  -setnumber h 0.025
```

The sbatch scripts check if the mesh exists before generating it:
```bash
if [ ! -f bp1/mesh/bp1_ss_100m.msh ]; then
    gmsh -2 bp1/mesh/bp1_selfsimilar.geo -o bp1/mesh/bp1_ss_100m.msh -setnumber h 0.1
fi
```

You can do the same check locally to avoid re-meshing.

---

### F.4 Running Locally (Laptop / Workstation)

Use `mpirun` (or `mpiexec`) with a small mesh for testing:

```bash
# Generate a coarse mesh first
gmsh -2 bp1/mesh/bp1.geo -o bp1/mesh/bp1_200m.msh -setnumber hf 0.2

# Run with 4 MPI ranks, coarse mesh, short time
mkdir -p bp1/results_test
mpirun -np 4 ./seas_bp1_bdrload \
    --mesh bp1/mesh/bp1_200m.msh \
    --ref-dir bp1/benchmark_data \
    --output-dir bp1/results_test \
    --checkpoint-interval 100
```

To see all command-line options:
```bash
./seas_bp1_bdrload --help
```

---

### F.5 Running on TACC (Frontera / Stampede)

#### Submitting a job

```bash
cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas
sbatch jobs/bp1_ss_100m.sbatch    # self-similar 100m mesh, 1 node, 32 ranks
sbatch jobs/bp1_ss_50m.sbatch     # self-similar 50m mesh,  2 nodes, 80 ranks
sbatch jobs/bp1_ss_25m.sbatch     # self-similar 25m mesh,  6 nodes, 300 ranks
sbatch jobs/bp1_bdrload_100m.sbatch  # unstructured 100m, 1 node, 10 ranks
sbatch jobs/bp1_bdrload_50m.sbatch   # unstructured 50m,  2 nodes, 40 ranks
sbatch jobs/bp1_bdrload_25m.sbatch   # unstructured 25m,  4 nodes, 100 ranks
```

#### Monitoring jobs

```bash
squeue -u $USER              # list your running/pending jobs
squeue -j <jobid>            # status of a specific job
scancel <jobid>              # cancel a job
```

#### Checking output while running

```bash
tail -f bp1_ss_100m_<jobid>.out    # follow stdout in real time
tail -f bp1_ss_100m_<jobid>.err    # follow stderr (errors appear here)
```

#### Key SLURM parameters in the sbatch files

| Parameter | Meaning |
|---|---|
| `#SBATCH -N <n>` | Number of nodes |
| `#SBATCH -n <n>` | Total MPI tasks across all nodes |
| `#SBATCH -t HH:MM:SS` | Wall-clock time limit |
| `#SBATCH -A EAR20006` | TACC allocation account |
| `ibrun ./executable` | TACC's MPI launcher (replaces `mpirun`) |

#### Module loads required on TACC

Every job script (and interactive session) must load these before running:
```bash
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}
```

The `LD_LIBRARY_PATH` line ensures the correct HYPRE version is found before any older bundled version.

---

### F.6 Running Unit Tests

#### Run all sequential unit tests

```bash
cd miniapps/seas
make test
```

This builds and runs:
- `seas_test_friction_law` — friction law (no MFEM mesh needed)
- `seas_test_state_evolution` — aging/slip law ODEs
- `seas_test_antiplane` — DG domain operator
- `seas_test_fault_operator` — fault geometry and rate-state RHS
- `seas_test_quasi_dynamic` — quasi-dynamic coupling
- `seas_test_bp2_short` — short BP2 simulation
- `seas_test_io` — output file writing
- `seas_test_checkpoint` — checkpoint save/restore

#### Run individual test suites

```bash
make test-friction         # friction law only
make test-state-evolution  # state evolution only
make test-antiplane        # DG operator only
make test-fault-operator   # fault operator only
make test-checkpoint       # checkpoint only
```

#### Run parallel unit tests

```bash
make test-parallel         # runs all parallel tests with 2 and 4 MPI ranks
make test-mpi-context      # MPI setup
make test-parallel-domain  # parallel domain operator
make test-br2-consistency  # BR2 stiffness matrix consistency
```

#### New Phase 1 tests (after Phase 1a/1b/1c implementation)

```bash
make test-seas-config           # TOML config parsing
make test-domain-factory        # factory dispatch
make test-fault-detection       # tag-based vs coordinate fault finding
make test-domain-interface      # NumSlipComponents, Dimension, GetFaultCoords2D
make test-config-integration    # end-to-end: toml → config → operator
```

#### Exit codes

Every test binary exits `0` on all-pass, `1` on any failure. The `make test` target stops on first failure. To run all tests regardless of failures:
```bash
make test || true     # continues even if some fail
```

---

### F.7 Viewing Output in ParaView

The simulation writes ParaView files to the `--output-dir` directory. To visualize:

1. Open ParaView
2. `File → Open` → navigate to `bp1/results_100m_bdrload/`
3. Open the `.pvd` file (collection of all time steps)
4. Click `Apply` in the Properties panel
5. Use `Warp By Scalar` filter to exaggerate displacement

For fault-only output (time series of slip velocity at fault nodes):
- Look for `mfem_bp1_*.txt` files in the output directory
- These are plain-text columns: `time, depth, slip_rate, slip, state`
- Use `bp1/visualize_results.py` or `bp1/mesh_convergence.py` for Python plots

---

### F.8 Common Issues and Fixes

| Symptom | Likely cause | Fix |
|---|---|---|
| `MFEM abort: Error during MUMPS numerical factorization` | Singular stiffness matrix (all-Natural BCs on a structured quad mesh) | Use `seas_bp1_bdrload` instead of `seas_bp1_full`; see Section A for diagnosis |
| `error: MFEM library is not built` | `make` cannot find `config/config.mk` | Check `MFEM_DIR` and `MFEM_BUILD_DIR` in the Makefile, or `cd` to the SEAS directory first |
| `gmsh: command not found` | Gmsh not in PATH | `module load gmsh` on TACC, or install locally |
| `ibrun: command not found` | Running on a login node | Submit via `sbatch` or use `idev` for interactive compute node |
| MUMPS memory error (`-8` or `-9`) | Not enough memory per MPI rank | Increase `#SBATCH -n` (more ranks = less memory each), or reduce `ICNTL(14)` |
| Test prints `FAILED:` but exits 0 | Test framework uses `num_failed` counter | Check if the test binary's final line says `0 failed` — if not, some tests silently passed despite logic errors; cross-check the exit code |

---

## Part G: Developer Manual

This section explains the internal architecture of the MFEM SEAS miniapp in enough detail that a developer unfamiliar with the codebase can locate any behavior, trace a bug, or add a new component without reading every file from scratch.

---

### G.1 Architecture Layers

The codebase is organized in five horizontal layers. Each layer depends only on the layers below it — no upward dependencies.

```
┌─────────────────────────────────────────────────────────────────┐
│  Layer 5 — Drivers & Verification                               │
│  pseas.cpp  bp1_bdrload.cpp  bp1_verification_full.cpp  …      │
│  Entry points: parse args, wire components, run time loop       │
├─────────────────────────────────────────────────────────────────┤
│  Layer 4 — Time Integration                                     │
│  solver/seas_operator.hpp          (SEAS coupling operator)     │
│  solver/time_stepper.hpp           (RK45 + adaptive dt)         │
│  solver/seas_bdrload_operator.hpp  (BP1 BC adapter)             │
├─────────────────────────────────────────────────────────────────┤
│  Layer 3 — Fault Physics                                        │
│  fault/rate_state_fault.hpp   (state ODE + stress balance)      │
│  fault/fault_geometry.hpp     (depths, a(z), MPI gather)        │
│  fault/fault_nodes.hpp        (DOF ↔ face mapping, mass matrix) │
│  friction/dieterich_ruina.hpp (friction law + Newton solver)    │
│  friction/state_evolution.hpp (aging/slip law, psi-space)       │
├─────────────────────────────────────────────────────────────────┤
│  Layer 2 — Domain PDE                                           │
│  domain/domain_operator.hpp         (abstract interface)        │
│  domain/antiplane_operator.hpp      (DG-BR2 Laplacian, full BC) │
│  domain/antiplane_bdrload_operator  (DG + far-field Dirichlet)  │
│  integrator/dg_br2_integrator.hpp   (BR2 lifting operator)      │
├─────────────────────────────────────────────────────────────────┤
│  Layer 1 — Infrastructure                                       │
│  config/bp1_params.hpp   config/bp2_params.hpp  (parameters)    │
│  domain/bp2_mesh.hpp     (mesh factory, boundary attributes)    │
│  io/benchmark_output.hpp  io/checkpoint.hpp  …  (output/restart)│
│  common/seas_types.hpp    common/mpi_context.hpp  (MPI, types)  │
│  MFEM library (mesh, FE space, bilinear forms, linear algebra)  │
└─────────────────────────────────────────────────────────────────┘
```

---

### G.2 Directory and File Reference

Every source file in `miniapps/seas/` and its exact role:

#### Top-level
| File | Purpose |
|---|---|
| `pseas.cpp` | Parallel BP2 driver: wires all layers, runs the time loop |
| `Makefile` | Build rules for all executables and test binaries |

#### `config/`
| File | Purpose |
|---|---|
| `bp1_params.hpp` | `MakeBP1Params()` — returns `BP2Params` with BP1 overrides (Dc=8mm, t_final=3000yr); `GetBP1ProbeDepths()` — 15 probe depths at 2.5 km spacing |
| `bp2_params.hpp` | `BP2Params` struct: all physical constants (mu, rho, cs, Vp, sigma_n, V0, f0, b, Dc, a0, amax, H, h, Wf, Lx, Lz, t_final); `a_of_z(z)` and `tau0()` methods |

#### `common/`
| File | Purpose |
|---|---|
| `seas_types.hpp` | `FESpaceForMesh<MeshType>`, `BilinearFormForMesh<MeshType>`, `GridFunctionForMesh<MeshType>` — conditional type aliases that expand to serial or parallel MFEM classes; `IsParallelMesh<T>` trait |
| `mpi_context.hpp` | `MPIContext` RAII class: init/finalize MPI, expose `Rank()`, `Size()`, `GlobalMax()`, `GlobalSum()`, `Barrier()`, `Bcast()` — all no-ops in serial build |
| `parallel_utils.hpp` | `GatherVectorToRoot()`, `ScatterVectorFromRoot()`, `BroadcastFromRoot()`, `AllRanksAgree()` — variable-length MPI collective helpers |

#### `domain/`
| File | Purpose |
|---|---|
| `domain_operator.hpp` | `DomainOperator<MeshType>` abstract base: declares `Solve()`, `ComputeTraction()`, `GetNumFaultDOFs()`, `GetFaultDepths()`, `GetFaultDOFs()` |
| `antiplane_operator.hpp` | `AntiplaneDomainOperator<MeshType>`: DG-BR2 (or IP) Laplacian with fault at x=0 as interior interface; `SetupFaultInfo()` finds fault faces by `|x| < tol`; assembles stiffness matrix once, solves with MUMPS |
| `antiplane_bdrload_operator.hpp` | `AntiplaneBdrLoadOperator<MeshType>`: same as above but adds Dirichlet BCs `u(±Lx) = ±Vp/2·t` on FARFIELD_LEFT/RIGHT, making matrix SPD |
| `bp2_mesh.hpp` | `BP2MeshGenerator` static factory: `Create()`, `CreateGraded()`, `LoadGmshMesh()`, `SaveVTK()`; `BP2BoundaryAttributes` constants (FARFIELD_LEFT=1, FARFIELD_RIGHT=2, FREE_SURFACE=3, BOTTOM=4) |

#### `fault/`
| File | Purpose |
|---|---|
| `fault_geometry.hpp` | `FaultGeometry<MeshType>`: extracts fault DOF positions, evaluates `a(z)`, `eta`, `sn`; `GatherToRoot()` / `GatherToRootDedup()` — MPI gather with deduplication of DG partition-boundary duplicates |
| `fault_nodes.hpp` | `FaultNodes<MeshType>`: maps fault interior face indices → local DOF pairs; precomputes 2×2 face mass-matrix inverses for L2 projection; `FaceToNodes()`, `FaceLength()`, `FaceMassInverse()` |
| `rate_state_fault.hpp` | `RateStateFaultOperator<MeshType>`: holds state `[slip₀, θ₀, slip₁, θ₁, …]`; `PreInit()`, `Init()` (stress balance at t=0), `ComputeRHS()` (ODE right-hand side per node), `VerifyStressEquilibrium()` |

#### `friction/`
| File | Purpose |
|---|---|
| `friction_law.hpp` | `FrictionLaw` abstract base: `FrictionCoefficient(V, θ, a)`, `SolveSlipRate(tau, sn, eta, a, theta)`, `InitialState(V)` |
| `dieterich_ruina.hpp` | `DieterichRuinaFriction`: regularized DR law; Newton-Raphson bracketing in `SolveSlipRate()`; `SolveSlipRatePsi()` for psi-space; `Constants` struct (V0, f0, b, Dc) |
| `state_evolution.hpp` | `AgingLaw`, `SlipLaw`, `AgingLawPsi`, `SlipLawPsi`: `Rate(V, theta)`, `SteadyState(V)`, `RateDerivativeV()`, `RateDerivativeTheta()` |

#### `integrator/`
| File | Purpose |
|---|---|
| `dg_br2_integrator.hpp` | `BR2InteriorFaceIntegrator`, `BR2BoundaryFaceIntegrator`: `AssembleFaceMatrix()` implements consistency + symmetry + lifting penalty; penalty σ = D+1 = 3 (2D) |

#### `solver/`
| File | Purpose |
|---|---|
| `seas_operator.hpp` | `SEASQuasiDynamicOperator<MeshType>`: `TimeDependentOperator` subclass; `Mult(state, dstate_dt)` extracts slip → domain solve → traction → fault RHS; `SetInitialCondition()` runs PreInit→Init→verify cycle |
| `time_stepper.hpp` | `DormandPrinceRK45`: 5(4) embedded RK, local error control; `AdaptiveTimeStepper`: slip-rate-based dt adjustment; matches Tandem's PETSc DOPRI5 settings |
| `seas_bdrload_operator.hpp` | `BdrLoadFaultAdapter<MeshType>`: wraps `RateStateFaultOperator` with BP1 behavior — below-Wf nodes are prescribed `V=Vp` (Dirichlet in fault sense) |

#### `io/`
| File | Purpose |
|---|---|
| `probe_output.hpp` | `ProbeOutput`: writes tab-separated time-series column file at a single depth |
| `benchmark_output.hpp` | `BenchmarkOutput<MeshType>`: manages all probe outputs, adaptive output frequency (based on `V_max`), SCEC column format, interpolation to probe depths from fault DOFs |
| `parallel_benchmark_output.hpp` | `ParallelBenchmarkOutput`: gather–deduplicate–write wrapper; synchronizes output decision across all ranks via `AllRanksAgree()` |
| `paraview_output.hpp` | `ParaviewOutput`: writes displacement, traction, slip, state fields as `.vtu`/`.pvd` for ParaView |
| `checkpoint.hpp` | `Checkpoint`: `Save(t, state)` and `Load(t, state)` — binary dump of the full ODE state vector and simulation time |

#### `tests/unit/`
| File | What it tests |
|---|---|
| `test_friction_law.cpp` | `SolveSlipRate()` accuracy, Newton convergence, stress balance at steady state |
| `test_state_evolution.cpp` | Aging/slip law `Rate()`, `SteadyState()`, psi-space consistency |
| `test_antiplane.cpp` | `AntiplaneDomainOperator` construction, fault face count, traction accuracy |
| `test_fault_operator.cpp` | `RateStateFaultOperator` PreInit/Init/RHS, stress equilibrium, pre-stress value |
| `test_psi_state.cpp` | Psi-space integration correctness |
| `test_quasi_dynamic.cpp` | Full coupling loop: domain + fault `Mult()` call |
| `test_bp2_short.cpp` | Short simulation (verify earthquake triggers) |
| `test_checkpoint.cpp` | Save/load roundtrip |
| `test_io.cpp` | `BenchmarkOutput` file writing format |

#### `tests/parallel/`
| File | What it tests |
|---|---|
| `test_mpi_context.cpp` | `GlobalMax`, `GlobalSum`, `Broadcast` correctness |
| `test_parallel_utils.cpp` | `GatherVectorToRoot`, `ScatterVectorFromRoot` |
| `test_parallel_domain.cpp` | Domain operator partition consistency |
| `test_parallel_fault.cpp` | `FaultGeometry` gather and deduplication |
| `test_br2_consistency.cpp` | BR2 vs IP equivalence in limit |
| `test_serial_parallel_consistency.cpp` | Serial and parallel tractions match |
| `mms_antiplane_parallel.cpp` | Convergence rate study (Method of Manufactured Solutions) |
| `test_scaling.cpp` | Weak/strong scaling timing |

#### `tests/verification/`
| File | What it tests |
|---|---|
| `bp2_serial_smoke.cpp` | Quick 1–2 earthquake cycle smoke test |
| `bp2_verification_first_cycle.cpp` | First BP2 earthquake cycle against Tandem reference |
| `bp2_verification_full.cpp` | Full 1200-year BP2 against SCEC benchmark data |
| `bp2_benchmark_parallel.cpp` | Parallel BP2 benchmark run |
| `bp1_verification_full.cpp` | Full 3000-year BP1 all-Natural BCs driver |
| `bp1_bdrload.cpp` | Full BP1 with far-field Dirichlet BCs driver |

---

### G.3 Class Hierarchy and Key Interfaces

```
mfem::TimeDependentOperator
└── SEASQuasiDynamicOperator<MeshType>
      owns: DomainOperator<MeshType>*     (solve PDE)
      owns: RateStateFaultOperator<MeshType>*  (ODE RHS)
      owns: BenchmarkOutput / ParallelBenchmarkOutput

DomainOperator<MeshType>  [abstract]
├── AntiplaneDomainOperator<MeshType>
│     has: BilinearFormForMesh<MeshType>  stiffness_
│     has: Array<int>                    fault_dofs_
│     has: MUMPSSolver                   mumps_
└── AntiplaneBdrLoadOperator<MeshType>
      has: BilinearFormForMesh<MeshType>  stiffness_
      has: LinearFormForMesh<MeshType>    rhs_     [time-dependent]
      has: MUMPSSolver                   mumps_

RateStateFaultOperator<MeshType>
      has: FaultGeometry<MeshType>        geom_
      has: FaultNodes<MeshType>           nodes_
      has: FrictionLaw*                   friction_
      has: StateEvolution*                evolution_
      has: Vector                         state_   [2 × num_nodes]
      has: Vector                         traction_ [num_nodes]

FaultGeometry<MeshType>
      has: Vector  depths_   [num_fault_dofs]
      has: Vector  a_vals_   [num_fault_dofs]
      has: Vector  eta_vals_ [num_fault_dofs]

FrictionLaw  [abstract]
└── DieterichRuinaFriction
      has: Constants  {V0, f0, b, Dc}

StateEvolution  [abstract]
├── AgingLaw
├── SlipLaw
├── AgingLawPsi
└── SlipLawPsi

mfem::BilinearFormIntegrator
├── BR2InteriorFaceIntegrator
└── BR2BoundaryFaceIntegrator
```

**Critical interface methods to know:**

```cpp
// DomainOperator — the only methods SEAS operator calls:
virtual void Solve(real_t t, const Vector& slip, GridFuncType& u) = 0;
virtual void ComputeTraction(const GridFuncType& u, const Vector& slip,
                              Vector& traction) = 0;
virtual int  GetNumFaultDOFs() const = 0;

// RateStateFaultOperator — ODE interface:
void PreInit(const DomainOpType& dom, real_t V_init);
void Init(const DomainOpType& dom, real_t V_init);
void ComputeRHS(const Vector& traction, Vector& dstate_dt);

// SEASQuasiDynamicOperator — MFEM ODE interface:
void Mult(const Vector& state, Vector& dstate_dt) override;
//  internally: GetSlip(state) → domain.Solve() → domain.ComputeTraction()
//              → fault.ComputeRHS()
```

---

### G.4 Data Flow: One Complete Time Step

This traces exactly what happens during a single call to `SEASQuasiDynamicOperator::Mult(state, dstate_dt)`:

```
state  [size = 2 × N_fault]
  layout: [slip_0, θ_0, slip_1, θ_1, ..., slip_{N-1}, θ_{N-1}]

Step 1 — Extract slip
  slip[i] = state[2*i]       (every other element)
  Result: slip  [size = N_fault]

Step 2 — Domain solve
  AntiplaneDomainOperator::Solve(t, slip, u)
    → Assembles RHS: b = b_static - K_slip * slip
        K_slip encodes the DG jump: [[u]] = slip on fault interior faces
    → Solves K * u = b  (MUMPS factorization, cached from setup)
    → Returns grid function u over the full domain

Step 3 — Compute traction
  AntiplaneDomainOperator::ComputeTraction(u, slip, traction)
    → For each fault interior face:
        tau_h = {{mu * grad(u)}} · n_fault  (DG average of normal derivative)
    → Projects onto fault DOFs via face mass-matrix L2 projection
    → Result: traction[i] = τ at fault DOF i  [size = N_fault]

Step 4 — Compute fault ODE RHS
  RateStateFaultOperator::ComputeRHS(traction, dstate_dt)
    → For each fault DOF i:
        tau_total = traction[i] + tau_pre[i]   (add pre-stress)
        V[i]      = friction.SolveSlipRate(tau_total, sn, eta, a[i], θ[i])
        dtheta_dt = evolution.Rate(V[i], θ[i])
        dstate_dt[2*i]   = V[i]       (slip rate)
        dstate_dt[2*i+1] = dtheta_dt  (state variable rate)

Output:
  dstate_dt  [size = 2 × N_fault]
    layout: [V_0, dθ/dt_0, V_1, dθ/dt_1, ..., V_{N-1}, dθ/dt_{N-1}]
```

The ODE solver (`DormandPrinceRK45`) calls `Mult()` 6 times per step (the 6 RK stages), then checks the local truncation error to accept or reject the step and choose the next `dt`.

---

### G.5 Initialization Sequence

Before the time loop starts, a precise initialization order is required (matches Tandem's `RateAndState::init()`):

```
1. domain.SetupFaultInfo()
      → finds fault interior faces by IsFaultFace()
      → sets up fault_dofs_ array
      → assembles stiffness matrix K (once — factored by MUMPS)

2. RateStateFaultOperator::PreInit(domain, V_init)
      → calls domain.GetFaultDepths() → builds FaultGeometry
      → sets state[2*i]   = 0         (zero initial slip)
      → sets state[2*i+1] = placeholder θ

3. domain.Solve(t=0, slip=0, u)
      → first domain solve with zero slip

4. domain.ComputeTraction(u, 0, traction)
      → initial traction from zero-slip state

5. RateStateFaultOperator::Init(domain, traction, V_init)
      → for each fault DOF: solve stress balance for θ_init such that
          sigma_n * f(V_init, θ_init) + eta * V_init = traction[i] + tau_pre[i]
      → sets state[2*i+1] = θ_init

6. VerifyStressEquilibrium()
      → checks |τ - (σ_n·f(V,θ) + η·V)| < tol for all nodes
      → aborts if not satisfied
```

---

### G.6 State Vector Layout

The ODE state vector has exactly `2 × N_fault_dofs` entries. Knowing this layout is essential for reading or writing state directly.

```
index  field         description
─────  ─────         ─────────────────────────────
0      slip₀         cumulative fault slip at DOF 0 (m)
1      θ₀            state variable at DOF 0 (s) — or ψ₀ in psi-space
2      slip₁         cumulative fault slip at DOF 1
3      θ₁            state variable at DOF 1
...
2i     slipᵢ         slip at DOF i
2i+1   θᵢ            state variable at DOF i
...
2(N-1) slip_{N-1}
2N-1   θ_{N-1}
```

**Psi-space mode** (`use_psi = true`): the state variable θ is replaced by `ψ = f₀ + b·ln(V₀·θ/Dc)`, which is the argument of the `asinh` in the friction law. This improves numerical stability for very slow or very fast slip but the layout is identical.

**BP5 state vector layout** (planned, `StatePerNode = 3`):

For BP5 (full 3D elasticity with 2 tangential slip components), the state vector expands to `3 × N_fault_dofs` entries:

```
index  field         description
─────  ─────         ─────────────────────────────
0      s2₀           along-strike slip at DOF 0 (m)
1      s3₀           along-dip slip at DOF 0 (m)
2      ψ₀            state variable at DOF 0 (dimensionless, psi-space)
3      s2₁           along-strike slip at DOF 1
4      s3₁           along-dip slip at DOF 1
5      ψ₁            state variable at DOF 1
...
3i     s2ᵢ           along-strike slip at DOF i
3i+1   s3ᵢ           along-dip slip at DOF i
3i+2   ψᵢ            state variable at DOF i
...
3(N-1) s2_{N-1}
3N-2   s3_{N-1}
3N-1   ψ_{N-1}
```

The friction law operates on the slip rate norm `V = sqrt(V2² + V3²)`, and the resolved component rates are `V2 = V · τ2/||τ||`, `V3 = V · τ3/||τ||` (following Tandem's `DieterichRuinaAgeing::slip_rate`).

**Fault DOF ordering**: DOFs are ordered by the MFEM mesh's interior face traversal order. After `FaultGeometry::GatherToRootDedup()`, the root rank has a deduplicated depth-sorted array for output. The order in `state[]` is NOT depth-sorted — it follows mesh traversal order.

---

### G.7 Template Machinery: Serial and Parallel Duality

All core classes are templated on `MeshType`:

```cpp
template <typename MeshType = mfem::Mesh>
class AntiplaneDomainOperator : public DomainOperator<MeshType> { ... };
```

The `common/seas_types.hpp` header provides type traits that select the right MFEM class:

```cpp
// Serial usage:
using SerialDomain = AntiplaneDomainOperator<mfem::Mesh>;
// FiniteElementSpace, BilinearForm, GridFunction used internally

// Parallel usage (when SEAS_USE_MPI is defined):
using ParDomain = AntiplaneDomainOperator<mfem::ParMesh>;
// ParFiniteElementSpace, ParBilinearForm, ParGridFunction used internally
```

The resolution happens through:
```cpp
// seas_types.hpp
template <typename MeshType>
using FESpaceForMesh = typename MeshConditional<MeshType,
    mfem::FiniteElementSpace,    // if Mesh
    mfem::ParFiniteElementSpace  // if ParMesh
>::type;
```

**Rule for developers**: never write `FiniteElementSpace` or `ParFiniteElementSpace` directly inside a templated class. Always use `FESpaceForMesh<MeshType>`. This ensures the code compiles in both serial and parallel builds.

---

### G.8 Tracing Specific Behaviors

A quick lookup table: "I want to find where X is implemented — go to file Y, method Z."

| Question | File | Class / Method |
|---|---|---|
| Where does the stiffness matrix get assembled? | `domain/antiplane_operator.hpp` | `AntiplaneDomainOperator::AssembleStiffness()` |
| Where does the fault slip enter as a boundary condition? | `domain/antiplane_operator.hpp` | `AntiplaneDomainOperator::Solve()` — assembles the slip RHS `b_slip` |
| Where is the DG traction `{{mu * grad u}} · n` computed? | `domain/antiplane_operator.hpp` | `AntiplaneDomainOperator::ComputeTraction()` |
| Where is the fault face identified (which faces are on the fault)? | `domain/antiplane_operator.hpp` | `AntiplaneDomainOperator::IsFaultFace()` and `SetupFaultInfo()` |
| Where is the friction law Newton solver? | `friction/dieterich_ruina.hpp` | `DieterichRuinaFriction::SolveSlipRate()` |
| Where is `dθ/dt = 1 - V·θ/Dc` evaluated? | `friction/state_evolution.hpp` | `AgingLaw::Rate()` |
| Where does the ODE RHS combine traction + friction? | `fault/rate_state_fault.hpp` | `RateStateFaultOperator::ComputeRHS()` |
| Where is the initial `θ₀` computed from stress balance? | `fault/rate_state_fault.hpp` | `RateStateFaultOperator::Init()` |
| Where is slip extracted from the ODE state vector? | `solver/seas_operator.hpp` | `SEASQuasiDynamicOperator::Mult()` — `GetSlip(state, slip)` call |
| Where is `dt` adjusted based on `V_max`? | `solver/time_stepper.hpp` | `AdaptiveTimeStepper::ComputeNewDt()` |
| Where is the RK45 step accepted/rejected? | `solver/time_stepper.hpp` | `DormandPrinceRK45::Step()` |
| Where are MPI fault DOFs deduplicated? | `fault/fault_geometry.hpp` | `FaultGeometry::GatherToRootDedup()` |
| Where is the probe output written to disk? | `io/benchmark_output.hpp` | `BenchmarkOutput::Write()` |
| Where is the parallel gather for output? | `io/parallel_benchmark_output.hpp` | `ParallelBenchmarkOutput::Write()` |
| Where is the `tau_pre` pre-stress set? | `fault/rate_state_fault.hpp` | `RateStateFaultOperator::InitPreStress()` |
| Where is the below-fault `V=Vp` prescription? | `solver/seas_bdrload_operator.hpp` | `BdrLoadFaultAdapter::ComputeRHS()` |
| Where is MUMPS configured as SPD or indefinite? | `domain/antiplane_operator.hpp` | `AntiplaneDomainOperator::AssembleStiffness()` — `mumps_->SetMatrixSymType(...)` |
| Where is the far-field Dirichlet loading applied? | `domain/antiplane_bdrload_operator.hpp` | `AntiplaneBdrLoadOperator::Solve()` — `DGDirichletLFIntegrator` added to RHS |
| Where are boundary attributes defined? | `domain/bp2_mesh.hpp` | `BP2BoundaryAttributes` struct |
| Where is the BR2 lifting penalty assembled? | `integrator/dg_br2_integrator.hpp` | `BR2InteriorFaceIntegrator::AssembleFaceMatrix()` |

---

### G.9 How to Add New Components

#### Add a new friction law

1. Create `friction/my_friction.hpp` with class `MyFriction : public FrictionLaw`.
2. Implement `FrictionCoefficient(V, theta, a)`, `SolveSlipRate(tau, sn, eta, a, theta)`, `InitialState(V)`.
3. In the driver (e.g., `pseas.cpp`), replace `DieterichRuinaFriction` with `MyFriction`.
4. Add `test_my_friction.cpp` following the pattern in `test_friction_law.cpp`.
5. No other files need to change — `FrictionLaw*` is stored by pointer.

#### Add a new domain operator (e.g., plane strain)

1. Create `domain/plane_strain_operator.hpp` with class `PlaneStrainOperator<MeshType> : public DomainOperator<MeshType>`.
2. Implement all pure-virtual methods: `Solve()`, `ComputeTraction()`, `GetNumFaultDOFs()`, `GetFaultDepths()`, `GetFaultDOFs()`.
3. Use `FESpaceForMesh<MeshType>` for FE space type — never write `FiniteElementSpace` directly.
4. In the driver, replace `AntiplaneDomainOperator<ParMesh>` with `PlaneStrainOperator<ParMesh>`.
5. `SEASQuasiDynamicOperator` does not need to change — it calls only the base interface.

#### Add a new output field

1. Locate `io/benchmark_output.hpp`, `BenchmarkOutput::Write()`.
2. The write method receives the fault DOF arrays (`slip`, `slip_rate`, `traction`, `state`). Add the new column to the output format.
3. Update the file header string to include the new column name.
4. If a new quantity requires a domain-level query (e.g., off-fault displacement), add a virtual method to `DomainOperator` with a default no-op, then implement it in the specific operator.

#### Add a new benchmark (e.g., BP3)

1. Create `bp3/bp3.geo` — Gmsh geometry for BP3 domain (60° dipping fault).
2. Create `config/bp3_params.hpp` — `BP3Params` struct with BP3-specific physical parameters.
3. Add `Physical Curve` tags for fault and BCs in the `.geo` file.
4. Create `tests/verification/bp3_verification.cpp` — driver following `bp1_bdrload.cpp` structure.
5. Key difference from BP1/BP2: `IsFaultFace()` must use the physical group tag, not a coordinate check, since the fault is dipping (see Part C for fault detection plan).

---

### G.10 Coordinate Conventions and Sign Conventions

**Spatial coordinates** (2D problems):
```
x  ∈ [-Lx, +Lx]   horizontal, fault-normal direction
z  ∈ [-Lz, 0]     vertical, z = 0 at free surface, z < 0 at depth
```
Fault is at `x = 0`. Left side: `x < 0`. Right side: `x > 0`.

**Depth convention**: depth is reported as `|z|` (positive), not `z` (negative). In `FaultGeometry`, `depths_[i] = std::abs(z_center_of_face)`.

**Traction sign**: traction `tau` is positive right-lateral (positive slip on positive side relative to negative side). Pre-stress `tau_pre` is **negative** (it opposes the applied shear), following Tandem's convention: the physical stress is `tau + tau_pre`.

**Slip sign**: `slip > 0` means the positive-x side moved in the positive-z direction relative to the negative-x side. This is the standard right-lateral convention for a vertical strike-slip fault.

**State variable**: `theta` has units of **seconds** (time). Steady state: `θ_ss = Dc / V`. Psi-space: `ψ = f₀ + b · ln(V₀·θ/Dc)`, dimensionless.

**Pre-stress `tau_pre`**: computed at `t=0` from `tau_pre = -(sigma_n * f(V_init, theta_ss) + eta * V_init)` so that at `t=0`, `traction + tau_pre = 0` when the domain solve returns `traction = tau_pre_magnitude`. In other words, `tau_pre` makes the initial net driving stress consistent with `V_init`.

---

### G.11 Key Numbers and Sanity Checks

When debugging, these numbers indicate the simulation is behaving correctly:

| Quantity | Expected value | Where to check |
|---|---|---|
| Initial `V_max` | ~`V_init = 1e-9 m/s` (BP1/BP2) | First line of output |
| Pre-stress `tau_pre` | ~26.546 MPa (BP2, shallow zone) | `test_fault_operator.cpp` TEST for this value |
| Time to first earthquake (BP2) | ~30–40 years | `bp2_verification_first_cycle.cpp` |
| Peak coseismic `V_max` | ~1–10 m/s | Output log during earthquake |
| Recurrence interval (BP2) | ~27–30 years | Full benchmark verification |
| Fault face count (BP2, 400m mesh) | ~200 faces | `test_antiplane.cpp` fault face count test |
| `theta_init` at seismogenic depth | `Dc / V_init = 4000 s` (BP2) | Computed in `Init()` |

If `tau_pre` is far from 26.546 MPa, the parameters or pre-stress computation are wrong. If the fault face count is zero, `IsFaultFace()` or mesh coordinates are wrong. If `V_max` explodes immediately, the initialization stress balance failed.

---

### G.12 Important Conventions Inherited from Tandem

The MFEM SEAS miniapp deliberately mirrors Tandem's interface and algorithmic choices to make porting verification easier. The key inherited conventions are:

| Convention | Tandem location | MFEM location |
|---|---|---|
| Psi-space state variable | `RateAndStateBase.h: psi_` | `rate_state_fault.hpp: use_psi_` |
| DG-BR2 lifting operator | `Elasticity.cpp: assemble_skeleton` | `dg_br2_integrator.hpp: AssembleFaceMatrix()` |
| Dirichlet-as-slip BC | `Poisson.cpp: rhs_integration` | `antiplane_operator.hpp: AssembleSlipRHS()` |
| RK45 step acceptance | `PETSc DOPRI5` | `time_stepper.hpp: DormandPrinceRK45::Step()` |
| Below-fault Vp prescription | `bp2.lua: boundary()` | `seas_bdrload_operator.hpp: BdrLoadFaultAdapter` |
| Fault physical group tags | `BC.h: enum BC::Fault = 3` | Planned in Part C (currently coordinate-based) |

When a result does not match Tandem, check this table first — the conventions may have diverged during porting.

---

## Part H: Prioritized Implementation Recommendations

Based on review of the plan against the actual codebase, the following prioritized action items are recommended. Items are ordered by impact on BP5 readiness.

| Priority | Action | Effort | Impact | Dependency |
|----------|--------|--------|--------|------------|
| 1 | **Implement Part C** (fault tag detection via Gmsh Physical Groups) | 2–3 days | Unblocks BP3/BP5 mesh design, removes fragile coordinate check | None |
| 2 | **Implement Phase 1b** (interface cleanup: `NumSlipComponents()`, generalize `SEASQuasiDynamicOperator` to use base `DomainOperator<MeshType>*`) | 1–2 days | Unblocks BP5 operator design, enables polymorphic domain dispatch | None |
| 3 | **Write `test_fault_detection.cpp`** (E.3) | 1 day | Validates Part C; use embedded `.msh` string for programmatic test mesh | Part C |
| 4 | **Write `test_domain_operator_interface.cpp`** (E.4) | 0.5 days | Validates Phase 1b interface changes | Phase 1b |
| 5 | **Extract shared BP1 driver code** into `bp1_common.hpp` (Phase 1c) | 1 day | Eliminates 95% duplication without TOML overhead | None |
| 6 | **Defer Part B** (TOML input system) until BP5 is working | — | Avoids premature abstraction; full config shape unknown until BP5 works | — |
| 7 | **Defer E.1, E.2, E.5** (TOML-dependent tests) with Phase 1d | — | Keeps test suite buildable; avoids broken targets | Phase 1d |
| 8 | **Defer E.6** (3D fault detection test) until Phase 3 starts | — | Test code that doesn't exist yet creates broken targets | Phase 3 |

### Key risk mitigations

1. **Backward compatibility**: `fault_tag = -1` default preserves existing coordinate-based behavior. No `.geo` changes required for existing jobs.
2. **Regression guard**: test E.3.5 (tag-based == coordinate-based on bp1_100m mesh) catches any accidental change to fault face counting.
3. **Incremental validation**: each phase has its own test file that must pass before the next phase starts.
4. **No premature abstraction**: the TOML config system is designed after BP5 reveals the full parameter space, not before.
