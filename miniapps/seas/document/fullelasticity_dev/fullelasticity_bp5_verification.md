# BP5 Benchmark Verification Plan

**Date**: 2026-03-04
**Target file**: `document/fullelasticity_dev/fullelasticity_bp5_verification.md`

## Context

All core physics components (Phases 1-4) are complete and verified: 3D DG elasticity,
vector fault operator, friction laws, SEAS coupling operator, adaptive RK45 time stepper,
and 7 passing integration tests. The remaining gap is the BP5 benchmark driver and output
infrastructure to run actual simulations and produce SCEC-format results.

This plan follows the **existing BP1 infrastructure patterns** exactly:
- Driver location: `tests/verification/bp5_verification_full.cpp` (like `bp1_verification_full.cpp`)
- Output classes: `BP5BenchmarkOutput` (like `BenchmarkOutput`), same ProbeOutput, same parallel wrapper pattern
- Directory structure: `bp5/mesh/`, `bp5/benchmark_data/` (like `bp1/mesh/`, `bp1/benchmark_data/`)
- Makefile: same pattern as `seas_bp1_full` target
- CLI args: same `--mesh`, `--mesh-scale`, `--output-dir`, `--tfinal`, `--checkpoint-interval`, `--restart`, `--comparison-only`

---

## 1. Readiness Assessment

### 1.1 Components Ready (Phases 1-4 Complete)

| Component | File | Status |
|-----------|------|--------|
| BP5 parameters | `config/bp5_params.hpp` | Complete (354 lines) |
| 3D DG elasticity | `domain/elasticity_operator.hpp` | Complete (~550 lines, BR2/IP) |
| Fault basis transform | `fault/fault_basis.hpp` | Complete (~270 lines) |
| Fault geometry (2D) | `fault/fault_geometry.hpp` | Complete (635 lines, BP5 constructor) |
| Vector fault operator | `fault/rate_state_fault.hpp` | Complete (804 lines, SlipComponents=2) |
| Dieterich-Ruina friction | `friction/dieterich_ruina.hpp` | Complete |
| Aging law (psi-space) | `friction/state_evolution.hpp` | Complete |
| SEAS operator | `solver/seas_operator.hpp` | Complete (239 lines, BP5SEASOp alias) |
| DormandPrince RK45 | `solver/time_stepper.hpp` | Complete (577 lines, SetStatePerNode(3)) |
| Probe output (generic) | `io/probe_output.hpp` | Complete (ProbeOutput + ProbeInterpolator) |
| Checkpoint I/O | `io/checkpoint.hpp` | Complete (size-agnostic) |
| MPI context | `common/mpi_context.hpp` | Complete |
| Integration tests | `tests/unit/test_bp5_integration.cpp` | Complete (7 tests passing) |

### 1.2 What's Missing

| Gap | Priority | Notes |
|-----|----------|-------|
| BP5 output infrastructure | **HIGH** | 8-column vector format, 2D probe matching |
| BP5 driver program | **HIGH** | `tests/verification/bp5_verification_full.cpp` |
| BP5 parallel output | **MEDIUM** | `io/bp5_parallel_output.hpp` (parallel gather wrapper) |
| BP5 mesh files | **MEDIUM** | Gmsh .geo + .msh for fault refinement |

---

## 2. BP5 Mesh Configurations

### 2.1 Coordinate Convention (MFEM ↔ SCEC ↔ Tandem)

| Direction | MFEM mesh | SCEC BP5 | Tandem |
|-----------|-----------|----------|--------|
| Fault-normal | x | x₁ | Y |
| Along-strike | y | x₂ | X |
| Depth (positive down) | z | x₃ | Z (negative) |

- Fault plane: x = 0 (internal interface)
- RSF zone: |y| ≤ lf/2 = 50 km, 0 ≤ z ≤ Wf = 40 km
- Free surface: z = 0 (natural BC = zero traction)
- Loading walls: y = ±Ly (Dirichlet: u_y = ±Vp·t/2)

### 2.2 Boundary Attribute Convention

| Attr | Face | BC | Description |
|------|------|----|-------------|
| 1 | x = -Lx | Natural | Far-field (zero traction) |
| 2 | x = +Lx | Natural | Far-field (zero traction) |
| 3 | y = +Ly | Dirichlet | Plate loading: u_y = +Vp·t/2 |
| 4 | y = -Ly | Dirichlet | Plate loading: u_y = -Vp·t/2 |
| 5 | z = 0 | Natural | Free surface (zero traction) |
| 6 | z = Lz | Natural | Deep boundary (zero traction) |

### 2.3 Inline Cartesian Mesh (for smoke tests)

All meshes use `Mesh::MakeCartesian3D(2*nx, 2*ny, nz, HEXAHEDRON, 2*Lx, 2*Ly, Lz)`
then shift to center: x∈[-Lx,Lx], y∈[-Ly,Ly], z∈[0,Lz].

| Config | nx | ny | nz | Elements | Cell size | Purpose |
|--------|----|----|-----|----------|-----------|---------|
| A Smoke | 2 | 2 | 1 | 16 | ~50 km | Code runs, no NaN |
| B Coarse | 4 | 4 | 2 | 128 | ~50 km | Init stability, output files |
| C Medium | 20 | 10 | 10 | 8000 | ~10 km | First earthquake |
| D Moderate | 40 | 20 | 20 | 64000 | ~5 km | Visible earthquake cycle |
| E Benchmark | 200 | 100 | 50 | 4M | ~1 km | SCEC benchmark quality |

Domain: [-200,200] × [-100,100] × [0,100] km (matching Tandem's 400×200×100 km).

### 2.4 Gmsh Mesh (bp5/mesh/bp5.geo)

Adapted from Tandem's `examples/tandem/3d/bp5.geo` with MFEM coordinate convention:
- Two half-volumes split at fault plane x=0
- OpenCASCADE BooleanFragments for fault surface embedding
- Refinement zones matching Tandem: fault res=1km, nucleation res=0.5km, far-field res=40km
- Physical surfaces tagged with MFEM boundary attributes (3,4=Dirichlet ±y)

### 2.5 BP5 Mesh Files (following BP1 directory structure: `bp1/mesh/`)

```
bp5/mesh/
├── bp5.geo              # Gmsh geometry (adapted from Tandem)
├── bp5_coarse.msh       # ~128 elements (for testing)
└── bp5_fine.msh         # ~64K+ elements (for production)
```

---

## 3. BP5 Output Infrastructure

### 3.1 Key Differences from BP1/BP2 Output

| Aspect | BP1/BP2 (`BenchmarkOutput`) | BP5 (`BP5BenchmarkOutput`) |
|--------|-----|-----|
| Probes | 1D depth (z) | 2D (x2, x3) = (strike, depth) |
| Interpolator | `ProbeInterpolator` (1D, z) | `Probe2DInterpolator` (2D, x2+x3) |
| Columns | 5 (scalar) | 8 (vector, 2 components) |
| Slip | 1 component | 2: slip_strike, slip_dip |
| Slip rate | 1 component | 2: V_strike, V_dip |
| Stress | 1: shear_stress | 2: tau_strike, tau_dip |
| State | log10(theta) | log10(theta) |
| Component swap | N/A | Internal [dip,strike] → Output [strike,dip] |
| Naming | `{prefix}_z{depth}km.txt` | `{prefix}_fltst_strk{x2}dp{x3}.txt` |

### 3.2 New File: `io/bp5_benchmark_output.hpp` (~400 lines)

**Pattern**: Follows `io/benchmark_output.hpp` exactly — same class structure, same
ProbeOutput usage, same adaptive output, same ForceWrite/Flush/Close interface.

#### Class: `Probe2DInterpolator`

2D version of `ProbeInterpolator` from `io/probe_output.hpp`:
- Takes fault DOF coordinates (x2, x3) from `ElasticityDomainOperator::GetFaultCoords2D()`
- Takes station list `{name, x2, x3}`
- Finds nearest DOF for each station (simple nearest-neighbor, no interpolation needed for DG)
- Reports matching distance for diagnostics

```cpp
Probe2DInterpolator(const Vector &fault_x2, const Vector &fault_x3,
                    const std::vector<Station> &stations);
int GetNearestDOF(int station_idx) const;
real_t GetMatchDistance(int station_idx) const;
int NumStations() const;
```

#### Class: `BP5BenchmarkOutput`

Same structure as `BenchmarkOutput<MeshType>`:

```cpp
template <typename MeshType = Mesh>
class BP5BenchmarkOutput
{
public:
   struct Station { std::string name; real_t x2; real_t x3; };

   BP5BenchmarkOutput(const std::string &prefix,
                      const BP5Params &params,
                      const std::vector<Station> &stations,
                      const Vector &fault_x2, const Vector &fault_x3);

   bool Write(real_t time, const Vector &state,
              const RateStateFaultOperator<MeshType, 2> &fault,
              const Vector &traction, real_t V_max);

   void ForceWrite(real_t time, const Vector &state,
                   const RateStateFaultOperator<MeshType, 2> &fault,
                   const Vector &traction, real_t V_max);

   void WriteFromGlobalData(real_t time,
                            const Vector &global_slip,
                            const Vector &global_theta,
                            const Vector &global_V,
                            const Vector &global_traction,
                            const BP5Params &params);

   void Flush();
   void Close();

   static std::vector<Station> DefaultStations();  // 10 SCEC on-fault stations
   static real_t OutputInterval(real_t V_max);      // Same as BenchmarkOutput
};
```

Output per station (8 columns, SCEC BP5 format):
```
# BP5-QD time series at fltst_strk+00dp+10 (x2=0.0km, x3=10.0km)
# Columns: t(s), slip_strike(m), slip_dip(m), log10(V_strike)(m/s),
#          log10(V_dip)(m/s), tau_strike(MPa), tau_dip(MPa), log10(state)(s)
```

Component ordering: Internal MFEM [dip=0, strike=1] → Output [strike, dip] **swap on output**.

Stress computation: `tau_i = tau0_i + traction_i` (total shear, divide by 1e6 for MPa).

Adaptive output intervals: Same `OutputInterval()` as BP1/BP2 (V-based phase detection).

#### Default Stations (10 SCEC BP5 on-fault stations)

```cpp
static std::vector<Station> DefaultStations()
{
   // SCEC BP5 stations: fltst_strk{x2}dp{x3}
   // x2 = along-strike (km), x3 = depth (km)
   return {
      {"fltst_strk-24dp+10", -24e3, 10e3},
      {"fltst_strk-16dp+10", -16e3, 10e3},
      {"fltst_strk+00dp+00",   0.0,  0.0},
      {"fltst_strk+00dp+05",   0.0,  5e3},
      {"fltst_strk+00dp+10",   0.0, 10e3},
      {"fltst_strk+00dp+15",   0.0, 15e3},
      {"fltst_strk+00dp+20",   0.0, 20e3},
      {"fltst_strk+00dp+30",   0.0, 30e3},
      {"fltst_strk+16dp+10",  16e3, 10e3},
      {"fltst_strk+24dp+10",  24e3, 10e3},
   };
}
```

#### Global Output

Following BP1 pattern (currently BP1 doesn't have a separate global output file,
but BP5 needs one per SCEC spec):

```
Columns: t(s), log10(Vmax)(m/s), moment_rate(N-m/s)
```

This can be a single ProbeOutput instance in the driver.

### 3.3 New File: `io/bp5_parallel_output.hpp` (~150 lines)

**Pattern**: Follows `io/parallel_benchmark_output.hpp` exactly.

```cpp
class ParallelBP5BenchmarkOutput
{
   FaultGeometry<ParMesh> &fault_geom_;
   MPIContext &mpi_ctx_;
   std::unique_ptr<BP5BenchmarkOutput<Mesh>> bench_out_;  // root only
   real_t last_write_time_;

public:
   ParallelBP5BenchmarkOutput(const std::string &prefix,
                               const BP5Params &params,
                               const std::vector<BP5BenchmarkOutput<Mesh>::Station> &stations,
                               FaultGeometry<ParMesh> &fault_geom,
                               MPIContext &mpi_ctx,
                               const Vector &dedup_x2, const Vector &dedup_x3);

   bool Write(real_t time, const Vector &state,
              const RateStateFaultOperator<ParMesh, 2> &fault,
              const Vector &traction, real_t global_V_max);
   void ForceWrite(...);
   void Flush();
   void Close();
};
```

Parallel gather: 7 scalar gathers per write (slip_dip, slip_strike, theta,
V_dip, V_strike, trac_dip, trac_strike) via `GatherFieldsToRootDedup`.

---

## 4. BP5 Driver Program

### 4.1 New File: `tests/verification/bp5_verification_full.cpp` (~900 lines)

**Pattern**: Follows `tests/verification/bp1_verification_full.cpp` exactly —
same CLI args, same component stack flow, same time loop, same checkpoint/restart,
same comparison infrastructure.

#### Structure (matching BP1 section by section)

```
1. Reference data loading utilities (adapted for 8-column BP5 format)
   - LoadBP5TimeSeriesFile(): 8 columns instead of 5
   - InterpolateOnto(): reused from BP1 (same function)
   - RelativeL2Error(): reused from BP1 (same function)
   - MakeFilename(): BP5 station naming (fltst_strk{x2}dp{x3}.txt)
   - RunComparison(): adapted for 8-column data

2. CLI argument parsing (same as BP1, plus BP5-specific)
   --mesh FILE                 Gmsh .msh file
   --mesh-scale S              km → m (default: 1000)
   --output-dir DIR            Output directory
   --output-prefix PFX         Output prefix (default: "bp5_full")
   --tfinal T                  Override t_final (default: 1800 years)
   --checkpoint-interval N     Steps between checkpoints (default: 5000)
   --restart PREFIX             Restart from checkpoint
   --ref-dir DIR               Reference data directory (default: bp5/benchmark_data)
   --comparison-only           Skip simulation, compare only
   --inline-mesh               Use inline Cartesian mesh (for smoke tests)
   --nx/--ny/--nz N            Inline mesh element counts
   --Lx/--Ly/--Lz L            Inline mesh domain half-sizes [m]
   --write-every-step          Write output at every accepted step

3. Mesh creation (matching BP1 pattern)
   If --mesh: Load Gmsh, apply mesh_scale
   If --inline-mesh: Create3DMesh(nx,ny,nz,Lx,Ly,Lz) with boundary attrs
   Default: require --mesh (like BP1)

4. BP5 component stack (matching BP1 stack, elasticity instead of antiplane)
   BP5Params params;
   params.Validate();
   ElasticityDomainOperator<ParMesh> domain(pmesh, order, params.lambda(),
      params.mu(), params.Vp, params.Wf, params.lf, DGMethod::BR2);
   FaultGeometry<ParMesh> fault_geom(domain, params, &mpi);
   DieterichRuinaFriction friction(fc);
   AgingLawPsi aging(params.b, params.V0, params.f0);
   RateStateFaultOperator<ParMesh, 2> fault_op(&fault_geom, &friction, &aging,
      params, &mpi, true);
   PBP5SEASOp seas_op(&domain, &fault_op, &mpi);

5. Initialization (same as BP1)
   Vector state(fault_op.StateSize());
   seas_op.SetInitialCondition(state);

6. I/O setup (matching BP1's parallel output pattern)
   // Gather global fault coords for probe matching
   Vector local_x2, local_x3;
   domain.GetFaultCoords2D(local_x2, local_x3);
   Vector dedup_x2, dedup_x3;
   // Gather and dedup x2, x3 coordinates
   ...
   auto stations = BP5BenchmarkOutput<Mesh>::DefaultStations();
   ParallelBP5BenchmarkOutput bench_out(prefix, params, stations,
      fault_geom, mpi, dedup_x2, dedup_x3);
   // Global output (root only)
   ProbeOutput global_out(prefix + "_global.txt",
      {"time(s)", "log10(Vmax)(m/s)"}, "BP5-QD global output");

7. Time integration (same as BP1, with BP5-specific settings)
   DormandPrinceRK45 rk45;
   rk45.SetMPIContext(&mpi);
   rk45.SetAbsTol(1e-7);
   rk45.SetRelTol(1e-50);
   rk45.SetDtMin(1e-6);
   rk45.SetDtMax(0.5 * seconds_per_year);
   rk45.SetStatePerNode(3);     // BP5: 3 values per node (slip_dip, slip_strike, psi)
   rk45.Init(seas_op);

8. Main loop (same structure as BP1)
   while (t < t_final && step < max_steps)
   {
      bool accepted = rk45.Step(seas_op, state, t, dt);
      if (!accepted) continue;
      step++;
      V_max = seas_op.GetMaxSlipRate();
      // NaN check (same as BP1)
      // Earthquake detection (same as BP1)
      // Adaptive output (same as BP1)
      bench_out.Write(t, state, fault_op, seas_op.GetTraction(), V_max);
      global_out.WriteStep({t, std::log10(V_max)});
      // Checkpoint (same as BP1)
      // Console output (same as BP1)
   }

9. Final checkpoint + close (same as BP1)

10. Post-simulation comparison (same pattern as BP1)
    if (mpi.IsRoot())
    {
       RunComparison(output_dir, output_prefix, ref_dir, stations, t_final);
    }
```

### 4.2 Inline Mesh Creation Helper

```cpp
/// Create a 3D hex mesh for BP5 with proper boundary attributes.
/// Domain: [-Lx,Lx] × [-Ly,Ly] × [0,Lz]
/// Fault plane at x=0 (internal interface between adjacent elements).
std::unique_ptr<Mesh> CreateBP5InlineMesh(
   int nx, int ny, int nz,
   real_t Lx, real_t Ly, real_t Lz);
```

This follows `bp2_mesh.hpp` pattern but for 3D. Shifts the Cartesian mesh so
x∈[-Lx,Lx], y∈[-Ly,Ly], z∈[0,Lz] and assigns boundary attributes 1-6.

### 4.3 Makefile Updates

Following the exact BP1 pattern in the existing Makefile:

```makefile
# BP5 full multi-cycle simulation
BP5_FULL_SRC = tests/verification/bp5_verification_full.cpp
BP5_FULL_OBJ = $(BP5_FULL_SRC:.cpp=.o)

# Add to PAR_LONG:
PAR_LONG = ... seas_bp5_full

# Build rule (same pattern as seas_bp1_full):
seas_bp5_full: $(BP5_FULL_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(BP5_FULL_OBJ) $(MFEM_LIBS)

$(BP5_FULL_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) \
                       $(PARALLEL_IO_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@
```

Add `io/bp5_benchmark_output.hpp` and `io/bp5_parallel_output.hpp` to headers.

---

## 5. Unit Tests for New Functions

### 5.1 New File: `tests/unit/test_bp5_output.cpp` (~300 lines)

Tests for new classes in `io/bp5_benchmark_output.hpp`.

#### Test: Probe2DInterpolator_NearestDOF
- Create synthetic fault DOF grid (10×8 = 80 DOFs) with known (x2, x3) coords
- Query 10 default SCEC stations
- Verify each station maps to the nearest DOF (check index and distance)
- Verify edge cases: station at exact DOF location → distance ≈ 0

#### Test: Probe2DInterpolator_NoFaultDOFs
- Empty fault DOF arrays → verify graceful handling (or MFEM_VERIFY)

#### Test: BP5BenchmarkOutput_FileCreation
- Construct BP5BenchmarkOutput with 3 test stations and synthetic fault coords
- Verify 3 output files are created with correct names (station-based naming)
- Verify each file has correct 8-column header

#### Test: BP5BenchmarkOutput_ComponentSwap
- Write one row of known data: slip=[dip=1.0, strike=2.0]
- Read output file
- Verify columns are: slip_strike=2.0, slip_dip=1.0 (swapped)

#### Test: BP5BenchmarkOutput_StressComputation
- Write data with known traction and tau0
- Verify output stress = (tau0 + traction) / 1e6 in MPa

#### Test: BP5BenchmarkOutput_AdaptiveOutput
- Call `OutputInterval()` with different V_max values
- Verify: V > 1e-3 → 0.001s, 1e-6 < V < 1e-3 → 0.1s, V < 1e-6 → 0.01 yr

#### Test: BP5BenchmarkOutput_DefaultStations
- Call `DefaultStations()`
- Verify 10 stations returned with correct (x2, x3) coordinates per SCEC spec

### 5.2 New File: `tests/unit/test_bp5_mesh.cpp` (~200 lines)

Tests for `CreateBP5InlineMesh` helper function.

#### Test: CreateBP5InlineMesh_ElementCount
- Create mesh with nx=2, ny=2, nz=1
- Verify element count = 4*4*1 = 16
- Verify vertex count matches hexahedral mesh

#### Test: CreateBP5InlineMesh_BoundaryAttributes
- Create mesh, iterate boundary elements
- Verify 6 distinct boundary attributes
- Verify attr 1,2 on x-faces, attr 3,4 on y-faces, attr 5,6 on z-faces

#### Test: CreateBP5InlineMesh_DomainExtent
- Create mesh with Lx=200e3, Ly=100e3, Lz=100e3
- Verify vertex coordinates span [-200e3, 200e3] × [-100e3, 100e3] × [0, 100e3]

#### Test: CreateBP5InlineMesh_FaultInterface
- Create mesh, check that internal faces at x≈0 exist
- Verify these can be identified as fault faces by `IsFaultFace3D()`

### 5.3 Makefile Updates for Unit Tests

```makefile
# BP5 Phase 5: Output tests
TEST_BP5_OUTPUT_SRC = tests/unit/test_bp5_output.cpp
TEST_BP5_OUTPUT_OBJ = $(TEST_BP5_OUTPUT_SRC:.cpp=.o)

# BP5 Phase 5: Mesh tests
TEST_BP5_MESH_SRC = tests/unit/test_bp5_mesh.cpp
TEST_BP5_MESH_OBJ = $(TEST_BP5_MESH_SRC:.cpp=.o)

seas_test_bp5_output: $(TEST_BP5_OUTPUT_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BP5_OUTPUT_OBJ) $(MFEM_LIBS)

seas_test_bp5_mesh: $(TEST_BP5_MESH_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BP5_MESH_OBJ) $(MFEM_LIBS)

# Add to SEQ_MINIAPPS and test target
```

---

## 6. Implementation Order

### Phase A: BP5 Output Infrastructure + Unit Tests (HIGH)

1. Create `io/bp5_benchmark_output.hpp`
   - `Probe2DInterpolator` (nearest-DOF matching in 2D)
   - `BP5BenchmarkOutput` (8-column writer, same interface as `BenchmarkOutput`)
   - Default 10 SCEC stations

2. Create `tests/unit/test_bp5_output.cpp`
   - 7 unit tests for Probe2DInterpolator and BP5BenchmarkOutput
   - Verify file creation, component swap, stress computation, adaptive output

3. Create `io/bp5_parallel_output.hpp`
   - `ParallelBP5BenchmarkOutput` (same pattern as `ParallelBenchmarkOutput`)

### Phase B: BP5 Mesh + Unit Tests (HIGH)

4. Implement `CreateBP5InlineMesh` helper (in driver or separate header)

5. Create `tests/unit/test_bp5_mesh.cpp`
   - 4 unit tests for mesh creation: element count, boundary attrs, domain extent, fault interface

### Phase C: BP5 Driver (HIGH)

6. Create `tests/verification/bp5_verification_full.cpp`
   - Full driver following BP1 pattern
   - Inline mesh support for smoke tests
   - Comparison/verification mode

7. Update `Makefile` (add `seas_bp5_full`, `seas_test_bp5_output`, `seas_test_bp5_mesh`)

### Phase D: Mesh + Smoke Test (MEDIUM)

8. Create `bp5/mesh/bp5.geo` (adapted from Tandem)
9. Generate meshes: `gmsh -3 bp5.geo -o bp5_coarse.msh`
10. Smoke test: `mpirun -np 2 ./seas_bp5_full --inline-mesh --nx 2 --ny 2 --nz 1`

### Phase E: Benchmark Runs (after code verified)

11. Run with Config C/D meshes for first earthquake detection
12. Compare with Tandem results if available

### Phase F: Deferred

13. Off-fault station output (9 bulk stations)
14. Moment rate computation
15. Earthquake catalog

---

## 7. Verification Strategy

### 7.1 Unit Tests (run immediately after implementation)
```bash
make seas_test_bp5_output && ./seas_test_bp5_output
make seas_test_bp5_mesh && ./seas_test_bp5_mesh
```
All 11 unit tests must pass before proceeding to integration testing.

### 7.2 Smoke Test (Config A, minutes)
```bash
mpirun -np 2 ./seas_bp5_full --inline-mesh --nx 2 --ny 2 --nz 1 \
    --Lx 100e3 --Ly 100e3 --Lz 50e3 --tfinal 3.15e10 --output-dir bp5_smoke
```
Verify: no NaN, 10 station files with 8 columns, global output.

### 7.3 Coarse Test (Config B, ~30 min)
```bash
mpirun -np 4 ./seas_bp5_full --inline-mesh --nx 4 --ny 4 --nz 2 \
    --Lx 200e3 --Ly 100e3 --Lz 100e3 --tfinal 6.3e10 --output-dir bp5_coarse
```
Verify: initialization (stress equilibrium < 1e-6), slip accumulates, V_max bounded.

### 7.4 Initialization Cross-Checks

| Check | Expected | Tolerance |
|-------|----------|-----------|
| V_max after init | ~V_nuc = 0.01 m/s | 50% relative |
| Stress equilibrium | < 1e-6 | Absolute |
| StateSize | 3 × num_fault_DOFs | Exact |
| SlipSize | 2 × num_fault_DOFs | Exact |
| Initial slip | 0.0 | Machine eps |
| Initial psi | f0 + b·ln(V0/V_init) ≈ 1.29 | 1% |
| Station files | 10 + 1 (global) | Exact |
| Columns per station | 8 | Exact |

### 7.5 Physics Sanity Checks

1. **Interseismic**: V_max ≈ V_init (1e-9 m/s), slowly increasing
2. **Nucleation**: V_max ramps in nucleation zone → 1e-3 m/s
3. **Coseismic**: V_max → ~1 m/s, rupture propagates
4. **Post-seismic**: V_max drops below 1e-3, then 1e-6
5. **Cycle repeat**: ~100-300 years later

---

## 8. Key Technical Notes

### 8.1 Reuse from BP1 Infrastructure

| Reused (unchanged) | Adapted for BP5 |
|---|---|
| `io/probe_output.hpp` (ProbeOutput) | New `Probe2DInterpolator` (2D nearest-DOF) |
| `io/checkpoint.hpp` (size-agnostic) | No changes needed |
| `solver/time_stepper.hpp` (SetStatePerNode(3)) | Already supports BP5 |
| `solver/seas_operator.hpp` (BP5SEASOp alias) | Already defined |
| `common/mpi_context.hpp` | Already supports BP5 |
| CLI arg parsing pattern | Same args + --inline-mesh |
| Earthquake detection logic | Same V thresholds |
| Checkpoint/restart logic | Same (size-agnostic) |

### 8.2 V_nuc Discrepancy
- SCEC spec: V_nuc = 0.03 m/s
- Tandem bp5.lua: V_nuc = 0.01 m/s
- `bp5_params.hpp` follows Tandem: V_nuc = 0.01 m/s

### 8.3 Component Ordering
- Internal MFEM: index 0 = dip, index 1 = strike
- SCEC output: strike first, dip second → **swap on output**

### 8.4 State Variable Convention
- BP5 uses psi-space: ψ = f₀ + b·ln(V₀·θ/Dc)
- Output: log10(θ) where θ = (Dc/V₀)·exp((ψ - f₀)/b)
- `fault.GetTheta()` returns physical theta (already converted)

### 8.5 Checkpoint Compatibility
Existing `io/checkpoint.hpp` is size-agnostic. BP5 state = 3N, displacement = 3×scalar_DOFs.
No changes needed.
