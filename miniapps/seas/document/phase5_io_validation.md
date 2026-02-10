# Phase 5: I/O and Serial Validation

## Overview

This phase implements the I/O infrastructure for SCEC-format output, probe-based time series output, ParaView visualization output, and performs a qualitative serial validation of the BP2 benchmark on a coarse mesh (800 m). The serial test verifies that the code runs stably for at least 5 full earthquake cycles (hundreds of years) and reproduces the fundamental seismic-aseismic cycling behavior. Quantitative comparison against SCEC benchmark data at tight tolerances is deferred to Phase 9, where refined meshes and MPI parallelism provide both the resolution and performance needed for an exact match.

## Dependencies

- **Phase 1-4**: Complete serial implementation
- **BP2 Benchmark Data**: SCEC reference results

## Reference

- **SCEC SEAS Benchmark**: https://strike.scec.org/cvws/seas/
- **BP2 Specification**: https://strike.scec.org/cvws/seas/download/SEAS_BP2_QD.pdf

## Tasks and Deliverables

| Task | Deliverable | Tests |
|------|-------------|-------|
| 5.1 SCEC output format | `io/benchmark_output.hpp` | `tests/unit/test_io.cpp` |
| 5.2 Probe output + interpolation | `io/probe_output.hpp` | `tests/unit/test_io.cpp` |
| 5.3 ParaView visualization output | `io/paraview_output.hpp` | `tests/unit/test_io.cpp` |
| 5.4 BP2 serial smoke test (800 m coarse mesh) | - | (planned) |
| 5.5 Comparison scripts | `compare_results.py` | (planned) |

## Validation Checkpoint

Serial code runs stably for at least 5 earthquake cycles (hundreds of years) on a coarse 800 m mesh, demonstrating seismic-aseismic phase transitions without crashes or solver failures. Qualitative behavior matches expected BP2 patterns.

---

## SCEC Output Format

### Time Series at Stations

Output quantities at each station:
1. Time (s)
2. Slip delta (m)
3. Slip rate V (log10 m/s)
4. Shear stress tau (MPa)
5. State theta (log10 s)

### Station Locations (BP2)

| Station | Depth (km) | Has Benchmark Data |
|---------|------------|-------------------|
| 1 | 0.0 | Yes |
| 2 | 2.4 | No |
| 3 | 4.8 | Yes |
| 4 | 7.2 | No |
| 5 | 9.6 | No |
| 6 | 12.0 | Yes |
| 7 | 14.4 | No |
| 8 | 16.8 | Yes |
| 9 | 19.2 | No |
| 10 | 24.0 | Yes |
| 11 | 28.8 | No |
| 12 | 36.0 | No |

### Output File Format

```
# BP2-QD time series at z = 0 km
# Columns: time(s), slip(m), log10(slip_rate)(m/s), shear_stress(MPa), log10(state)(s)
0.000000e+00  0.000000e+00  -9.000000e+00  2.654600e+01  3.602000e+00
1.000000e+06  1.000000e-03  -9.000000e+00  2.654650e+01  3.602100e+00
...
```

### SCEC Naming Convention

Output files follow SCEC naming convention:

```
{code}_{benchmark}_z{depth}km.txt

Examples:
- mfem_bp2qd_z0km.txt      (z = 0 km)
- mfem_bp2qd_z4.8km.txt    (z = 4.8 km)
- mfem_bp2qd_z12km.txt     (z = 12 km)
```

### Output Frequency

| Phase | Output Interval | Notes |
|-------|-----------------|-------|
| Interseismic | Every 0.01 years | Coarse sampling |
| Nucleation | Every 0.1 seconds | Fine sampling |
| Coseismic | Every 0.001 seconds | Very fine sampling |

Adaptive output based on slip rate:
- V < 10⁻⁶ m/s: Interseismic
- 10⁻⁶ < V < 10⁻³ m/s: Nucleation
- V > 10⁻³ m/s: Coseismic

---

## Detailed Component Design

### 5.1 Benchmark Output

```cpp
// io/benchmark_output.hpp

template <typename MeshType = Mesh>
class BenchmarkOutput {
public:
    /// Constructor: prefix, params, probe depths (m), fault DOF depths (from domain)
    BenchmarkOutput(const std::string &prefix,
                    const BP2Params &params,
                    const std::vector<real_t> &probe_depths,
                    const Vector &fault_depths);

    /// Write if adaptive schedule allows (returns true if written)
    /// Takes fault operator + traction separately (not SEAS operator)
    bool Write(real_t time, const Vector &state,
               const RateStateFaultOperator<MeshType> &fault,
               const Vector &traction);

    /// Force write regardless of schedule
    void ForceWrite(real_t time, const Vector &state,
                    const RateStateFaultOperator<MeshType> &fault,
                    const Vector &traction);

    void Flush();
    void Close();
    int NumProbes() const;

    /// Adaptive output interval based on V_max
    static real_t OutputInterval(real_t V_max);

private:
    ProbeInterpolator interpolator_;  // Linear interpolation to probe depths
    std::vector<std::unique_ptr<ProbeOutput>> probes_;
    real_t last_write_time_;

    /// SCEC filename: "{prefix}_z{depth}km.txt"
    static std::string MakeFilename(const std::string &prefix, real_t depth_m);
};
```

Key design points:
- Template class `<MeshType>` for serial/parallel support
- Takes `RateStateFaultOperator` + `traction` directly (not `SEASQuasiDynamicOperator`)
- Uses `ProbeInterpolator` for linear interpolation between fault DOFs
- One `ProbeOutput` file per probe depth
- `tau_total = tau0 + probe_tau(p)` computed internally
- Values clamped before log10: `V = max(V, 1e-30)`, `theta = max(theta, 1e-30)`

### 5.2 Probe Output and Interpolation

Both `ProbeOutput` and `ProbeInterpolator` are in `io/probe_output.hpp`.

```cpp
// io/probe_output.hpp

/// Single-file ASCII writer for time series data
class ProbeOutput {
public:
    ProbeOutput(const std::string &filename,
                const std::vector<std::string> &column_names);
    ~ProbeOutput();

    /// Write one row of data
    void WriteStep(const std::vector<real_t> &values);

    void Flush();
    void Close();
    bool Good() const;

    // Non-copyable, movable
};

/// Linear interpolation from fault DOFs to probe locations
class ProbeInterpolator {
public:
    /// Precomputes bracket indices and weights from fault_depths to probe_depths
    ProbeInterpolator(const Vector &fault_depths,
                      const std::vector<real_t> &probe_depths);

    /// Interpolate full field to all probes
    void Interpolate(const Vector &field, Vector &probe_values) const;

    /// Interpolate single probe
    real_t InterpolateOne(const Vector &field, int probe_idx) const;

    int NumProbes() const;

private:
    std::vector<int> idx_lo_, idx_hi_;   // Bracket DOF indices
    std::vector<real_t> weight_;         // Interpolation weight (0=lo, 1=hi)
};
```

Key design points:
- `ProbeOutput` writes scientific notation with 12 digits of precision
- Header line: `# Columns: time(s), slip(m), ...`
- `ProbeInterpolator` uses linear interpolation between bracketing DOFs; falls back to nearest-neighbor when probes are outside the DOF range
- Takes `Vector &fault_depths` directly (not `FaultGeometry&`)

### 5.3 ParaView Visualization Output

MFEM provides `ParaViewDataCollection` which writes PVD/VTU files natively supported by ParaView. This enables full-field visualization of the domain displacement, fault slip, stress, and state variable evolution over time.

#### Output Fields

| Field Name | Type | Description |
|------------|------|-------------|
| `displacement` | Scalar GridFunction | Anti-plane displacement u(x,z,t) on the 2D domain |
| `slip` | Scalar GridFunction | Fault slip delta(z,t) on the fault interface |
| `slip_rate` | Scalar GridFunction | Fault slip rate V(z,t) on the fault interface |
| `shear_stress` | Scalar GridFunction | Shear stress tau(z,t) on the fault interface |
| `state` | Scalar GridFunction | State variable theta(z,t) on the fault interface |

#### File Format

ParaView output uses VTK XML formats:
- **`.pvd`** (ParaView Data): Top-level XML file indexing all timesteps, serves as the entry point for ParaView
- **`.vtu`** (VTK Unstructured Grid): Per-timestep mesh and field data in XML format
- Directory structure: `<prefix>/Cycle<NNN>/*.vtu`

These formats are natively supported by ParaView without any plugins or conversion.

#### Output Frequency

ParaView output uses a separate, coarser output schedule than the probe time series to manage file sizes:

| Phase | Output Interval | Notes |
|-------|-----------------|-------|
| Interseismic | Every 1.0 years | Coarse snapshots of slow deformation |
| Nucleation | Every 1.0 seconds | Capture nucleation zone development |
| Coseismic | Every 0.01 seconds | Capture rupture propagation |

Adaptive output based on maximum slip rate (same thresholds as probe output).

#### Class Design

```cpp
// io/paraview_output.hpp

template <typename MeshType = Mesh>
class ParaViewOutput {
public:
    /// Output interval in time steps (write every N steps)
    int output_every_n_steps = 100;

    ParaViewOutput(const std::string &prefix,
                   MeshType &mesh,
                   int order);

    /// Register a domain field (displacement) in the single collection
    void RegisterDomainField(const std::string &name, GridFunction *gf);

    /// Initialize fault output using L2 p=0 on the domain mesh
    /// Maps fault DOF values to adjacent domain elements
    void InitFaultOutput(const Array<int> &fault_interior_faces);

    /// Update fault field data before saving
    void UpdateFaultFields(const Vector &slip, const Vector &slip_rate,
                           const Vector &stress, const Vector &state,
                           int dofs_per_face = 1);

    /// Write if step count matches interval
    bool Save(int cycle, real_t time);

    /// Write based on adaptive time schedule (V-based)
    bool Save(int cycle, real_t time, real_t V_max);

    /// Force write
    void ForceSave(int cycle, real_t time);

    void SetDataFormat(VTKFormat fmt);
    void SetHighOrderOutput(bool enable);
    int GetNumFaultFaces() const;
    bool HasFaultOutput() const;

    static real_t OutputInterval(real_t V_max);

private:
    ParaViewDataCollection pv_;   // Single collection on domain mesh

    // L2 p=0 GridFunctions for fault fields (on domain mesh)
    std::unique_ptr<GridFunction> fault_slip_, fault_slip_rate_;
    std::unique_ptr<GridFunction> fault_stress_, fault_state_;
    std::vector<int> fault_face_elem1_, fault_face_elem2_;
};
```

#### Key Design Decision: Single Collection with L2 p=0

The actual implementation uses **one** `ParaViewDataCollection` on the 2D domain mesh, not two separate collections. Fault fields are stored as **piecewise-constant (L2 p=0) GridFunctions** on the domain mesh:
- Elements adjacent to fault faces carry fault values
- All other elements are zero
- Each fault face maps to two adjacent elements (L2 p=0: DOF index = element index)

This avoids the complexity of a separate 1D fault submesh while still allowing ParaView visualization.

#### Usage Example

```cpp
ParaViewOutput<Mesh> pv("ParaView/bp2", *mesh, order);
pv.RegisterDomainField("displacement", &u_gf);

// Initialize fault output from domain operator's fault faces
pv.InitFaultOutput(domain.GetFaultInteriorFaces());

// In time loop:
pv.UpdateFaultFields(slip, slip_rate, stress, state, /*dofs_per_face=*/1);
pv.Save(cycle, time);  // or pv.ForceSave(cycle, time)
```

#### Implementation Notes

1. **Single data collection**: All fields (domain + fault) go into one `ParaViewDataCollection`. Fault fields use L2 p=0 FE space on the same mesh.

2. **Default settings**: Binary format, high-order output, levels of detail matching FE order.

3. **Parallel support**: `ParaViewDataCollection` works transparently in parallel. No code changes needed for `ParMesh`.

4. **Viewing in ParaView**: Open the `.pvd` file. Apply "Threshold" filter to see only fault-adjacent elements for fault fields.

---

### Implementation Details

```cpp
template <typename MeshType>
bool BenchmarkOutput<MeshType>::Write(
    real_t time, const Vector &state,
    const RateStateFaultOperator<MeshType> &fault,
    const Vector &traction)
{
    // Check adaptive schedule
    real_t V_max = fault.GetMaxSlipRate();
    real_t dt_out = OutputInterval(V_max);
    if (time - last_write_time_ < dt_out * 0.99) { return false; }

    // Extract fault quantities
    Vector slip, theta;
    fault.GetSlip(state, slip);
    fault.GetTheta(state, theta);
    const Vector &slip_rate = fault.GetSlipRate();
    real_t tau0 = fault.GetTau0();

    // Interpolate to all probes
    Vector probe_slip, probe_theta, probe_V, probe_tau;
    interpolator_.Interpolate(slip, probe_slip);
    interpolator_.Interpolate(theta, probe_theta);
    interpolator_.Interpolate(slip_rate, probe_V);
    interpolator_.Interpolate(traction, probe_tau);

    // Write one row per probe
    for (int p = 0; p < interpolator_.NumProbes(); p++)
    {
        real_t V = std::max(probe_V(p), 1e-30);
        real_t th = std::max(probe_theta(p), 1e-30);
        real_t tau_total = tau0 + probe_tau(p);

        probes_[p]->WriteStep({time, probe_slip(p),
                               std::log10(V), tau_total / 1e6,
                               std::log10(th)});
    }

    last_write_time_ = time;
    return true;
}
```

---

## Serial BP2 Smoke Test (Coarse Mesh)

### Mesh and Scope

The serial BP2 test uses the **coarse 800 m mesh** from the available benchmark data. At this resolution the simulation will not match SCEC reference results exactly — the purpose is to verify:

1. **Stability**: The code runs without crashes, solver failures, or NaN/Inf values for the entire simulation span (hundreds of simulated years).
2. **Seismic-aseismic cycling**: The simulation transitions between interseismic (V ~ 10⁻⁹ m/s) and coseismic (V > 10⁻³ m/s) phases, producing recognizable earthquake cycles.
3. **Qualitative behavior**: Slip accumulates, stress builds and drops during earthquakes, and the state variable evolves as expected.
4. **I/O correctness**: All output files (SCEC-format probes, ParaView) are written correctly.

Quantitative comparison against SCEC benchmark data at tight tolerances requires **refined meshes** (e.g., 100-200 m element size), which in turn requires **MPI parallelism** for practical runtimes. This is addressed in **Phase 9**.

### 5.4 BP2 Serial Smoke Test (Planned)

**Status:** Not yet implemented. The `bp2_serial_smoke.cpp` test is planned for when the full simulation pipeline is validated.

The serial smoke test will:
1. Run BP2 on a coarse 800 m mesh for 600+ simulated years
2. Verify at least 5 earthquake cycles (V_max > 10⁻³ m/s)
3. Verify interseismic recovery (V drops below 10⁻⁶ between events)
4. Verify output files are written correctly

### 5.5 Comparison Script (Planned)

**Status:** `compare_results.py` is planned for use with the serial smoke test. It will generate comparison plots of slip, slip rate, shear stress, and state variable time series against SCEC benchmark data.

### 5.6 I/O Unit Tests (Implemented)

All I/O tests are in `tests/unit/test_io.cpp` (7 test functions):

1. **TestProbeOutput** — File creation, header, 3 data rows, cleanup
2. **TestProbeInterpolator** — Linear interpolation at exact and intermediate depths
3. **TestBenchmarkOutputFilenames** — Adaptive output interval computation (interseismic/nucleation/coseismic)
4. **TestBenchmarkOutputIntegration** — Full integration: creates SCEC files, verifies content (5 columns, initial values)
5. **TestParaViewOutput** — Creates PVD file with displacement field
6. **TestParaViewCombinedOutput** — Domain + fault fields in single collection, verifies VTU contains all 5 fields
7. **TestParaViewOutputInterval** — Adaptive ParaView output intervals

---

## Validation Criteria (Coarse 800 m Mesh)

These criteria are intentionally qualitative. The coarse mesh will not match benchmark data exactly, but it must demonstrate correct physical behavior.

### Initial Values (Mesh-Independent)

| Metric | Expected Value | Tolerance |
|--------|----------------|-----------|
| Initial tau | 26.546 MPa | +/- 0.01 MPa |
| Initial theta (z=0) | 10^3.602 s | +/- 10% |
| Initial V | 10^-9 m/s | exact |

### Stability and Cycling (Primary Goal)

| Metric | Criterion |
|--------|-----------|
| Simulation duration | Runs for 600+ simulated years without crash |
| Solver convergence | No solver failures or divergence |
| Numerical health | No NaN or Inf values at any time step |
| Earthquake detection | At least 5 seismic events (V_max > 10⁻³ m/s) |
| Interseismic recovery | V_max drops below 10⁻⁶ m/s between events |

### Qualitative Time Evolution (Loose Tolerances)

These values are approximate on the coarse mesh. Exact match is **not** required.

| Metric | Expected Range | Notes |
|--------|----------------|-------|
| First earthquake time | 50-150 years | Coarse mesh may shift timing |
| Earthquake recurrence | 50-150 years | Qualitative periodicity sufficient |
| Maximum slip rate | 10⁻⁴ to 10¹ m/s | Must reach coseismic rates |
| Slip per event | 0.5-10 m | Order-of-magnitude check |

### Quantitative Comparison (Deferred to Phase 9)

Tight quantitative comparison against SCEC benchmark data requires refined meshes (100-200 m element size) and MPI parallelism. The following metrics will be evaluated in Phase 9:

| Quantity | Error Metric | Tolerance | Required Mesh |
|----------|--------------|-----------|---------------|
| Slip | Relative error | < 5% | Refined (MPI) |
| log10(V) | Absolute error | < 0.5 | Refined (MPI) |
| Shear stress | Absolute error | < 0.5 MPa | Refined (MPI) |
| log10(theta) | Absolute error | < 0.5 | Refined (MPI) |

---

## File Organization

```
io/
├── benchmark_output.hpp         # SCEC-format output (template, uses ProbeInterpolator)
├── probe_output.hpp             # ProbeOutput (ASCII writer) + ProbeInterpolator (linear interp)
└── paraview_output.hpp          # ParaView PVD/VTU output (single collection, L2 p=0 fault fields)

tests/unit/
└── test_io.cpp                  # 7 I/O unit tests
```

---

## Acceptance Criteria

### I/O Unit Tests (Implemented)

All 7 tests in `tests/unit/test_io.cpp` must pass:

| Test | Criterion |
|------|-----------|
| TestProbeOutput | File created, header written, 3 data rows with 12-digit precision |
| TestProbeInterpolator | Exact values at node depths, correct linear interpolation between nodes |
| TestBenchmarkOutputFilenames | Correct adaptive output interval (interseismic/nucleation/coseismic thresholds) |
| TestBenchmarkOutputIntegration | SCEC probe files created with 5-column format and correct initial values |
| TestParaViewOutput | PVD file created with displacement field |
| TestParaViewCombinedOutput | Single VTU contains domain displacement + 4 fault fields (slip, slip_rate, stress, state) |
| TestParaViewOutputInterval | Adaptive ParaView output intervals respond to V_max thresholds |

### BP2 Serial Smoke Test (Planned)

| Metric | Criterion |
|--------|-----------|
| Initial tau | 26.546 MPa +/- 0.01 MPa |
| Initial theta (z=0) | 10^3.602 s +/- 10% |
| Initial V | 10⁻⁹ m/s (exact) |
| Long-term stability | Runs 600+ years without crash, NaN, or solver failure |
| Seismic-aseismic cycling | At least 5 earthquake events detected |
| Interseismic recovery | Slip rate returns to < 10⁻⁶ m/s between events |
| I/O correctness | SCEC probe files and ParaView PVD/VTU files written |

### Phase 9 (Parallel, Refined Mesh — Quantitative Match)

| Metric | Criterion |
|--------|-----------|
| Slip time series | Relative error < 5% vs SCEC benchmark |
| log10(V) time series | Absolute error < 0.5 vs SCEC benchmark |
| Shear stress time series | Absolute error < 0.5 MPa vs SCEC benchmark |
| log10(theta) time series | Absolute error < 0.5 vs SCEC benchmark |
| Serial-parallel consistency | Parallel matches serial within 1e-10 |
