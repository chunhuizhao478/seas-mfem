# Phase 9: Parallel Verification and Scaling

## Overview

This phase verifies parallel correctness through serial-parallel consistency tests, evaluates parallel performance through scaling studies, and performs the **quantitative validation** of BP2 benchmark results on refined meshes. Phase 5 validated the serial code qualitatively on a coarse 800 m mesh (stability, seismic-aseismic cycling). Here we use MPI parallelism to run on refined meshes (100-200 m) where tight comparison against SCEC benchmark data is meaningful.

## Dependencies

- **Phase 1-5**: Complete serial implementation (reference for comparison)
- **Phase 6-8**: Complete parallel implementation

## Reference

- **Phase 5** (`phase5_io_validation.md`): Serial validation criteria (same for parallel)

## Tasks and Deliverables

| Task | Deliverable | Tests |
|------|-------------|-------|
| 9.1 Serial-parallel consistency | - | `test_serial_parallel_consistency.cpp` |
| 9.2 BP2 parallel validation (refined mesh) | - | `bp2_benchmark_parallel.cpp` |
| 9.3 Strong scaling test | - | `test_scaling.cpp` |
| 9.4 Weak scaling test | - | Included above |
| 9.5 Communication profiling | - | Timing reports |

## Verification Checkpoint

- Parallel results match serial within tolerance
- Reasonable parallel efficiency (>70% at 8 ranks)
- **Quantitative match** against SCEC BP2 benchmark data on refined mesh (100-200 m element size)

---

## Serial-Parallel Consistency

### Consistency Requirements

The parallel implementation must produce results identical (within numerical tolerance) to the serial implementation:

| Component | Consistency Requirement |
|-----------|------------------------|
| Domain solution | Parallel matches serial within 1e-10 |
| Fault traction | Parallel matches serial within 1e-10 |
| Slip rate | Parallel matches serial within 1e-10 |
| State evolution | Identical after N time steps |
| Output | Probe values identical |

### Consistency Test

```cpp
// tests/parallel/test_serial_parallel_consistency.cpp

class SerialParallelConsistencyTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
        MPI_Comm_size(MPI_COMM_WORLD, &size_);
    }
    int rank_, size_;
};

/// Test 1: Domain solution consistency
TEST_F(SerialParallelConsistencyTest, DomainSolution) {
    BP2Params params;
    int N = 32;

    // Compute serial solution on rank 0
    Vector serial_u;
    if (rank_ == 0) {
        Mesh mesh = CreateBP2Mesh(params, N);
        AntiplaneDomainOperator<Mesh> domain(mesh, 1, params.mu);

        Vector slip(domain.NumFaultDOFs());
        slip = 0.001;  // Some test slip

        GridFunction u(&domain.GetFESpace());
        domain.Solve(0.0, slip, u);
        serial_u = u;
    }

    // Compute parallel solution
    Mesh serial_mesh = CreateBP2Mesh(params, N);
    ParMesh pmesh(MPI_COMM_WORLD, serial_mesh);
    AntiplaneDomainOperator<ParMesh> par_domain(pmesh, 1, params.mu);

    Vector slip(par_domain.NumLocalFaultDOFs());
    slip = 0.001;

    ParGridFunction par_u(&par_domain.GetFESpace());
    par_domain.Solve(0.0, slip, par_u);

    // Gather parallel solution to root and compare
    Vector gathered_u;
    par_domain.GatherSolutionToRoot(par_u, gathered_u);

    if (rank_ == 0) {
        real_t max_diff = 0.0;
        for (int i = 0; i < serial_u.Size(); i++) {
            real_t diff = std::abs(serial_u(i) - gathered_u(i));
            max_diff = std::max(max_diff, diff);
        }
        EXPECT_LT(max_diff, 1e-10 * serial_u.Norml2());
    }
}

/// Test 2: Fault state consistency after time stepping
TEST_F(SerialParallelConsistencyTest, FaultState) {
    BP2Params params;
    params.t_final = 1e6;  // Short run

    // Run serial simulation (on rank 0)
    Vector serial_state;
    if (rank_ == 0) {
        // ... run serial ...
        serial_state = final_state;
    }

    // Run parallel simulation
    Vector parallel_state = RunParallelSimulation(params);

    // Gather and compare
    Vector gathered_state;
    GatherStateToRoot(parallel_state, gathered_state);

    if (rank_ == 0) {
        for (int i = 0; i < serial_state.Size(); i++) {
            EXPECT_NEAR(serial_state(i), gathered_state(i),
                        1e-10 * std::abs(serial_state(i)));
        }
    }
}

/// Test 3: Output probe values consistency
TEST_F(SerialParallelConsistencyTest, ProbeOutput) {
    // Run both serial and parallel for same time period
    // Compare probe output values

    Vector serial_probes = RunSerialAndGetProbes();
    Vector parallel_probes = RunParallelAndGetProbes();

    if (rank_ == 0) {
        for (int i = 0; i < serial_probes.Size(); i++) {
            EXPECT_NEAR(serial_probes(i), parallel_probes(i), 1e-12);
        }
    }
}

/// Test 4: Reproducibility across different rank counts
TEST_F(SerialParallelConsistencyTest, Reproducibility) {
    // Results should be identical regardless of number of ranks used
    // This is verified by running same problem with np=1,2,4,8
}
```

---

## BP2 Parallel Validation (Refined Mesh — Quantitative Comparison)

### Context

Phase 5 validated the serial implementation qualitatively on a coarse 800 m mesh, confirming stability and seismic-aseismic cycling. Here in Phase 9 we perform the **quantitative comparison** against SCEC BP2 benchmark data using a **refined mesh** (100-200 m element size). The refined mesh requires MPI parallelism for practical runtimes.

### Mesh Resolution and Convergence

| Mesh Size | Fault DOFs | Purpose |
|-----------|-----------|---------|
| 800 m | ~50 | Phase 5 serial smoke test (qualitative) |
| 200 m | ~200 | Intermediate convergence check |
| 100 m | ~400 | Quantitative benchmark comparison target |

A mesh convergence study should be performed to confirm that the 100 m result is in the asymptotic regime.

### Validation Against SCEC Benchmark

```cpp
// tests/verification/bp2_benchmark_parallel.cpp

/// Quantitative comparison on refined mesh against SCEC benchmark data.
/// This test requires MPI and runs on a refined mesh (100-200 m element size).
TEST(BP2ParallelBenchmark, QuantitativeMatchRefinedMesh) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Run full BP2 simulation on refined mesh in parallel
    BP2Params params;
    params.mesh_size = 100.0;       // 100 m refined mesh
    params.t_final = 1e10;          // ~317 years
    auto result = RunBP2Parallel(params);

    // Load benchmark data and compare (on root only)
    if (rank == 0) {
        auto bench_z0 = LoadBenchmarkData("bp2-qd-z0km-res.txt");
        auto bench_z4_8 = LoadBenchmarkData("bp2-qd-z4.8km-res.txt");
        auto bench_z12 = LoadBenchmarkData("bp2-qd-z12km-res.txt");
        auto bench_z16_8 = LoadBenchmarkData("bp2-qd-z16.8km-res.txt");
        auto bench_z24 = LoadBenchmarkData("bp2-qd-z24km-res.txt");

        // Compare at selected times
        std::vector<real_t> compare_times = {1e6, 1e7, 1e8, 1e9, 1e10};

        for (real_t t : compare_times) {
            // Compare slip (relative error < 5%)
            real_t slip_sim = InterpolateSlip(result, t, 0.0);
            real_t slip_bench = InterpolateSlip(bench_z0, t);
            EXPECT_NEAR(slip_sim, slip_bench, std::abs(slip_bench) * 0.05);

            // Compare slip rate (within half decade)
            real_t V_sim = InterpolateSlipRate(result, t, 0.0);
            real_t V_bench = InterpolateSlipRate(bench_z0, t);
            EXPECT_NEAR(std::log10(V_sim), std::log10(V_bench), 0.5);

            // Compare stress (within 0.5 MPa)
            real_t tau_sim = InterpolateStress(result, t, 0.0);
            real_t tau_bench = InterpolateStress(bench_z0, t);
            EXPECT_NEAR(tau_sim, tau_bench, 0.5);
        }
    }
}

/// Verify parallel results match serial BP2 results (on same mesh)
TEST(BP2ParallelBenchmark, MatchSerialResults) {
    // Run serial BP2 on rank 0 (coarse mesh for tractability)
    // Run parallel BP2 on all ranks (same coarse mesh)
    // Compare final state and time series — must match within 1e-10
}

/// Mesh convergence study: verify results converge as mesh is refined
TEST(BP2ParallelBenchmark, MeshConvergence) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Run at multiple resolutions
    std::vector<real_t> mesh_sizes = {800.0, 400.0, 200.0, 100.0};
    std::vector<real_t> first_eq_times;

    for (real_t h : mesh_sizes) {
        BP2Params params;
        params.mesh_size = h;
        params.t_final = 1e10;
        auto result = RunBP2Parallel(params);

        if (rank == 0) {
            first_eq_times.push_back(result.first_earthquake_time);
        }
    }

    // Verify convergence: successive refinements should yield
    // smaller changes in first earthquake time
    if (rank == 0) {
        for (size_t i = 2; i < first_eq_times.size(); i++) {
            real_t diff_coarse = std::abs(first_eq_times[i-1] - first_eq_times[i-2]);
            real_t diff_fine   = std::abs(first_eq_times[i]   - first_eq_times[i-1]);
            EXPECT_LT(diff_fine, diff_coarse);  // Converging
        }
    }
}
```

---

## Scaling Tests

### Strong Scaling

**Definition**: Fixed problem size, increasing processor count.

**Goal**: Time decreases proportionally with processor count.

```cpp
// tests/parallel/test_scaling.cpp

/// Strong scaling test: fixed problem size, increasing processors
TEST(Scaling, StrongScaling) {
    // Problem: 128x128 mesh, 1-year simulation
    // Run with np = 1, 2, 4, 8, 16

    BP2Params params;
    params.nx = 128;
    params.nz = 128;
    params.t_final = 3.15e7;  // 1 year

    // Measure wall time
    double start = MPI_Wtime();
    RunSimulation(params);
    double end = MPI_Wtime();
    double elapsed = end - start;

    // Collect times from all ranks
    double max_time;
    MPI_Reduce(&elapsed, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // On root, compute and report speedup
    if (rank == 0) {
        // Reference: time with 1 process
        double T1 = reference_serial_time;
        double Tn = max_time;
        double speedup = T1 / Tn;
        double efficiency = speedup / size;

        std::cout << "np=" << size
                  << " time=" << Tn
                  << " speedup=" << speedup
                  << " efficiency=" << efficiency << std::endl;

        // Target: >70% efficiency at 8 processes
        if (size == 8) {
            EXPECT_GT(efficiency, 0.70);
        }
    }
}
```

### Expected Strong Scaling Results

| Ranks | Time (s) | Speedup | Efficiency |
|-------|----------|---------|------------|
| 1 | T1 | 1.0 | 100% |
| 2 | T1/1.8 | 1.8 | 90% |
| 4 | T1/3.2 | 3.2 | 80% |
| 8 | T1/5.6 | 5.6 | 70% |
| 16 | T1/9.6 | 9.6 | 60% |

### Weak Scaling

**Definition**: Problem size scales with processor count (fixed work per processor).

**Goal**: Time remains constant as problem and processors scale together.

```cpp
/// Weak scaling test: problem size scales with processors
TEST(Scaling, WeakScaling) {
    // Each rank handles 32x32 mesh elements
    // np = 1:  32x32
    // np = 4:  64x64
    // np = 16: 128x128

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Scale mesh with rank count
    int N = 32 * std::sqrt(size);

    BP2Params params;
    params.nx = N;
    params.nz = N;

    // Run and measure time
    double start = MPI_Wtime();
    RunSimulation(params);
    double end = MPI_Wtime();
    double elapsed = end - start;

    // Collect times
    double max_time;
    MPI_Reduce(&elapsed, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "np=" << size
                  << " N=" << N
                  << " time=" << max_time << std::endl;

        // Target: <1.5x increase in time at 16 processes
        if (size == 16 && reference_time_np1 > 0) {
            EXPECT_LT(max_time, reference_time_np1 * 1.5);
        }
    }
}
```

### Expected Weak Scaling Results

| Ranks | Problem Size | Time (s) | Overhead |
|-------|-------------|----------|----------|
| 1 | 32x32 | T | 0% |
| 4 | 64x64 | T*1.1 | 10% |
| 16 | 128x128 | T*1.3 | 30% |

### Communication Overhead Test

```cpp
/// Communication overhead test
TEST(Scaling, CommunicationOverhead) {
    // Measure time spent in MPI operations

    double mpi_time = 0.0;
    double total_time = 0.0;

    // Instrument key MPI calls
    double start = MPI_Wtime();

    for (int step = 0; step < num_steps; step++) {
        // Domain solve (includes MPI in HYPRE)
        double mpi_start = MPI_Wtime();
        // ... MPI operations ...
        mpi_time += MPI_Wtime() - mpi_start;
    }

    total_time = MPI_Wtime() - start;

    double comm_fraction = mpi_time / total_time;

    if (rank == 0) {
        std::cout << "Communication fraction: " << comm_fraction * 100 << "%" << std::endl;

        // Target: <10% of total time
        EXPECT_LT(comm_fraction, 0.10);
    }
}
```

---

## Profiling and Timing

### Key Metrics to Measure

| Metric | Description | Tool |
|--------|-------------|------|
| Wall time | Total execution time | MPI_Wtime |
| Solver time | Time in linear solver | Profiling |
| MPI time | Time in MPI calls | MPI profiling |
| Memory | Memory per rank | System tools |

### Profiling Commands

```bash
# Run with MPI profiling
mpirun -np 4 ./pseas -bp bp2 -tf 1e10 --profile

# Use Score-P or similar
scorep mpirun -np 4 ./pseas -bp bp2

# Memory profiling
valgrind --tool=massif mpirun -np 2 ./pseas -bp bp2 -tf 1e8
```

---

## Acceptance Criteria

### Parallel Unit Tests

| Test | Pass Criteria |
|------|---------------|
| Mesh distribution | All ranks have elements, global count preserved |
| DOF distribution | Unique DOF ownership, no gaps or duplicates |
| Parallel solver | Converges to same tolerance as serial |
| Communication | Gather/scatter operations preserve data |

### Parallel Consistency Tests

| Test | Pass Criteria |
|------|---------------|
| Solution match | Parallel matches serial within 1e-10 relative error |
| Traction match | Parallel matches serial within 1e-10 |
| State evolution | Parallel state matches serial after N time steps |
| Output match | Probe values identical between serial and parallel |

### Parallel Scaling Criteria

| Metric | Target |
|--------|--------|
| Strong scaling efficiency (8 ranks) | > 70% |
| Weak scaling overhead (16 ranks) | < 50% increase in time |
| Communication overhead | < 10% of total time |
| Memory per rank | Decreases linearly with rank count |

### Parallel Validation Criteria

| Metric | Pass Criteria |
|--------|---------------|
| BP2 parallel results | Match serial results within 1e-10 (same mesh) |
| Reproducibility | Same result with different rank counts (1,2,4,8) |

### Quantitative BP2 Benchmark Comparison (Refined Mesh)

These tight tolerances apply to the refined mesh (100-200 m) run with MPI. This is the definitive validation that was deferred from Phase 5 (coarse 800 m serial smoke test).

| Quantity | Error Metric | Tolerance |
|----------|--------------|-----------|
| Slip | Relative error | < 5% |
| log10(V) | Absolute error | < 0.5 |
| Shear stress | Absolute error | < 0.5 MPa |
| log10(theta) | Absolute error | < 0.5 |
| Mesh convergence | First earthquake time converges | Successive diffs decrease |

---

## File Organization

```
tests/
├── parallel/
│   ├── test_parallel_domain.cpp
│   ├── test_parallel_fault.cpp
│   ├── test_parallel_comm.cpp
│   ├── test_serial_parallel_consistency.cpp
│   └── test_scaling.cpp
└── verification/
    ├── mms_antiplane_parallel.cpp
    └── bp2_benchmark_parallel.cpp
```
