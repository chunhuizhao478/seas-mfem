# Phase 8: Parallel Fault and I/O

## Overview

This phase extends the fault operator and I/O components to support parallel (MPI) execution, including distributed fault data management and parallel probe output.

## Dependencies

- **Phase 3** (`phase3_fault_operator.md`): Serial fault operator
- **Phase 5** (`phase5_io_validation.md`): Serial I/O and output format
- **Phase 6** (`phase6_parallel_infrastructure.md`): MPI context and utilities
- **Phase 7** (`phase7_parallel_domain.md`): Parallel domain operator

## Tasks and Deliverables

| Task | Deliverable | Tests |
|------|-------------|-------|
| 8.1 Parallel fault data | `fault/parallel_fault_data.hpp` | `test_parallel_fault.cpp` |
| 8.2 Distributed state mgmt | In `rate_state_fault.hpp` | - |
| 8.3 Parallel time stepper | `solver/parallel_time_stepper.hpp` | - |
| 8.4 Parallel probe output | `io/parallel_probe_output.hpp` | `test_parallel_comm.cpp` |
| 8.5 Parallel main driver | `pseas.cpp` | - |

## Verification Checkpoint

Parallel code runs on multiple ranks.

---

## Fault DOF Distribution

In parallel execution, fault DOFs are distributed based on mesh partitioning:

```
+-------------------------------------------------------------------------+
|                    Fault DOF Distribution Example                        |
|                    (Mesh partitioned by depth)                           |
+-------------------------------------------------------------------------+
|                                                                          |
|  Depth (z)    0 km                                            40 km     |
|     |         |                                                 |        |
|     v         v                                                 v        |
|  +----------+----------+----------+----------+----------+----------+    |
|  |  Rank 0  |  Rank 1  |  Rank 2  |  Rank 3  |  Rank 4  |  Rank 5  |    |
|  | (0-7km)  | (7-14km) |(14-20km) |(20-26km) |(26-33km) |(33-40km) |    |
|  |  ~10 DOFs|  ~10 DOFs|  ~10 DOFs|  ~10 DOFs|  ~10 DOFs|  ~10 DOFs|    |
|  +----------+----------+----------+----------+----------+----------+    |
|                                                                          |
|  Key Points:                                                             |
|  - Each rank owns fault DOFs on its local mesh portion                   |
|  - State variables (theta) are computed locally (no communication)       |
|  - Slip rates (V) are computed locally from local traction               |
|  - Global max V requires MPI_Allreduce for time step control             |
|  - Output stations may span ranks (need gather for I/O)                  |
|                                                                          |
|  Output Stations (may not align with rank boundaries):                   |
|  z = 0 km    -> Rank 0 (local)                                          |
|  z = 4.8 km  -> Rank 0 (local)                                          |
|  z = 12 km   -> Rank 1 (local)                                          |
|  z = 16.8 km -> Rank 2 (local)                                          |
|  z = 24 km   -> Rank 3 (local)                                          |
|                                                                          |
+--------------------------------------------------------------------------+
```

---

## Detailed Component Design

### 8.1 Parallel Fault Data

```cpp
// fault/parallel_fault_data.hpp

namespace mfem {
namespace seas {

/// Manages distribution of fault data across MPI ranks
class ParallelFaultData {
public:
    /// Construct from parallel mesh and fault boundary attribute
    ParallelFaultData(ParMesh &pmesh, int fault_attr);

    /// Get number of fault DOFs owned by this rank
    int NumLocalDOFs() const { return num_local_; }

    /// Get total number of fault DOFs globally
    int NumGlobalDOFs() const { return num_global_; }

    /// Get local-to-global DOF mapping
    const Array<int> &LocalToGlobal() const { return l2g_; }

    /// Get global-to-local DOF mapping (-1 if not local)
    int GlobalToLocal(int global_dof) const;

    /// Check if a global DOF is owned by this rank
    bool IsOwned(int global_dof) const;

    /// Get depth (z-coordinate) for each local fault DOF
    const Vector &GetLocalDepths() const { return local_depths_; }

    /// Gather local data to global array on root rank
    /// Used for I/O operations
    void GatherToRoot(const Vector &local_data, Vector &global_data) const;

    /// Scatter global data from root to local arrays
    void ScatterFromRoot(const Vector &global_data, Vector &local_data) const;

    /// Interpolate global data to specific probe depths (on root only)
    void InterpolateToProbes(const Vector &global_data,
                              const Array<real_t> &probe_depths,
                              Vector &probe_values) const;

private:
    ParMesh &pmesh_;
    int fault_attr_;

    int num_local_;
    int num_global_;
    Array<int> l2g_;         // Local to global mapping
    Array<int> g2l_;         // Global to local (-1 if not owned)
    Array<int> owner_;       // Which rank owns each global DOF
    Vector local_depths_;

    // For MPI communication
    Array<int> recv_counts_;
    Array<int> recv_displs_;

    void SetupDistribution();
};

} // namespace seas
} // namespace mfem
```

### 8.2 Parallel Fault Operator Extension

```cpp
// fault/rate_state_fault.hpp (with parallel support)

class RateStateFaultOperator {
public:
    RateStateFaultOperator(FaultGeometry *geom,
                           FrictionLaw *friction,
                           StateEvolution *state_law,
                           const BP2Params &params,
                           const MPIContext *mpi_ctx = nullptr);

    /// Number of LOCAL state DOFs (slip + state variable)
    /// In parallel, this is the local portion only
    int LocalStateSize() const { return num_local_fault_dofs_ * 2; }

    /// Number of GLOBAL state DOFs
    int GlobalStateSize() const { return num_global_fault_dofs_ * 2; }

    /// Get LOCAL maximum slip rate
    real_t GetLocalMaxSlipRate() const { return V_max_local_; }

    /// Get GLOBAL maximum slip rate (requires MPI reduction in parallel)
    real_t GetGlobalMaxSlipRate() const {
        if (mpi_ctx_) {
            return mpi_ctx_->GlobalMax(V_max_local_);
        }
        return V_max_local_;
    }

    /// Access local fault DOF information
    int NumLocalDOFs() const { return num_local_fault_dofs_; }

    /// Get local-to-global mapping for fault DOFs (parallel only)
    const Array<int> &GetLocalToGlobalMap() const { return local_to_global_; }

private:
    const MPIContext *mpi_ctx_;

    int num_local_fault_dofs_;   // Local fault DOFs on this rank
    int num_global_fault_dofs_;  // Total fault DOFs across all ranks
    real_t V_max_local_;

    // Parallel data distribution
    Array<int> local_to_global_;
    Array<int> fault_dof_owner_;
};
```

### 8.3 Parallel Time Stepper

```cpp
// solver/parallel_time_stepper.hpp

namespace mfem {
namespace seas {

/// Adaptive time stepper aware of parallel execution
/// Synchronizes time step across all MPI ranks
class ParallelAdaptiveTimeStepper {
public:
    ParallelAdaptiveTimeStepper(const MPIContext *mpi_ctx = nullptr);

    /// Configuration
    void SetDtMin(real_t dt_min) { dt_min_ = dt_min; }
    void SetDtMax(real_t dt_max) { dt_max_ = dt_max; }
    void SetTargetSlipRateMax(real_t V_target) { V_target_ = V_target; }

    /// Compute new time step based on slip rate
    /// In parallel, uses global max slip rate
    real_t ComputeNewDt(real_t V_max_local, real_t dt_current) const {
        // Get global max slip rate
        real_t V_max = mpi_ctx_ ? mpi_ctx_->GlobalMax(V_max_local)
                                : V_max_local;

        // Adaptive formula based on slip rate
        real_t dt_new;
        if (V_max < 1e-12) {
            dt_new = dt_max_;
        } else if (V_max > V_target_) {
            dt_new = dt_current * V_target_ / V_max;
        } else {
            dt_new = std::min(dt_current * 1.5, dt_max_);
        }

        return std::max(dt_min_, std::min(dt_max_, dt_new));
    }

    /// Get current time step (same on all ranks)
    real_t GetDt() const { return dt_; }

    /// Update time step (must be called synchronously by all ranks)
    void SetDt(real_t dt) { dt_ = dt; }

private:
    const MPIContext *mpi_ctx_;
    real_t dt_min_ = 1e-6;
    real_t dt_max_ = 3.15e7;
    real_t V_target_ = 1e-6;
    real_t dt_ = 1e3;
};

} // namespace seas
} // namespace mfem
```

### 8.4 Parallel Probe Output

```cpp
// io/parallel_probe_output.hpp

namespace mfem {
namespace seas {

/// Parallel-aware probe output
/// Gathers data from multiple ranks before writing
class ParallelProbeOutput {
public:
    ParallelProbeOutput(const std::string &filename,
                        const std::vector<real_t> &probe_depths,
                        const ParallelFaultData &fault_data,
                        const MPIContext &mpi_ctx);

    /// Write output at current time
    /// Only root rank actually writes to file
    void Write(real_t time, const Vector &local_state,
               const Vector &local_slip_rate,
               const Vector &local_stress);

private:
    std::ofstream file_;  // Only valid on root
    std::vector<real_t> probe_depths_;
    const ParallelFaultData &fault_data_;
    const MPIContext &mpi_ctx_;

    // Work arrays for gathering
    Vector global_state_;
    Vector global_slip_rate_;
    Vector global_stress_;

    /// Gather and interpolate to probe locations
    void GatherAndInterpolate(const Vector &local_data,
                               Vector &probe_values);
};

} // namespace seas
} // namespace mfem
```

### 8.5 Parallel Main Driver

```cpp
// pseas.cpp

#include "common/seas_types.hpp"
#include "common/mpi_context.hpp"
#include "domain/antiplane_operator.hpp"
#include "fault/rate_state_fault.hpp"
#include "solver/seas_operator.hpp"
#include "io/parallel_probe_output.hpp"

int main(int argc, char *argv[]) {
    // Initialize MPI context
    mfem::seas::MPIContext mpi_ctx(&argc, &argv);

    // Parse command line
    mfem::OptionsParser args(argc, argv);
    // ... add options ...

    // Create serial mesh on rank 0, then distribute
    Mesh serial_mesh;
    if (mpi_ctx.IsRoot()) {
        serial_mesh = CreateBP2Mesh(params);
    }
    ParMesh pmesh(MPI_COMM_WORLD, serial_mesh);

    // Create operators using parallel types
    AntiplaneDomainOperator<ParMesh> domain(pmesh, order, params.mu);

    RateStateFaultOperator fault(&geom, &friction, &aging_law, params, &mpi_ctx);

    SEASQuasiDynamicOperator<ParMesh> seas_op(&domain, &fault, &mpi_ctx);

    // Time integration
    ParallelAdaptiveTimeStepper stepper(&mpi_ctx);
    RK4Solver ode_solver;
    ode_solver.Init(seas_op);

    Vector state(fault.LocalStateSize());
    seas_op.SetInitialCondition(state);

    // Output (parallel-aware)
    ParallelProbeOutput output("bp2_qd", probe_depths, fault_data, mpi_ctx);

    real_t t = 0.0, t_final = params.t_final;
    while (t < t_final) {
        real_t dt = stepper.GetDt();
        ode_solver.Step(state, t, dt);
        t += dt;

        // Update time step based on global max slip rate
        real_t V_max = seas_op.GetMaxSlipRate();  // Uses MPI reduction
        stepper.SetDt(stepper.ComputeNewDt(V_max, dt));

        // Output (gathers data from all ranks, writes on root)
        output.Write(t, state, fault.GetSlipRate(), fault.GetShearStress());
    }

    return 0;
}
```

---

## Parallel Unit Tests

### Parallel Fault Tests

```cpp
// tests/parallel/test_parallel_fault.cpp

/// Test 1: Fault DOF distribution
TEST_F(ParallelFaultTest, DOFDistribution) {
    ParMesh pmesh(MPI_COMM_WORLD, serial_mesh);
    ParallelFaultData fault_data(pmesh, fault_attr);

    // Verify global count equals sum of local counts
    int local_dofs = fault_data.NumLocalDOFs();
    int global_dofs;
    MPI_Allreduce(&local_dofs, &global_dofs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    EXPECT_EQ(global_dofs, fault_data.NumGlobalDOFs());
}

/// Test 2: Local-to-global mapping is unique
TEST_F(ParallelFaultTest, UniqueMapping) {
    ParallelFaultData fault_data(pmesh, fault_attr);

    // Gather all local-to-global maps to root
    std::vector<int> all_global_dofs;
    // ... gather ...

    if (rank_ == 0) {
        std::set<int> unique_dofs(all_global_dofs.begin(), all_global_dofs.end());
        EXPECT_EQ(unique_dofs.size(), fault_data.NumGlobalDOFs());
    }
}

/// Test 3: Global max slip rate reduction
TEST_F(ParallelFaultTest, GlobalMaxSlipRate) {
    // Set known slip rates on each rank
    real_t local_V_max;
    if (rank_ == 0) local_V_max = 1e-9;
    else if (rank_ == 1) local_V_max = 1e-6;  // Should be global max
    else local_V_max = 1e-8;

    real_t global_V_max = mpi_ctx.GlobalMax(local_V_max);

    EXPECT_NEAR(global_V_max, 1e-6, 1e-12);
}

/// Test 4: Gather to root
TEST_F(ParallelFaultTest, GatherToRoot) {
    ParallelFaultData fault_data(pmesh, fault_attr);

    // Create local data = global DOF index
    Vector local_data(fault_data.NumLocalDOFs());
    for (int i = 0; i < local_data.Size(); i++) {
        local_data(i) = fault_data.LocalToGlobal()[i];
    }

    // Gather to root
    Vector global_data;
    fault_data.GatherToRoot(local_data, global_data);

    // On root, verify values
    if (rank_ == 0) {
        EXPECT_EQ(global_data.Size(), fault_data.NumGlobalDOFs());
        for (int i = 0; i < global_data.Size(); i++) {
            EXPECT_NEAR(global_data(i), i, 1e-12);
        }
    }
}
```

### Parallel Communication Tests

```cpp
// tests/parallel/test_parallel_comm.cpp

/// Test 1: Gather operation
TEST_F(ParallelCommTest, GatherVector) {
    // Test GatherVectorToRoot
}

/// Test 2: Scatter operation
TEST_F(ParallelCommTest, ScatterVector) {
    // Test ScatterVectorFromRoot
}

/// Test 3: Probe interpolation in parallel
TEST_F(ParallelCommTest, ProbeInterpolation) {
    // Verify that output at probe stations is correctly gathered
    // from potentially different ranks
}
```

---

## Communication Patterns

### Required MPI Communications

| Operation | Communication | Frequency |
|-----------|--------------|-----------|
| Time step control | MPI_Allreduce (max V) | Every time step |
| Probe output | MPI_Gatherv | Output frequency |
| Convergence check | MPI_Allreduce (residual) | Solver iterations |

### Communication Overhead Minimization

1. **Fault state evolution**: No communication (local computation)
2. **Slip rate computation**: No communication (local from local traction)
3. **Time step**: Single scalar reduction per step
4. **Output**: Only at specified output times

---

## File Organization

```
fault/
├── fault_geometry.hpp           # Fault surface geometry
├── rate_state_fault.hpp         # Rate-state fault operator
├── fault_output.hpp             # Fault output handling
└── parallel_fault_data.hpp      # Distributed fault data management

solver/
├── seas_operator.hpp            # TimeDependentOperator for SEAS
├── time_stepper.hpp             # Adaptive time stepping
└── parallel_time_stepper.hpp    # Parallel-aware time control

io/
├── benchmark_output.hpp         # SCEC-format output
├── probe_output.hpp             # Serial probe output
└── parallel_probe_output.hpp    # MPI-aware probe gathering
```

---

## Acceptance Criteria

| Test | Pass Criteria |
|------|---------------|
| DOF distribution | Global count = sum of local counts |
| Unique ownership | Each global DOF owned by exactly one rank |
| Global max reduction | Returns correct maximum across ranks |
| Gather operation | Data correctly collected on root |
| Time step sync | All ranks use same dt |
| Output | Probe values match serial execution |
