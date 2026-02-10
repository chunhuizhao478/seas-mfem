# Phase 6: Parallel Infrastructure

## Overview

This phase establishes the foundational infrastructure for parallel (MPI) execution, including type abstractions, MPI context wrappers, and CMake configuration for parallel builds.

## Dependencies

- **Phase 1-5**: Complete serial implementation (must pass all serial tests first)

## Dependent Phases

- **Phase 7** (`phase7_parallel_domain.md`): Uses type abstractions
- **Phase 8** (`phase8_parallel_fault_io.md`): Uses MPI context and parallel utilities
- **Phase 9** (`phase9_parallel_verification.md`): Uses all parallel components

## Tasks and Deliverables

| Task | Deliverable | Tests |
|------|-------------|-------|
| 6.1 Type abstraction | `common/seas_types.hpp` | - |
| 6.2 MPI context wrapper | `common/mpi_context.hpp` | `test_mpi_context.cpp` |
| 6.3 Parallel utilities | `common/parallel_utils.hpp` | `test_parallel_utils.cpp` |
| 6.4 CMake parallel config | `CMakeLists.txt` updates | Build test |

## Verification Checkpoint

Code compiles with and without `MFEM_USE_MPI`.

---

## Design Philosophy

### Serial/Parallel Design Strategy

**Key Principle**: Design for MPI parallelism from the start, but implement serial version first for concept verification.

**Implementation Strategy**:
1. **Template-based abstraction**: Use template parameters to switch between serial and parallel MFEM classes
2. **Conditional compilation**: Use `#ifdef MFEM_USE_MPI` for parallel-specific code paths
3. **Serial-first verification**: Complete all unit tests with serial version before parallel implementation
4. **Incremental parallelization**: Parallelize components one at a time with verification at each step

---

## Detailed Component Design

### 6.1 Serial/Parallel Type Abstraction

```cpp
// common/seas_types.hpp

#ifndef SEAS_TYPES_HPP
#define SEAS_TYPES_HPP

#include "mfem.hpp"

namespace mfem {
namespace seas {

#ifdef SEAS_USE_MPI

// Parallel type aliases
using SEASMesh = ParMesh;
using SEASFiniteElementSpace = ParFiniteElementSpace;
using SEASBilinearForm = ParBilinearForm;
using SEASLinearForm = ParLinearForm;
using SEASGridFunction = ParGridFunction;

// Parallel solver types
using SEASSolver = HyprePCG;
using SEASPreconditioner = HypreBoomerAMG;

// MPI communicator wrapper
inline MPI_Comm GetSEASComm() { return MPI_COMM_WORLD; }
inline int GetSEASRank() {
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank); return rank;
}
inline int GetSEASSize() {
    int size; MPI_Comm_size(MPI_COMM_WORLD, &size); return size;
}
inline bool IsSEASRoot() { return GetSEASRank() == 0; }

#else

// Serial type aliases
using SEASMesh = Mesh;
using SEASFiniteElementSpace = FiniteElementSpace;
using SEASBilinearForm = BilinearForm;
using SEASLinearForm = LinearForm;
using SEASGridFunction = GridFunction;

// Serial solver types
using SEASSolver = CGSolver;
using SEASPreconditioner = GSSmoother;

// Serial stubs for MPI functions
inline int GetSEASRank() { return 0; }
inline int GetSEASSize() { return 1; }
inline bool IsSEASRoot() { return true; }

#endif

} // namespace seas
} // namespace mfem

#endif // SEAS_TYPES_HPP
```

### MFEM Serial/Parallel Class Mapping

| Serial Class | Parallel Class | Notes |
|--------------|----------------|-------|
| `Mesh` | `ParMesh` | Distributed mesh with ghost elements |
| `FiniteElementSpace` | `ParFiniteElementSpace` | True DOFs vs local DOFs |
| `BilinearForm` | `ParBilinearForm` | Shared face assembly |
| `LinearForm` | `ParLinearForm` | Parallel vector assembly |
| `GridFunction` | `ParGridFunction` | Distributed solution vector |
| `Vector` | `HypreParVector` | Hypre-compatible parallel vectors |
| `CGSolver` | `HyprePCG` | Parallel Krylov solver |
| `GSSmoother` | `HypreBoomerAMG` | Algebraic multigrid preconditioner |

### 6.2 MPI Context Wrapper

```cpp
// common/mpi_context.hpp

#ifndef MPI_CONTEXT_HPP
#define MPI_CONTEXT_HPP

#include "mfem.hpp"

namespace mfem {
namespace seas {

/// RAII wrapper for MPI initialization (no-op in serial)
class MPIContext {
public:
#ifdef SEAS_USE_MPI
    MPIContext(int *argc, char ***argv) {
        MPI_Init(argc, argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
        MPI_Comm_size(MPI_COMM_WORLD, &size_);
    }
    ~MPIContext() { MPI_Finalize(); }
#else
    MPIContext(int *argc, char ***argv) : rank_(0), size_(1) {}
    ~MPIContext() = default;
#endif

    int Rank() const { return rank_; }
    int Size() const { return size_; }
    bool IsRoot() const { return rank_ == 0; }

    /// Global reduction for maximum value (for CFL condition)
    real_t GlobalMax(real_t local_val) const {
#ifdef SEAS_USE_MPI
        real_t global_val;
        MPI_Allreduce(&local_val, &global_val, 1, MPI_DOUBLE, MPI_MAX,
                      MPI_COMM_WORLD);
        return global_val;
#else
        return local_val;
#endif
    }

    /// Global reduction for minimum value
    real_t GlobalMin(real_t local_val) const {
#ifdef SEAS_USE_MPI
        real_t global_val;
        MPI_Allreduce(&local_val, &global_val, 1, MPI_DOUBLE, MPI_MIN,
                      MPI_COMM_WORLD);
        return global_val;
#else
        return local_val;
#endif
    }

    /// Global sum reduction
    real_t GlobalSum(real_t local_val) const {
#ifdef SEAS_USE_MPI
        real_t global_val;
        MPI_Allreduce(&local_val, &global_val, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        return global_val;
#else
        return local_val;
#endif
    }

    /// Global sum reduction for integers
    int GlobalSumInt(int local_val) const {
#ifdef SEAS_USE_MPI
        int global_val;
        MPI_Allreduce(&local_val, &global_val, 1, MPI_INT, MPI_SUM,
                      MPI_COMM_WORLD);
        return global_val;
#else
        return local_val;
#endif
    }

    /// Barrier synchronization
    void Barrier() const {
#ifdef SEAS_USE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }

private:
    int rank_;
    int size_;
};

} // namespace seas
} // namespace mfem

#endif // MPI_CONTEXT_HPP
```

### 6.3 Parallel Utilities

```cpp
// common/parallel_utils.hpp

#ifndef PARALLEL_UTILS_HPP
#define PARALLEL_UTILS_HPP

#include "mfem.hpp"
#include "mpi_context.hpp"

namespace mfem {
namespace seas {

#ifdef SEAS_USE_MPI

/// Gather vector data from all ranks to root
void GatherVectorToRoot(const Vector &local_data,
                        Vector &global_data,
                        const Array<int> &recv_counts,
                        const Array<int> &displacements,
                        MPI_Comm comm = MPI_COMM_WORLD);

/// Scatter vector data from root to all ranks
void ScatterVectorFromRoot(const Vector &global_data,
                           Vector &local_data,
                           const Array<int> &send_counts,
                           const Array<int> &displacements,
                           MPI_Comm comm = MPI_COMM_WORLD);

/// Broadcast scalar value from root
template <typename T>
void BroadcastFromRoot(T &value, MPI_Comm comm = MPI_COMM_WORLD) {
    MPI_Bcast(&value, sizeof(T), MPI_BYTE, 0, comm);
}

/// Check if all ranks have the same value (for debugging)
template <typename T>
bool AllRanksAgree(T local_value, MPI_Comm comm = MPI_COMM_WORLD) {
    T min_val, max_val;
    MPI_Allreduce(&local_value, &min_val, 1, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(&local_value, &max_val, 1, MPI_DOUBLE, MPI_MAX, comm);
    return (min_val == max_val);
}

#else

// Serial stubs
inline void GatherVectorToRoot(const Vector &local_data,
                               Vector &global_data,
                               const Array<int> &,
                               const Array<int> &,
                               int = 0) {
    global_data = local_data;
}

inline void ScatterVectorFromRoot(const Vector &global_data,
                                  Vector &local_data,
                                  const Array<int> &,
                                  const Array<int> &,
                                  int = 0) {
    local_data = global_data;
}

template <typename T>
void BroadcastFromRoot(T &, int = 0) {}

template <typename T>
bool AllRanksAgree(T, int = 0) { return true; }

#endif

} // namespace seas
} // namespace mfem

#endif // PARALLEL_UTILS_HPP
```

### 6.4 CMake Configuration

```cmake
# CMakeLists.txt

cmake_minimum_required(VERSION 3.10)
project(SEAS LANGUAGES CXX)

# Find MFEM
find_package(MFEM REQUIRED)

# Common source files (shared between serial and parallel)
set(SEAS_COMMON_SOURCES
    common/seas_types.hpp
    common/mpi_context.hpp
    common/parallel_utils.hpp
    config/bp2_params.hpp
    friction/friction_law.hpp
    friction/dieterich_ruina.hpp
    friction/state_evolution.hpp
    domain/domain_operator.hpp
    domain/antiplane_operator.hpp
    fault/fault_geometry.hpp
    fault/rate_state_fault.hpp
    solver/seas_operator.hpp
    solver/time_stepper.hpp
    io/benchmark_output.hpp
    io/probe_output.hpp
)

# Always build serial version
add_executable(seas seas.cpp ${SEAS_COMMON_SOURCES})
target_link_libraries(seas mfem)
target_include_directories(seas PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Build parallel version if MPI is available
if(MFEM_USE_MPI)
    add_executable(pseas pseas.cpp ${SEAS_COMMON_SOURCES})
    target_link_libraries(pseas mfem)
    target_include_directories(pseas PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_compile_definitions(pseas PRIVATE SEAS_USE_MPI)

    # Parallel-specific source files
    target_sources(pseas PRIVATE
        solver/parallel_time_stepper.hpp
        io/parallel_probe_output.hpp
        fault/parallel_fault_data.hpp
    )
endif()

# Unit tests
enable_testing()

# Serial tests
add_executable(test_friction tests/unit/test_friction_law.cpp)
target_link_libraries(test_friction mfem gtest gtest_main)
add_test(NAME friction_law COMMAND test_friction)

add_executable(test_state_evolution tests/unit/test_state_evolution.cpp)
target_link_libraries(test_state_evolution mfem gtest gtest_main)
add_test(NAME state_evolution COMMAND test_state_evolution)

add_executable(test_antiplane tests/unit/test_antiplane.cpp)
target_link_libraries(test_antiplane mfem gtest gtest_main)
add_test(NAME antiplane COMMAND test_antiplane)

# Parallel tests (if MPI available)
if(MFEM_USE_MPI)
    add_executable(test_parallel_domain tests/parallel/test_parallel_domain.cpp)
    target_link_libraries(test_parallel_domain mfem gtest gtest_main)
    target_compile_definitions(test_parallel_domain PRIVATE SEAS_USE_MPI)

    # MPI test runner
    add_test(NAME parallel_domain_np2
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                     $<TARGET_FILE:test_parallel_domain>)
    add_test(NAME parallel_domain_np4
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                     $<TARGET_FILE:test_parallel_domain>)

    add_executable(test_parallel_fault tests/parallel/test_parallel_fault.cpp)
    target_link_libraries(test_parallel_fault mfem gtest gtest_main)
    target_compile_definitions(test_parallel_fault PRIVATE SEAS_USE_MPI)

    add_test(NAME parallel_fault_np4
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                     $<TARGET_FILE:test_parallel_fault>)

    add_executable(test_scaling tests/parallel/test_scaling.cpp)
    target_link_libraries(test_scaling mfem gtest gtest_main)
    target_compile_definitions(test_scaling PRIVATE SEAS_USE_MPI)
endif()
```

---

## Parallel Execution Model

```
+-----------------------------------------------------------------------+
|                         MPI Parallel Execution                         |
+-----------------------------------------------------------------------+
|                                                                        |
|  Rank 0              Rank 1              Rank 2              Rank N    |
|  +-----------+      +-----------+      +-----------+      +-----------+|
|  | ParMesh   |      | ParMesh   |      | ParMesh   |      | ParMesh   ||
|  | (local)   |      | (local)   |      | (local)   |      | (local)   ||
|  | + ghosts  |      | + ghosts  |      | + ghosts  |      | + ghosts  ||
|  +-----+-----+      +-----+-----+      +-----+-----+      +-----+-----+|
|        |                  |                  |                  |      |
|  +-----+-----+      +-----+-----+      +-----+-----+      +-----+-----+|
|  |Local Fault|      |Local Fault|      |Local Fault|      |Local Fault||
|  |DOFs + theta|     |DOFs + theta|     |DOFs + theta|     |DOFs + theta|
|  +-----+-----+      +-----+-----+      +-----+-----+      +-----+-----+|
|        |                  |                  |                  |      |
|        +------------------+------------------+------------------+      |
|                                 |                                      |
|                    +------------+------------+                         |
|                    |   HyprePCG + BoomerAMG  |                         |
|                    |   (Parallel solve)      |                         |
|                    +-------------------------+                         |
|                                                                        |
|  Communication Points:                                                 |
|  1. Domain solve: HYPRE handles parallel linear algebra                |
|  2. Fault state: No communication needed (local computation)           |
|  3. Output: MPI_Gather for probe stations spanning ranks               |
|  4. Time stepping: MPI_Allreduce for global max slip rate (CFL)        |
|                                                                        |
+------------------------------------------------------------------------+
```

---

## Build and Test Commands

```bash
# Serial build
mkdir build-serial && cd build-serial
cmake .. -DMFEM_USE_MPI=OFF
make seas

# Parallel build
mkdir build-parallel && cd build-parallel
cmake .. -DMFEM_USE_MPI=ON
make pseas

# Run serial
./seas -bp bp2 -tf 1e10

# Run parallel (4 ranks)
mpirun -np 4 ./pseas -bp bp2 -tf 1e10

# Run tests
ctest --output-on-failure
```

---

## Unit Tests

### MPI Context Tests

```cpp
// tests/parallel/test_mpi_context.cpp

TEST(MPIContext, Initialization) {
    MPIContext ctx(nullptr, nullptr);

    // Rank should be in [0, size-1]
    EXPECT_GE(ctx.Rank(), 0);
    EXPECT_LT(ctx.Rank(), ctx.Size());

    // Size should be positive
    EXPECT_GT(ctx.Size(), 0);
}

TEST(MPIContext, GlobalMax) {
    MPIContext ctx(nullptr, nullptr);

    // Each rank contributes its rank number
    real_t local_val = ctx.Rank();
    real_t global_max = ctx.GlobalMax(local_val);

    EXPECT_EQ(global_max, ctx.Size() - 1);
}

TEST(MPIContext, GlobalSum) {
    MPIContext ctx(nullptr, nullptr);

    // Each rank contributes 1
    int local_val = 1;
    int global_sum = ctx.GlobalSumInt(local_val);

    EXPECT_EQ(global_sum, ctx.Size());
}
```

---

## File Organization

```
common/
├── seas_types.hpp           # Type aliases for serial/parallel switching
├── mpi_context.hpp          # MPI wrapper (no-op in serial)
└── parallel_utils.hpp       # MPI helper functions
```

---

## Acceptance Criteria

| Test | Pass Criteria |
|------|---------------|
| Serial compile | Code compiles without SEAS_USE_MPI defined |
| Parallel compile | Code compiles with SEAS_USE_MPI defined |
| Type aliases | Same code works with serial and parallel types |
| MPI context | Initialization, rank/size queries work correctly |
| Global reductions | GlobalMax, GlobalMin, GlobalSum return correct values |
