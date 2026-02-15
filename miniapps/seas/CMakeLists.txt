# Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
# at the Lawrence Livermore National Laboratory. All Rights reserved. See files
# LICENSE and NOTICE for details. LLNL-CODE-806117.
#
# This file is part of the MFEM library. For more information and source code
# availability visit https://mfem.org.
#
# MFEM is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.

# SEAS (Sequences of Earthquakes and Aseismic Slip) miniapp
# Implements SCEC SEAS benchmark problems (BP1, BP2, etc.)

# =============================================================================
# Header files - Phase 1: Friction and State Evolution
# =============================================================================

list(APPEND SEAS_HEADERS
    friction/friction_law.hpp
    friction/state_evolution.hpp
    friction/dieterich_ruina.hpp
    config/bp2_params.hpp
    common/seas_types.hpp
    domain/domain_operator.hpp
    domain/bp2_mesh.hpp
    domain/antiplane_operator.hpp
    integrator/dg_br2_integrator.hpp)

# Phase 3: Fault operators
list(APPEND SEAS_HEADERS
    fault/fault_geometry.hpp
    fault/rate_state_fault.hpp)

# Phase 4: SEAS operators and time stepping
# list(APPEND SEAS_HEADERS
#     solver/seas_operator.hpp
#     solver/quasi_dynamic.hpp
#     solver/time_stepper.hpp)

# Phase 5: I/O
# list(APPEND SEAS_HEADERS
#     io/benchmark_output.hpp
#     io/probe_output.hpp)

# Phase 6: Parallel infrastructure
list(APPEND SEAS_HEADERS
    common/mpi_context.hpp
    common/parallel_utils.hpp)

convert_filenames_to_full_paths(SEAS_HEADERS)

set(SEAS_COMMON_FILES EXTRA_HEADERS ${SEAS_HEADERS})

# =============================================================================
# Phase 1 Unit Tests
# =============================================================================

# Friction law unit tests
add_mfem_miniapp(seas_test_friction_law
  MAIN tests/unit/test_friction_law.cpp
  ${SEAS_COMMON_FILES}
  LIBRARIES mfem)

# State evolution unit tests
add_mfem_miniapp(seas_test_state_evolution
  MAIN tests/unit/test_state_evolution.cpp
  ${SEAS_COMMON_FILES}
  LIBRARIES mfem)

# Antiplane domain operator unit tests (Phase 2)
add_mfem_miniapp(seas_test_antiplane
  MAIN tests/unit/test_antiplane.cpp
  ${SEAS_COMMON_FILES}
  LIBRARIES mfem)

# Fault operator unit tests (Phase 3)
add_mfem_miniapp(seas_test_fault_operator
  MAIN tests/unit/test_fault_operator.cpp
  ${SEAS_COMMON_FILES}
  LIBRARIES mfem)

# Add tests if testing is enabled
if (MFEM_ENABLE_TESTING)
  add_test(NAME seas_test_friction_law
    COMMAND $<TARGET_FILE:seas_test_friction_law>)
  add_test(NAME seas_test_state_evolution
    COMMAND $<TARGET_FILE:seas_test_state_evolution>)
  add_test(NAME seas_test_antiplane
    COMMAND $<TARGET_FILE:seas_test_antiplane>)
  add_test(NAME seas_test_fault_operator
    COMMAND $<TARGET_FILE:seas_test_fault_operator>)
endif()

# =============================================================================
# Phase 6: Parallel Unit Tests (requires MPI)
# =============================================================================

if (MFEM_USE_MPI)
  add_mfem_miniapp(seas_test_mpi_context
    MAIN tests/parallel/test_mpi_context.cpp
    ${SEAS_COMMON_FILES}
    LIBRARIES mfem)
  target_compile_definitions(seas_test_mpi_context PRIVATE SEAS_USE_MPI)

  add_mfem_miniapp(seas_test_parallel_utils
    MAIN tests/parallel/test_parallel_utils.cpp
    ${SEAS_COMMON_FILES}
    LIBRARIES mfem)
  target_compile_definitions(seas_test_parallel_utils PRIVATE SEAS_USE_MPI)

  # Phase 7: Parallel domain operator tests
  add_mfem_miniapp(seas_test_parallel_domain
    MAIN tests/parallel/test_parallel_domain.cpp
    ${SEAS_COMMON_FILES}
    LIBRARIES mfem)
  target_compile_definitions(seas_test_parallel_domain PRIVATE SEAS_USE_MPI)

  add_mfem_miniapp(seas_mms_antiplane_parallel
    MAIN tests/parallel/mms_antiplane_parallel.cpp
    ${SEAS_COMMON_FILES}
    LIBRARIES mfem)
  target_compile_definitions(seas_mms_antiplane_parallel PRIVATE SEAS_USE_MPI)

  if (MFEM_ENABLE_TESTING)
    add_test(NAME seas_test_mpi_context_np2
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
              $<TARGET_FILE:seas_test_mpi_context>)
    add_test(NAME seas_test_mpi_context_np4
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
              $<TARGET_FILE:seas_test_mpi_context>)
    add_test(NAME seas_test_parallel_utils_np2
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
              $<TARGET_FILE:seas_test_parallel_utils>)
    add_test(NAME seas_test_parallel_utils_np4
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
              $<TARGET_FILE:seas_test_parallel_utils>)
    add_test(NAME seas_test_parallel_domain_np2
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
              $<TARGET_FILE:seas_test_parallel_domain>)
    add_test(NAME seas_test_parallel_domain_np4
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
              $<TARGET_FILE:seas_test_parallel_domain>)
    add_test(NAME seas_mms_antiplane_parallel_np2
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
              $<TARGET_FILE:seas_mms_antiplane_parallel>)
  endif()
endif()

# =============================================================================
# Future: Main driver (uncomment when implemented)
# =============================================================================

# Serial driver
# add_mfem_miniapp(seas
#   MAIN seas.cpp
#   ${SEAS_COMMON_FILES}
#   LIBRARIES mfem)

# Parallel driver (if MPI is enabled)
# if (MFEM_USE_MPI)
#   add_mfem_miniapp(pseas
#     MAIN pseas.cpp
#     ${SEAS_COMMON_FILES}
#     LIBRARIES mfem)
# endif()

# =============================================================================
# Future: BP2 example (uncomment when implemented)
# =============================================================================

# add_mfem_miniapp(seas_bp2
#   MAIN examples/bp2_qd.cpp
#   ${SEAS_COMMON_FILES}
#   LIBRARIES mfem)
