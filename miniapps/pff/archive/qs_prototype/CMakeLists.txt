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

# Phase Field Fracture (PFF) miniapp
# Implements phase field method for brittle fracture simulation

# Header files - Phase 1: Materials
list(APPEND PFF_HEADERS
    materials/degradation_function.hpp
    materials/spectral_decomposition.hpp
    materials/pff_material.hpp)

# Header files - Phase 2: Integrators
list(APPEND PFF_HEADERS
    integrators/degraded_elasticity_integrator.hpp
    integrators/damage_diffusion_integrator.hpp
    integrators/damage_source_integrator.hpp
    operators/elasticity_operator.hpp
    operators/damage_operator.hpp
    pff_solver.hpp)

convert_filenames_to_full_paths(PFF_HEADERS)

set(PFF_COMMON_FILES EXTRA_HEADERS ${PFF_HEADERS})

# =============================================================================
# Unit Tests
# =============================================================================

# Phase 1 unit tests (materials)
add_mfem_miniapp(pff_unit_tests_phase1
  MAIN tests/unit_tests_phase1.cpp
  ${PFF_COMMON_FILES}
  LIBRARIES mfem)

# Phase 2 unit tests (integrators)
add_mfem_miniapp(pff_unit_tests_phase2
  MAIN tests/unit_tests_phase2.cpp
  ${PFF_COMMON_FILES}
  LIBRARIES mfem)

# Phase 3 unit tests (operators)
add_mfem_miniapp(pff_unit_tests_phase3
  MAIN tests/unit_tests_phase3.cpp
  ${PFF_COMMON_FILES}
  LIBRARIES mfem)

# Phase 4 unit tests (staggered solver)
add_mfem_miniapp(pff_unit_tests_phase4
  MAIN tests/unit_tests_phase4.cpp
  ${PFF_COMMON_FILES}
  LIBRARIES mfem)

# Add tests if testing is enabled
if (MFEM_ENABLE_TESTING)
  add_test(NAME pff_unit_tests_phase1
    COMMAND $<TARGET_FILE:pff_unit_tests_phase1>)
  add_test(NAME pff_unit_tests_phase2
    COMMAND $<TARGET_FILE:pff_unit_tests_phase2>)
  add_test(NAME pff_unit_tests_phase3
    COMMAND $<TARGET_FILE:pff_unit_tests_phase3>)
  add_test(NAME pff_unit_tests_phase4
    COMMAND $<TARGET_FILE:pff_unit_tests_phase4>)
endif()

# =============================================================================
# Phase 5: Main driver
# =============================================================================

add_mfem_miniapp(pff
  MAIN pff.cpp
  ${PFF_COMMON_FILES}
  LIBRARIES mfem)

# =============================================================================
# Examples
# =============================================================================

# Mode 2 Pure Shear with Gmsh mesh
add_mfem_miniapp(pff_mode2_gmsh
  MAIN examples/mode2_pure_shear/code/pff.cpp
  ${PFF_COMMON_FILES}
  LIBRARIES mfem)
