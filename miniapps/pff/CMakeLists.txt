# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC. Produced
# at the Lawrence Livermore National Laboratory. All Rights reserved. See files
# LICENSE and NOTICE for details. LLNL-CODE-806117.
#
# This file is part of the MFEM library. For more information and source code
# availability visit https://mfem.org.
#
# MFEM is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.

# pff: dynamic phase-field fracture (port of farms@cdms spe_hftc and
# cmame_revision). The primary build is the Makefile (`make test`); this file
# builds the same unit tests under CMake. The quasi-static prototype in
# archive/qs_prototype is not built.

set(PFF_UNIT_TESTS
    test_phase_field_functions
    test_strain_split
    test_poro_damage
    test_pulse_load
    test_config
    test_explicit_operator
    test_damage_step)

foreach(test ${PFF_UNIT_TESTS})
  add_mfem_miniapp(pff_${test}
    MAIN tests/unit/${test}.cpp
    LIBRARIES mfem)
  target_include_directories(pff_${test} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/extern/toml11)
  if (MFEM_ENABLE_TESTING)
    # test_config reads <model>/configs/ relative to the miniapp directory.
    add_test(NAME pff_${test} COMMAND $<TARGET_FILE:pff_${test}>
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  endif()
endforeach()

# Drivers
add_mfem_miniapp(pff_qs_driver
  MAIN drivers/qs_driver.cpp
  LIBRARIES mfem)
target_include_directories(pff_qs_driver PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/toml11)

add_mfem_miniapp(pff_driver
  MAIN drivers/pff_driver.cpp
  LIBRARIES mfem)
target_include_directories(pff_driver PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/toml11)
