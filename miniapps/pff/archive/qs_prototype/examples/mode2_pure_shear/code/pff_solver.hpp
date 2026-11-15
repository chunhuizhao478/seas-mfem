// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// Wrapper header for pff_solver.hpp that adjusts include paths for the example
// directory structure.
//
// This example is located in: miniapps/pff/examples/mode2_pure_shear/code/
// The main pff sources are in: miniapps/pff/

#ifndef MFEM_PFF_EXAMPLES_MODE2_SOLVER_HPP
#define MFEM_PFF_EXAMPLES_MODE2_SOLVER_HPP

// Include the main pff_solver.hpp from the parent pff directory
// The relative path from examples/mode2_pure_shear/code/ to pff/ is ../../../
#include "../../../pff_solver.hpp"

#endif // MFEM_PFF_EXAMPLES_MODE2_SOLVER_HPP
