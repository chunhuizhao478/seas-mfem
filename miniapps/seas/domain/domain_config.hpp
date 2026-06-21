// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef MFEM_SEAS_DOMAIN_CONFIG_HPP
#define MFEM_SEAS_DOMAIN_CONFIG_HPP

#include "mfem.hpp"

namespace mfem
{
namespace seas
{

// Forward-declare enums (defined in antiplane_operator.hpp and
// elasticity_operator.hpp respectively). We re-declare them here
// so domain_config.hpp can be included without pulling in the full
// operator headers. The actual enum definitions remain authoritative.
//
// NOTE: These are NOT redefinitions — they reference the existing enums.
// DGMethod is in antiplane_operator.hpp, SolverType in elasticity_operator.hpp.

/// @brief Discretization and solver configuration for the domain operator.
///
/// Bundles all non-physics parameters that were previously passed as
/// constructor arguments or set via post-construction setters.
/// Default values match the current production configuration.
struct DomainConfig
{
   int face_basis_type = BasisType::GaussLobatto;
   real_t penalty_factor = 1.0;
   real_t blr_tol = 1e-12;
   bool check_residual = false;
   bool match_quad_order = false;
   real_t cfl_factor = 0.5;  ///< CFL stability factor for explicit time stepping

   // Krylov/AMG knobs for the iterative solver paths (CG_AMG / GMRES_AMG),
   // Phase 4.  Defaults reproduce the pre-Phase-4 hardcoded dispatch values so
   // existing callers are unchanged; the quasi-dynamic driver overrides them
   // from the [solver] config.  (MUMPS/SuperLU/STRUMPACK ignore these.)
   real_t ksp_rtol = 1e-10;        ///< Krylov relative tolerance
   real_t ksp_atol = 0.0;          ///< Krylov absolute tolerance
   int    ksp_maxit = 10000;       ///< Krylov max iterations
   int    amg_print_level = 0;     ///< AMG / Krylov print verbosity
   /// CG_AMG only: call HypreBoomerAMG::SetElasticityOptions (CFEM rigid-body
   /// near-null-space).  Default true reproduces the pre-existing dispatch, but
   /// it is UNPROVEN on DG (R-006) and bloats the AMG hierarchy there (slow
   /// V-cycles); set false to use SetSystemsOptions only (the GMRES_AMG path,
   /// which the iterative-vs-direct test verifies converges).
   bool   amg_elasticity_options = true;
   /// BoomerAMG tuning (CG_AMG/GMRES_AMG).  relax 8 (l1-symmetric-GS) is fully
   /// parallel (the HYPRE default hybrid-GS stalls at np>1).  Aggressive
   /// coarsening is OFF by default: on the ill-conditioned DG-elasticity operator
   /// it makes CG STAGNATE (~1e-7, below ksp_rtol) on a loaded solve -> a
   /// near-null-space-contaminated solution -> garbage traction -> NaN.
   int    amg_relax_type = 8;          ///< HYPRE relax type (8/16/3/18); -1 = HYPRE default
   int    amg_aggressive_levels = 0;   ///< aggressive-coarsening levels (0 = off, DEFAULT)
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_DOMAIN_CONFIG_HPP
