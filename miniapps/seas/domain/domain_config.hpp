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

#include <memory>

namespace mfem
{
namespace seas
{

// fault/fault_frame_field.hpp (MPI builds).
class FaultFrameField;
struct FaultFrameConfig;

/// Where an ElasticityDomainOperator takes its fault orientation from
/// (2D in-plane design rev 3, phase A: one orientation system for the
/// quasi-dynamic and the dynamic operators).
///
/// Both empty (the default): the operator's historical constants, reference
/// normal (0, -1, 0) and up (0, 0, 1), in the fault frames, the
/// prescribed-slip sign and the BR2 slip sign (serial meshes, legacy drivers).
/// `field`: read the orientation from this field; the operator checks that
/// its own fault and interior Dirichlet face lists are the field's.
/// `config`: build the field from the operator's own face lists (ParMesh
/// only); FrameField() then hands it to the other operators of the run.
struct FaultFrameSource
{
   std::shared_ptr<const FaultFrameField> field;
   std::shared_ptr<const FaultFrameConfig> config;
   bool Empty() const { return !field && !config; }
};

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
   FaultFrameSource fault_frames;  ///< fault orientation (default: the constants)
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_DOMAIN_CONFIG_HPP
