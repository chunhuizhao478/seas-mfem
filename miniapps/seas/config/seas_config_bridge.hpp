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

#ifndef MFEM_SEAS_CONFIG_BRIDGE_HPP
#define MFEM_SEAS_CONFIG_BRIDGE_HPP

#include "seas_config.hpp"
#include "bp5_params.hpp"
#include "../domain/boundary_config.hpp"
#include "../domain/domain_config.hpp"
#include "../domain/elasticity_operator.hpp"

namespace mfem
{
namespace seas
{

/// Build BP5Params from SEASConfig fields.
/// Fields not in SEASConfig (V_zero, seconds_per_year) keep BP5Params defaults.
inline BP5Params BuildBP5Params(const SEASConfig &config)
{
   BP5Params p;
   p.rho = config.material.density;
   p.cs = config.material.cs;
   p.nu = config.material.nu;
   p.V0 = config.friction.V0;
   p.f0 = config.friction.f0;
   p.b = config.friction.b;
   p.L0 = config.friction.L0;
   p.L_nuc = config.friction.L_nuc;
   p.a0 = config.friction.a0;
   p.amax = config.friction.amax;
   p.sigma_n = config.friction.sigma_n;
   p.Vp = config.loading.Vp;
   p.V_init = config.loading.V_init;
   p.V_nuc = config.loading.V_nuc;
   p.delta_tau_factor = config.loading.delta_tau_factor;
   p.nucleation_eps = config.loading.nucleation_eps;
   p.smooth_nucleation = config.loading.smooth_nucleation;
   p.Wf = config.fault_geom.Wf;
   p.lf = config.fault_geom.lf;
   p.hs = config.fault_geom.hs;
   p.ht = config.fault_geom.ht;
   p.H = config.fault_geom.H;
   p.l_vw = config.fault_geom.l_vw;
   p.w_nuc = config.fault_geom.w_nuc;
   p.t_final = config.time.t_final;
   return p;
}

/// Parse DG method string to enum.
/// Accepts only values that SEASConfigParser::Validate() allows.
inline DGMethod ParseDGMethod(const std::string &s)
{
   if (s == "IP" || s == "ip") { return DGMethod::IP; }
   if (s == "BR2" || s == "br2") { return DGMethod::BR2; }
   MFEM_ABORT("Unknown DG method: " << s);
   return DGMethod::IP;
}

/// Parse solver type string to enum.
inline SolverType ParseSolverType(const std::string &s)
{
   if (s == "cg" || s == "CG") { return SolverType::CG_AMG; }
   if (s == "mumps" || s == "MUMPS") { return SolverType::MUMPS; }
   if (s == "mumps-blr" || s == "MUMPS-BLR") { return SolverType::MUMPS_BLR; }
   if (s == "superlu") { return SolverType::SUPERLU; }
   if (s == "strumpack") { return SolverType::STRUMPACK; }
   if (s == "gmres") { return SolverType::GMRES_AMG; }
   MFEM_ABORT("Unknown solver type: " << s);
   return SolverType::CG_AMG;
}

/// Build BoundaryConfig from TOML boundary section + plate rate.
inline BoundaryConfig BuildBoundaryConfig(const BoundaryTomlConfig &bdr,
                                          real_t Vp)
{
   BoundaryConfig bc;
   bc.dirichlet_attrs = bdr.dirichlet_attrs;
   bc.natural_attrs = bdr.natural_attrs;
   bc.fault_attr = bdr.fault_attr;
   bc.default_dirichlet_func = MakeBP5DirichletFunc(Vp);
   return bc;
}

/// Build DomainConfig from solver config section.
inline DomainConfig BuildDomainConfig(const SolverConfig &solver)
{
   DomainConfig dc;
   dc.face_basis_type = solver.face_basis_type;
   dc.penalty_factor = solver.penalty_factor;
   dc.blr_tol = solver.blr_tol;
   dc.check_residual = solver.check_residual;
   return dc;
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CONFIG_BRIDGE_HPP
