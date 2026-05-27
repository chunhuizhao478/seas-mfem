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

#ifndef MFEM_SEAS_CONFIG_HPP
#define MFEM_SEAS_CONFIG_HPP

#include "mfem.hpp"
#include "../domain/domain_config.hpp"
#include <string>
#include <set>
#include <vector>

namespace mfem
{
namespace seas
{

/// Mesh configuration
struct MeshConfig
{
   std::string file;           ///< Path to mesh file (.msh)
   real_t scale = 1000.0;      ///< Coordinate scale factor (default: km→m)
   int order = 1;              ///< Polynomial order for DG
};

/// Material configuration
struct MaterialConfig
{
   real_t density = 2670.0;    ///< Density rho [kg/m³]
   real_t cs = 3464.0;         ///< S-wave velocity [m/s]
   real_t nu = 0.25;           ///< Poisson's ratio [-]

   real_t mu() const { return density * cs * cs; }
   real_t lambda() const { return 2.0 * nu * mu() / (1.0 - 2.0 * nu); }
};

/// Rate-and-state friction configuration
struct FrictionConfig
{
   real_t V0 = 1.0e-6;        ///< Reference slip rate [m/s]
   real_t f0 = 0.6;           ///< Reference friction coefficient [-]
   real_t b = 0.03;           ///< State evolution parameter [-]
   real_t L0 = 0.14;          ///< Critical slip distance [m]
   real_t L_nuc = 0.13;       ///< Critical slip distance in nucleation [m]
   real_t a0 = 0.004;         ///< Direct effect (VW zone) [-]
   real_t amax = 0.04;        ///< Direct effect (VS zone) [-]
   real_t sigma_n = 25.0e6;   ///< Effective normal stress [Pa]
};

/// Loading configuration
struct LoadingConfig
{
   real_t Vp = 1.0e-9;        ///< Plate rate [m/s]
   real_t V_init = 1.0e-9;    ///< Initial slip rate [m/s]
   real_t V_nuc = 0.01;       ///< Nucleation slip rate [m/s]
   real_t delta_tau_factor = 0.0;  ///< Nucleation pre-stress multiplier
   real_t nucleation_eps = 1.0e-3; ///< Nucleation zone tolerance [m]
   bool smooth_nucleation = false; ///< Gaussian taper for nucleation
};

/// Fault geometry configuration
struct FaultGeomConfig
{
   real_t Wf = 40.0e3;        ///< Fault depth [m]
   real_t lf = 100.0e3;       ///< Fault length [m]
   real_t hs = 2.0e3;         ///< Shallow VS zone width [m]
   real_t ht = 2.0e3;         ///< VW-VS transition width [m]
   real_t H = 12.0e3;         ///< Uniform VW region width [m]
   real_t l_vw = 60.0e3;      ///< VW region length [m]
   real_t w_nuc = 12.0e3;     ///< Nucleation zone width [m]
};

/// Boundary configuration (from TOML)
struct BoundaryTomlConfig
{
   std::set<int> dirichlet_attrs;    ///< User must specify (or use benchmark preset)
   std::set<int> natural_attrs;      ///< User must specify (or use benchmark preset)
   int fault_attr = 0;               ///< 0 = not set (user must specify or use preset)
};

/// Solver configuration
struct SolverConfig
{
   std::string dg_method = "IP";        ///< "IP" or "BR2"
   std::string solver_type = "mumps";   ///< "cg", "mumps", "mumps-blr"
   real_t penalty_factor = 1.0;
   real_t blr_tol = 1e-12;
   int face_basis_type = 1;  // BasisType::GaussLobatto = 1
   bool check_residual = false;
};

/// Time stepping configuration
struct TimeConfig
{
   real_t t_final = 56844000000.0;  ///< Final time [s] (~1800 yr)
   real_t atol = 1e-7;
   real_t rtol = 1e-50;
   int max_steps = 10000000;
   int checkpoint_interval = 5000;
   bool tandem_time_stepping = true;
   bool use_petsc_ts = true;        ///< Use PETSc TS RK45 (default, matches Tandem)
   std::string petsc_ts_options = "tests/verification/petsc_ts_rk45_tandem.cfg";
};

/// Output configuration
struct OutputConfig
{
   std::string output_dir = ".";
   std::string output_prefix = "bp5";
   bool write_every_step = false;
   std::string ref_dir;
   real_t regression_tolerance = -1.0;  ///< Negative = no regression check
};

/// Simulation mode configuration
struct SimulationConfig
{
   std::string mode = "qd";   ///< "qd" (quasi-dynamic), "dynamic", "hybrid"
};

/// SAFS sidecar stress configuration (Phase 6 §7 of PLAN_onfaultstress.md).
///
/// Strict opt-in: `use_sidecar = false` (the default) preserves every
/// existing BP5 / BP2 / TPV102 / TPV205 code path bit-exact.  When set
/// to true, the driver must:
///   1. Instantiate `StressField3D(sidecar_path)`.
///   2. Call `FaultGeometry::ComputeParams(field, P_p_pa, P_p_grad_pa_per_m, min_sigma_n_pa)`.
///   3. Toggle `RateStateFaultOperator::SetSAFSMode(true, ...)` with
///      pointers to FaultGeometry's per-DOF pre-stress / sigma_n.
struct StressConfig
{
   bool use_sidecar = false;            ///< false (default) → BP5 path
   std::string sidecar_path;            ///< Path to stress_safs.h5 (required if use_sidecar)
   real_t P_p_pa = 0.0;                 ///< Constant pore-pressure offset [Pa]
   real_t P_p_grad_pa_per_m = 0.0;      ///< Depth gradient of pore pressure [Pa/m]
   real_t min_sigma_n_pa = 0.0;         ///< Optional Pa-valued floor on σ_n; 0 = no clamp
};

/// @brief Top-level SEAS simulation configuration.
///
/// Pure data struct — no factory methods, no logic.
/// Populated from TOML file by SEASConfigParser.
struct SEASConfig
{
   std::string benchmark;     ///< Optional preset: "bp5", "bp1", "bp2"
   MeshConfig mesh;
   MaterialConfig material;
   FrictionConfig friction;
   LoadingConfig loading;
   FaultGeomConfig fault_geom;
   BoundaryTomlConfig boundary;
   SolverConfig solver;
   TimeConfig time;
   OutputConfig output;
   SimulationConfig simulation;
   StressConfig stress;       ///< Phase 6 §7: SAFS sidecar (opt-in)
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CONFIG_HPP
