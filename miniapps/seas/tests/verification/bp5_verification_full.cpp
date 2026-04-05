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

// BP5 Full Simulation: 3D Earthquake Cycles
//
// Runs a parallel BP5-QD simulation with DG elasticity on a 3D hex mesh.
// Outputs SCEC-format time series at 10 on-fault stations with
// checkpoint/restart.
//
// Usage:
//   mpirun -np N ./seas_bp5_full [options]
//
// Options:
//   --mesh FILE                Load mesh from Gmsh .msh file
//   --mesh-scale S             Coordinate scale factor (default: 1000 = km->m)
//   --inline-mesh              Use inline Cartesian mesh (for smoke tests)
//   --nx/--ny/--nz N           Inline mesh element counts per half-domain
//   --Lx/--Ly/--Lz L           Inline mesh domain half-sizes [m]
//   --output-dir DIR           Output directory (default: .)
//   --output-prefix PFX        Output prefix (default: "bp5_full")
//   --tfinal T                 Override final time in seconds
//   --checkpoint-interval N    Write checkpoint every N steps (default: 5000)
//   --restart PREFIX           Restart from checkpoint
//   --ref-dir DIR              Reference data directory
//   --comparison-only          Skip simulation, compare only
//   --write-every-step         Write output at every accepted step
//   --V-nuc VAL                Nucleation slip rate [m/s] (default: 0.03)
//   --delta-tau-factor VAL     Delta-tau multiplier (default: 1.0, Tandem: 0.0)
//   --solver mumps|cg           Linear solver (default: cg)
//   --check-residual           Warn if post-solve ||K*x-b||/||b|| > 1e-8
//   --monitor-traction N       Log tau_pre/traction/total every N RHS evals
//   --dump-bdr-vtk             Output boundary attributes to VTK
//   --psi-clamp                Enable legacy post-step psi clamping (diagnostic)
//   --diag-psi-clamp           Print/summary diagnostics for psi clamp activity
//   --tandem-time-stepping     Use Tandem-style startup/acceptance defaults
//   --tandem-dt-init DT        Initial dt [s] for Tandem-style startup
//   --petsc-ts                 Use PETSc TS RK45 path (exact Tandem framework)
//   --petsc-ts-options FILE    PETSc options file (default: built-in Tandem rk45)

#include "mfem.hpp"
#ifdef MFEM_USE_PETSC
#include "petsc.h"
#if PETSC_VERSION_LT(3,19,0)
#define PETSC_SUCCESS 0
#endif
#endif
#include "../../solver/seas_operator.hpp"
#include "../../solver/time_stepper.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../fault/fault_geometry.hpp"
#include "../../fault/rate_state_fault.hpp"
#include "../../friction/dieterich_ruina.hpp"
#include "../../friction/state_evolution.hpp"
#include "../../config/bp5_params.hpp"
#include "../../io/bp5_parallel_output.hpp"
#include "../../io/probe_output.hpp"
#include "../../io/checkpoint.hpp"
#include "../../common/mpi_context.hpp"
#include "../../trace/face_trace_logger.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <memory>
#include <vector>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <set>
#include <map>

using namespace mfem;
using namespace mfem::seas;

// ============================================================================
// Inline mesh creation for smoke tests
// ============================================================================

/// Create a 3D hex mesh for BP5 with Tandem boundary attributes.
/// Domain: [-Lx,Lx] x [-Ly,Ly] x [-Lz,0]
/// Boundary attributes (Tandem tags):
///   1 = Natural (z=0 top, z=-Lz bottom)
///   5 = Dirichlet (x=±Lx, y=±Ly far-field)
std::unique_ptr<Mesh> CreateBP5InlineMesh(
   int nx, int ny, int nz,
   real_t Lx, real_t Ly, real_t Lz)
{
   auto mesh = std::make_unique<Mesh>(
      Mesh::MakeCartesian3D(2 * nx, 2 * ny, nz,
                            Element::HEXAHEDRON,
                            2.0 * Lx, 2.0 * Ly, Lz));

   for (int i = 0; i < mesh->GetNV(); i++)
   {
      real_t *v = mesh->GetVertex(i);
      v[0] -= Lx;
      v[1] -= Ly;
      v[2] -= Lz;  // Z ranges [-Lz, 0]
   }

   const real_t tol = 1e-6 * std::max({Lx, Ly, Lz});

   for (int i = 0; i < mesh->GetNBE(); i++)
   {
      Array<int> vertices;
      mesh->GetBdrElementVertices(i, vertices);

      real_t cx = 0.0, cy = 0.0, cz = 0.0;
      for (int j = 0; j < vertices.Size(); j++)
      {
         const real_t *v = mesh->GetVertex(vertices[j]);
         cx += v[0]; cy += v[1]; cz += v[2];
      }
      cx /= vertices.Size();
      cy /= vertices.Size();
      cz /= vertices.Size();

      int attr;
      if (std::abs(cz) < tol || std::abs(cz + Lz) < tol)
      {
         attr = 1;  // Natural (top z=0, bottom z=-Lz)
      }
      else
      {
         attr = 5;  // Dirichlet (far-field x=±Lx, y=±Ly)
      }

      mesh->SetBdrAttribute(i, attr);
   }

   mesh->SetAttributes();
   return mesh;
}

// ============================================================================
// Reference data loading and comparison (adapted for 8-column BP5 format)
// ============================================================================

struct BP5TimeSeriesData
{
   std::vector<double> time;
   std::vector<double> slip_strike;
   std::vector<double> slip_dip;
   std::vector<double> log10_V_strike;
   std::vector<double> log10_V_dip;
   std::vector<double> tau_strike;
   std::vector<double> tau_dip;
   std::vector<double> log10_state;

   int Size() const { return static_cast<int>(time.size()); }
};

/// Load BP5 SCEC-format time series file (8 columns, # comment lines).
bool LoadBP5TimeSeriesFile(const std::string &filename,
                            BP5TimeSeriesData &data,
                            double t_max = 1e30)
{
   std::ifstream file(filename);
   if (!file.is_open()) { return false; }

   std::string line;
   while (std::getline(file, line))
   {
      if (line.empty() || line[0] == '#') { continue; }
      if (line.find("slip") != std::string::npos) { continue; }

      std::istringstream iss(line);
      double t, ss, sd, vs, vd, ts, td, st;
      if (!(iss >> t >> ss >> sd >> vs >> vd >> ts >> td >> st)) { continue; }
      if (t > t_max) { break; }

      data.time.push_back(t);
      data.slip_strike.push_back(ss);
      data.slip_dip.push_back(sd);
      data.log10_V_strike.push_back(vs);
      data.log10_V_dip.push_back(vd);
      data.tau_strike.push_back(ts);
      data.tau_dip.push_back(td);
      data.log10_state.push_back(st);
   }
   return data.Size() > 0;
}

/// Linear interpolation of y_ref at x_ref onto x_target grid.
std::vector<double> InterpolateOnto(const std::vector<double> &x_ref,
                                     const std::vector<double> &y_ref,
                                     const std::vector<double> &x_target)
{
   std::vector<double> y_interp(x_target.size());
   int n_ref = static_cast<int>(x_ref.size());
   int j = 0;

   for (size_t i = 0; i < x_target.size(); i++)
   {
      double xt = x_target[i];
      if (xt <= x_ref.front()) { y_interp[i] = y_ref.front(); continue; }
      if (xt >= x_ref.back())  { y_interp[i] = y_ref.back(); continue; }

      while (j < n_ref - 2 && x_ref[j + 1] < xt) { j++; }
      double frac = (xt - x_ref[j]) / (x_ref[j + 1] - x_ref[j]);
      y_interp[i] = y_ref[j] + frac * (y_ref[j + 1] - y_ref[j]);
   }
   return y_interp;
}

/// Compute relative L2 error: ||sim - ref||_2 / ||ref||_2
double RelativeL2Error(const std::vector<double> &sim,
                       const std::vector<double> &ref)
{
   if (sim.size() != ref.size() || sim.empty()) { return 1e30; }

   double num = 0.0, den = 0.0;
   for (size_t i = 0; i < sim.size(); i++)
   {
      double diff = sim[i] - ref[i];
      num += diff * diff;
      den += ref[i] * ref[i];
   }
   if (den < 1e-30) { return (num < 1e-30) ? 0.0 : 1e30; }
   return std::sqrt(num / den);
}

/// Run comparison for all stations (placeholder — needs reference data).
bool RunComparison(const std::string &output_dir,
                   const std::string &output_prefix,
                   const std::string &ref_dir,
                   const std::vector<Probe2DInterpolator::Station> &stations,
                   double t_final_s)
{
   std::cout << "\n" << std::string(80, '=') << "\n";
   std::cout << "BP5 Full Simulation Verification: Comparison Summary\n";
   std::cout << std::string(80, '=') << "\n\n";

   int num_loaded = 0;
   for (const auto &st : stations)
   {
      std::string sim_file = output_dir + "/" + output_prefix + "_"
                             + st.name + ".txt";
      BP5TimeSeriesData sim_data;
      bool loaded = LoadBP5TimeSeriesFile(sim_file, sim_data, t_final_s);

      std::cout << "  Station " << std::setw(24) << st.name << ": ";
      if (loaded)
      {
         std::cout << sim_data.Size() << " data points";
         num_loaded++;
      }
      else
      {
         std::cout << "MISSING";
      }
      std::cout << "\n";
   }

   // Check for reference data
   for (const auto &st : stations)
   {
      std::string ref_file = ref_dir + "/bp5-qd-" + st.name + ".txt";
      BP5TimeSeriesData ref_data;
      bool ref_loaded = LoadBP5TimeSeriesFile(ref_file, ref_data, t_final_s);
      if (ref_loaded)
      {
         std::cout << "  Reference data found for " << st.name
                   << ": " << ref_data.Size() << " points\n";
      }
   }

   std::cout << "\n=== Verification Summary ===\n";
   std::cout << "  Stations with output: " << num_loaded << " / "
             << stations.size() << "\n";
   std::cout << std::string(80, '=') << "\n";

   return num_loaded > 0;
}

// ============================================================================
// Main driver
// ============================================================================

int main(int argc, char *argv[])
{
   MPIContext mpi(&argc, &argv);

   // =========================================================================
   // Parse command-line arguments
   // =========================================================================
   std::string mesh_file;
   real_t mesh_scale = 1000.0;
   bool inline_mesh = false;
   int nx = 2, ny = 2, nz = 1;
   real_t Lx = 200e3, Ly = 100e3, Lz = 100e3;

   std::string output_dir = ".";
   std::string output_prefix = "bp5_full";
   std::string ref_dir = "bp5/benchmark_data";
   bool comparison_only = false;
   double tfinal_override = -1.0;
   int checkpoint_interval = 5000;
   std::string restart_prefix;
   bool write_every_step = false;
   std::string solver_str = "mumps-blr";
   bool check_residual = false;
   int monitor_traction = 0;
   std::string dg_method_str = "BR2";
   double V_nuc_override = 0.0;
   double delta_tau_factor_override = -1.0;
   double nucleation_eps_override = -1.0;
   bool dump_bdr_vtk = false;
   std::string bc_mode_str = "far-field";
   std::string psi_init_mode_str = "tandem";  // "tandem" (default) or "scec"
   int order = 1;
   double blr_tol = 1e-10;  // MUMPS-BLR tolerance (default 1e-10)

   // v49 Phase 1 flags
   bool smooth_nucleation = false;
   bool match_quad_order = false;
   bool zero_dip_traction = false;   // v51: zero tau_dip after ComputeTraction
   // v54: elastic sigma_n ON by default (matches Tandem DieterichRuinaBase.h:87)
   bool elastic_sigma_n = true;
   // v49 Phase 2: CFL-aware dt and V guard
   real_t dt_init_override = -1.0;  // Manual dt_init override (negative = auto)
   real_t v_guard_factor = -1.0;    // V guard threshold factor (negative = use default 100)
   bool no_v_guard = false;         // v50: disable V guard (for testing only)
   bool tandem_time_stepping = false; // Use Tandem-style startup/acceptance policy
   real_t tandem_dt_init = 0.01;      // Default Tandem-style startup dt [s]
   // v50a: penalty scaling factor (1.0 = default, <1.0 = reduced penalty)
   real_t penalty_factor = 1.0;
   bool no_psi_clamp = true;           // Default OFF: match Tandem (no post-step psi clamp)
   bool dt_init_explicit = false;
   bool v_guard_explicit = false;
   bool psi_clamp_explicit = false;
   bool use_petsc_ts = false;          // Exact Tandem framework: PETSc TS
   std::string petsc_ts_options_file;  // Optional PETSc options file
   bool petsc_initialized = false;
   // v50g: face DOF node type (GaussLobatto has cond(M)=2901 at p=4, ClosedUniform=58)
   int face_basis_type = BasisType::GaussLobatto;
   std::string face_basis_str = "GaussLobatto";

   for (int i = 1; i < argc; i++)
   {
      std::string arg(argv[i]);
      if (arg == "--mesh" && i + 1 < argc) { mesh_file = argv[++i]; }
      if (arg == "--mesh-scale" && i + 1 < argc)
      {
         mesh_scale = std::atof(argv[++i]);
      }
      if (arg == "--inline-mesh") { inline_mesh = true; }
      if (arg == "--nx" && i + 1 < argc) { nx = std::atoi(argv[++i]); }
      if (arg == "--ny" && i + 1 < argc) { ny = std::atoi(argv[++i]); }
      if (arg == "--nz" && i + 1 < argc) { nz = std::atoi(argv[++i]); }
      if (arg == "--Lx" && i + 1 < argc) { Lx = std::atof(argv[++i]); }
      if (arg == "--Ly" && i + 1 < argc) { Ly = std::atof(argv[++i]); }
      if (arg == "--Lz" && i + 1 < argc) { Lz = std::atof(argv[++i]); }
      if (arg == "--output-dir" && i + 1 < argc)
      {
         output_dir = argv[++i];
      }
      if (arg == "--output-prefix" && i + 1 < argc)
      {
         output_prefix = argv[++i];
      }
      if (arg == "--ref-dir" && i + 1 < argc) { ref_dir = argv[++i]; }
      if (arg == "--comparison-only") { comparison_only = true; }
      if (arg == "--tfinal" && i + 1 < argc)
      {
         tfinal_override = std::atof(argv[++i]);
      }
      if (arg == "--checkpoint-interval" && i + 1 < argc)
      {
         checkpoint_interval = std::atoi(argv[++i]);
      }
      if (arg == "--restart" && i + 1 < argc)
      {
         restart_prefix = argv[++i];
      }
      if (arg == "--write-every-step") { write_every_step = true; }
      if (arg == "--mumps") { solver_str = "mumps"; }
      if (arg == "--solver" && i + 1 < argc) { solver_str = argv[++i]; }
      if (arg == "--check-residual") { check_residual = true; }
      if (arg == "--monitor-traction" && i + 1 < argc)
      {
         monitor_traction = std::atoi(argv[++i]);
      }
      if (arg == "--dg-method" && i + 1 < argc) { dg_method_str = argv[++i]; }
      if (arg == "--V-nuc" && i + 1 < argc)
      {
         V_nuc_override = std::atof(argv[++i]);
      }
      if (arg == "--delta-tau-factor" && i + 1 < argc)
      {
         delta_tau_factor_override = std::atof(argv[++i]);
      }
      if (arg == "--nucleation-eps" && i + 1 < argc)
      {
         nucleation_eps_override = std::atof(argv[++i]);
      }
      if (arg == "--dump-bdr-vtk") { dump_bdr_vtk = true; }
      if (arg == "--bc-mode" && i + 1 < argc) { bc_mode_str = argv[++i]; }
      if (arg == "--psi-init-mode" && i + 1 < argc)
      {
         psi_init_mode_str = argv[++i];
      }
      if (arg == "--order" && i + 1 < argc) { order = std::atoi(argv[++i]); }
      if (arg == "--blr-tol" && i + 1 < argc) { blr_tol = std::atof(argv[++i]); }
      if (arg == "--smooth-nucleation") { smooth_nucleation = true; }
      if (arg == "--match-quad-order") { match_quad_order = true; }
      if (arg == "--zero-dip-traction") { zero_dip_traction = true; }
      if (arg == "--elastic-sigma-n") { elastic_sigma_n = true; }
      if (arg == "--no-elastic-sigma-n") { elastic_sigma_n = false; }
      // v49 Phase 2: CFL fix and V guard
      if (arg == "--dt-init" && i + 1 < argc)
      {
         dt_init_override = std::atof(argv[++i]);
         dt_init_explicit = true;
      }
      if (arg == "--v-guard" && i + 1 < argc)
      {
         v_guard_factor = std::atof(argv[++i]);
         v_guard_explicit = true;
      }
      if (arg == "--no-v-guard")
      {
         no_v_guard = true;
         v_guard_explicit = true;
      }
      if (arg == "--tandem-time-stepping") { tandem_time_stepping = true; }
      if (arg == "--tandem-dt-init" && i + 1 < argc)
      {
         tandem_dt_init = std::atof(argv[++i]);
      }
      if (arg == "--petsc-ts") { use_petsc_ts = true; }
      if (arg == "--petsc-ts-options" && i + 1 < argc)
      {
         petsc_ts_options_file = argv[++i];
      }
      if (arg == "--penalty-factor" && i + 1 < argc) { penalty_factor = std::atof(argv[++i]); }
      if (arg == "--psi-clamp")
      {
         no_psi_clamp = false;
         psi_clamp_explicit = true;
      }
      if (arg == "--no-psi-clamp")
      {
         no_psi_clamp = true;
         psi_clamp_explicit = true;
      }
      if (arg == "--face-basis-type" && i + 1 < argc)
      {
         face_basis_str = argv[++i];
         if (face_basis_str == "GaussLobatto" || face_basis_str == "gl")
         { face_basis_type = BasisType::GaussLobatto; face_basis_str = "GaussLobatto"; }
         else if (face_basis_str == "ClosedUniform" || face_basis_str == "equi")
         { face_basis_type = BasisType::ClosedUniform; face_basis_str = "ClosedUniform"; }
         else if (face_basis_str == "ClosedGL" || face_basis_str == "cgl")
         { face_basis_type = BasisType::ClosedGL; face_basis_str = "ClosedGL"; }
         else
         { MFEM_ABORT("Unknown face-basis-type: " << face_basis_str); }
      }
   }

   if (tandem_time_stepping)
   {
      if (!dt_init_explicit) { dt_init_override = tandem_dt_init; }
      if (!v_guard_explicit) { no_v_guard = true; }
      if (!psi_clamp_explicit) { no_psi_clamp = true; }
   }

   if (use_petsc_ts)
   {
      no_v_guard = true;
      no_psi_clamp = true;
   }

#ifndef MFEM_USE_PETSC
   if (use_petsc_ts)
   {
      if (mpi.IsRoot())
      {
         std::cerr << "ERROR: --petsc-ts requires MFEM built with PETSc "
                   << "(current config has MFEM_USE_PETSC=NO).\n";
      }
      return 2;
   }
#else
   if (use_petsc_ts)
   {
      if (petsc_ts_options_file.empty())
      {
         petsc_ts_options_file = "tests/verification/petsc_ts_rk45_tandem.cfg";
      }
      MFEMInitializePetsc(&argc, &argv,
                          petsc_ts_options_file.c_str(), NULL);
      petsc_initialized = true;
   }
#endif

   // PETSc TS does not yet serialize its internal state, so
   // checkpoint/restart is not supported. Reject the combination
   // early so users don't discover it mid-run.
   if (use_petsc_ts && !restart_prefix.empty())
   {
      if (mpi.IsRoot())
      {
         std::cerr << "ERROR: --restart is not supported with --petsc-ts "
                   << "(PETSc TS state is not serialized in checkpoints).\n";
      }
      return 2;
   }

   // Parse DG method
   DGMethod dg_method = DGMethod::BR2;
   if (dg_method_str == "IP" || dg_method_str == "ip")
   {
      dg_method = DGMethod::IP;
   }
   else if (dg_method_str == "BR2" || dg_method_str == "br2")
   {
      dg_method = DGMethod::BR2;
   }

   // Process --solver flag
   SolverType solver_type = SolverType::MUMPS_BLR;
   if (solver_str == "mumps" || solver_str == "MUMPS")
   {
      solver_type = SolverType::MUMPS;
   }
   else if (solver_str == "mumps-blr" || solver_str == "MUMPS-BLR")
   {
      solver_type = SolverType::MUMPS_BLR;
   }
   else if (solver_str == "superlu" || solver_str == "SUPERLU")
   {
      solver_type = SolverType::SUPERLU;
   }
   else if (solver_str == "strumpack" || solver_str == "STRUMPACK")
   {
      solver_type = SolverType::STRUMPACK;
   }
   else if (solver_str == "gmres" || solver_str == "GMRES")
   {
      solver_type = SolverType::GMRES_BlockILU;
   }
   else if (solver_str == "gmres-amg" || solver_str == "GMRES-AMG")
   {
      solver_type = SolverType::GMRES_AMG;
   }
   else if (solver_str == "cg" || solver_str == "CG")
   {
      solver_type = SolverType::CG_AMG;
   }

   // Verify requested solver is actually available at compile time.
   // Without this check, MUMPS/SuperLU/STRUMPACK requests silently
   // fall through to CG+AMG, which may converge differently.
#ifndef MFEM_USE_MUMPS
   if (solver_type == SolverType::MUMPS || solver_type == SolverType::MUMPS_BLR)
   {
      if (mpi.IsRoot())
      {
         std::cerr << "ERROR: --solver " << solver_str
                   << " requires MFEM built with MFEM_USE_MUMPS=YES.\n"
                   << "  Available solvers: cg, gmres, gmres-amg";
#ifdef MFEM_USE_SUPERLU
         std::cerr << ", superlu";
#endif
#ifdef MFEM_USE_STRUMPACK
         std::cerr << ", strumpack";
#endif
         std::cerr << "\n";
      }
      return 2;
   }
#endif
#ifndef MFEM_USE_SUPERLU
   if (solver_type == SolverType::SUPERLU)
   {
      if (mpi.IsRoot())
      {
         std::cerr << "ERROR: --solver superlu requires MFEM built with "
                   << "MFEM_USE_SUPERLU=YES.\n";
      }
      return 2;
   }
#endif
#ifndef MFEM_USE_STRUMPACK
   if (solver_type == SolverType::STRUMPACK)
   {
      if (mpi.IsRoot())
      {
         std::cerr << "ERROR: --solver strumpack requires MFEM built with "
                   << "MFEM_USE_STRUMPACK=YES.\n";
      }
      return 2;
   }
#endif

   // Process --bc-mode flag
   BCMode bc_mode = BCMode::FarField;
   if (bc_mode_str == "far-field" || bc_mode_str == "farfield" ||
       bc_mode_str == "FarField")
   {
      bc_mode = BCMode::FarField;
   }
   else if (bc_mode_str == "x-only" || bc_mode_str == "xonly" ||
            bc_mode_str == "XOnly")
   {
      bc_mode = BCMode::XOnly;
   }
   else if (bc_mode_str == "all-dirichlet" || bc_mode_str == "AllDirichlet" ||
            bc_mode_str == "all")
   {
      bc_mode = BCMode::AllDirichlet;
   }

   // Default stations
   auto stations = BP5BenchmarkOutput<Mesh>::DefaultStations();

   // BP5 parameters
   BP5Params params;
   if (V_nuc_override > 0.0) { params.V_nuc = V_nuc_override; }
   if (delta_tau_factor_override >= 0.0)
   {
      params.delta_tau_factor = delta_tau_factor_override;
   }
   if (nucleation_eps_override >= 0.0)
   {
      params.nucleation_eps = nucleation_eps_override;
   }
   if (smooth_nucleation) { params.smooth_nucleation = true; }
   params.Validate();
   double t_final = params.t_final;
   if (tfinal_override >= 0.0) { t_final = tfinal_override; }
   params.t_final = t_final;

   std::string full_prefix = output_dir + "/" + output_prefix;

   // =========================================================================
   // Comparison-only mode
   // =========================================================================
   if (comparison_only)
   {
      if (mpi.IsRoot())
      {
         RunComparison(output_dir, output_prefix, ref_dir,
                       stations, t_final);
      }
#ifdef MFEM_USE_PETSC
      if (petsc_initialized) { MFEMFinalizePetsc(); }
#endif
      return 0;
   }

   // =========================================================================
   // Load or create mesh
   // =========================================================================
   std::unique_ptr<Mesh> serial_mesh;

   if (inline_mesh)
   {
      serial_mesh = CreateBP5InlineMesh(nx, ny, nz, Lx, Ly, Lz);
      if (mpi.IsRoot())
      {
         std::cout << "BP5 Full Simulation: Inline Mesh\n";
         std::cout << "================================\n";
         std::cout << "  nx=" << nx << " ny=" << ny << " nz=" << nz << "\n";
         std::cout << "  Lx=" << Lx/1e3 << "km Ly=" << Ly/1e3
                   << "km Lz=" << Lz/1e3 << "km\n";
         std::cout << "  Elements: " << serial_mesh->GetNE() << "\n";
      }
   }
   else if (!mesh_file.empty())
   {
      // Load Gmsh mesh
      serial_mesh = std::make_unique<Mesh>(mesh_file.c_str(), 1, 1);
      if (mesh_scale != 1.0)
      {
         for (int i = 0; i < serial_mesh->GetNV(); i++)
         {
            real_t *v = serial_mesh->GetVertex(i);
            v[0] *= mesh_scale;
            v[1] *= mesh_scale;
            v[2] *= mesh_scale;
         }
         serial_mesh->SetAttributes();
      }
      if (mpi.IsRoot())
      {
         std::cout << "BP5 Full Simulation: Gmsh Mesh\n";
         std::cout << "==============================\n";
         std::cout << "  Mesh file: " << mesh_file << "\n";
         std::cout << "  Mesh scale: " << mesh_scale << "\n";
         std::cout << "  Elements: " << serial_mesh->GetNE() << "\n";
      }
   }
   else
   {
      if (mpi.IsRoot())
      {
         std::cerr << "ERROR: --mesh <file.msh> or --inline-mesh required.\n";
      }
#ifdef MFEM_USE_PETSC
      if (petsc_initialized) { MFEMFinalizePetsc(); }
#endif
      return 1;
   }

   if (mpi.IsRoot())
   {
      std::cout << "  Ranks: " << mpi.Size() << "\n";
      std::cout << "  DG order: " << order << "\n";
      std::cout << "  DG method: " << dg_method_str << "\n";
      std::string solver_desc = "CG+AMG (iterative)";
      if (solver_type == SolverType::MUMPS) solver_desc = "MUMPS (direct)";
      else if (solver_type == SolverType::MUMPS_BLR) solver_desc = "MUMPS BLR (approximate direct)";
      else if (solver_type == SolverType::SUPERLU) solver_desc = "SuperLU_DIST (direct)";
      else if (solver_type == SolverType::STRUMPACK) solver_desc = "STRUMPACK BLR (approximate direct)";
      else if (solver_type == SolverType::GMRES_BlockILU) solver_desc = "GMRES+BlockILU (iterative)";
      std::cout << "  Solver: " << solver_desc << "\n";
      std::string bc_desc = "FarField (attrs 1-4 Dirichlet, 5-6 Natural)";
      if (bc_mode == BCMode::XOnly) bc_desc = "XOnly (attrs 1-2 Dirichlet, 3-6 Natural)";
      else if (bc_mode == BCMode::AllDirichlet) bc_desc = "AllDirichlet (all attrs Dirichlet, legacy)";
      std::cout << "  BC mode: " << bc_desc << "\n";
      std::cout << "  Psi clamp: "
                << (no_psi_clamp ? "OFF (default, Tandem-style)"
                                 : "ON [-5, 3] (--psi-clamp)")
                << "\n";
      if (tandem_time_stepping)
      {
         real_t ts_dt = (dt_init_override > 0.0) ? dt_init_override : tandem_dt_init;
         std::cout << "  Time stepping policy: Tandem-style"
                   << " (dt_init=" << ts_dt
                   << " s, no V-guard, no psi clamp)\n";
      }
      if (use_petsc_ts)
      {
         std::cout << "  Time stepping policy: PETSc TS RK45"
                   << " (exact Tandem framework";
         if (!petsc_ts_options_file.empty())
         {
            std::cout << ", options=" << petsc_ts_options_file;
         }
         std::cout << ")\n";
      }
      std::cout << "  t_final: " << t_final / BP5Params::seconds_per_year
                << " years\n";
      std::cout << "  Output prefix: " << full_prefix << "\n";
      params.Print();
      if (smooth_nucleation) { std::cout << "  Smooth nucleation: ON\n"; }
      if (match_quad_order) { std::cout << "  Match quad order (2p): ON\n"; }
      std::cout << "\n";
   }

   ParMesh pmesh(mpi.GetComm(), *serial_mesh);
   serial_mesh.reset();

   long long global_ne = pmesh.GetGlobalNE();

   // v50: Compute minimum element size for CFL-aware dt
   real_t h_min, h_max, kappa_min, kappa_max;
   pmesh.GetCharacteristics(h_min, h_max, kappa_min, kappa_max);

   if (mpi.IsRoot())
   {
      std::cout << "  ParMesh: " << global_ne << " global elements\n";
      std::cout << "  h_min = " << h_min << " m, h_max = " << h_max << " m\n";
   }

   // Dump boundary attributes to VTK for visual verification
   if (dump_bdr_vtk)
   {
      pmesh.PrintBdrVTU(output_dir + "/boundary_attributes");
      if (mpi.IsRoot())
      {
         std::cout << "  Wrote boundary VTK: " << output_dir
                   << "/boundary_attributes\n";
      }
   }

   // =========================================================================
   // Domain operator: 3D DG elasticity
   // =========================================================================
   ElasticityDomainOperator<ParMesh> domain(
      pmesh, order, params.lambda(), params.mu(),
      params.Vp, params.Wf, params.lf, dg_method, solver_type, bc_mode,
      face_basis_type);

   if (check_residual) { domain.SetCheckResidual(true); }
   if (blr_tol != 1e-10) { domain.SetBLRTol(blr_tol); }
   if (match_quad_order) { domain.SetMatchQuadOrder(true); }
   if (penalty_factor != 1.0)
   {
      domain.SetPenaltyFactor(penalty_factor);
      if (mpi.IsRoot())
      {
         std::cout << "  [v50a] penalty_factor = " << penalty_factor << "\n";
      }
   }

   if (face_basis_type != BasisType::GaussLobatto && mpi.IsRoot())
   {
      std::cout << "  [v50g] face-basis-type: " << face_basis_str << "\n";
   }
   if (mpi.IsRoot())
   {
      std::cout << "  Local fault DOFs: " << domain.GetNumFaultDOFs()
                << " (rank 0)\n";
   }

   // =========================================================================
   // Fault components
   // =========================================================================
   FaultGeometry<ParMesh> fault_geom(domain, params, &mpi);

   if (mpi.IsRoot())
   {
      std::cout << "  Global fault DOFs: " << fault_geom.NumGlobalFaultDOFs()
                << "\n";
      fault_geom.Print();
   }

   DieterichRuinaFriction::Constants fc;
   fc.V0 = params.V0;
   fc.f0 = params.f0;
   fc.b = params.b;
   fc.Dc = params.L0;  // Default Dc (per-DOF Dc from FaultGeometry overrides)
   DieterichRuinaFriction friction(fc);
   AgingLawPsi aging(params.b, params.V0, params.f0);

   RateStateFaultOperator<ParMesh, 2> fault_op(
      &fault_geom, &friction, &aging, params, &mpi);

   // Set psi initialization mode
   // Default: Tandem-style (InitialStatePsi from stress equilibrium)
   // Override: --psi-init-mode scec for SCEC fixed psi with delta_tau overstress
   if (psi_init_mode_str == "scec" || psi_init_mode_str == "SCEC")
   {
      fault_op.SetScecPsiInit(true);
      if (mpi.IsRoot())
      {
         std::cout << "  Psi init mode: SCEC (genuine delta_tau overstress)\n";
      }
   }
   else
   {
      // Default: Tandem-style (scec_psi_init_ already false by default)
      if (mpi.IsRoot())
      {
         std::cout << "  Psi init mode: Tandem (equilibrium psi from stress)\n";
      }
   }

   if (monitor_traction > 0)
   {
      fault_op.SetTractionMonitoring(monitor_traction);
   }

   // =========================================================================
   // SEAS quasi-dynamic operator
   // =========================================================================
   PBP5SEASOp seas_op(&domain, &fault_op, &mpi);

   // v51: dip traction diagnostics
   if (zero_dip_traction)
   {
      seas_op.SetZeroDipTraction(true);
      if (mpi.IsRoot()) { std::cout << "  [v51] zero-dip-traction: ON\n"; }
   }
   // v54: elastic sigma_n ON by default (Tandem DieterichRuinaBase.h:87)
   seas_op.SetElasticSigmaN(elastic_sigma_n);
   if (mpi.IsRoot())
   {
      std::cout << "  [v54] elastic-sigma-n: "
                << (elastic_sigma_n ? "ON (default, Tandem)" : "OFF (--no-elastic-sigma-n)")
                << "\n";
   }

   // Face tracer for per-face diagnostics (rank-local, no MPI)
   seas::TraceConfig trace_cfg;
   trace_cfg.explicit_rank = 96;
   trace_cfg.explicit_fi = {21, 28, 29, 31, 33, 37};
   trace_cfg.use_coord_window = true;
   trace_cfg.x2_min = -45e3; trace_cfg.x2_max = -25e3;
   trace_cfg.x3_min = -40e3; trace_cfg.x3_max = -35e3;
   trace_cfg.num_control_faces = 2;
   trace_cfg.output_dir = output_dir;
   seas::FaceTraceLogger<ParMesh> face_tracer(trace_cfg, mpi.Rank());
   face_tracer.SelectFaces(domain, fault_geom);
   seas_op.SetFaceTracer(&face_tracer);

   Vector state(fault_op.StateSize());
   seas_op.SetInitialCondition(state);

   real_t V_init = seas_op.GetMaxSlipRate();
   if (mpi.IsRoot())
   {
      std::cout << "\nInitial conditions:\n";
      std::cout << "  V_max = " << V_init << " m/s\n";
      std::cout << "  StateSize = " << fault_op.StateSize() << "\n";
      std::cout << "  SlipSize = " << fault_op.SlipSize() << "\n";
   }

   // =========================================================================
   // I/O: distributed probe output (Tandem-style)
   // =========================================================================
   // Each rank uses the owned local fault coords for probe location.
   // This follows Tandem's owned-state + ghost-read layout.
   Vector local_x2 = fault_geom.GetCoordsX2();
   Vector local_x3 = fault_geom.GetCoordsX3();

   int N_local = fault_geom.NumLocalFaultDOFs();
   Vector local_tp_dip(N_local), local_tp_strike(N_local);
   const Vector &local_tau_pre = fault_geom.GetTauPre();
   for (int i = 0; i < N_local; i++)
   {
      local_tp_dip(i) = local_tau_pre(2 * i);
      local_tp_strike(i) = local_tau_pre(2 * i + 1);
   }

   // Clean up old output files
   if (mpi.IsRoot())
   {
      for (const auto &st : stations)
      {
         std::string fn = full_prefix + "_" + st.name + ".txt";
         std::remove(fn.c_str());
      }
      std::remove((full_prefix + "_global.txt").c_str());
   }
   mpi.Barrier();

   ParallelBP5BenchmarkOutput bench_out(
      full_prefix, params, stations, fault_geom, mpi,
      local_x2, local_x3, local_tp_dip, local_tp_strike,
      domain.GetNbfPerFace(), face_basis_type);

   if (mpi.IsRoot())
   {
      bench_out.PrintDiagnostics(local_x2, local_x3);
   }

   // Global output (root only)
   std::unique_ptr<ProbeOutput> global_out;
   if (mpi.IsRoot())
   {
      global_out = std::make_unique<ProbeOutput>(
         full_prefix + "_global.txt",
         std::vector<std::string>{"time(s)", "log10(Vmax)(m/s)"},
         "BP5-QD global output");
   }

   // Write initial state
   bench_out.ForceWrite(0.0, state, fault_op, seas_op.GetTraction(), V_init);
   bench_out.Flush();
   if (mpi.IsRoot() && global_out)
   {
      global_out->WriteStep({0.0, V_init > 0.0 ? std::log10(V_init) : -300.0});
   }

   // =========================================================================
   // Time integration (Dormand-Prince RK45)
   // =========================================================================
   DormandPrinceRK45 ode_solver;
   ode_solver.SetMPIContext(&mpi);
   ode_solver.SetAbsTol(1e-7);
   ode_solver.SetRelTol(1e-50);
   ode_solver.SetDtMin(1e-6);
   ode_solver.SetDtMax(0.1 * BP5Params::seconds_per_year);

   // v50: CFL-aware initial dt selection.
   // The IP penalty creates a CFL-like stability constraint: dt must satisfy
   //   z = lambda_eff * dt < z_crit  where lambda_eff = beta * mu / (eta * h)
   // Exceeding z_crit causes exponential RK stage amplification (cascade).
   // See bp5_debug_v49.md Section 4 for full derivation.
   //
   // Two dt limits:
   //   dt_V   = 0.01 * Dc / V_max        (physics: slip per step << Dc)
   //   dt_CFL = C * eta * h_min / (beta(p) * mu)  (stability: z < z_crit)
   // with C=2.0 (safety below empirical z_crit≈2.5-3.0)
   //
   // v50+: beta scales with c_N_1 = p*(p+dim-1)/dim (penalty coefficient).
   //   beta(p) = beta_ref * c_N_1(p) / c_N_1(p_ref)
   // where beta_ref=4.0 was calibrated at p_ref=2.
   //   p=1: c_N_1=1.0   -> beta=1.5
   //   p=2: c_N_1=2.67  -> beta=4.0
   //   p=4: c_N_1=8.0   -> beta=12.0
   //   p=6: c_N_1=16.0  -> beta=24.0
   int dim = 3;
   real_t c_N_1 = order * (order + dim - 1.0) / dim;
   real_t c_N_1_ref = 2.0 * (2.0 + dim - 1.0) / dim;  // c_N_1 at p=2 = 8/3
   real_t beta = 4.0 * c_N_1 / c_N_1_ref;
   real_t V_max_init = std::max(V_init, params.V_nuc);
   real_t dt_V = std::min(1e3, 0.01 * params.L_nuc /
                          std::max(V_max_init, 1e-20));
   real_t dt_CFL = 2.0 * params.eta() * h_min / (beta * params.mu());
   real_t dt_init = std::min(dt_V, dt_CFL);

   // Manual override if --dt-init flag provided (for testing)
   if (dt_init_override > 0)
   {
      dt_init = dt_init_override;
      if (mpi.IsRoot())
      {
         std::cout << "  [override] dt_init = " << dt_init << " s\n";
      }
   }
   real_t current_dt = dt_init;
   int step_rejections = 0;
   int max_steps = 10000000;
   Vector empty_k0;

   if (!use_petsc_ts)
   {
      ode_solver.SetDt(dt_init);
      if (mpi.IsRoot())
      {
         std::cout << "  CFL: c_N_1=" << c_N_1 << " beta=" << beta
                   << " (order=" << order << ")\n";
         std::cout << "  dt_V = " << dt_V << " s, dt_CFL = " << dt_CFL << " s\n";
         std::cout << "  Initial dt: " << dt_init << " s"
                   << (dt_init <= dt_CFL ? " (CFL-limited)" : " (V-limited)")
                   << "\n";
      }
      ode_solver.SetStatePerNode(3);  // BP5: [slip_dip, slip_strike, psi]

      // v50: V-guard ON by default (factor=100). Prevents RK cascade overflow.
      // Use --v-guard <factor> to change threshold, --no-v-guard to disable.
      if (!no_v_guard)
      {
         real_t factor = (v_guard_factor > 0) ? v_guard_factor : 100.0;
         ode_solver.SetVGuard(factor);
         if (mpi.IsRoot())
         {
            std::cout << "  V-guard: ON (factor=" << factor << ")\n";
         }
      }
      else if (mpi.IsRoot())
      {
         std::cout << "  V-guard: OFF (--no-v-guard)\n";
      }
      ode_solver.Init(seas_op);
   }
#ifdef MFEM_USE_PETSC
   std::unique_ptr<PetscODESolver> petsc_ode;
   // Persistent PETSc vector for state — avoids MFEM's PlaceMemory/
   // ResetMemory cycle which can crash with certain PETSc versions.
   std::unique_ptr<PetscParVector> petsc_state;
   if (use_petsc_ts)
   {
      // No prefix: options file uses unprefixed names (-ts_type, etc.)
      // so the solver must also be unprefixed.
      petsc_ode = std::make_unique<PetscODESolver>(mpi.GetComm(), "");
      petsc_ode->SetAbsTol(1e-7);
      petsc_ode->SetRelTol(1e-50);
      petsc_ode->SetMaxIter(max_steps);
      petsc_ode->Init(seas_op, PetscODESolver::ODE_SOLVER_GENERAL);
      petsc::TS ts = *petsc_ode;

      // MFEM's PetscODESolver constructor disables TS adaptivity
      // (TSAdaptSetType(TSADAPTNONE) in petsc.cpp:4194). Re-enable
      // adaptive stepping so PETSc's RK45 error control works.
      // TSSetFromOptions below will pick up -ts_adapt_* from the cfg.
      {
         TSAdapt tsad;
         PetscErrorCode ierr2 = TSGetAdapt(ts, &tsad);
         MFEM_VERIFY(ierr2 == PETSC_SUCCESS, "TSGetAdapt failed");
         ierr2 = TSAdaptSetType(tsad, TSADAPTBASIC);
         MFEM_VERIFY(ierr2 == PETSC_SUCCESS, "TSAdaptSetType(BASIC) failed");
      }

      // Apply all PETSc options from the cfg file (-ts_type, -ts_rk_type, etc.)
      PetscErrorCode ierr = TSSetFromOptions(ts);
      MFEM_VERIFY(ierr == PETSC_SUCCESS, "TSSetFromOptions failed");

      ierr = TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP);
      MFEM_VERIFY(ierr == PETSC_SUCCESS, "TSSetExactFinalTime(MATCHSTEP) failed");
      ierr = TSSetMaxTime(ts, t_final);
      MFEM_VERIFY(ierr == PETSC_SUCCESS, "TSSetMaxTime() failed");

      // Create a persistent PETSc-managed copy of the state vector.
      // PetscODESolver::Step uses PlaceMemory/ResetMemory which fails
      // with PETSc 3.15. Instead, we copy state ↔ petsc_state at each
      // step and call TSStep directly.
      petsc_state = std::make_unique<PetscParVector>(
         mpi.GetComm(), state, true);
      ierr = TSSetSolution(ts, *petsc_state);
      MFEM_VERIFY(ierr == PETSC_SUCCESS, "TSSetSolution() failed");

      current_dt = dt_init;
      if (mpi.IsRoot())
      {
         std::cout << "  PETSc TS initial dt: "
                   << ((dt_init_override > 0.0) ? dt_init : 0.1)
                   << " s";
         if (dt_init_override <= 0.0)
         {
            std::cout << " (PETSc default if not overridden)";
         }
         std::cout << "\n";
      }
   }
#endif

   real_t t = 0.0;
   int step = 0;

   // Earthquake detection
   bool in_seismic_event = false;
   int num_seismic_events = 0;
   real_t V_threshold_seismic = 1e-3;
   real_t V_threshold_interseismic = 1e-6;
   long long psi_clamp_hits_hi = 0;
   long long psi_clamp_hits_lo = 0;
   int psi_clamp_first_step = -1;
   real_t psi_clamp_first_time = -1.0;
   real_t psi_clamp_max_hi = -std::numeric_limits<real_t>::infinity();
   real_t psi_clamp_min_lo = std::numeric_limits<real_t>::infinity();

   // =========================================================================
   // Restart from checkpoint (if requested)
   // =========================================================================
   if (!restart_prefix.empty())
   {
      real_t restart_dt;
      Vector restart_disp, restart_traction, restart_slip_rate, restart_k0;
      bool restart_fsal;

      bool ok = ReadCheckpoint(restart_prefix, t, restart_dt,
                               step, num_seismic_events, in_seismic_event,
                               state, restart_disp, restart_traction,
                               restart_slip_rate, restart_fsal, restart_k0,
                               &mpi);
      MFEM_VERIFY(ok, "Failed to load checkpoint: " << restart_prefix);

      current_dt = restart_dt;
      seas_op.SetDisplacement(restart_disp);
      fault_op.SetSlipRate(restart_slip_rate);
      if (!use_petsc_ts)
      {
         ode_solver.SetDt(restart_dt);
         if (restart_fsal && restart_k0.Size() > 0)
         {
            ode_solver.RestoreFSAL(restart_k0);
         }
      }

      if (mpi.IsRoot())
      {
         std::cout << "\nRestarted from checkpoint:\n";
         std::cout << "  Time: " << t / BP5Params::seconds_per_year << " yr\n";
         std::cout << "  Step: " << step << "\n";
         std::cout << "  dt: " << restart_dt << " s\n";
      }
   }

   if (mpi.IsRoot())
   {
      std::cout << "\n" << std::setw(10) << "Step"
                << std::setw(16) << "Time [yr]"
                << std::setw(14) << "dt [s]"
                << std::setw(16) << "V_max [m/s]"
                << std::setw(8) << "EQs"
                << "\n";
      std::cout << std::string(64, '-') << "\n";
   }

   int print_step_interval = 10;

   // =========================================================================
   // Main time-stepping loop
   // =========================================================================
   while (t < t_final && step < max_steps)
   {
      if (!use_petsc_ts && t + ode_solver.GetDt() > t_final)
      {
         ode_solver.SetDt(t_final - t);
      }

      real_t dt;
      bool accepted = true;
      if (!use_petsc_ts)
      {
         accepted = ode_solver.Step(seas_op, state, t, dt);
         if (!accepted) { continue; }
         current_dt = ode_solver.GetDt();
      }
#ifdef MFEM_USE_PETSC
      else
      {
         MFEM_VERIFY(petsc_ode && petsc_state,
                     "PETSc TS solver was not initialized");
         petsc::TS ts = *petsc_ode;
         PetscErrorCode ierr;

         // Copy state → PETSc Vec (bypass MFEM PlaceMemory)
         {
            PetscScalar *arr;
            ierr = VecGetArray(*petsc_state, &arr);
            MFEM_VERIFY(ierr == PETSC_SUCCESS, "VecGetArray failed");
            for (int i = 0; i < state.Size(); i++) { arr[i] = state(i); }
            ierr = VecRestoreArray(*petsc_state, &arr);
            MFEM_VERIFY(ierr == PETSC_SUCCESS, "VecRestoreArray failed");
         }

         // Take one PETSc TS step directly
         dt = current_dt;
         ierr = TSSetTime(ts, t);
         MFEM_VERIFY(ierr == PETSC_SUCCESS, "TSSetTime failed");
         ierr = TSSetTimeStep(ts, dt);
         MFEM_VERIFY(ierr == PETSC_SUCCESS, "TSSetTimeStep failed");
         ierr = TSStep(ts);
         MFEM_VERIFY(ierr == PETSC_SUCCESS, "TSStep failed");

         // Read back updated time and state
         PetscReal pt;
         ierr = TSGetTime(ts, &pt);
         MFEM_VERIFY(ierr == PETSC_SUCCESS, "TSGetTime failed");
         dt = pt - t;
         t = pt;

         // Copy PETSc Vec → state
         {
            const PetscScalar *arr;
            ierr = VecGetArrayRead(*petsc_state, &arr);
            MFEM_VERIFY(ierr == PETSC_SUCCESS, "VecGetArrayRead failed");
            for (int i = 0; i < state.Size(); i++) { state(i) = arr[i]; }
            ierr = VecRestoreArrayRead(*petsc_state, &arr);
            MFEM_VERIFY(ierr == PETSC_SUCCESS, "VecRestoreArrayRead failed");
         }

         PetscReal next_dt = 0.0;
         ierr = TSGetTimeStep(ts, &next_dt);
         MFEM_VERIFY(ierr == PETSC_SUCCESS, "TSGetTimeStep() failed");
         current_dt = next_dt;
         PetscInt rejects = 0;
         ierr = TSGetStepRejections(ts, &rejects);
         MFEM_VERIFY(ierr == PETSC_SUCCESS, "TSGetStepRejections() failed");
         step_rejections = static_cast<int>(rejects);
      }
#endif
      step++;

      // Post-step psi clamping is MFEM-specific. Tandem does not do this, so
      // keep it switchable and instrumented while debugging dynamic mismatch.
      if (!no_psi_clamp)
      {
         const int spn = 3;  // BP5: [slip_dip, slip_strike, psi]
         const int psi_idx = 2;
         const real_t psi_max = 3.0;
         const real_t psi_min = -5.0;
         int local_hi = 0;
         int local_lo = 0;
         real_t local_max_hi = -std::numeric_limits<real_t>::infinity();
         real_t local_min_lo = std::numeric_limits<real_t>::infinity();
         int n_nodes = state.Size() / spn;
         for (int i = 0; i < n_nodes; i++)
         {
            real_t &psi = state(i * spn + psi_idx);
            if (psi > psi_max)
            {
               local_hi++;
               local_max_hi = std::max(local_max_hi, psi);
               psi = psi_max;
            }
            else if (psi < psi_min)
            {
               local_lo++;
               local_min_lo = std::min(local_min_lo, psi);
               psi = psi_min;
            }
         }

         int global_hi = mpi.GlobalSumInt(local_hi);
         int global_lo = mpi.GlobalSumInt(local_lo);
         real_t global_max_hi = (global_hi > 0)
                                  ? mpi.GlobalMax(local_max_hi)
                                  : -std::numeric_limits<real_t>::infinity();
         real_t global_min_lo = (global_lo > 0)
                                  ? mpi.GlobalMin(local_min_lo)
                                  : std::numeric_limits<real_t>::infinity();

         if (global_hi > 0 || global_lo > 0)
         {
            psi_clamp_hits_hi += global_hi;
            psi_clamp_hits_lo += global_lo;
            if (global_hi > 0)
            {
               psi_clamp_max_hi = std::max(psi_clamp_max_hi, global_max_hi);
            }
            if (global_lo > 0)
            {
               psi_clamp_min_lo = std::min(psi_clamp_min_lo, global_min_lo);
            }
            if (psi_clamp_first_step < 0)
            {
               psi_clamp_first_step = step + 1;
               psi_clamp_first_time = t;
            }
         }
      }

      // Face tracer: commit after psi clamp so psi reflects the carried state.
      // Note: traction/rate/correction were staged in Mult() before clamping,
      // so rows are only fully self-consistent when psi clamp is off (default).
      if (face_tracer.IsActive())
      {
         face_tracer.CommitStep(step, t, state, false);
      }

      real_t V_max = seas_op.GetMaxSlipRate();

      // NaN/Inf check
      bool has_nan = !std::isfinite(V_max);
      if (!has_nan)
      {
         for (int i = 0; i < std::min(state.Size(), 100); i++)
         {
            if (!std::isfinite(state(i))) { has_nan = true; break; }
         }
      }
      {
         int local_nan = has_nan ? 1 : 0;
         int global_nan = mpi.GlobalSumInt(local_nan);
         has_nan = (global_nan > 0);
      }
      if (has_nan)
      {
         if (mpi.IsRoot())
         {
            std::cerr << "NaN/Inf detected at step " << step
                      << ", t = " << t / BP5Params::seconds_per_year
                      << " yr, V_max = " << V_max << "\n";
         }
         break;
      }

      // Earthquake detection
      if (!in_seismic_event && V_max > V_threshold_seismic)
      {
         in_seismic_event = true;
         num_seismic_events++;
         if (mpi.IsRoot())
         {
            std::cout << "  *** EARTHQUAKE #" << num_seismic_events
                      << " at t = " << std::fixed << std::setprecision(1)
                      << t / BP5Params::seconds_per_year << " yr"
                      << ", V_max = " << std::scientific
                      << std::setprecision(2) << V_max << " m/s ***\n";
         }
      }
      else if (in_seismic_event && V_max < V_threshold_interseismic)
      {
         in_seismic_event = false;
         if (mpi.IsRoot())
         {
            std::cout << "  Earthquake #" << num_seismic_events
                      << " resolved at t = "
                      << std::fixed << std::setprecision(1)
                      << t / BP5Params::seconds_per_year << " yr\n";
         }
      }

      // I/O: adaptive schedule or every step
      if (write_every_step)
      {
         bench_out.ForceWrite(t, state, fault_op, seas_op.GetTraction(),
                              V_max);

         bench_out.Flush();
      }
      else if (bench_out.Write(t, state, fault_op, seas_op.GetTraction(),
                               V_max))
      {

         bench_out.Flush();
      }

      // Global output
      if (mpi.IsRoot() && global_out)
      {
         global_out->WriteStep(
            {t, V_max > 0.0 ? std::log10(V_max) : -300.0});
      }

      // Checkpoint
      if (checkpoint_interval > 0 && step % checkpoint_interval == 0)
      {
         Vector u_vec;
         u_vec = seas_op.GetDisplacement();
         WriteCheckpoint(full_prefix, t, current_dt,
                         step, num_seismic_events, in_seismic_event,
                         state, u_vec,
                         seas_op.GetTraction(), fault_op.GetSlipRate(),
                         use_petsc_ts ? false : ode_solver.IsInitialized(),
                         use_petsc_ts ? empty_k0 : ode_solver.GetK0(),
                         &mpi);
      }

      // Console output
      if (mpi.IsRoot() &&
          (step % print_step_interval == 0 || V_max > V_threshold_seismic))
      {
         std::cout << std::setw(10) << step
                   << std::setw(16) << std::scientific << std::setprecision(6)
                   << t / BP5Params::seconds_per_year
                   << std::setw(14) << std::scientific << std::setprecision(3)
                   << dt
                   << std::setw(16) << std::scientific << std::setprecision(3)
                   << V_max
                   << std::setw(8) << num_seismic_events
                   << "\n";
         std::cout.flush();
      }
   }

   // =========================================================================
   // Final checkpoint and close
   // =========================================================================
   if (checkpoint_interval > 0)
   {
      Vector u_vec;
      u_vec = seas_op.GetDisplacement();
      WriteCheckpoint(full_prefix, t, current_dt,
                      step, num_seismic_events, in_seismic_event,
                      state, u_vec,
                      seas_op.GetTraction(), fault_op.GetSlipRate(),
                      use_petsc_ts ? false : ode_solver.IsInitialized(),
                      use_petsc_ts ? empty_k0 : ode_solver.GetK0(),
                      &mpi);
   }

   // Face tracer: finalize (commit final step + write summary)
   if (face_tracer.IsActive())
   {
      face_tracer.CommitStep(step, t, state, true);
   }
   face_tracer.Finalize();

   bench_out.ForceWrite(t, state, fault_op, seas_op.GetTraction(),
                        seas_op.GetMaxSlipRate());
   bench_out.Close();
   if (mpi.IsRoot() && global_out)
   {
      global_out->Flush();
      global_out->Close();
   }

   if (mpi.IsRoot())
   {
      std::cout << "\n=== BP5 Simulation Summary ===\n";
      std::cout << "  Ranks: " << mpi.Size() << "\n";
      std::cout << "  Final time: " << t / BP5Params::seconds_per_year
                << " years\n";
      std::cout << "  Total steps: " << step << "\n";
      std::cout << "  Step rejections: "
                << (use_petsc_ts ? step_rejections
                                 : ode_solver.GetTotalRejections()) << "\n";
      std::cout << "  Seismic events: " << num_seismic_events << "\n\n";
      if (!no_psi_clamp)
      {
         std::cout << "  Psi clamp hits (high): " << psi_clamp_hits_hi << "\n";
         std::cout << "  Psi clamp hits (low):  " << psi_clamp_hits_lo << "\n";
         if (psi_clamp_first_step >= 0)
         {
            std::cout << "  First psi clamp: step " << psi_clamp_first_step
                      << ", t = " << psi_clamp_first_time << " s\n";
         }
         else
         {
            std::cout << "  First psi clamp: none\n";
         }
         if (psi_clamp_hits_hi > 0)
         {
            std::cout << "  Max unclamped psi: " << psi_clamp_max_hi << "\n";
         }
         if (psi_clamp_hits_lo > 0)
         {
            std::cout << "  Min unclamped psi: " << psi_clamp_min_lo << "\n";
         }
         std::cout << "\n";
      }
   }

   // =========================================================================
   // Post-simulation comparison
   // =========================================================================
   if (mpi.IsRoot())
   {
      RunComparison(output_dir, output_prefix, ref_dir,
                    stations, t_final);
   }

#ifdef MFEM_USE_PETSC
   if (petsc_initialized)
   {
      MFEMFinalizePetsc();
   }
#endif

   return 0;
}
