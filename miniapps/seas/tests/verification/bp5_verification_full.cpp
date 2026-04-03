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
//   --diag-vtk                 Output diagnostic VTK: displacement, fault a,
//                              tau_pre, V_init, and boundary attributes
//   --diag-traction-decomp     Print traction decomposition (stress vs penalty
//                              correction vs jump) for each fault DOF
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
   bool diag_vtk = false;
   bool diag_traction_decomp = false;
   std::string bc_mode_str = "far-field";
   std::string psi_init_mode_str = "tandem";  // "tandem" (default) or "scec"
   int order = 1;
   double blr_tol = 1e-10;  // MUMPS-BLR tolerance (default 1e-10)

   // v49 Phase 1 diagnostic flags
   bool smooth_nucleation = false;
   bool match_quad_order = false;
   bool diag_normals = false;
   bool diag_first_traction = false;
   bool diag_rk_stages = false;
   bool diag_dip_traction = false;   // v51: per-component traction diagnostic
   bool zero_dip_traction = false;   // v51: zero tau_dip after ComputeTraction
   bool diag_coseismic_dip = false;  // v51: dump dip/strike ratio during coseismic
   real_t coseismic_dip_threshold = 0.1;  // v51: V_max threshold for coseismic dump
   bool diag_uz_fault = false;       // v51: dump u_z at fault faces after solve
   // v54: elastic sigma_n ON by default (matches Tandem DieterichRuinaBase.h:87)
   bool elastic_sigma_n = true;
   bool diag_traction_coherence = false; // v52: solve vs traction coherence test
   bool diag_rhs_z = false;              // v52: dump f_z components of RHS
   // v49 Phase 2: CFL-aware dt and V guard
   real_t dt_init_override = -1.0;  // Manual dt_init override (negative = auto)
   real_t v_guard_factor = -1.0;    // V guard threshold factor (negative = use default 100)
   bool no_v_guard = false;         // v50: disable V guard (for testing only)
   bool tandem_time_stepping = false; // Use Tandem-style startup/acceptance policy
   real_t tandem_dt_init = 0.01;      // Default Tandem-style startup dt [s]
   // v50a: penalty scaling factor (1.0 = default, <1.0 = reduced penalty)
   real_t penalty_factor = 1.0;
   bool diag_station_traction_decomp = false; // Write station-level stress/correction traction
   bool diag_station_jump_residual = false;   // Write station-level [[u]]-delta residual
   bool no_psi_clamp = true;           // Default OFF: match Tandem (no post-step psi clamp)
   bool diag_psi_clamp = false;        // Report psi clamp activation statistics
   bool dt_init_explicit = false;
   bool v_guard_explicit = false;
   bool psi_clamp_explicit = false;
   bool use_petsc_ts = false;          // Exact Tandem framework: PETSc TS
   bool diag_tip_step1 = false;        // v58: one-step tip reproducer then exit
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
      if (arg == "--diag-traction-decomp") { diag_traction_decomp = true; }
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
      if (arg == "--diag-vtk") { diag_vtk = true; }
      if (arg == "--bc-mode" && i + 1 < argc) { bc_mode_str = argv[++i]; }
      if (arg == "--psi-init-mode" && i + 1 < argc)
      {
         psi_init_mode_str = argv[++i];
      }
      if (arg == "--order" && i + 1 < argc) { order = std::atoi(argv[++i]); }
      if (arg == "--blr-tol" && i + 1 < argc) { blr_tol = std::atof(argv[++i]); }
      // v49 Phase 1 diagnostic flags
      if (arg == "--smooth-nucleation") { smooth_nucleation = true; }
      if (arg == "--match-quad-order") { match_quad_order = true; }
      if (arg == "--diag-normals") { diag_normals = true; }
      if (arg == "--diag-first-traction") { diag_first_traction = true; }
      if (arg == "--diag-rk-stages") { diag_rk_stages = true; }
      if (arg == "--diag-dip-traction") { diag_dip_traction = true; }
      if (arg == "--zero-dip-traction") { zero_dip_traction = true; }
      if (arg == "--diag-coseismic-dip") { diag_coseismic_dip = true; }
      if (arg == "--diag-coseismic-dip-threshold" && i + 1 < argc)
      {
         diag_coseismic_dip = true;
         coseismic_dip_threshold = std::atof(argv[++i]);
      }
      if (arg == "--diag-uz-fault") { diag_uz_fault = true; }
      if (arg == "--elastic-sigma-n") { elastic_sigma_n = true; }
      if (arg == "--no-elastic-sigma-n") { elastic_sigma_n = false; }
      if (arg == "--diag-traction-coherence") { diag_traction_coherence = true; }
      if (arg == "--diag-rhs-z") { diag_rhs_z = true; }
      if (arg == "--diag-tip-step1") { diag_tip_step1 = true; }
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
      if (arg == "--diag-station-traction-decomp")
      {
         diag_station_traction_decomp = true;
      }
      if (arg == "--diag-station-jump-residual")
      {
         diag_station_jump_residual = true;
      }
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
      if (arg == "--diag-psi-clamp") { diag_psi_clamp = true; }
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
                << (diag_psi_clamp ? " [diagnostic]" : "") << "\n";
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
      // v49 Phase 1 diagnostic flags
      if (smooth_nucleation) { std::cout << "  [v49] Smooth nucleation: ON\n"; }
      if (match_quad_order) { std::cout << "  [v49] Match quad order (2p): ON\n"; }
      if (diag_normals) { std::cout << "  [v49] Diag normals: ON\n"; }
      if (diag_first_traction) { std::cout << "  [v49] Diag first traction: ON\n"; }
      if (diag_rk_stages) { std::cout << "  [v49] Diag RK stages: ON\n"; }
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
   if (diag_traction_decomp) { domain.SetDiagTractionDecomp(true); }
   if (blr_tol != 1e-10) { domain.SetBLRTol(blr_tol); }
   // v49 Phase 1 flags
   if (match_quad_order) { domain.SetMatchQuadOrder(true); }
   if (diag_normals) { domain.SetDiagNormals(true); }
   if (diag_first_traction) { domain.SetDiagFirstTraction(true); }
   if (diag_dip_traction) { domain.SetDiagDipTraction(true); }
   if (diag_uz_fault) { domain.SetDiagUzFault(true); }
   if (diag_traction_coherence) { domain.SetDiagTractionCoherence(true); }
   if (diag_rhs_z) { domain.SetDiagRhsZ(true); }
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
   if (diag_coseismic_dip)
   {
      seas_op.SetDiagCoseismicDip(true, coseismic_dip_threshold);
      if (mpi.IsRoot()) { std::cout << "  [v51] diag-coseismic-dip: ON (threshold=" << coseismic_dip_threshold << " m/s)\n"; }
   }
   // v54: elastic sigma_n ON by default (Tandem DieterichRuinaBase.h:87)
   seas_op.SetElasticSigmaN(elastic_sigma_n);
   if (mpi.IsRoot())
   {
      std::cout << "  [v54] elastic-sigma-n: "
                << (elastic_sigma_n ? "ON (default, Tandem)" : "OFF (--no-elastic-sigma-n)")
                << "\n";
   }

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

   if (diag_station_traction_decomp)
   {
      bench_out.EnableTractionDecompositionOutput();
      if (mpi.IsRoot())
      {
         std::cout << "  Diagnostic station traction decomposition: ON\n";
      }
   }
   if (diag_station_jump_residual)
   {
      bench_out.EnableJumpResidualOutput();
      if (mpi.IsRoot())
      {
         std::cout << "  Diagnostic station jump residual: ON\n";
      }
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
   // Diagnostic VTK output
   // =========================================================================
   if (diag_vtk)
   {
      if (mpi.IsRoot())
      {
         std::cout << "\n=== Diagnostic VTK Output ===\n";
      }

      // Solve domain at t=1yr with zero slip to show Dirichlet BC pattern
      real_t diag_t = BP5Params::seconds_per_year;
      int N_loc = fault_geom.NumLocalFaultDOFs();
      Vector zero_slip(2 * N_loc);
      zero_slip = 0.0;

      ParGridFunction u_diag(&domain.GetFESpace());
      u_diag = 0.0;
      domain.Solve(diag_t, zero_slip, u_diag);

      if (mpi.IsRoot())
      {
         std::cout << "  Dirichlet BC solve at t=1yr: |u| = "
                   << u_diag.Norml2() << "\n";
      }

      // === Test 2 diagnostic: boundary loading traction verification ===
      // Compute traction at fault from boundary-only solve (zero slip)
      {
         Vector diag_traction(2 * N_loc);
         domain.ComputeTraction(u_diag, zero_slip, diag_traction);

         // Gather to root for analysis
         Vector local_trac_dip(N_loc), local_trac_strike(N_loc);
         for (int i = 0; i < N_loc; i++)
         {
            local_trac_dip(i) = diag_traction(2 * i);
            local_trac_strike(i) = diag_traction(2 * i + 1);
         }

         Vector global_trac_dip, global_trac_strike;
         fault_geom.GatherToRoot(local_trac_dip, global_trac_dip);
         fault_geom.GatherToRoot(local_trac_strike, global_trac_strike);

         // Gather coords for this diagnostic only
         Vector diag_gx2, diag_gx3;
         fault_geom.GatherToRoot(local_x2, diag_gx2);
         fault_geom.GatherToRoot(local_x3, diag_gx3);

         if (mpi.IsRoot())
         {
            // Analytical estimate: tau = mu*Vp*1yr/(2*Lx)
            real_t mu = params.mu();
            real_t tau_analytical = mu * params.Vp * BP5Params::seconds_per_year
                                    / (2.0 * 100e3);
            std::cout << "\n  === Boundary Loading Traction Verification ===\n";
            std::cout << "  Analytical estimate (2D antiplane): "
                      << tau_analytical << " Pa = "
                      << tau_analytical / 1e6 << " MPa\n";

            int M = global_trac_strike.Size();
            // Find station nearest to (x2=0, x3=10km) and (x2=0, x3=22km)
            for (int target_z : {0, 10000, 22000})
            {
               int best = -1;
               real_t best_dist = 1e30;
               for (int j = 0; j < M; j++)
               {
                  real_t d = std::abs(diag_gx2(j)) +
                             std::abs(diag_gx3(j) - target_z);
                  if (d < best_dist) { best_dist = d; best = j; }
               }
               if (best >= 0)
               {
                  real_t ts = global_trac_strike(best);
                  real_t td = global_trac_dip(best);
                  std::cout << "  z=" << target_z/1000 << "km: trac_strike="
                            << ts << " Pa (" << ts/1e6 << " MPa)"
                            << "  trac_dip=" << td << " Pa (" << td/1e6 << " MPa)"
                            << "  ratio_to_analytical=" << ts/tau_analytical
                            << "\n";
               }
            }
            std::cout << "\n";
         }
      }

      // Create L2 p=0 fields for fault parameter visualization
      L2_FECollection l2_fec(0, 3);
      ParFiniteElementSpace l2_fes(&pmesh, &l2_fec);

      // --- Centroid-based analytic fields (kept from original) ---
      ParGridFunction a_field(&l2_fes);
      a_field = 0.0;
      ParGridFunction tau_pre_mag(&l2_fes);
      tau_pre_mag = 0.0;
      ParGridFunction V_init_field(&l2_fes);
      V_init_field = 0.0;

      // Map fault parameters to elements adjacent to fault plane (x=0)
      int n_fault_elem = 0;
      for (int i = 0; i < pmesh.GetNE(); i++)
      {
         Array<int> verts;
         pmesh.GetElementVertices(i, verts);
         real_t min_x = 1e30, max_x = -1e30;
         Vector center(3);
         center = 0.0;
         for (int v = 0; v < verts.Size(); v++)
         {
            const real_t *coords = pmesh.GetVertex(verts[v]);
            min_x = std::min(min_x, coords[0]);
            max_x = std::max(max_x, coords[0]);
            for (int d = 0; d < 3; d++) { center(d) += coords[d]; }
         }
         center /= verts.Size();

         if (min_x > 0.0 || max_x < 0.0) { continue; }

         real_t y = center(1);
         real_t z = center(2);

         if (std::abs(y) > params.lf / 2.0 || z > params.Wf) { continue; }

         a_field(i) = params.a_of_x2_x3(y, z);

         real_t tau[2];
         params.tau0_vec(y, z, tau);
         tau_pre_mag(i) = std::sqrt(tau[0] * tau[0] + tau[1] * tau[1]);

         real_t V[2];
         params.V_init_vec(y, z, V);
         V_init_field(i) = std::sqrt(V[0] * V[0] + V[1] * V[1]);

         n_fault_elem++;
      }

      if (mpi.IsRoot())
      {
         std::cout << "  Fault-adjacent elements (centroid, rank 0): "
                   << n_fault_elem << "\n";
      }

      // --- State-based fields via fault face → element mapping ---
      const Array<int> &fault_int_faces = domain.GetFaultInteriorFaces();
      int nf_int = fault_int_faces.Size();

      std::vector<int> face_elem1(nf_int), face_elem2(nf_int);
      for (int i = 0; i < nf_int; i++)
      {
         FaceElementTransformations *FTr =
            pmesh.GetInteriorFaceTransformations(fault_int_faces[i]);
         face_elem1[i] = FTr->Elem1No;
         face_elem2[i] = FTr->Elem2No;
      }

      // Extract fault quantities from state
      Vector slip, theta;
      fault_op.GetSlip(state, slip);
      fault_op.GetTheta(state, theta);
      const Vector &V_rate = fault_op.GetSlipRate();
      const Vector &traction = seas_op.GetTraction();
      const Vector &tau_pre = fault_geom.GetTauPre();
      const Vector &a_vals = fault_geom.GetAValues();
      const Vector &dc_vals = fault_geom.GetDcValues();

      // Create L2 p=0 fields for state-based quantities
      ParGridFunction slip_dip_f(&l2_fes);    slip_dip_f = 0.0;
      ParGridFunction slip_strike_f(&l2_fes); slip_strike_f = 0.0;
      ParGridFunction V_dip_f(&l2_fes);       V_dip_f = 0.0;
      ParGridFunction V_strike_f(&l2_fes);    V_strike_f = 0.0;
      ParGridFunction V_mag_f(&l2_fes);       V_mag_f = 0.0;
      ParGridFunction tau_dip_f(&l2_fes);     tau_dip_f = 0.0;
      ParGridFunction tau_strike_f(&l2_fes);  tau_strike_f = 0.0;
      ParGridFunction tau_mag_f(&l2_fes);     tau_mag_f = 0.0;
      ParGridFunction psi_f(&l2_fes);         psi_f = 0.0;
      ParGridFunction a_state_f(&l2_fes);     a_state_f = 0.0;
      ParGridFunction dc_f(&l2_fes);          dc_f = 0.0;

      // Map fault DOFs to adjacent elements
      // DOF ordering: interior faces first (0..nf_int-1), then shared faces
      for (int i = 0; i < nf_int; i++)
      {
         int e1 = face_elem1[i];
         int e2 = face_elem2[i];

         real_t sd = slip(2 * i);
         real_t ss = slip(2 * i + 1);
         slip_dip_f(e1) = sd;    slip_dip_f(e2) = sd;
         slip_strike_f(e1) = ss; slip_strike_f(e2) = ss;

         real_t vd = V_rate(2 * i);
         real_t vs = V_rate(2 * i + 1);
         real_t vm = std::sqrt(vd * vd + vs * vs);
         V_dip_f(e1) = vd;    V_dip_f(e2) = vd;
         V_strike_f(e1) = vs; V_strike_f(e2) = vs;
         V_mag_f(e1) = vm;    V_mag_f(e2) = vm;

         // Total stress = pre-stress + elastic traction
         real_t td = tau_pre(2 * i) + traction(2 * i);
         real_t ts = tau_pre(2 * i + 1) + traction(2 * i + 1);
         real_t tm = std::sqrt(td * td + ts * ts);
         tau_dip_f(e1) = td;    tau_dip_f(e2) = td;
         tau_strike_f(e1) = ts; tau_strike_f(e2) = ts;
         tau_mag_f(e1) = tm;    tau_mag_f(e2) = tm;

         psi_f(e1) = theta(i); psi_f(e2) = theta(i);
         a_state_f(e1) = a_vals(i); a_state_f(e2) = a_vals(i);
         dc_f(e1) = dc_vals(i); dc_f(e2) = dc_vals(i);
      }

      if (mpi.IsRoot())
      {
         std::cout << "  Fault interior faces mapped: " << nf_int << "\n";
      }

      // Save via ParaViewDataCollection
      ParaViewDataCollection pv("bp5_diag", &pmesh);
      pv.SetPrefixPath(output_dir);
      pv.SetDataFormat(VTKFormat::ASCII);
      // Volume displacement
      pv.RegisterField("displacement", &u_diag);
      // Centroid-based analytic fields
      pv.RegisterField("fault_a", &a_field);
      pv.RegisterField("tau_pre_magnitude", &tau_pre_mag);
      pv.RegisterField("V_init_magnitude", &V_init_field);
      // State-based fields
      pv.RegisterField("slip_dip", &slip_dip_f);
      pv.RegisterField("slip_strike", &slip_strike_f);
      pv.RegisterField("V_dip", &V_dip_f);
      pv.RegisterField("V_strike", &V_strike_f);
      pv.RegisterField("V_magnitude", &V_mag_f);
      pv.RegisterField("tau_dip", &tau_dip_f);
      pv.RegisterField("tau_strike", &tau_strike_f);
      pv.RegisterField("tau_magnitude", &tau_mag_f);
      pv.RegisterField("state_psi", &psi_f);
      pv.RegisterField("fault_a_state", &a_state_f);
      pv.RegisterField("fault_dc", &dc_f);
      pv.SetCycle(0);
      pv.SetTime(0.0);
      pv.Save();

      // Also output boundary attributes
      pmesh.PrintBdrVTU(output_dir + "/boundary_attributes");

      if (mpi.IsRoot())
      {
         std::cout << "  ParaView output: " << output_dir << "/bp5_diag/\n";
         std::cout << "  Fields: displacement, fault_a, tau_pre_magnitude, "
                   << "V_init_magnitude\n";
         std::cout << "  State fields: slip_dip/strike, V_dip/strike/magnitude, "
                   << "tau_dip/strike/magnitude, state_psi, fault_a_state, "
                   << "fault_dc\n";
         std::cout << "  Boundary VTK: " << output_dir
                   << "/boundary_attributes\n";
         std::cout << "  Open in ParaView: File > Open > bp5_diag.pvd\n";
         std::cout << "  To see fault 'a': Threshold filter on fault_a > 0\n";
         std::cout << "=== Diagnostic VTK Complete ===\n\n";
      }
   }

   auto write_station_fault_diagnostics = [&](real_t time_now)
   {
      if (!diag_station_traction_decomp && !diag_station_jump_residual)
      {
         return;
      }

      Vector slip_diag;
      fault_op.GetSlip(state, slip_diag);

      Vector traction_diag, traction_stress_diag, traction_corr_diag;
      Vector jump_residual_diag;
      if (diag_station_jump_residual)
      {
         domain.ComputeTractionDiagnostics(seas_op.GetDisplacement(),
                                           slip_diag,
                                           traction_diag,
                                           traction_stress_diag,
                                           traction_corr_diag,
                                           jump_residual_diag);
      }
      else
      {
         domain.ComputeTractionComponents(seas_op.GetDisplacement(),
                                          slip_diag,
                                          traction_diag,
                                          traction_stress_diag,
                                          traction_corr_diag);
      }

      if (diag_station_traction_decomp)
      {
         bench_out.WriteTractionDecomposition(time_now,
                                              traction_stress_diag,
                                              traction_corr_diag);
      }
      if (diag_station_jump_residual)
      {
         bench_out.WriteJumpResidual(time_now, jump_residual_diag);
      }
   };

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
      if (diag_rk_stages) { ode_solver.SetDiagRKStages(true); }

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
   // v58: One-step tip reproducer — isolate first-solve asymmetry
   // =========================================================================
   if (diag_tip_step1)
   {
      if (mpi.IsRoot())
      {
         std::cout << "\n[TIP-STEP1] One-step reproducer: applying first nonzero slip...\n";
      }

      // Apply first nonzero slip: V_init * dt with dt=0.01s (Tandem default)
      // This mimics what the first RK stage does.
      real_t dt_step1 = 0.01;
      int num_nodes = fault_op.NumNodes();
      const int spn = 3;  // BP5: [slip_dip, slip_strike, psi]

      // Create a state with first-step slip increment
      Vector state_step1 = state;  // copy initial state
      for (int i = 0; i < num_nodes; i++)
      {
         // dslip/dt = V_init, so slip += V_init * dt * RK_coeff
         // For RK45 stage 2: coeff = a21 = 1/5
         real_t coeff = 0.2;
         const Vector &V = fault_op.GetSlipRate();
         state_step1(i * spn + 0) += V(2*i)   * dt_step1 * coeff;
         state_step1(i * spn + 1) += V(2*i+1) * dt_step1 * coeff;
      }

      // Extract slip from state
      Vector slip_step1(fault_op.SlipSize());
      fault_op.GetSlip(state_step1, slip_step1);

      // Expand to local and solve
      Vector local_slip_step1;
      domain.ExpandOwnedToLocalFault(slip_step1, local_slip_step1,
                                      domain.NumSlipComponents());

      // === Solve K*u = b(slip) ===
      ParGridFunction u_step1(&domain.GetFESpace());
      domain.Solve(0.0, local_slip_step1, u_step1);

      // === Compute traction with decomposition ===
      Vector trac_step1, trac_stress_step1, trac_corr_step1;
      domain.ComputeTractionComponents(u_step1, local_slip_step1,
                                        trac_step1, nullptr,
                                        &trac_stress_step1,
                                        &trac_corr_step1, nullptr);

      // === Dump mirror tip faces ===
      const auto *geom = fault_op.GetGeometry();
      const Vector &x2 = geom->GetCoordsX2();
      const Vector &x3 = geom->GetCoordsX3();

      for (int i = 0; i < num_nodes; i++)
      {
         // Mirror tip DOFs: |x2| > 49km AND depth < 2.5km
         if (std::abs(x2(i)) > 49000.0 && x3(i) < 2500.0)
         {
            real_t tau_stress_d = trac_stress_step1(2*i);
            real_t tau_stress_s = trac_stress_step1(2*i+1);
            real_t tau_corr_d = trac_corr_step1(2*i);
            real_t tau_corr_s = trac_corr_step1(2*i+1);
            real_t tau_d = trac_step1(2*i);
            real_t tau_s = trac_step1(2*i+1);
            real_t slip_d = slip_step1(2*i);
            real_t slip_s = slip_step1(2*i+1);

            mfem::out << std::scientific << std::setprecision(8)
               << "[TIP-STEP1] rank=" << mpi.Rank()
               << " dof=" << i
               << " x2=" << x2(i) << " x3=" << x3(i)
               << "\n  slip=(" << slip_d << "," << slip_s << ")"
               << "\n  tau_total=(" << tau_d << "," << tau_s << ")"
               << "  |tau|=" << std::sqrt(tau_d*tau_d + tau_s*tau_s)
               << "\n  tau_stress=(" << tau_stress_d << "," << tau_stress_s << ")"
               << "  |stress|=" << std::sqrt(tau_stress_d*tau_stress_d +
                                              tau_stress_s*tau_stress_s)
               << "\n  tau_corr=(" << tau_corr_d << "," << tau_corr_s << ")"
               << "  |corr|=" << std::sqrt(tau_corr_d*tau_corr_d +
                                            tau_corr_s*tau_corr_s)
               << "\n";
         }
      }

      if (mpi.IsRoot())
      {
         std::cout << "[TIP-STEP1] Done. Exiting.\n";
      }
      MPI_Finalize();
      return 0;
   }

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
               if (diag_psi_clamp && mpi.IsRoot())
               {
                  std::cout << "  [psi-clamp] first activation at accepted step "
                            << psi_clamp_first_step
                            << ", t=" << std::scientific
                            << std::setprecision(6) << psi_clamp_first_time
                            << " s, hi_hits=" << global_hi
                            << ", lo_hits=" << global_lo;
                  if (global_hi > 0)
                  {
                     std::cout << ", max_hi=" << global_max_hi;
                  }
                  if (global_lo > 0)
                  {
                     std::cout << ", min_lo=" << global_min_lo;
                  }
                  std::cout << "\n";
               }
            }
         }
      }

      real_t V_max = seas_op.GetMaxSlipRate();

      // v58 tip DOF monitor: sparse-to-dense schedule
      // Early: powers of 2; mid-run: every 1000; dense near failure: every 10
      {
         real_t t_yr = t / BP5Params::seconds_per_year;
         bool should_log = false;
         if (step <= 64 && (step & (step - 1)) == 0) { should_log = true; }  // powers of 2
         else if (step % 1000 == 0) { should_log = true; }                    // every 1000
         else if (t_yr > 0.45 && step % 10 == 0) { should_log = true; }      // dense near failure

         if (should_log)
         {
            const auto *geom = fault_op.GetGeometry();
            if (geom)
            {
               const Vector &x2 = geom->GetCoordsX2();
               const Vector &x3 = geom->GetCoordsX3();
               const Vector &slip_rate = fault_op.GetSlipRate();
               const Vector &traction = seas_op.GetTraction();
               int nn = fault_op.NumNodes();
               int spn = 3;  // BP5: slip_dip, slip_strike, psi

               for (int i = 0; i < nn; i++)
               {
                  // Monitor exact crash DOFs: |x2| > 48km AND x3 < 2.5km
                  if (std::abs(x2(i)) > 48000.0 && x3(i) < 2500.0)
                  {
                     real_t psi = state(i * spn + 2);
                     real_t slip_d = state(i * spn + 0);
                     real_t slip_s = state(i * spn + 1);
                     real_t V_d = slip_rate(2*i);
                     real_t V_s = slip_rate(2*i+1);
                     real_t V_abs = std::sqrt(V_d*V_d + V_s*V_s);
                     real_t tau_d = traction(2*i);
                     real_t tau_s = traction(2*i+1);
                     real_t tau_abs = std::sqrt(tau_d*tau_d + tau_s*tau_s);
                     mfem::out << "[TIP-MON] step=" << step
                               << " t_yr=" << std::scientific << std::setprecision(6)
                               << t_yr
                               << " dt=" << dt
                               << " rank=" << mpi.Rank()
                               << " dof=" << i
                               << " x2=" << x2(i)
                               << " x3=" << x3(i)
                               << " psi=" << psi
                               << " |V|=" << V_abs
                               << " |tau|=" << tau_abs
                               << " slip=(" << slip_d << "," << slip_s << ")"
                               << "\n";
                  }
               }
            }
         }
      }

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
         write_station_fault_diagnostics(t);
         bench_out.Flush();
      }
      else if (bench_out.Write(t, state, fault_op, seas_op.GetTraction(),
                               V_max))
      {
         write_station_fault_diagnostics(t);
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
