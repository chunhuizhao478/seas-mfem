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

#include "mfem.hpp"
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

/// Create a 3D hex mesh for BP5 with proper boundary attributes.
/// Domain: [-Lx,Lx] x [-Ly,Ly] x [0,Lz]
/// Boundary attributes: 1=x-, 2=x+, 3=y+, 4=y-, 5=z=0, 6=z=Lz
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
      if (std::abs(cx - (-Lx)) < tol)      { attr = 1; }
      else if (std::abs(cx - Lx) < tol)     { attr = 2; }
      else if (std::abs(cy - Ly) < tol)     { attr = 3; }
      else if (std::abs(cy - (-Ly)) < tol)  { attr = 4; }
      else if (std::abs(cz) < tol)          { attr = 5; }
      else if (std::abs(cz - Lz) < tol)     { attr = 6; }
      else                                   { attr = 1; }

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
   bool dump_bdr_vtk = false;
   bool diag_vtk = false;
   std::string bc_mode_str = "far-field";
   std::string psi_init_mode_str = "scec";  // "scec" or "tandem"

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
      if (arg == "--dump-bdr-vtk") { dump_bdr_vtk = true; }
      if (arg == "--diag-vtk") { diag_vtk = true; }
      if (arg == "--bc-mode" && i + 1 < argc) { bc_mode_str = argv[++i]; }
      if (arg == "--psi-init-mode" && i + 1 < argc)
      {
         psi_init_mode_str = argv[++i];
      }
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
   else if (solver_str == "cg" || solver_str == "CG")
   {
      solver_type = SolverType::CG_AMG;
   }

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
      return 1;
   }

   if (mpi.IsRoot())
   {
      std::cout << "  Ranks: " << mpi.Size() << "\n";
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
      std::cout << "  t_final: " << t_final / BP5Params::seconds_per_year
                << " years\n";
      std::cout << "  Output prefix: " << full_prefix << "\n";
      params.Print();
      std::cout << "\n";
   }

   ParMesh pmesh(mpi.GetComm(), *serial_mesh);
   serial_mesh.reset();

   long long global_ne = pmesh.GetGlobalNE();
   if (mpi.IsRoot())
   {
      std::cout << "  ParMesh: " << global_ne << " global elements\n";
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
   // Domain operator: 3D DG elasticity, order 1
   // =========================================================================
   int order = 1;
   ElasticityDomainOperator<ParMesh> domain(
      pmesh, order, params.lambda(), params.mu(),
      params.Vp, params.Wf, params.lf, dg_method, solver_type, bc_mode);

   if (check_residual) { domain.SetCheckResidual(true); }

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
   if (psi_init_mode_str == "tandem" || psi_init_mode_str == "Tandem")
   {
      fault_op.SetScecPsiInit(false);
      if (mpi.IsRoot())
      {
         std::cout << "  Psi init mode: Tandem (absorb delta_tau into psi)\n";
      }
   }
   else
   {
      if (mpi.IsRoot())
      {
         std::cout << "  Psi init mode: SCEC (genuine delta_tau overstress)\n";
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
   // I/O: gather global fault coords and tau_pre for parallel output
   // =========================================================================
   Vector local_x2, local_x3;
   domain.GetFaultCoords2D(local_x2, local_x3);

   // Gather coordinates and tau_pre to root
   Vector global_x2, global_x3;
   fault_geom.GatherToRoot(local_x2, global_x2);
   fault_geom.GatherToRoot(local_x3, global_x3);

   // Split and gather tau_pre components
   int N_local = fault_geom.NumLocalFaultDOFs();
   Vector local_tp_dip(N_local), local_tp_strike(N_local);
   const Vector &local_tau_pre = fault_geom.GetTauPre();
   for (int i = 0; i < N_local; i++)
   {
      local_tp_dip(i) = local_tau_pre(2 * i);
      local_tp_strike(i) = local_tau_pre(2 * i + 1);
   }

   Vector global_tp_dip, global_tp_strike;
   fault_geom.GatherToRoot(local_tp_dip, global_tp_dip);
   fault_geom.GatherToRoot(local_tp_strike, global_tp_strike);

   if (mpi.IsRoot())
   {
      std::cout << "  Global fault coords gathered: "
                << global_x2.Size() << " DOFs\n";
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
      global_x2, global_x3, global_tp_dip, global_tp_strike);

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
                  real_t d = std::abs(global_x2(j)) +
                             std::abs(global_x3(j) - target_z);
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

   // =========================================================================
   // Time integration (Dormand-Prince RK45)
   // =========================================================================
   DormandPrinceRK45 ode_solver;
   ode_solver.SetMPIContext(&mpi);
   ode_solver.SetAbsTol(1e-7);
   ode_solver.SetRelTol(1e-50);
   ode_solver.SetDtMin(1e-6);
   ode_solver.SetDtMax(0.1 * BP5Params::seconds_per_year);

   // Initial dt must be small enough for the nucleation zone dynamics.
   // The nucleation zone has V_nuc = 0.03 m/s with overstressing that
   // drives acceleration. Using a conservative multiplier (0.01) ensures
   // the RK45 adaptive controller can resolve the initial transient
   // without producing intermediate states that blow up. The controller
   // will quickly ramp up dt during the interseismic period.
   real_t dt_init = std::min(1e3, 0.01 * params.L_nuc /
                             std::max(V_init, 1e-20));
   ode_solver.SetDt(dt_init);
   if (mpi.IsRoot())
   {
      std::cout << "  Initial dt: " << dt_init << " s"
                << " (V_init_max = " << V_init << ")\n";
   }
   ode_solver.SetStatePerNode(3);  // BP5: [slip_dip, slip_strike, psi]
   ode_solver.Init(seas_op);

   real_t t = 0.0;
   int step = 0;
   int max_steps = 10000000;

   // Earthquake detection
   bool in_seismic_event = false;
   int num_seismic_events = 0;
   real_t V_threshold_seismic = 1e-3;
   real_t V_threshold_interseismic = 1e-6;

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

      ode_solver.SetDt(restart_dt);
      seas_op.SetDisplacement(restart_disp);
      fault_op.SetSlipRate(restart_slip_rate);
      if (restart_fsal && restart_k0.Size() > 0)
      {
         ode_solver.RestoreFSAL(restart_k0);
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
      if (t + ode_solver.GetDt() > t_final)
      {
         ode_solver.SetDt(t_final - t);
      }

      real_t dt;
      bool accepted = ode_solver.Step(seas_op, state, t, dt);
      if (!accepted) { continue; }
      step++;

      // Post-step psi clamping: prevent unphysical state variable values.
      // The explicit RK45 can overshoot psi during post-earthquake healing
      // (stiff exp((f0-psi)/b) term). Clamp psi to a physically reasonable
      // range to prevent the fault from getting trapped at V ≈ 0.
      {
         const int spn = 3;  // BP5: [slip_dip, slip_strike, psi]
         const int psi_idx = 2;
         // psi_max: steady-state at V = 1e-20 m/s with generous margin
         // psi_ss(1e-20) = f0 + b*ln(V0/1e-20) = 0.6 + 0.03*32.2 ≈ 1.57
         const real_t psi_max = 3.0;   // well above any physical steady state
         const real_t psi_min = -5.0;  // generous lower bound
         int n_nodes = state.Size() / spn;
         for (int i = 0; i < n_nodes; i++)
         {
            real_t &psi = state(i * spn + psi_idx);
            if (psi > psi_max) { psi = psi_max; }
            else if (psi < psi_min) { psi = psi_min; }
         }
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
         WriteCheckpoint(full_prefix, t, ode_solver.GetDt(),
                         step, num_seismic_events, in_seismic_event,
                         state, u_vec,
                         seas_op.GetTraction(), fault_op.GetSlipRate(),
                         ode_solver.IsInitialized(), ode_solver.GetK0(),
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
                   << ode_solver.GetDt()
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
      WriteCheckpoint(full_prefix, t, ode_solver.GetDt(),
                      step, num_seismic_events, in_seismic_event,
                      state, u_vec,
                      seas_op.GetTraction(), fault_op.GetSlipRate(),
                      ode_solver.IsInitialized(), ode_solver.GetK0(),
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
      std::cout << "  Seismic events: " << num_seismic_events << "\n\n";
   }

   // =========================================================================
   // Post-simulation comparison
   // =========================================================================
   if (mpi.IsRoot())
   {
      RunComparison(output_dir, output_prefix, ref_dir,
                    stations, t_final);
   }

   return 0;
}
