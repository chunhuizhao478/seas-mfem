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

// BP2 Verification: First Earthquake Cycle
//
// Runs a parallel BP2-QD simulation at ~200m resolution, outputs SCEC-format
// time series at 12 probe depths, then compares against Erickson's 50m
// finite-difference reference solution for the first earthquake cycle (~209 yr).
//
// Usage:
//   mpirun -np 10 ./seas_bp2_verify [options]
//
// Options:
//   --ref-dir DIR              Reference data directory (default: bp2/benchmark_data_200m)
//   --output-dir DIR           Output directory for simulation files (default: .)
//   --comparison-only          Skip simulation, only compare existing output vs reference
//   --tfinal T                 Override final time in seconds
//   --checkpoint-interval N    Write checkpoint every N steps (default: 5000, 0=off)
//   --restart PREFIX           Restart from checkpoint files PREFIX_checkpoint_r*.txt

#include "mfem.hpp"
#include "../../solver/seas_operator.hpp"
#include "../../solver/time_stepper.hpp"
#include "../../domain/antiplane_operator.hpp"
#include "../../fault/fault_geometry.hpp"
#include "../../fault/rate_state_fault.hpp"
#include "../../friction/dieterich_ruina.hpp"
#include "../../friction/state_evolution.hpp"
#include "../../config/bp2_params.hpp"
#include "../../domain/bp2_mesh.hpp"
#include "../../io/parallel_benchmark_output.hpp"
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
#include <numeric>

using namespace mfem;
using namespace mfem::seas;

// ============================================================================
// Reference data loading and comparison utilities
// ============================================================================

struct TimeSeriesData
{
   std::vector<double> time;
   std::vector<double> slip;
   std::vector<double> log10_slip_rate;
   std::vector<double> shear_stress;
   std::vector<double> log10_state;

   int Size() const { return static_cast<int>(time.size()); }
};

/// Load SCEC-format time series file (5 columns, # comment lines).
bool LoadTimeSeriesFile(const std::string &filename, TimeSeriesData &data,
                        double t_max = 1e30)
{
   std::ifstream file(filename);
   if (!file.is_open()) { return false; }

   std::string line;
   while (std::getline(file, line))
   {
      // Skip comment lines and header lines
      if (line.empty() || line[0] == '#') { continue; }
      // Skip the header line with column names
      if (line.find("slip") != std::string::npos) { continue; }

      std::istringstream iss(line);
      double t, s, vlog, tau, thlog;
      if (!(iss >> t >> s >> vlog >> tau >> thlog)) { continue; }

      if (t > t_max) { break; }

      data.time.push_back(t);
      data.slip.push_back(s);
      data.log10_slip_rate.push_back(vlog);
      data.shear_stress.push_back(tau);
      data.log10_state.push_back(thlog);
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

      // Clamp to reference range
      if (xt <= x_ref.front())
      {
         y_interp[i] = y_ref.front();
         continue;
      }
      if (xt >= x_ref.back())
      {
         y_interp[i] = y_ref.back();
         continue;
      }

      // Advance j so that x_ref[j] <= xt < x_ref[j+1]
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

/// Generate SCEC-convention filename for a given depth.
std::string MakeFilename(const std::string &prefix, double depth_m)
{
   double depth_km = std::abs(depth_m) / 1000.0;
   std::ostringstream oss;
   oss << prefix << "_z";
   if (std::abs(depth_km - std::round(depth_km)) < 1e-6)
   {
      oss << static_cast<int>(std::round(depth_km));
   }
   else
   {
      oss << std::fixed << std::setprecision(1) << depth_km;
   }
   oss << "km.txt";
   return oss.str();
}

/// Generate Erickson reference filename for a given depth.
std::string MakeRefFilename(const std::string &ref_dir, double depth_m)
{
   double depth_km = std::abs(depth_m) / 1000.0;
   std::ostringstream oss;
   oss << ref_dir << "/bp2-qd-erickson-z";
   if (std::abs(depth_km - std::round(depth_km)) < 1e-6)
   {
      oss << static_cast<int>(std::round(depth_km));
   }
   else
   {
      oss << std::fixed << std::setprecision(1) << depth_km;
   }
   oss << "km-res.txt";
   return oss.str();
}

struct ProbeComparison
{
   double depth_km;
   double slip_l2_error;
   double slip_rate_l2_error;
   double stress_l2_error;
   double state_l2_error;
   double nucleation_time_sim_yr;
   double nucleation_time_ref_yr;
   bool sim_loaded;
   bool ref_loaded;
};

/// Run comparison for all 12 probe depths.
/// Returns true if comparison passes basic sanity checks.
bool RunComparison(const std::string &output_dir,
                   const std::string &output_prefix,
                   const std::string &ref_dir,
                   const std::vector<double> &probe_depths_m,
                   double t_final_s)
{
   const double seconds_per_year = 365.25 * 24.0 * 3600.0;
   std::vector<ProbeComparison> results;

   std::cout << "\n" << std::string(80, '=') << "\n";
   std::cout << "BP2 First Cycle Verification: Comparison Summary\n";
   std::cout << std::string(80, '=') << "\n\n";

   for (double depth_m : probe_depths_m)
   {
      ProbeComparison pc;
      pc.depth_km = std::abs(depth_m) / 1000.0;
      pc.sim_loaded = false;
      pc.ref_loaded = false;
      pc.nucleation_time_sim_yr = -1.0;
      pc.nucleation_time_ref_yr = -1.0;

      // Build filenames
      std::string sim_file = output_dir + "/" +
         MakeFilename(output_prefix, depth_m);
      std::string ref_file = MakeRefFilename(ref_dir, depth_m);

      // Load data
      TimeSeriesData sim_data, ref_data;
      pc.sim_loaded = LoadTimeSeriesFile(sim_file, sim_data, t_final_s);
      pc.ref_loaded = LoadTimeSeriesFile(ref_file, ref_data, t_final_s);

      if (!pc.sim_loaded)
      {
         std::cerr << "  WARNING: Could not load simulation data: "
                   << sim_file << "\n";
      }
      if (!pc.ref_loaded)
      {
         std::cerr << "  WARNING: Could not load reference data: "
                   << ref_file << "\n";
      }

      if (pc.sim_loaded && pc.ref_loaded)
      {
         // Interpolate reference onto simulation time grid
         auto ref_slip = InterpolateOnto(ref_data.time, ref_data.slip,
                                          sim_data.time);
         auto ref_vlog = InterpolateOnto(ref_data.time,
                                          ref_data.log10_slip_rate,
                                          sim_data.time);
         auto ref_tau = InterpolateOnto(ref_data.time,
                                         ref_data.shear_stress,
                                         sim_data.time);
         auto ref_thlog = InterpolateOnto(ref_data.time,
                                           ref_data.log10_state,
                                           sim_data.time);

         pc.slip_l2_error = RelativeL2Error(sim_data.slip, ref_slip);
         pc.slip_rate_l2_error = RelativeL2Error(sim_data.log10_slip_rate,
                                                  ref_vlog);
         pc.stress_l2_error = RelativeL2Error(sim_data.shear_stress, ref_tau);
         pc.state_l2_error = RelativeL2Error(sim_data.log10_state, ref_thlog);

         // Extract nucleation time (first time log10(V) > -3, i.e. V > 1e-3)
         for (int i = 0; i < sim_data.Size(); i++)
         {
            if (sim_data.log10_slip_rate[i] > -3.0)
            {
               pc.nucleation_time_sim_yr = sim_data.time[i] / seconds_per_year;
               break;
            }
         }
         for (size_t i = 0; i < ref_data.time.size(); i++)
         {
            if (ref_data.log10_slip_rate[i] > -3.0)
            {
               pc.nucleation_time_ref_yr = ref_data.time[i] / seconds_per_year;
               break;
            }
         }
      }
      else
      {
         pc.slip_l2_error = -1.0;
         pc.slip_rate_l2_error = -1.0;
         pc.stress_l2_error = -1.0;
         pc.state_l2_error = -1.0;
      }

      results.push_back(pc);
   }

   // Print summary table
   std::cout << std::setw(10) << "Depth(km)"
             << std::setw(14) << "Slip L2"
             << std::setw(14) << "SlipRate L2"
             << std::setw(14) << "Stress L2"
             << std::setw(14) << "State L2"
             << std::setw(14) << "Tnuc_sim(yr)"
             << std::setw(14) << "Tnuc_ref(yr)"
             << "\n";
   std::cout << std::string(94, '-') << "\n";

   for (const auto &pc : results)
   {
      std::cout << std::setw(10) << std::fixed << std::setprecision(1)
                << pc.depth_km;
      if (pc.sim_loaded && pc.ref_loaded)
      {
         std::cout << std::setw(14) << std::scientific << std::setprecision(3)
                   << pc.slip_l2_error
                   << std::setw(14) << pc.slip_rate_l2_error
                   << std::setw(14) << pc.stress_l2_error
                   << std::setw(14) << pc.state_l2_error;
         if (pc.nucleation_time_sim_yr > 0)
         {
            std::cout << std::setw(14) << std::fixed << std::setprecision(1)
                      << pc.nucleation_time_sim_yr;
         }
         else
         {
            std::cout << std::setw(14) << "N/A";
         }
         if (pc.nucleation_time_ref_yr > 0)
         {
            std::cout << std::setw(14) << std::fixed << std::setprecision(1)
                      << pc.nucleation_time_ref_yr;
         }
         else
         {
            std::cout << std::setw(14) << "N/A";
         }
      }
      else
      {
         std::cout << std::setw(14) << "MISSING"
                   << std::setw(14) << "MISSING"
                   << std::setw(14) << "MISSING"
                   << std::setw(14) << "MISSING"
                   << std::setw(14) << "N/A"
                   << std::setw(14) << "N/A";
      }
      std::cout << "\n";
   }

   // Write CSV summary
   std::string csv_file = output_dir + "/bp2_verify_comparison.csv";
   std::ofstream csv(csv_file);
   if (csv.is_open())
   {
      csv << "depth_km,slip_l2,slip_rate_l2,stress_l2,state_l2,"
          << "nucleation_time_sim_yr,nucleation_time_ref_yr\n";
      for (const auto &pc : results)
      {
         csv << std::fixed << std::setprecision(1) << pc.depth_km << ","
             << std::scientific << std::setprecision(6)
             << pc.slip_l2_error << ","
             << pc.slip_rate_l2_error << ","
             << pc.stress_l2_error << ","
             << pc.state_l2_error << ","
             << std::fixed << std::setprecision(2)
             << pc.nucleation_time_sim_yr << ","
             << pc.nucleation_time_ref_yr << "\n";
      }
      csv.close();
      std::cout << "\nCSV summary written to: " << csv_file << "\n";
   }

   // Check pass/fail criteria
   bool pass = true;
   int num_with_data = 0;
   int num_earthquakes_detected = 0;
   bool has_nan = false;

   for (const auto &pc : results)
   {
      if (!pc.sim_loaded || !pc.ref_loaded) { continue; }
      num_with_data++;

      if (std::isnan(pc.slip_l2_error) || std::isinf(pc.slip_l2_error) ||
          std::isnan(pc.stress_l2_error) || std::isinf(pc.stress_l2_error))
      {
         has_nan = true;
      }

      if (pc.nucleation_time_sim_yr > 0) { num_earthquakes_detected++; }
   }

   // Check nucleation time relative error for probes that detected earthquakes
   for (const auto &pc : results)
   {
      if (pc.nucleation_time_sim_yr > 0 && pc.nucleation_time_ref_yr > 0)
      {
         double rel_err = std::abs(pc.nucleation_time_sim_yr -
                                    pc.nucleation_time_ref_yr) /
                          pc.nucleation_time_ref_yr;
         if (rel_err > 0.10)
         {
            std::cout << "\n  WARNING: Nucleation time error " << rel_err * 100
                      << "% exceeds 10% threshold at z="
                      << pc.depth_km << "km\n";
         }
      }
   }

   std::cout << "\n=== Verification Summary ===\n";
   std::cout << "  Probes with data: " << num_with_data << " / "
             << results.size() << "\n";
   std::cout << "  Probes detecting earthquake: " << num_earthquakes_detected
             << "\n";
   std::cout << "  NaN/Inf detected: " << (has_nan ? "YES" : "NO") << "\n";

   if (num_with_data == 0) { pass = false; }
   if (num_earthquakes_detected == 0) { pass = false; }
   if (has_nan) { pass = false; }

   std::cout << "  Overall: " << (pass ? "PASS" : "FAIL") << "\n";
   std::cout << std::string(80, '=') << "\n";

   return pass;
}

// ============================================================================
// Main driver
// ============================================================================

int main(int argc, char *argv[])
{
   // Initialize MPI + Hypre via MPIContext
   MPIContext mpi(&argc, &argv);

   // Parse command-line arguments
   std::string ref_dir = "bp2/benchmark_data_200m";
   std::string output_dir = ".";
   std::string output_prefix = "bp2_verify";
   bool comparison_only = false;
   double tfinal_override = 0.0;
   int checkpoint_interval = 5000;   // steps between checkpoints (0 = disabled)
   std::string restart_prefix;       // non-empty = restart from checkpoint

   for (int i = 1; i < argc; i++)
   {
      std::string arg(argv[i]);
      if (arg == "--ref-dir" && i + 1 < argc) { ref_dir = argv[++i]; }
      if (arg == "--output-dir" && i + 1 < argc) { output_dir = argv[++i]; }
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
   }

   // Standard 12 probe depths [m]
   std::vector<double> probe_depths_m = {
      0.0, -2400.0, -4800.0, -7200.0, -9600.0,
      -12000.0, -14400.0, -16800.0, -19200.0,
      -24000.0, -28800.0, -36000.0
   };

   // 250 years in seconds
   double t_final = 250.0 * BP2Params::seconds_per_year;
   if (tfinal_override > 0.0) { t_final = tfinal_override; }

   // Prepend output_dir to prefix
   std::string full_prefix = output_dir + "/" + output_prefix;

   // =========================================================================
   // Comparison-only mode: skip simulation
   // =========================================================================
   if (comparison_only)
   {
      if (mpi.IsRoot())
      {
         RunComparison(output_dir, output_prefix, ref_dir,
                       probe_depths_m, t_final);
      }
      return 0;
   }

   // =========================================================================
   // Simulation parameters: ~200m resolution matching Tandem configuration
   // =========================================================================
   BP2Params params;
   params.t_final = t_final;

   int mesh_nx = 25;
   int mesh_nz = 250;
   real_t grading_x = 7.0;
   real_t grading_z = 4.0;
   real_t Lx = 400.0e3;   // 400 km
   real_t Lz = 400.0e3;   // 400 km

   if (mpi.IsRoot())
   {
      std::cout << "BP2 Verification: First Earthquake Cycle\n";
      std::cout << "========================================\n";
      std::cout << "  Ranks: " << mpi.Size() << "\n";
      std::cout << "  Resolution: ~200m (nx=" << mesh_nx
                << ", nz=" << mesh_nz << ")\n";
      std::cout << "  Domain: Lx=" << Lx / 1e3 << "km, Lz="
                << Lz / 1e3 << "km\n";
      std::cout << "  Grading: x=" << grading_x << ", z=" << grading_z << "\n";
      std::cout << "  t_final: " << t_final / BP2Params::seconds_per_year
                << " years\n";
      std::cout << "  Reference dir: " << ref_dir << "\n";
      std::cout << "  Output prefix: " << full_prefix << "\n\n";
   }

   // =========================================================================
   // Mesh: create serial mesh on all ranks, then distribute
   // =========================================================================
   BP2MeshGenerator::Parameters mesh_params;
   mesh_params.Lx = Lx;
   mesh_params.Lz = Lz;
   mesh_params.Wf = params.Wf;
   mesh_params.nx = mesh_nx;
   mesh_params.nz = mesh_nz;
   mesh_params.grading_x = grading_x;
   mesh_params.grading_z = grading_z;

   auto serial_mesh = BP2MeshGenerator::CreateGraded(mesh_params);
   ParMesh pmesh(mpi.GetComm(), *serial_mesh);
   serial_mesh.reset();

   // GetGlobalNE() is collective (MPI_Allreduce) — call on all ranks
   long long global_ne = pmesh.GetGlobalNE();
   if (mpi.IsRoot())
   {
      std::cout << "  ParMesh: " << global_ne << " global elements\n";
   }

   // =========================================================================
   // Domain operator: DG order 1
   // =========================================================================
   int order = 1;
   AntiplaneDomainOperator<ParMesh> domain(
      pmesh, order, params.mu(), params.Vp, params.Wf, DGMethod::BR2);

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
   fc.V0 = params.V0;  fc.f0 = params.f0;
   fc.b = params.b;     fc.Dc = params.Dc;
   DieterichRuinaFriction friction(fc);
   AgingLaw aging;

   RateStateFaultOperator<ParMesh> fault_op(
      &fault_geom, &friction, &aging, params, &mpi);

   // =========================================================================
   // SEAS operator
   // =========================================================================
   SEASQuasiDynamicOperator<ParMesh> seas_op(&domain, &fault_op, &mpi);

   Vector state(fault_op.StateSize());
   seas_op.SetInitialCondition(state);

   real_t V_init = seas_op.GetMaxSlipRate();
   if (mpi.IsRoot())
   {
      std::cout << "\nInitial conditions:\n";
      std::cout << "  tau0 = " << fault_op.GetTau0() / 1e6 << " MPa\n";
      std::cout << "  V_max = " << V_init << " m/s\n";
   }

   // =========================================================================
   // I/O: gather global fault depths, create parallel benchmark output
   // =========================================================================
   Vector local_fault_depths;
   domain.GetFaultDepths(local_fault_depths);

   Vector dedup_fault_depths;
   {
      Vector dedup_data;
      fault_geom.GatherToRootDedup(local_fault_depths, dedup_data,
                                    dedup_fault_depths);
   }

   if (mpi.IsRoot())
   {
      std::cout << "  Deduplicated global fault DOFs: "
                << dedup_fault_depths.Size() << "\n";
   }

   // Clean up old output files
   if (mpi.IsRoot())
   {
      std::string rm_cmd = "rm -f " + full_prefix + "_z*.txt";
      std::system(rm_cmd.c_str());
   }
   mpi.Barrier();

   std::vector<real_t> probe_depths_rt(probe_depths_m.begin(),
                                        probe_depths_m.end());
   ParallelBenchmarkOutput bench_out(
      full_prefix, params, probe_depths_rt, fault_geom, mpi,
      dedup_fault_depths);

   // Write initial state
   bench_out.ForceWrite(0.0, state, fault_op, seas_op.GetTraction(), V_init);

   // =========================================================================
   // Time integration (Dormand-Prince RK45)
   // =========================================================================
   DormandPrinceRK45 ode_solver;
   ode_solver.SetMPIContext(&mpi);
   ode_solver.SetAbsTol(1e-7);
   ode_solver.SetRelTol(1e-50);  // Match Tandem/PETSc: pure absolute tolerance
   ode_solver.SetDtMin(1e-6);
   ode_solver.SetDtMax(0.1 * BP2Params::seconds_per_year);
   ode_solver.SetDt(1e3);
   ode_solver.Init(seas_op);

   real_t t = 0.0;
   int step = 0;
   int max_steps = 10000000;

   // Earthquake detection
   bool in_seismic_event = false;
   int num_seismic_events = 0;
   real_t V_threshold_seismic = 1e-3;
   real_t V_threshold_interseismic = 1e-6;

   // Early termination: after first EQ resolves + 40 yr postseismic
   real_t first_eq_end_time = -1.0;
   real_t postseismic_duration = 40.0 * BP2Params::seconds_per_year;

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
      fault_op.InitPreStress();  // Must set tau0_ before ComputeRHS
      seas_op.SetDisplacement(restart_disp);
      fault_op.SetSlipRate(restart_slip_rate);
      if (restart_fsal && restart_k0.Size() > 0)
      {
         ode_solver.RestoreFSAL(restart_k0);
      }

      if (mpi.IsRoot())
      {
         std::cout << "\nRestarted from checkpoint:\n";
         std::cout << "  Time: " << t / BP2Params::seconds_per_year << " yr\n";
         std::cout << "  Step: " << step << "\n";
         std::cout << "  dt: " << restart_dt << " s\n";
         std::cout << "  Seismic events: " << num_seismic_events << "\n";
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

   real_t next_print_time = 0.0;
   real_t print_interval = 10.0 * BP2Params::seconds_per_year;

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

      real_t V_max = seas_op.GetMaxSlipRate();

      // Check for NaN/Inf
      if (std::isnan(V_max) || std::isinf(V_max))
      {
         if (mpi.IsRoot())
         {
            std::cerr << "NaN/Inf detected at step " << step
                      << ", t = " << t / BP2Params::seconds_per_year
                      << " yr\n";
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
                      << t / BP2Params::seconds_per_year << " yr"
                      << ", V_max = " << std::scientific
                      << std::setprecision(2) << V_max << " m/s ***\n";
         }
      }
      else if (in_seismic_event && V_max < V_threshold_interseismic)
      {
         in_seismic_event = false;
         if (num_seismic_events == 1 && first_eq_end_time < 0)
         {
            first_eq_end_time = t;
            if (mpi.IsRoot())
            {
               std::cout << "  First earthquake resolved at t = "
                         << std::fixed << std::setprecision(1)
                         << t / BP2Params::seconds_per_year << " yr\n";
               std::cout << "  Will continue for "
                         << postseismic_duration / BP2Params::seconds_per_year
                         << " yr postseismic\n";
            }
         }
      }

      // Early termination after first earthquake + postseismic
      if (first_eq_end_time > 0 && t > first_eq_end_time + postseismic_duration)
      {
         if (mpi.IsRoot())
         {
            std::cout << "  Early termination: postseismic period complete\n";
         }
         break;
      }

      // I/O
      if (bench_out.Write(t, state, fault_op, seas_op.GetTraction(), V_max))
      {
         bench_out.Flush();
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

      // Periodic console output
      if (mpi.IsRoot() && (t >= next_print_time || V_max > V_threshold_seismic))
      {
         std::cout << std::setw(10) << step
                   << std::setw(16) << std::fixed << std::setprecision(2)
                   << t / BP2Params::seconds_per_year
                   << std::setw(14) << std::scientific << std::setprecision(3)
                   << ode_solver.GetDt()
                   << std::setw(16) << std::scientific << std::setprecision(3)
                   << V_max
                   << std::setw(8) << num_seismic_events
                   << "\n";
         next_print_time = t + print_interval;
      }
   }

   // Final checkpoint
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

   // Final I/O
   bench_out.ForceWrite(t, state, fault_op, seas_op.GetTraction(),
                        seas_op.GetMaxSlipRate());
   bench_out.Close();

   if (mpi.IsRoot())
   {
      std::cout << "\n=== Simulation Summary ===\n";
      std::cout << "  Ranks: " << mpi.Size() << "\n";
      std::cout << "  Final time: " << t / BP2Params::seconds_per_year
                << " years\n";
      std::cout << "  Total steps: " << step << "\n";
      std::cout << "  Seismic events: " << num_seismic_events << "\n\n";
   }

   // =========================================================================
   // Post-simulation comparison (root only)
   // =========================================================================
   if (mpi.IsRoot())
   {
      RunComparison(output_dir, output_prefix, ref_dir,
                    probe_depths_m, t_final);
   }

   return 0;
}
