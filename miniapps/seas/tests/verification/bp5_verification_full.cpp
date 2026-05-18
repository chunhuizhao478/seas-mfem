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
//   --regression-tolerance TOL Regression test: fail if any per-field relative
//                              L2 error exceeds TOL (negative = skip check)
//   --ref-prefix PFX           Reference file prefix (auto-detected if omitted)

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
#include "../../io/petsc_ts_checkpoint.hpp"   // V2 PETSc-TS restart support
#include "../../io/paraview_output.hpp"
#include "../../io/hdf5_error_filter.hpp"   // Suppress dual-HDF5 noise on Frontera
#include "../../common/mpi_context.hpp"
#include "../../trace/face_trace_logger.hpp"
#include "../../config/bp5_mesh_utils.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <filesystem>   // weakly_canonical for the --restart / --output-dir safety check
#include <sstream>
#include <cmath>
#include <cctype>
#include <memory>
#include <vector>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <set>
#include <map>
#include <dirent.h>

using namespace mfem;
using namespace mfem::seas;

// CreateBP5InlineMesh() moved to config/bp5_mesh_utils.hpp

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
/// Precondition: x_ref and x_target must be sorted in ascending order.
std::vector<double> InterpolateOnto(const std::vector<double> &x_ref,
                                     const std::vector<double> &y_ref,
                                     const std::vector<double> &x_target)
{
   std::vector<double> y_interp(x_target.size(), 0.0);
   if (x_ref.empty() || y_ref.empty()) { return y_interp; }
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
   double norm_num = std::sqrt(num);
   double norm_den = std::sqrt(den);
   if (norm_den < 1e-30) { return (norm_num < 1e-30) ? 0.0 : 1e30; }
   return norm_num / norm_den;
}

/// Check if string `s` ends with `suffix`.
static bool EndsWith(const std::string &s, const std::string &suffix)
{
   return s.size() >= suffix.size() &&
          s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

/// Auto-detect reference file prefix by scanning ref_dir for *_global.txt
/// or *_fltst_strk+00dp+00.txt. Returns empty string if no match found.
/// Deterministic: sorts matches alphabetically, warns on ambiguity.
std::string DetectRefPrefix(const std::string &ref_dir)
{
   const std::string suffixes[] = {
      "_global.txt",
      "_fltst_strk+00dp+00.txt"
   };

   DIR *dir = opendir(ref_dir.c_str());
   if (!dir) { return ""; }

   // Collect all matching prefixes per suffix
   for (const auto &suffix : suffixes)
   {
      std::vector<std::string> prefixes;
      rewinddir(dir);
      struct dirent *entry;
      while ((entry = readdir(dir)) != nullptr)
      {
         std::string fname(entry->d_name);
         if (EndsWith(fname, suffix))
         {
            prefixes.push_back(fname.substr(0, fname.size() - suffix.size()));
         }
      }
      if (!prefixes.empty())
      {
         std::sort(prefixes.begin(), prefixes.end());
         if (prefixes.size() > 1)
         {
            std::cout << "WARNING: Multiple prefixes in " << ref_dir
                      << ", using '" << prefixes[0] << "'\n";
         }
         closedir(dir);
         return prefixes[0];
      }
   }
   closedir(dir);
   return "";
}

/// Extract a single field from BP5TimeSeriesData by column index (0-6).
/// Order: slip_strike, slip_dip, log10_V_strike, log10_V_dip,
///        tau_strike, tau_dip, log10_state
static const char *BP5FieldNames[7] = {
   "slip_strike", "slip_dip", "log10_V_strike", "log10_V_dip",
   "tau_strike", "tau_dip", "log10_state"
};

const std::vector<double> &GetBP5Field(const BP5TimeSeriesData &d, int idx)
{
   switch (idx)
   {
      case 0: return d.slip_strike;
      case 1: return d.slip_dip;
      case 2: return d.log10_V_strike;
      case 3: return d.log10_V_dip;
      case 4: return d.tau_strike;
      case 5: return d.tau_dip;
      case 6: return d.log10_state;
      default: MFEM_ABORT("Invalid BP5 field index"); return d.slip_strike;
   }
}

/// Load BP5 global output file (2 columns: time, log10(Vmax)).
bool LoadBP5GlobalFile(const std::string &filename,
                       std::vector<double> &time,
                       std::vector<double> &log10_vmax,
                       double t_max = 1e30)
{
   std::ifstream file(filename);
   if (!file.is_open()) { return false; }

   std::string line;
   while (std::getline(file, line))
   {
      if (line.empty() || line[0] == '#') { continue; }
      std::istringstream iss(line);
      double t, v;
      if (!(iss >> t >> v)) { continue; }
      if (t > t_max) { break; }
      time.push_back(t);
      log10_vmax.push_back(v);
   }
   return !time.empty();
}

/// Count seismic events: upward crossings of threshold in log10(Vmax).
int CountEvents(const std::vector<double> &log10_vmax,
                double threshold_log10 = -3.0)
{
   int count = 0;
   bool above = false;
   for (double v : log10_vmax)
   {
      if (v > threshold_log10 && !above)
      {
         count++;
         above = true;
      }
      else if (v <= threshold_log10)
      {
         above = false;
      }
   }
   return count;
}

/// Run comparison for all stations + global output.
/// Returns 0 on pass (all errors within tolerance), nonzero on fail.
/// If regression_tol < 0, comparison is informational only (always returns 0).
int RunComparison(const std::string &output_dir,
                  const std::string &output_prefix,
                  const std::string &ref_dir,
                  const std::string &ref_prefix,
                  const std::vector<Probe2DInterpolator::Station> &stations,
                  double t_final_s,
                  double regression_tol)
{
   std::cout << "\n" << std::string(80, '=') << "\n";
   std::cout << "BP5 Verification: Regression Comparison\n";
   std::cout << std::string(80, '=') << "\n";
   if (regression_tol >= 0.0)
   {
      std::cout << "  Regression tolerance: " << std::scientific
                << std::setprecision(2) << regression_tol << "\n";
   }
   else
   {
      std::cout << "  Mode: informational (no tolerance check)\n";
   }
   std::cout << "\n";

   bool any_fail = false;
   double max_error = 0.0;
   std::string worst_station, worst_field;
   int stations_compared = 0;

   // --- Station comparison (7 SCEC columns) ---
   for (const auto &st : stations)
   {
      std::string sim_file = output_dir + "/" + output_prefix + "_"
                             + st.name + ".txt";
      std::string ref_file = ref_dir + "/" + ref_prefix + "_"
                             + st.name + ".txt";

      BP5TimeSeriesData sim_data, ref_data;
      bool sim_ok = LoadBP5TimeSeriesFile(sim_file, sim_data, t_final_s);
      bool ref_ok = LoadBP5TimeSeriesFile(ref_file, ref_data, t_final_s);

      if (!sim_ok)
      {
         std::cout << "  Station " << std::setw(24) << st.name
                   << ": sim MISSING\n";
         if (regression_tol >= 0.0) { any_fail = true; }
         continue;
      }
      if (!ref_ok)
      {
         std::cout << "  Station " << std::setw(24) << st.name
                   << ": ref MISSING (" << ref_file << ")\n";
         continue;  // Missing ref is not a failure — may be informational
      }

      stations_compared++;
      std::cout << "  Station " << std::setw(24) << st.name
                << " (" << sim_data.Size() << " vs " << ref_data.Size()
                << " pts):\n";

      for (int f = 0; f < 7; f++)
      {
         // Interpolate reference onto simulation time grid
         const auto &ref_field = GetBP5Field(ref_data, f);
         auto ref_interp = InterpolateOnto(ref_data.time, ref_field,
                                           sim_data.time);
         const auto &sim_field = GetBP5Field(sim_data, f);
         double err = RelativeL2Error(sim_field, ref_interp);

         bool field_fail = (regression_tol >= 0.0 && err > regression_tol);
         const char *status = field_fail ? "FAIL" : "ok";

         std::cout << "    " << std::setw(16) << BP5FieldNames[f]
                   << ": L2_rel = " << std::scientific << std::setprecision(6)
                   << err << "  [" << status << "]\n";

         if (field_fail) { any_fail = true; }
         if (err > max_error)
         {
            max_error = err;
            worst_station = st.name;
            worst_field = BP5FieldNames[f];
         }
      }
   }

   // --- Global output comparison (log10(Vmax)) ---
   std::string sim_global = output_dir + "/" + output_prefix + "_global.txt";
   std::string ref_global = ref_dir + "/" + ref_prefix + "_global.txt";
   std::vector<double> sim_gt, sim_gv, ref_gt, ref_gv;
   bool sim_g_ok = LoadBP5GlobalFile(sim_global, sim_gt, sim_gv, t_final_s);
   bool ref_g_ok = LoadBP5GlobalFile(ref_global, ref_gt, ref_gv, t_final_s);

   if (sim_g_ok && ref_g_ok)
   {
      auto ref_gv_interp = InterpolateOnto(ref_gt, ref_gv, sim_gt);
      double err = RelativeL2Error(sim_gv, ref_gv_interp);
      bool g_fail = (regression_tol >= 0.0 && err > regression_tol);
      const char *status = g_fail ? "FAIL" : "ok";

      std::cout << "  Global log10(Vmax): L2_rel = " << std::scientific
                << std::setprecision(6) << err << "  [" << status << "]\n";
      if (g_fail) { any_fail = true; }
      if (err > max_error)
      {
         max_error = err;
         worst_station = "global";
         worst_field = "log10_Vmax";
      }

      // Event count comparison
      int sim_events = CountEvents(sim_gv);
      int ref_events = CountEvents(ref_gv);
      if (sim_events != ref_events)
      {
         std::cout << "  Event count: " << sim_events << " vs "
                   << ref_events << " [MISMATCH]\n";
         if (regression_tol >= 0.0) { any_fail = true; }
      }
      else
      {
         std::cout << "  Event count: " << sim_events << " [ok]\n";
      }
   }
   else
   {
      if (!sim_g_ok)
      {
         std::cout << "  Global: sim file missing (" << sim_global << ")\n";
      }
      if (!ref_g_ok)
      {
         std::cout << "  Global: ref file missing (" << ref_global << ")\n";
      }
   }

   // --- Summary ---
   std::cout << "\n" << std::string(80, '-') << "\n";
   std::cout << "  Stations compared: " << stations_compared << " / "
             << stations.size() << "\n";
   std::cout << "  Max relative L2 error: " << std::scientific
             << std::setprecision(6) << max_error;
   if (!worst_station.empty())
   {
      std::cout << " (" << worst_station << " / " << worst_field << ")";
   }
   std::cout << "\n";

   if (regression_tol >= 0.0)
   {
      if (stations_compared == 0)
      {
         std::cout << "  REGRESSION CHECK: FAIL (no stations compared — "
                      "check --ref-dir and --ref-prefix)\n";
         any_fail = true;
      }
      else if (any_fail)
      {
         std::cout << "  REGRESSION CHECK: FAIL (tolerance "
                   << regression_tol << " exceeded)\n";
      }
      else
      {
         std::cout << "  REGRESSION CHECK: PASS\n";
      }
   }
   std::cout << std::string(80, '=') << "\n";

   return any_fail ? 1 : 0;
}

// ============================================================================
// TSSolve monitor context and callback (Tandem-style time stepping)
// ============================================================================
#ifdef MFEM_USE_PETSC

/// Context passed to PETSc TSMonitor callback.  Holds references to all
/// objects needed for per-step monitoring (I/O, earthquake detection,
/// console output, checkpoint, face tracer).
struct BP5MonitorCtx
{
   // Core objects (non-owning)
   MPIContext *mpi;
   Vector *state;                                       // shares memory with PETSc Vec
   ParallelBP5BenchmarkOutput *bench_out;
   seas::FaceTraceLogger<ParMesh> *face_tracer;
   RateStateFaultOperator<ParMesh, 2> *fault_op;
   PBP5SEASOp *seas_op;
   ProbeOutput *global_out;                              // may be nullptr (non-root)

   // Settings
   bool write_every_step;
   int  print_step_interval;
   int  checkpoint_interval;
   std::string full_prefix;
   bool debug_first_step_dump;

   // Mutable per-step state
   int  num_seismic_events;
   bool in_seismic_event;
   real_t V_threshold_seismic;
   real_t V_threshold_interseismic;

   // Current dt (read from TS after each step)
   real_t current_dt;

   // ParaView output (may be nullptr if --paraview not set)
   std::function<void(int, real_t, real_t)> paraview_write_fn;

   // V2 PETSc-TS restart support (R-001 / R-005).
   // The monitor callback runs in a context where `pv_out`,
   // `use_petsc_ts`, `petsc_ode`, and `restart_rejections_carryover`
   // (all local to main) are OUT OF SCOPE.  Plumb them through here so
   // the monitor-site WritePetscTSCheckpoint snippet can read them
   // (R-001) and so the cumulative-rejection accumulator works across
   // restart chains (R-005).
   seas::ParaViewOutput<ParMesh> *pv_out = nullptr;     // may be nullptr
   int restart_rejections_carryover = 0;                 // R-303 / R-005
};

/// PETSc TSMonitor callback — called after every accepted step inside TSSolve.
/// Matches Tandem's MonitorQD::monitor pattern: read state, do I/O.
static PetscErrorCode bp5_ts_monitor_callback(
   TS ts, PetscInt step, PetscReal time, Vec u, void *ctx)
{
   PetscFunctionBeginUser;
   auto *mon = static_cast<BP5MonitorCtx*>(ctx);

   // state already reflects the accepted step (PlaceMemory shares memory).
   Vector &state  = *mon->state;
   auto   &mpi    = *mon->mpi;

   // Get dt for console output
   PetscReal dt;
   TSGetTimeStep(ts, &dt);
   mon->current_dt = dt;

   // Face tracer
   if (mon->face_tracer->IsActive())
   {
      mon->face_tracer->CommitStep(static_cast<int>(step), time, state, false);
   }

   // V_max & NaN check
   real_t V_max = mon->seas_op->GetMaxSlipRate();
   {
      bool has_nan = !std::isfinite(V_max);
      if (!has_nan)
      {
         for (int i = 0; i < std::min(state.Size(), 100); i++)
         {
            if (!std::isfinite(state(i))) { has_nan = true; break; }
         }
      }
      int local_nan = has_nan ? 1 : 0;
      int global_nan = mpi.GlobalSumInt(local_nan);
      has_nan = (global_nan > 0);
      if (has_nan)
      {
         if (mpi.IsRoot())
         {
            std::cerr << "NaN/Inf detected at step " << step
                      << ", t = " << time / BP5Params::seconds_per_year
                      << " yr, V_max = " << V_max << "\n";
         }
         TSSetConvergedReason(ts, TS_DIVERGED_STEP_REJECTED);
         PetscFunctionReturn(PETSC_SUCCESS);
      }
   }

   // Earthquake detection
   if (!mon->in_seismic_event && V_max > mon->V_threshold_seismic)
   {
      mon->in_seismic_event = true;
      mon->num_seismic_events++;
      if (mpi.IsRoot())
      {
         std::cout << "  *** EARTHQUAKE #" << mon->num_seismic_events
                   << " at t = " << std::fixed << std::setprecision(1)
                   << time / BP5Params::seconds_per_year << " yr"
                   << ", V_max = " << std::scientific
                   << std::setprecision(2) << V_max << " m/s ***\n";
      }
   }
   else if (mon->in_seismic_event && V_max < mon->V_threshold_interseismic)
   {
      mon->in_seismic_event = false;
      if (mpi.IsRoot())
      {
         std::cout << "  Earthquake #" << mon->num_seismic_events
                   << " resolved at t = "
                   << std::fixed << std::setprecision(1)
                   << time / BP5Params::seconds_per_year << " yr\n";
      }
   }

   // I/O: adaptive schedule or every step
   if (mon->write_every_step)
   {
      mon->bench_out->ForceWrite(time, state, *mon->fault_op,
                                 mon->seas_op->GetTraction(), V_max);
      mon->bench_out->Flush();
   }
   else if (mon->bench_out->Write(time, state, *mon->fault_op,
                                  mon->seas_op->GetTraction(), V_max))
   {
      mon->bench_out->Flush();
   }

   // ParaView output (adaptive schedule, MPI-collective)
   if (mon->paraview_write_fn)
   {
      mon->paraview_write_fn(static_cast<int>(step), time, V_max);
   }

   // Global output
   if (mpi.IsRoot() && mon->global_out)
   {
      mon->global_out->WriteStep(
         {time, V_max > 0.0 ? std::log10(V_max) : -300.0});
   }

   // Checkpoint
   int istep = static_cast<int>(step);
   if (mon->checkpoint_interval > 0 && istep % mon->checkpoint_interval == 0)
   {
      Vector u_vec;
      u_vec = mon->seas_op->GetDisplacement();
      Vector empty_k0;
      WriteCheckpoint(mon->full_prefix, time, mon->current_dt,
                      istep, mon->num_seismic_events, mon->in_seismic_event,
                      state, u_vec,
                      mon->seas_op->GetTraction(), mon->fault_op->GetSlipRate(),
                      false, empty_k0, mon->mpi);

      // V2 PETSc-TS trailing block (plan §5 monitor site).
      //
      // The monitor is installed ONLY when --petsc-ts is active, so the
      // use_petsc_ts gate is structurally true and elided.  `ts` is the
      // callback's TS parameter — no *petsc_ode indirection (petsc_ode
      // is local to main and OUT OF SCOPE here; R-001).
      PetscReal ts_dt_next_q;
      PetscInt  ts_step_q, ts_rejections_q;
      TSGetTimeStep(ts, &ts_dt_next_q);
      TSGetStepNumber(ts, &ts_step_q);
      TSGetStepRejections(ts, &ts_rejections_q);
      const int    pv_snap          = mon->pv_out
                                      ? mon->pv_out->GetTotalSnapshotsWritten()
                                      : 0;
      const real_t pv_last_write    = mon->pv_out
                                      ? mon->pv_out->GetLastWriteTime()
                                      : -1e30;
      const real_t pv_last_vmax     = mon->pv_out
                                      ? mon->pv_out->GetLastVMax()
                                      :  0.0;
      const int    pv_regime        = mon->pv_out
                                      ? mon->pv_out->GetCurrentRegime()
                                      :  0;
      const int    pv_last_commit   = mon->pv_out
                                      ? mon->pv_out->GetLastCommittedCycle()
                                      : std::numeric_limits<int>::min();   // R-004
      const real_t pv_last_vol_time = mon->pv_out
                                      ? mon->pv_out->GetLastVolumeWriteTime()
                                      : -1e30;                              // R-006
      // R-005: save CUMULATIVE rejection count (carryover + this-run),
      // not this-run alone, so chained restarts preserve prior counts.
      const int    cum_rejects      = mon->restart_rejections_carryover
                                      + static_cast<int>(ts_rejections_q);
      seas::WritePetscTSCheckpoint(mon->full_prefix, time, ts_dt_next_q,
                                   static_cast<int>(ts_step_q),
                                   cum_rejects,
                                   pv_snap, pv_last_write, pv_last_vmax,
                                   pv_regime, pv_last_commit,
                                   pv_last_vol_time,
                                   mon->mpi);
   }

   // Console output
   if (mpi.IsRoot() &&
       (istep % mon->print_step_interval == 0 ||
        V_max > mon->V_threshold_seismic))
   {
      std::cout << std::setw(10) << istep
                << std::setw(16) << std::scientific << std::setprecision(6)
                << time / BP5Params::seconds_per_year
                << std::setw(14) << std::scientific << std::setprecision(3)
                << static_cast<double>(dt)
                << std::setw(16) << std::scientific << std::setprecision(3)
                << V_max
                << std::setw(8) << mon->num_seismic_events
                << "\n";
      std::cout.flush();
   }

   // Debug: stop after first accepted step
   if (mon->debug_first_step_dump)
   {
      if (mpi.IsRoot())
      {
         std::cout << "  [debug] stopping after first accepted step dump\n";
      }
      TSSetConvergedReason(ts, TS_CONVERGED_USER);
   }

   PetscFunctionReturn(PETSC_SUCCESS);
}

#endif // MFEM_USE_PETSC

// ============================================================================
// Main driver
// ============================================================================

int main(int argc, char *argv[])
{
   MPIContext mpi(&argc, &argv);

   // Suppress HDF5 auto-print of internal error stacks.  On the
   // Frontera build, PETSc 3.15 pulls in HDF5 1.10 (libhdf5.so.200)
   // while seas/MFEM uses HDF5 1.14 (libhdf5.so.310); each instance
   // has its own ID table and cross-instance closes produce noisy
   // "can't locate ID (already closed?)" stacks on every ParaView
   // write.  See debug_document/paraview_output_debug_document/
   //   hdf5_diag_noise_2026-05-17.md
   mfem::seas::InstallHdf5ErrorFilter();

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
   bool use_paraview = false;          // Enable ParaView PVD/VTU output
   int  paraview_step_interval = 0;   // 0 = adaptive/time schedule, >0 = every N steps
   real_t paraview_dt = 0.0;          // >0 = fixed time interval (seconds) between writes
   // Adaptive-schedule CLI overrides (negative = inherit schedule default).
   real_t pv_dt_co    = -1.0;         // coseismic interval override (s)
   real_t pv_dt_nu    = -1.0;         // nucleation interval override (s)
   real_t pv_dt_inter = -1.0;         // interseismic interval override (s)
   real_t pv_v_co     = -1.0;         // coseismic V threshold override (m/s)
   real_t pv_v_nu     = -1.0;         // nucleation V threshold override (m/s)
   real_t pv_hyst     = -1.0;         // hysteresis factor override (must be >= 1)
   bool   pv_fault_only = false;      // if true, skip volume PVD Save()
   // Phase 4 (volume PV decouple) — `--paraview-fault-only` is the
   // BP5 driver's pre-existing equivalent of `--no-volume-pv`; the
   // tpv-style names are accepted as aliases for cross-driver parity.
   real_t pv_volume_pv_dt = 0.0;      // --volume-pv-dt X (overrides off)
   // Phase 2b — fault back-end selectors.
   bool   pv_force_vtu          = false;
   bool   pv_force_hdf5         = false;
   bool   pv_legacy_ascii_vtu   = false;
   // Phase 2d.3 — chunk filter selectors (HDF5 only).  -1 / 0 = unset.
   real_t pv_fault_zfp_tol      = 0.0;
   int    pv_fault_deflate_level = -1;
   real_t pv_bulk_zfp_tol       = 0.0;
   int    pv_bulk_deflate_level = -1;
   // Phase 6.3 / 6.3a — volume back-end + primary-collection compression.
   bool   pv_volume_force_vtu   = false;
   bool   pv_volume_force_hdf5  = false;
   real_t pv_volume_zfp_tol     = 0.0;
   int    pv_volume_deflate_level = -1;
   // Phase 3 — snapshot cap.
   int    pv_max_snapshots      = 0;
   // Fault-surface VTU field filter.  Empty ⇒ emit all 12 standard
   // fields + the 5 _k4 fields when stage-4 buffers are present
   // (backward-compatible default).  Non-empty ⇒ only the named
   // CellData arrays are emitted.  Parsed from --paraview-fields.
   std::set<std::string> pv_fault_fields;
   int  max_steps = 10000000;         // Maximum number of time steps
   // v50g: face DOF node type (GaussLobatto has cond(M)=2901 at p=4, ClosedUniform=58)
   int face_basis_type = BasisType::GaussLobatto;
   std::string face_basis_str = "GaussLobatto";
   bool verify_parallel = false;  // Run production-mesh verification diagnostics
   double regression_tolerance = -1.0;  // Negative = no regression check
   std::string ref_prefix;               // Reference file prefix (auto-detect if empty)

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
      if (arg == "--max-steps" && i + 1 < argc)
      {
         max_steps = std::atoi(argv[++i]);
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
      if (arg == "--verify") { verify_parallel = true; }
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
      if (arg == "--paraview") { use_paraview = true; }
      if (arg == "--paraview-every" && i + 1 < argc)
      {
         use_paraview = true;
         paraview_step_interval = std::atoi(argv[++i]);
      }
      if (arg == "--paraview-dt" && i + 1 < argc)
      {
         use_paraview = true;
         paraview_dt = std::atof(argv[++i]);
      }
      if (arg == "--paraview-adaptive") { use_paraview = true; }
      if (arg == "--paraview-dt-co" && i + 1 < argc)
      {
         use_paraview = true;
         pv_dt_co = std::atof(argv[++i]);
      }
      if (arg == "--paraview-dt-nu" && i + 1 < argc)
      {
         use_paraview = true;
         pv_dt_nu = std::atof(argv[++i]);
      }
      if (arg == "--paraview-dt-inter-yr" && i + 1 < argc)
      {
         use_paraview = true;
         pv_dt_inter = std::atof(argv[++i]) * BP5Params::seconds_per_year;
      }
      if (arg == "--paraview-v-co" && i + 1 < argc)
      {
         use_paraview = true;
         pv_v_co = std::atof(argv[++i]);
      }
      if (arg == "--paraview-v-nu" && i + 1 < argc)
      {
         use_paraview = true;
         pv_v_nu = std::atof(argv[++i]);
      }
      if (arg == "--paraview-hyst" && i + 1 < argc)
      {
         use_paraview = true;
         pv_hyst = std::atof(argv[++i]);
      }
      if (arg == "--paraview-fault-only")
      {
         use_paraview = true;
         pv_fault_only = true;
      }
      // Phase 4 — `--no-volume-pv` is the canonical (tpv-driver) name.
      // `--no-domain-pv` is preserved as a deprecated alias.  Both map
      // to the BP5 driver's `pv_fault_only`.
      if (arg == "--no-volume-pv" || arg == "--no-domain-pv")
      {
         use_paraview = true;
         pv_fault_only = true;
      }
      if (arg == "--volume-pv-dt" && i + 1 < argc)
      {
         use_paraview = true;
         pv_volume_pv_dt = std::atof(argv[++i]);
      }
      // Phase 2b — fault back-end selector.
      if (arg == "--paraview-fault-vtu")
      {
         use_paraview = true;
         pv_force_vtu = true;
      }
      if (arg == "--paraview-fault-hdf5")
      {
         use_paraview = true;
         pv_force_hdf5 = true;
      }
      if (arg == "--paraview-fault-legacy-ascii")
      {
         use_paraview = true;
         pv_legacy_ascii_vtu = true;
      }
      // Phase 2d.3 — chunk filter selectors.
      if (arg == "--paraview-fault-zfp-tol" && i + 1 < argc)
      {
         use_paraview = true;
         pv_fault_zfp_tol = std::atof(argv[++i]);
      }
      if (arg == "--paraview-fault-deflate-level" && i + 1 < argc)
      {
         use_paraview = true;
         pv_fault_deflate_level = std::atoi(argv[++i]);
      }
      if (arg == "--paraview-bulk-zfp-tol" && i + 1 < argc)
      {
         use_paraview = true;
         pv_bulk_zfp_tol = std::atof(argv[++i]);
      }
      if (arg == "--paraview-bulk-deflate-level" && i + 1 < argc)
      {
         use_paraview = true;
         pv_bulk_deflate_level = std::atoi(argv[++i]);
      }
      // Phase 6.3 — volume back-end selectors.
      if (arg == "--paraview-volume-vtu")
      {
         use_paraview = true;
         pv_volume_force_vtu = true;
      }
      if (arg == "--paraview-volume-hdf5")
      {
         use_paraview = true;
         pv_volume_force_hdf5 = true;
      }
      // Phase 6.3a — volume PRIMARY-collection compression.
      if (arg == "--paraview-volume-zfp-tol" && i + 1 < argc)
      {
         use_paraview = true;
         pv_volume_zfp_tol = std::atof(argv[++i]);
      }
      if (arg == "--paraview-volume-deflate-level" && i + 1 < argc)
      {
         use_paraview = true;
         pv_volume_deflate_level = std::atoi(argv[++i]);
      }
      // Phase 3 — snapshot cap.
      if (arg == "--paraview-max-snapshots" && i + 1 < argc)
      {
         use_paraview = true;
         pv_max_snapshots = std::atoi(argv[++i]);
      }
      if (arg == "--paraview-fields" && i + 1 < argc)
      {
         // Comma-separated list of CellData field names to emit.
         // Empty list (not providing the flag) ⇒ emit all fields.
         // Example: --paraview-fields slip_rate_strike,normal_stress
         // Silent handling of unknown names; they simply do not
         // appear in the VTU.
         use_paraview = true;
         std::string list = argv[++i];
         std::string token;
         size_t pos = 0;
         while (pos <= list.size()) {
            size_t comma = list.find(',', pos);
            size_t end   = (comma == std::string::npos) ? list.size() : comma;
            token = list.substr(pos, end - pos);
            // strip leading/trailing whitespace
            size_t a = token.find_first_not_of(" \t");
            size_t b = token.find_last_not_of(" \t");
            if (a != std::string::npos && b != std::string::npos)
            {
               pv_fault_fields.insert(token.substr(a, b - a + 1));
            }
            if (comma == std::string::npos) { break; }
            pos = comma + 1;
         }
      }
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
      if (arg == "--regression-tolerance" && i + 1 < argc)
      {
         regression_tolerance = std::atof(argv[++i]);
      }
      if (arg == "--ref-prefix" && i + 1 < argc)
      {
         ref_prefix = argv[++i];
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

   // PETSc TS restart: the V2 trailing-block checkpoint format
   // (`miniapps/seas/io/petsc_ts_checkpoint.hpp`) carries the PETSc TS
   // internal state needed to resume `--restart` + `--petsc-ts`.  The
   // V2 read happens in the restart block below (see "Restart from
   // checkpoint"); the V2 write happens in the monitor + final
   // checkpoint sites.  The old early-abort at this location is
   // intentionally removed.

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

   // -------------------------------------------------------------------
   // RESTART / OUTPUT-DIR collision safety check.
   //
   // When --restart PREFIX is supplied, the driver READS from
   // <dirname(PREFIX)>/...  and WRITES new output (fault.vtkhdf,
   // volume.vtkhdf, *_checkpoint_*.txt, *_fltst_*.txt, *_global.txt,
   // ...) into <output-dir>/...  If those two directories are the
   // same, the new run truncate-overwrites the previous run's
   // outputs — the prior fault.vtkhdf is lost, the prior checkpoint
   // is replaced with the post-restart-final version, etc.  The
   // V2 PETSc-TS restart machinery (plan §"Phase 1") preserves
   // SCHEDULE STATE across the seam, but the on-disk output files
   // do NOT continue (plan §"Out of scope": "The VTKHDF writer
   // currently overwrites").
   //
   // Refuse to start when those paths resolve to the same directory.
   // Compare CANONICAL paths via std::filesystem::weakly_canonical so
   // we catch trailing slashes, "./", relative-vs-absolute, etc.
   // weakly_canonical (vs canonical) handles output_dir not existing
   // yet (we may be about to mkdir it).
   if (!restart_prefix.empty())
   {
      namespace fs = std::filesystem;
      try
      {
         const fs::path restart_path(restart_prefix);
         fs::path restart_dir_path = restart_path.parent_path();
         if (restart_dir_path.empty()) { restart_dir_path = "."; }

         const fs::path restart_canonical =
            fs::weakly_canonical(restart_dir_path);
         const fs::path output_canonical =
            fs::weakly_canonical(fs::path(output_dir));

         if (restart_canonical == output_canonical)
         {
            if (mpi.IsRoot())
            {
               std::cerr
                  << "ERROR: --output-dir (" << output_dir
                  << ") resolves to the SAME directory as the parent "
                  "of --restart (" << restart_dir_path.string()
                  << ").\n"
                  "       Continuing would clobber the previous run's "
                  "outputs (fault.vtkhdf, volume.vtkhdf, "
                  "*_checkpoint_r*.txt, probe CSVs, ...).\n"
                  "       The V2 PETSc-TS restart preserves SCHEDULE "
                  "STATE across the seam, but the on-disk output "
                  "files do NOT continue (the VTKHDF writer "
                  "truncate-overwrites on construction).\n"
                  "       Pick a DIFFERENT --output-dir for the "
                  "restarted run.  Recommended chained-restart "
                  "pattern (used by "
                  "jobs/bp5/bp5_restart_test_v2_dev_2hr.sbatch — "
                  "single base dir + segment_NNN subdirs):\n"
                  "         --output-dir <BASE>/segment_001  (initial run)\n"
                  "         --output-dir <BASE>/segment_002  (restart 1)\n"
                  "         --output-dir <BASE>/segment_003  (restart 2)\n"
                  "         ...                              (increment "
                  "for each link in the chain)\n";
            }
            return 3;
         }

         // R-006 (REVIEW.md round 6): soft warning for parent/child
         // path relationships.  Strict equality (above) catches the
         // most common misuse, but the user can still cause partial
         // clobber by pointing --output-dir at a parent or subdir of
         // the restart's directory.  Print a warning so the operator
         // can decide whether the layout is intentional.
         {
            const std::string r = restart_canonical.string();
            const std::string o = output_canonical.string();
            const bool r_is_parent_of_o =
               (o.size() > r.size())
               && (o.compare(0, r.size(), r) == 0)
               && (o[r.size()] == '/');
            const bool o_is_parent_of_r =
               (r.size() > o.size())
               && (r.compare(0, o.size(), o) == 0)
               && (r[o.size()] == '/');
            if ((r_is_parent_of_o || o_is_parent_of_r)
                && mpi.IsRoot())
            {
               std::cerr
                  << "WARNING: --output-dir and --restart have a "
                  "parent/child directory relationship\n"
                  "         (restart=" << r << ",\n"
                  "          output =" << o << ").\n"
                  "         Phase B's output may partially overlap with "
                  "Phase A's if file names collide.  Consider distinct "
                  "sibling dirs.\n";
            }
         }
      }
      catch (const fs::filesystem_error &e)
      {
         if (mpi.IsRoot())
         {
            std::cerr
               << "ERROR: failed to canonicalise --restart / "
               "--output-dir paths: " << e.what() << "\n"
               "       restart_prefix = " << restart_prefix << "\n"
               "       output_dir     = " << output_dir << "\n";
         }
         return 3;
      }
   }

   std::string full_prefix = output_dir + "/" + output_prefix;

   // =========================================================================
   // Comparison-only mode
   // =========================================================================
   if (comparison_only)
   {
      int rc = 0;
      if (mpi.IsRoot())
      {
         std::string rp = ref_prefix;
         if (rp.empty()) { rp = DetectRefPrefix(ref_dir); }
         if (rp.empty()) { rp = output_prefix; }
         rc = RunComparison(output_dir, output_prefix, ref_dir, rp,
                            stations, t_final, regression_tolerance);
      }
      mpi.Bcast(rc);
#ifdef MFEM_USE_PETSC
      if (petsc_initialized) { MFEMFinalizePetsc(); }
#endif
      return rc;
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
      std::string bc_desc = "FarField";
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
   // Coordinate window selects faces in the region where prior blowups occurred.
   // All ranks participate; max_traced_faces limits output per rank.
   auto env_truthy = [](const char *name) -> bool
   {
      const char *v = std::getenv(name);
      if (!v || !*v) { return false; }
      std::string s(v);
      std::transform(s.begin(), s.end(), s.begin(),
                     [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
      return !(s == "0" || s == "false" || s == "off" || s == "no");
   };
   auto env_int = [](const char *name, int default_value) -> int
   {
      const char *v = std::getenv(name);
      if (!v || !*v) { return default_value; }
      return std::atoi(v);
   };
   const bool debug_first_step_dump = env_truthy("SEAS_DEBUG_FIRST_STEP_DUMP");
   const int debug_first_step_rank =
      env_int("SEAS_DEBUG_FIRST_STEP_TARGET_RANK", 96);
   // Debug gate for per-face trace CSVs (trace_faces / trace_timeseries /
   // trace_events / trace_summary).  These were diagnostic outputs used
   // during v59 fault-tip blowup investigations and are NOT needed for
   // production runs.  Default OFF; opt in with `export
   // SEAS_DEBUG_FACE_TRACE=1` when re-debugging.
   const bool debug_face_trace = env_truthy("SEAS_DEBUG_FACE_TRACE");

   seas::TraceConfig trace_cfg;
   if (debug_face_trace)
   {
      trace_cfg.use_coord_window = true;
      trace_cfg.x2_min = -45e3; trace_cfg.x2_max = -25e3;
      trace_cfg.x3_min = 35e3; trace_cfg.x3_max = 40e3;
      trace_cfg.num_control_faces = 2;
      trace_cfg.max_traced_faces = 50;
   }
   trace_cfg.output_dir = output_dir;
   // SelectFaces with `use_coord_window = false` and `explicit_rank = -1`
   // (defaults) selects no faces, so `active_` stays false and no
   // trace_*.csv files are opened on any rank.
   seas::FaceTraceLogger<ParMesh> face_tracer(trace_cfg, mpi.Rank());
   face_tracer.SelectFaces(domain, fault_geom);
   seas_op.SetFaceTracer(&face_tracer);

   if (debug_first_step_dump)
   {
      ElasticityDomainOperator<ParMesh>::FirstStepDebugConfig dbg_cfg;
      dbg_cfg.enabled = true;
      dbg_cfg.target_rank = debug_first_step_rank;
      dbg_cfg.output_dir = output_dir;
      // Fault-tip faces on rank 96 (1000m mesh, 400 ranks):
      //   fi=4:  near x2=-37 km, depth~39 km (interior, just above tip)
      //   fi=28: x2=-37 km, depth=40 km (tip corner, v59 blowup DOF)
      //   fi=29: x2=-36 km, depth=40 km (adjacent tip)
      //   fi=33: x2=-37 km, depth~39 km (interior control)
      dbg_cfg.target_fault_faces = {4, 28, 29, 33};
      domain.SetFirstStepDebugConfig(dbg_cfg);
      if (mpi.IsRoot())
      {
         std::cout << "  [debug] first-step IP dump: ON"
                   << " (target rank " << debug_first_step_rank
                   << ", fault faces 4/28/29/33)\n";
      }
   }

   // =========================================================================
   // Production-mesh verification: Phase 1 (before K assembly)
   // =========================================================================
   if (verify_parallel)
   {
      if (mpi.IsRoot())
      {
         std::cout << "\n=== Production-Mesh Verification ===\n";
      }
      // Test 1: Ghost DOF communication (no K needed)
      real_t ghost_err = domain.VerifyGhostDOFCommunication();
      (void)ghost_err;
   }

   // R-003 (REVIEW.md 2026-05-16): PetscParVector::PlaceMemory /
   // ResetMemory cannot tolerate a NULL-backed Memory.  On ranks with
   // zero owned fault DOFs (a strike-slip fault embedded in a 3D box
   // leaves many ranks touching no fault face when partitioned over
   // hundreds of MPI tasks), `fault_op.StateSize() == 0`.  The naive
   // `Vector state(0)` does NOT allocate (vector.hpp:574-582: the
   // constructor's `data.New(s)` is gated on `s > 0`), so
   // `state.GetMemory().Empty() == true` (h_ptr == NULL).
   // PlaceMemory aliases the NULL pointer happily, but ResetMemory at
   // the end of TSSolve checks `MFEM_VERIFY(!pdata.Empty(),...)` and
   // aborts (petsc.cpp:899).  This crashed BP5 v62/v63-style runs
   // late in TSSolve on Frontera 8N×400r.  The workaround below pads
   // to at least one element so the underlying Memory is allocated,
   // then shrinks back to the actual size: `SetSize` only
   // re-allocates when `new_size > capacity` (vector.hpp:584-602),
   // so the 1-element allocation survives the shrink.
   //
   // R-105 (REVIEW.md 2026-05-16 round 2) — FRAGILE DEPENDENCY:
   // this fix DEPENDS on MFEM `Vector::SetSize` preserving the
   // underlying Memory allocation on shrink (vector.hpp:584-602:
   // re-allocate only when `new_size > capacity`).  If a future MFEM
   // upgrade changes that contract — e.g., to "always reallocate" or
   // "free memory when shrinking to 0" — the `SetSize(actual)` shrink
   // would free the `SetSize(padded)` allocation and the workaround
   // would silently regress to the original NULL-h_ptr crash.  The
   // companion unit test `seas_test_bp5_petsc_ts_zero_fault_rank`
   // pins this invariant; RUN IT after any MFEM bump (`make
   // test-bp5-petsc-ts-zero-fault-rank` from miniapps/seas).
   //
   // Long-term fix: land the MFEM-side patch documented in the
   // previous REVIEW.md round 1 (R-003 alternative — allow
   // zero-length aliases in `PetscParVector::ResetMemory`), then
   // delete this driver-side workaround.
   Vector state;
   {
      const int actual = fault_op.StateSize();
      const int padded = std::max(actual, 1);
      state.SetSize(padded);
      state.SetSize(actual);  // shrink back; allocation is preserved
   }
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
   // Production-mesh verification: Phase 2 (after K assembly)
   // =========================================================================
   if (verify_parallel)
   {
      // K is now assembled (triggered by SetInitialCondition → Solve)

      // Test 2: Dirichlet skip-set audit (root cause diagnostic)
      domain.VerifyDirichletSkipSets();

      Vector zero_slip(2 * domain.GetNumFaultDOFs());
      zero_slip = 0.0;

      // Test 3: RHS norms at t=0 (slip only, Dirichlet=0)
      domain.VerifyRHSNorms(0.0, zero_slip);

      // Test 4: RHS norms at t=1yr (Dirichlet active)
      domain.VerifyRHSNorms(3.15576e7, zero_slip);

      // Test 5: Per-face shared Dirichlet diagnostic
      domain.VerifySharedDirichletPerFace(3.15576e7);

      if (mpi.IsRoot())
      {
         std::cout << "=== Verification Complete ===\n\n"
                   << "Compare ||b_dir|| and ||b_total|| between 1-rank and\n"
                   << "N-rank runs. Mismatch -> shared-face assembly bug.\n"
                   << "If SIGN_SAME on any face -> dir_sign does not flip.\n\n";
      }
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

   // =========================================================================
   // ParaView output (displacement + fault fields)
   // =========================================================================
   std::unique_ptr<seas::ParaViewOutput<ParMesh>> pv_out;
   // Fault field scratch vectors (local = all faces on this rank)
   Vector pv_local_slip, pv_local_slip_rate, pv_local_traction, pv_local_state;
   Vector pv_local_normal_stress;
   if (use_paraview)
   {
      // Phase 6.3 / 6.3a: validate volume CLI flags + select volume back end.
      if (pv_volume_force_vtu && pv_volume_force_hdf5)
      {
         MFEM_ABORT("--paraview-volume-vtu and --paraview-volume-hdf5 are "
                    "mutually exclusive.");
      }
#ifndef MFEM_USE_HDF5
      if (pv_volume_force_hdf5)
      {
         MFEM_ABORT("--paraview-volume-hdf5 requires the seas-mfem build to "
                    "define MFEM_USE_HDF5=YES; current build has it disabled.");
      }
      if (pv_volume_zfp_tol > 0.0)
      {
         MFEM_ABORT("--paraview-volume-zfp-tol requires MFEM_USE_HDF5=YES; "
                    "current build has it disabled.");
      }
      if (pv_volume_deflate_level >= 0)
      {
         MFEM_ABORT("--paraview-volume-deflate-level requires MFEM_USE_HDF5=YES; "
                    "current build has it disabled.");
      }
#endif
#ifndef MFEM_USE_H5Z_ZFP
      if (pv_volume_zfp_tol > 0.0)
      {
         MFEM_ABORT("--paraview-volume-zfp-tol requires MFEM_USE_H5Z_ZFP=YES; "
                    "current build has it disabled.");
      }
#endif
      if (pv_volume_zfp_tol > 0.0 && pv_volume_deflate_level >= 0)
      {
         MFEM_ABORT("--paraview-volume-zfp-tol and "
                    "--paraview-volume-deflate-level are mutually exclusive — "
                    "choose ZFP-accuracy OR deflate, not both.");
      }
      // Phase 6.4: BP5 has no secondary collection; --paraview-bulk-* is a
      // no-op (R-310 migration note).  Warn so the user knows.
      if (pv_bulk_zfp_tol > 0.0 || pv_bulk_deflate_level >= 0)
      {
         if (mpi.IsRoot())
         {
            mfem::out
               << "warning: --paraview-bulk-zfp-tol / "
                  "--paraview-bulk-deflate-level set on BP5 driver; "
                  "BP5 has no secondary 'bulk' collection — the flags "
                  "have no effect.  To compress the volume PV use "
                  "--paraview-volume-zfp-tol / "
                  "--paraview-volume-deflate-level instead.\n";
         }
      }

      auto volume_mode =
         seas::ParaViewOutput<ParMesh>::DefaultVolumeOutputMode();
      if (pv_volume_force_vtu)
      { volume_mode = seas::ParaViewOutput<ParMesh>::VolumeOutputMode::Vtu; }
      if (pv_volume_force_hdf5)
      { volume_mode = seas::ParaViewOutput<ParMesh>::VolumeOutputMode::Hdf5; }
      // Renamed from "volume" to "kinematics" per
      // PLAN_split_bulk_solutions_2026-05-12.  On-disk filename is
      // <output>/ParaView/kinematics.vtkhdf.  BP5 has no stress
      // collection (quasi-dynamic — no wavefield).
      pv_out = std::make_unique<seas::ParaViewOutput<ParMesh>>(
         output_dir + "/ParaView", pmesh, order,
         /*collection_name=*/"kinematics", volume_mode);
      // Displacement: non-owning pointer to the live ParGridFunction in seas_op
      pv_out->RegisterDomainField("displacement",
         const_cast<ParGridFunction*>(
            &seas_op.GetDisplacement()));

      // MPI rank field: L2 order-0 (constant per element) for partition visualization
      {
         auto *l2_fec = new L2_FECollection(0, 3);
         auto *l2_fes = new ParFiniteElementSpace(&pmesh, l2_fec);
         auto *rank_gf = new ParGridFunction(l2_fes);
         *rank_gf = static_cast<real_t>(mpi.Rank());
         rank_gf->MakeOwner(l2_fec);  // GF owns FEC+FES, freed on destruction
         pv_out->RegisterDomainField("mpi_rank", rank_gf);
      }

      // Fault fields: L2-p0 projection onto volume elements
      pv_out->InitFaultOutputBP5(
         domain.GetFaultInteriorFaces(),
         domain.GetFaultSharedFaces(),
         domain.GetNbfPerFace());

      // Apply optional CellData field filter from --paraview-fields.
      // Empty set ⇒ emit all fields (default).  See ParaViewOutput::
      // SetFaultVTUFields.
      if (!pv_fault_fields.empty())
      {
         pv_out->SetFaultVTUFields(pv_fault_fields);
         if (mpi.IsRoot())
         {
            std::cout << "  Fault-VTU field filter (--paraview-fields):";
            for (const auto &f : pv_fault_fields)
            {
               std::cout << " " << f;
            }
            std::cout << "\n";
         }
      }

      // Pre-allocate local fault vectors
      const int n_local_dofs = domain.GetNumFaultDOFs();
      pv_local_slip.SetSize(2 * n_local_dofs);
      pv_local_slip_rate.SetSize(2 * n_local_dofs);
      pv_local_traction.SetSize(2 * n_local_dofs);
      pv_local_state.SetSize(n_local_dofs);
      pv_local_normal_stress.SetSize(n_local_dofs);

      // Set static friction parameters and fault coordinates for visualization.
      // Expand owned → local so they map to the same faces as the dynamic fields.
      {
         const auto &geom = fault_geom;
         Vector local_a, local_Dc, local_x2, local_x3;
         domain.ExpandOwnedToLocalFault(geom.GetAValues(), local_a, 1);
         domain.ExpandOwnedToLocalFault(geom.GetDcValues(), local_Dc, 1);
         domain.ExpandOwnedToLocalFault(geom.GetCoordsX2(), local_x2, 1);
         domain.ExpandOwnedToLocalFault(geom.GetCoordsX3(), local_x3, 1);
         pv_out->SetFaultParamsBP5(local_a, local_Dc, local_x2, local_x3);
      }

      if (paraview_step_interval > 0)
      {
         pv_out->output_every_n_steps = paraview_step_interval;
      }
      if (paraview_dt > 0.0)
      {
         pv_out->fixed_dt = paraview_dt;
      }

      // Apply CLI overrides to the adaptive schedule (inside the use_paraview
      // guard — pv_out is nullptr otherwise).  Validate() runs once at
      // configuration time; a bad user input aborts BEFORE the first
      // time step instead of inside the hot-path NextRegime.
      {
         auto &sched = pv_out->GetSchedule();
         if (pv_dt_co    > 0) { sched.dt_coseismic      = pv_dt_co; }
         if (pv_dt_nu    > 0) { sched.dt_nucleation     = pv_dt_nu; }
         if (pv_dt_inter > 0) { sched.dt_interseismic   = pv_dt_inter; }
         if (pv_v_co     > 0) { sched.v_coseismic       = pv_v_co; }
         if (pv_v_nu     > 0) { sched.v_nucleation      = pv_v_nu; }
         if (pv_hyst     > 0) { sched.hysteresis_factor = pv_hyst; }
         // Phase 3 — snapshot cap.
         if (pv_max_snapshots > 0)
         { sched.max_total_snapshots = pv_max_snapshots; }
         sched.Validate();
      }

      // Phase 3 — total run time required for cap projection.  Set
      // unconditionally so the cap engages as soon as the user opts in.
      pv_out->SetTotalRunTime(t_final);

      // Phase 4 — volume-PV decouple.  `--volume-pv-dt X` (when > 0)
      // re-enables the volume save at an independent cadence even if
      // `--paraview-fault-only` / `--no-volume-pv` is also set ("the
      // explicit dt wins").  Otherwise `pv_fault_only` translates
      // directly to `SetVolumeSaveEnabled(false)`.
      const bool volume_save_enabled =
         (pv_volume_pv_dt > 0.0) || !pv_fault_only;
      pv_out->SetVolumeSaveEnabled(volume_save_enabled);
      if (pv_volume_pv_dt > 0.0)
      { pv_out->SetVolumePVDt(pv_volume_pv_dt); }

      // Phase 2b — fault back-end selector.
      if (pv_force_vtu && pv_force_hdf5)
      {
         MFEM_ABORT("--paraview-fault-vtu and --paraview-fault-hdf5 are "
                    "mutually exclusive.");
      }
      if (pv_legacy_ascii_vtu && pv_force_hdf5)
      {
         MFEM_ABORT("--paraview-fault-legacy-ascii implies the binary "
                    "VTU back end and is incompatible with "
                    "--paraview-fault-hdf5.");
      }
#ifndef MFEM_USE_HDF5
      if (pv_force_hdf5)
      {
         MFEM_ABORT("--paraview-fault-hdf5 requires the seas-mfem build "
                    "to define MFEM_USE_HDF5=YES; current build has it "
                    "disabled.");
      }
      if (pv_fault_deflate_level >= 0)
      {
         MFEM_ABORT("--paraview-fault-deflate-level requires the seas-mfem "
                    "build to define MFEM_USE_HDF5=YES; current build has "
                    "it disabled.");
      }
      if (pv_bulk_deflate_level >= 0)
      {
         MFEM_ABORT("--paraview-bulk-deflate-level requires the seas-mfem "
                    "build to define MFEM_USE_HDF5=YES; current build has "
                    "it disabled.");
      }
#endif
#ifndef MFEM_USE_H5Z_ZFP
      if (pv_fault_zfp_tol > 0.0)
      {
         MFEM_ABORT("--paraview-fault-zfp-tol requires the seas-mfem "
                    "build to define MFEM_USE_H5Z_ZFP=YES; current build "
                    "has it disabled.");
      }
      if (pv_bulk_zfp_tol > 0.0)
      {
         MFEM_ABORT("--paraview-bulk-zfp-tol requires the seas-mfem "
                    "build to define MFEM_USE_H5Z_ZFP=YES; current build "
                    "has it disabled.");
      }
#endif
      if (pv_fault_zfp_tol > 0.0 && pv_fault_deflate_level >= 0)
      {
         MFEM_ABORT("--paraview-fault-zfp-tol and "
                    "--paraview-fault-deflate-level are mutually "
                    "exclusive — choose ZFP-accuracy OR deflate, not both.");
      }
      if (pv_bulk_zfp_tol > 0.0 && pv_bulk_deflate_level >= 0)
      {
         MFEM_ABORT("--paraview-bulk-zfp-tol and "
                    "--paraview-bulk-deflate-level are mutually "
                    "exclusive — choose ZFP-accuracy OR deflate, not both.");
      }
      if (pv_force_vtu)
      { pv_out->SetFaultOutputMode(seas::ParaViewOutput<ParMesh>::FaultOutputMode::Vtu); }
      if (pv_force_hdf5)
      { pv_out->SetFaultOutputMode(seas::ParaViewOutput<ParMesh>::FaultOutputMode::Hdf5); }
      if (pv_legacy_ascii_vtu)
      {
         pv_out->SetFaultOutputMode(seas::ParaViewOutput<ParMesh>::FaultOutputMode::Vtu);
         pv_out->SetLegacyAsciiVTU(true);
      }
#ifdef MFEM_USE_HDF5
      // Phase 2d.3 — chunk filter (fault).
      if (pv_fault_zfp_tol > 0.0)
      {
         pv_out->SetFaultHDFCompression(
            mfem::ParaViewHDFDataCollection::HDFCompression::ZfpAccuracy,
            pv_fault_zfp_tol);
      }
      else if (pv_fault_deflate_level >= 0)
      {
         pv_out->SetFaultHDFCompression(
            mfem::ParaViewHDFDataCollection::HDFCompression::Deflate,
            static_cast<double>(pv_fault_deflate_level));
      }
      // Phase 6.3a (R-301): --paraview-volume-* controls primary `pv_out`
      // compression.  BP5 has no secondary `pv_bulk_out` so
      // --paraview-bulk-* flags do not route here — see the
      // "no effect" warning above.
      if (pv_volume_zfp_tol > 0.0)
      {
         pv_out->SetVolumeHDFCompression(
            mfem::ParaViewHDFDataCollection::HDFCompression::ZfpAccuracy,
            pv_volume_zfp_tol);
      }
      else if (pv_volume_deflate_level >= 0)
      {
         pv_out->SetVolumeHDFCompression(
            mfem::ParaViewHDFDataCollection::HDFCompression::Deflate,
            static_cast<double>(pv_volume_deflate_level));
      }
#endif

      if (mpi.IsRoot())
      {
         if (paraview_step_interval > 0)
         {
            std::cout << "  ParaView output: ON (every "
                      << paraview_step_interval << " steps)\n";
         }
         else if (paraview_dt > 0.0)
         {
            std::cout << "  ParaView output: ON (every "
                      << paraview_dt / BP5Params::seconds_per_year
                      << " yr)\n";
         }
         else
         {
            const auto &s = pv_out->GetSchedule();
            std::cout << "  ParaView output: ON (adaptive schedule"
                      << (pv_fault_only ? ", fault-only" : "")
                      << ")\n"
                      << "    V_co=" << s.v_coseismic
                      << "  V_nu=" << s.v_nucleation
                      << "  hyst=" << s.hysteresis_factor
                      << "  dt_co=" << s.dt_coseismic << "s"
                      << "  dt_nu=" << s.dt_nucleation << "s"
                      << "  dt_inter="
                      << s.dt_interseismic / BP5Params::seconds_per_year
                      << "yr\n";
         }
      }
   }

   // =========================================================================
   // Fault DOF point cloud (CSV) — diagnostic for coordinate validation
   // =========================================================================
   // Writes one CSV per rank with the EXACT owned-DOF coordinates used by
   // the friction parameter computation. If the points form a clean fault
   // rectangle, the coordinates are correct and any ParaView scatter is a
   // projection artifact. If points are scattered here too, the coordinates
   // are wrong and that's the root cause.
   //
   // Gated behind `SEAS_DEBUG_FAULT_DOF_COORDS`.  Default OFF to keep
   // production output directories clean.  Set the env var to 1 when
   // debugging coordinate / projection issues.
   //
   // NOTE (R-106 REVIEW.md 2026-05-16 round 2): the per-face trace
   // CSVs (trace_faces_*, trace_timeseries_*, trace_events_*,
   // trace_summary_*) are controlled INDEPENDENTLY by
   // `SEAS_DEBUG_FACE_TRACE` (see line ~1493).  Set BOTH env vars to
   // recover the pre-fix "always-on" behaviour for the full debug
   // bundle.  Earlier revisions tied fault_dof_coords to either gate,
   // which was confusing.
   const bool debug_fault_dof_coords =
      env_truthy("SEAS_DEBUG_FAULT_DOF_COORDS");
   if (use_paraview && debug_fault_dof_coords)
   {
      const int n_owned = fault_geom.NumFaultDOFs();  // owned count
      const Vector &x2 = fault_geom.GetCoordsX2();
      const Vector &x3 = fault_geom.GetCoordsX3();
      const Vector &a_vals = fault_geom.GetAValues();
      const Vector &dc_vals = fault_geom.GetDcValues();

      std::string vtp_file = output_dir + "/fault_dof_coords_r"
                           + std::to_string(mpi.Rank()) + ".csv";
      std::ofstream vtp(vtp_file);
      vtp << "owned_dof,x2,x3,param_a,param_Dc,psi_init\n";
      vtp << std::setprecision(10);
      for (int i = 0; i < n_owned; i++)
      {
         real_t psi_i = state(i * 3 + 2);  // BP5: [slip_dip, slip_strike, psi]
         vtp << i << "," << x2(i) << "," << x3(i)
             << "," << a_vals(i) << "," << dc_vals(i)
             << "," << psi_i << "\n";
      }
      vtp.close();
      if (mpi.IsRoot())
      {
         std::cout << "  Fault DOF point cloud: " << vtp_file
                   << " (" << n_owned << " owned DOFs)\n";
      }
   }

   // Helper lambda: update + save ParaView output at a given time step.
   // All ranks must call collectively (ParaViewDataCollection::Save is
   // MPI-collective).  V_max must already be globally reduced.
   auto paraview_write = [&](int step_num, real_t time, real_t V_max)
   {
      if (!pv_out) { return; }

      // -------- Fault-only path (skips volume PVD) --------
      // Gated by PeekShouldWrite (const), populates local vectors, commits
      // the schedule (two-arg: advances last_write_time_ AND current_regime_
      // in lockstep so hysteresis keeps working), then writes the fault
      // surface VTU.  We skip pv_out->UpdateFaultFieldsBP5 (volume-PVD
      // GridFunction update) and pv_out->Save (volume PVD write) — neither
      // is needed when the volume mesh is not being output.
      //
      // Phase 4: query the library's `GetVolumeSaveEnabled()` rather
      // than the driver-local `pv_fault_only`, so `--volume-pv-dt X`
      // (which re-enables the volume save with an independent cadence
      // even when `--no-volume-pv` is also set) takes the volume-PVD
      // path instead of being silently dropped.  Mirrors the
      // tpv102/104/205 driver pattern.
      if (!pv_out->GetVolumeSaveEnabled())
      {
         if (!pv_out->PeekShouldWrite(step_num, time, V_max)) { return; }

         // Populate pv_local_* for WriteFaultSurfaceVTU.  Same five
         // expansions as the volume path below, MINUS UpdateFaultFieldsBP5.
         Vector owned_slip;
         fault_op.GetSlip(state, owned_slip);
         domain.ExpandOwnedToLocalFault(owned_slip, pv_local_slip, 2);
         domain.ExpandOwnedToLocalFault(fault_op.GetSlipRate(),
                                        pv_local_slip_rate, 2);
         domain.ExpandOwnedToLocalFault(seas_op.GetTraction(),
                                        pv_local_traction, 2);
         {
            const int spn = 3;  // BP5: [slip_dip, slip_strike, psi]
            const int n_owned = fault_op.NumNodes();
            Vector owned_psi(n_owned);
            for (int i = 0; i < n_owned; i++)
            {
               owned_psi(i) = state(i * spn + 2);
            }
            domain.ExpandOwnedToLocalFault(owned_psi, pv_local_state, 1);
         }
         if (seas_op.ElasticSigmaNEnabled() &&
             seas_op.GetNormalTraction().Size() > 0)
         {
            domain.ExpandOwnedToLocalFault(seas_op.GetNormalTraction(),
                                           pv_local_normal_stress, 1);
         }
         else
         {
            pv_local_normal_stress = 0.0;
         }

         // Advance BOTH last_write_time_ AND current_regime_ in lockstep.
         // This is the two-arg overload — the one-arg shim would leave
         // current_regime_ pinned at 0, silently defeating Phase 3 hysteresis.
         pv_out->CommitSchedule(time, V_max);

         Vector local_a, local_Dc, local_x2, local_x3;
         domain.ExpandOwnedToLocalFault(fault_geom.GetAValues(),  local_a,  1);
         domain.ExpandOwnedToLocalFault(fault_geom.GetDcValues(), local_Dc, 1);
         domain.ExpandOwnedToLocalFault(fault_geom.GetCoordsX2(), local_x2, 1);
         domain.ExpandOwnedToLocalFault(fault_geom.GetCoordsX3(), local_x3, 1);

         pv_out->WriteFaultSurfaceVTU(
            output_dir, step_num, time, mpi.Rank(), mpi.Size(),
            pv_local_slip, pv_local_slip_rate, pv_local_traction,
            pv_local_state, pv_local_normal_stress,
            local_a, local_Dc, local_x2, local_x3);
         return;
      }

      // -------- Volume-PVD path (default) --------
      // Expand owned fault vectors to local (all faces) for visualization
      Vector owned_slip;
      fault_op.GetSlip(state, owned_slip);
      domain.ExpandOwnedToLocalFault(owned_slip, pv_local_slip, 2);
      domain.ExpandOwnedToLocalFault(fault_op.GetSlipRate(),
                                     pv_local_slip_rate, 2);
      domain.ExpandOwnedToLocalFault(seas_op.GetTraction(),
                                     pv_local_traction, 2);
      // State (psi): extract from ODE state vector, 1 comp per DOF
      {
         const int spn = 3;  // BP5: [slip_dip, slip_strike, psi]
         const int n_owned = fault_op.NumNodes();
         Vector owned_psi(n_owned);
         for (int i = 0; i < n_owned; i++)
         {
            owned_psi(i) = state(i * spn + 2);
         }
         domain.ExpandOwnedToLocalFault(owned_psi, pv_local_state, 1);
      }
      // Normal stress (1 comp per DOF, may be empty if disabled)
      if (seas_op.ElasticSigmaNEnabled() &&
          seas_op.GetNormalTraction().Size() > 0)
      {
         domain.ExpandOwnedToLocalFault(seas_op.GetNormalTraction(),
                                        pv_local_normal_stress, 1);
      }
      else
      {
         pv_local_normal_stress = 0.0;
      }
      pv_out->UpdateFaultFieldsBP5(pv_local_slip, pv_local_slip_rate,
                                   pv_local_traction, pv_local_state,
                                   pv_local_normal_stress);
      bool wrote = pv_out->Save(step_num, time, V_max);

      // Fault surface VTU (proper triangle geometry, no L2-p0 artifacts)
      // Only write when Save() schedule triggers (same gating).
      if (wrote)
      {
         Vector local_a, local_Dc, local_x2, local_x3;
         domain.ExpandOwnedToLocalFault(fault_geom.GetAValues(), local_a, 1);
         domain.ExpandOwnedToLocalFault(fault_geom.GetDcValues(), local_Dc, 1);
         domain.ExpandOwnedToLocalFault(fault_geom.GetCoordsX2(), local_x2, 1);
         domain.ExpandOwnedToLocalFault(fault_geom.GetCoordsX3(), local_x3, 1);
         pv_out->WriteFaultSurfaceVTU(
            output_dir, step_num, time, mpi.Rank(), mpi.Size(),
            pv_local_slip, pv_local_slip_rate, pv_local_traction,
            pv_local_state, pv_local_normal_stress,
            local_a, local_Dc, local_x2, local_x3);
      }
   };

   // Write initial state
   bench_out.ForceWrite(0.0, state, fault_op, seas_op.GetTraction(), V_init);
   bench_out.Flush();
   // t=0 IC write: always triggers because last_write_time_ starts at -1e30.
   paraview_write(0, 0.0, V_init);
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
   // R-005: pre-restart cumulative rejection count, loaded from the V2
   // checkpoint when --restart is supplied with --petsc-ts.  Stays 0 on
   // fresh runs.  Added to TSGetStepRejections at every WritePetscTS
   // call site (so chained restarts don't lose prior runs' counts)
   // and at the post-Run end-of-summary accumulation.
   int restart_rejections_carryover = 0;
   // R-008 (REVIEW.md round 4): captured copy of the V2-authoritative
   // dt set by the V2 restart block.  Asserted equal to `current_dt`
   // at the Run() call site below, so a future CFL clamp / dt_init
   // override inserted between the V2 block and Run() trips a clear
   // failure rather than silently clobbering the restart-state dt.
   // Initialised to -1.0 (sentinel for "V2 block did not run").
   real_t v2_authoritative_dt = -1.0;
   int print_step_interval = 10;
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
   BP5MonitorCtx petsc_mon_ctx{};  // zero-initialized
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

      // Set initial dt once — MFEM's Run() will pass this to TSSolve.
      // Matches Tandem: set dt once, then TSSolve handles everything.
      ierr = TSSetTimeStep(ts, dt_init);
      MFEM_VERIFY(ierr == PETSC_SUCCESS, "TSSetTimeStep(dt_init) failed");

      // Register monitor callback for per-step I/O, earthquake detection,
      // console output, checkpointing — exactly Tandem's TSMonitorSet pattern.
      // MFEM's Run() calls TSSolve; the monitor fires after each accepted step.
      petsc_mon_ctx.mpi = &mpi;
      petsc_mon_ctx.state = &state;
      petsc_mon_ctx.bench_out = &bench_out;
      petsc_mon_ctx.face_tracer = &face_tracer;
      petsc_mon_ctx.fault_op = &fault_op;
      petsc_mon_ctx.seas_op = &seas_op;
      petsc_mon_ctx.global_out = global_out.get();
      petsc_mon_ctx.write_every_step = write_every_step;
      petsc_mon_ctx.print_step_interval = print_step_interval;
      petsc_mon_ctx.checkpoint_interval = checkpoint_interval;
      petsc_mon_ctx.full_prefix = full_prefix;
      petsc_mon_ctx.debug_first_step_dump = debug_first_step_dump;
      petsc_mon_ctx.num_seismic_events = 0;
      petsc_mon_ctx.in_seismic_event = false;
      petsc_mon_ctx.V_threshold_seismic = 1e-3;
      petsc_mon_ctx.V_threshold_interseismic = 1e-6;
      petsc_mon_ctx.paraview_write_fn = paraview_write;
      petsc_mon_ctx.current_dt = dt_init;
      // R-001 / R-005: thread `pv_out` and the rejection carryover
      // through to the monitor callback so its V2 WritePetscTSCheckpoint
      // snippet can reach them (they are local to main and otherwise
      // out of scope inside the static callback).  The carryover is 0
      // on a fresh run; the V2 restart block (further down) will reload
      // it from the checkpoint and re-propagate to petsc_mon_ctx.
      petsc_mon_ctx.pv_out = pv_out.get();
      petsc_mon_ctx.restart_rejections_carryover =
         restart_rejections_carryover;

      ierr = TSMonitorSet(ts, bp5_ts_monitor_callback, &petsc_mon_ctx,
                          nullptr);
      MFEM_VERIFY(ierr == PETSC_SUCCESS, "TSMonitorSet failed");

      current_dt = dt_init;
      if (mpi.IsRoot())
      {
         std::cout << "  PETSc TS initial dt: " << dt_init << " s\n";
         std::cout << "  Time stepping: TSSolve (Tandem-style, no manual TSStep loop)\n";
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
      // R-002 (REVIEW.md round 4): `restart_traction` is intentionally
      // NOT restored — `seas_op` has no `SetTraction` setter, and the
      // traction field is recomputed on the next `ComputeTraction`
      // call from the just-restored displacement + slip.  This means
      // the FIRST post-restart step starts with a freshly-recomputed
      // traction, NOT the byte-identical pre-checkpoint traction.
      // For Phase-1's "tolerance-correct" restart contract this is
      // acceptable: the recomputation differs from the pre-checkpoint
      // traction by at most ~atol (1e-7), well below the trajectory
      // tolerance bound of `atol + rtol*|y|`.  Phase 3 (bit-exact
      // restart, plan §"Phase 3") would require either a SetTraction
      // setter on seas_op + saving traction in V2, OR proving that the
      // recomputation is bit-deterministic from (displacement, slip).
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

#ifdef MFEM_USE_PETSC
   // =========================================================================
   // V2 PETSc-TS restart (plan §"Phase 1" §4).
   //
   // Placed IMMEDIATELY AFTER the V1 restart block (R-302) — the cross-check
   // below depends on `t` having been populated by the V1 `ReadCheckpoint`
   // above.  Do NOT place this inside the PetscTS init block (~line 2269);
   // `t` is still 0 there and the cross-check would always fail.
   // =========================================================================
   if (use_petsc_ts && !restart_prefix.empty())
   {
      real_t ts_t = 0.0, ts_dt_next = 0.0;
      int ts_step = 0, ts_rejections = 0, pv_snapshots = 0;
      real_t ts_last_write_time = -1e30;                          // R-304
      real_t ts_last_v_max       = 0.0;                            // R-304
      int    ts_current_regime   = 0;                              // R-304
      int    ts_last_committed_cycle =
                std::numeric_limits<int>::min();                   // R-004
      real_t ts_last_volume_write_time = -1e30;                    // R-006

      const bool have_ts_state = seas::ReadPetscTSCheckpoint(
         restart_prefix, ts_t, ts_dt_next, ts_step, ts_rejections,
         pv_snapshots, ts_last_write_time, ts_last_v_max,
         ts_current_regime, ts_last_committed_cycle,
         ts_last_volume_write_time, &mpi);
      MFEM_VERIFY(have_ts_state,
                  "--restart with --petsc-ts requires a V2 checkpoint "
                  "(file with a PETSC_TS_V2 trailing block, written by a "
                  "build that includes io/petsc_ts_checkpoint.hpp).  V1 "
                  "checkpoints do not contain PETSc TS state; cannot "
                  "continue.  Re-write with the current build or restart "
                  "on the MFEM time-stepper (--no-petsc-ts).");

      // V1 ReadCheckpoint already set `t = ts_t_v1`.  Cross-check vs. V2.
      //
      // R-003 (REVIEW.md round 4): use max(|t|, 1.0) as the scale so the
      // tolerance does not degenerate to 0 when t == 0 (e.g., a debug
      // checkpoint taken before the first accepted TS step, or a test
      // fixture with a t=0 prefix).  The plan-text format round-trips
      // real_t bit-exactly via 17-digit scientific, so the legitimate
      // diff is 0 in practice; the tolerance only guards against file
      // corruption.
      const real_t cross_check_scale = std::max(std::abs(t), real_t(1.0));
      MFEM_VERIFY(std::abs(t - ts_t) < 1e-12 * cross_check_scale,
                  "Checkpoint inconsistency: V1 time=" << t
                  << " differs from V2 time=" << ts_t);

      petsc::TS ts = *petsc_ode;
      PetscErrorCode ierr;
      // R-002: PetscODESolver::Run() unconditionally calls
      // TSSetTime(ts, t) and TSSetTimeStep(ts, dt) on entry
      // (linalg/petsc.cpp:4362-4363).  An explicit TSSetTime /
      // TSSetTimeStep call HERE would be silently overwritten on the
      // next `petsc_ode->Run(state, t, current_dt, t_final)`.  Instead
      // update the C++ `t` and `current_dt` variables that Run() reads
      // on entry — those are the load-bearing ones.  Only
      // TSSetStepNumber survives Run() (Run never resets the step
      // counter), so it stays.
      ierr = TSSetStepNumber(ts, static_cast<PetscInt>(ts_step));   // R-010
      // The plan suggested PCHKERRQ here, but this driver does not
      // pull in the PETSc private header that defines it.  Match the
      // existing convention at lines 2279/2286 instead.
      MFEM_VERIFY(ierr == PETSC_SUCCESS,
                  "TSSetStepNumber(ts_step=" << ts_step << ") failed");

      // Make V2 authoritative for `t` and `current_dt`.  V1 ReadCheckpoint
      // already set `t = ts_t_v1`; the cross-check above guarantees
      // ts_t == t, so the reassignment is a no-op today.  But
      // `current_dt` was set to V1's restart_dt, which may diverge from
      // ts_dt_next in any future change that adds a CFL clamp or
      // dt_init override.
      t          = ts_t;
      current_dt = ts_dt_next;

      // R-003: TSGetTimeStep can return 0 if the checkpoint was
      // written before the first accepted step or just after a
      // TSSetConvergedReason(TS_DIVERGED_*).  Fall back to dt_init
      // and log on rank 0 so PETSc has a non-zero starting dt.
      if (current_dt <= 0.0)
      {
         if (mpi.IsRoot())
         {
            std::cout << "PETSc TS restart: V2 ts_dt_next was "
                      << current_dt << " <= 0; falling back to "
                      << "dt_init = " << dt_init << " s\n";
         }
         current_dt = dt_init;
      }

      // R-008 (REVIEW.md round 4): snapshot the V2-authoritative dt so
      // the assertion at the Run() call site can catch any subsequent
      // overwrite (CFL clamp, dt_init override).  Includes the R-003
      // fallback so the asserted invariant is "current_dt at Run() ==
      // current_dt at end of V2 block", not "current_dt == ts_dt_next".
      v2_authoritative_dt = current_dt;

      if (pv_out)
      {
         pv_out->SetTotalSnapshotsWritten(pv_snapshots);
         // R-304 + R-006 + R-007: restore schedule state across the
         // seam.  Without this the first ShouldWrite after restart
         // fires unconditionally, the regime state machine resets to
         // interseismic, and (with --volume-pv-dt) the first
         // ForceSaveImpl emits a spurious volume snapshot.
         pv_out->RestoreScheduleState(ts_last_write_time,
                                      ts_last_v_max,
                                      ts_current_regime,
                                      ts_last_volume_write_time);
         // R-004: restore the dedup cycle key so the first
         // post-restart CommitSchedule does not over-bump
         // total_snapshots_written_ by 1.
         pv_out->SetLastCommittedCycle(ts_last_committed_cycle);
      }

      // R-303 / R-005: PRE-restart rejection count.  The post-Run code
      // at the end of TSSolve will OVERWRITE `step_rejections` with
      // the THIS-run count from TSGetStepRejections; we accumulate
      // with the carryover so the summary line reports the sum across
      // the restart seam (and the next checkpoint's V2 block stores
      // the cumulative value, not just this-run's).
      restart_rejections_carryover = ts_rejections;
      // Thread the carryover into the monitor too, so monitor-site
      // V2 writes save the correct cumulative count (R-001 + R-005).
      //
      // R-009 (REVIEW.md round 4): this thread-back is structurally
      // unreachable on the `--no-petsc-ts` path because the V2 restart
      // block entry gate (`use_petsc_ts && !restart_prefix.empty()`)
      // at the top of this block prevents entry.  Therefore
      // `petsc_mon_ctx` is guaranteed to have been initialised by the
      // PetscTS init block earlier, and accessing its members here is
      // safe.  Do NOT remove the gate without re-thinking this.
      petsc_mon_ctx.restart_rejections_carryover =
         restart_rejections_carryover;

      if (mpi.IsRoot())
      {
         std::cout << "PETSc TS restart (V2): t=" << ts_t << " s, "
                   << "dt_next=" << ts_dt_next << " s, step=" << ts_step
                   << ", cumulative_rejections=" << ts_rejections
                   << ", paraview_snapshots=" << pv_snapshots << "\n";
      }
   }
#endif

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

   // =========================================================================
   // Main time-stepping loop
   // =========================================================================
#ifdef MFEM_USE_PETSC
   if (use_petsc_ts)
   {
      // ---- Tandem-style TSSolve: single call, PETSc manages everything ----
      // MFEM's Run() maps state memory into PETSc Vec via PlaceMemory,
      // calls TSSolve, then copies back final time and dt.
      // Per-step monitoring (I/O, earthquake detection, console output)
      // is handled by the bp5_ts_monitor_callback registered above.
      MFEM_VERIFY(petsc_ode, "PETSc TS solver was not initialized");

      if (mpi.IsRoot())
      {
         std::cout << "  Entering TSSolve (t=" << t << " → " << t_final
                   << " s) ...\n";
         std::cout.flush();
      }

      // R-008 (REVIEW.md round 4): if the V2 restart block ran, assert
      // that nothing between the V2 block and here has clobbered the
      // restart-state dt.  `v2_authoritative_dt` is -1.0 on fresh runs
      // (sentinel for "V2 block did not run"), so the check is gated
      // on `v2_authoritative_dt > 0`.  Catches future CFL clamps /
      // dt_init overrides inserted between the V2 block and Run() that
      // would silently re-introduce the R-002 failure mode this whole
      // V2 machinery was built to prevent.
      if (v2_authoritative_dt > 0.0)
      {
         // R-006 (REVIEW.md round 5): bit-exact `==` is intentional.
         // Any modification of current_dt between the V2 block and
         // here — even a value-preserving one like
         // `current_dt = std::min(current_dt, dt_max)` — may change
         // the bit pattern under some compilers and fire this
         // assertion.  That's by design: any insertion here deserves
         // a deliberate re-examination of whether V2 is still the
         // authoritative source of post-restart dt.  If you
         // legitimately need to clamp post-V2 dt, update
         // `v2_authoritative_dt` in the same statement, OR widen
         // this check to a relative tolerance with a documented
         // bound.
         MFEM_VERIFY(current_dt == v2_authoritative_dt,
                     "R-008: current_dt (" << current_dt
                     << ") was modified between the V2 restart block "
                     "and the Run() call (V2 set it to "
                     << v2_authoritative_dt << ").  This breaks the "
                     "R-002 contract that V2 is the source of truth "
                     "for the post-restart dt.  Check for a CFL clamp "
                     "or dt_init override that should be gated on "
                     "restart_prefix.empty().");
      }

      petsc_ode->Run(state, t, current_dt, t_final);

      // Retrieve final step count and rejection count from PETSc
      {
         petsc::TS ts = *petsc_ode;
         PetscInt ts_steps = 0;
         TSGetStepNumber(ts, &ts_steps);
         step = static_cast<int>(ts_steps);
         PetscInt rejects = 0;
         TSGetStepRejections(ts, &rejects);
         // R-303 / R-005: TSGetStepRejections returns THIS-Run's
         // rejections only (it is NOT pre-populated from the V2
         // checkpoint; the checkpointed value lives in
         // restart_rejections_carryover).  Accumulate so the summary
         // line reports the cumulative count across the restart seam.
         step_rejections = restart_rejections_carryover
                         + static_cast<int>(rejects);
      }

      // Copy monitor state back for summary
      num_seismic_events = petsc_mon_ctx.num_seismic_events;
      in_seismic_event = petsc_mon_ctx.in_seismic_event;
   }
   else
#endif
   {
      // ---- Original MFEM Dormand-Prince loop (non-PETSc path) ----
      while (t < t_final && step < max_steps)
      {
         if (t + ode_solver.GetDt() > t_final)
         {
            ode_solver.SetDt(t_final - t);
         }

         real_t dt;
         bool accepted = ode_solver.Step(seas_op, state, t, dt);
         if (!accepted) { continue; }
         current_dt = ode_solver.GetDt();
         step++;

         // Post-step psi clamping is MFEM-specific. Tandem does not do this,
         // so keep it switchable and instrumented while debugging.
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

         // ParaView output (adaptive schedule, MPI-collective)
         paraview_write(step, t, V_max);

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
                      << dt
                      << std::setw(16) << std::scientific << std::setprecision(3)
                      << V_max
                      << std::setw(8) << num_seismic_events
                      << "\n";
            std::cout.flush();
         }

         // Intentional early exit: dump only the first accepted step for
         // cross-verification with Tandem (v59 diagnostic).
         if (debug_first_step_dump)
         {
            if (mpi.IsRoot())
            {
               std::cout << "  [debug] stopping after first accepted step dump\n";
            }
            break;
         }
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

#ifdef MFEM_USE_PETSC
      // V2 PETSc-TS trailing block (plan §5 final site).
      // Final-site is in main, so petsc_ode, use_petsc_ts, pv_out, and
      // restart_rejections_carryover are all in scope as locals.
      if (use_petsc_ts && petsc_ode)
      {
         petsc::TS ts = *petsc_ode;
         PetscReal ts_dt_next_q;
         PetscInt  ts_step_q, ts_rejections_q;
         TSGetTimeStep(ts, &ts_dt_next_q);
         TSGetStepNumber(ts, &ts_step_q);
         TSGetStepRejections(ts, &ts_rejections_q);
         const int    pv_snap          = pv_out
                                         ? pv_out->GetTotalSnapshotsWritten()
                                         : 0;
         const real_t pv_last_write    = pv_out
                                         ? pv_out->GetLastWriteTime()
                                         : -1e30;
         const real_t pv_last_vmax     = pv_out
                                         ? pv_out->GetLastVMax()
                                         :  0.0;
         const int    pv_regime        = pv_out
                                         ? pv_out->GetCurrentRegime()
                                         :  0;
         const int    pv_last_commit   = pv_out
                                         ? pv_out->GetLastCommittedCycle()
                                         : std::numeric_limits<int>::min();  // R-004
         const real_t pv_last_vol_time = pv_out
                                         ? pv_out->GetLastVolumeWriteTime()
                                         : -1e30;                             // R-006
         // R-005: cumulative rejection count, not this-run alone.
         const int    cum_rejects      = restart_rejections_carryover
                                         + static_cast<int>(ts_rejections_q);
         seas::WritePetscTSCheckpoint(full_prefix, t, ts_dt_next_q,
                                      static_cast<int>(ts_step_q),
                                      cum_rejects,
                                      pv_snap, pv_last_write, pv_last_vmax,
                                      pv_regime, pv_last_commit,
                                      pv_last_vol_time, &mpi);
      }
#endif
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
   int regression_rc = 0;
   if (mpi.IsRoot())
   {
      std::string rp = ref_prefix;
      if (rp.empty()) { rp = DetectRefPrefix(ref_dir); }
      if (rp.empty()) { rp = output_prefix; }
      regression_rc = RunComparison(output_dir, output_prefix, ref_dir, rp,
                                    stations, t_final, regression_tolerance);
   }
   mpi.Bcast(regression_rc);

#ifdef MFEM_USE_PETSC
   if (petsc_initialized)
   {
      MFEMFinalizePetsc();
   }
#endif

   return regression_rc;
}
