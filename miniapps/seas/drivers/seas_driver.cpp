// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// SEAS TOML-based driver: reads configuration from a TOML file and runs
// a complete BP5 quasi-dynamic earthquake cycle simulation using all
// new code paths (ConstitutiveModel, BoundaryConfig, DomainConfig).
//
// Usage:
//   mpirun -np N ./seas_driver config.toml [--override key=value ...]
//                              [--max-steps N] [--write-every-step]
//
// For BP5: use config/bp5_example.toml as a starting point.

#include "mfem.hpp"
#ifdef MFEM_USE_PETSC
#include "petsc.h"
#if PETSC_VERSION_LT(3,19,0)
#define PETSC_SUCCESS 0
#endif
#endif
#include "../solver/seas_operator.hpp"
#include "../solver/time_stepper.hpp"
#include "../domain/elasticity_operator.hpp"
#include "../domain/boundary_config.hpp"
#include "../fault/fault_geometry.hpp"
#include "../fault/rate_state_fault.hpp"
#include "../friction/dieterich_ruina.hpp"
#include "../friction/state_evolution.hpp"
#include "../config/bp5_params.hpp"
#include "../config/seas_config.hpp"
#include "../config/seas_config_parser.hpp"
#include "../config/seas_config_bridge.hpp"
#include "../config/bp5_mesh_utils.hpp"
#include "../io/bp5_parallel_output.hpp"
#include "../io/probe_output.hpp"
// #include "../io/checkpoint.hpp"  // Full checkpoint needs displacement/traction vectors
#include "../common/mpi_context.hpp"
#include "../constitutive/linear_elastic.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <filesystem>

using namespace mfem;
using namespace mfem::seas;

// ============================================================================
// PETSc TS monitor context and callback
// ============================================================================
#ifdef MFEM_USE_PETSC

/// Context passed to PETSc TSMonitor callback.
struct DriverMonitorCtx
{
   MPIContext *mpi;
   Vector *state;
   ParallelBP5BenchmarkOutput *bench_out;
   RateStateFaultOperator<ParMesh, 2> *fault_op;
   PBP5SEASOp *seas_op;
   ProbeOutput *global_out;  // may be nullptr (non-root)

   bool write_every_step;
   int  print_step_interval;
   int  checkpoint_interval;
   std::string full_prefix;

   int  num_seismic_events;
   bool in_seismic_event;
   real_t V_threshold_seismic;
   real_t V_threshold_interseismic;
   real_t current_dt;
};

/// PETSc TSMonitor callback — called after every accepted step inside TSSolve.
static PetscErrorCode driver_ts_monitor_callback(
   TS ts, PetscInt step, PetscReal time, Vec u, void *ctx)
{
   PetscFunctionBeginUser;
   auto *mon = static_cast<DriverMonitorCtx*>(ctx);

   Vector &state  = *mon->state;
   auto   &mpi    = *mon->mpi;

   PetscReal dt;
   TSGetTimeStep(ts, &dt);
   mon->current_dt = dt;

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

   // Global output
   if (mpi.IsRoot() && mon->global_out)
   {
      mon->global_out->WriteStep(
         {time, V_max > 0.0 ? std::log10(V_max) : -300.0});
   }

   // Console output
   int istep = static_cast<int>(step);
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

   PetscFunctionReturn(PETSC_SUCCESS);
}

#endif // MFEM_USE_PETSC

int main(int argc, char *argv[])
{
   MPIContext mpi(&argc, &argv);

   // =========================================================================
   // Stage 1: Parse TOML config + CLI overrides
   // =========================================================================
   if (argc < 2)
   {
      if (mpi.IsRoot())
      {
         std::cerr << "Usage: " << argv[0]
                   << " config.toml [--override key=value ...] "
                      "[--max-steps N] [--write-every-step]\n";
      }
      return 1;
   }

   std::string config_file = argv[1];
   std::vector<std::string> overrides;
   bool write_every_step = false;

   for (int i = 2; i < argc; i++)
   {
      std::string arg(argv[i]);
      if (arg == "--override" && i + 1 < argc)
      {
         overrides.push_back(argv[++i]);
      }
      else if (arg == "--max-steps" && i + 1 < argc)
      {
         overrides.push_back("time.max_steps=" + std::string(argv[++i]));
      }
      else if (arg == "--write-every-step")
      {
         write_every_step = true;
      }
   }

   SEASConfig config = SEASConfigParser::ParseFile(config_file);
   SEASConfigParser::ApplyCLIOverrides(config, overrides);
   SEASConfigParser::Validate(config);

   // PETSc initialization (must happen before any PETSc calls)
   bool petsc_initialized = false;
   bool use_petsc_ts = config.time.use_petsc_ts;
#ifndef MFEM_USE_PETSC
   if (use_petsc_ts)
   {
      if (mpi.IsRoot())
      {
         std::cerr << "WARNING: use_petsc_ts=true but MFEM built without PETSc. "
                   << "Falling back to DormandPrince RK45.\n";
      }
      use_petsc_ts = false;
   }
#else
   if (use_petsc_ts)
   {
      MFEMInitializePetsc(&argc, &argv,
                          config.time.petsc_ts_options.c_str(), NULL);
      petsc_initialized = true;
   }
#endif

   // Build bridge objects
   BP5Params params = BuildBP5Params(config);
   params.Validate();
   DGMethod dg_method = ParseDGMethod(config.solver.dg_method);
   SolverType solver_type = ParseSolverType(config.solver.solver_type);
   BoundaryConfig bdr_config = BuildBoundaryConfig(config.boundary, config.loading.Vp);
   DomainConfig domain_config = BuildDomainConfig(config.solver);

   if (mpi.IsRoot())
   {
      std::cout << "SEAS Driver\n"
                << "================================\n"
                << "  Config:    " << config_file << "\n"
                << "  Benchmark: "
                << (config.benchmark.empty() ? "(none)" : config.benchmark) << "\n"
                << "  DG method: " << config.solver.dg_method << "\n"
                << "  Solver:    " << config.solver.solver_type << "\n"
                << "  Stepper:   " << (use_petsc_ts ? "PETSc TS RK45" : "DormandPrince RK45") << "\n"
                << "  Ranks:     " << mpi.Size() << "\n\n";
   }

   // =========================================================================
   // Stage 2: Load mesh
   // =========================================================================
   std::unique_ptr<Mesh> serial_mesh;
   bool is_inline = (config.mesh.file == "__inline__");

   if (is_inline)
   {
      serial_mesh = CreateBP5InlineMesh(2, 2, 1, 200e3, 100e3, 100e3);
      if (mpi.IsRoot())
      {
         std::cout << "  Mesh: inline (" << serial_mesh->GetNE()
                   << " elements)\n";
      }
   }
   else
   {
      serial_mesh = std::make_unique<Mesh>(config.mesh.file.c_str(), 1, 1);
      if (config.mesh.scale != 1.0)
      {
         for (int i = 0; i < serial_mesh->GetNV(); i++)
         {
            real_t *v = serial_mesh->GetVertex(i);
            v[0] *= config.mesh.scale;
            v[1] *= config.mesh.scale;
            v[2] *= config.mesh.scale;
         }
      }
      if (mpi.IsRoot())
      {
         std::cout << "  Mesh: " << config.mesh.file << " ("
                   << serial_mesh->GetNE() << " elements)\n";
      }
   }

   ParMesh pmesh(mpi.GetComm(), *serial_mesh);
   serial_mesh.reset();

   // GetGlobalNE() does MPI_Allreduce — must be called by ALL ranks.
   long long global_ne = pmesh.GetGlobalNE();

   real_t h_min, h_max, kappa_min, kappa_max;
   pmesh.GetCharacteristics(h_min, h_max, kappa_min, kappa_max);

   if (mpi.IsRoot())
   {
      std::cout << "  Global elements: " << global_ne << "\n"
                << "  h_min = " << h_min << " m, h_max = " << h_max << " m\n\n";
   }

   // =========================================================================
   // Stage 3: Domain operator (NEW constructor path)
   // =========================================================================
   LinearElastic material(params.lambda(), params.mu());

   ElasticityDomainOperator<ParMesh> domain(
      pmesh, config.mesh.order,
      material,
      params.Vp, params.Wf, params.lf,
      bdr_config,
      dg_method, solver_type, domain_config);

   if (mpi.IsRoot())
   {
      std::cout << "  Fault DOFs (global): "
                << domain.GetNumFaultDOFs() << "\n";
   }

   // =========================================================================
   // Stage 4: Fault components
   // =========================================================================
   FaultGeometry<ParMesh> fault_geom(domain, params, &mpi);

   DieterichRuinaFriction::Constants fc;
   fc.V0 = params.V0;
   fc.f0 = params.f0;
   fc.b = params.b;
   fc.Dc = params.L0;
   DieterichRuinaFriction friction(fc);
   AgingLawPsi aging(params.b, params.V0, params.f0);

   RateStateFaultOperator<ParMesh, 2> fault_op(
      &fault_geom, &friction, &aging, params, &mpi);

   // =========================================================================
   // Stage 5: SEAS operator + initial condition
   // =========================================================================
   PBP5SEASOp seas_op(&domain, &fault_op, &mpi);
   seas_op.SetElasticSigmaN(true);

   Vector state(fault_op.StateSize());
   seas_op.SetInitialCondition(state);

   real_t V_init = seas_op.GetMaxSlipRate();
   if (mpi.IsRoot())
   {
      std::cout << "\n  Initial V_max = " << V_init << " m/s\n"
                << "  State size:    " << fault_op.StateSize() << "\n\n";
   }

   // =========================================================================
   // Stage 6: I/O setup
   // =========================================================================
   std::string output_dir = config.output.output_dir;
   std::string output_prefix = config.output.output_prefix;
   std::string full_prefix = output_dir + "/" + output_prefix;

   auto stations = BP5BenchmarkOutput<ParMesh>::DefaultStations();

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

   // Ensure the output directory exists before any writer opens a file there.
   // ProbeOutput / ParallelBP5BenchmarkOutput open <output_dir>/<prefix>_*.txt
   // directly; without this the run aborts ("Cannot open probe output file")
   // whenever output_dir does not already exist.  Matches spatial_dyn_driver.
   if (mpi.IsRoot() && !output_dir.empty() && output_dir != ".")
   {
      std::error_code ec;
      std::filesystem::create_directories(output_dir, ec);
      if (ec)
      {
         std::cerr << "ERROR: could not create output_dir '" << output_dir
                   << "': " << ec.message() << "\n";
      }
   }
   mpi.Barrier();

   // Clean up old output files
   if (mpi.IsRoot())
   {
      for (const auto &st : stations)
      {
         std::remove((full_prefix + "_" + st.name + ".txt").c_str());
      }
      std::remove((full_prefix + "_global.txt").c_str());
   }
   mpi.Barrier();

   ParallelBP5BenchmarkOutput bench_out(
      full_prefix, params, stations, fault_geom, mpi,
      local_x2, local_x3, local_tp_dip, local_tp_strike,
      domain.GetNbfPerFace(), config.solver.face_basis_type);

   std::unique_ptr<ProbeOutput> global_out;
   if (mpi.IsRoot())
   {
      global_out = std::make_unique<ProbeOutput>(
         full_prefix + "_global.txt",
         std::vector<std::string>{"time(s)", "log10(Vmax)(m/s)"},
         "BP5-QD global output");
   }

   // =========================================================================
   // Stage 7: Time stepper setup
   // =========================================================================

   // CFL-aware initial dt (used by both PETSc and native paths)
   int dim = 3;
   real_t c_N_1 = config.mesh.order * (config.mesh.order + dim - 1.0) / dim;
   real_t c_N_1_ref = 2.0 * (2.0 + dim - 1.0) / dim;
   real_t beta = 4.0 * c_N_1 / c_N_1_ref;
   real_t V_max_init = std::max(V_init, params.V_nuc);
   real_t dt_V = std::min(1e3, 0.01 * params.L_nuc / std::max(V_max_init, 1e-20));
   real_t dt_CFL = (h_min > 0) ? 2.0 * params.eta() * h_min / (beta * params.mu()) : dt_V;
   real_t dt_init = std::min(dt_V, dt_CFL);

   if (config.time.tandem_time_stepping)
   {
      dt_init = 0.01;
   }

   // Native DormandPrince solver (used when PETSc-TS is not active)
   DormandPrinceRK45 ode_solver;
   if (!use_petsc_ts)
   {
      ode_solver.SetAbsTol(config.time.atol);
      ode_solver.SetRelTol(config.time.rtol);
      ode_solver.SetDtMin(1e-6);
      ode_solver.SetDtMax(0.1 * BP5Params::seconds_per_year);
      ode_solver.SetMPIContext(&mpi);

      if (!config.time.tandem_time_stepping)
      {
         ode_solver.SetVGuard(100.0);
      }

      ode_solver.SetDt(dt_init);
      ode_solver.Init(seas_op);
   }

#ifdef MFEM_USE_PETSC
   std::unique_ptr<PetscODESolver> petsc_ode;
   DriverMonitorCtx petsc_mon_ctx{};
   if (use_petsc_ts)
   {
      petsc_ode = std::make_unique<PetscODESolver>(mpi.GetComm(), "");
      petsc_ode->SetAbsTol(config.time.atol);
      petsc_ode->SetRelTol(config.time.rtol);
      petsc_ode->SetMaxIter(config.time.max_steps);
      petsc_ode->Init(seas_op, PetscODESolver::ODE_SOLVER_GENERAL);
      petsc::TS ts = *petsc_ode;

      // Re-enable adaptive stepping (MFEM disables it by default)
      {
         TSAdapt tsad;
         PetscErrorCode ierr2 = TSGetAdapt(ts, &tsad);
         MFEM_VERIFY(ierr2 == PETSC_SUCCESS, "TSGetAdapt failed");
         ierr2 = TSAdaptSetType(tsad, TSADAPTBASIC);
         MFEM_VERIFY(ierr2 == PETSC_SUCCESS, "TSAdaptSetType(BASIC) failed");
      }

      PetscErrorCode ierr = TSSetFromOptions(ts);
      MFEM_VERIFY(ierr == PETSC_SUCCESS, "TSSetFromOptions failed");

      ierr = TSSetTimeStep(ts, dt_init);
      MFEM_VERIFY(ierr == PETSC_SUCCESS, "TSSetTimeStep(dt_init) failed");

      // Register monitor callback
      petsc_mon_ctx.mpi = &mpi;
      petsc_mon_ctx.state = &state;
      petsc_mon_ctx.bench_out = &bench_out;
      petsc_mon_ctx.fault_op = &fault_op;
      petsc_mon_ctx.seas_op = &seas_op;
      petsc_mon_ctx.global_out = global_out.get();
      petsc_mon_ctx.write_every_step = write_every_step;
      petsc_mon_ctx.print_step_interval = 10;
      petsc_mon_ctx.checkpoint_interval = config.time.checkpoint_interval;
      petsc_mon_ctx.full_prefix = full_prefix;
      petsc_mon_ctx.num_seismic_events = 0;
      petsc_mon_ctx.in_seismic_event = false;
      petsc_mon_ctx.V_threshold_seismic = 1e-3;
      petsc_mon_ctx.V_threshold_interseismic = 1e-6;
      petsc_mon_ctx.current_dt = dt_init;

      ierr = TSMonitorSet(ts, driver_ts_monitor_callback, &petsc_mon_ctx,
                          nullptr);
      MFEM_VERIFY(ierr == PETSC_SUCCESS, "TSMonitorSet failed");
   }
#endif

   if (mpi.IsRoot())
   {
      std::cout << "  dt_init = " << dt_init << " s"
                << (config.time.tandem_time_stepping ? " (Tandem-style)" : " (CFL)")
                << "\n\n";
   }

   // Write initial state (t=0)
   bench_out.ForceWrite(0.0, state, fault_op, seas_op.GetTraction(), V_init);
   bench_out.Flush();
   if (global_out)
   {
      real_t log10_V = (V_init > 0) ? std::log10(V_init) : -300.0;
      global_out->WriteStep({0.0, log10_V});
   }

   // =========================================================================
   // Stage 8: Time loop
   // =========================================================================
   real_t t = 0.0;
   real_t t_final = config.time.t_final;
   int max_steps = config.time.max_steps;
   int step = 0;
   int eq_count = 0;
   bool in_event = false;
   real_t V_max = V_init;
   int step_rejections = 0;
   real_t current_dt = dt_init;

   if (mpi.IsRoot())
   {
      std::cout << std::setw(10) << "Step"
                << std::setw(16) << "Time [yr]"
                << std::setw(14) << "dt [s]"
                << std::setw(16) << "V_max [m/s]"
                << std::setw(8) << "EQs" << "\n"
                << std::string(64, '-') << "\n";
   }

#ifdef MFEM_USE_PETSC
   if (use_petsc_ts)
   {
      // ---- PETSc TSSolve: single call, PETSc manages everything ----
      MFEM_VERIFY(petsc_ode, "PETSc TS solver was not initialized");

      if (mpi.IsRoot())
      {
         std::cout << "  Entering TSSolve (t=" << t << " → " << t_final
                   << " s) ...\n";
         std::cout.flush();
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
         step_rejections = static_cast<int>(rejects);
      }

      eq_count = petsc_mon_ctx.num_seismic_events;
      in_event = petsc_mon_ctx.in_seismic_event;
      V_max = seas_op.GetMaxSlipRate();
   }
   else
#endif
   {
      // ---- Native DormandPrince RK45 loop ----
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

         V_max = seas_op.GetMaxSlipRate();

         // NaN/Inf check
         if (std::isnan(V_max) || std::isinf(V_max))
         {
            if (mpi.IsRoot())
            {
               std::cerr << "\n[FATAL] NaN/Inf detected at step " << step
                         << ", t = " << t << " s. Aborting.\n";
            }
            return 2;
         }

         // Earthquake detection
         if (!in_event && V_max > 1e-3)
         {
            eq_count++;
            in_event = true;
         }
         else if (in_event && V_max < 1e-6)
         {
            in_event = false;
         }

         // I/O: adaptive schedule or every step
         if (write_every_step)
         {
            bench_out.ForceWrite(t, state, fault_op, seas_op.GetTraction(), V_max);
            bench_out.Flush();
         }
         else if (bench_out.Write(t, state, fault_op, seas_op.GetTraction(), V_max))
         {
            bench_out.Flush();
         }

         // Global output (every accepted step)
         if (global_out)
         {
            real_t log10_V = (V_max > 0) ? std::log10(V_max) : -300.0;
            global_out->WriteStep({t, log10_V});
         }

         // Console output
         if (step % 10 == 0 && mpi.IsRoot())
         {
            real_t t_yr = t / BP5Params::seconds_per_year;
            std::cout << std::setw(10) << step
                      << std::scientific << std::setprecision(6)
                      << std::setw(16) << t_yr
                      << std::setw(14) << dt
                      << std::setw(16) << V_max
                      << std::setw(8) << eq_count
                      << "\n" << std::flush;
         }

         // Checkpoint
         if (config.time.checkpoint_interval > 0 &&
             step % config.time.checkpoint_interval == 0 && mpi.IsRoot())
         {
            std::cout << "Checkpoint at step " << step << ", t=" << t << " s\n";
         }
      }
      step_rejections = ode_solver.GetTotalRejections();
   }

   // =========================================================================
   // Stage 9: Final output
   // =========================================================================
   bench_out.ForceWrite(t, state, fault_op, seas_op.GetTraction(), V_max);
   bench_out.Close();

   if (mpi.IsRoot())
   {
      std::cout << "\n=== SEAS Simulation Complete ===\n"
                << "  Final time: " << t / BP5Params::seconds_per_year << " years\n"
                << "  Steps:      " << step << "\n"
                << "  Rejections: " << step_rejections << "\n"
                << "  Events:     " << eq_count << "\n"
                << "  Stepper:    " << (use_petsc_ts ? "PETSc TS" : "DormandPrince RK45") << "\n"
                << "================================\n";
   }

#ifdef MFEM_USE_PETSC
   if (petsc_initialized) { MFEMFinalizePetsc(); }
#endif

   return 0;
}
