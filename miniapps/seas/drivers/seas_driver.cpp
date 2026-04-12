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

using namespace mfem;
using namespace mfem::seas;

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
   DormandPrinceRK45 ode_solver;
   ode_solver.SetAbsTol(config.time.atol);
   ode_solver.SetRelTol(config.time.rtol);
   ode_solver.SetDtMin(1e-6);
   ode_solver.SetDtMax(0.1 * BP5Params::seconds_per_year);
   ode_solver.SetMPIContext(&mpi);

   // CFL-aware initial dt
   int dim = 3;
   real_t c_N_1 = config.mesh.order * (config.mesh.order + dim - 1.0) / dim;
   real_t c_N_1_ref = 2.0 * (2.0 + dim - 1.0) / dim;
   real_t beta = 4.0 * c_N_1 / c_N_1_ref;
   real_t V_max_init = std::max(V_init, params.V_nuc);
   real_t dt_V = std::min(1e3, 0.01 * params.L_nuc / std::max(V_max_init, 1e-20));
   real_t dt_CFL = (h_min > 0) ? 2.0 * params.eta() * h_min / (beta * params.mu()) : dt_V;
   real_t dt_init = std::min(dt_V, dt_CFL);

   // Tandem-style time stepping overrides
   if (config.time.tandem_time_stepping)
   {
      dt_init = 0.01;
      // Tandem-style: no V-guard rejection (matches old driver line 882)
   }
   else
   {
      // Default: V-guard ON (factor=100) — reject RK stages where V > 100*V_stage0
      // Matches old driver (bp5_verification_full.cpp:1712-1715)
      ode_solver.SetVGuard(100.0);
   }

   ode_solver.SetDt(dt_init);
   ode_solver.Init(seas_op);

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

   const real_t V_seismic = 1e-3;
   const real_t V_interseismic = 1e-6;

   if (mpi.IsRoot())
   {
      std::cout << std::setw(10) << "Step"
                << std::setw(16) << "Time [yr]"
                << std::setw(14) << "dt [s]"
                << std::setw(16) << "V_max [m/s]"
                << std::setw(8) << "EQs" << "\n"
                << std::string(64, '-') << "\n";
   }

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
      if (!in_event && V_max > V_seismic)
      {
         eq_count++;
         in_event = true;
      }
      else if (in_event && V_max < V_interseismic)
      {
         in_event = false;
      }

      // I/O: adaptive schedule (SaveScheduler) or every step
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

      // Checkpoint (periodic log only — full checkpoint requires displacement/traction)
      if (config.time.checkpoint_interval > 0 &&
          step % config.time.checkpoint_interval == 0 && mpi.IsRoot())
      {
         std::cout << "Checkpoint at step " << step << ", t=" << t << " s\n";
      }
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
                << "  Rejections: " << ode_solver.GetTotalRejections() << "\n"
                << "  Events:     " << eq_count << "\n"
                << "================================\n";
   }

   return 0;
}
