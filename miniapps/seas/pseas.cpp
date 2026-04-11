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

// Parallel SEAS BP2 driver
//
// Runs the SEAS quasi-dynamic simulation in parallel using MPI.
// Usage: mpirun -np N ./pseas [--quick] [--nz NZ] [--tfinal T] [--prefix P]

#include "mfem.hpp"
#include "solver/seas_operator.hpp"
#include "solver/time_stepper.hpp"
#include "domain/antiplane_operator.hpp"
#include "fault/fault_geometry.hpp"
#include "fault/rate_state_fault.hpp"
#include "friction/dieterich_ruina.hpp"
#include "friction/state_evolution.hpp"
#include "config/bp2_params.hpp"
#include "domain/bp2_mesh.hpp"
#include "io/parallel_benchmark_output.hpp"
#include "common/mpi_context.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <memory>
#include <vector>
#include <string>
#include <cstdlib>

using namespace mfem;
using namespace mfem::seas;

int main(int argc, char *argv[])
{
   // Initialize MPI
   MPIContext mpi(&argc, &argv);

   // Parse command-line arguments
   bool quick_mode = false;
   int nz_override = 0;
   double tfinal_override = 0.0;
   std::string prefix_override;
   int order = 1;
   for (int i = 1; i < argc; i++)
   {
      std::string arg(argv[i]);
      if (arg == "--quick" || arg == "-q") { quick_mode = true; }
      if (arg == "--nz" && i + 1 < argc) { nz_override = std::atoi(argv[++i]); }
      if (arg == "--prefix" && i + 1 < argc) { prefix_override = argv[++i]; }
      if (arg == "--tfinal" && i + 1 < argc) { tfinal_override = std::atof(argv[++i]); }
      if (arg == "--order" && i + 1 < argc) { order = std::atoi(argv[++i]); }
   }

   // Simulation parameters
   BP2Params params;

   int mesh_nx = 5;
   int mesh_nz = 125;
   real_t grading_x = 3.0;
   real_t grading_z = 1.0;
   real_t Lx = 50.0e3;
   real_t Lz = 100.0e3;
   real_t t_final;

   if (quick_mode)
   {
      t_final = 1e10;  // ~317 years
   }
   else
   {
      t_final = 1.2e10;  // ~380 years
   }

   if (nz_override > 0) { mesh_nz = nz_override; }
   if (tfinal_override > 0.0) { t_final = tfinal_override; }
   std::string output_prefix = prefix_override.empty() ?
      "par_bp2qd" : prefix_override;

   params.t_final = t_final;

   // =========================================================================
   // Mesh: create serial mesh on all ranks, then distribute
   // =========================================================================
   if (mpi.IsRoot())
   {
      std::cout << "SEAS Parallel BP2 Driver (np=" << mpi.Size() << ")\n";
      std::cout << "  Mesh: " << 2*mesh_nx << "x" << mesh_nz
                << " graded, Lx=" << Lx/1e3 << "km, Lz=" << Lz/1e3 << "km\n";
   }

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
   serial_mesh.reset();  // Free serial mesh

   // GetGlobalNE() is collective (MPI_Allreduce) — call on all ranks
   long long global_ne = pmesh.GetGlobalNE();
   if (mpi.IsRoot())
   {
      std::cout << "  ParMesh: " << global_ne << " global elements\n";
   }

   // =========================================================================
   // Domain operator: DG
   // =========================================================================
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

   // Gather and deduplicate fault depths (handles DG partition-boundary dupes)
   Vector dedup_fault_depths;
   {
      // Use a dummy field (depths themselves) to get dedup_depths
      Vector dedup_data;
      fault_geom.GatherToRootDedup(local_fault_depths, dedup_data,
                                    dedup_fault_depths);
   }

   if (mpi.IsRoot())
   {
      std::cout << "  Deduplicated global fault DOFs: "
                << dedup_fault_depths.Size() << "\n";
   }

   std::vector<real_t> probe_depths = {
      0.0, -2400.0, -4800.0, -7200.0, -9600.0,
      -12000.0, -14400.0, -16800.0, -19200.0,
      -24000.0, -28800.0, -36000.0
   };

   // Clean up old output files
   if (mpi.IsRoot())
   {
      std::string rm_cmd = "rm -f " + output_prefix + "_z*.txt";
      std::system(rm_cmd.c_str());
   }
   mpi.Barrier();

   ParallelBenchmarkOutput bench_out(
      output_prefix, params, probe_depths, fault_geom, mpi,
      dedup_fault_depths);

   // Write initial state
   bench_out.ForceWrite(0.0, state, fault_op, seas_op.GetTraction(),
                        V_init);

   // =========================================================================
   // Time integration (Dormand-Prince RK45)
   // =========================================================================
   DormandPrinceRK45 ode_solver;
   ode_solver.SetMPIContext(&mpi);
   ode_solver.SetAbsTol(1e-7);
   ode_solver.SetRelTol(1e-50);  // Match Tandem/PETSc: pure absolute tolerance
   ode_solver.SetDtMin(1e-6);
   ode_solver.SetDtMax(0.5 * BP2Params::seconds_per_year);
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

      real_t V_max = seas_op.GetMaxSlipRate();  // Already global

      // Check for NaN/Inf
      if (std::isnan(V_max) || std::isinf(V_max))
      {
         if (mpi.IsRoot())
         {
            std::cerr << "NaN/Inf detected at step " << step
                      << ", t = " << t / BP2Params::seconds_per_year << " yr\n";
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
                      << ", V_max = " << std::scientific << std::setprecision(2)
                      << V_max << " m/s ***\n";
         }
      }
      else if (in_seismic_event && V_max < V_threshold_interseismic)
      {
         in_seismic_event = false;
      }

      // I/O
      if (bench_out.Write(t, state, fault_op, seas_op.GetTraction(), V_max))
      {
         bench_out.Flush();
      }

      // Periodic console output (root only)
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

   // Final I/O
   bench_out.ForceWrite(t, state, fault_op, seas_op.GetTraction(),
                        seas_op.GetMaxSlipRate());
   bench_out.Close();

   if (mpi.IsRoot())
   {
      std::cout << "\n=== Parallel Simulation Summary ===\n";
      std::cout << "  Ranks: " << mpi.Size() << "\n";
      std::cout << "  Final time: " << t / BP2Params::seconds_per_year
                << " years\n";
      std::cout << "  Total steps: " << step << "\n";
      std::cout << "  Seismic events: " << num_seismic_events << "\n";
   }

   return 0;
}
