// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// fd_cycle_driver.cpp — fully dynamic BP5 cycle driver, seas_fd_cycle_driver
// (plan rev 3, item j).  Quasi-dynamic between earthquakes, Nitsche
// elastodynamics during them, switched on the slip rate.
//
// Phase P1: the QD regime only ([cycle].switching = false).
//
//   cfg   = CycleConfig::Load(deck)                  # [cycle] + shared sections
//   mesh  = LoadScaledMesh([mesh].path, [cycle].mesh_scale);  ParMesh
//   qd    = QDPhase(mesh, cfg)                       # domain + geometry (H1a ctor)
//   par   = BuildFaultParams(qd.Table(), [friction.rate_state], [cycle.fields])
//   pre   = BuildPrestress(par, ic set)
//   qd.Setup(par, pre)                               # H1a, H1b, kernel, IC
//   out   = FdCycleOutput(cfg.output, qd.LocateStations())
//   restart: Checkpoint read, files truncated to the checkpoint, TS step set
//   qd.Run(monitor):  at the start step and after every accepted step
//       Vmax = global max V;  NaN check
//       out.UpdateEvents;  out.Write (cadence);  out.WriteGlobal
//       every N steps: Checkpoint, and the next step restarts its RK stage
//   end: Checkpoint (if N > 0), out.ForceWrite (seas_bp5_full's final row)
//
// The per-step work and its order follow seas_bp5_full's TS monitor
// (tests/verification/bp5_verification_full.cpp:529-718), so with the gold
// deck the station and global files match seas_bp5_full's.
//
// Usage (from miniapps/seas, built with make -f cycle/cycle.mk):
//   mpirun -np N ./seas_fd_cycle_driver --deck config/fd_cycle/bp5_gold.toml
//          [--max-steps 50] [--output-dir DIR] [--checkpoint-every N]
//          [--restart <checkpoint prefix>] [--allow-other-commit]
// A restart from a checkpoint written at a multiple of N repeats the
// uninterrupted run bit for bit (test_cycle_restart); one from the final
// checkpoint of a run stopped elsewhere (max steps) matches it to round-off.

#include "mfem.hpp"
#include "../cycle/build_info.hpp"
#include "../cycle/scaled_mesh.hpp"
#include "../cycle/cycle_config.hpp"
#include "../cycle/fault_fields.hpp"
#include "../cycle/prestress.hpp"
#include "../cycle/fd_cycle_output.hpp"
#include "../cycle/checkpoint.hpp"
#include "../cycle/qd_phase.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

using namespace mfem;
using namespace mfem::seas;

int main(int argc, char *argv[])
{
   MPIContext mpi(&argc, &argv);
   const bool root = mpi.IsRoot();

   const char *deck = "config/fd_cycle/bp5_gold.toml";
   int max_steps = 10000000;
   const char *restart = "";
   bool allow_other_commit = false;
   const char *output_dir = "";
   int checkpoint_every = -1;

   OptionsParser args(argc, argv);
   args.AddOption(&deck, "-d", "--deck", "Cycle deck (TOML): [cycle] and the shared sections.");
   args.AddOption(&max_steps, "-n", "--max-steps",
                  "Stop after this TS step number (counted from the start of the run "
                  "chain, as seas_bp5_full --max-steps).");
   args.AddOption(&restart, "-r", "--restart", "Checkpoint prefix to continue from.");
   args.AddOption(&allow_other_commit, "-aoc", "--allow-other-commit", "-no-aoc",
                  "--no-allow-other-commit",
                  "Accept a checkpoint written by another commit (needs approval).");
   args.AddOption(&output_dir, "-o", "--output-dir",
                  "Run directory of the output and checkpoint files (overrides "
                  "[output].output_dir).");
   args.AddOption(&checkpoint_every, "-ce", "--checkpoint-every",
                  "Checkpoint every N accepted QD steps, 0 = none (overrides "
                  "[cycle.checkpoint].every_qd_steps; -1 keeps the deck value).");
   args.Parse();
   if (!args.Good())
   {
      if (root) { args.PrintUsage(std::cout); }
      return 1;
   }

   const std::string commit = cycle::BuildCommit();
   cycle::CycleConfig cfg = cycle::CycleConfig::Load(deck);
   // Run-time overrides: where the files go and how often to checkpoint.
   if (output_dir[0] != '\0') { cfg.output.dir = output_dir; }
   if (checkpoint_every >= 0) { cfg.checkpoint_every = checkpoint_every; }
   const std::string restart_prefix = restart;
   const std::string ckpt_prefix = cfg.output.dir + "/" + cfg.output.prefix;

   if (root)
   {
      std::cout << "seas_fd_cycle_driver (P1: QD regime)\n"
                << "  commit: " << commit
                << (cycle::BuildIsClean() ? "" : "   [not a clean build of a commit]") << "\n"
                << "  ranks: " << mpi.Size() << "\n"
                << "  deck: " << deck << "\n"
                << "  IC set: " << cfg.ic << "\n"
                << "  mesh: " << cfg.shared.mesh.path << " x " << cfg.mesh_scale
                << ", p = " << cfg.shared.mesh.order << "\n"
                << "  QD solver: " << cfg.qd.solver << ", PETSc options: "
                << cfg.qd.petsc_options << "\n"
                << "  t_final: " << cfg.t_final << " s, dt_init: " << cfg.dt_init << " s\n"
                << "  output: " << ckpt_prefix << "_*\n"
                << "  checkpoint every: " << cfg.checkpoint_every << " QD steps"
                << (cfg.checkpoint_every > 0 ? "" : " (none)") << "\n";
   }

   if (root) { std::filesystem::create_directories(cfg.output.dir); }
   mpi.Barrier();

#ifndef MFEM_USE_PETSC
   if (root)
   {
      std::cerr << "ERROR: the QD regime steps with PETSc TS (rk 5dp), as seas_bp5_full "
                   "--petsc-ts; this MFEM was built without PETSc (MFEM_USE_PETSC=NO).\n";
   }
   return 2;
#else
   {
      std::ifstream probe(cfg.qd.petsc_options);
      MFEM_VERIFY(probe.good(), "[cycle.qd].petsc_options file " << cfg.qd.petsc_options
                  << " not found (paths are relative to the run directory)");
   }
   MFEMInitializePetsc(&argc, &argv, cfg.qd.petsc_options.c_str(), NULL);
   int rc = 0;
   {
      // ---- mesh ----
      std::unique_ptr<Mesh> smesh = cycle::LoadScaledMesh(cfg.shared.mesh.path, cfg.mesh_scale);
      if (root)
      {
         std::cout << "  elements: " << smesh->GetNE() << ", vertices: " << smesh->GetNV()
                   << ", boundary faces: " << smesh->GetNBE() << "\n";
      }
      ParMesh pmesh(mpi.GetComm(), *smesh);
      smesh.reset();

      // ---- QD regime ----
      cycle::QDPhase qd(pmesh, mpi, cfg);
      const spatial::SpatialFrictionConfig &s = cfg.shared;
      const MaterialField material = MaterialField::MakeConstant(
         s.material_fallback.lambda, s.material_fallback.mu, s.material_fallback.rho);
      // sigma_n_default is already effective: no pore pressure here (Load
      // refuses a non-zero [pore_pressure]).
      const spatial::PorePressureSpec no_pp;
      const cycle::FaultParams par = cycle::BuildFaultParams(
         qd.Table(), *s.rate_state, cfg.fields, material, pmesh, no_pp);
      const cycle::Prestress pre = cycle::BuildPrestress(par, cfg.prestress);
      qd.Setup(par, pre);

      // psi0 of the kernel against BuildPrestress's (a check value).
      {
         const Vector &st = qd.State();
         real_t dmax = 0.0;
         for (int i = 0; i < pre.psi0.Size(); i++)
         {
            dmax = std::max(dmax, std::abs(st(3 * i + 2) - pre.psi0(i)));
         }
         dmax = mpi.GlobalMax(dmax);
         if (root)
         {
            std::cout << "  fault DOFs: " << qd.NumGlobalFaultDOFs() << " (global)\n"
                      << "  initial V_max: " << qd.MaxV() << " m/s\n"
                      << "  max |psi0(kernel) - psi0(BuildPrestress)|: " << dmax << "\n";
         }
      }

      // ---- restart or fresh output ----
      real_t t = 0.0;
      real_t dt = cfg.dt_init;
      int rejections_before = 0;
      std::unique_ptr<cycle::FdCycleOutput> out;
      const std::vector<cycle::StationLocation> locs = qd.LocateStations(cfg.output);
      qd.InitTS(dt);
      if (!restart_prefix.empty())
      {
         cycle::CheckpointInfo info;
         cycle::CheckpointData data;
         cycle::ReadCycleCheckpoint(restart_prefix, commit, allow_other_commit, info, data, mpi);
         MFEM_VERIFY(info.ic == cfg.ic, "restart: the checkpoint has IC set " << info.ic
                     << ", the deck " << cfg.ic);
         MFEM_VERIFY(data.state.Size() == qd.State().Size(),
                     "restart: state size " << data.state.Size() << " vs "
                     << qd.State().Size() << " (same mesh, order and rank count needed)");
         qd.State() = data.state;
         qd.Operator().SetDisplacement(data.displacement);
         qd.Fault().SetSlipRate(data.slip_rate);
         t = info.t;
         dt = info.dt_next;
         rejections_before = info.rejections;
         qd.RestoreTS(info.step);
         out = std::make_unique<cycle::FdCycleOutput>(cfg.output, mpi.GetComm(), locs,
                                                      commit, &info.out);
         if (root)
         {
            std::cout << "  restarted from " << restart_prefix << ": step " << info.step
                      << ", t = " << t << " s, next dt = " << dt << " s, written by "
                      << info.commit << "\n";
         }
      }
      else
      {
         out = std::make_unique<cycle::FdCycleOutput>(cfg.output, mpi.GetComm(), locs, commit);
      }

      auto fill = [&qd](cycle::StationFields &f) { qd.FillStationFields(f); };
      auto checkpoint = [&](int step, real_t tt)
      {
         cycle::CheckpointInfo info;
         info.t = tt;
         info.dt_next = qd.NextDt();
         info.step = step;
         info.rejections = rejections_before + qd.Rejections();
         info.commit = commit;
         info.ic = cfg.ic;
         info.out = out->GetState();
         cycle::CheckpointData d;
         d.state = qd.State();
         d.displacement = qd.Operator().GetDisplacement();
         d.traction = qd.Operator().GetTraction();
         d.slip_rate = qd.Fault().GetSlipRate();
         cycle::WriteCycleCheckpoint(ckpt_prefix, info, d, mpi);
      };

      bool diverged = false;
      if (root)
      {
         std::cout << "\n" << std::setw(10) << "Step" << std::setw(16) << "Time [yr]"
                   << std::setw(14) << "dt [s]" << std::setw(16) << "V_max [m/s]"
                   << std::setw(8) << "EQs" << "\n" << std::string(64, '-') << "\n";
      }
      const real_t yr = BP5Params::seconds_per_year;
      auto monitor = [&](int step, real_t tt) -> bool
      {
         const real_t Vmax = qd.MaxV();
         bool bad = !std::isfinite(Vmax);
         const Vector &st = qd.State();
         for (int i = 0; i < std::min(st.Size(), 100) && !bad; i++)
         {
            if (!std::isfinite(st(i))) { bad = true; }
         }
         if (mpi.GlobalSumInt(bad ? 1 : 0) > 0)
         {
            if (root)
            {
               std::cerr << "NaN/Inf at step " << step << ", t = " << tt / yr
                         << " yr, V_max = " << Vmax << "\n";
            }
            diverged = true;
            return false;
         }
         out->UpdateEvents(tt, Vmax);
         out->Write(tt, Vmax, fill);
         out->WriteGlobal(tt, Vmax);
         if (cfg.checkpoint_every > 0 && step % cfg.checkpoint_every == 0)
         {
            checkpoint(step, tt);
            qd.RestartStage();
         }
         if (root && (step % 10 == 0 || Vmax > cfg.output.event_V))
         {
            std::cout << std::setw(10) << step << std::setw(16) << std::scientific
                      << std::setprecision(6) << tt / yr << std::setw(14)
                      << std::setprecision(3) << qd.NextDt() << std::setw(16)
                      << std::setprecision(3) << Vmax << std::setw(8) << out->NumEvents()
                      << "\n" << std::flush;
         }
         return true;
      };

      qd.Run(t, dt, cfg.t_final, max_steps, monitor);
      const int steps = qd.StepNumber();

      if (diverged)
      {
         rc = 3;
      }
      else
      {
         if (cfg.checkpoint_every > 0) { checkpoint(steps, t); }
         out->ForceWrite(t, qd.MaxV(), fill);
      }
      out->Flush();
      if (root)
      {
         std::cout << "\n=== fd_cycle summary ===\n"
                   << "  commit: " << commit << "\n"
                   << "  final time: " << t << " s (" << t / yr << " yr)\n"
                   << "  TS steps: " << steps << "\n"
                   << "  step rejections: " << rejections_before + qd.Rejections() << "\n"
                   << "  events: " << out->NumEvents() << "\n"
                   << "  status: " << (diverged ? "DIVERGED" : "done") << "\n";
      }
   }
   MFEMFinalizePetsc();
   return rc;
#endif
}
