// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// fd_cycle_driver.cpp — fully dynamic BP5 cycle driver, seas_fd_cycle_driver
// (plan rev 3, item j).  Quasi-dynamic between earthquakes, Nitsche
// elastodynamics during them, switched on the slip rate.
//
// Phase P1: the QD regime only ([cycle].switching = false).
// Phase P2: the first BP5-FD event ([cycle].switching = true, ic "scec_fd",
// [cycle.switch].max_dynamic_phases = 1), see RunFirstDynamicPhase.
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
//
// P2 (switching on):
//   mpirun -np N ./seas_fd_cycle_driver --deck config/fd_cycle/bp5_fd_first_event.toml
//          [--output-dir DIR] [--max-dyn-steps N] [--log-every N]
// --max-dyn-steps stops the dynamic phase early (timing and smoke runs).

#include "mfem.hpp"
#include "../cycle/build_info.hpp"
#include "../cycle/scaled_mesh.hpp"
#include "../cycle/cycle_config.hpp"
#include "../cycle/fault_fields.hpp"
#include "../cycle/prestress.hpp"
#include "../cycle/fd_cycle_output.hpp"
#include "../cycle/checkpoint.hpp"
#include "../cycle/qd_phase.hpp"
#include "../cycle/regime_controller.hpp"
#include "../cycle/state_bridge.hpp"
#include "../cycle/dyn_phase.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>

using namespace mfem;
using namespace mfem::seas;

namespace
{

/// Plan P2: the first BP5-FD event (plan 3.1 without the switch back).
///
///   dyn    = DynPhase(mesh, cfg, table, par, pre)      # W, H2, dt, damping
///   Vmax   = qd.Solve(0)                 # static u and V from the QD balance
///   rows at t = 0 from the QD fields (as the P1 monitor's start step)
///   regime = ctl.Start(Vmax)             # dynamic if Vmax > V_on
///   dyn.Enter(Bridge::ToNitsche(qd), t = 0)   # U = u_s, U' = 0, (d, V, psi)
///   until the controller ends the phase or T_max:
///      dyn.Step();  rupture time, peak V per DOF
///      out.UpdateEvents;  out.Write (cadence) and a global row with it
///   end: a station row, the fault table, the switch log, the summary
///
/// A QD start (Vmax <= V_on at t = 0) needs the QD -> dynamic switch at
/// t > 0, which is P4, so it stops here.  Returns 0 at the end of the phase
/// or at --max-dyn-steps, 4 when the phase reached T_max.
int RunFirstDynamicPhase(MPIContext &mpi, const cycle::CycleConfig &cfg, ParMesh &pmesh,
                         cycle::QDPhase &qd, const cycle::FaultParams &par,
                         const cycle::Prestress &pre, cycle::FdCycleOutput &out,
                         const std::string &commit, long max_dyn_steps, int log_every)
{
   const bool root = mpi.IsRoot();
   const cycle::SwitchConfig &sw = cfg.switch_cfg;
   const cycle::DynamicConfig &dc = cfg.dyn;
   const std::string prefix = cfg.output.dir + "/" + cfg.output.prefix;

   const double wall_setup = MPI_Wtime();
   cycle::DynPhase dyn(pmesh, cfg, qd.Table(), par, pre);
   const double setup_s = MPI_Wtime() - wall_setup;
   const int p = cfg.shared.mesh.order;
   if (root)
   {
      std::cout << std::setprecision(6)
                << "\n  dynamic regime (Nitsche, explicit central difference):\n"
                << "    dt = " << dyn.Dt() << " s (cfl " << dc.cfl << " x C(" << p << ") "
                << cycle::DynPhase::kCeig[p] << "), set up in " << setup_s << " s\n"
                << "    Kelvin-Voigt: ";
      if (dyn.HasKelvinVoigt())
      {
         std::cout << "eta_bar " << dc.kv_eta << ", scale " << dc.kv_scale << " (h_e "
                   << dc.kv_length << "), beta " << dyn.KelvinVoigtBetaMin() << " .. "
                   << dyn.KelvinVoigtBetaMax() << " s, " << dyn.KelvinVoigtElements()
                   << " elements\n";
      }
      else { std::cout << "off\n"; }
      std::cout << "    Prakash-Clifton: ";
      if (dyn.HasPrakashClifton())
      {
         std::cout << "v_star " << dc.pc_v_star << " m/s, length " << dc.pc_length << " m\n";
      }
      else { std::cout << "off\n"; }
      std::cout << "    shared fault faces: " << dyn.NumSharedFaultFaces()
                << "; frame checks: ghost copies " << dyn.GhostFrameMismatch()
                << ", eps_f n_k vs oriented normal " << dyn.NormalMismatch() << "\n"
                << "    switch: V_on " << sw.V_on << ", V_off " << sw.V_off << " m/s, T_dwell "
                << sw.T_dwell << " s, T_min " << sw.T_min << " s, T_max " << sw.T_max << " s\n";
   }

   // ---- start: the QD balance at t = 0 ----
   real_t t = 0.0;
   const real_t Vstart = qd.Solve(t);
   auto fill_qd = [&qd](cycle::StationFields &f) { qd.FillStationFields(f); };
   out.UpdateEvents(t, Vstart);
   out.Write(t, Vstart, fill_qd);
   out.WriteGlobal(t, Vstart);

   std::unique_ptr<std::ofstream> swlog;
   if (root)
   {
      swlog = std::make_unique<std::ofstream>(prefix + "_switch.txt");
      *swlog << "# fd_cycle switch log, commit " << commit << "\n"
             << "# Columns: time(s), Vmax(m/s), from, to, reason\n"
             << std::scientific << std::uppercase << std::setprecision(12);
   }
   cycle::RegimeController ctl(sw);
   if (ctl.Start(Vstart) != cycle::Regime::Dynamic)
   {
      if (root)
      {
         std::cerr << "ERROR: max V = " << Vstart << " m/s at t = 0 is not above V_on = "
                   << sw.V_on << " m/s, so the run would start in the QD regime; the "
                   "QD -> dynamic switch at t > 0 is plan phase P4, not in this build.\n";
      }
      return 5;
   }

   // ---- hand-over (plan 6.2): U = u_s, U' = 0, (d, V, psi), sigma* ----
   Vector sk, psi_qd, sn, d, V, psi;
   qd.Slip(sk);
   qd.Psi(psi_qd);
   qd.SigmaN(sn);
   cycle::bridge::ToNitsche(dyn.Frames(), sk, qd.SlipRate(), psi_qd, d, V, psi);
   dyn.Enter(qd.Displacement(), d, V, psi, sn, t);
   ctl.EnterDynamic();
   if (swlog) { *swlog << t << " " << Vstart << " QD dynamic start\n" << std::flush; }

   const int no = dyn.NumOwnedFaultDOFs();
   const real_t V_rupt = 1.0e-3;   // SCEC: first V >= 1 mm/s (BP5 section 4.6)
   Vector trupt(no), peakV(no), vabs;
   trupt = -1.0;
   peakV = 0.0;
   auto track = [&](real_t tt)
   {
      dyn.SlipRateMagnitude(vabs);
      for (int j = 0; j < no; j++)
      {
         if (trupt(j) < 0.0 && vabs(j) >= V_rupt) { trupt(j) = tt; }
         peakV(j) = std::max(peakV(j), vabs(j));
      }
   };
   track(t);

   auto fill_dyn = [&dyn](cycle::StationFields &f) { dyn.FillStationFields(f); };
   if (root)
   {
      std::cout << "\n  initial max V (QD balance) = " << Vstart << " m/s; dynamic phase from t = "
                << t << " s\n\n" << std::setw(10) << "step" << std::setw(14) << "tau [s]"
                << std::setw(14) << "max V" << std::setw(14) << "max VdT/Dc" << std::setw(28)
                << "sigma_n_eff [MPa]" << std::setw(26) << "sigma_PC [MPa]"
                << std::setw(14) << "s/step" << "\n" << std::string(120, '-') << "\n";
   }

   cycle::DynVerdict verdict = cycle::DynVerdict::Stay;
   bool step_cap = false;
   double step_s = 0.0, wall0 = MPI_Wtime();
   cycle::DynPhase::StepInfo first{}, worst{};
   worst.min_sigma_eff = std::numeric_limits<real_t>::infinity();
   worst.min_sigma_pc = std::numeric_limits<real_t>::infinity();
   real_t Vmax = Vstart;
   while (true)
   {
      const double w = MPI_Wtime();
      dyn.Step();
      step_s += MPI_Wtime() - w;
      t = dyn.Time();
      const cycle::DynPhase::StepInfo r = dyn.LastStepInfo();
      Vmax = r.max_V;
      if (dyn.Steps() == 1) { first = r; }
      worst.max_V = std::max(worst.max_V, r.max_V);
      worst.max_V_dt_over_Dc = std::max(worst.max_V_dt_over_Dc, r.max_V_dt_over_Dc);
      worst.min_sigma_eff = std::min(worst.min_sigma_eff, r.min_sigma_eff);
      worst.max_sigma_eff = std::max(worst.max_sigma_eff, r.max_sigma_eff);
      worst.min_sigma_pc = std::min(worst.min_sigma_pc, r.min_sigma_pc);
      worst.max_sigma_pc = std::max(worst.max_sigma_pc, r.max_sigma_pc);
      track(t);
      out.UpdateEvents(t, Vmax);
      if (out.Write(t, Vmax, fill_dyn)) { out.WriteGlobal(t, Vmax); }
      verdict = ctl.AfterDynStep(Vmax, dyn.Tau(), dyn.Dt());
      step_cap = (max_dyn_steps >= 0 && dyn.Steps() >= max_dyn_steps);
      const bool last = (verdict != cycle::DynVerdict::Stay) || step_cap;
      if (root && (dyn.Steps() % log_every == 0 || dyn.Steps() <= 3 || last))
      {
         std::cout << std::setw(10) << dyn.Steps() << std::scientific << std::setprecision(5)
                   << std::setw(14) << dyn.Tau() << std::setw(14) << r.max_V << std::setw(14)
                   << r.max_V_dt_over_Dc << std::setw(14) << r.min_sigma_eff / 1e6
                   << std::setw(14) << r.max_sigma_eff / 1e6 << std::setw(13)
                   << r.min_sigma_pc / 1e6 << std::setw(13) << r.max_sigma_pc / 1e6
                   << std::setw(14) << step_s / dyn.Steps() << "\n" << std::flush;
      }
      if (last) { break; }
   }
   const double loop_s = MPI_Wtime() - wall0;

   out.ForceWrite(t, Vmax, fill_dyn);
   out.WriteGlobal(t, Vmax);
   out.Flush();
   const char *reason = step_cap && verdict == cycle::DynVerdict::Stay ? "max_dyn_steps"
                        : verdict == cycle::DynVerdict::EndOfPhase  ? "end_of_phase"
                                                                     : "T_max";
   if (swlog) { *swlog << t << " " << Vmax << " dynamic stop " << reason << "\n" << std::flush; }

   // One file for the fault: rupture time, final slip, peak V, final psi.
   {
      Vector sk_end;
      dyn.SlipKernel(sk_end);
      const Vector &X = dyn.FaultPositions();
      const int nc = 8;
      Vector tab(nc * no);
      for (int j = 0; j < no; j++)
      {
         real_t *r = &tab(nc * j);
         r[0] = X(3 * j);
         r[1] = X(3 * j + 1);
         r[2] = X(3 * j + 2);
         r[3] = trupt(j);
         r[4] = sk_end(2 * j + 1);
         r[5] = sk_end(2 * j);
         r[6] = peakV(j);
         r[7] = dyn.Psi()(j);
      }
      cycle::WriteFaultTable(prefix + "_fault.txt", mpi.GetComm(),
                             "fd_cycle fault table at t = " + std::to_string(t) +
                             " s, end of the first dynamic phase (" + reason + "); t_rupture: "
                             "first V >= 1e-3 m/s, -1 = never; commit " + commit,
                             {"x(m)", "y(m)", "z(m)", "t_rupture(s)", "slip_strike(m)",
                              "slip_dip(m)", "peak_V(m/s)", "psi"}, tab);
   }

   // Step cost: rank-0 average and the slowest rank's.
   const real_t step_avg = dyn.Steps() > 0 ? step_s / dyn.Steps() : 0.0;
   const real_t step_max = mpi.GlobalMax(step_avg);
   if (root)
   {
      std::cout << std::setprecision(6) << "\n=== fd_cycle P2 summary ===\n"
                << "  commit: " << commit << "\n"
                << "  stop: " << reason << " at t = " << t << " s (phase length "
                << dyn.Tau() << " s, " << dyn.Steps() << " steps of " << dyn.Dt() << " s)\n"
                << "  max V: at t = 0 (QD balance) " << Vstart << ", after step 1 "
                << first.max_V << ", over the phase " << worst.max_V << " m/s\n"
                << "  max |V| dt / Dc: " << worst.max_V_dt_over_Dc << " (the operator aborts above 2)\n"
                << "  sigma_n_eff range: [" << worst.min_sigma_eff / 1e6 << ", "
                << worst.max_sigma_eff / 1e6 << "] MPa";
      if (dyn.HasPrakashClifton())
      {
         std::cout << "; sigma_PC range: [" << worst.min_sigma_pc / 1e6 << ", "
                   << worst.max_sigma_pc / 1e6 << "] MPa";
      }
      std::cout << "\n  time per dynamic step: " << step_avg << " s (rank 0), " << step_max
                << " s (slowest rank), on " << mpi.Size() << " ranks; loop with output "
                << loop_s << " s\n"
                << "  events: " << out.NumEvents() << (out.InEvent() ? " (still open)" : "")
                << "\n  files: " << prefix << "_{<stations>,global,catalog,switch,fault}.txt\n";
   }
   return verdict == cycle::DynVerdict::TooLong ? 4 : 0;
}

} // namespace

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
   int max_dyn_steps = -1;
   int log_every = 1000;

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
   args.AddOption(&max_dyn_steps, "-nd", "--max-dyn-steps",
                  "Switching on: stop the dynamic phase after this many steps "
                  "(-1 = until the controller ends it).");
   args.AddOption(&log_every, "-le", "--log-every",
                  "Switching on: a log line every this many dynamic steps.");
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
   MFEM_VERIFY(!cfg.switching || restart_prefix.empty(),
               "--restart with [cycle].switching = true: checkpoints come with the "
               "switch back to QD (plan P3); P2 runs one dynamic phase from t = 0");
   const std::string ckpt_prefix = cfg.output.dir + "/" + cfg.output.prefix;

   if (root)
   {
      std::cout << "seas_fd_cycle_driver ("
                << (cfg.switching ? "P2: first dynamic phase" : "P1: QD regime") << ")\n"
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
         const real_t Vmax0 = qd.MaxV();   // collective: every rank calls it
         if (root)
         {
            std::cout << "  fault DOFs: " << qd.NumGlobalFaultDOFs() << " (global)\n"
                      << "  initial V_max: " << Vmax0 << " m/s\n"
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

      if (cfg.switching)
      {
         // P2: the first dynamic phase from t = 0, then stop.
         rc = RunFirstDynamicPhase(mpi, cfg, pmesh, qd, par, pre, *out, commit,
                                   max_dyn_steps, log_every);
      }
      else
      {
         // P1: the QD regime only.
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
   }
   MFEMFinalizePetsc();
   return rc;
#endif
}
