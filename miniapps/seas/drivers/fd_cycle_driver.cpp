// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// fd_cycle_driver.cpp — fully dynamic BP5 cycle driver, seas_fd_cycle_driver
// (plan rev 3, item j).  Quasi-dynamic between earthquakes, Nitsche
// elastodynamics during them, switched on the slip rate.
//
// Phase P1: the QD regime only ([cycle].switching = false).
// Phase P2: the first BP5-FD event ([cycle].switching = true, ic "scec_fd",
// [cycle.switch].max_dynamic_phases = 1), see RunFirstDynamicPhase.
// Phase P3: the switch back ([cycle.switch].max_dynamic_phases = 0): at the
// end of the first dynamic phase the QD regime takes the slip and psi over
// (plan 6.3) and runs to t_final, with a checkpoint at the switch.  A QD
// max V above V_on after that is the QD -> dynamic switch of P4, which
// stops the run (exit 6).  [cycle.probe] adds the switch probe.
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
// Switching on (P2 deck: stop after the first event; P3 deck: back to QD):
//   mpirun -np N ./seas_fd_cycle_driver --deck config/fd_cycle/bp5_fd_first_event.toml
//          [--output-dir DIR] [--max-dyn-steps N] [--log-every N]
//          [--max-steps N] [--checkpoint-every N] [--restart <prefix>]
// --max-dyn-steps stops the dynamic phase early (timing and smoke runs);
// --end-phase-after N ends it after N steps and switches back (smoke runs).
// Checkpoints exist only in the QD regime, so --restart continues the QD
// part of a switched run (from the switch, or every N QD steps after it).
// Exit codes: 0 done, 2 no PETSc, 3 non-finite, 4 a dynamic phase reached
// T_max, 5 max V <= V_on at t = 0, 6 QD -> dynamic at t > 0 (P4), 7 the
// guard: the QD balance after the switch back is above V_on.

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
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>

using namespace mfem;
using namespace mfem::seas;

namespace
{

/// Where the first dynamic phase leaves the run.
struct PhaseEnd
{
   int rc = 0;               ///< non-zero: stop with this exit code
   bool back_to_qd = false;  ///< the QD regime continues from t
   real_t t = 0.0;
};

/// The global max of `v` over all ranks and the position (3 values) of a
/// DOF that attains it.  Collective.
real_t GlobalMaxAt(const MPIContext &mpi, real_t v, const real_t *x_local, real_t x[3])
{
   struct { double v; int r; } in{static_cast<double>(v), mpi.Rank()}, out{0.0, 0};
   MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, mpi.GetComm());
   double buf[3] = {0.0, 0.0, 0.0};
   if (mpi.Rank() == out.r && x_local)
   {
      for (int d = 0; d < 3; d++) { buf[d] = x_local[d]; }
   }
   MPI_Bcast(buf, 3, MPI_DOUBLE, out.r, mpi.GetComm());
   for (int d = 0; d < 3; d++) { x[d] = static_cast<real_t>(buf[d]); }
   return static_cast<real_t>(out.v);
}

/// The switch probe of plan P3.  At the switch time of one (V_off, T_dwell)
/// pair: the QD balance from the dynamic regime's (slip, psi) at t
/// (Bridge::ToQD, one QD solve), and the station rows of both regimes from
/// the same state.  The root writes one row per station to `rows` and one
/// summary row to `summary`.  Collective; leaves qd.State() at that state.
void ProbeSwitch(const MPIContext &mpi, const cycle::CycleConfig &cfg, cycle::QDPhase &qd,
                 cycle::DynPhase &dyn, cycle::FdCycleOutput &out, real_t V_off,
                 real_t T_dwell, real_t Vmax_dyn, std::ostream *rows, std::ostream *summary)
{
   const real_t t = dyn.Time();
   auto fill_dyn = [&dyn](cycle::StationFields &f) { dyn.FillStationFields(f); };
   auto fill_qd = [&qd](cycle::StationFields &f) { qd.FillStationFields(f); };
   const std::vector<real_t> rd = out.Sample(t, fill_dyn);
   cycle::StationFields fd, fq;
   dyn.FillStationFields(fd);

   Vector sk, psi;
   cycle::bridge::ToQD(dyn.Frames(), dyn.Slip(), dyn.Psi(), sk, psi);
   qd.SetFaultState(sk, psi);
   const real_t Vmax_qd = qd.Solve(t);
   const std::vector<real_t> rq = out.Sample(t, fill_qd);
   qd.FillStationFields(fq);

   // Fault-wide: the largest change of the shear traction vector, with the
   // station formula tau_c = -(tau_a_c + tau_b_c + eta V_c) of each regime.
   const int N = dyn.NumOwnedFaultDOFs();
   real_t dmax = 0.0;
   int jmax = -1;
   for (int j = 0; j < N; j++)
   {
      const real_t dd = -(fq.tau_a_dip(j) + fq.tau_b_dip(j) + fq.eta * fq.V_dip(j))
                        + (fd.tau_a_dip(j) + fd.tau_b_dip(j) + fd.eta * fd.V_dip(j));
      const real_t ds = -(fq.tau_a_strike(j) + fq.tau_b_strike(j) + fq.eta * fq.V_strike(j))
                        + (fd.tau_a_strike(j) + fd.tau_b_strike(j) + fd.eta * fd.V_strike(j));
      const real_t d = std::sqrt(dd * dd + ds * ds);
      if (d > dmax) { dmax = d; jmax = j; }
   }
   real_t xmax[3];
   const real_t dtau = GlobalMaxAt(mpi, dmax, jmax >= 0 ? &dyn.FaultPositions()(3 * jmax)
                                                         : nullptr, xmax);
   if (!mpi.IsRoot()) { return; }

   const int nc = cycle::FdCycleOutput::kStationColumns;
   const std::vector<cycle::StationSpec> &st = cfg.output.stations;
   for (size_t s = 0; s < st.size(); s++)
   {
      const real_t *a = &rd[nc * s], *b = &rq[nc * s];
      *rows << V_off << " " << T_dwell << " " << t << " " << dyn.Tau() << " " << dyn.Steps()
            << " " << st[s].name << " " << a[5] << " " << b[5] << " " << a[6] << " " << b[6]
            << " " << a[3] << " " << b[3] << " " << a[4] << " " << b[4] << "\n";
   }
   rows->flush();
   *summary << V_off << " " << T_dwell << " " << t << " " << dyn.Tau() << " " << dyn.Steps()
            << " " << Vmax_dyn << " " << Vmax_qd << " " << dtau / 1e6 << " " << xmax[0] << " "
            << xmax[1] << " " << xmax[2] << "\n" << std::flush;
}

/// Plan P2 and P3: the first BP5-FD event and the switch back to QD.
///
///   dyn    = DynPhase(mesh, cfg, table, par, pre)      # W, H2, dt, damping
///   Vmax   = qd.Solve(0)                 # static u and V from the QD balance
///   rows at t = 0 from the QD fields (as the P1 monitor's start step)
///   regime = ctl.Start(Vmax)             # dynamic if Vmax > V_on
///   dyn.Enter(Bridge::ToNitsche(qd), t = 0)   # U = u_s, U' = 0, (d, V, psi)
///   until the controller ends the phase or T_max:
///      dyn.Step();  rupture time, peak V per DOF
///      out.UpdateEvents;  out.Write (cadence) and a global row with it
///   end of the phase: a station row, the fault table, the switch log, the
///      summary; then (plan 6.3, unless max_dynamic_phases = 1)
///      (s, psi) = Bridge::ToQD(dyn);  Vqd = qd.Reenter(s, psi, t)
///      checkpoint;  Vqd >= V_on: stop (the guard of plan 6.3)
///
/// With [cycle.probe] the phase runs on after the deck's switch until every
/// probe pair has switched (or T_max) for the measurement alone: the files
/// stop at the deck's switch, and the switch back uses the state saved
/// there.  At each pair's switch ProbeSwitch writes the jumps.
///
/// A QD start (Vmax <= V_on at t = 0) needs the QD -> dynamic switch at
/// t > 0, which is P4, so it stops here.  Exit codes: 0 at the end of the
/// phase or at --max-dyn-steps, 3 a non-finite QD balance after the switch,
/// 4 the phase reached T_max, 5 a QD start, 7 the guard.
PhaseEnd RunFirstDynamicPhase(MPIContext &mpi, const cycle::CycleConfig &cfg,
                              ParMesh &pmesh, cycle::QDPhase &qd,
                              const cycle::FaultParams &par, const cycle::Prestress &pre,
                              cycle::FdCycleOutput &out, cycle::RegimeController &ctl,
                              const std::string &commit, long max_dyn_steps,
                              long end_phase_after, int log_every, int &rejections_before,
                              const std::function<void(int, real_t)> &checkpoint)
{
   PhaseEnd end;
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
                << sw.T_dwell << " s, T_min " << sw.T_min << " s, T_max " << sw.T_max
                << " s, max_dynamic_phases " << sw.max_dynamic_phases << "\n";
   }

   // ---- start: the QD balance at t = 0 ----
   real_t t = 0.0;
   const real_t Vstart = qd.Solve(t);
   auto fill_qd = [&qd](cycle::StationFields &f) { qd.FillStationFields(f); };
   out.UpdateEvents(t, Vstart);
   out.Write(t, Vstart, fill_qd);
   out.WriteGlobal(t, Vstart);

   if (ctl.Start(Vstart) != cycle::Regime::Dynamic)
   {
      if (root)
      {
         std::cerr << "ERROR: max V = " << Vstart << " m/s at t = 0 is not above V_on = "
                   << sw.V_on << " m/s, so the run would start in the QD regime; the "
                   "QD -> dynamic switch at t > 0 is plan phase P4, not in this build.\n";
      }
      end.rc = 5;
      return end;
   }

   // ---- hand-over (plan 6.2): U = u_s, U' = 0, (d, V, psi), sigma* ----
   {
      Vector sk, psi_qd, sn, d, V, psi;
      qd.Slip(sk);
      qd.Psi(psi_qd);
      qd.SigmaN(sn);
      cycle::bridge::ToNitsche(dyn.Frames(), sk, qd.SlipRate(), psi_qd, d, V, psi);
      dyn.Enter(qd.Displacement(), d, V, psi, sn, t);
   }
   ctl.EnterDynamic(t);
   out.SwitchLog(t, Vstart, "QD", "dynamic", "start");

   // ---- the switch probe (plan P3): one controller per pair ----
   std::vector<cycle::RegimeController> probes;
   std::vector<char> probe_done;
   std::unique_ptr<std::ofstream> probe_rows, probe_sum;
   if (cfg.probe.Enabled())
   {
      for (real_t voff : cfg.probe.V_off)
      {
         for (real_t tdw : cfg.probe.T_dwell)
         {
            cycle::SwitchConfig s = sw;
            s.V_off = voff;
            s.T_dwell = tdw;
            s.max_dynamic_phases = 0;
            probes.emplace_back(s);
            probes.back().EnterDynamic(t);
         }
      }
      probe_done.assign(probes.size(), 0);
      if (root)
      {
         const std::string head = "fd_cycle switch probe (plan P3), commit " + commit +
                                  "; deck pair V_off " + std::to_string(sw.V_off) +
                                  " m/s, T_dwell " + std::to_string(sw.T_dwell) + " s, V_on " +
                                  std::to_string(sw.V_on) + " m/s";
         probe_rows = std::make_unique<std::ofstream>(prefix + "_switch_probe.txt");
         *probe_rows << "# " << head << "\n# At each (V_off, T_dwell) pair's switch time: "
                     "the dynamic regime's station values (end of step: X* + gamma V) and "
                     "the QD balance's from the same slip, psi and t.  tau in MPa, V in "
                     "m/s, station formulas of the station files.\n# Columns: V_off(m/s), "
                     "T_dwell(s), t_switch(s), phase_length(s), dyn_steps, station, "
                     "tau_s_dyn, tau_s_qd, tau_d_dyn, tau_d_qd, log10V_s_dyn, log10V_s_qd, "
                     "log10V_d_dyn, log10V_d_qd\n"
                     << std::scientific << std::uppercase << std::setprecision(12);
         probe_sum = std::make_unique<std::ofstream>(prefix + "_switch_probe_summary.txt");
         *probe_sum << "# " << head << "\n# One row per pair, in switch order; a pair that "
                    "did not switch before T_max or --max-dyn-steps has t_switch = -1 and "
                    "the phase time it reached.  dtau_max: the largest change of the shear "
                    "traction vector over the fault [MPa], at (x, y, z).\n# Columns: "
                    "V_off(m/s), T_dwell(s), t_switch(s), phase_length(s), dyn_steps, "
                    "Vmax_dyn(m/s), Vmax_qd(m/s), dtau_max(MPa), x(m), y(m), z(m)\n"
                    << std::scientific << std::uppercase << std::setprecision(12);
      }
   }

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
                << std::setw(14) << "s/step" << std::setw(28) << "max V at x, y, z [km]"
                << "\n" << std::string(148, '-') << "\n";
   }

   cycle::DynVerdict verdict = cycle::DynVerdict::Stay;
   bool step_cap = false, phase_over = false;
   double step_s = 0.0, wall0 = MPI_Wtime();
   cycle::DynPhase::StepInfo first{}, worst{};
   worst.min_sigma_eff = std::numeric_limits<real_t>::infinity();
   worst.min_sigma_pc = std::numeric_limits<real_t>::infinity();
   real_t Vmax = Vstart;
   // The deck pair's switch state (plan 6.3 hands over slip and psi only).
   Vector snap_sk, snap_psi;
   real_t snap_Vmax = 0.0;

   // End of the phase for the files: the last station row, the switch log,
   // the fault table and the summary.  Runs once, at the deck pair's end
   // (or at T_max, --max-dyn-steps).
   auto finish_phase = [&](const char *reason, bool to_qd)
   {
      out.ForceWrite(t, Vmax, fill_dyn);
      out.WriteGlobal(t, Vmax);
      out.Flush();
      out.SwitchLog(t, Vmax, "dynamic", to_qd ? "QD" : "stop", reason);

      // One file for the fault: rupture time, final slip, peak V, final psi.
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

      // Step cost: rank-0 average and the slowest rank's.
      const double loop_s = MPI_Wtime() - wall0;
      const real_t step_avg = dyn.Steps() > 0 ? step_s / dyn.Steps() : 0.0;
      const real_t step_max = mpi.GlobalMax(step_avg);
      if (root)
      {
         std::cout << std::setprecision(6) << "\n=== fd_cycle dynamic phase 1 summary ===\n"
                   << "  commit: " << commit << "\n"
                   << "  end: " << reason << " at t = " << t << " s (phase length "
                   << dyn.Tau() << " s, " << dyn.Steps() << " steps of " << dyn.Dt() << " s)"
                   << (to_qd ? "; back to QD" : "; the run stops") << "\n"
                   << "  max V: at t = 0 (QD balance) " << Vstart << ", after step 1 "
                   << first.max_V << ", over the phase " << worst.max_V << " m/s\n"
                   << "  max |V| dt / Dc: " << worst.max_V_dt_over_Dc
                   << " (the operator aborts above 2)\n"
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
                   << "\n  files: " << prefix << "_{<stations>,global,catalog,switch,fault}.txt\n"
                   << std::flush;
      }
      phase_over = true;
   };

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
      if (!phase_over)
      {
         track(t);
         out.UpdateEvents(t, Vmax);
         if (out.Write(t, Vmax, fill_dyn)) { out.WriteGlobal(t, Vmax); }
         verdict = ctl.AfterDynStep(Vmax, t, dyn.Tau());
      }
      for (size_t k = 0; k < probes.size(); k++)
      {
         if (probe_done[k]) { continue; }
         const cycle::DynVerdict v = probes[k].AfterDynStep(Vmax, t, dyn.Tau());
         if (v == cycle::DynVerdict::EndOfPhase)
         {
            ProbeSwitch(mpi, cfg, qd, dyn, out, probes[k].Config().V_off,
                        probes[k].Config().T_dwell, Vmax, probe_rows.get(), probe_sum.get());
            probe_done[k] = 1;
         }
         else if (v == cycle::DynVerdict::TooLong) { probe_done[k] = 2; }
      }
      step_cap = (max_dyn_steps >= 0 && dyn.Steps() >= max_dyn_steps);
      const bool probes_open = std::find(probe_done.begin(), probe_done.end(), 0)
                               != probe_done.end();
      const bool forced = !phase_over && end_phase_after >= 0 &&
                          dyn.Steps() >= end_phase_after && verdict == cycle::DynVerdict::Stay;
      if (!phase_over && (verdict == cycle::DynVerdict::EndOfPhase || forced))
      {
         const bool to_qd = !ctl.LastPhase();
         verdict = cycle::DynVerdict::EndOfPhase;
         finish_phase(forced ? "end_phase_after" : "end_of_phase", to_qd);
         cycle::bridge::ToQD(dyn.Frames(), dyn.Slip(), dyn.Psi(), snap_sk, snap_psi);
         snap_Vmax = Vmax;
         end.t = t;
         end.back_to_qd = to_qd;
      }
      const bool last = (phase_over && !probes_open) || step_cap ||
                        (!phase_over && verdict == cycle::DynVerdict::TooLong);
      if (dyn.Steps() % log_every == 0 || dyn.Steps() <= 3 || last)
      {
         // Where max V is (collective, so outside the root-only print).
         dyn.SlipRateMagnitude(vabs);
         int jv = -1;
         for (int j = 0; j < no; j++) { if (jv < 0 || vabs(j) > vabs(jv)) { jv = j; } }
         real_t xv[3];
         GlobalMaxAt(mpi, jv >= 0 ? vabs(jv) : 0.0,
                     jv >= 0 ? &dyn.FaultPositions()(3 * jv) : nullptr, xv);
         if (root)
         {
            std::cout << std::setw(10) << dyn.Steps() << std::scientific << std::setprecision(5)
                      << std::setw(14) << dyn.Tau() << std::setw(14) << r.max_V << std::setw(14)
                      << r.max_V_dt_over_Dc << std::setw(14) << r.min_sigma_eff / 1e6
                      << std::setw(14) << r.max_sigma_eff / 1e6 << std::setw(13)
                      << r.min_sigma_pc / 1e6 << std::setw(13) << r.max_sigma_pc / 1e6
                      << std::setw(14) << step_s / dyn.Steps() << std::fixed
                      << std::setprecision(2) << std::setw(10) << xv[0] / 1e3 << std::setw(9)
                      << xv[1] / 1e3 << std::setw(9) << xv[2] / 1e3 << "\n" << std::flush;
         }
      }
      if (last) { break; }
   }

   // Probe pairs that never switched.
   if (probe_sum)
   {
      for (size_t k = 0; k < probes.size(); k++)
      {
         if (probe_done[k] == 1) { continue; }
         *probe_sum << probes[k].Config().V_off << " " << probes[k].Config().T_dwell
                    << " -1 " << dyn.Tau() << " " << dyn.Steps()
                    << " 0 0 0 0 0 0\n";
      }
      probe_sum->flush();
   }
   if (root && !probes.empty())
   {
      const long n_ok = std::count(probe_done.begin(), probe_done.end(), 1);
      std::cout << "\n  switch probe: " << n_ok << " of " << probes.size()
                << " pairs switched by tau = " << dyn.Tau() << " s; files " << prefix
                << "_switch_probe{,_summary}.txt\n";
   }

   if (!phase_over)
   {
      // T_max or --max-dyn-steps before the deck pair's end: stop here.
      const bool too_long = (verdict == cycle::DynVerdict::TooLong);
      finish_phase(too_long ? "T_max" : "max_dyn_steps", false);
      end.rc = too_long ? 4 : 0;
      return end;
   }
   if (!end.back_to_qd) { return end; }   // max_dynamic_phases reached

   // ---- switch back (plan 6.3): (s, psi) at t = t0 + tau, one QD solve, new TS ----
   const int step0 = qd.StepNumber();
   rejections_before += qd.Rejections();
   const real_t Vqd = qd.Reenter(snap_sk, snap_psi, end.t, step0, cfg.dt_init);
   if (root)
   {
      std::cout << std::setprecision(6) << "\n  back to QD at t = " << end.t
                << " s (TS step " << step0 << ", first dt " << cfg.dt_init
                << " s): max V dynamic " << snap_Vmax << " m/s, QD balance " << Vqd
                << " m/s\n" << std::flush;
   }
   if (!std::isfinite(Vqd))
   {
      if (root) { std::cerr << "ERROR: the QD balance after the switch back is not finite.\n"; }
      out.SwitchLog(end.t, Vqd, "QD", "stop", "non_finite");
      end.rc = 3;
      end.back_to_qd = false;
      return end;
   }
   if (cfg.checkpoint_every > 0) { checkpoint(step0, end.t); }
   if (Vqd >= sw.V_on)
   {
      if (root)
      {
         std::cerr << "STOP: max V of the QD balance after the switch back, " << Vqd
                   << " m/s, is not below V_on = " << sw.V_on << " m/s (the guard of plan "
                   "6.3; a chatter policy is decided after the P3 data).\n";
      }
      out.SwitchLog(end.t, Vqd, "QD", "stop", "guard_V_on");
      end.rc = 7;
      end.back_to_qd = false;
      return end;
   }
   return end;
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
   int end_phase_after = -1;
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
   args.AddOption(&end_phase_after, "-epa", "--end-phase-after",
                  "Switching on, smoke tests only: end the first dynamic phase after "
                  "this many steps as if the controller had ended it (-1 = off).");
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
   const std::string ckpt_prefix = cfg.output.dir + "/" + cfg.output.prefix;

   if (root)
   {
      std::cout << "seas_fd_cycle_driver ("
                << (cfg.switching ? "P3: the first event dynamic, then QD" : "P1: QD regime")
                << ")\n"
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
      cycle::RegimeController ctl(cfg.switch_cfg);
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
         ctl.Restore(info.dyn_phases);
         qd.RestoreTS(info.step);
         out = std::make_unique<cycle::FdCycleOutput>(cfg.output, mpi.GetComm(), locs,
                                                      commit, &info.out);
         if (root)
         {
            std::cout << "  restarted from " << restart_prefix << ": step " << info.step
                      << ", t = " << t << " s, next dt = " << dt << " s, dynamic phases "
                      << info.dyn_phases << ", written by " << info.commit << "\n";
         }
      }
      else
      {
         out = std::make_unique<cycle::FdCycleOutput>(cfg.output, mpi.GetComm(), locs, commit);
      }
      if (cfg.switching) { out->OpenSwitchLog(commit); }

      auto checkpoint = [&](int step, real_t tt)
      {
         cycle::CheckpointInfo info;
         info.t = tt;
         info.dt_next = qd.NextDt();
         info.step = step;
         info.rejections = rejections_before + qd.Rejections();
         info.commit = commit;
         info.ic = cfg.ic;
         info.dyn_phases = ctl.Phases();
         info.out = out->GetState();
         cycle::CheckpointData d;
         d.state = qd.State();
         d.displacement = qd.Operator().GetDisplacement();
         d.traction = qd.Operator().GetTraction();
         d.slip_rate = qd.Fault().GetSlipRate();
         cycle::WriteCycleCheckpoint(ckpt_prefix, info, d, mpi);
      };

      // P1, a restart, or after the switch back: the QD regime.
      bool run_qd = !cfg.switching || !restart_prefix.empty();
      if (cfg.switching && restart_prefix.empty())
      {
         // P2 and P3: the first event dynamic from t = 0, then back to QD.
         const PhaseEnd end = RunFirstDynamicPhase(mpi, cfg, pmesh, qd, par, pre, *out, ctl,
                                                   commit, max_dyn_steps, end_phase_after,
                                                   log_every, rejections_before, checkpoint);
         rc = end.rc;
         run_qd = end.back_to_qd;
         t = end.t;
         dt = cfg.dt_init;
      }
      if (run_qd)
      {
         auto fill = [&qd](cycle::StationFields &f) { qd.FillStationFields(f); };
         bool diverged = false, to_dynamic = false;
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
            if (cfg.switching && ctl.GoDynamic(Vmax))
            {
               // A new event: the QD -> dynamic switch at t > 0 is P4.
               if (root)
               {
                  std::cerr << "STOP: max V = " << Vmax << " m/s > V_on = "
                            << cfg.switch_cfg.V_on << " m/s at step " << step << ", t = "
                            << tt / yr << " yr: the QD -> dynamic switch at t > 0 is plan "
                            "phase P4, not in this build.\n";
               }
               out->SwitchLog(tt, Vmax, "QD", "stop", "QD_to_dynamic_is_P4");
               to_dynamic = true;
               return false;
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
         else if (to_dynamic)
         {
            rc = 6;
            if (cfg.checkpoint_every > 0) { checkpoint(steps, t); }
            out->ForceWrite(t, qd.MaxV(), fill);
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
                      << "  status: " << (diverged     ? "DIVERGED"
                                          : to_dynamic ? "stopped at a QD -> dynamic switch (P4)"
                                                       : "done") << "\n";
         }
      }
   }
   MFEMFinalizePetsc();
   return rc;
#endif
}
