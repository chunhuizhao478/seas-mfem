// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// fd_cycle_driver.cpp — fully dynamic BP5 cycle driver, seas_fd_cycle_driver
// (plan rev 3, item j).  Quasi-dynamic between earthquakes, Nitsche
// elastodynamics during them, switched on the slip rate.
//
// Phase P1: the QD regime only ([cycle].switching = false).
// Phase P2: the first BP5-FD event ([cycle].switching = true, ic "scec_fd",
// [cycle.switch].max_dynamic_phases = 1), see RunDynamicPhase.
// Phase P3: the switch back: at the end of a dynamic phase the QD regime
// takes the slip and psi over (plan 6.3) and continues, with a checkpoint
// at the switch.  [cycle.probe] adds the switch probe to the first phase.
// Phase P4: the QD -> dynamic switch at t > 0: when the QD regime's max V
// passes V_on, the dynamic regime takes over from the QD balance at that
// step (U = u_s, U' = 0), so the regimes alternate to t_final
// ([cycle.switch].max_dynamic_phases = 0) or to the end of the N-th
// dynamic phase (N > 0).  A run whose QD balance at t = 0 is not above V_on
// starts in QD.
// The first dynamic phase starts from u = u' = 0 and the QD balance at t = 0:
// 0.015 m/s in the BP5 patch with the scec_fd prestress (eq. 23, dtau = 0),
// V_i = 0.03 m/s with [cycle.prestress].eta_V_init = true (tau0 += eta V_i).
//
//   cfg   = CycleConfig::Load(deck)                  # [cycle] + shared sections
//   mesh  = LoadScaledMesh([mesh].path, [cycle].mesh_scale);  ParMesh
//   qd    = QDPhase(mesh, cfg)                       # domain + geometry (H1a ctor)
//   par   = BuildFaultParams(qd.Table(), [friction.rate_state], [cycle.fields])
//   pre   = BuildPrestress(par, ic set, [cycle.prestress].eta_V_init, .a)
//           (a = "amax": a_max = GlobalMaxA(par), one allreduce)
//   qd.Setup(par, pre)                               # H1a, H1b, kernel, IC
//   out   = FdCycleOutput(cfg.output, qd.LocateStations())
//   restart: Checkpoint read, files truncated to the checkpoint, TS step set
//   switching: dyn = DynPhase(mesh, cfg, table, par, pre)   # one W for the run
//       fresh run: Vmax = qd.Solve(0); Vmax > V_on: rows at t = 0 and
//       RunDynamicPhase(0) (the first event), back to QD at its end
//   loop:
//     qd.Run(monitor):  at the start step and after every accepted step
//         Vmax = global max V;  NaN check
//         out.UpdateEvents;  out.Write (cadence);  out.WriteGlobal
//         every N steps: Checkpoint, and the next step restarts its RK stage
//         switching and Vmax > V_on: stop the QD run at this step
//     stopped for a switch: Vqd = qd.Solve(t); RunDynamicPhase(t)
//         (back to QD at its end, unless it stops the run); next loop
//   end: Checkpoint (if N > 0), out.ForceWrite (seas_bp5_full's final row)
//   snapshots ([cycle.snapshots], cycle/snapshots.hpp): with the station
//   rows, at the QD and dynamic marks, at every switch and at the end
//   --wall-limit S: at the first QD step after S seconds of wall time,
//   a checkpoint and a clean stop (exit 8), for chained jobs
//   [cycle.benchmark] (cycle/bench_output.hpp): the SCEC BP5 outputs besides
//   the station files, from the same steps: source and off-fault rows with
//   every station row, profiles on their own cadence, the event rows and
//   rupture times at each event's end (Track at every step, every fault
//   step under LTS)
//
// Speed-ups, all deck options and off by default (each repeats a trunk
// path; the cycle calls the same functions):
//   QD regime:  [cycle.qd].operators = "matrix" (stage operators, QD speed
//               item 1), [cycle.qd.green] (Green's function, dense or
//               H-matrix, items 3-4; --qd-green-cache / --qd-green-hcache
//               override its files), the options file (rk 5dp, atol)
//   dynamic:    [cycle.lts] mode = "rate2" (local time stepping; the mesh
//               is cut with BuildLtsAwarePartition for both regimes)
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
// Switching on (P2 deck: stop after the first event; P3 deck: back to QD and
// the next events dynamic):
//   mpirun -np N ./seas_fd_cycle_driver --deck config/fd_cycle/bp5_fd_first_event.toml
//          [--output-dir DIR] [--max-dyn-steps N] [--log-every N]
//          [--max-steps N] [--checkpoint-every N] [--restart <prefix>]
//          [--qd-green-cache FILE] [--qd-green-hcache FILE]
// --max-dyn-steps stops the run in a dynamic phase after N steps (timing and
// smoke runs); --end-phase-after N ends every dynamic phase after N steps
// and switches back (smoke runs).  Checkpoints exist only in the QD regime,
// so --restart continues the QD part of a switched run (from a switch back,
// or every N QD steps after it), with the dynamic phases so far counted.
// Exit codes: 0 done (or --max-dyn-steps), 2 no PETSc, 3 non-finite, 4 a
// dynamic phase reached T_max, 7 the guard: the QD balance after a switch
// back is above V_on, 8 stopped at --wall-limit (checkpoint written; a
// --restart from it continues the run and its files).

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
#include "../cycle/snapshots.hpp"
#include "../cycle/bench_output.hpp"

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

/// Where a dynamic phase leaves the run.
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

   const int nc = out.NumStationColumns();
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

/// The dynamic regime's settings, once per run (main builds one DynPhase and
/// every dynamic phase re-enters it).
void PrintDynamicSetup(const MPIContext &mpi, const cycle::CycleConfig &cfg,
                       cycle::DynPhase &dyn, double setup_s)
{
   const cycle::SwitchConfig &sw = cfg.switch_cfg;
   const cycle::DynamicConfig &dc = cfg.dyn;
   const int p = cfg.shared.mesh.order;
   if (!mpi.IsRoot()) { return; }
   std::cout << std::setprecision(6)
             << "\n  dynamic regime (Nitsche, explicit central difference):\n"
             << "    dt = " << dyn.Dt() << " s (cfl " << dc.cfl << " x C(" << p << ") "
             << dyn.StableConstant() << "), set up in " << setup_s << " s\n"
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
   std::cout << "    absorbing boundary: ";
   if (cfg.absorbing.On())
   {
      const cycle::OuterFaceTags &ab = dyn.AbsorbingFaces();
      std::cout << cfg.absorbing.faces << ", " << ab.natural
                << " traction-free faces (attribute " << ab.natural_to << ")";
      if (cfg.absorbing.Walls())
      {
         std::cout << ", " << ab.dirichlet << " wall faces (attribute " << ab.dirichlet_to
                   << ", frozen load)";
      }
      std::cout << "; " << ab.free_surface << " free-surface faces stay traction-free\n";
   }
   else { std::cout << "off (reflecting)\n"; }
   std::cout << "    shared fault faces: " << dyn.NumSharedFaultFaces()
             << "; frame checks: ghost copies " << dyn.GhostFrameMismatch()
             << ", eps_f n_k vs oriented normal " << dyn.NormalMismatch() << "\n"
             << "    switch: V_on " << sw.V_on << ", V_off " << sw.V_off << " m/s, T_dwell "
             << sw.T_dwell << " s, T_min " << sw.T_min << " s, T_max " << sw.T_max
             << " s, max_dynamic_phases " << sw.max_dynamic_phases << "\n";
}

/// One dynamic phase (plan 6.2, 6.3): from the QD balance at t_s (max V =
/// Vstart > V_on, just solved) to the switch back.  Phase 1 at t = 0 is the
/// first BP5-FD event (P2, P3); a later phase starts where the QD regime's
/// max V passed V_on (P4).
///
///   dyn.Enter(Bridge::ToNitsche(qd), t_s)   # U = u_s, U' = 0, (d, V, psi)
///   until the controller ends the phase or T_max:
///      dyn.Step();  rupture time, peak V per DOF (after every fault step)
///      out.UpdateEvents;  out.Write (cadence) and a global row with it
///      (with LTS a Step ends at a sync point, so all of these act there)
///   end of the phase: a station row, the fault table, the switch log, the
///      summary; then (plan 6.3, unless it is the deck's last phase)
///      (s, psi) = Bridge::ToQD(dyn);  Vqd = qd.Reenter(s, psi, t)
///      checkpoint;  Vqd >= V_on: stop (the guard of plan 6.3)
///
/// With [cycle.probe] phase 1 runs on after the deck's switch until every
/// probe pair has switched (or T_max) for the measurement alone: the files
/// stop at the deck's switch, and the switch back uses the state saved
/// there.  At each pair's switch ProbeSwitch writes the jumps.  Later phases
/// have no probe.  Exit codes: 0 at the end of the phase or at
/// --max-dyn-steps, 3 a non-finite QD balance after the switch, 4 the phase
/// reached T_max, 7 the guard.
PhaseEnd RunDynamicPhase(MPIContext &mpi, const cycle::CycleConfig &cfg, cycle::QDPhase &qd,
                         cycle::DynPhase &dyn, cycle::FdCycleOutput &out,
                         cycle::BenchmarkOutput *bench,
                         cycle::Snapshots &snaps, cycle::RegimeController &ctl,
                         const std::string &commit, real_t t_s,
                         real_t Vstart, long max_dyn_steps, long end_phase_after,
                         int log_every, int &rejections_before,
                         const std::function<void(int, real_t)> &checkpoint)
{
   PhaseEnd end;
   const bool root = mpi.IsRoot();
   const cycle::SwitchConfig &sw = cfg.switch_cfg;
   const std::string prefix = cfg.output.dir + "/" + cfg.output.prefix;
   real_t t = t_s;

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
   const int phase = ctl.Phases();
   out.SwitchLog(t, Vstart, "QD", "dynamic", t == 0.0 ? "start" : "V_on");
   snaps.EnterDynamic();

   // ---- the switch probe (plan P3): one controller per pair ----
   std::vector<cycle::RegimeController> probes;
   std::vector<char> probe_done;
   std::unique_ptr<std::ofstream> probe_rows, probe_sum;
   if (cfg.probe.Enabled() && phase == 1)
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
   Vector trupt(no), peakV(no), vabs, vabs_hook;
   trupt = -1.0;
   peakV = 0.0;
   bool phase_over = false;
   auto track = [&](real_t tt, const Vector &va)
   {
      for (int j = 0; j < no; j++)
      {
         if (trupt(j) < 0.0 && va(j) >= V_rupt) { trupt(j) = tt; }
         peakV(j) = std::max(peakV(j), va(j));
      }
      if (bench) { bench->Track(tt, va); }
   };
   dyn.SlipRateMagnitude(vabs);
   track(t, vabs);
   // After every fault step (every step without LTS), until the deck's
   // switch: the files stop there even when the probe runs on.
   struct HookGuard
   {
      cycle::DynPhase &d;
      ~HookGuard() { d.SetFaultStepHook(nullptr); }
   } hook_guard{dyn};
   dyn.SetFaultStepHook([&](real_t tt)
   {
      if (phase_over) { return; }
      dyn.SlipRateMagnitude(vabs_hook);
      track(tt, vabs_hook);
   });

   auto fill_dyn = [&dyn](cycle::StationFields &f) { dyn.FillStationFields(f); };
   auto speed_dyn = [&dyn](Vector &v) { dyn.SlipRateMagnitude(v); };
   cycle::VolumeFieldsFn vol_dyn = [&dyn](const Vector *&u, const Vector *&v)
   {
      u = &dyn.U();
      v = &dyn.Ud();
   };
   if (root)
   {
      if (t == 0.0)
      {
         std::cout << "\n  initial max V (QD balance) = " << Vstart << " m/s; dynamic phase "
                   << phase << " from t = " << t << " s\n\n";
      }
      else
      {
         std::cout << "\n  QD -> dynamic at t = " << t << " s (" << t / BP5Params::seconds_per_year
                   << " yr): max V (QD balance) = " << Vstart << " m/s > V_on; dynamic phase "
                   << phase << "\n\n";
      }
      std::cout << std::setw(10) << "step" << std::setw(14) << "tau [s]"
                << std::setw(14) << "max V" << std::setw(14) << "max VdT/Dc" << std::setw(28)
                << "sigma_n_eff [MPa]" << std::setw(26) << "sigma_PC [MPa]"
                << std::setw(14) << "s/step" << std::setw(28) << "max V at x, y, z [km]"
                << "\n" << std::string(148, '-') << "\n";
   }

   cycle::DynVerdict verdict = cycle::DynVerdict::Stay;
   bool step_cap = false;
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
      if (bench)
      {
         bench->Rows(t, fill_dyn, vol_dyn);
         bench->Flush();
      }
      out.Flush();
      snaps.Dyn(t, dyn.Tau(), dyn.Steps(), phase, fill_dyn, dyn.U(), dyn.Ud(), true);
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
      const std::string fault_file =
         prefix + (phase == 1 ? std::string("_fault.txt")
                              : "_fault_phase" + std::to_string(phase) + ".txt");
      cycle::WriteFaultTable(fault_file, mpi.GetComm(),
                             "fd_cycle fault table at t = " + std::to_string(t) +
                             " s, end of dynamic phase " + std::to_string(phase) + " (from t = " +
                             std::to_string(t_s) + " s; " + reason + "); t_rupture: "
                             "first V >= 1e-3 m/s, -1 = never; commit " + commit,
                             {"x(m)", "y(m)", "z(m)", "t_rupture(s)", "slip_strike(m)",
                              "slip_dip(m)", "peak_V(m/s)", "psi"}, tab);

      // Step cost: rank-0 average and the slowest rank's.
      const double loop_s = MPI_Wtime() - wall0;
      const real_t step_avg = dyn.Steps() > 0 ? step_s / dyn.Steps() : 0.0;
      const real_t step_max = mpi.GlobalMax(step_avg);
      if (root)
      {
         std::cout << std::setprecision(6) << "\n=== fd_cycle dynamic phase " << phase
                   << " summary ===\n"
                   << "  commit: " << commit << "\n"
                   << "  end: " << reason << " at t = " << t << " s (phase length "
                   << dyn.Tau() << " s, " << dyn.Steps() << " steps of " << dyn.StepDt() << " s"
                   << (dyn.Lts() ? ", LTS sync intervals; " + std::to_string(dyn.FaultSteps()) +
                       " fault steps of " + std::to_string(dyn.Dt()) + " s" : std::string())
                   << ")"
                   << (to_qd ? "; back to QD" : "; the run stops") << "\n"
                   << "  max V: at the start (QD balance) " << Vstart << ", after step 1 "
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
                   << "\n  files: " << prefix << "_{<stations>,global,catalog,switch}.txt, "
                   << fault_file << "\n"
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
         const cycle::EventUpdate ev = out.UpdateEvents(t, Vmax);
         if (bench) { bench->OnEvents(ev, t, fill_dyn, speed_dyn); }
         if (out.Write(t, Vmax, fill_dyn))
         {
            out.WriteGlobal(t, Vmax);
            if (bench) { bench->Rows(t, fill_dyn, vol_dyn); }
         }
         if (bench) { bench->Profiles(t, Vmax, fill_dyn); }
         snaps.Dyn(t, dyn.Tau(), dyn.Steps(), phase, fill_dyn, dyn.U(), dyn.Ud());
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
      end.t = t;
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
   double wall_limit = -1.0;
   const char *green_cache = "";
   const char *green_hcache = "";

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
                  "Switching on: stop the run in a dynamic phase after this many of its "
                  "steps (-1 = until the controller ends it).");
   args.AddOption(&end_phase_after, "-epa", "--end-phase-after",
                  "Switching on, smoke tests only: end every dynamic phase after "
                  "this many steps as if the controller had ended it (-1 = off).");
   args.AddOption(&log_every, "-le", "--log-every",
                  "Switching on: a log line every this many dynamic steps.");
   args.AddOption(&wall_limit, "-wl", "--wall-limit",
                  "Stop with a checkpoint (exit 8) at the first QD step after this many "
                  "seconds of wall time since the start (-1 = off; needs checkpoints).");
   args.AddOption(&green_cache, "-gc", "--qd-green-cache",
                  "[cycle.qd.green]: the dense G file (overrides the deck's cache).");
   args.AddOption(&green_hcache, "-ghc", "--qd-green-hcache",
                  "[cycle.qd.green]: the H-matrix file (overrides the deck's hcache).");
   args.Parse();
   if (!args.Good())
   {
      if (root) { args.PrintUsage(std::cout); }
      return 1;
   }

   const double wall_start = MPI_Wtime();
   const std::string commit = cycle::BuildCommit();
   cycle::CycleConfig cfg = cycle::CycleConfig::Load(deck);
   // Run-time overrides: where the files go and how often to checkpoint.
   if (output_dir[0] != '\0') { cfg.output.dir = output_dir; }
   if (checkpoint_every >= 0) { cfg.checkpoint_every = checkpoint_every; }
   if (green_cache[0] != '\0' || green_hcache[0] != '\0')
   {
      MFEM_VERIFY(cfg.qd.green.enabled, "--qd-green-cache / --qd-green-hcache need "
                  "[cycle.qd.green] in the deck");
      if (green_cache[0] != '\0') { cfg.qd.green.cache = green_cache; }
      if (green_hcache[0] != '\0')
      {
         MFEM_VERIFY(cfg.qd.green.compress > 0.0, "--qd-green-hcache needs "
                     "[cycle.qd.green].compress > 0");
         cfg.qd.green.hcache = green_hcache;
      }
   }
   MFEM_VERIFY(wall_limit < 0.0 || cfg.checkpoint_every > 0,
               "--wall-limit needs checkpoints ([cycle.checkpoint].every_qd_steps or "
               "--checkpoint-every > 0): the stop writes one to continue from");
   const std::string restart_prefix = restart;
   const std::string ckpt_prefix = cfg.output.dir + "/" + cfg.output.prefix;

   if (root)
   {
      std::cout << "seas_fd_cycle_driver ("
                << (cfg.switching ? "P4: QD and dynamic regimes, switched on max V"
                                  : "P1: QD regime")
                << ")\n"
                << "  commit: " << commit
                << (cycle::BuildIsClean() ? "" : "   [not a clean build of a commit]") << "\n"
                << "  ranks: " << mpi.Size() << "\n"
                << "  deck: " << deck << "\n"
                << "  IC set: " << cfg.IcLabel()
                << (cfg.prestress.include_eta && !cfg.IsGold()
                    ? " (eta |V_init| in tau0, [cycle.prestress].eta_V_init)" : "") << "\n"
                << "  mesh: " << cfg.shared.mesh.path << " x " << cfg.mesh_scale
                << ", p = " << cfg.shared.mesh.order << "\n"
                << "  QD solver: " << cfg.qd.solver << ", PETSc options: "
                << cfg.qd.petsc_options << " (atol " << cfg.qd.ts_atol << ")\n"
                << "  QD stages: " << cfg.qd.operators;
      if (cfg.qd.green.enabled)
      {
         const cycle::QDGreenDeck &g = cfg.qd.green;
         std::cout << ", Green's function " << (g.compress > 0.0 ? "H-matrix eps " +
                   std::to_string(g.compress) : std::string("dense")) << " (cache '"
                   << g.cache << "', hcache '" << g.hcache << "', check every "
                   << g.check_every << " stages)";
      }
      std::cout << "\n  dynamic stepping: " << (cfg.lts.On() ? "LTS rate2, cluster cap " +
                std::to_string(cfg.lts.cluster_cap) + ", fault buffer " +
                std::to_string(cfg.lts.fault_buffer) : std::string("global (central difference)"))
                << "\n  benchmark output: " << (cfg.bench.enabled ? std::to_string(
                cfg.bench.volume_stations.size()) + " off-fault stations, " +
                std::to_string(cfg.bench.profiles.size()) + " profiles" : std::string("off"))
                << "\n"
                << "  t_final: " << cfg.t_final << " s, dt_init: " << cfg.dt_init << " s\n"
                << "  output: " << ckpt_prefix << "_*\n"
                << "  checkpoint every: " << cfg.checkpoint_every << " QD steps"
                << (cfg.checkpoint_every > 0 ? "" : " (none)") << "\n"
                << "  wall limit: "
                << (wall_limit > 0.0 ? std::to_string(wall_limit) + " s" : std::string("none"))
                << "\n";
      const cycle::SnapshotConfig &sc = cfg.snapshots;
      std::cout << "  snapshots: fault "
                << (sc.fault.enabled ? "dt_qd " + std::to_string(sc.fault.dt_qd) + " s, dt_dyn " +
                                       std::to_string(sc.fault.dt_dyn) + " s"
                                     : std::string("off"))
                << "; surface "
                << (sc.surface.enabled ? "dt_qd " + std::to_string(sc.surface.dt_qd) +
                                         " s, dt_dyn " + std::to_string(sc.surface.dt_dyn) + " s"
                                       : std::string("off"))
                << "\n";
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
      cfg.CheckDimension(smesh->Dimension());   // 2D in-plane rev 3, section 5.5
      if (root)
      {
         std::cout << "  elements: " << smesh->GetNE() << ", vertices: " << smesh->GetNV()
                   << ", boundary faces: " << smesh->GetNBE() << "\n";
      }
      // [cycle.lts] rate2: clusters and the LTS-aware cut, then the ParMesh.
      cycle::LtsLayout lts = cycle::DynPhase::PlanLts(*smesh, cfg, mpi.GetComm());
      ParMesh pmesh(mpi.GetComm(), *smesh, lts.part.empty() ? nullptr : lts.part.data());
      cycle::DynPhase::MapLts(lts, *smesh, pmesh);
      if (root && lts.on) { std::cout << "  LTS: " << lts.report << "\n"; }
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
      // [cycle.prestress].a = "amax" (BP3): a_max over every rank's DOFs.
      cycle::PrestressConfig pcfg = cfg.prestress;
      if (pcfg.AMax())
      {
         pcfg.a_max = cycle::GlobalMaxA(par, mpi.GetComm());
         if (root) { std::cout << "  prestress: a = \"amax\", a_max = " << pcfg.a_max << "\n"; }
      }
      const cycle::Prestress pre = cycle::BuildPrestress(par, pcfg);
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
      std::vector<cycle::SnapshotFileState> snap_restart;
      cycle::BenchmarkOutput::State bench_restart;
      bool have_bench_restart = false;
      const std::vector<cycle::StationLocation> locs = qd.LocateStations(cfg.output);
      qd.InitTS(dt);
      if (!restart_prefix.empty())
      {
         cycle::CheckpointInfo info;
         cycle::CheckpointData data;
         cycle::ReadCycleCheckpoint(restart_prefix, commit, allow_other_commit, info, data, mpi);
         MFEM_VERIFY(info.ic == cfg.IcLabel(), "restart: the checkpoint has IC set "
                     << info.ic << ", the deck " << cfg.IcLabel());
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
         MFEM_VERIFY(info.bench == cfg.bench.enabled, "restart: the checkpoint was written "
                     "with [cycle.benchmark] " << (info.bench ? "on" : "off") << ", the deck "
                     "has it " << (cfg.bench.enabled ? "on" : "off"));
         if (!std::isnan(info.fault_digest))
         {
            const real_t dg = cycle::FaultDigest(qd.Table().coords);
            MFEM_VERIFY(std::abs(dg - info.fault_digest) <= 1e-12 * (std::abs(dg) + 1.0),
                        "restart: rank " << mpi.Rank() << " owns other fault DOFs than when "
                        "the checkpoint was written (digest " << dg << " vs "
                        << info.fault_digest << "): the partition differs");
         }
         out = std::make_unique<cycle::FdCycleOutput>(cfg.output, mpi.GetComm(), locs,
                                                      commit, &info.out);
         bench_restart = info.bench_state;
         have_bench_restart = info.bench;
         snap_restart = info.snapshots;
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
      std::unique_ptr<cycle::BenchmarkOutput> bench;
      if (cfg.bench.enabled)
      {
         const cycle::BenchFaultGeometry bg = qd.BenchGeometry(cfg);
         bench = std::make_unique<cycle::BenchmarkOutput>(
            cfg, mpi.GetComm(), pmesh, qd.Domain().GetFESpace(), bg, commit,
            have_bench_restart ? &bench_restart : nullptr);
         if (root)
         {
            std::cout << "  benchmark output: fault area " << bench->FaultArea()
                      << " m^2, domain A " << bench->AreaA() << " m^2; off-fault stations";
            for (int s = 0; s < bench->NumVolumeStations(); s++)
            {
               std::cout << " " << cfg.bench.volume_stations[s].name << " (rank "
                         << bench->VolumeOwner(s) << ")";
            }
            std::cout << "\n";
         }
      }
      cycle::Snapshots snaps(cfg, mpi.GetComm(), pmesh, qd.Domain().GetFESpace(), qd.Table(),
                             qd.Domain().GetNbfPerFace(), commit,
                             restart_prefix.empty() ? nullptr : &snap_restart);
      if (root && snaps.Any())
      {
         if (cfg.snapshots.fault.enabled)
         {
            std::cout << "  fault snapshots: " << snaps.FaultPath() << " (" << snaps.FaultPoints()
                      << " points, " << snaps.FaultSteps() << " snapshots so far)\n";
         }
         if (cfg.snapshots.surface.enabled)
         {
            std::cout << "  surface snapshots: " << snaps.SurfacePath() << " ("
                      << snaps.SurfacePoints() << " points, " << snaps.SurfaceSteps()
                      << " snapshots so far)\n";
         }
      }

      auto checkpoint = [&](int step, real_t tt)
      {
         cycle::CheckpointInfo info;
         info.t = tt;
         info.dt_next = qd.NextDt();
         info.step = step;
         info.rejections = rejections_before + qd.Rejections();
         info.commit = commit;
         info.ic = cfg.IcLabel();
         info.dyn_phases = ctl.Phases();
         info.out = out->GetState();
         info.snapshots = snaps.GetState();
         info.fault_digest = cycle::FaultDigest(qd.Table().coords);
         info.bench = (bench != nullptr);
         if (bench) { info.bench_state = bench->GetState(); }
         cycle::CheckpointData d;
         d.state = qd.State();
         d.displacement = qd.Displacement();   // solved here with the Green's function
         d.traction = qd.Operator().GetTraction();
         d.slip_rate = qd.Fault().GetSlipRate();
         cycle::WriteCycleCheckpoint(ckpt_prefix, info, d, mpi);
      };

      // One dynamic regime for the whole run; every dynamic phase re-enters it.
      std::unique_ptr<cycle::DynPhase> dyn;
      if (cfg.switching)
      {
         const double w = MPI_Wtime();
         // Phase A (2D in-plane design rev 3): the dynamic regime reads the
         // QD regime's fault orientation (one system for both).
         dyn = std::make_unique<cycle::DynPhase>(pmesh, cfg, qd.Table(), par, pre, 1.0,
                                                 lts.on ? &lts : nullptr,
                                                 qd.Domain().FrameField());
         PrintDynamicSetup(mpi, cfg, *dyn, MPI_Wtime() - w);
      }
      auto dynamic_phase = [&](real_t ts, real_t Vs)
      {
         return RunDynamicPhase(mpi, cfg, qd, *dyn, *out, bench.get(), snaps, ctl, commit, ts, Vs,
                                max_dyn_steps, end_phase_after, log_every, rejections_before,
                                checkpoint);
      };

      // The QD regime's volume fields for the snapshots and the off-fault
      // rows: the velocity is a solve, done at most once per monitored state
      // (qd_vel_ok is cleared whenever the state changes).
      Vector qd_vel;
      bool qd_vel_ok = false;
      cycle::VolumeFieldsFn vol_qd = [&](const Vector *&u, const Vector *&v)
      {
         if (!qd_vel_ok)
         {
            qd.Velocity(qd_vel);
            qd_vel_ok = true;
         }
         u = &qd.Displacement();
         v = &qd_vel;
      };
      auto speed_qd = [&qd](Vector &v) { qd.SlipRateMagnitude(v); };

      // A fresh switched run starts from the QD balance at t = 0: dynamic if
      // its max V is above V_on (the BP5-FD first event), else in QD.  P1
      // and a restart start in QD.
      bool run_qd = true;
      if (cfg.switching && restart_prefix.empty())
      {
         const real_t Vstart = qd.Solve(t);
         if (ctl.Start(Vstart) == cycle::Regime::Dynamic)
         {
            // The rows at t = 0 from the QD fields (as the P1 monitor's start step).
            auto fill0 = [&qd](cycle::StationFields &f) { qd.FillStationFields(f); };
            qd_vel_ok = false;
            const cycle::EventUpdate ev = out->UpdateEvents(t, Vstart);
            if (bench) { bench->OnEvents(ev, t, fill0, speed_qd); }
            const bool row = out->Write(t, Vstart, fill0);
            out->WriteGlobal(t, Vstart);
            if (row && bench) { bench->Rows(t, fill0, vol_qd); }
            if (bench) { bench->Profiles(t, Vstart, fill0); }
            snaps.QD(t, qd.StepNumber(), 0, cycle::SnapReason::Start, fill0, vol_qd);
            const PhaseEnd end = dynamic_phase(t, Vstart);
            rc = end.rc;
            run_qd = end.back_to_qd;
            t = end.t;
            dt = cfg.dt_init;
         }
         else if (root)
         {
            std::cout << "\n  initial max V (QD balance) = " << Vstart << " m/s, not above V_on "
                      << "= " << cfg.switch_cfg.V_on << " m/s: the run starts in the QD regime\n";
         }
      }

      // The QD regime until t_final, --max-steps, a NaN or the next
      // QD -> dynamic switch (P4); after each dynamic phase it continues
      // from the switch back.
      const real_t yr = BP5Params::seconds_per_year;
      auto fill = [&qd](cycle::StationFields &f) { qd.FillStationFields(f); };
      bool diverged = false, to_dynamic = false, qd_ran = false, qd_done = false;
      bool wall_stop = false;
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
         qd_vel_ok = false;
         const cycle::EventUpdate ev = out->UpdateEvents(tt, Vmax);
         if (bench) { bench->OnEvents(ev, tt, fill, speed_qd); }
         const bool row = out->Write(tt, Vmax, fill);
         out->WriteGlobal(tt, Vmax);
         if (row && bench) { bench->Rows(tt, fill, vol_qd); }
         if (bench) { bench->Profiles(tt, Vmax, fill); }
         const bool go = cfg.switching && ctl.GoDynamic(Vmax);
         snaps.QD(tt, step, ctl.Phases(),
                  go ? cycle::SnapReason::ToDynamic : cycle::SnapReason::Cadence, fill,
                  vol_qd);
         if (go)
         {
            // A new event: the dynamic regime takes over at this step (P4).
            // No checkpoint here: a restart from it would switch one step
            // later (its start step is not monitored).
            to_dynamic = true;
            return false;
         }
         const bool ckpt = cfg.checkpoint_every > 0 && step % cfg.checkpoint_every == 0;
         if (ckpt)
         {
            checkpoint(step, tt);
            qd.RestartStage();
         }
         if (wall_limit > 0.0)
         {
            // Rank 0's clock decides, so every rank stops at the same step.
            int over = (MPI_Wtime() - wall_start >= wall_limit) ? 1 : 0;
            MPI_Bcast(&over, 1, MPI_INT, 0, mpi.GetComm());
            if (over)
            {
               if (!ckpt) { checkpoint(step, tt); }
               if (root)
               {
                  std::cout << "  wall limit " << wall_limit << " s reached at step " << step
                            << ", t = " << tt / yr << " yr: checkpoint written, stopping\n";
               }
               wall_stop = true;
               return false;
            }
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
      while (run_qd)
      {
         if (root)
         {
            std::cout << "\n" << std::setw(10) << "Step" << std::setw(16) << "Time [yr]"
                      << std::setw(14) << "dt [s]" << std::setw(16) << "V_max [m/s]"
                      << std::setw(8) << "EQs" << "\n" << std::string(64, '-') << "\n";
         }
         to_dynamic = false;
         qd_ran = true;
         qd.Run(t, dt, cfg.t_final, max_steps, monitor);
         if (diverged) { rc = 3; break; }
         if (wall_stop) { rc = 8; break; }
         if (!to_dynamic) { qd_done = true; break; }   // t_final or --max-steps

         // ---- QD -> dynamic (P4): the QD balance at the accepted step ----
         const real_t Vqd = qd.Solve(t);
         const PhaseEnd end = dynamic_phase(t, Vqd);
         rc = end.rc;
         t = end.t;
         if (!end.back_to_qd) { break; }
         dt = cfg.dt_init;
      }

      if (qd_ran)
      {
         const int steps = qd.StepNumber();
         if (qd_done)
         {
            if (cfg.checkpoint_every > 0) { checkpoint(steps, t); }
            qd_vel_ok = false;
            const real_t Vend = qd.MaxV();
            out->ForceWrite(t, Vend, fill);
            if (bench)
            {
               bench->Rows(t, fill, vol_qd);
               bench->Profiles(t, Vend, fill, true);
            }
            snaps.QD(t, steps, ctl.Phases(), cycle::SnapReason::End, fill, vol_qd);
         }
         out->Flush();
         if (bench) { bench->Flush(); }
         qd.Operator().PrintQDGreenCheckSummary(std::cout);
         if (root)
         {
            // rc 7 (the guard) and 3 (a non-finite balance) stop at the switch back.
            const std::string ph = std::to_string(ctl.Phases());
            const std::string stopped =
               (rc == 8)              ? std::string("stopped at the wall limit")
               : (rc == 7 || rc == 3) ? "stopped at the switch back after dynamic phase " + ph
                                      : "stopped in dynamic phase " + ph;
            std::cout << "\n=== fd_cycle summary ===\n"
                      << "  commit: " << commit << "\n"
                      << "  final time: " << t << " s (" << t / yr << " yr)\n"
                      << "  TS steps: " << steps << "\n"
                      << "  step rejections: " << rejections_before + qd.Rejections() << "\n"
                      << "  dynamic phases: " << ctl.Phases() << "\n"
                      << "  events: " << out->NumEvents() << "\n"
                      << "  snapshots: fault " << snaps.FaultSteps() << ", surface "
                      << snaps.SurfaceSteps() << "\n"
                      << "  status: " << (diverged  ? "DIVERGED"
                                          : qd_done ? "done"
                                                    : stopped + " (exit " +
                                                      std::to_string(rc) + ")")
                      << "\n";
         }
      }
   }
   MFEMFinalizePetsc();
   return rc;
#endif
}
