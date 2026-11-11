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
//           [stress] kind = "constant_tensor": BuildTensorPrestress(par,
//           qd.ProjectTensor(tensor), ...) instead (cycle/prestress.hpp)
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
//               [cycle.quench] on (quench plan rev1, Chunhui Zhao 2026-11-06):
//               once a phase's max V < V_q the damping C_q = d M_rho goes on
//               (DynPhase::SetQuench), off again above V_q or at the phase
//               end; "quench_on" / "quench_off" lines in the switch log
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
//          [--restart-drop-relaxation]
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
#include "../cycle/pml_partition.hpp"
#include "../cycle/snapshots.hpp"
#include "../cycle/bench_output.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
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

/// [cycle.output].handback_table: the jump of plan 6.3 at one switch back,
/// on every fault DOF.  fd: the dynamic regime's station fields at its
/// switch; qd: just re-entered at the same slip, psi and t (Reenter's solve
/// is its last evaluation).  The root writes <prefix>_handback_phase<N>.txt
/// (rank 0's owned DOFs first, then rank 1's, ...) and prints where the
/// shear traction, sigma_n' and V of the two regimes differ most.  tau is
/// the station files' -(tau_a + tau_b + eta V) of each regime.  Collective.
void WriteHandbackTable(const MPIContext &mpi, const cycle::CycleConfig &cfg,
                        cycle::QDPhase &qd, const cycle::StationFields &fd, const Vector &X,
                        int phase, real_t t, const std::string &commit)
{
   cycle::StationFields fq;
   qd.FillStationFields(fq);
   const int N = fd.V_strike.Size();
   MFEM_VERIFY(fq.V_strike.Size() == N && X.Size() == 3 * N,
               "handback table: " << N << " dynamic DOFs, " << fq.V_strike.Size()
               << " QD DOFs and " << X.Size() / 3 << " positions on rank " << mpi.Rank());
   auto tau = [](const cycle::StationFields &f, int j, bool strike)
   {
      return strike ? -(f.tau_a_strike(j) + f.tau_b_strike(j) + f.eta * f.V_strike(j))
                    : -(f.tau_a_dip(j) + f.tau_b_dip(j) + f.eta * f.V_dip(j));
   };
   // [cycle.opening]: two more columns, the open flags of both regimes; mode
   // contact (B4): two more again, the gap of both regimes
   const bool opening = cfg.opening.On();
   const bool contact = cfg.opening.Contact();
   MFEM_VERIFY(!opening || (fd.open.Size() == N && fq.open.Size() == N),
               "handback table: [cycle.opening] without open flags");
   MFEM_VERIFY(!contact || (fd.gap.Size() == N && fq.gap.Size() == N),
               "handback table: mode contact without the gap");
   const int nc = contact ? 19 : opening ? 17 : 15;
   std::vector<double> rows(static_cast<size_t>(nc) * N);
   for (int j = 0; j < N; j++)
   {
      double *r = &rows[static_cast<size_t>(nc) * j];
      r[0] = X(3 * j);
      r[1] = X(3 * j + 1);
      r[2] = X(3 * j + 2);
      r[3] = fd.slip_strike(j);
      r[4] = fd.slip_dip(j);
      r[5] = fd.V_strike(j);
      r[6] = fq.V_strike(j);
      r[7] = fd.V_dip(j);
      r[8] = fq.V_dip(j);
      r[9] = tau(fd, j, true) / 1e6;
      r[10] = tau(fq, j, true) / 1e6;
      r[11] = tau(fd, j, false) / 1e6;
      r[12] = tau(fq, j, false) / 1e6;
      r[13] = fd.sigma_n(j) / 1e6;
      r[14] = fq.sigma_n(j) / 1e6;
      if (opening)
      {
         r[15] = fd.open(j);
         r[16] = fq.open(j);
      }
      if (contact)
      {
         r[17] = fd.gap(j);
         r[18] = fq.gap(j);
      }
   }
   const int np = mpi.Size(), n_loc = nc * N;
   std::vector<int> counts(np), displs(np);
   MPI_Gather(&n_loc, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, mpi.GetComm());
   std::vector<double> all;
   if (mpi.IsRoot())
   {
      int s = 0;
      for (int r = 0; r < np; r++)
      {
         displs[r] = s;
         s += counts[r];
      }
      all.resize(s);
   }
   MPI_Gatherv(rows.data(), n_loc, MPI_DOUBLE, all.data(), counts.data(), displs.data(),
               MPI_DOUBLE, 0, mpi.GetComm());
   if (!mpi.IsRoot()) { return; }

   const std::string path = cfg.output.dir + "/" + cfg.output.prefix + "_handback_phase" +
                            std::to_string(phase) + ".txt";
   std::ofstream f(path);
   MFEM_VERIFY(f, "handback table: cannot open " << path);
   f << "# fd_cycle hand-back table, dynamic phase " << phase << " -> QD at t = "
     << std::setprecision(17) << t << " s, commit " << commit << "\n# Both regimes from "
     "the same slip, psi and t (plan 6.3): the dynamic regime's fields at its switch, the QD "
     "balance's after the switch back (one QD solve).  tau = -(tau_a + tau_b + eta V) of "
     "each regime (the station files' shear stress), sigma_n' compression positive.  Rows: "
     "the owned fault DOFs of rank 0, then rank 1, ...\n# Columns: x(m), y(m), z(m), "
     "slip_strike(m), slip_dip(m), V_strike_dyn(m/s), V_strike_qd(m/s), V_dip_dyn(m/s), "
     "V_dip_qd(m/s), tau_strike_dyn(MPa), tau_strike_qd(MPa), tau_dip_dyn(MPa), "
     "tau_dip_qd(MPa), sigma_n_dyn(MPa), sigma_n_qd(MPa)"
     << (contact ? ", open_dyn, open_qd ([cycle.opening]: 1 open, contact)"
         : opening ? ", open_dyn, open_qd ([cycle.opening]: 1 open, shear-free)" : "")
     << (contact ? ", gap_dyn(m), gap_qd(m)" : "") << "\n"
     << std::scientific << std::uppercase << std::setprecision(12);
   const size_t n = all.size() / nc;
   size_t jt = 0, js = 0, jv = 0;
   double dt_max = -1.0, ds_max = -1.0, v_max = -1.0;
   for (size_t j = 0; j < n; j++)
   {
      const double *r = &all[nc * j];
      for (int k = 0; k < nc; k++) { f << r[k] << (k + 1 < nc ? " " : "\n"); }
      const double dtau = std::hypot(r[10] - r[9], r[12] - r[11]);
      const double dsig = std::abs(r[14] - r[13]);
      const double vq = std::hypot(r[6], r[8]);
      if (dtau > dt_max) { dt_max = dtau; jt = j; }
      if (dsig > ds_max) { ds_max = dsig; js = j; }
      if (vq > v_max) { v_max = vq; jv = j; }
   }
   f.close();
   auto at = [&](size_t j)
   {
      std::ostringstream s;
      s << std::fixed << std::setprecision(2) << "(" << all[nc * j] / 1e3 << ", "
        << all[nc * j + 1] / 1e3 << ", " << all[nc * j + 2] / 1e3 << ") km";
      return s.str();
   };
   std::ostringstream s;
   s << std::scientific << std::setprecision(6) << "  hand-back table " << path << " (" << n
     << " DOFs): largest |tau_qd - tau_dyn| " << dt_max << " MPa at " << at(jt)
     << "; largest |sigma_n_qd - sigma_n_dyn| " << ds_max << " MPa at " << at(js)
     << "; largest V_qd " << v_max << " m/s at " << at(jv) << ", V_dyn there "
     << std::hypot(all[nc * jv + 5], all[nc * jv + 7]) << " m/s\n";
   std::cout << s.str() << std::flush;
}

/// [cycle.qd] integrator = "esdirk" (option I, phase I5): what the implicit
/// solve did in the QD regime since the last call, printed on root.
void PrintQDImplicit(const MPIContext &mpi, cycle::QDPhase &qd, const std::string &when)
{
#ifdef MFEM_USE_PETSC
   if (!qd.Implicit()) { return; }
   const cycle::QDPhase::ImplicitPhaseInfo f = qd.ImplicitPhase();
   if (!mpi.IsRoot()) { return; }
   const double yr = BP5Params::seconds_per_year;
   const double n = static_cast<double>(std::max(f.steps, 1L));
   std::ostringstream s;
   s << std::setprecision(4) << "  esdirk in the QD regime " << when << ": " << f.steps
     << " steps over " << (f.t1 - f.t0) / yr << " yr; per step " << f.newton_its / n
     << " Newton iterations and " << f.factors / n << " LU factors; " << f.failed_solves
     << " failed stage solves; " << f.columns << " new gap columns; most open points in a "
     << "solve so far " << f.max_open << "; QD loop " << f.wall << " s (" << f.wall / n
     << " s per step, output included)\n";
   std::cout << s.str() << std::flush;
#endif
}

/// [cycle].sigma_n_floor_pa: where the floor acted in the QD regime since the
/// last call, printed on root; the count starts again.  Collective.
void PrintQDFloor(const MPIContext &mpi, const cycle::CycleConfig &cfg, cycle::QDPhase &qd,
                  const std::string &when)
{
   if (cfg.sigma_n_floor <= 0.0) { return; }
   const cycle::QDPhase::FloorInfo f = qd.FloorCount();
   qd.ResetFloorCount();
   if (!mpi.IsRoot()) { return; }
   // its own stream: std::cout's precision stays as the lines after it expect
   std::ostringstream s;
   s << std::setprecision(6) << "  sigma_n' floor in the QD regime " << when << ": ";
   if (f.calls == 0) { s << "never acted\n"; }
   else
   {
      s << "acted in " << f.calls << " RHS evaluations (the most on one rank), up to "
        << f.max_nodes << " DOFs at once on one rank; lowest true sigma_n' "
        << f.min_sigma / 1e6 << " MPa at (" << f.x[0] / 1e3 << ", " << f.x[1] / 1e3
        << ", " << f.x[2] / 1e3 << ") km\n";
   }
   std::cout << s.str();
}

/// [cycle.relaxation]: r and the lowest sigma_n' the law sees before the
/// floor (sigma0 + dsig - r) at the relaxation's time, printed on root.
/// Collective.
void PrintRelaxation(const MPIContext &mpi, const cycle::CycleConfig &cfg,
                     const cycle::QDPhase &qd, const std::string &when)
{
   if (!cfg.Relaxation()) { return; }
   const cycle::QDPhase::RelaxInfo r = qd.RelaxationSummary();
   if (!mpi.IsRoot()) { return; }
   std::ostringstream s;
   s << std::setprecision(6) << "  sigma_n' relaxation " << when << " (t = "
     << qd.Relaxation()->Time() / BP5Params::seconds_per_year << " yr): r in ["
     << r.r_min / 1e6 << ", " << r.r_max / 1e6 << "] MPa; lowest sigma0 + dsig - r "
     << r.min_sigma / 1e6 << " MPa at (" << r.x[0] / 1e3 << ", " << r.x[1] / 1e3 << ", "
     << r.x[2] / 1e3 << ") km\n";
   std::cout << s.str();
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
             << dyn.StableConstant() << "), set up in " << setup_s << " s\n";
   if (cfg.finite_width.On())
   {
      std::cout << "    2.5D finite width: K + M_beta, substrate load at t0 + tau; the CFL "
                << "step times (1 - " << std::scientific << std::setprecision(3)
                << 1.0 - dyn.FiniteWidthDtFactor() << std::defaultfloat << std::setprecision(6)
                << ") for beta / rho; Kelvin-Voigt on K alone\n";
   }
   if (dyn.HasPml())
   {
      std::cout << "    PML band (dynamic): K_p (penalty x kappa), C_d and the memory terms; "
                << "the CFL step divided by sqrt(kappa_max = " << dyn.PmlKappaMax() << ")\n";
   }
   if (dyn.HasSponge())
   {
      std::cout << "    sponge band (comparator): rho sum_a d_a U' on the band elements, "
                << "K and the CFL step unchanged\n";
   }
   std::cout << "    Kelvin-Voigt: ";
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
   if (cfg.quench.On())
   {
      std::cout << "    quench ([cycle.quench]): C_q = d M_rho on every element while max V < V_q = "
                << cfg.quench.V_q << " m/s; d = " << dyn.QuenchRate() << " 1/s";
      if (cfg.quench.AutoRate())
      {
         std::cout << " (\"auto\" = 2 pi c_s,min / L_max: c_s,min " << dyn.QuenchShearSpeedMin()
                   << " m/s, L_max " << dyn.QuenchBoxLength() / 1e3 << " km)";
      }
      std::cout << "; decay time 2 / d = " << 2.0 / dyn.QuenchRate() << " s\n";
   }
}

/// --rank-costs: each rank's work and wait in one dynamic phase (speed-up
/// options, D0 design rev1).  Collective.
///
/// Input: the phase's loop wall time, this rank's time in dyn.Step() and in
/// the fault work (FaultSeconds), and the sums of common/rank_costs.hpp.
/// Output: one screen row per rank (ms per step) and the same numbers in
/// <prefix>_rank_costs_phase<k>.txt.
///
///   wait  = w_sync + w_pml + w_stiff + w_fault + w_out   (blocking MPI)
///   busy  = loop - wait
///   pml   = PML Advance + AddForce without w_pml
///   stiff = the dense stiffness product without w_stiff
///   fault = the fault advance and slip loads without w_fault
///   rhs, damp, mass = wall loads, damping blocks, inverse mass
///   other = busy - pml - stiff - fault - rhs - damp - mass  (vector
///           updates, output, log)
void ReportRankCosts(MPIContext &mpi, cycle::DynPhase &dyn, const std::string &prefix,
                     int phase, const std::string &commit, double loop_s, double step_s,
                     double fault_s)
{
   const double *S = RankCosts::Sum();
   ParMesh &pm = dyn.Op().GetMesh();
   const int NT = 15, NC = NT + 4;
   double mine[NC];
   const double wait = S[RankCosts::WSync] + S[RankCosts::WPml] + S[RankCosts::WStiff] +
                       S[RankCosts::WFault] + S[RankCosts::WOut];
   const double pml = S[RankCosts::Pml] - S[RankCosts::WPml];
   const double stiff = S[RankCosts::Stiff] - S[RankCosts::WStiff];
   const double fault = fault_s - S[RankCosts::WFault];
   mine[0] = loop_s;
   mine[1] = step_s;
   mine[2] = S[RankCosts::WSync];
   mine[3] = S[RankCosts::WPml];
   mine[4] = S[RankCosts::WStiff];
   mine[5] = S[RankCosts::WFault];
   mine[6] = S[RankCosts::WOut];
   mine[7] = loop_s - wait;
   mine[8] = pml;
   mine[9] = stiff;
   mine[10] = fault;
   mine[11] = S[RankCosts::Rhs];
   mine[12] = S[RankCosts::Damp];
   mine[13] = S[RankCosts::Mass];
   mine[14] = mine[7] - pml - stiff - fault - mine[11] - mine[12] - mine[13];
   mine[15] = pm.GetNE();
   mine[16] = dyn.Op().Pml() ? dyn.Op().Pml()->NumBandElements() : 0;
   mine[17] = dyn.Slip().Size() / 3;
   mine[18] = pm.GetNSharedFaces();
   const int np = mpi.Size();
   std::vector<double> all(mpi.IsRoot() ? NC * np : 0);
   MPI_Gather(mine, NC, MPI_DOUBLE, all.data(), NC, MPI_DOUBLE, 0, mpi.GetComm());
   if (!mpi.IsRoot()) { return; }
   const long steps = std::max(1L, static_cast<long>(dyn.Steps()));
   const char *names[NC] = {"loop", "step", "w_sync", "w_pml", "w_stiff", "w_fault", "w_out",
                            "busy", "pml", "stiff", "fault", "rhs", "damp", "mass", "other",
                            "n_elem", "n_band", "n_fault_dofs", "n_shared_faces"
                           };
   const std::string file = prefix + "_rank_costs_phase" + std::to_string(phase) + ".txt";
   std::ofstream f(file);
   f << "# fd_cycle --rank-costs, dynamic phase " << phase << ", " << steps << " steps on "
     << np << " ranks; commit " << commit << "\n# times in s summed over the phase; wait = "
     "w_sync + w_pml + w_stiff + w_fault + w_out; busy = loop - wait; other = busy - pml - "
     "stiff - fault - rhs - damp - mass\n# rank";
   for (int c = 0; c < NC; c++) { f << " " << names[c]; }
   f << "\n" << std::setprecision(9);
   std::cout << "  rank costs, ms per step (wait = blocked in MPI; busy = loop - wait):\n"
             << "  rank    loop    step |  w_sync   w_pml w_stiff w_fault   w_out |    busy"
             "     pml   stiff   fault     rhs    damp    mass   other |   n_elem   n_band"
             "  n_fdof n_shared\n";
   for (int r = 0; r < np; r++)
   {
      const double *a = &all[NC * r];
      f << r;
      for (int c = 0; c < NT; c++) { f << " " << a[c]; }
      for (int c = NT; c < NC; c++) { f << " " << static_cast<long>(a[c]); }
      f << "\n";
      std::cout << std::fixed << std::setprecision(3) << std::setw(6) << r;
      for (int c = 0; c < NT; c++)
      {
         if (c == 2 || c == 7) { std::cout << " |"; }
         std::cout << std::setw(8) << 1e3 * a[c] / steps;
      }
      std::cout << " |" << std::setw(9) << static_cast<long>(a[NT]) << std::setw(9)
                << static_cast<long>(a[NT + 1]) << std::setw(8) << static_cast<long>(a[NT + 2])
                << std::setw(9) << static_cast<long>(a[NT + 3]) << "\n";
   }
   std::cout << "  rank costs file: " << file << "\n" << std::flush;
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
      Vector sk, psi_qd, sn, d, V, psi, gk;
      qd.Slip(sk);
      qd.Psi(psi_qd);
      qd.SigmaN(sn);
      if (cfg.opening.Contact()) { qd.Gap(gk); }   // contact (B4): g = d . n
      cycle::bridge::ToNitsche(dyn.Frames(), sk, qd.SlipRate(), psi_qd, d, V, psi,
                               cfg.opening.Contact() ? &gk : nullptr);
      dyn.Enter(qd.Displacement(), d, V, psi, sn, t);
      if (cfg.opening.On()) { dyn.SetOpenFlags(qd.OpenFlags()); }   // [cycle.opening]
   }
   if (cfg.opening.On())
   {
      const int n_open = qd.Opening()->NumOpenGlobal();
      const real_t gmax = cfg.opening.Contact() ? qd.MaxGapGlobal() : 0.0;   // collective
      if (root)
      {
         if (cfg.opening.Contact())
         {
            std::cout << std::setprecision(6) << "  [cycle.opening] " << n_open
                      << " open DOFs at the switch to the dynamic regime (contact, left "
                      << "out of max V; largest gap " << gmax << " m)\n";
         }
         else
         {
            std::cout << "  [cycle.opening] " << n_open << " open DOFs at the switch to the "
                      << "dynamic regime (shear-free, left out of max V)\n";
         }
      }
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
   PrintQDFloor(mpi, cfg, qd, "before dynamic phase " + std::to_string(phase));
   PrintQDImplicit(mpi, qd, "before dynamic phase " + std::to_string(phase));
   PrintRelaxation(mpi, cfg, qd, "before dynamic phase " + std::to_string(phase));
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
   RankCosts::Reset();   // --rank-costs: this phase's sums
   const double fault_s0 = dyn.Op().FaultSeconds();
   // [cycle.quench] (quench plan rev1, Q4): on at the first step with max V
   // below V_q, off if it rises above V_q again; the counts for the summary
   const cycle::QuenchConfig &qc = cfg.quench;
   long quench_steps = 0, quench_switches = 0;
   double quench_s = 0.0;
   cycle::DynPhase::StepInfo first{}, worst{};
   worst.min_sigma_eff = std::numeric_limits<real_t>::infinity();
   worst.min_sigma_pc = std::numeric_limits<real_t>::infinity();
   real_t Vmax = Vstart;
   // The deck pair's switch state (plan 6.3 hands over slip and psi only).
   Vector snap_sk, snap_psi;
   std::vector<char> snap_open;   // [cycle.opening]: the open flags there
   Vector snap_gap;               // mode contact: the gap there (g = d . n)
   real_t snap_Vmax = 0.0;
   cycle::StationFields snap_fd;   // [cycle.output].handback_table

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
      long floor_hits = 0, floor_dofs = 0;
      if (cfg.sigma_n_floor > 0.0) { dyn.FloorCount(floor_hits, floor_dofs); }
      const long open_end = cfg.opening.On() ? dyn.OpenCount() : 0;   // collective
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
         if (cfg.sigma_n_floor > 0.0)
         {
            std::cout << "\n  sigma_n' floor " << cfg.sigma_n_floor / 1e6 << " MPa: the law saw it "
                      << "at " << floor_hits << " DOF-steps on " << floor_dofs << " DOFs";
         }
         if (cfg.opening.On())
         {
            std::cout << "\n  [cycle.opening] open DOFs at the last step: " << open_end
                      << " (max V above is over the closed DOFs)";
         }
         if (qc.On())
         {
            std::cout << "\n  quench: on for " << quench_steps << " of " << dyn.Steps()
                      << " steps (" << quench_steps * dyn.StepDt() << " s of phase time), "
                      << "switched on " << quench_switches << " time(s), d = "
                      << dyn.QuenchRate() << " 1/s; " << (quench_steps > 0 ? quench_s / quench_steps : 0.0)
                      << " s per quenched step (rank 0)";
         }
         std::cout << "\n  time per dynamic step: " << step_avg << " s (rank 0), " << step_max
                   << " s (slowest rank), on " << mpi.Size() << " ranks; loop with output "
                   << loop_s << " s\n"
                   << "  events: " << out.NumEvents() << (out.InEvent() ? " (still open)" : "")
                   << "\n  files: " << prefix << "_{<stations>,global,catalog,switch}.txt, "
                   << fault_file << "\n"
                   << std::flush;
      }
      if (RankCosts::On())
      {
         ReportRankCosts(mpi, dyn, prefix, phase, commit, loop_s, step_s,
                         dyn.Op().FaultSeconds() - fault_s0);
      }
      phase_over = true;
   };

   while (true)
   {
      const bool quenched_step = dyn.Quenched();
      const double w = MPI_Wtime();
      dyn.Step();
      const double w_step = MPI_Wtime() - w;
      step_s += w_step;
      if (quenched_step)
      {
         quench_s += w_step;
         quench_steps++;
      }
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
      if (qc.On())
      {
         // [cycle.quench]: on below V_q, off above it (a second rupture in
         // the phase runs undamped); the probe pairs below read the same field
         const bool want = (Vmax < qc.V_q);
         if (want != dyn.Quenched())
         {
            if (want)
            {
               dyn.SetQuench();
               quench_switches++;
            }
            else { dyn.ClearQuench(); }
            out.SwitchLog(t, Vmax, "dynamic", "dynamic", want ? "quench_on" : "quench_off");
            if (root)
            {
               std::cout << std::scientific << std::setprecision(3) << "  quench "
                         << (want ? "on" : "off") << " at tau = " << dyn.Tau() << " s (step "
                         << dyn.Steps() << ", max V = " << Vmax << (want ? " < " : " >= ")
                         << "V_q = " << qc.V_q << " m/s), d = " << dyn.QuenchRate() << " 1/s\n"
                         << std::flush;
            }
         }
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
         cycle::bridge::ToQD(dyn.Frames(), dyn.Slip(), dyn.Psi(), snap_sk, snap_psi,
                             cfg.opening.Contact() ? &snap_gap : nullptr);
         if (cfg.opening.On()) { snap_open = dyn.OpenFlags(); }
         if (cfg.output.handback_table && to_qd) { dyn.FillStationFields(snap_fd); }
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
   long open_dyn = 0;
   if (cfg.opening.On())
   {
      open_dyn = std::count(snap_open.begin(), snap_open.end(), char(1));
      MPI_Allreduce(MPI_IN_PLACE, &open_dyn, 1, MPI_LONG, MPI_SUM, mpi.GetComm());
   }
   const real_t Vqd = qd.Reenter(snap_sk, snap_psi, end.t, step0, cfg.dt_init,
                                 cfg.opening.On() ? &snap_open : nullptr,
                                 cfg.opening.Contact() ? &snap_gap : nullptr);
   real_t open_res = 0.0, gmax_qd = 0.0;
   if (cfg.opening.On()) { open_res = qd.OpenShearResidual(); }   // collective
   if (cfg.opening.Contact()) { gmax_qd = qd.MaxGapGlobal(); }    // collective
   if (root)
   {
      std::cout << std::setprecision(6) << "\n  back to QD at t = " << end.t
                << " s (TS step " << step0 << ", first dt " << cfg.dt_init
                << " s): max V dynamic " << snap_Vmax << " m/s, QD balance " << Vqd
                << " m/s" << (cfg.opening.On() ? " (closed DOFs)" : "") << "\n";
      if (cfg.opening.On())
      {
         std::cout << "  [cycle.opening] open DOFs: " << open_dyn << " at the dynamic "
                   << "switch, " << qd.Opening()->NumOpenGlobal()
                   << (cfg.qd.Esdirk() ? " after the consistent start (" : " after the QD settle (")
                   << qd.ReenterFlips() << " flips); max |tau| at open DOFs "
                   << open_res / 1e6 << " MPa\n";
         if (cfg.opening.Contact())
         {
            std::cout << "  [cycle.opening] largest gap after the "
                      << (cfg.qd.Esdirk() ? "consistent start " : "settle ") << gmax_qd
                      << " m\n";
         }
      }
      std::cout << std::flush;
   }
   if (cfg.output.handback_table)
   {
      WriteHandbackTable(mpi, cfg, qd, snap_fd, dyn.FaultPositions(), phase, end.t, commit);
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
   bool restart_drop_relaxation = false;
   const char *output_dir = "";
   int checkpoint_every = -1;
   int max_dyn_steps = -1;
   int end_phase_after = -1;
   int log_every = 1000;
   double wall_limit = -1.0;
   const char *green_cache = "";
   const char *green_hcache = "";
   bool rank_costs = false;
   bool early_sync = false;

   OptionsParser args(argc, argv);
   args.AddOption(&deck, "-d", "--deck", "Cycle deck (TOML): [cycle] and the shared sections.");
   args.AddOption(&max_steps, "-n", "--max-steps",
                  "Stop after this TS step number (counted from the start of the run "
                  "chain, as seas_bp5_full --max-steps).");
   args.AddOption(&restart, "-r", "--restart", "Checkpoint prefix to continue from.");
   args.AddOption(&allow_other_commit, "-aoc", "--allow-other-commit", "-no-aoc",
                  "--no-allow-other-commit",
                  "Accept a checkpoint written by another commit (needs approval).");
   args.AddOption(&restart_drop_relaxation, "-rdr", "--restart-drop-relaxation", "-no-rdr",
                  "--no-restart-drop-relaxation",
                  "Continue a checkpoint written with [cycle.relaxation] on with a deck "
                  "that has it off: the checkpoint's r is not read, so the law's normal "
                  "stress changes by r at the restart (a model change; needs approval).");
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
   args.AddOption(&rank_costs, "-rc", "--rank-costs", "-no-rc", "--no-rank-costs",
                  "Switching on: at the end of every dynamic phase print each rank's "
                  "work and wait per step and write <prefix>_rank_costs_phase<k>.txt "
                  "(measures only; common/rank_costs.hpp).");
   args.AddOption(&early_sync, "-es", "--early-sync", "-no-es", "--no-early-sync",
                  "Switching on, global stepping: exchange the stiffness halo at the start "
                  "of a dynamic step and reduce the fault record at its end, so the ranks "
                  "do not wait for each other in between (same results bit for bit; "
                  "NitscheWaveOperator::SetEarlySync).");
   args.Parse();
   if (!args.Good())
   {
      if (root) { args.PrintUsage(std::cout); }
      return 1;
   }
   RankCosts::On() = rank_costs;

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
   MFEM_VERIFY(!restart_drop_relaxation || !restart_prefix.empty(),
               "--restart-drop-relaxation needs --restart");
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
                << cfg.qd.petsc_options << " (atol " << cfg.qd.ts_atol << ")"
                << (cfg.qd.Esdirk() ? ", integrator esdirk (option I: implicit fault "
                    "solve, arkimex " + cfg.qd.esdirk_type + "; the options file's "
                    "-ts_type is replaced)" : "")
                << (cfg.qd.exact_retry ? ", exact retry ([cycle.qd].exact_retry: a "
                    "rejected step retries from a copy of its start state; the whole "
                    "state is scanned for NaN after each step)" : "") << "\n"
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
                std::to_string(cfg.lts.fault_buffer) +
                (cfg.lts.static_split ? ", static split" : "")
                : std::string("global (central difference)"))
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
                                       : std::string("off"));
      if (sc.volume.enabled)
      {
         std::cout << "; volume dt_qd " << std::to_string(sc.volume.dt_qd) << " s, dt_dyn "
                   << std::to_string(sc.volume.dt_dyn) << " s";
         if (sc.volume_box)
         {
            std::cout << ", box [" << sc.volume_box_min[0] << ", " << sc.volume_box_min[1] << ", "
                      << sc.volume_box_min[2] << "] .. [" << sc.volume_box_max[0] << ", "
                      << sc.volume_box_max[1] << ", " << sc.volume_box_max[2] << "] m";
         }
      }
      std::cout << "\n";
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
      // [cycle.pml] partition_weight > 0: a cut that shares the band's work
      // (cycle/pml_partition.hpp; never with [cycle.lts], which [cycle.pml]
      // refuses).  Otherwise MFEM's default cut.
      const cycle::PmlPartition pcut = cycle::PlanPmlPartition(*smesh, cfg, mpi.GetComm());
      MFEM_VERIFY(lts.part.empty() || pcut.part.empty(), "two mesh cuts ([cycle.lts] and "
                  "[cycle.pml] partition_weight)");
      const int *cut = !lts.part.empty() ? lts.part.data()
                       : (!pcut.part.empty() ? pcut.part.data() : nullptr);
      ParMesh pmesh(mpi.GetComm(), *smesh, cut);
      if (root && !pcut.report.empty()) { std::cout << "  " << pcut.report << "\n"; }
      cycle::DynPhase::MapLts(lts, *smesh, pmesh);
      if (root && lts.on) { std::cout << "  LTS: " << lts.report << "\n"; }
      smesh.reset();

      // ---- QD regime ----
      cycle::QDPhase qd(pmesh, mpi, cfg);
      if (root && cfg.finite_width.On())
      {
         std::cout << "  2.5D finite width ([cycle.finite_width]): "
                   << cycle::FiniteWidthReport(cfg) << "\n";
      }
      const spatial::SpatialFrictionConfig &s = cfg.shared;
      const MaterialField material = MaterialField::MakeConstant(
         s.material_fallback.lambda, s.material_fallback.mu, s.material_fallback.rho);
      // sigma_n_default is already effective: no pore pressure here (Load
      // refuses a non-zero [pore_pressure]).
      const spatial::PorePressureSpec no_pp;
      cycle::FaultParams par = cycle::BuildFaultParams(
         qd.Table(), *s.rate_state, cfg.fields, material, pmesh, no_pp, cfg.fault_attr);
      if (root && !par.a_report.empty())
      {
         std::cout << "  [cycle.fields.a] tip zone: " << par.a_report;
      }
      // [cycle.prestress].a = "amax" (BP3): a_max over every rank's DOFs.
      cycle::PrestressConfig pcfg = cfg.prestress;
      if (pcfg.AMax())
      {
         pcfg.a_max = cycle::GlobalMaxA(par, mpi.GetComm());
         if (root) { std::cout << "  prestress: a = \"amax\", a_max = " << pcfg.a_max << "\n"; }
      }
      cycle::Prestress pre;
      if (cfg.TensorPrestress())
      {
         // [stress] kind = "constant_tensor": the tensor's traction on each
         // fault DOF, the dynamic operator's projection (shared function).
         Vector K;
         qd.ProjectTensor(cfg.stress_tensor, K);
         pre = cycle::BuildTensorPrestress(par, K, pcfg, qd.Table());
         real_t r[6] = {std::numeric_limits<real_t>::max(), 0.0,
                        std::numeric_limits<real_t>::max(), 0.0,
                        std::numeric_limits<real_t>::max(), 0.0};
         for (int i = 0; i < pre.sigma0.Size(); i++)
         {
            const real_t sn = pre.sigma0(i);
            const real_t tau = std::sqrt(pre.tau0(2 * i) * pre.tau0(2 * i) +
                                         pre.tau0(2 * i + 1) * pre.tau0(2 * i + 1));
            const real_t v[3] = {sn, tau, tau / sn};
            for (int k = 0; k < 3; k++)
            {
               r[2 * k] = std::min(r[2 * k], v[k]);
               r[2 * k + 1] = std::max(r[2 * k + 1], v[k]);
            }
         }
         for (int k = 0; k < 3; k++)
         {
            r[2 * k] = -mpi.GlobalMax(-r[2 * k]);
            r[2 * k + 1] = mpi.GlobalMax(r[2 * k + 1]);
         }
         if (root)
         {
            const std::array<real_t, 9> &S = cfg.stress_tensor;
            std::cout << "  prestress: [stress] constant_tensor sigma_xx yy zz xy yz xz = "
                      << S[0] << " " << S[4] << " " << S[8] << " " << S[1] << " " << S[5]
                      << " " << S[2] << " Pa (tension positive)\n"
                      << "    sigma_n [" << r[0] / 1e6 << ", " << r[1] / 1e6 << "] MPa, |tau0| ["
                      << r[2] / 1e6 << ", " << r[3] / 1e6 << "] MPa, |tau0| / sigma_n ["
                      << r[4] << ", " << r[5] << "]; V_init along -tau0 with the deck's |V_init|\n";
         }
      }
      else { pre = cycle::BuildPrestress(par, pcfg); }
      if (cfg.sigma_n_floor > 0.0)
      {
         // [cycle].sigma_n_floor_pa must not act at t = 0, where psi0 comes
         // from the true sigma_n'
         real_t smin = std::numeric_limits<real_t>::max();
         for (int i = 0; i < pre.sigma0.Size(); i++) { smin = std::min(smin, pre.sigma0(i)); }
         smin = -mpi.GlobalMax(-smin);
         MFEM_VERIFY(smin >= cfg.sigma_n_floor, "[cycle].sigma_n_floor_pa = "
                     << cfg.sigma_n_floor << " Pa is above the prestress sigma_n' (lowest "
                     << smin << " Pa): the floor would act at t = 0");
         if (root)
         {
            std::cout << "  sigma_n' floor: " << cfg.sigma_n_floor / 1e6 << " MPa in the "
                      << "friction law of both regimes ([cycle].sigma_n_floor_pa); the "
                      << "prestress sigma_n' is >= " << smin / 1e6 << " MPa\n";
         }
      }
      if (cfg.Relaxation() && root)
      {
         std::cout << "  sigma_n' relaxation: t_r = " << cfg.relax_t / BP5Params::seconds_per_year
                   << " yr ([cycle.relaxation]); the law sees sigma0 + dsig - r, dr/dt = "
                   << "(dsig - r) / t_r in the QD regime, r frozen in the dynamic regime\n";
      }
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
         // --restart-drop-relaxation: only checkpoint on, deck off.  The
         // checkpoint's r (what the law's normal stress was reduced by) is
         // not carried: from the restart the law sees sigma0 + dsig.
         const bool drop_relax = restart_drop_relaxation && info.relax && !cfg.Relaxation();
         MFEM_VERIFY(!restart_drop_relaxation || drop_relax, "--restart-drop-relaxation needs "
                     "a checkpoint written with [cycle.relaxation] on and a deck with it off "
                     "(checkpoint " << (info.relax ? "on" : "off") << ", deck "
                     << (cfg.Relaxation() ? "on" : "off") << ")");
         MFEM_VERIFY(info.relax == cfg.Relaxation() || drop_relax,
                     "restart: the checkpoint was written "
                     "with [cycle.relaxation] " << (info.relax ? "on" : "off") << ", the deck "
                     "has it " << (cfg.Relaxation() ? "on" : "off"));
         if (drop_relax)
         {
            const real_t inf = std::numeric_limits<real_t>::infinity();
            real_t rlo = info.relax_r.Size() ? info.relax_r.Min() : inf;
            real_t rhi = info.relax_r.Size() ? info.relax_r.Max() : -inf;
            rlo = -mpi.GlobalMax(-rlo);
            rhi = mpi.GlobalMax(rhi);
            if (root)
            {
               std::cout << "  restart: --restart-drop-relaxation: the checkpoint's sigma_n' "
                         << "relaxation r (" << rlo / 1e6 << " to " << rhi / 1e6
                         << " MPa over the fault) is dropped; the law's normal stress "
                         << "changes by r at this restart\n";
            }
         }
         else if (info.relax)
         {
            MFEM_VERIFY(info.relax_t_relax == cfg.relax_t && info.relax_t_n == info.t,
                        "restart: the checkpoint's relaxation has t_relax "
                        << info.relax_t_relax << " s at t " << info.relax_t_n
                        << " s; the deck has " << cfg.relax_t << " s, the checkpoint t "
                        << info.t << " s");
            qd.SetRelaxationState(info.relax_t_n, info.relax_r, info.relax_d);
         }
         MFEM_VERIFY(info.opening == cfg.opening.On(), "restart: the checkpoint was written "
                     "with [cycle.opening] " << (info.opening ? "on" : "off") << ", the deck "
                     "has it " << (cfg.opening.On() ? "on" : "off"));
         if (info.opening)
         {
            MFEM_VERIFY(static_cast<int>(info.open.size()) == qd.Table().Size(),
                        "restart: " << info.open.size() << " open flags for "
                        << qd.Table().Size() << " owned fault DOFs");
            qd.SetOpenFlags(info.open);
         }
         MFEM_VERIFY(info.contact == cfg.opening.Contact(), "restart: the checkpoint "
                     "was written with [cycle.opening] mode "
                     << (info.contact ? "\"contact\"" : "without the gap block")
                     << ", the deck has mode "
                     << (cfg.opening.Contact() ? "\"contact\"" : "\""
                         + cfg.opening.mode + "\""));
         if (info.contact) { qd.SetGap(info.gap); }   // after the flags (B4)
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
      // vertex-normal plan rev1, P3: the fault's vertices with their turn
      // angles and normals, for checking the averaged vertex normals
      if (cfg.dyn.averaged_vertex_normals)
      {
         Vector vrows;
         cycle::FaultVertexRows(*qd.Domain().FrameField(), mpi.GetComm(), vrows);
         std::ostringstream vd;
         vd << "fd_cycle fault vertices: [fault_geometry] vertex_normals = \"averaged\", "
               "vertex_normal_max_turn_deg = " << cfg.dyn.vertex_normal_max_turn_deg
            << "; kind 0 tip, 1 smooth (n_v averaged), 2 kink (each face keeps its normal, "
               "n_v = 0); commit " << commit;
         cycle::WriteFaultTable(cfg.output.dir + "/" + cfg.output.prefix + "_vertices.txt",
                                mpi.GetComm(), vd.str(), cycle::FaultVertexColumns(), vrows);
      }
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
         if (cfg.snapshots.volume.enabled)
         {
            // Two 3-vectors in float64 per point.
            std::cout << "  volume snapshots: " << snaps.VolumePath() << " ("
                      << snaps.VolumePoints() << " points, " << snaps.VolumeCells()
                      << " cells, " << 48.0 * static_cast<double>(snaps.VolumePoints()) / 1e6
                      << " MB per snapshot, " << snaps.VolumeSteps()
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
         if (const cycle::NormalStressRelaxation *rx = qd.Relaxation())
         {
            info.relax = true;
            info.relax_t_relax = rx->TRelax();
            info.relax_t_n = rx->Time();
            info.relax_r = rx->R();
            info.relax_d = rx->Delta();
         }
         if (cfg.opening.On())
         {
            info.opening = true;
            info.open = qd.OpenFlags();
            if (cfg.opening.Contact())   // V7 (B4): the gap rides along
            {
               info.contact = true;
               qd.Gap(info.gap);
            }
         }
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
                                                 qd.Domain().FrameField(),
                                                 qd.SubstrateShape(),
                                                 qd.PmlBandPtr(), qd.PmlWallSpeedFactor());
         if (cfg.Relaxation()) { dyn->SetNormalStressRelaxation(qd.Relaxation()->R()); }
         if (early_sync)
         {
            MFEM_VERIFY(!lts.on, "--early-sync with [cycle.lts]: the LTS run has its own "
                        "halo and record");
            dyn->Op().SetEarlySync(true);
            if (root)
            {
               std::cout << "  --early-sync: stiffness halo at the start of a dynamic step, "
                         "fault record reduced at its end\n";
            }
         }
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
      // [cycle.opening]: the open set of the start state (a restart has the
      // checkpoint's)
      if (cfg.opening.On() && restart_prefix.empty())
      {
         const long flips = qd.SettleOpening(t);
         if (root)
         {
            std::cout << "  [cycle.opening] mode \"" << cfg.opening.mode << "\": "
                      << qd.Opening()->NumOpenGlobal() << " open DOFs at t = " << t
                      << " s (" << flips << " flips; " << qd.Opening()->UnitProducts()
                      << " Green's-function columns)\n";
         }
      }
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
         // cycle/qd_timers.hpp: the pieces of the monitor, for -log_view
         cycle::QDTimer tm_scan("QDMonScan");
         const real_t Vmax = qd.MaxV();
         bool bad = !std::isfinite(Vmax);
         const Vector &st = qd.State();
         // the first 100 entries per rank, or all with [cycle.qd].exact_retry
         // (Vista 1127480: the bad entry lay further down its rank's list)
         const int nscan = cfg.qd.exact_retry ? st.Size() : std::min(st.Size(), 100);
         for (int i = 0; i < nscan && !bad; i++)
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
         tm_scan.Stop();
         cycle::QDTimer tm_rows("QDMonRows");
         qd_vel_ok = false;
         const cycle::EventUpdate ev = out->UpdateEvents(tt, Vmax);
         if (bench) { bench->OnEvents(ev, tt, fill, speed_qd); }
         const bool row = out->Write(tt, Vmax, fill);
         out->WriteGlobal(tt, Vmax);
         if (row && bench) { bench->Rows(tt, fill, vol_qd); }
         if (bench) { bench->Profiles(tt, Vmax, fill); }
         tm_rows.Stop();
         cycle::QDTimer tm_snaps("QDMonSnaps");
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
         tm_snaps.Stop();
         cycle::QDTimer tm_ckpt("QDMonCkpt");
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
         PrintQDFloor(mpi, cfg, qd, "since the last switch");
         PrintQDImplicit(mpi, qd, "since the last switch");
         PrintRelaxation(mpi, cfg, qd, "at the end");
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
                      << (cfg.qd.exact_retry
                          ? "  exact retries: " + std::to_string(qd.ExactRetries()) + " ("
                            + std::to_string(qd.ExactRetriesNonFinite())
                            + " after a rollback with non-finite entries)\n"
                          : std::string())
                      << (qd.Implicit()
                          ? "  esdirk (option I): Newton iterations "
                            + std::to_string(qd.Implicit()->GetCounts().newton_its)
                            + ", LU factors "
                            + std::to_string(qd.Implicit()->GetCounts().factors)
                            + ", failed stage solves "
                            + std::to_string(qd.Implicit()->GetCounts().newton_fails)
                            + ", gap columns "
                            + std::to_string(qd.Implicit()->GetCounts().columns)
                            + ", most open points in a solve "
                            + std::to_string(qd.Implicit()->GetCounts().max_open) + "\n"
                          : std::string())
                      << "  dynamic phases: " << ctl.Phases() << "\n"
                      << "  events: " << out->NumEvents() << "\n"
                      << "  snapshots: fault " << snaps.FaultSteps() << ", surface "
                      << snaps.SurfaceSteps()
                      << (cfg.snapshots.volume.enabled
                          ? ", volume " + std::to_string(snaps.VolumeSteps()) : std::string())
                      << "\n"
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
