// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_cycle_p3.cpp — fully dynamic BP5 cycle, P3 tests (plan rev 3, phase
// P3): the switch back from the dynamic regime to QD (plan 6.3).
//
//   bridge_roundtrip   QD -> DYN -> QD with no time advance.  A QD state with
//                      slip and a perturbed psi at t_s = 1e7 s, one QD solve,
//                      Bridge::ToNitsche and DynPhase::Enter, then
//                      Bridge::ToQD from W's slip and psi and
//                      QDPhase::Reenter at W's time: slip, psi and the whole
//                      QD state come back bit for bit, the time is t_s
//                      exactly, and the re-entry solve repeats the first
//                      (max V, slip rate, traction, displacement) bit for bit.
//                      The fault is the plane y = 0, so the frame map is a
//                      signed permutation and the round trip is exact; on a
//                      curved fault it would hold to round-off.
//   qd_restart_after_switch
//                      A: QD from the IC set for K steps, then a dynamic
//                      phase of n steps from that state (DynPhase::Enter at
//                      t_K), Bridge::ToQD and Reenter (TS step counter K),
//                      a checkpoint, then M QD steps.
//                      B: a new QDPhase (its own setup and factorization)
//                      Reenter'ed with the same slip, psi, t and step: a
//                      fresh QD run from that state.  M steps.
//                      C: a new QDPhase restarted from A's checkpoint the way
//                      the driver restarts (InitTS, state, RestoreTS).
//                      B and C must repeat A bit for bit: state, time, next
//                      dt, step number and rejections.
//   catalog_matches_controller
//                      Synthetic max V histories (no solver) through a
//                      RegimeController and FdCycleOutput's event catalog
//                      with the same threshold and dwell: the phase must end
//                      on the step the catalog closes the event, and the
//                      row's t_end must be the first time below.  Cases: the
//                      Mac run's dt with a dip that comes back above and a
//                      step exactly at the threshold, a dt that divides the
//                      dwell (10,000 steps), T_dwell = 0, and t0 = 100 yr.
//
// Bit for bit holds on 1 rank.  On more ranks MUMPS's parallel solve is not
// reproducible to the last bit from one call to the next (the same solve
// twice differs by about 1e-16 in u, 1e-13 in the traction; bridge_roundtrip
// prints this control), so there the re-entry solve is compared to 1e-12
// and B and C to 1e-10 relative (the next dt to 1e-8, see DtTol); slip,
// psi and the state copies stay exact.
// `make test-cycle-p3` runs 1 rank (exact) and 2 ranks.
//
// Usage (from miniapps/seas; PETSc needed):
//   mpirun -np N ./seas_test_cycle_p3 [--deck config/fd_cycle/bp5_fd_first_event.toml]
//          [--qd-steps-before 5] [--dyn-steps 20] [--qd-steps-after 10]
//          [--work-dir DIR]
// Exit code 0 = PASSED.

#include "mfem.hpp"
#include "../../cycle/scaled_mesh.hpp"
#include "../../cycle/cycle_config.hpp"
#include "../../cycle/fault_fields.hpp"
#include "../../cycle/prestress.hpp"
#include "../../cycle/checkpoint.hpp"
#include "../../cycle/qd_phase.hpp"
#include "../../cycle/dyn_phase.hpp"
#include "../../cycle/state_bridge.hpp"
#include "../../cycle/regime_controller.hpp"
#include "../../cycle/fd_cycle_output.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

using namespace mfem;
using namespace mfem::seas;

namespace
{

struct Checker
{
   const MPIContext &mpi;
   int failures = 0;
   /// Tolerances of what goes through a parallel MUMPS solve (file comment).
   real_t SolveTol() const { return mpi.Size() == 1 ? 0.0 : 1e-12; }
   real_t RunTol() const { return mpi.Size() == 1 ? 0.0 : 1e-10; }
   /// The RK45 controller's next dt: its error estimate is a difference of
   /// nearly equal stages, so the solve's round-off grows in it (2 Mac
   /// ranks: 2e-10 and 4e-10 with the state at 1e-13).
   real_t DtTol() const { return mpi.Size() == 1 ? 0.0 : 1e-8; }

   /// Pass when m <= tol (m and tol global); prints both.
   void Report(const std::string &name, real_t m, real_t tol, const std::string &what)
   {
      const bool ok = (m <= tol);
      if (!ok) { failures++; }
      if (mpi.IsRoot())
      {
         std::cout << "  " << (ok ? "ok  " : "FAIL") << " " << std::left << std::setw(44)
                   << name << std::right << " " << what << " " << std::scientific
                   << std::setprecision(3) << m << " (tolerance " << tol << ")\n";
      }
   }
   void True(const std::string &name, bool ok_local, const std::string &detail)
   {
      const bool ok = mpi.GlobalMinInt(ok_local ? 1 : 0) == 1;
      if (!ok) { failures++; }
      if (mpi.IsRoot())
      {
         std::cout << "  " << (ok ? "ok  " : "FAIL") << " " << std::left << std::setw(44)
                   << name << std::right << " " << detail << "\n";
      }
   }
   /// max |x - y| / max |y| over all ranks (0 when both are zero).
   real_t RelInf(const Vector &x, const Vector &y) const
   {
      MFEM_VERIFY(x.Size() == y.Size(), "RelInf: sizes " << x.Size() << " vs " << y.Size());
      real_t d = 0.0, s = 0.0;
      for (int i = 0; i < x.Size(); i++)
      {
         d = std::max(d, std::abs(x(i) - y(i)));
         s = std::max(s, std::abs(y(i)));
      }
      d = mpi.GlobalMax(d);
      s = mpi.GlobalMax(s);
      return (s > 0.0) ? d / s : d;
   }
};

/// A smooth deterministic field of position, partition independent.
real_t Field(const real_t *x, int c)
{
   return std::sin(1.3e-4 * x[0] + 0.7 * c) * std::cos(2.1e-4 * x[2] - 0.3 * c)
          + 0.25 * std::sin(3.7e-5 * x[1] + 1.1 * c);
}

struct Setup
{
   const cycle::CycleConfig &cfg;
   ParMesh &pmesh;
   MPIContext &mpi;
   const cycle::FaultParams &par;
   const cycle::Prestress &pre;

   std::unique_ptr<cycle::QDPhase> NewQD() const
   {
      auto qd = std::make_unique<cycle::QDPhase>(pmesh, mpi, cfg);
      qd->Setup(par, pre);
      return qd;
   }
};

// ------------------------------------------------------------- bridge_roundtrip

void TestBridgeRoundtrip(Checker &c, const Setup &S, cycle::QDPhase &qd, cycle::DynPhase &dyn)
{
   const MPIContext &mpi = c.mpi;
   const real_t ts = 1.0e7;   // the creep step at the fault edge stays 1 cm (test_cycle_p2)
   if (mpi.IsRoot()) { std::cout << "\nbridge_roundtrip (t_s = " << ts << " s)\n"; }

   // Slip up to 2 m, zero at the fault's edges; psi0 plus up to 0.01.
   const Vector &X = dyn.FaultPositions();
   const int N = dyn.NumOwnedFaultDOFs();
   Vector s0(2 * N), psi0(N);
   for (int i = 0; i < N; i++)
   {
      const real_t x = X(3 * i) / 50.0e3, z = X(3 * i + 2) / 40.0e3;
      const real_t w = std::max(0.0, 1.0 - x * x) * std::max(0.0, 1.0 - z * z);
      s0(2 * i) = 2.0 * w * Field(&X(3 * i), 0);
      s0(2 * i + 1) = 2.0 * w * Field(&X(3 * i), 1);
      psi0(i) = S.pre.psi0(i) + 0.01 * Field(&X(3 * i), 2);
   }
   qd.SetFaultState(s0, psi0);
   const Vector st0 = qd.State();
   const real_t V0 = qd.Solve(ts);
   const Vector V0k = qd.SlipRate(), tr0 = qd.Traction(), sn0 = qd.NormalTraction(),
                u0 = qd.Displacement();
   // Control: the same solve again.  The QD right-hand side is a function
   // of (t, state) alone, so any difference here is the linear solver's own.
   qd.Solve(ts);
   const real_t ctl_u = c.RelInf(qd.Displacement(), u0), ctl_tr = c.RelInf(qd.Traction(), tr0);
   if (mpi.IsRoot())
   {
      std::cout << std::scientific << std::setprecision(3)
                << "  control: the same QD solve twice: displacement " << ctl_u
                << ", traction " << ctl_tr << " (max rel diff)\n";
   }

   // QD -> DYN.
   Vector sk, psi_qd, sn, d, V, psi;
   qd.Slip(sk);
   qd.Psi(psi_qd);
   qd.SigmaN(sn);
   cycle::bridge::ToNitsche(dyn.Frames(), sk, qd.SlipRate(), psi_qd, d, V, psi);
   dyn.Enter(qd.Displacement(), d, V, psi, sn, ts);
   c.True("W's time after Enter is t_s", dyn.Time() == ts && dyn.Tau() == 0.0,
          "t0 + tau with tau = 0");

   // DYN -> QD with no step.
   Vector s1, psi1;
   cycle::bridge::ToQD(dyn.Frames(), dyn.Slip(), dyn.Psi(), s1, psi1);
   c.Report("slip: QD -> W -> QD", c.RelInf(s1, sk), 0.0, "max rel diff");
   c.Report("psi: QD -> W -> QD", c.RelInf(psi1, psi_qd), 0.0, "max rel diff");

   const real_t V1 = qd.Reenter(s1, psi1, dyn.Time(), 0, S.cfg.dt_init);
   c.Report("QD state after Reenter vs before", c.RelInf(qd.State(), st0), 0.0,
            "max rel diff");
   const real_t ts_ = c.SolveTol();
   c.Report("re-entry solve: max V", std::abs(V1 - V0) / V0, ts_, "rel diff");
   c.Report("re-entry solve: slip rate", c.RelInf(qd.SlipRate(), V0k), ts_, "max rel diff");
   c.Report("re-entry solve: traction", c.RelInf(qd.Traction(), tr0), ts_, "max rel diff");
   c.Report("re-entry solve: normal traction", c.RelInf(qd.NormalTraction(), sn0), ts_,
            "max rel diff");
   c.Report("re-entry solve: displacement", c.RelInf(qd.Displacement(), u0), ts_,
            "max rel diff");
   c.True("TS after Reenter: step counter and dt",
          qd.StepNumber() == 0 && qd.NextDt() == S.cfg.dt_init, "step 0, dt = dt_init");
   if (mpi.IsRoot())
   {
      std::cout << std::scientific << std::setprecision(6) << "  max V " << V0
                << " m/s; max |slip| 2 m\n";
   }
}

// ------------------------------------------------- qd_restart_after_switch

struct RunResult
{
   Vector state;
   real_t t = 0.0, dt = 0.0;
   int step = 0, rejections = 0;
};

/// M accepted steps from the current TS (t, dt0) on; no output.
RunResult RunSteps(cycle::QDPhase &qd, real_t t, real_t dt0, int last_step)
{
   RunResult r;
   real_t dt = dt0;
   qd.Run(t, dt, 1.0e30, last_step, [](int, real_t) { return true; });
   r.state = qd.State();
   r.t = t;
   r.dt = qd.NextDt();
   r.step = qd.StepNumber();
   r.rejections = qd.Rejections();
   return r;
}

void Compare(Checker &c, const std::string &who, const RunResult &x, const RunResult &a)
{
   const real_t tol = c.RunTol();
   c.Report(who + ": state vs A", c.RelInf(x.state, a.state), tol, "max rel diff");
   c.Report(who + ": t vs A", std::abs(x.t - a.t) / std::abs(a.t), tol, "rel diff");
   c.Report(who + ": next dt vs A", std::abs(x.dt - a.dt) / std::abs(a.dt), c.DtTol(),
            "rel diff");
   c.True(who + ": step and rejections vs A", x.step == a.step && x.rejections == a.rejections,
          "step " + std::to_string(x.step) + ", rejections " + std::to_string(x.rejections));
}

void TestQDRestartAfterSwitch(Checker &c, const Setup &S, const cycle::CycleConfig &cfg0,
                              int K, int n_dyn, int M, const std::string &work)
{
   const MPIContext &mpi = c.mpi;
   const bool root = mpi.IsRoot();
   if (root)
   {
      std::cout << "\nqd_restart_after_switch (K = " << K << " QD steps, " << n_dyn
                << " dynamic steps, M = " << M << " QD steps)\n";
   }
   Vector snap_s, snap_psi;
   real_t t_sw = 0.0;
   RunResult A;
   const std::string prefix = work + "/p3_switch";
   {
      // A: QD for K steps from the IC set at t = 0.
      std::unique_ptr<cycle::QDPhase> qd = S.NewQD();
      qd->InitTS(S.cfg.dt_init);
      real_t t = 0.0, dt = S.cfg.dt_init;
      qd->Run(t, dt, 1.0e30, K, [](int, real_t) { return true; });
      MFEM_VERIFY(qd->StepNumber() == K, "A: " << qd->StepNumber() << " QD steps, not " << K);

      // A dynamic phase from that state.
      const real_t VK = qd->Solve(t);
      cycle::DynPhase dyn(S.pmesh, cfg0, qd->Table(), S.par, S.pre);
      Vector sk, psi_qd, sn, d, V, psi;
      qd->Slip(sk);
      qd->Psi(psi_qd);
      qd->SigmaN(sn);
      cycle::bridge::ToNitsche(dyn.Frames(), sk, qd->SlipRate(), psi_qd, d, V, psi);
      dyn.Enter(qd->Displacement(), d, V, psi, sn, t);
      for (int k = 0; k < n_dyn; k++) { dyn.Step(); }
      t_sw = dyn.Time();
      cycle::bridge::ToQD(dyn.Frames(), dyn.Slip(), dyn.Psi(), snap_s, snap_psi);
      const real_t Vdyn = dyn.MaxV();

      // Back to QD (the driver's order): Reenter, checkpoint, M steps.
      const real_t Vq = qd->Reenter(snap_s, snap_psi, t_sw, K, S.cfg.dt_init);
      cycle::CheckpointInfo info;
      info.t = t_sw;
      info.dt_next = qd->NextDt();
      info.step = K;
      info.rejections = qd->Rejections();
      info.commit = SEAS_CYCLE_COMMIT;
      info.ic = S.cfg.ic;
      info.dyn_phases = 1;
      cycle::CheckpointData cd;
      cd.state = qd->State();
      cd.displacement = qd->Operator().GetDisplacement();
      cd.traction = qd->Operator().GetTraction();
      cd.slip_rate = qd->Fault().GetSlipRate();
      cycle::WriteCycleCheckpoint(prefix, info, cd, mpi);
      A = RunSteps(*qd, t_sw, S.cfg.dt_init, K + M);
      if (root)
      {
         std::cout << std::scientific << std::setprecision(6) << "  A: max V " << VK
                   << " m/s at t_K = " << t_sw - n_dyn * dyn.Dt() << " s; dynamic to t = "
                   << t_sw << " s (max V " << Vdyn << "); QD balance " << Vq << " m/s; "
                   << M << " QD steps to t = " << A.t << " s\n";
      }
      c.True("A: M QD steps after the switch", A.step == K + M,
             "step " + std::to_string(A.step));
   }
   {
      // B: a fresh QD run from the switch state.
      std::unique_ptr<cycle::QDPhase> qd = S.NewQD();
      qd->InitTS(S.cfg.dt_init);
      qd->Reenter(snap_s, snap_psi, t_sw, K, S.cfg.dt_init);
      const RunResult B = RunSteps(*qd, t_sw, S.cfg.dt_init, K + M);
      Compare(c, "B (fresh QD run)", B, A);
   }
   {
      // C: the driver's restart from A's switch checkpoint.
      std::unique_ptr<cycle::QDPhase> qd = S.NewQD();
      qd->InitTS(S.cfg.dt_init);
      cycle::CheckpointInfo info;
      cycle::CheckpointData cd;
      cycle::ReadCycleCheckpoint(prefix, SEAS_CYCLE_COMMIT, false, info, cd, mpi);
      c.True("C: checkpoint block", info.step == K && info.t == t_sw && info.dyn_phases == 1,
             "step, t and dynamic phases read back");
      qd->State() = cd.state;
      qd->Operator().SetDisplacement(cd.displacement);
      qd->Fault().SetSlipRate(cd.slip_rate);
      qd->RestoreTS(info.step);
      const RunResult C = RunSteps(*qd, info.t, info.dt_next, K + M);
      Compare(c, "C (restart from the switch checkpoint)", C, A);
   }
}

} // namespace

// ------------------------------------------------- catalog_matches_controller

/// One piece of a synthetic history: `steps` steps at max V = `V`.
struct Piece
{
   int steps;
   real_t V;
};

/// The last data row of a root-written row file (empty when none).
std::vector<real_t> LastRow(const std::string &path, int &rows)
{
   std::ifstream in(path);
   std::string line;
   std::vector<real_t> last;
   rows = 0;
   while (std::getline(in, line))
   {
      if (line.empty() || line[0] == '#') { continue; }
      std::istringstream ls(line);
      last.clear();
      real_t v;
      while (ls >> v) { last.push_back(v); }
      rows++;
   }
   return last;
}

void CatalogCase(Checker &c, const std::string &work, const std::string &name, real_t t0,
                 real_t dt, real_t dwell, const std::vector<Piece> &history)
{
   const real_t V_off = 1.0e-3;
   cycle::OutputConfig oc;
   oc.dir = work;
   oc.prefix = "catalog_" + name;
   oc.event_V = V_off;
   oc.event_dwell = dwell;
   cycle::SwitchConfig sw;
   sw.V_off = V_off;
   sw.T_dwell = dwell;
   sw.T_max = 1e30;
   cycle::RegimeController ctl(sw);
   int end_step = -1, close_step = -1, n = 0;
   real_t t = t0, tau = 0.0, dip_start = -1.0, prev_V = 1.0, t_end_phase = -1.0;
   {
      cycle::FdCycleOutput out(oc, c.mpi.GetComm(), {}, "test");
      // As the driver: the QD balance at t0 opens the event, then the phase.
      out.UpdateEvents(t, 1.5e-2);
      ctl.EnterDynamic(t);
      for (const Piece &p : history)
      {
         for (int k = 0; k < p.steps && (end_step < 0 || close_step < 0); k++)
         {
            n++;
            tau += dt;
            t = t0 + tau;
            if (p.V < V_off && prev_V >= V_off) { dip_start = t; }
            prev_V = p.V;
            out.UpdateEvents(t, p.V);
            const cycle::DynVerdict v = ctl.AfterDynStep(p.V, t, tau);
            if (end_step < 0 && v == cycle::DynVerdict::EndOfPhase) { end_step = n; t_end_phase = t; }
            if (close_step < 0 && !out.InEvent()) { close_step = n; }
         }
      }
      out.Flush();
   }
   std::ostringstream d;
   d << "phase ends at step " << end_step << " (t = " << std::setprecision(12) << t_end_phase
     << " s), catalog closes at step " << close_step;
   c.True("catalog_matches_controller " + name, end_step > 0 && close_step == end_step, d.str());
   int rows = 0;
   std::vector<real_t> row;
   bool ok = true;
   real_t rel = 0.0;
   if (c.mpi.IsRoot())
   {
      row = LastRow(work + "/" + oc.prefix + "_catalog.txt", rows);
      ok = (rows == 1 && row.size() == 4);
      if (ok) { rel = std::abs(row[2] - dip_start) / std::abs(dip_start); }
      ok = ok && rel <= 1e-11 && row[0] == 1.0 && row[1] == t0;
   }
   std::ostringstream e;
   e << "rows " << rows << ", t_end vs first time below: relative " << std::scientific
     << std::setprecision(2) << rel;
   c.True("catalog_matches_controller " + name + " row", ok, e.str());
}

void TestCatalogMatchesController(Checker &c, const std::string &work)
{
   if (c.mpi.IsRoot()) { std::cout << "\ncatalog_matches_controller\n"; }
   const real_t dt_mac = 1.719631e-3;
   // 3.4 s of rupture, 5.2 s below, back above for 0.9 s (one step exactly at
   // the threshold), then 14 s below: only the last dip ends the phase.
   CatalogCase(c, work, "mac_dt", 0.0, dt_mac, 10.0,
               {{2000, 3.0}, {3000, 5e-4}, {250, 2e-3}, {1, 1e-3}, {250, 2e-3},
                {8000, 8.6e-4}});
   // 10 s is exactly 10,000 steps of 1 ms (the tie the old code missed by
   // one step either way).
   CatalogCase(c, work, "round_dt", 0.0, 1e-3, 10.0, {{100, 1.0}, {12000, 9e-4}});
   CatalogCase(c, work, "dwell_0", 0.0, 1e-3, 0.0, {{50, 1.0}, {10, 9e-4}});
   // A later event: t0 = 100 yr, where t carries fewer digits of tau.
   CatalogCase(c, work, "t0_100yr", 100.0 * 365.25 * 86400.0, dt_mac, 10.0,
               {{2000, 3.0}, {8000, 8.6e-4}});
}

int main(int argc, char *argv[])
{
   MPIContext mpi(&argc, &argv);
   const char *deck = "config/fd_cycle/bp5_fd_first_event.toml";
   int K = 5, n_dyn = 20, M = 10;
   const char *work = ".";
   bool bridge_only = false;
   OptionsParser args(argc, argv);
   args.AddOption(&deck, "-d", "--deck", "Switching deck (ic scec_fd).");
   args.AddOption(&K, "-k", "--qd-steps-before", "QD steps before the dynamic phase.");
   args.AddOption(&n_dyn, "-nd", "--dyn-steps", "Dynamic steps before the switch back.");
   args.AddOption(&M, "-m", "--qd-steps-after", "QD steps after the switch back.");
   args.AddOption(&work, "-w", "--work-dir", "Directory of the test checkpoint.");
   args.AddOption(&bridge_only, "-bo", "--bridge-only", "-all", "--all",
                  "Only bridge_roundtrip.");
   args.Parse();
   if (!args.Good())
   {
      if (mpi.IsRoot()) { args.PrintUsage(std::cout); }
      return 1;
   }
#ifndef MFEM_USE_PETSC
   if (mpi.IsRoot()) { std::cerr << "test_cycle_p3 needs PETSc (the QD TS).\n"; }
   return 2;
#else
   const cycle::CycleConfig cfg = cycle::CycleConfig::Load(deck);
   MFEM_VERIFY(cfg.switching, deck << ": [cycle] switching is false; the P3 tests need "
               "the dynamic regime's sections");
   MFEM_VERIFY(K >= 1 && n_dyn >= 1 && M >= 1, "K, n and M must be >= 1");
   // The tests' W: damping and the normal-stress regularization off.
   cycle::CycleConfig cfg0 = cycle::CycleConfig::Load(deck);
   cfg0.dyn.kv_eta = 0.0;
   cfg0.dyn.pc = false;
   MFEMInitializePetsc(&argc, &argv, cfg.qd.petsc_options.c_str(), NULL);
   int failures = 0;
   {
      if (mpi.IsRoot()) { std::filesystem::create_directories(work); }
      mpi.Barrier();
      std::unique_ptr<Mesh> smesh = cycle::LoadScaledMesh(cfg.shared.mesh.path, cfg.mesh_scale);
      ParMesh pmesh(mpi.GetComm(), *smesh);
      smesh.reset();
      if (mpi.IsRoot())
      {
         std::cout << "test_cycle_p3: " << deck << ", " << cfg.shared.mesh.path << " x "
                   << cfg.mesh_scale << ", p = " << cfg.shared.mesh.order << ", ranks "
                   << mpi.Size() << ", commit " << SEAS_CYCLE_COMMIT << "\n";
      }

      // One QD regime at a time (each holds a factorization).
      auto qd = std::make_unique<cycle::QDPhase>(pmesh, mpi, cfg);
      const spatial::SpatialFrictionConfig &s = cfg.shared;
      const MaterialField material = MaterialField::MakeConstant(
         s.material_fallback.lambda, s.material_fallback.mu, s.material_fallback.rho);
      const spatial::PorePressureSpec no_pp;
      const cycle::FaultParams par = cycle::BuildFaultParams(
         qd->Table(), *s.rate_state, cfg.fields, material, pmesh, no_pp);
      const cycle::Prestress pre = cycle::BuildPrestress(par, cfg.prestress);
      Setup S{cfg, pmesh, mpi, par, pre};
      Checker c{mpi};
      TestCatalogMatchesController(c, work);
      {
         qd->Setup(par, pre);
         cycle::DynPhase dyn(pmesh, cfg0, qd->Table(), par, pre);
         TestBridgeRoundtrip(c, S, *qd, dyn);
      }
      qd.reset();
      if (!bridge_only) { TestQDRestartAfterSwitch(c, S, cfg0, K, n_dyn, M, work); }
      failures = c.failures;
      if (mpi.IsRoot())
      {
         std::cout << "\ntest_cycle_p3: " << (failures ? "FAILED" : "PASSED") << " ("
                   << failures << " failing checks)\n";
      }
   }
   MFEMFinalizePetsc();
   return failures ? 1 : 0;
#endif
}
