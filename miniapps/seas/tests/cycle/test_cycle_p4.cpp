// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_cycle_p4.cpp — fully dynamic BP5 cycle, P4 tests (plan rev 3, phase
// P4): the QD -> dynamic switch at t > 0 (plan 6.2), with one dynamic
// regime for the whole run.
//
//   handoff_closure    The QD balance at t_s = 1e7 s with slip (up to 2 m)
//                      and a perturbed psi, then the driver's hand-over
//                      (Bridge::ToNitsche, DynPhase::Enter: U = u_s, U' = 0).
//                      At every fault DOF W's closure X - gamma V = F(V, psi)
//                      with X = T_h(u_s, delta_s) + S0 (the traction with the
//                      bulk at rest) and gamma = eta returns the QD slip
//                      rate: the QD balance tau0 + dtau - eta V = F is W's
//                      closure with the radiation-damping impedance.  That
//                      is the rule a switch at t > 0 and the eta_V_init start
//                      at t = 0 share.  The first step's own gamma
//                      (sigma_F dt / 2, plus the Kelvin-Voigt slope) is not
//                      eta, so |V^1| / |V_qd| is printed only.
//   reenter_repeat     The driver keeps one DynPhase and re-enters it at
//                      every switch.  Enter state A, n steps; Enter state B
//                      (other slip, psi and t), n steps; Enter A again from
//                      the same inputs, n steps.  The third phase must repeat
//                      the first bit for bit: U, U', slip, slip rate, psi,
//                      sigma* (Prakash-Clifton), time and steps, and the
//                      station fields right after each Enter of A (before the
//                      first step) and after its n steps.  With the deck's
//                      damping.  No solve between, so bit for bit on any
//                      number of ranks.
//
// Usage (from miniapps/seas; PETSc needed):
//   mpirun -np N ./seas_test_cycle_p4 [--deck config/fd_cycle/bp5_fd_first_event.toml]
//          [--dyn-steps 20]
// Exit code 0 = PASSED.

#include "mfem.hpp"
#include "../../cycle/scaled_mesh.hpp"
#include "../../cycle/cycle_config.hpp"
#include "../../cycle/fault_fields.hpp"
#include "../../cycle/prestress.hpp"
#include "../../cycle/qd_phase.hpp"
#include "../../cycle/dyn_phase.hpp"
#include "../../cycle/state_bridge.hpp"
#include "../../cycle/fd_cycle_output.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace
{

struct Checker
{
   const MPIContext &mpi;
   int failures = 0;

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

/// A QD fault state on owned DOFs: slip (kernel pairs) of up to `amp` m that
/// vanishes at the fault's edges (x = +-50 km, z = -40 km), and psi0 plus
/// `dpsi` times a smooth field.  `seed` picks the fields.
void QDFaultState(const cycle::DynPhase &dyn, const cycle::Prestress &pre, real_t amp,
                  real_t dpsi, int seed, Vector &s, Vector &psi)
{
   const Vector &X = dyn.FaultPositions();
   const int N = dyn.NumOwnedFaultDOFs();
   s.SetSize(2 * N);
   psi.SetSize(N);
   for (int i = 0; i < N; i++)
   {
      const real_t x = X(3 * i) / 50.0e3, z = X(3 * i + 2) / 40.0e3;
      const real_t w = std::max(0.0, 1.0 - x * x) * std::max(0.0, 1.0 - z * z);
      s(2 * i) = amp * w * Field(&X(3 * i), seed);
      s(2 * i + 1) = amp * w * Field(&X(3 * i), seed + 1);
      psi(i) = pre.psi0(i) + dpsi * Field(&X(3 * i), seed + 2);
   }
}

/// The driver's hand-over inputs from the QD balance at t (plan 6.2).
struct Handover
{
   Vector u, d, V, psi, sn;
   real_t t = 0.0, Vmax = 0.0;
};

Handover QDBalance(cycle::QDPhase &qd, const cycle::DynPhase &dyn, const Vector &s,
                   const Vector &psi0, real_t t)
{
   Handover h;
   qd.SetFaultState(s, psi0);
   h.Vmax = qd.Solve(t);
   Vector sk, psi_qd;
   qd.Slip(sk);
   qd.Psi(psi_qd);
   qd.SigmaN(h.sn);
   cycle::bridge::ToNitsche(dyn.Frames(), sk, qd.SlipRate(), psi_qd, h.d, h.V, h.psi);
   h.u = qd.Displacement();
   h.t = t;
   return h;
}

void Enter(cycle::DynPhase &dyn, const Handover &h)
{
   dyn.Enter(h.u, h.d, h.V, h.psi, h.sn, h.t);
}

// ------------------------------------------------------------ handoff_closure

void TestHandoffClosure(Checker &c, cycle::QDPhase &qd, cycle::DynPhase &dyn,
                        const cycle::FaultParams &par, const cycle::Prestress &pre)
{
   const MPIContext &mpi = c.mpi;
   NitscheWaveOperator &W = dyn.Op();
   // 1e7 s keeps the slip step at the fault's edge at 1 cm (test_cycle_p2).
   const real_t ts = 1.0e7;
   if (mpi.IsRoot()) { std::cout << "\nhandoff_closure (t_s = " << ts << " s, slip != 0)\n"; }
   Vector s0, psi0;
   QDFaultState(dyn, pre, 2.0, 0.01, 0, s0, psi0);
   const Handover h = QDBalance(qd, dyn, s0, psi0, ts);
   Enter(dyn, h);

   // X = T_h(u_s, delta_s) + S0 at every owned fault DOF, then the closure.
   Vector T;
   W.TrialTraction(h.u, h.d, T);
   const Vector &S0 = W.PrestressTraction();
   const int N = dyn.NumOwnedFaultDOFs();
   real_t rel = 0.0, cosmin = 1.0, vmin = std::numeric_limits<real_t>::infinity(), vmax = 0.0;
   for (int j = 0; j < N; j++)
   {
      const real_t *n = &W.FaultFrames()(9 * j), *dv = &W.FaultFrames()(9 * j + 3),
                    *sv = &W.FaultFrames()(9 * j + 6);
      real_t TS[3], X[2], s, Vt[2], vg[3];
      for (int k = 0; k < 3; k++) { TS[k] = T(3 * j + k) + S0(3 * j + k); }
      fault_frame::ToLocal(TS, 0.0, n, dv, sv, X[0], X[1], s);
      FaultPointState st;
      st.psi = h.psi(j);
      st.a = W.FaultA()(j);
      st.Dc = W.FaultDc()(j);
      W.Closure().SolveSlipRate(X, par.rs.eta(j), s, st, Vt);
      fault_frame::Embed(Vt[0], Vt[1], dv, sv, vg);
      const real_t *vq = &h.V(3 * j);
      const real_t nq = std::sqrt(fault_frame::Dot(vq, vq)),
                   ng = std::sqrt(fault_frame::Dot(vg, vg));
      real_t dd = 0.0;
      for (int k = 0; k < 3; k++) { dd = std::max(dd, std::abs(vg[k] - vq[k])); }
      rel = std::max(rel, dd / nq);
      cosmin = std::min(cosmin, fault_frame::Dot(vq, vg) / (nq * ng));
      vmin = std::min(vmin, nq);
      vmax = std::max(vmax, nq);
   }
   vmin = mpi.GlobalMin(vmin);
   vmax = mpi.GlobalMax(vmax);
   if (mpi.IsRoot())
   {
      std::cout << std::scientific << std::setprecision(3) << "  QD balance: |V| from " << vmin
                << " to " << vmax << " m/s\n";
   }
   c.Report("closure(T_h(u_s, d_s) + S0, psi, gamma eta) vs V_qd", mpi.GlobalMax(rel), 1e-8,
            "max |dV| / |V_qd| per DOF");
   c.Report("  direction: 1 - min cos(V, V_qd)", 1.0 - mpi.GlobalMin(cosmin), 1e-10, "");

   // One step with W's own gamma: printed only.
   dyn.Step();
   Vector v1;
   dyn.SlipRateMagnitude(v1);
   real_t rmin = std::numeric_limits<real_t>::infinity(), rmax = 0.0;
   for (int j = 0; j < N; j++)
   {
      const real_t nq = std::sqrt(fault_frame::Dot(&h.V(3 * j), &h.V(3 * j)));
      rmin = std::min(rmin, v1(j) / nq);
      rmax = std::max(rmax, v1(j) / nq);
   }
   rmin = mpi.GlobalMin(rmin);
   rmax = mpi.GlobalMax(rmax);
   if (mpi.IsRoot())
   {
      std::cout << std::fixed << std::setprecision(4) << "  after one step (damping off): "
                << "|V^1| / |V_qd| from " << rmin << " to " << rmax << " (printed only)\n";
   }
}

// ------------------------------------------------------------- reenter_repeat

struct DynSnapshot
{
   Vector U, Ud, slip, V, psi, sigpc;
   cycle::StationFields f;
   real_t time = 0.0, tau = 0.0;
   long steps = 0;
};

DynSnapshot Snap(cycle::DynPhase &dyn)
{
   DynSnapshot s;
   s.U = dyn.U();
   s.Ud = dyn.Ud();
   s.slip = dyn.Slip();
   s.V = dyn.Op().SlipRate();
   s.psi = dyn.Psi();
   if (dyn.HasPrakashClifton()) { s.sigpc = dyn.Op().PrakashCliftonSigma(); }
   dyn.FillStationFields(s.f);
   s.time = dyn.Time();
   s.tau = dyn.Tau();
   s.steps = dyn.Steps();
   return s;
}

void Same(Checker &c, const std::string &when, const DynSnapshot &x, const DynSnapshot &a,
          bool pc)
{
   c.Report(when + ": U", c.RelInf(x.U, a.U), 0.0, "max rel diff");
   c.Report(when + ": U'", c.RelInf(x.Ud, a.Ud), 0.0, "max rel diff");
   c.Report(when + ": slip", c.RelInf(x.slip, a.slip), 0.0, "max rel diff");
   c.Report(when + ": slip rate", c.RelInf(x.V, a.V), 0.0, "max rel diff");
   c.Report(when + ": psi", c.RelInf(x.psi, a.psi), 0.0, "max rel diff");
   if (pc) { c.Report(when + ": sigma* (Prakash-Clifton)", c.RelInf(x.sigpc, a.sigpc), 0.0,
                      "max rel diff"); }
   real_t fd = 0.0;
   const std::vector<const Vector *> fx = {&x.f.slip_dip, &x.f.slip_strike, &x.f.V_dip,
                                           &x.f.V_strike, &x.f.tau_a_dip, &x.f.tau_a_strike,
                                           &x.f.theta};
   const std::vector<const Vector *> fa = {&a.f.slip_dip, &a.f.slip_strike, &a.f.V_dip,
                                           &a.f.V_strike, &a.f.tau_a_dip, &a.f.tau_a_strike,
                                           &a.f.theta};
   for (size_t k = 0; k < fx.size(); k++) { fd = std::max(fd, c.RelInf(*fx[k], *fa[k])); }
   c.Report(when + ": station fields", fd, 0.0, "max rel diff (7 fields)");
   c.True(when + ": time, tau, steps", x.time == a.time && x.tau == a.tau && x.steps == a.steps,
          "t " + std::to_string(x.time) + " s, " + std::to_string(x.steps) + " steps");
}

void TestReenterRepeat(Checker &c, cycle::QDPhase &qd, cycle::DynPhase &dyn,
                       const cycle::Prestress &pre, int n)
{
   const MPIContext &mpi = c.mpi;
   if (mpi.IsRoot())
   {
      std::cout << "\nreenter_repeat (the deck's damping; " << n << " steps per phase)\n";
   }
   Vector sA, pA, sB, pB;
   QDFaultState(dyn, pre, 2.0, 0.01, 0, sA, pA);
   QDFaultState(dyn, pre, -1.0, -0.01, 3, sB, pB);
   const Handover A = QDBalance(qd, dyn, sA, pA, 1.0e7);
   const Handover B = QDBalance(qd, dyn, sB, pB, 2.0e7);

   Enter(dyn, A);
   const DynSnapshot A0 = Snap(dyn);
   for (int k = 0; k < n; k++) { dyn.Step(); }
   const DynSnapshot A1 = Snap(dyn);
   Enter(dyn, B);
   for (int k = 0; k < n; k++) { dyn.Step(); }
   const DynSnapshot B1 = Snap(dyn);
   Enter(dyn, A);
   const DynSnapshot C0 = Snap(dyn);
   for (int k = 0; k < n; k++) { dyn.Step(); }
   const DynSnapshot C1 = Snap(dyn);

   // Control: B's phase is a different one, so equal results are not trivial.
   const real_t dU = c.RelInf(B1.U, A1.U), dV = c.RelInf(B1.V, A1.V);
   if (mpi.IsRoot())
   {
      std::cout << std::scientific << std::setprecision(3) << "  max V at the hand-over: A "
                << A.Vmax << ", B " << B.Vmax << " m/s; control, B vs A after " << n
                << " steps: U " << dU << ", slip rate " << dV << " (max rel diff)\n";
   }
   c.True("control: B differs from A", dU > 1e-6 && dV > 1e-6, "U and slip rate");
   const bool pc = dyn.HasPrakashClifton();
   Same(c, "A again, at Enter", C0, A0, pc);
   Same(c, "A again, after the steps", C1, A1, pc);
}

} // namespace

int main(int argc, char *argv[])
{
   MPIContext mpi(&argc, &argv);
   const char *deck = "config/fd_cycle/bp5_fd_first_event.toml";
   int n_dyn = 20;
   OptionsParser args(argc, argv);
   args.AddOption(&deck, "-d", "--deck", "Switching deck.");
   args.AddOption(&n_dyn, "-nd", "--dyn-steps", "Dynamic steps per phase (reenter_repeat).");
   args.Parse();
   if (!args.Good())
   {
      if (mpi.IsRoot()) { args.PrintUsage(std::cout); }
      return 1;
   }
#ifndef MFEM_USE_PETSC
   if (mpi.IsRoot()) { std::cerr << "test_cycle_p4 needs PETSc (the QD TS).\n"; }
   return 2;
#else
   const cycle::CycleConfig cfg = cycle::CycleConfig::Load(deck);
   MFEM_VERIFY(cfg.switching, deck << ": [cycle] switching is false; the P4 tests need "
               "the dynamic regime's sections");
   MFEM_VERIFY(n_dyn >= 1, "--dyn-steps must be >= 1");
   // handoff_closure's W: damping and the normal-stress regularization off.
   cycle::CycleConfig cfg0 = cycle::CycleConfig::Load(deck);
   cfg0.dyn.kv_eta = 0.0;
   cfg0.dyn.pc = false;
   MFEMInitializePetsc(&argc, &argv, cfg.qd.petsc_options.c_str(), NULL);
   int failures = 0;
   {
      std::unique_ptr<Mesh> smesh = cycle::LoadScaledMesh(cfg.shared.mesh.path, cfg.mesh_scale);
      ParMesh pmesh(mpi.GetComm(), *smesh);
      smesh.reset();
      if (mpi.IsRoot())
      {
         std::cout << "test_cycle_p4: " << deck << ", " << cfg.shared.mesh.path << " x "
                   << cfg.mesh_scale << ", p = " << cfg.shared.mesh.order << ", ranks "
                   << mpi.Size() << ", IC set " << cfg.IcLabel() << ", commit "
                   << SEAS_CYCLE_COMMIT << "\n";
      }
      cycle::QDPhase qd(pmesh, mpi, cfg);
      const spatial::SpatialFrictionConfig &s = cfg.shared;
      const MaterialField material = MaterialField::MakeConstant(
         s.material_fallback.lambda, s.material_fallback.mu, s.material_fallback.rho);
      const spatial::PorePressureSpec no_pp;
      const cycle::FaultParams par = cycle::BuildFaultParams(
         qd.Table(), *s.rate_state, cfg.fields, material, pmesh, no_pp);
      const cycle::Prestress pre = cycle::BuildPrestress(par, cfg.prestress);
      qd.Setup(par, pre);
      Checker c{mpi};
      {
         cycle::DynPhase dyn(pmesh, cfg0, qd.Table(), par, pre);
         TestHandoffClosure(c, qd, dyn, par, pre);
      }
      {
         cycle::DynPhase dyn(pmesh, cfg, qd.Table(), par, pre);
         TestReenterRepeat(c, qd, dyn, pre, n_dyn);
      }
      failures = c.failures;
      if (mpi.IsRoot())
      {
         std::cout << "\ntest_cycle_p4: " << (failures ? "FAILED" : "PASSED") << " ("
                   << failures << " failing checks)\n";
      }
   }
   MFEMFinalizePetsc();
   return failures ? 1 : 0;
#endif
}
