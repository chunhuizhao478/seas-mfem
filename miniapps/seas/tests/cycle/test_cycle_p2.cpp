// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_cycle_p2.cpp — fully dynamic BP5 cycle, P2 tests (plan rev 3,
// phase P2): the dynamic regime W (DynPhase, a NitscheWaveOperator) against
// the QD regime on the same mesh, boundaries, material and penalty.
//
//   frame_map       owned fault DOFs at the same points in both operators
//                   (DynPhase verifies it bitwise), the two kernel bases
//                   equal, FrameMap round trips (slip and traction), no
//                   normal component, eps_f n_k = n, and the ghost copies
//                   (R8): eps_f t_k of every local copy equals the owner's,
//                   and the slip load G delta through the QD path (kernel
//                   pairs copied to the ghost copies) equals W's (global
//                   vectors copied, each copy's own triad).
//   same_K          K u for a fixed u, F_D(t) at several t (side walls and
//                   the y = 0 creep plane, eq. 15) and the DOF coordinates of
//                   the two finite element spaces.
//   prestress_sign  W.PrestressTraction() = FrameMap(tau0, sigma0) exactly
//                   (hook H2); then the physics: W's own closure, fed the
//                   background traction alone with gamma = 0 and psi0,
//                   returns the QD V_init (mapped to global) in magnitude and
//                   direction, with sigma_eff = sigma0.  A flipped shear sign
//                   reverses V; a flipped normal makes sigma_eff negative.
//   load_clock      after Enter at t_s = 3.1e9 s W's clock reads 0, and W's
//                   acceleration from its own loads at clock tau equals
//                   M^-1 F_D of the QD operator at t_s + tau, bit for bit;
//                   F_D at tau alone differs (the control).
//   static_handoff  a QD static state with slip at t_s = 1e7 s: W's
//                   residual r = F_D - K u_s + G delta_s is at the QD solve
//                   tolerance, and W's traction T_h(u_s, delta_s) in the
//                   kernel frame equals the QD traction and normal traction.
//   energy          (--energy T > 0) Kelvin-Voigt and Prakash-Clifton off,
//                   from the first event's start (at t = 0 the walls and
//                   the slip are zero, so the QD balance is u = 0, delta = 0
//                   with V_init and psi0, and no QD solve is needed): the
//                   discrete energy
//                   balance of nitsche_energy.hpp over T at dt and dt/2.  The
//                   accumulated residual is O(dt^2) for this scheme (T-15 on
//                   the TPV102 box: ratio 3.5-4.5 per dt halving), so the
//                   check is that ratio, plus |sum rho| / work as a report.
//   timing          (--time-steps N > 0) N steps with the deck's damping from
//                   the background state; time per step.  No other test.
//
// --frames-only runs frame_map alone (np = 1, 2, 4, 8 in `make
// test-cycle-p2`); --split-y partitions the mesh so the two sides of the
// fault sit on different ranks and every fault face is shared (needs an even
// rank count).  The default partition is the drivers' (METIS).
//
// Usage (from miniapps/seas):
//   mpirun -np N ./seas_test_cycle_p2 [--deck config/fd_cycle/bp5_fd_first_event.toml]
//          [--frames-only] [--split-y] [--energy 2.0] [--time-steps 40]
// --energy and --time-steps each run alone.
// Exit code 0 = PASSED.

#include "mfem.hpp"
#include "../../cycle/scaled_mesh.hpp"
#include "../../cycle/cycle_config.hpp"
#include "../../cycle/fault_fields.hpp"
#include "../../cycle/prestress.hpp"
#include "../../cycle/qd_phase.hpp"
#include "../../cycle/dyn_phase.hpp"
#include "../../cycle/state_bridge.hpp"
#include "../../domain/nitsche_energy.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
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
   real_t MaxAbs(const Vector &x) const { return mpi.GlobalMax(x.Normlinf()); }
};

/// A smooth deterministic field of position, partition independent.
real_t Field(const real_t *x, int c)
{
   return std::sin(1.3e-4 * x[0] + 0.7 * c) * std::cos(2.1e-4 * x[2] - 0.3 * c)
          + 0.25 * std::sin(3.7e-5 * x[1] + 1.1 * c);
}

/// Kernel pairs (2 N) from Field on the owned fault DOF positions.
void KernelField(const Vector &X, real_t scale, Vector &k)
{
   const int N = X.Size() / 3;
   k.SetSize(2 * N);
   for (int j = 0; j < N; j++)
   {
      k(2 * j) = scale * Field(&X(3 * j), 0);
      k(2 * j + 1) = scale * Field(&X(3 * j), 1);
   }
}

/// Elements with centroid y < 0 on even ranks, y > 0 on odd ranks, np / 2
/// slabs in x: the fault plane y = 0 lies between ranks everywhere.
std::unique_ptr<int[]> SplitY(Mesh &m, int np)
{
   MFEM_VERIFY(np % 2 == 0, "--split-y needs an even rank count, not " << np);
   const int ns = np / 2;
   real_t xmin = 1e300, xmax = -1e300;
   for (int v = 0; v < m.GetNV(); v++)
   {
      xmin = std::min(xmin, m.GetVertex(v)[0]);
      xmax = std::max(xmax, m.GetVertex(v)[0]);
   }
   std::unique_ptr<int[]> part(new int[m.GetNE()]);
   Vector c(3);
   for (int e = 0; e < m.GetNE(); e++)
   {
      m.GetElementCenter(e, c);
      const int slab = std::min(ns - 1, static_cast<int>((c(0) - xmin) / (xmax - xmin) * ns));
      part[e] = 2 * slab + (c(1) > 0.0 ? 1 : 0);
   }
   return part;
}

/// The background state of the first event: zero displacement and slip, V
/// and psi from the IC set, at t = 0 (timing and prestress tests).
void EnterBackground(cycle::DynPhase &dyn, const cycle::FaultParams &par,
                     const cycle::Prestress &pre)
{
   Vector d(3 * dyn.NumOwnedFaultDOFs()), V, U(dyn.Op().Height());
   d = 0.0;
   U = 0.0;
   dyn.Frames().ToGlobal(par.V_init, V);
   dyn.Enter(U, d, V, pre.psi0, pre.sigma0, 0.0);
}

// ------------------------------------------------------------- frame_map

void TestFrameMap(Checker &c, cycle::QDPhase &qd, cycle::DynPhase &dyn)
{
   const MPIContext &mpi = c.mpi;
   NitscheWaveOperator &W = dyn.Op();
   ElasticityDomainOperator<ParMesh> &Eq = qd.Domain();
   ElasticityDomainOperator<ParMesh> &Ew = W.Elasticity();
   const cycle::FrameMap &f = dyn.Frames();
   const int N = dyn.NumOwnedFaultDOFs();
   if (mpi.IsRoot())
   {
      std::cout << "\nframe_map: " << qd.NumGlobalFaultDOFs() << " fault DOFs, "
                << dyn.NumSharedFaultFaces() << " shared fault faces\n";
   }

   // Same points (DynPhase verified it bitwise; repeated here as a report).
   c.Report("owned DOF positions, W vs QD table", c.RelInf(W.FaultPositions(), qd.Table().coords),
            0.0, "max rel diff");

   // The two kernel bases on every local fault DOF.
   DenseMatrix Bq, Bw;
   Eq.GetFaultDOFBasis(Bq);
   Ew.GetFaultDOFBasis(Bw);
   c.True("kernel basis sizes", Bq.Height() == Bw.Height() && Bq.Width() == Bw.Width(),
          std::to_string(Bw.Height()) + " x local fault DOFs");
   {
      Vector vq(Bq.GetData(), Bq.Height() * Bq.Width()), vw(Bw.GetData(), Bw.Height() * Bw.Width());
      c.Report("kernel basis (n_k, t1_k, t2_k), W vs QD", c.RelInf(vw, vq), 0.0, "max rel diff");
   }
   c.Report("ghost copies: eps_f t_k vs owner's (R8)", dyn.GhostFrameMismatch(), 1e-14, "max");
   c.Report("eps_f n_k vs oriented normal", dyn.NormalMismatch(), 1e-14, "max");

   // Round trips.
   Vector k, g, k2;
   KernelField(W.FaultPositions(), 1.0, k);
   f.ToGlobal(k, g);
   f.ToKernel(g, k2);
   c.Report("kernel -> global -> kernel", c.RelInf(k2, k), 1e-14, "max rel diff");
   real_t nmax = 0.0;
   for (int j = 0; j < N; j++)
   {
      nmax = std::max(nmax, std::abs(fault_frame::Dot(&g(3 * j), &W.FaultFrames()(9 * j))));
   }
   c.Report("normal component of global slip", mpi.GlobalMax(nmax) / c.MaxAbs(k), 1e-14,
            "max |delta.n| / max |s|");
   Vector sn(N), tau2, sn2, T;
   for (int j = 0; j < N; j++) { sn(j) = 5.0e7 + 1.0e6 * Field(&W.FaultPositions()(3 * j), 2); }
   Vector tau = k;
   tau *= 3.0e7;
   f.TractionToGlobal(tau, sn, T);
   f.TractionToKernel(T, tau2, sn2);
   c.Report("traction kernel -> global -> kernel (shear)", c.RelInf(tau2, tau), 1e-14,
            "max rel diff");
   c.Report("traction kernel -> global -> kernel (sigma_n)", c.RelInf(sn2, sn), 1e-14,
            "max rel diff");

   // R8 through the solvers: the slip load of the same slip, QD path vs W's.
   Vector s, sloc, Gq, Gw, d;
   KernelField(W.FaultPositions(), 2.0, s);
   Eq.ExpandOwnedToLocalFault(s, sloc, 2);
   Eq.AssembleSlipOnlyRHS(Gq, sloc);
   f.ToGlobal(s, d);
   W.SlipLoad(d, Gw);
   c.Report("slip load G delta, W vs QD (ghost copies)", c.RelInf(Gw, Gq), 1e-12,
            "max rel diff");
}

// --------------------------------------------------------------- same_K

void TestSameK(Checker &c, cycle::QDPhase &qd, cycle::DynPhase &dyn)
{
   const MPIContext &mpi = c.mpi;
   ElasticityDomainOperator<ParMesh> &Eq = qd.Domain();
   ElasticityDomainOperator<ParMesh> &Ew = dyn.Op().Elasticity();
   if (mpi.IsRoot()) { std::cout << "\nsame_K\n"; }

   // DOF coordinates: x projected into both spaces.
   ParFiniteElementSpace &fq = Eq.GetFESpace(), &fw = Ew.GetFESpace();
   c.True("vector sizes", fq.GetVSize() == fw.GetVSize() && fq.GetTrueVSize() == fw.GetTrueVSize()
          && fq.GetOrdering() == fw.GetOrdering(),
          "V-size " + std::to_string(fw.GetVSize()) + " on rank 0");
   {
      VectorFunctionCoefficient xc(3, [](const Vector &x, Vector &y) { y = x; });
      ParGridFunction xq(&fq), xw(&fw);
      xq.ProjectCoefficient(xc);
      xw.ProjectCoefficient(xc);
      c.Report("DOF coordinates", c.RelInf(xw, xq), 0.0, "max rel diff");
   }

   // K u for a fixed u (true-vector layout of both).
   const int n = fw.GetTrueVSize();
   Vector u(n), Kq, Kw;
   const HYPRE_BigInt off = fw.GetMyTDofOffset();
   for (int i = 0; i < n; i++) { u(i) = std::sin(0.731 * static_cast<real_t>(off + i)); }
   Eq.ApplyStiffness(u, Kq);
   Ew.ApplyStiffness(u, Kw);
   c.Report("K u, W vs QD", c.RelInf(Kw, Kq), 1e-12, "max rel diff");

   // F_D(t): W's clock is t - t0 (t0 = 0 before Enter).
   for (const real_t t : {0.0, 1.0, 3.15576e7, 1.0e10, 4.5e10})
   {
      Vector Fq, Fw;
      Eq.AssembleDirichletOnlyRHS(Fq, t);
      Ew.AssembleDirichletOnlyRHS(Fw, t - dyn.T0());
      std::ostringstream nm;
      nm << "F_D(t = " << std::setprecision(6) << t << " s), W vs QD";
      c.Report(nm.str(), c.RelInf(Fw, Fq), 1e-12, "max rel diff");
   }
}

// -------------------------------------------------------- prestress_sign

void TestPrestressSign(Checker &c, cycle::DynPhase &dyn, const cycle::FaultParams &par,
                       const cycle::Prestress &pre)
{
   const MPIContext &mpi = c.mpi;
   NitscheWaveOperator &W = dyn.Op();
   const int N = dyn.NumOwnedFaultDOFs();
   if (mpi.IsRoot()) { std::cout << "\nprestress_sign\n"; }

   Vector S0;
   dyn.Frames().TractionToGlobal(pre.tau0, pre.sigma0, S0);
   c.Report("H2: PrestressTraction() vs FrameMap(tau0, sigma0)",
            c.RelInf(W.PrestressTraction(), S0), 0.0, "max rel diff");

   // The closure, gamma = 0, on the background traction alone.
   Vector Vexp, Vgot(3 * N), sig(N);
   dyn.Frames().ToGlobal(par.V_init, Vexp);
   real_t vrel = 0.0, smis = 0.0, cosmin = 1.0;
   for (int j = 0; j < N; j++)
   {
      const real_t *n = &W.FaultFrames()(9 * j), *dv = &W.FaultFrames()(9 * j + 3),
                    *sv = &W.FaultFrames()(9 * j + 6);
      real_t X[2], s;
      fault_frame::ToLocal(&W.PrestressTraction()(3 * j), 0.0, n, dv, sv, X[0], X[1], s);
      FaultPointState st;
      st.psi = pre.psi0(j);
      st.a = W.FaultA()(j);
      st.Dc = W.FaultDc()(j);
      real_t Vt[2];
      W.Closure().SolveSlipRate(X, 0.0, s, st, Vt);
      fault_frame::Embed(Vt[0], Vt[1], dv, sv, &Vgot(3 * j));
      sig(j) = s;
      const real_t *ve = &Vexp(3 * j), *vg = &Vgot(3 * j);
      const real_t ne = std::sqrt(fault_frame::Dot(ve, ve)), ng = std::sqrt(fault_frame::Dot(vg, vg));
      real_t dd = 0.0;
      for (int i = 0; i < 3; i++) { dd = std::max(dd, std::abs(vg[i] - ve[i])); }
      vrel = std::max(vrel, dd / ne);
      cosmin = std::min(cosmin, fault_frame::Dot(ve, vg) / (ne * ng));
      smis = std::max(smis, std::abs(s - pre.sigma0(j)) / pre.sigma0(j));
   }
   c.Report("closure(S0, psi0, gamma 0): V vs V_init", mpi.GlobalMax(vrel), 1e-9,
            "max |dV| / |V_init|");
   c.Report("  direction: 1 - min cos(V, V_init)", 1.0 - mpi.GlobalMin(cosmin), 1e-12, "");
   c.Report("closure: sigma_eff vs sigma0", mpi.GlobalMax(smis), 1e-14, "max rel diff");
   (void)sig;
}

// ------------------------------------------------------------ load_clock

void TestLoadClock(Checker &c, cycle::QDPhase &qd, cycle::DynPhase &dyn,
                   const cycle::Prestress &pre)
{
   const MPIContext &mpi = c.mpi;
   NitscheWaveOperator &W = dyn.Op();
   const real_t ts = 3.1e9;
   if (mpi.IsRoot()) { std::cout << "\nload_clock (t_s = " << ts << " s)\n"; }
   const int n = W.Height(), N = dyn.NumOwnedFaultDOFs();
   Vector U(n), Z(n), d(3 * N), V(3 * N);
   U = 0.0;
   Z = 0.0;
   d = 0.0;
   V = 0.0;
   dyn.Enter(U, d, V, pre.psi0, pre.sigma0, ts);
   c.True("Enter: W's clock is 0", W.GetTime() == 0.0 && dyn.T0() == ts,
          "clock " + std::to_string(W.GetTime()));
   for (const real_t tau : {0.0, 1.0e-3, 0.37, 17.0})
   {
      W.SetTime(tau);
      Vector aw(n), aq(n), Fq, F0;
      W.Mult(U, Z, aw);   // M^-1 (F_D(clock) - K 0 + G 0)
      qd.Domain().AssembleDirichletOnlyRHS(Fq, ts + tau);
      W.Mass().Apply(Fq, aq);
      std::ostringstream nm;
      nm << "W at tau = " << tau << " vs QD at t_s + tau";
      c.Report(nm.str(), c.RelInf(aw, aq), 0.0, "max rel diff");
      if (tau == 1.0e-3)
      {
         qd.Domain().AssembleDirichletOnlyRHS(F0, tau);
         Vector a0(n);
         W.Mass().Apply(F0, a0);
         const real_t r = c.RelInf(aw, a0);
         c.True("control: W at tau vs QD at tau alone differs", r > 0.5,
                "max rel diff " + std::to_string(r));
      }
   }
   W.SetTime(0.0);
}

// -------------------------------------------------------- static_handoff

void TestStaticHandoff(Checker &c, cycle::QDPhase &qd, cycle::DynPhase &dyn)
{
   const MPIContext &mpi = c.mpi;
   NitscheWaveOperator &W = dyn.Op();
   ElasticityDomainOperator<ParMesh> &Ew = W.Elasticity();
   // The creeping region outside the fault slips Vp t_s while the fault
   // holds the field below, so the fault edge carries a slip step of Vp t_s:
   // at t_s = 1e9 s (1 m) the corner (-49.3 km, 39.3 km) opened.  1e7 s keeps
   // the step at 1 cm; load_clock covers large t.
   const real_t ts = 1.0e7;
   if (mpi.IsRoot()) { std::cout << "\nstatic_handoff (t_s = " << ts << " s, slip != 0)\n"; }

   // A QD state with slip: a smooth field of up to 2 m that vanishes at the
   // fault's edges (x = +-50 km, z = -40 km), then one solve.
   const int N = dyn.NumOwnedFaultDOFs();
   Vector s;
   KernelField(W.FaultPositions(), 2.0, s);
   for (int i = 0; i < N; i++)
   {
      const real_t x = W.FaultPositions()(3 * i) / 50.0e3, z = W.FaultPositions()(3 * i + 2) / 40.0e3;
      const real_t w = std::max(0.0, 1.0 - x * x) * std::max(0.0, 1.0 - z * z);
      s(2 * i) *= w;
      s(2 * i + 1) *= w;
   }
   Vector &st = qd.State();
   using FO = cycle::QDPhase::FaultOp;
   for (int i = 0; i < N; i++)
   {
      st(FO::StatePerNode * i + FO::SlipIndex) = s(2 * i);
      st(FO::StatePerNode * i + FO::SlipIndex + 1) = s(2 * i + 1);
   }
   qd.Solve(ts);

   Vector sk, psi_qd, sn, d, V, psi;
   qd.Slip(sk);
   qd.Psi(psi_qd);
   qd.SigmaN(sn);
   c.Report("QD slip read back", c.RelInf(sk, s), 0.0, "max rel diff");
   cycle::bridge::ToNitsche(dyn.Frames(), sk, qd.SlipRate(), psi_qd, d, V, psi);
   const Vector &us = qd.Displacement();
   dyn.Enter(us, d, V, psi, sn, ts);

   // r = F_D(t_s) - K u_s + G delta_s, and a = M^-1 r through W.Mult.
   const int n = W.Height();
   Vector Ku, F, Gd, r(n), Z(n), a(n), aK(n);
   Z = 0.0;
   Ew.ApplyStiffness(us, Ku);
   Ew.AssembleDirichletOnlyRHS(F, W.GetTime());
   W.SlipLoad(d, Gd);
   for (int i = 0; i < n; i++) { r(i) = F(i) - Ku(i) + Gd(i); }
   W.Mult(us, Z, a);
   W.Mass().Apply(Ku, aK);
   // Collective norms first: every rank calls MaxAbs.
   const real_t nKu = c.MaxAbs(Ku), nF = c.MaxAbs(F), nGd = c.MaxAbs(Gd), nr = c.MaxAbs(r),
                na = c.MaxAbs(a), naK = c.MaxAbs(aK);
   if (mpi.IsRoot())
   {
      std::cout << std::scientific << std::setprecision(3) << "  max |K u_s| = " << nKu
                << ", max |F_D| = " << nF << ", max |G delta| = " << nGd << "; max |M^-1 r| = "
                << na << " m/s^2 (max |M^-1 K u_s| = " << naK << ")\n";
   }
   c.Report("residual: max |r| / max |K u_s|", nr / nKu, 1e-10, "");
   c.Report("acceleration: max |M^-1 r| / max |M^-1 K u_s|", na / naK, 1e-10, "");

   // Traction: W's T_h(u_s, delta_s) in the kernel frame vs QD's.
   Vector T, tau_k, sn_k;
   W.TrialTraction(us, d, T);
   dyn.Frames().TractionToKernel(T, tau_k, sn_k);
   c.Report("traction (kernel shear), W vs QD", c.RelInf(tau_k, qd.Traction()), 1e-12,
            "max rel diff");
   c.Report("normal traction, W vs QD", c.RelInf(sn_k, qd.NormalTraction()), 1e-12,
            "max rel diff");
   const real_t ntau = c.MaxAbs(qd.Traction()), nsn = c.MaxAbs(qd.NormalTraction());
   if (mpi.IsRoot())
   {
      std::cout << std::scientific << std::setprecision(3) << "  max |tau| = " << ntau
                << " Pa, max |sigma_n change| = " << nsn << " Pa\n";
   }
}

// ---------------------------------------------------------------- energy

struct EnergyRun
{
   real_t dt = 0.0, sum_rho = 0.0, sum_abs = 0.0, Emax = 0.0;
   real_t W_ext = 0.0, W_pre = 0.0, W_fric = 0.0;
   NitscheEnergyRecord last;
   long steps = 0;
   double step_s = 0.0;
};

EnergyRun RunEnergy(ParMesh &pmesh, const cycle::CycleConfig &cfg, cycle::QDPhase &qd,
                    const cycle::FaultParams &par, const cycle::Prestress &pre,
                    real_t T, real_t dt_scale)
{
   EnergyRun R;
   cycle::DynPhase dyn(pmesh, cfg, qd.Table(), par, pre, dt_scale);
   EnterBackground(dyn, par, pre);
   NitscheWaveOperator &W = dyn.Op();
   NitscheEnergyMonitor mon(W);
   const int n = W.Height();
   Vector A0(n);
   W.Mult(dyn.U(), dyn.Ud(), A0);
   NitscheEnergyRecord r = mon.Evaluate(dyn.U(), dyn.Ud(), A0, 0.0);
   mon.Start(r);
   NitscheEnergyRecord prev = r;
   R.dt = dyn.Dt();
   R.steps = std::lround(std::ceil(T / R.dt - 1e-9));
   const double w0 = MPI_Wtime();
   for (long k = 0; k < R.steps; k++)
   {
      dyn.Step();
      r = mon.Evaluate(dyn.U(), dyn.Ud(), W.LastAcceleration(), dyn.Tau());
      mon.Step(r);
      R.W_ext += 0.5 * R.dt * (r.P_ext + prev.P_ext);
      R.W_pre += 0.5 * R.dt * (r.P_pre + prev.P_pre);
      R.W_fric += 0.5 * R.dt * (r.P_fric + prev.P_fric);
      R.Emax = std::max(R.Emax, std::abs(r.E_h));
      prev = r;
   }
   R.step_s = (MPI_Wtime() - w0) / std::max(R.steps, 1L);
   R.sum_rho = mon.SumRho();
   R.sum_abs = mon.SumAbsRho();
   R.last = r;
   return R;
}

void TestEnergy(Checker &c, ParMesh &pmesh, const cycle::CycleConfig &cfg0, cycle::QDPhase &qd,
                const cycle::FaultParams &par, const cycle::Prestress &pre, real_t T)
{
   const MPIContext &mpi = c.mpi;
   if (mpi.IsRoot())
   {
      std::cout << "\nenergy (Kelvin-Voigt and Prakash-Clifton off, " << T
                << " s from the first event's start)\n";
   }
   const EnergyRun a = RunEnergy(pmesh, cfg0, qd, par, pre, T, 1.0);
   const EnergyRun b = RunEnergy(pmesh, cfg0, qd, par, pre, T, 0.5);
   const real_t ratio = std::abs(a.sum_rho / b.sum_rho);
   const real_t scale = std::max({std::abs(a.W_pre), std::abs(a.W_fric), a.Emax});
   if (mpi.IsRoot())
   {
      std::cout << std::scientific << std::setprecision(6);
      for (const EnergyRun *e : {&a, &b})
      {
         std::cout << "  dt = " << e->dt << " s (" << e->steps << " steps, "
                   << std::setprecision(3) << e->step_s << " s/step with the monitor): E_h(T) = "
                   << std::setprecision(6) << e->last.E_h << " J [kin " << e->last.E_kin
                   << ", el " << e->last.E_el << "]\n     work: loads " << e->W_ext
                   << ", prestress " << e->W_pre << ", friction " << e->W_fric
                   << "; sum rho = " << e->sum_rho << " (sum |rho| " << e->sum_abs << ")\n";
      }
      std::cout << "  |sum rho| / work scale = " << std::abs(a.sum_rho) / scale << " at dt, "
                << std::abs(b.sum_rho) / scale << " at dt/2\n";
   }
   c.Report("accumulated residual ratio dt / (dt/2), 4 - r", std::abs(4.0 - ratio), 0.5,
            "|4 - " + std::to_string(ratio) + "| =");
}

// ---------------------------------------------------------------- timing

void Timing(const MPIContext &mpi, cycle::DynPhase &dyn, const cycle::FaultParams &par,
            const cycle::Prestress &pre, int nsteps)
{
   EnterBackground(dyn, par, pre);
   dyn.Step();   // the first step assembles K (and K_eta)
   const double w0 = MPI_Wtime();
   for (int k = 1; k < nsteps; k++) { dyn.Step(); }
   const double per = (MPI_Wtime() - w0) / std::max(nsteps - 1, 1);
   const real_t pmax = mpi.GlobalMax(per);
   if (mpi.IsRoot())
   {
      std::cout << "\ntiming: " << nsteps << " steps of " << dyn.Dt() << " s, Kelvin-Voigt "
                << (dyn.HasKelvinVoigt() ? "on" : "off") << ", Prakash-Clifton "
                << (dyn.HasPrakashClifton() ? "on" : "off") << ": " << per
                << " s/step (rank 0), " << pmax << " s/step (slowest rank), on " << mpi.Size()
                << " ranks; max V " << dyn.MaxV() << " m/s\n";
   }
}

} // namespace

int main(int argc, char *argv[])
{
   MPIContext mpi(&argc, &argv);
   const char *deck = "config/fd_cycle/bp5_fd_first_event.toml";
   bool frames_only = false, split_y = false;
   real_t energy_T = 0.0;
   int time_steps = 0;
   OptionsParser args(argc, argv);
   args.AddOption(&deck, "-d", "--deck", "P2 deck (switching on).");
   args.AddOption(&frames_only, "-fo", "--frames-only", "-all", "--all",
                  "Only frame_map (the multi-rank check).");
   args.AddOption(&split_y, "-sy", "--split-y", "-metis", "--metis",
                  "Partition y < 0 and y > 0 onto different ranks (every fault face shared).");
   args.AddOption(&energy_T, "-e", "--energy", "Energy balance over this time [s] (0 = skip).");
   args.AddOption(&time_steps, "-ts", "--time-steps",
                  "Only time this many steps with the deck's damping (0 = skip).");
   args.Parse();
   if (!args.Good())
   {
      if (mpi.IsRoot()) { args.PrintUsage(std::cout); }
      return 1;
   }

   const cycle::CycleConfig cfg = cycle::CycleConfig::Load(deck);
   MFEM_VERIFY(cfg.switching, deck << ": [cycle] switching is false; the P2 tests need "
               "the dynamic regime's sections");
   // The tests' W: damping and the normal-stress regularization off.
   cycle::CycleConfig cfg0 = cycle::CycleConfig::Load(deck);
   cfg0.dyn.kv_eta = 0.0;
   cfg0.dyn.pc = false;

   std::unique_ptr<Mesh> smesh = cycle::LoadScaledMesh(cfg.shared.mesh.path, cfg.mesh_scale);
   std::unique_ptr<int[]> part;
   if (split_y) { part = SplitY(*smesh, mpi.Size()); }
   ParMesh pmesh(mpi.GetComm(), *smesh, part.get());
   smesh.reset();
   if (mpi.IsRoot())
   {
      std::cout << "test_cycle_p2: " << deck << ", " << cfg.shared.mesh.path << " x "
                << cfg.mesh_scale << ", p = " << cfg.shared.mesh.order << ", ranks "
                << mpi.Size() << ", partition " << (split_y ? "split at y = 0" : "METIS")
                << ", commit " << SEAS_CYCLE_COMMIT << "\n";
   }

   cycle::QDPhase qd(pmesh, mpi, cfg);
   const spatial::SpatialFrictionConfig &s = cfg.shared;
   const MaterialField material = MaterialField::MakeConstant(
      s.material_fallback.lambda, s.material_fallback.mu, s.material_fallback.rho);
   const spatial::PorePressureSpec no_pp;
   const cycle::FaultParams par = cycle::BuildFaultParams(
      qd.Table(), *s.rate_state, cfg.fields, material, pmesh, no_pp);
   const cycle::Prestress pre = cycle::BuildPrestress(par, cfg.prestress);

   if (time_steps > 0)
   {
      cycle::DynPhase dyn(pmesh, cfg, qd.Table(), par, pre);
      Timing(mpi, dyn, par, pre, time_steps);
      return 0;
   }

   Checker c{mpi};
   if (energy_T > 0.0)
   {
      TestEnergy(c, pmesh, cfg0, qd, par, pre, energy_T);
   }
   else
   {
      {
         cycle::DynPhase dyn(pmesh, cfg0, qd.Table(), par, pre);
         TestFrameMap(c, qd, dyn);
         if (!frames_only)
         {
            TestSameK(c, qd, dyn);
            TestPrestressSign(c, dyn, par, pre);
            TestLoadClock(c, qd, dyn, pre);
         }
      }
      if (!frames_only)
      {
         qd.Setup(par, pre);   // the initial QD solve (factorization)
         cycle::DynPhase dyn(pmesh, cfg0, qd.Table(), par, pre);
         TestStaticHandoff(c, qd, dyn);
      }
   }

   if (mpi.IsRoot())
   {
      std::cout << "\ntest_cycle_p2" << (frames_only ? " (frame_map)" : "") << ": "
                << (c.failures ? "FAILED" : "PASSED") << " (" << c.failures
                << " failing checks)\n";
   }
   return c.failures ? 1 : 0;
}
