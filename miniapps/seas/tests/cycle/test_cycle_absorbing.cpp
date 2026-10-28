// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_cycle_absorbing.cpp — fully dynamic BP5 cycle, P5 tests (plan rev 3,
// phase P5, Sections 6.4 and 11.2): the dynamic regime's absorbing outer
// faces ([cycle.absorbing], cycle/absorbing_faces.hpp, hook H3).
//
// One QD balance at t_s = 3e9 s is handed over (DynPhase::Enter) to three
// dynamic regimes
// built in turn from the deck, with its damping: reflecting (faces = "off",
// the reference), "outer" and "bottom".  The method is DG, so "the wall
// elements" (the elements with a wall face) own their DOFs, and the walls'
// Nitsche terms and the dashpots touch only those.  Accelerations are
// measured against a_K = max |M^-1 K u_s| of the reflecting W.
//
//   loaded_state       The fault state: the walls and the creep plane at
//                      +-Vp t_s / 2 (1.5 m); strike slip s_p (1 - w / 2),
//                      w = (1 - x^2 / 50^2 km^2)(1 - z^2 / 40^2 km^2) (a locked
//                      patch with half the plate slip, the slip continuous
//                      with the creep plane's at the fault's edges), dip slip
//                      0.1 m w f(x), psi0 + 0.01 f(x).  s_p = +-Vp t_s: the
//                      sign for which a uniform strike slip s_p with psi0 (the
//                      two halves moved with the walls) leaves the balance's
//                      max V as at t = 0 (the IC) to 1e-6; exactly one sign
//                      must (the kernel strike's sense depends on the deck's
//                      ref_normal).
//   retag              "outer" on the deck's mesh: the re-tagged face
//                      counts (on bp5_tandem_exact.msh 138 bottom, 240 wall
//                      and 1876 free-surface faces), the QD mesh's boundary
//                      attributes unchanged, W's absorbing attributes.
//   frozen_load        b = K_abs u_s - G d_s - F_c(t_s) on the wall elements
//                      is not zero there and is zero elsewhere; it is the
//                      walls' reaction, so its net force sum_i b_i,c in each
//                      direction c is zero (the body is in equilibrium:
//                      traction-free top and bottom, the fault and the creep
//                      plane internal) up to the QD solve residual:
//                      |sum b_c| / sum |b_c| <= 1e-6.
//   switch_equilibrium The first acceleration A^0 of the "outer" W, with b:
//                      off the wall elements it equals the reflecting W's
//                      (with Kelvin-Voigt both hold M^-1 G_eta V_s; plan
//                      11.1), on them it is 0, both to 1e-12 a_K.  The
//                      reflecting A^0 on them, M^-1 r (the QD residual the
//                      frozen load leaves out), is printed.  Control: without
//                      b, A^0 on the wall elements differs from the
//                      reflecting one by about max |M^-1 b| (> 0.5 of it),
//                      the static load the walls would drop.
//   first_steps        n steps of each: the slip, slip rate and psi agree
//                      with the reflecting run to 1e-12 (no wave has reached
//                      a boundary).  U is printed only: the reflecting walls
//                      follow their load, Vp tau / 2, the frozen ones do not.
//   bottom             "bottom": no wall faces and no frozen load, A^0 equal
//                      to the reflecting one to 1e-13 a_K (C U' = 0), the
//                      slip rate and U after n steps to 1e-12.
//
// Usage (from miniapps/seas; PETSc needed):
//   mpirun -np N ./seas_test_cycle_absorbing [--deck config/fd_cycle/bp5_fd_first_event.toml]
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

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <set>
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

   void Report(const std::string &name, real_t m, real_t tol, const std::string &what)
   {
      const bool ok = (m <= tol);
      if (!ok) { failures++; }
      if (mpi.IsRoot())
      {
         std::cout << "  " << (ok ? "ok  " : "FAIL") << " " << std::left << std::setw(46)
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
         std::cout << "  " << (ok ? "ok  " : "FAIL") << " " << std::left << std::setw(46)
                   << name << std::right << " " << detail << "\n";
      }
   }
   real_t MaxAbs(const Vector &x) const
   {
      real_t m = 0.0;
      for (int i = 0; i < x.Size(); i++) { m = std::max(m, std::abs(x(i))); }
      return mpi.GlobalMax(m);
   }
   /// max |x - y| over all ranks.
   real_t DiffInf(const Vector &x, const Vector &y) const
   {
      MFEM_VERIFY(x.Size() == y.Size(), "DiffInf: sizes " << x.Size() << " vs " << y.Size());
      real_t d = 0.0;
      for (int i = 0; i < x.Size(); i++) { d = std::max(d, std::abs(x(i) - y(i))); }
      return mpi.GlobalMax(d);
   }
   /// max |x_i| over the local i with mask[i] == want, all ranks.
   real_t MaxAbsOn(const Vector &x, const std::vector<char> &mask, char want) const
   {
      real_t m = 0.0;
      for (int i = 0; i < x.Size(); i++)
      {
         if (mask[i] == want) { m = std::max(m, std::abs(x(i))); }
      }
      return mpi.GlobalMax(m);
   }
   /// max |x_i - y_i| over the local i with mask[i] == want, all ranks.
   real_t DiffInfOn(const Vector &x, const Vector &y, const std::vector<char> &mask,
                    char want) const
   {
      MFEM_VERIFY(x.Size() == y.Size(), "DiffInfOn: sizes " << x.Size() << " vs " << y.Size());
      real_t d = 0.0;
      for (int i = 0; i < x.Size(); i++)
      {
         if (mask[i] == want) { d = std::max(d, std::abs(x(i) - y(i))); }
      }
      return mpi.GlobalMax(d);
   }
   /// max |x - y| / max |y| over all ranks (0 when both are zero).
   real_t RelInf(const Vector &x, const Vector &y) const
   {
      const real_t d = DiffInf(x, y), s = MaxAbs(y);
      return (s > 0.0) ? d / s : d;
   }
};

/// A smooth deterministic field of position, partition independent.
real_t Field(const real_t *x, int c)
{
   return std::sin(1.3e-4 * x[0] + 0.7 * c) * std::cos(2.1e-4 * x[2] - 0.3 * c)
          + 0.25 * std::sin(3.7e-5 * x[1] + 1.1 * c);
}

/// The loaded fault state (see loaded_state): strike slip S (1 - f w), dip slip
/// `dip` w f_1(x), psi0 + dpsi f_2(x); w = 0 at the fault's edges (x = +-50 km,
/// z = -40 km).
void LoadedFaultState(const cycle::QDPhase &qd, const cycle::Prestress &pre, real_t S,
                      real_t f, real_t dip, real_t dpsi, Vector &s, Vector &psi)
{
   const Vector &X = qd.Table().coords;
   const int N = qd.Table().Size();
   s.SetSize(2 * N);
   psi.SetSize(N);
   for (int i = 0; i < N; i++)
   {
      const real_t x = X(3 * i) / 50.0e3, z = X(3 * i + 2) / 40.0e3;
      const real_t w = std::max(0.0, 1.0 - x * x) * std::max(0.0, 1.0 - z * z);
      s(2 * i) = dip * w * Field(&X(3 * i), 1);      // kernel pairs: (dip, strike)
      s(2 * i + 1) = S * (1.0 - f * w);
      psi(i) = pre.psi0(i) + dpsi * Field(&X(3 * i), 2);
   }
}

/// The QD balance, solved once; the bridge needs a DynPhase's frames.
struct Balance
{
   Vector u, sk, Vqd, psi, sn;
   real_t t = 0.0, Vmax = 0.0;
};

struct Handover
{
   Vector u, d, V, psi, sn;
   real_t t = 0.0;
};

Handover ToNitsche(const Balance &b, const cycle::DynPhase &dyn)
{
   Handover h;
   cycle::bridge::ToNitsche(dyn.Frames(), b.sk, b.Vqd, b.psi, h.d, h.V, h.psi);
   h.u = b.u;
   h.sn = b.sn;
   h.t = b.t;
   return h;
}

/// What a run leaves for the comparisons.
struct Run
{
   Vector A0, U, slip, V, psi;
   Vector Minv_b;   ///< M^-1 b ("outer" only)
   real_t tau = 0.0;   ///< W's clock after the steps (time since Enter)
};

/// Enter, A^0 = W.Mult(u_s, 0), n steps.
Run Phase(cycle::DynPhase &dyn, const Balance &b, int n)
{
   Run r;
   const Handover h = ToNitsche(b, dyn);
   dyn.Enter(h.u, h.d, h.V, h.psi, h.sn, h.t);
   Vector zero(h.u.Size());
   zero = 0.0;
   dyn.Op().Mult(h.u, zero, r.A0);
   if (dyn.FrozenWallLoad().Size() > 0)
   {
      dyn.Op().Mass().Apply(dyn.FrozenWallLoad(), r.Minv_b);
   }
   for (int k = 0; k < n; k++) { dyn.Step(); }
   r.U = dyn.U();
   r.slip = dyn.Slip();
   r.V = dyn.Op().SlipRate();
   r.psi = dyn.Psi();
   r.tau = dyn.Op().GetTime();
   return r;
}

/// 1 on the local DOFs of the elements with a boundary face of attribute
/// `attr` (DG: the element's own DOFs), 0 elsewhere.
std::vector<char> ElementDofMask(NitscheWaveOperator &W, int attr)
{
   ParMesh &m = W.GetMesh();
   const ParFiniteElementSpace &fes = W.Elasticity().GetFESpace();
   std::vector<char> on(fes.GetVSize(), 0);
   Array<int> vd;
   for (int be = 0; be < m.GetNBE(); be++)
   {
      if (m.GetBdrAttribute(be) != attr) { continue; }
      int e = -1, info = 0;
      m.GetBdrElementAdjacentElement(be, e, info);
      fes.GetElementVDofs(e, vd);
      for (int i = 0; i < vd.Size(); i++) { on[vd[i] >= 0 ? vd[i] : -1 - vd[i]] = 1; }
   }
   return on;
}

std::vector<int> BdrAttributes(const ParMesh &m)
{
   std::vector<int> a(m.GetNBE());
   for (int be = 0; be < m.GetNBE(); be++) { a[be] = m.GetBdrAttribute(be); }
   return a;
}

} // namespace

int main(int argc, char *argv[])
{
   MPIContext mpi(&argc, &argv);
   const char *deck = "config/fd_cycle/bp5_fd_first_event.toml";
   int n_dyn = 20;
   OptionsParser args(argc, argv);
   args.AddOption(&deck, "-d", "--deck", "Switching deck.");
   args.AddOption(&n_dyn, "-nd", "--dyn-steps", "Dynamic steps per run.");
   args.Parse();
   if (!args.Good())
   {
      if (mpi.IsRoot()) { args.PrintUsage(std::cout); }
      return 1;
   }
#ifndef MFEM_USE_PETSC
   if (mpi.IsRoot()) { std::cerr << "test_cycle_absorbing needs PETSc (the QD TS).\n"; }
   return 2;
#else
   const cycle::CycleConfig cfg = cycle::CycleConfig::Load(deck);
   MFEM_VERIFY(cfg.switching, deck << ": [cycle] switching is false; the P5 tests need "
               "the dynamic regime's sections");
   MFEM_VERIFY(!cfg.lts.On(), deck << ": the P5 tests compare global stepping runs");
   MFEM_VERIFY(n_dyn >= 1, "--dyn-steps must be >= 1");
   cycle::CycleConfig cfg_off = cycle::CycleConfig::Load(deck);
   cfg_off.absorbing.faces = "off";
   cycle::CycleConfig cfg_outer = cfg_off, cfg_bottom = cfg_off;
   cfg_outer.absorbing.faces = "outer";
   cfg_bottom.absorbing.faces = "bottom";
   MFEMInitializePetsc(&argc, &argv, cfg.qd.petsc_options.c_str(), NULL);
   int failures = 0;
   {
      std::unique_ptr<Mesh> smesh = cycle::LoadScaledMesh(cfg.shared.mesh.path, cfg.mesh_scale);
      ParMesh pmesh(mpi.GetComm(), *smesh);
      smesh.reset();
      if (mpi.IsRoot())
      {
         std::cout << "test_cycle_absorbing: " << deck << ", " << cfg.shared.mesh.path << " x "
                   << cfg.mesh_scale << ", p = " << cfg.shared.mesh.order << ", ranks "
                   << mpi.Size() << ", commit " << SEAS_CYCLE_COMMIT << "\n";
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

      // loaded_state: the slip's sign and size, then the QD balance at t_s.
      const real_t ts = 3.0e9;
      real_t Sp = 0.0;
      {
         Vector s0, psi0;
         LoadedFaultState(qd, pre, 0.0, 0.0, 0.0, 0.0, s0, psi0);
         qd.SetFaultState(s0, psi0);
         const real_t V0 = qd.Solve(0.0);
         if (mpi.IsRoot())
         {
            std::cout << std::scientific << std::setprecision(6) << "\nloaded_state: max V "
                      << V0 << " m/s at t = 0 (slip 0)\n";
         }
         int matches = 0;
         for (const real_t sg : {1.0, -1.0})
         {
            const real_t S = sg * cfg.plate_rate * ts;
            LoadedFaultState(qd, pre, S, 0.0, 0.0, 0.0, s0, psi0);
            qd.SetFaultState(s0, psi0);
            const real_t Vu = qd.Solve(ts), rel = std::abs(Vu / V0 - 1.0);
            if (mpi.IsRoot())
            {
               std::cout << "  uniform strike slip " << S << " m at t_s: max V " << Vu
                         << " m/s, |ratio - 1| " << rel << "\n";
            }
            if (rel <= 1e-6) { matches++; Sp = S; }
         }
         c.True("loaded_state: one sign of uniform slip is stress free", matches == 1,
                "|max V ratio - 1| <= 1e-6 for exactly one sign");
         if (matches != 1) { Sp = cfg.plate_rate * ts; }
      }
      Balance bal;
      {
         Vector s0, psi0;
         LoadedFaultState(qd, pre, Sp, 0.5, 0.1, 0.01, s0, psi0);
         qd.SetFaultState(s0, psi0);
         bal.Vmax = qd.Solve(ts);
         qd.Slip(bal.sk);
         qd.Psi(bal.psi);
         qd.SigmaN(bal.sn);
         bal.Vqd = qd.SlipRate();
         bal.u = qd.Displacement();
         bal.t = ts;
      }
      const real_t mus = c.MaxAbs(bal.u);   // collective: outside the root's print
      if (mpi.IsRoot())
      {
         std::cout << std::scientific << std::setprecision(3) << "\nQD balance at t_s = " << ts
                   << " s: walls at +-" << 0.5 * cfg.plate_rate * ts << " m, max V "
                   << bal.Vmax << " m/s, max |u_s| " << mus << " m\n";
      }
      const std::vector<int> qd_attrs = BdrAttributes(pmesh);

      // ---- reflecting (the reference) ----
      Run refl;
      real_t aK = 0.0;   // max |M^-1 K u_s|, the reflecting W
      {
         cycle::DynPhase dyn(pmesh, cfg_off, qd.Table(), par, pre);
         refl = Phase(dyn, bal, n_dyn);
         c.True("off: no frozen load", dyn.FrozenWallLoad().Size() == 0, "");
         Vector Ku, MKu;
         dyn.Op().StiffnessMult(bal.u, Ku);
         dyn.Op().Mass().Apply(Ku, MKu);
         aK = c.MaxAbs(MKu);
         if (mpi.IsRoot())
         {
            std::cout << std::scientific << std::setprecision(3)
                      << "  reflecting W: a_K = max |M^-1 K u_s| = " << aK << " m/s^2\n";
         }
      }

      // ---- "outer" ----
      if (mpi.IsRoot()) { std::cout << "\nretag, frozen_load, switch_equilibrium (\"outer\")\n"; }
      Run outer;
      {
         cycle::DynPhase dyn(pmesh, cfg_outer, qd.Table(), par, pre);
         const cycle::OuterFaceTags &t = dyn.AbsorbingFaces();
         if (mpi.IsRoot())
         {
            std::cout << "  re-tagged: " << t.natural << " traction-free faces -> attribute "
                      << t.natural_to << ", " << t.dirichlet << " wall faces -> attribute "
                      << t.dirichlet_to << "; " << t.free_surface << " free-surface faces kept\n";
         }
         const std::string mp = cfg.shared.mesh.path;
         if (mp.find("bp5_tandem_exact") != std::string::npos)
         {
            c.True("retag: BP5 mesh counts", t.natural == 138 && t.dirichlet == 240
                   && t.free_surface == 1876, "138 bottom, 240 wall, 1876 free-surface faces");
         }
         c.True("retag: the QD mesh is unchanged", BdrAttributes(pmesh) == qd_attrs,
                "every boundary attribute of the QD regime's mesh");
         const std::set<int> &abs = dyn.Boundaries().absorbing_attrs;
         c.True("retag: W absorbs on both new attributes",
                abs.count(t.natural_to) == 1 && abs.count(t.dirichlet_to) == 1
                && abs.size() == 2, "");

         outer = Phase(dyn, bal, 0);   // A^0 with b, no step yet

         // frozen_load
         NitscheWaveOperator &W = dyn.Op();
         const std::vector<char> wall = ElementDofMask(W, t.dirichlet_to);
         const Vector &b = dyn.FrozenWallLoad();
         const real_t bw = c.MaxAbsOn(b, wall, 1), bo = c.MaxAbsOn(b, wall, 0);
         const ParFiniteElementSpace &fes = W.Elasticity().GetFESpace();
         real_t fsum[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};   // sum b_c, sum |b_c|
         for (int i = 0; i < fes.GetNDofs(); i++)
         {
            for (int d = 0; d < 3; d++)
            {
               const real_t v = b(fes.DofToVDof(i, d));
               fsum[d] += v;
               fsum[3 + d] += std::abs(v);
            }
         }
         MPI_Allreduce(MPI_IN_PLACE, fsum, 6, MPITypeMap<real_t>::mpi_type, MPI_SUM,
                       mpi.GetComm());
         real_t net = 0.0;
         for (int d = 0; d < 3; d++)
         {
            net = std::max(net, fsum[3 + d] > 0.0 ? std::abs(fsum[d]) / fsum[3 + d] : 0.0);
         }
         if (mpi.IsRoot())
         {
            std::cout << std::scientific << std::setprecision(3)
                      << "  frozen load: max |b| on the wall elements " << bw << " N, elsewhere "
                      << bo << " N; net force (x, y, z) " << fsum[0] << ", " << fsum[1] << ", "
                      << fsum[2] << " N against sum |b_c| " << fsum[3] << ", " << fsum[4] << ", "
                      << fsum[5] << " N\n";
         }
         c.True("frozen_load: the walls carry a load", bw > 0.0, "max |b| on the wall elements > 0");
         c.Report("frozen_load: b off the wall elements", bo, 0.0, "max |b| (N)");
         c.Report("frozen_load: net force of b", net, 1e-6, "max_c |sum b_c| / sum |b_c|");

         // switch_equilibrium
         const real_t mb = c.MaxAbs(outer.Minv_b);
         c.Report("switch_equilibrium: A0 off the wall elements",
                  c.DiffInfOn(outer.A0, refl.A0, wall, 0) / aK, 1e-12,
                  "max |A0_outer - A0_refl| / a_K");
         c.Report("switch_equilibrium: A0 on the wall elements",
                  c.MaxAbsOn(outer.A0, wall, 1) / aK, 1e-12, "max |A0_outer| / a_K");
         const real_t rw = c.MaxAbsOn(refl.A0, wall, 1) / aK;
         Vector none, a_nob, zero(bal.u.Size());
         zero = 0.0;
         W.SetConstantLoad(none);
         W.Mult(bal.u, zero, a_nob);
         const real_t ctl = c.DiffInfOn(a_nob, refl.A0, wall, 1) / mb;
         if (mpi.IsRoot())
         {
            std::cout << std::scientific << std::setprecision(3)
                      << "  reflecting A0 on the wall elements (M^-1 r, left out of b): " << rw
                      << " a_K; max |M^-1 b| = " << mb / aK << " a_K\n"
                      << "  control, without b: max |A0 - A0_refl| on the wall elements = "
                      << ctl << " max |M^-1 b|\n";
         }
         c.True("control: without b the walls unload", ctl > 0.5,
                "max |A0 - A0_refl| / max |M^-1 b| > 0.5 on the wall elements");

         // first_steps: enter again (b set again), n steps
         outer = Phase(dyn, bal, n_dyn);
      }
      if (mpi.IsRoot()) { std::cout << "\nfirst_steps (" << n_dyn << " steps)\n"; }
      c.Report("first_steps outer: slip", c.RelInf(outer.slip, refl.slip), 1e-12, "max rel diff");
      c.Report("first_steps outer: slip rate", c.RelInf(outer.V, refl.V), 1e-12, "max rel diff");
      c.Report("first_steps outer: psi", c.RelInf(outer.psi, refl.psi), 1e-12, "max rel diff");
      {
         Vector du = refl.U;
         du -= bal.u;
         const real_t dU = c.DiffInf(outer.U, refl.U), mU = c.MaxAbs(du);
         if (mpi.IsRoot())
         {
            std::cout << std::scientific << std::setprecision(3) << "  U: max |U_outer - U_refl| "
                      << dU << " m against max |U_refl - u_s| " << mU
                      << " m (printed only: the reflecting walls move Vp tau / 2 = "
                      << 0.5 * cfg.plate_rate * refl.tau << " m in tau = " << refl.tau
                      << " s)\n";
         }
      }

      // ---- "bottom" ----
      if (mpi.IsRoot()) { std::cout << "\nbottom\n"; }
      {
         cycle::DynPhase dyn(pmesh, cfg_bottom, qd.Table(), par, pre);
         const cycle::OuterFaceTags &t = dyn.AbsorbingFaces();
         c.True("bottom: no wall faces, no frozen load", t.dirichlet == 0 && t.natural > 0,
                std::to_string(t.natural) + " traction-free faces absorb");
         const Run bot = Phase(dyn, bal, n_dyn);
         c.True("bottom: no frozen load after Enter", dyn.FrozenWallLoad().Size() == 0, "");
         c.Report("bottom: A0 vs reflecting", c.DiffInf(bot.A0, refl.A0) / aK, 1e-13,
                  "max |dA0| / a_K");
         c.Report("bottom: slip rate", c.RelInf(bot.V, refl.V), 1e-12, "max rel diff");
         c.Report("bottom: U", c.RelInf(bot.U, refl.U), 1e-12, "max rel diff");
      }

      failures = c.failures;
      if (mpi.IsRoot())
      {
         std::cout << "\ntest_cycle_absorbing: " << (failures ? "FAILED" : "PASSED") << " ("
                   << failures << " failing checks)\n";
      }
   }
   MFEMFinalizePetsc();
   return failures ? 1 : 0;
#endif
}
