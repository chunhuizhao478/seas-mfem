// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// tests/2d_inplane/test_fault_opening.cpp — [cycle.opening] mode =
// "shear_free", version A of the fault opening design rev 1.1 (Chunhui Zhao,
// 2026-11-04: "start version A", "build and test here";
// cycle/fault_opening.hpp).  Operator-level checks on two simple-case decks
// (their mesh and parameters; the Green's function and the opening are
// switched on here, no deck is run):
//
//   T1  static crack (D0, straight fault y = 0, |x| <= 15 km, box +-240 km):
//       a uniform shear pre-stress tau_inf, the open DOFs' slip from the QD
//       open-set solve against plane-strain mode II,
//         delta(x) = -2 (1 - nu) tau_inf sqrt(a^2 - x^2) / mu
//       (the kernel's slip is anti-parallel to tau).  a) the whole fault
//       open (a = 15 km, the fault tips): error at |x| <= 0.8 a below 3% of
//       delta(0); b) a patch |x| <= 7.5 km open in a closed, unslipped fault:
//       delta(0) between the values for a = 7.5 km and a + h.  Both: shear
//       at the open DOFs 0 to 1e-9 tau_inf.
//   T2  open and re-close (D9f arc, the 2.5 km bend of 0.8 km offset): a
//       tapered slip A w(x) on the closed fault, A = 1.5 x the amplitude
//       that first makes sigma_n' tensile.  a) the settle opens DOFs and ends
//       consistent: every open DOF sigma < 0, every closed sigma >= 0; shear
//       0 at the open DOFs; b) re-closing them all at the same slip leaves
//       their shear 0 (no traction jump at re-closure); c) the reversed slip
//       -A w closes DOFs of a) and opens others; d) a short QD run from a):
//       psi of the open DOFs unchanged bit for bit, their V = 0, max V = the
//       closed DOFs' max, shear 0 at the open DOFs after every accepted step;
//       e) restart: from a copy of (state, flags) at step k the next steps
//       repeat bit for bit.
//   T4  switch (D9f, from T2 a): dynamic steps with the QD open flags handed
//       over: at every DOF open in a step the end-of-step shear X + gamma V
//       is 0, max V of the record is the closed DOFs' max, psi of DOFs open
//       in every step unchanged; back to QD (Reenter with the dynamic flags):
//       shear 0 at the QD open DOFs, so the jump there is 0; the QD max V
//       (the guard's) is the closed DOFs' max.
//   C1  checkpoint V6 round trip: the open flags, with and without the
//       relaxation record.
// Refusals (--refuse n, each must abort with its message; the deck is D0
// with a [cycle.opening] section appended, written to
// test_fault_opening_refuse.toml, which --clean deletes): 1 mode "contact"
// with D0's [cycle.qd.green].check_every (the solve-path check carries no
// gap), 2 an unknown mode, 3 max_passes 0, 4 no [cycle.qd.green], 5 a 3-D
// mesh (the deck: config/fd_cycle/bp5_scec_fd.toml with
// [cycle.qd].operators = "matrix" and [cycle.qd.green] on), 6
// [friction.prakash_clifton], 7 [cycle.probe].  Contact (version B, plan
// rev 1): 8 the H-matrix path (compress), 9 a contact knob (floor_taper)
// under mode "shear_free", 10 floor_taper without [cycle].sigma_n_floor_pa,
// 11 sigma_tol_pa < 0.
//
// Deck mode: --load-deck d.toml only reads a deck (CycleConfig::Load) and
// prints its [cycle.opening] (the committed opening decks, no solve).
//
// Usage (from miniapps/seas; PETSc needed): mpirun -np N ./seas_test_fault_opening
//        [--only T1,T2,T4,C1] [--refuse n] [--clean] [--load-deck d.toml]

#include "mfem.hpp"
#include "../../cycle/scaled_mesh.hpp"
#include "../../cycle/cycle_config.hpp"
#include "../../cycle/fault_fields.hpp"
#include "../../cycle/prestress.hpp"
#include "../../cycle/qd_phase.hpp"
#include "../../cycle/dyn_phase.hpp"
#include "../../cycle/state_bridge.hpp"
#include "../../cycle/checkpoint.hpp"
#include "../../common/mpi_context.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace
{
int fails = 0;
void Check(MPIContext &mpi, bool ok, const std::string &what, real_t v)
{
   if (mpi.IsRoot()) { std::printf("  [%s] %s  (%.4e)\n", ok ? "PASS" : "FAIL", what.c_str(), v); }
   if (!ok) { fails++; }
}
real_t GMax(MPIContext &mpi, real_t v)
{
   MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, mpi.GetComm());
   return v;
}
real_t GMin(MPIContext &mpi, real_t v) { return -GMax(mpi, -v); }
long GSum(MPIContext &mpi, long v)
{
   MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_LONG, MPI_SUM, mpi.GetComm());
   return v;
}
bool Wanted(const std::string &only, const char *name)
{
   return only.empty() || ("," + only + ",").find(std::string(",") + name + ",") != std::string::npos;
}

/// The QD stack of a deck with the Green's function and the opening on.
struct Stack
{
   cycle::CycleConfig cfg;
   std::unique_ptr<ParMesh> pmesh;
   std::unique_ptr<cycle::QDPhase> qd;
   cycle::FaultParams par;
   cycle::Prestress pre;

   Stack(MPIContext &mpi, const std::string &deck, bool opening,
         const std::function<void(cycle::Prestress &, const cycle::FaultDOFTable &)> &edit = nullptr)
      : cfg(cycle::CycleConfig::Load(deck))
   {
      cfg.qd.operators = "matrix";
      cfg.qd.green.enabled = true;
      cfg.qd.green.cache.clear();
      cfg.qd.green.check_every = 0;
      if (opening) { cfg.opening.mode = "shear_free"; }
      std::unique_ptr<Mesh> smesh = cycle::LoadScaledMesh(cfg.shared.mesh.path, cfg.mesh_scale);
      pmesh = std::make_unique<ParMesh>(mpi.GetComm(), *smesh);
      smesh.reset();
      cfg.CheckDimension(pmesh->Dimension());
      qd = std::make_unique<cycle::QDPhase>(*pmesh, mpi, cfg);
      const spatial::SpatialFrictionConfig &s = cfg.shared;
      const MaterialField material = MaterialField::MakeConstant(
         s.material_fallback.lambda, s.material_fallback.mu, s.material_fallback.rho);
      const spatial::PorePressureSpec no_pp;
      par = cycle::BuildFaultParams(qd->Table(), *s.rate_state, cfg.fields, material, *pmesh,
                                    no_pp);
      cycle::PrestressConfig pc = cfg.prestress;
      pc.a_max = cycle::GlobalMaxA(par, mpi.GetComm());
      pre = cycle::BuildPrestress(par, pc);
      if (edit) { edit(pre, qd->Table()); }
      qd->Setup(par, pre);
   }
   int N() const { return qd->Table().Size(); }
   real_t X(int j) const { return qd->Table().coords(3 * j); }
};

/// Largest |tau_pre + T| (in-plane slot) over the DOFs with flag 1.
real_t ShearAt(MPIContext &mpi, cycle::QDPhase &qd, const std::vector<char> &which)
{
   const Vector &T = qd.Traction();
   const Vector &tp = qd.Fault().TauPre();
   real_t m = 0.0;
   for (size_t j = 0; j < which.size(); j++)
   {
      if (!which[j]) { continue; }
      m = std::max(m, std::hypot(tp(2 * j) + T(2 * j), tp(2 * j + 1) + T(2 * j + 1)));
   }
   return GMax(mpi, m);
}

// ------------------------------------------------------------------ T1
void TestT1(MPIContext &mpi)
{
   if (mpi.IsRoot()) { std::cout << "\nT1: static mode II crack (D0)\n"; }
   const real_t tau_inf = 5.0e6;
   int slot = -1;
   Stack st(mpi, "simple_cases/configs/D0_straight.toml", true,
            [&](cycle::Prestress &pre, const cycle::FaultDOFTable &)
   {
      // uniform shear in the in-plane slot (map view: strike, 1; checked
      // below), none out of the plane
      for (int i = 0; i < pre.tau0.Size() / 2; i++)
      {
         pre.tau0(2 * i) = 0.0;
         pre.tau0(2 * i + 1) = tau_inf;
      }
   });
   cycle::QDPhase &qd = *st.qd;
   slot = qd.Opening()->Slot();
   MFEM_VERIFY(slot == 1, "T1: the planar slot is " << slot << ", the test set tau_pre in 1");
   const int N = st.N();
   const auto &m = st.cfg.shared.material_fallback;
   const real_t nu = m.lambda / (2.0 * (m.lambda + m.mu));
   const real_t C = 2.0 * (1.0 - nu) / m.mu;
   real_t xlo = 1e30, xhi = -1e30;
   for (int j = 0; j < N; j++) { xlo = std::min(xlo, st.X(j)); xhi = std::max(xhi, st.X(j)); }
   xlo = GMin(mpi, xlo);
   xhi = GMax(mpi, xhi);
   const real_t a_f = 0.5 * (xhi - xlo), xc = 0.5 * (xhi + xlo);
   real_t tp_dev = 0.0;
   for (int j = 0; j < N; j++)
   {
      tp_dev = std::max(tp_dev, std::abs(qd.Fault().TauPre()(2 * j + slot) - tau_inf));
   }
   tp_dev = GMax(mpi, tp_dev);
   Check(mpi, tp_dev == 0.0, "the kernel's tau_pre is the uniform tau_inf", tp_dev);
   if (mpi.IsRoot())
   {
      std::printf("  fault x in [%.1f, %.1f] km, %d global DOFs, nu %.3f, slot %d\n", xlo / 1e3,
                  xhi / 1e3, qd.NumGlobalFaultDOFs(), nu, slot);
   }
   Vector zero_slip(2 * N), psi;
   zero_slip = 0.0;
   qd.Psi(psi);

   auto run = [&](real_t a_open, const char *name, bool whole)
   {
      std::vector<char> f(N, 0);
      for (int j = 0; j < N; j++) { f[j] = (std::abs(st.X(j) - xc) <= a_open + 1e-6) ? 1 : 0; }
      qd.SetFaultState(zero_slip, psi);
      qd.SetOpenFlags(f);
      qd.Solve(0.0);
      const Vector &s = qd.Operator().GetStageSlip();
      const real_t res = ShearAt(mpi, qd, f);
      Check(mpi, res <= 1e-9 * tau_inf, std::string(name) + ": shear at the open DOFs <= 1e-9 tau_inf",
            res / tau_inf);
      // centre slip and the profile error
      const real_t d0 = C * tau_inf * a_open;
      real_t err = 0.0, sgn_bad = 0.0, center = 0.0, cdist = 1e30;
      for (int j = 0; j < N; j++)
      {
         if (!f[j]) { continue; }
         const real_t x = st.X(j) - xc;
         const real_t num = s(2 * j + slot);
         if (num * tau_inf > 0.0) { sgn_bad = std::max(sgn_bad, std::abs(num)); }
         if (std::abs(x) < cdist) { cdist = std::abs(x); center = num; }
         if (std::abs(x) <= 0.8 * a_open)
         {
            const real_t ex = -C * tau_inf * std::sqrt(a_open * a_open - x * x);
            err = std::max(err, std::abs(num - ex));
         }
      }
      err = GMax(mpi, err);
      sgn_bad = GMax(mpi, sgn_bad);
      // the centre value from the rank holding the DOF nearest x = 0
      struct { double d; int r; } in{static_cast<double>(cdist), mpi.Rank()}, out;
      MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, mpi.GetComm());
      double cv = center;
      MPI_Bcast(&cv, 1, MPI_DOUBLE, out.r, mpi.GetComm());
      Check(mpi, sgn_bad == 0.0, std::string(name) + ": open slip anti-parallel to tau_pre "
            "(the kernel's convention)", sgn_bad);
      if (mpi.IsRoot())
      {
         std::printf("  %s: %d open DOFs, delta(0) = %.6e m (mode II, a = %.2f km: %.6e m; "
                     "a + h: %.6e m), max error at |x| <= 0.8 a %.4e of delta(0)\n", name,
                     qd.Opening()->NumOpenGlobal(), cv, a_open / 1e3, -d0,
                     -C * tau_inf * (a_open + 500.0), err / d0);
      }
      if (whole)
      {
         Check(mpi, err <= 0.03 * d0, std::string(name) + ": |delta - mode II| at |x| <= 0.8 a "
               "<= 3% of delta(0)", err / d0);
      }
      else
      {
         const real_t lo = C * tau_inf * a_open, hi = C * tau_inf * (a_open + 500.0);
         Check(mpi, -cv >= lo && -cv <= hi, std::string(name) + ": |delta(0)| between mode II "
               "for a and a + h", -cv / lo);
      }
      // the normal traction of the open slip on a straight fault (the
      // continuum: 0 by antisymmetry; the mesh is not symmetric, and the
      // crack-tip field is singular): over all DOFs and over |x| <= 0.8 a
      real_t dn = 0.0, dn_in = 0.0, x_dn = 0.0;
      const Vector &Tn = qd.NormalTraction();
      for (int j = 0; j < N; j++)
      {
         const real_t x = st.X(j) - xc;
         if (std::abs(Tn(j)) > dn) { dn = std::abs(Tn(j)); x_dn = x; }
         if (std::abs(x) <= 0.8 * a_open) { dn_in = std::max(dn_in, std::abs(Tn(j))); }
      }
      struct { double v; int r; } ti{static_cast<double>(dn), mpi.Rank()}, to;
      MPI_Allreduce(&ti, &to, 1, MPI_DOUBLE_INT, MPI_MAXLOC, mpi.GetComm());
      double xd = x_dn;
      MPI_Bcast(&xd, 1, MPI_DOUBLE, to.r, mpi.GetComm());
      dn = to.v;
      dn_in = GMax(mpi, dn_in);
      if (mpi.IsRoot())
      {
         std::printf("  %s: |T_n| of the open slip on the straight fault: max %.3e tau_inf at "
                     "x = %.2f km, max over |x| <= 0.8 a %.3e tau_inf (information)\n", name,
                     dn / tau_inf, xd / 1e3, dn_in / tau_inf);
      }
   };
   run(a_f, "T1a (whole fault open)", true);
   run(7.5e3, "T1b (patch |x| <= 7.5 km open)", false);
   if (mpi.IsRoot())
   {
      std::printf("  unit products (Green's-function columns fetched) %ld\n",
                  qd.Opening()->UnitProducts());
   }
}

// ------------------------------------------------------------------ T2, T4
void TestT2T4(MPIContext &mpi, bool do_t2, bool do_t4)
{
   if (mpi.IsRoot()) { std::cout << "\nT2/T4: open, re-close and switch (D9f arc)\n"; }
   Stack st(mpi, "simple_cases/configs/D9f_arc08L2p5_h250_avg_floor10.toml", true);
   cycle::QDPhase &qd = *st.qd;
   const int N = st.N();
   const int slot = qd.Opening()->Slot();
   real_t xlo = 1e30, xhi = -1e30;
   for (int j = 0; j < N; j++) { xlo = std::min(xlo, st.X(j)); xhi = std::max(xhi, st.X(j)); }
   xlo = GMin(mpi, xlo);
   xhi = GMax(mpi, xhi);
   const real_t xc = 0.5 * (xlo + xhi), L = 0.45 * (xhi - xlo);
   Vector w(2 * N), psi0;
   w = 0.0;
   for (int j = 0; j < N; j++)
   {
      const real_t u = (st.X(j) - xc) / L;
      w(2 * j + slot) = std::max(0.0, 1.0 - u * u);
   }
   qd.Psi(psi0);
   const std::vector<char> none(N, 0);

   // sigma at zero slip and at slip w (all closed): the amplitude of first tension
   Vector zero(2 * N), sn0, sn1;
   zero = 0.0;
   qd.SetOpenFlags(none);
   qd.SetFaultState(zero, psi0);
   qd.Solve(0.0);
   qd.SigmaN(sn0);
   qd.SetFaultState(w, psi0);
   qd.Solve(0.0);
   qd.SigmaN(sn1);
   // by linearity sigma(A) = sn0 + A (sn1 - sn0): the amplitude of the
   // first tension for +A w (Ap) and for -A w (Am); A opens both sides
   real_t Ap = 1e30, Am = 1e30;
   for (int j = 0; j < N; j++)
   {
      const real_t d = sn1(j) - sn0(j);
      if (d < 0.0) { Ap = std::min(Ap, sn0(j) / -d); }
      if (d > 0.0) { Am = std::min(Am, sn0(j) / d); }
   }
   Ap = GMin(mpi, Ap);
   Am = GMin(mpi, Am);
   MFEM_VERIFY(Ap < 1e29 && Am < 1e29, "T2: the slip w makes no tension on one side "
               "(Ap " << Ap << ", Am " << Am << ")");
   const real_t A = 1.5 * std::max(Ap, Am);
   if (mpi.IsRoot())
   {
      std::printf("  fault x in [%.2f, %.2f] km, %d global DOFs, slot %d; taper half-length "
                  "%.2f km; first tension at A = %.4f m (+w), %.4f m (-w); A = %.4f m\n",
                  xlo / 1e3, xhi / 1e3, qd.NumGlobalFaultDOFs(), slot, L / 1e3, Ap, Am, A);
   }
   Vector sA(w);
   sA *= A;

   // a) settle at +A w
   qd.SetFaultState(sA, psi0);
   const long flips_a = qd.SettleOpening(0.0);
   const std::vector<char> Oa = qd.OpenFlags();
   const long na = GSum(mpi, std::count(Oa.begin(), Oa.end(), char(1)));
   Vector sn;
   qd.SigmaN(sn);
   real_t open_max = -1e30, closed_min = 1e30;
   for (int j = 0; j < N; j++)
   {
      if (Oa[j]) { open_max = std::max(open_max, sn(j)); }
      else { closed_min = std::min(closed_min, sn(j)); }
   }
   open_max = GMax(mpi, open_max);
   closed_min = GMin(mpi, closed_min);
   const real_t tau_scale = 1.0e6;
   const real_t res_a = ShearAt(mpi, qd, Oa);
   if (mpi.IsRoot())
   {
      std::printf("  a) +A w: %ld open DOFs (%ld flips), max sigma open %.4f MPa, min sigma "
                  "closed %.4f MPa, shear at open %.3e Pa\n", na, flips_a, open_max / 1e6,
                  closed_min / 1e6, res_a);
   }
   if (do_t2)
   {
      Check(mpi, na > 0, "T2a: +A w opens DOFs", static_cast<real_t>(na));
      Check(mpi, open_max < 0.0 && closed_min >= 0.0, "T2a: the settled set is consistent "
            "(open sigma < 0, closed sigma >= 0)", open_max);
      Check(mpi, res_a <= 1e-9 * tau_scale, "T2a: shear at the open DOFs <= 1e-3 Pa", res_a);
   }
   // the settled state, for d), e) and T4
   Vector state_a(qd.State());

   if (do_t2)
   {
      // b) re-close all at the same slip
      qd.SetOpenFlags(none);
      qd.Solve(0.0);
      const real_t res_b = ShearAt(mpi, qd, Oa);
      Check(mpi, res_b <= 1e-9 * tau_scale, "T2b: re-closed at the same slip: shear at those "
            "DOFs stays 0 (no jump)", res_b);

      // c) the reversed slip
      Vector sm(w);
      sm *= -A;
      qd.SetFaultState(sm, psi0);
      qd.SettleOpening(0.0);
      const std::vector<char> Oc = qd.OpenFlags();
      long closed_again = 0, newly = 0;
      for (int j = 0; j < N; j++)
      {
         if (Oa[j] && !Oc[j]) { closed_again++; }
         if (!Oa[j] && Oc[j]) { newly++; }
      }
      closed_again = GSum(mpi, closed_again);
      newly = GSum(mpi, newly);
      if (mpi.IsRoot())
      {
         std::printf("  c) -A w: %ld of a)'s open DOFs closed, %ld others opened\n",
                     closed_again, newly);
      }
      Check(mpi, closed_again == na, "T2c: -A w re-closes every DOF a) opened",
            static_cast<real_t>(closed_again));
      Check(mpi, newly > 0, "T2c: -A w opens the other side", static_cast<real_t>(newly));

      // d) a short QD run from a)
      qd.State() = state_a;
      qd.SetOpenFlags(Oa);
      qd.Solve(0.0);
      qd.InitTS(1.0e-3);
      Vector psi_start;
      qd.Psi(psi_start);
      real_t worst_res = 0.0, psi_dev = 0.0, vopen = 0.0, vmax_err = 0.0;
      std::vector<char> still_open(Oa);   // open in every accepted step so far
      long flips_run = 0;
      const long flips0 = qd.Opening()->Flips();
      int steps = 0;
      std::vector<Vector> rec_a;
      const int k_ckpt = 6, n_more = 6;
      Vector state_k, u_k, V_k;
      std::vector<char> flags_k;
      real_t t_k = 0.0, dt_k = 0.0;
      int step_k = 0;
      long flips_prev = flips0, flip_steps = 0, dt_bad = 0;
      auto mon = [&](int step, real_t tt) -> bool
      {
         const std::vector<char> &O = qd.OpenFlags();
         if (qd.Opening()->Flips() != flips_prev)
         {
            // a step with flips: the next one starts from [time].dt_initial
            flips_prev = qd.Opening()->Flips();
            flip_steps++;
            if (qd.NextDt() != st.cfg.dt_init) { dt_bad++; }
         }
         worst_res = std::max(worst_res, ShearAt(mpi, qd, O));
         Vector p, vabs;
         qd.Psi(p);
         qd.SlipRateMagnitude(vabs);
         real_t vc = 0.0, vo = 0.0, dp = 0.0;
         for (int j = 0; j < N; j++)
         {
            if (!O[j]) { still_open[j] = 0; }
            if (O[j])
            {
               vo = std::max(vo, vabs(j));
               if (still_open[j]) { dp = std::max(dp, std::abs(p(j) - psi_start(j))); }
            }
            else { vc = std::max(vc, vabs(j)); }
         }
         vc = GMax(mpi, vc);
         vopen = std::max(vopen, GMax(mpi, vo));
         psi_dev = std::max(psi_dev, GMax(mpi, dp));
         vmax_err = std::max(vmax_err, std::abs(qd.MaxV() - vc) / vc);
         steps = step;
         if (step == k_ckpt)
         {
            state_k = qd.State();
            u_k = qd.Displacement();
            V_k = qd.SlipRate();
            flags_k = O;
            t_k = tt;
            dt_k = qd.NextDt();
            step_k = step;
            qd.RestartStage();   // as a checkpoint step does
         }
         if (step > k_ckpt) { rec_a.push_back(qd.State()); }
         return true;
      };
      real_t t = 0.0, dt = 1.0e-3;
      qd.Run(t, dt, 1.0e9, k_ckpt + n_more, mon);
      flips_run = qd.Opening()->Flips() - flips0;
      long n_still = 0;
      for (int j = 0; j < N; j++) { n_still += still_open[j]; }
      n_still = GSum(mpi, n_still);
      if (mpi.IsRoot())
      {
         std::printf("  d) QD run: %d steps to t = %.4e s, %ld flips, %ld DOFs open in every "
                     "step; max shear at open %.3e Pa, max |psi - psi0| there %.3e, max V at "
                     "open %.3e, |MaxV - closed max| / closed max %.3e\n", steps, t, flips_run,
                     n_still, worst_res, psi_dev, vopen, vmax_err);
      }
      Check(mpi, steps == k_ckpt + n_more, "T2d: the QD run took its steps", steps);
      Check(mpi, worst_res <= 1e-9 * tau_scale, "T2d: shear at the open DOFs 0 after every "
            "accepted step", worst_res);
      Check(mpi, psi_dev == 0.0, "T2d: psi frozen at open DOFs (bit for bit)", psi_dev);
      Check(mpi, vopen == 0.0, "T2d: V = 0 at open DOFs", vopen);
      Check(mpi, n_still > 0, "T2d: DOFs open in every step", static_cast<real_t>(n_still));
      Check(mpi, vmax_err <= 1e-14, "T2d: max V = the closed DOFs' max (to rounding)",
            vmax_err);
      if (mpi.IsRoot())
      {
         std::printf("  d) %ld steps with flips; next dt after each: [time].dt_initial = %.3e s "
                     "(%ld not)\n", flip_steps, st.cfg.dt_init, dt_bad);
      }
      Check(mpi, flip_steps > 0 && dt_bad == 0, "T2d: after a step with flips the next dt is "
            "[time].dt_initial", static_cast<real_t>(dt_bad));

      // e) restart from step k: the next steps repeat bit for bit
      std::vector<Vector> rec_b;
      qd.State() = state_k;
      qd.Operator().SetDisplacement(u_k);   // as the driver's restart
      qd.Fault().SetSlipRate(V_k);
      qd.SetOpenFlags(flags_k);
      qd.RestoreTS(step_k);
      auto mon_b = [&](int step, real_t) -> bool
      {
         if (step > k_ckpt) { rec_b.push_back(qd.State()); }
         return true;
      };
      t = t_k;
      dt = dt_k;
      qd.Run(t, dt, 1.0e9, k_ckpt + n_more, mon_b);
      bool same = rec_a.size() == rec_b.size() && !rec_a.empty();
      for (size_t i = 0; same && i < rec_a.size(); i++)
      {
         same = (rec_a[i].Size() == rec_b[i].Size()) &&
                std::memcmp(rec_a[i].GetData(), rec_b[i].GetData(),
                            sizeof(real_t) * rec_a[i].Size()) == 0;
      }
      int all_same = same ? 1 : 0;
      MPI_Allreduce(MPI_IN_PLACE, &all_same, 1, MPI_INT, MPI_MIN, mpi.GetComm());
      Check(mpi, all_same == 1, "T2e: restart from (state, flags) at step k repeats the next "
            "steps bit for bit", static_cast<real_t>(rec_b.size()));
   }

   if (!do_t4) { return; }
   // ---- T4: the switch from a) ----
   qd.State() = state_a;
   qd.SetOpenFlags(Oa);
   const real_t ts = 0.0;
   qd.Solve(ts);
   Vector skq, psiq, snq;
   qd.Slip(skq);
   qd.Psi(psiq);
   qd.SigmaN(snq);
   const Vector Vqd = qd.SlipRate();
   const Vector us = qd.Displacement();
   cycle::DynPhase dyn(*st.pmesh, st.cfg, qd.Table(), st.par, st.pre, 1.0, nullptr,
                       qd.Domain().FrameField(), qd.SubstrateShape(), qd.PmlBandPtr(),
                       qd.PmlWallSpeedFactor());
   Vector d, V, psin;
   cycle::bridge::ToNitsche(dyn.Frames(), skq, Vqd, psiq, d, V, psin);
   dyn.Enter(us, d, V, psin, snq, ts);
   dyn.SetOpenFlags(Oa);
   Vector psi_enter(dyn.Psi());
   std::vector<char> open_all(N, 1);
   real_t worst_tau = 0.0, vmax_err = 0.0, vopen_dyn = 0.0;
   long open_steps = 0;
   const int n_dyn = 60;
   cycle::StationFields f;
   for (int k = 0; k < n_dyn; k++)
   {
      dyn.Step();
      const std::vector<char> &O = dyn.OpenFlags();
      dyn.FillStationFields(f);
      Vector vabs;
      dyn.SlipRateMagnitude(vabs);
      real_t tm = 0.0, vc = 0.0;
      long no = 0;
      for (int j = 0; j < N; j++)
      {
         if (O[j])
         {
            no++;
            tm = std::max(tm, std::hypot(f.tau_a_dip(j), f.tau_a_strike(j)));
            vopen_dyn = std::max(vopen_dyn, vabs(j));
         }
         else
         {
            vc = std::max(vc, vabs(j));
            open_all[j] = 0;
         }
      }
      worst_tau = std::max(worst_tau, GMax(mpi, tm));
      const real_t vcg = GMax(mpi, vc);
      vmax_err = std::max(vmax_err, std::abs(dyn.LastStepInfo().max_V - vcg) / vcg);
      open_steps += GSum(mpi, no);
   }
   real_t psi_dev = 0.0;
   long n_open_all = 0;
   for (int j = 0; j < N; j++)
   {
      if (!open_all[j]) { continue; }
      n_open_all++;
      psi_dev = std::max(psi_dev, std::abs(dyn.Psi()(j) - psi_enter(j)));
   }
   psi_dev = GMax(mpi, psi_dev);
   n_open_all = GSum(mpi, n_open_all);
   vopen_dyn = GMax(mpi, vopen_dyn);
   const long n_open_end = dyn.OpenCount();
   if (mpi.IsRoot())
   {
      std::printf("  T4 dynamic: %d steps (tau = %.4e s), open DOF-steps %ld, open at the end "
                  "%ld, open in every step %ld; max |X + gamma V| at open %.3e Pa, |max V - "
                  "closed max| / closed max %.3e, psi change there %.3e; max V at open DOFs "
                  "%.4e m/s (left out of max V)\n", n_dyn, dyn.Tau(), open_steps,
                  n_open_end, n_open_all, worst_tau, vmax_err, psi_dev, vopen_dyn);
   }
   Check(mpi, open_steps > 0, "T4: open DOFs in the dynamic steps", static_cast<real_t>(open_steps));
   Check(mpi, worst_tau <= 1e-9 * tau_scale, "T4: end-of-step shear X + gamma V = 0 at open DOFs",
         worst_tau);
   Check(mpi, vmax_err <= 1e-14, "T4: the record's max V = the closed DOFs' max (to rounding)",
         vmax_err);
   Check(mpi, psi_dev == 0.0, "T4: psi unchanged at DOFs open in every step", psi_dev);

   // back to QD with the dynamic flags
   Vector sk2, psi2;
   cycle::bridge::ToQD(dyn.Frames(), dyn.Slip(), dyn.Psi(), sk2, psi2);
   std::vector<char> Od = dyn.OpenFlags();
   dyn.FillStationFields(f);
   const real_t Vg = qd.Reenter(sk2, psi2, dyn.Time(), 100, 1.0e-3, &Od);
   const std::vector<char> &Oq = qd.OpenFlags();
   const real_t res_q = ShearAt(mpi, qd, Oq);
   real_t jump = 0.0, vc = 0.0;
   long both = 0;
   cycle::StationFields fq;
   qd.FillStationFields(fq);
   Vector vq;
   qd.SlipRateMagnitude(vq);
   for (int j = 0; j < N; j++)
   {
      if (!Oq[j]) { vc = std::max(vc, vq(j)); }
      if (Od[j] && Oq[j])
      {
         both++;
         const real_t td = std::hypot(f.tau_a_dip(j), f.tau_a_strike(j));
         const real_t tq = std::hypot(fq.tau_a_dip(j) + fq.tau_b_dip(j),
                                      fq.tau_a_strike(j) + fq.tau_b_strike(j));
         jump = std::max(jump, std::abs(tq - td));
      }
   }
   jump = GMax(mpi, jump);
   vc = GMax(mpi, vc);
   both = GSum(mpi, both);
   if (mpi.IsRoot())
   {
      std::printf("  T4 back to QD: %d open (%ld flips in the settle), %ld open in both; shear "
                  "at QD open %.3e Pa, max jump there %.3e Pa; guard max V %.4e (closed max "
                  "%.4e)\n", qd.Opening()->NumOpenGlobal(), qd.ReenterFlips(), both, res_q,
                  jump, Vg, vc);
   }
   Check(mpi, res_q <= 1e-9 * tau_scale, "T4: shear 0 at the QD open DOFs after the switch back",
         res_q);
   Check(mpi, both > 0 && jump <= 1e-9 * tau_scale, "T4: the hand-back jump at DOFs open in both "
         "regimes is 0", jump);
   Check(mpi, std::abs(Vg - vc) <= 1e-14 * vc, "T4: the guard's max V is the closed DOFs' max",
         std::abs(Vg - vc) / vc);
}

// ------------------------------------------------------------------ C1
void TestC1(MPIContext &mpi)
{
   if (mpi.IsRoot()) { std::cout << "\nC1: checkpoint V6 round trip\n"; }
   for (int relax = 0; relax <= 1; relax++)
   {
      cycle::CheckpointInfo w;
      w.t = 1.25e9;
      w.dt_next = 3.5;
      w.step = 42;
      w.commit = "opening-test";
      w.ic = "scec_fd";
      w.fault_digest = 1.0;
      const int n = 7 + mpi.Rank();
      w.opening = true;
      w.open.resize(n);
      for (int i = 0; i < n; i++) { w.open[i] = (i % 3 == 1) ? 1 : 0; }
      w.relax = (relax == 1);
      if (w.relax)
      {
         w.relax_t_relax = 1.5e11;
         w.relax_t_n = w.t;
         w.relax_r.SetSize(n);
         w.relax_d.SetSize(n);
         for (int i = 0; i < n; i++) { w.relax_r(i) = 1e5 * i; w.relax_d(i) = -2e5 * i; }
      }
      cycle::CheckpointData d;
      d.state.SetSize(3 * n);
      d.state = 0.5;
      d.displacement.SetSize(4);
      d.displacement = 0.0;
      d.traction.SetSize(2 * n);
      d.traction = 1.0;
      d.slip_rate.SetSize(2 * n);
      d.slip_rate = 1e-9;
      const std::string prefix = "test_fault_opening_ckpt" + std::to_string(relax);
      cycle::WriteCycleCheckpoint(prefix, w, d, mpi);
      cycle::CheckpointInfo r;
      cycle::CheckpointData e;
      cycle::ReadCycleCheckpoint(prefix, "opening-test", false, r, e, mpi);
      bool ok = r.opening && r.open == w.open && r.relax == w.relax;
      if (w.relax)
      {
         for (int i = 0; ok && i < n; i++)
         {
            ok = r.relax_r(i) == w.relax_r(i) && r.relax_d(i) == w.relax_d(i);
         }
      }
      int all = ok ? 1 : 0;
      MPI_Allreduce(MPI_IN_PLACE, &all, 1, MPI_INT, MPI_MIN, mpi.GetComm());
      Check(mpi, all == 1, std::string("C1: V6 flags read back") + (relax ? " (with the "
            "relaxation record)" : " (no relaxation)"), 0.0);
      std::remove(CheckpointFilename(prefix, mpi.Rank()).c_str());
   }
}

/// Refusal n: D0 (5: BP5) with [cycle.opening] (and what n needs), read: it
/// must abort.
int Refuse(int n)
{
   const char *base = (n == 5) ? "config/fd_cycle/bp5_scec_fd.toml"
                               : "simple_cases/configs/D0_straight.toml";
   std::ifstream in(base);
   MFEM_VERIFY(in.good(), "cannot read " << base);
   std::stringstream ss;
   ss << in.rdbuf();
   std::string deck = ss.str();
   std::string mode = "shear_free", extra, open_extra;
   int passes = 20;
   auto erase_green = [&deck]()
   {
      const std::string g = "[cycle.qd.green]\ncheck_every = 10000\n";
      const size_t at = deck.find(g);
      MFEM_VERIFY(at != std::string::npos, "D0 has no '[cycle.qd.green] check_every' block");
      deck.erase(at, g.size());
   };
   // 1: D0 keeps check_every = 10000, which contact refuses (the solve-path
   // check carries no gap)
   if (n == 1) { mode = "contact"; }
   if (n == 2) { mode = "gap"; }
   if (n == 3) { passes = 0; }
   if (n == 4) { erase_green(); }
   if (n == 8)
   {
      mode = "contact";
      erase_green();
      extra = "\n[cycle.qd.green]\ncompress = 0.5\n";
   }
   if (n == 9) { open_extra = "floor_taper = true\n"; }
   if (n == 10)
   {
      // floor_taper defaults true and D0 has no [cycle].sigma_n_floor_pa
      mode = "contact";
      erase_green();
      extra = "\n[cycle.qd.green]\nblock = 8\n";
   }
   if (n == 11)
   {
      mode = "contact";
      erase_green();
      extra = "\n[cycle.qd.green]\nblock = 8\n";
      const std::string cy = "[cycle]\n";
      const size_t at = deck.find(cy);
      MFEM_VERIFY(at != std::string::npos, "D0 has no [cycle] section");
      deck.insert(at + cy.size(), "sigma_n_floor_pa = 10.0e6\n");
      open_extra = "sigma_tol_pa = -1.0\n";
   }
   if (n == 5)
   {
      const std::string q = "[cycle.qd]\n";
      const size_t at = deck.find(q);
      MFEM_VERIFY(at != std::string::npos && deck.find("operators", at) > deck.find("\n[", at + 1),
                  base << ": no [cycle.qd] section, or it has operators already");
      deck.insert(at + q.size(), "operators = \"matrix\"\n");
      extra = "\n[cycle.qd.green]\nenabled = true\n";
   }
   if (n == 6) { extra = "\n[friction.prakash_clifton]\nv_star = 1.0\nlength = 100.0\n"; }
   if (n == 7) { extra = "\n[cycle.probe]\nV_off = [1.0e-3]\nT_dwell = [10.0]\n"; }
   deck += "\n[cycle.opening]\nmode = \"" + mode + "\"\nmax_passes = " +
           std::to_string(passes) + "\n" + open_extra + extra;
   const std::string fn = "test_fault_opening_refuse.toml";
   {
      std::ofstream out(fn);
      out << deck;
   }
   cycle::CycleConfig c = cycle::CycleConfig::Load(fn);
   if (n == 5) { c.CheckDimension(3); }
   std::remove(fn.c_str());
   std::cerr << "refusal " << n << " did not abort\n";
   return 3;
}

}  // namespace

int main(int argc, char *argv[])
{
   MPIContext mpi(&argc, &argv);
   const char *only = "";
   int refuse = 0;
   bool clean = false;
   const char *load_deck = "";
   OptionsParser args(argc, argv);
   args.AddOption(&load_deck, "-ld", "--load-deck", "Only read this deck and print [cycle.opening].");
   args.AddOption(&clean, "-c", "--clean", "-no-c", "--no-clean",
                  "Delete the refusals' deck file and stop.");
   args.AddOption(&only, "-o", "--only", "Comma list of T1, T2, T4, C1 (default: all).");
   args.AddOption(&refuse, "-r", "--refuse", "Run refusal n (it must abort).");
   args.Parse();
   if (!args.Good())
   {
      if (mpi.IsRoot()) { args.PrintUsage(std::cout); }
      return 1;
   }
   if (clean) { std::remove("test_fault_opening_refuse.toml"); return 0; }
   if (load_deck[0])
   {
      const cycle::CycleConfig c = cycle::CycleConfig::Load(load_deck);
      if (mpi.IsRoot())
      {
         std::cout << "[cycle.opening] mode " << c.opening.mode << ", max_passes "
                   << c.opening.max_passes << "; [cycle.qd.green] " << (c.qd.green.enabled ? 1 : 0)
                   << ", switching " << (c.switching ? 1 : 0) << ", relaxation "
                   << (c.Relaxation() ? 1 : 0) << ", floor " << c.sigma_n_floor / 1e6
                   << " MPa\n";
      }
      return 0;
   }
   if (refuse > 0) { return Refuse(refuse); }
#ifndef MFEM_USE_PETSC
   if (mpi.IsRoot()) { std::cerr << "test_fault_opening needs PETSc.\n"; }
   return 2;
#else
   const std::string o = only;
   {
      const cycle::CycleConfig c0 = cycle::CycleConfig::Load("simple_cases/configs/D0_straight.toml");
      MFEMInitializePetsc(&argc, &argv, c0.qd.petsc_options.c_str(), NULL);
   }
   if (mpi.IsRoot()) { std::cout << "test_fault_opening on " << mpi.Size() << " ranks\n"; }
   if (Wanted(o, "T1")) { TestT1(mpi); }
   if (Wanted(o, "T2") || Wanted(o, "T4")) { TestT2T4(mpi, Wanted(o, "T2"), Wanted(o, "T4")); }
   if (Wanted(o, "C1")) { TestC1(mpi); }
   if (mpi.IsRoot()) { std::cout << (fails ? "FAILED: " : "ALL PASS: ") << fails << "\n"; }
   MFEMFinalizePetsc();
   return fails ? 1 : 0;
#endif
}
