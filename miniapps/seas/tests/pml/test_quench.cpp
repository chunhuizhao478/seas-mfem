// tests/pml/test_quench.cpp - quench plan rev1 (Chunhui Zhao 2026-11-06),
// T-Q2/T-Q3 at the operator level: the dynamic relaxation C_q = d M_rho of
// NitscheWaveOperator::SetQuench (Section 3 of the plan) on a box [-H, H]^2
// at rest, every wall held (attribute 5, g = 0), under a body load b
// switched on at t = 0 (a Gaussian patch; the settle check of the hand-back
// offset explainer).  The static answer is K u = b, with the band K0 u = b
// (the band's static operator).  d = 2 pi c_s / (2 H), the plan's "auto".
//
// Parts (2-D, p = 1, h = 10 km, H = 240 km):
//   E  the discrete energy identity (plain walls, no Kelvin-Voigt, the
//      quench on from t = 0).  With v^{n+1/2} = (U^{n+1} - U^n) / dt,
//        E^{n+1/2} = 1/2 v^{n+1/2}.M v^{n+1/2} + 1/2 U^n.K U^{n+1}
//                    - 1/2 b.(U^n + U^{n+1}),
//      the central-difference step with the implicit damping (the operator's
//      damped branch) gives exactly
//        E^{n+1/2} - E^{n-1/2} = -dt vbar.C vbar,  vbar = (U^{n+1} - U^{n-1}) / (2 dt),
//      checked to 1e-10 of the energy scale at every step; the dissipation
//      is never negative (C_q is positive definite: the quench only removes
//      energy, plan Section 3.2)
//   D  decay (plain walls, Kelvin-Voigt 0.1 as the decks): the interior
//      misfit vs the static solution against the modal envelope
//      2 (1 + d t / 2) exp(-d t / 2) (every mode of M u'' + d M u' + K u = b
//      decays at least that fast, plan Section 3.3): misfit <= 10 x envelope
//      from 2 tau_1 on (tau_1 = 2 / d), and below 1 % at 5 tau_1; the
//      unquenched run (to 620 s) keeps a misfit above 1 % over 500-620 s,
//      and the quenched misfit at about 600 s is below 1e-3 of it
//   S  the static end state: the quenched misfit at 1200 s (27 tau_1,
//      envelope 1e-10) below 1e-8; with the CFS band (f_alpha 0.01 Hz, L 150
//      km, as the decks) the quenched field settles to the K0 static solution
//      below 1e-6 at 1200 s
//   R  the round trip: SetQuench + ClearQuench after 20 steps, then the run
//      equals the plain run bit for bit (the damped mass and the blocks are
//      put back as they were)
//
// Usage (from miniapps/seas): mpirun -np N ./seas_test_quench [--only E,D,S,R]

#include "mfem.hpp"
#include "../../domain/nitsche_wave_operator.hpp"
#include "../../solver/nitsche_cfl.hpp"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace
{

int g_rank = 0, num_tests = 0, num_failed = 0;
std::set<std::string> g_only;

void Check(bool ok, const std::string &what, real_t v)
{
   num_tests++;
   if (!ok) { num_failed++; }
   if (g_rank == 0) { std::printf("  [%s] %s  (%.4e)\n", ok ? "PASS" : "FAIL", what.c_str(), v); }
}
bool Want(const char *p) { return g_only.empty() || g_only.count(p); }

const real_t kLam = 3.2e10, kMu = 3.2e10, kRho = 2670.0;
real_t Cp() { return std::sqrt((kLam + 2 * kMu) / kRho); }
real_t Cs() { return std::sqrt(kMu / kRho); }

real_t Dot(const Vector &a, const Vector &b, MPI_Comm c)
{
   real_t l = a * b, g = 0;
   MPI_Allreduce(&l, &g, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, c);
   return g;
}
real_t Norm(const Vector &a, MPI_Comm c) { return std::sqrt(Dot(a, a, c)); }
real_t GMax(real_t v, MPI_Comm c)
{
   MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, c);
   return v;
}
bool GSame(const Vector &a, const Vector &b, MPI_Comm c)
{
   int eq = (a.Size() == b.Size()) ? 1 : 0;
   for (int i = 0; eq && i < a.Size(); i++) { if (a(i) != b(i)) { eq = 0; } }
   MPI_Allreduce(MPI_IN_PLACE, &eq, 1, MPI_INT, MPI_MIN, c);
   return eq == 1;
}

/// [-H, H]^2, n cells per side, every wall attribute 5.
Mesh Box(int n, real_t H)
{
   Mesh m = Mesh::MakeCartesian2D(n, n, Element::TRIANGLE, true, 2 * H, 2 * H);
   for (int v = 0; v < m.GetNV(); v++)
   {
      real_t *x = m.GetVertex(v);
      x[0] -= H;
      x[1] -= H;
   }
   for (int b = 0; b < m.GetNBE(); b++) { m.SetBdrAttribute(b, 5); }
   m.SetAttributes();
   return m;
}

PmlBandParams Params(real_t L, real_t f_alpha)
{
   PmlBandParams p;
   p.thickness = L;
   p.power = 2;
   p.R0 = 1e-3;
   p.c_ref = Cp();
   p.f_alpha = f_alpha;
   return p;
}

BoundaryConfig BC()
{
   BoundaryConfig bc;
   bc.fault_attr = 3;   // absent
   bc.dirichlet_attrs = {5};
   bc.default_dirichlet_func = [](const Vector &, real_t, Vector &u)
   {
      u.SetSize(2);
      u = 0.0;
   };
   return bc;
}

std::unique_ptr<NitscheWaveOperator> Op(ParMesh &pm, int p)
{
   return std::make_unique<NitscheWaveOperator>(
      pm, p, MaterialField::MakeConstant(kLam, kMu, kRho), BC());
}

/// omega_max^2 of M^{-1} K (power iteration on the operator's own K and M).
real_t OmegaMax2(NitscheWaveOperator &op)
{
   MPI_Comm c = op.GetMesh().GetComm();
   const int n = op.Height();
   Vector v(n), Kv(n), w(n);
   v.Randomize(4242 + g_rank);
   v -= 0.5;
   v /= Norm(v, c);
   real_t lam = 0, old = 0;
   for (int it = 0; it < 3000; it++)
   {
      op.Elasticity().ApplyStiffness(v, Kv);
      op.Mass().Apply(Kv, w);
      lam = Norm(w, c);
      if (it > 10 && std::abs(lam - old) < 1e-7 * lam) { break; }
      old = lam;
      v = w;
      v /= lam;
   }
   return lam;
}

/// The settle problem: the box, the band (null without), the body load b
/// (x direction, Gaussian patch) and the static solution K u = b (K0 with
/// the band) of the operator opS.
struct Settle
{
   ParMesh *pm = nullptr;
   int p = 1;
   std::shared_ptr<const PmlBand> band;
   Vector b, ust;
   std::vector<int> mask;   ///< interior DOFs (element centroids |x|, |y| < Hin)

   void Build(ParMesh &mesh, int order, real_t L, real_t f_alpha, real_t Hin)
   {
      pm = &mesh;
      p = order;
      MPI_Comm c = pm->GetComm();
      if (L > 0.0)
      {
         band = std::make_shared<PmlBand>(*pm, std::vector<int>{5}, Params(L, f_alpha), c);
      }
      auto opS = Op(*pm, p);
      if (band) { opS->Elasticity().SetPmlBand(band, false); }
      const ParFiniteElementSpace &fes = opS->Elasticity().GetFESpace();
      const int N = opS->Height();
      b.SetSize(N);
      {
         ParLinearForm lf(const_cast<ParFiniteElementSpace *>(&fes));
         VectorFunctionCoefficient fc(2, [](const Vector &x, Vector &f)
         {
            const real_t r2 = x(0) * x(0) + x(1) * x(1);
            f = 0.0;
            f(0) = 1.0 * std::exp(-r2 / (10e3 * 10e3));   // N/m^3
         });
         lf.AddDomainIntegrator(new VectorDomainLFIntegrator(fc));
         lf.Assemble();
         b = lf;
      }
      ust.SetSize(N);
      ust = 0.0;
      {
         HypreParMatrix &K = opS->Elasticity().ParStiffness();
         HypreBoomerAMG amg(K);
         amg.SetPrintLevel(0);
         amg.SetSystemsOptions(2, true);
         CGSolver cg(c);
         cg.SetOperator(K);
         cg.SetPreconditioner(amg);
         cg.SetRelTol(1e-14);
         cg.SetAbsTol(0.0);
         cg.SetMaxIter(20000);
         cg.SetPrintLevel(0);
         cg.Mult(b, ust);
         if (g_rank == 0)
         {
            std::printf("  static solve (%s): CG %d iterations\n", band ? "K0, band" : "plain K",
                        cg.GetNumIterations());
         }
      }
      mask.clear();
      Array<int> vd;
      Vector x(2);
      IntegrationPoint ip;
      ip.Set2(1.0 / 3, 1.0 / 3);
      for (int e = 0; e < pm->GetNE(); e++)
      {
         pm->GetElementTransformation(e)->Transform(ip, x);
         if (std::abs(x(0)) >= Hin || std::abs(x(1)) >= Hin) { continue; }
         fes.GetElementVDofs(e, vd);
         for (int i = 0; i < vd.Size(); i++) { mask.push_back(vd[i]); }
      }
   }
   /// The relative interior misfit |U - u_st| / |u_st| on the mask.
   real_t Misfit(const Vector &U) const
   {
      real_t d = 0, nn = 0;
      for (int i : mask)
      {
         d += (U(i) - ust(i)) * (U(i) - ust(i));
         nn += ust(i) * ust(i);
      }
      real_t a[2] = {d, nn};
      MPI_Allreduce(MPI_IN_PLACE, a, 2, MPITypeMap<real_t>::mpi_type, MPI_SUM, pm->GetComm());
      return std::sqrt(a[0] / a[1]);
   }
   /// A dynamic operator for the problem: the band (if any), dt = 0.5 x the
   /// stable step, Kelvin-Voigt kv (0: off), the load b, from rest.
   std::unique_ptr<NitscheWaveOperator> Dyn(real_t kv, real_t &dt) const
   {
      auto op = Op(*pm, p);
      if (band) { op->SetPml(band); }
      dt = 0.5 * 2.0 / std::sqrt(OmegaMax2(*op));
      op->SetTimeStep(dt);
      if (kv > 0) { op->SetKelvinVoigt(kv, NitscheWaveOperator::KVScale::Node, -1.0); }
      op->SetTime(0.0);
      op->SetConstantLoad(b);
      if (band)
      {
         Vector us(op->Height());
         us = 0.0;
         op->PmlEnter(us);
      }
      return op;
   }
};

/// The modal envelope of plan Section 3.3 for the rate d: every mode of
/// M u'' + d M u' + K u = b is within 2 (1 + d t / 2) exp(-d t / 2) of its
/// static value, relative to its start (the slowest mode critically damped,
/// the others underdamped at exp(-d t / 2) with the sqrt(1 + (d / 2 w)^2) <=
/// sqrt 2 factor of their phase).
real_t Envelope(real_t d, real_t t) { return 2.0 * (1.0 + 0.5 * d * t) * std::exp(-0.5 * d * t); }

} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   const char *only = "";
   OptionsParser args(argc, argv);
   args.AddOption(&only, "-o", "--only", "Comma-separated parts to run (E, D, S, R).");
   args.Parse();
   if (!args.Good())
   {
      if (g_rank == 0) { args.PrintUsage(std::cout); }
      return 1;
   }
   {
      std::stringstream ss(only);
      std::string tok;
      while (std::getline(ss, tok, ',')) { if (!tok.empty()) { g_only.insert(tok); } }
   }

   const real_t H = 240e3, h = 10e3, Hin = 80e3, L = 150e3, fa = 0.01, kv_deck = 0.1;
   const int p = 1, n = static_cast<int>(std::lround(2 * H / h));
   const real_t d = 2.0 * M_PI * Cs() / (2.0 * H);   // the plan's "auto"
   const real_t tau1 = 2.0 / d;
   Mesh sm = Box(n, H);
   ParMesh pm(MPI_COMM_WORLD, sm);
   MPI_Comm c = pm.GetComm();
   if (g_rank == 0)
   {
      std::printf("test_quench: box +-%.0f km, h %.0f km, p %d, %d ranks; d = 2 pi c_s / L = "
                  "%.4e 1/s, tau_1 = 2/d = %.1f s\n", H / 1e3, h / 1e3, p, Mpi::WorldSize(), d,
                  tau1);
   }

   Settle plain;
   plain.Build(pm, p, 0.0, 0.0, Hin);

   // ---- E: the discrete energy identity (plain walls, no KV, the quench on)
   if (Want("E"))
   {
      real_t dt = 0.0;
      auto op = plain.Dyn(0.0, dt);
      op->SetQuench(d);
      const ParFiniteElementSpace &fes = op->Elasticity().GetFESpace();
      ParBilinearForm mf(const_cast<ParFiniteElementSpace *>(&fes));
      ConstantCoefficient rho(kRho);
      mf.AddDomainIntegrator(new VectorMassIntegrator(rho));
      mf.Assemble();
      mf.Finalize();
      const SparseMatrix &M = mf.SpMat();   // block diagonal (DG): local
      const int N = op->Height();
      Vector U(N), V(N), Up(N), Um(N), vh(N), Mv(N), KU(N), vbar(N), Cv(N);
      U = 0.0;
      V = 0.0;
      CentralDifferenceSolver ode;
      ode.Init(*op);
      real_t t = 0.0;
      const int nsteps = static_cast<int>(std::lround(300.0 / dt));
      // E^{n+1/2} from U^n (Um) and U^{n+1} (U)
      auto energy = [&](const Vector &Un, const Vector &Un1)
      {
         vh = Un1;
         vh -= Un;
         vh /= dt;
         M.Mult(vh, Mv);
         op->Elasticity().ApplyStiffness(Un1, KU);
         const real_t kin = 0.5 * Dot(vh, Mv, c), pot = 0.5 * Dot(Un, KU, c);
         const real_t work = 0.5 * (Dot(plain.b, Un, c) + Dot(plain.b, Un1, c));
         return kin + pot - work;
      };
      const real_t escale = 0.5 * std::abs(Dot(plain.b, plain.ust, c));
      real_t worst = 0.0, diss_min = 0.0, e_prev = 0.0;
      Um = U;
      for (int s = 1; s <= nsteps; s++)
      {
         Up = U;   // U^{n-1} before the step (n = s - 1)
         real_t dts = dt;
         ode.Step(U, V, t, dts);   // U^{n+1}
         // after step s: Up = U^{s-1}, Um = U^{s-1}?  keep the last three
         if (s == 1)
         {
            e_prev = energy(Up, U);   // E^{1/2}
            Um = Up;                  // U^0
            continue;
         }
         // Um = U^{s-2}, Up = U^{s-1}, U = U^s: E^{s-1/2} vs E^{s-3/2}
         const real_t e_now = energy(Up, U);
         vbar = U;
         vbar -= Um;
         vbar /= (2.0 * dt);
         op->ApplyDamping(vbar, Cv);
         const real_t diss = dt * Dot(vbar, Cv, c);
         const real_t res = std::abs(e_now - e_prev + diss);
         const real_t scale = std::max({std::abs(e_now), std::abs(e_prev), diss, escale});
         worst = std::max(worst, res / scale);
         diss_min = std::min(diss_min, diss / scale);
         e_prev = e_now;
         Um = Up;
      }
      if (g_rank == 0)
      {
         std::printf("  E: %d steps of %.4e s; identity residual max %.3e of the energy scale, "
                     "smallest dissipation %.3e\n", nsteps, dt, worst, diss_min);
      }
      Check(worst <= 1e-10, "E: E^{n+1/2} - E^{n-1/2} = -dt vbar.C vbar at every step (1e-10)",
            worst);
      Check(diss_min >= 0.0, "E: the quench never adds energy", diss_min);
   }

   // ---- D and S (plain walls, KV as the decks): decay and the static end state
   if (Want("D") || Want("S"))
   {
      const real_t T = 1200.0;
      real_t dt = 0.0;
      auto opq = plain.Dyn(kv_deck, dt);
      opq->SetQuench(d);
      real_t dt0 = 0.0;
      auto op0 = plain.Dyn(kv_deck, dt0);
      const int N = opq->Height();
      Vector Uq(N), Vq(N), U0(N), V0(N);
      Uq = 0.0; Vq = 0.0; U0 = 0.0; V0 = 0.0;
      CentralDifferenceSolver odeq, ode0;
      odeq.Init(*opq);
      ode0.Init(*op0);
      real_t tq = 0.0, t0 = 0.0;
      const int nsteps = static_cast<int>(std::lround(T / dt));
      const int every = std::max(1, static_cast<int>(std::lround(20.0 / dt)));
      // the unquenched comparison run goes to T0 (its last report, at about
      // 600 s, gives m600q and m600_0)
      const real_t T0 = 620.0;
      real_t worst_env = 0.0, m5 = -1.0, m600q = -1.0, m600_0 = -1.0, last0 = 0.0, lastq = 0.0;
      if (g_rank == 0)
      {
         std::printf("  D/S: plain walls, KV %.2f, dt %.4e s, %d steps to %.0f s\n", kv_deck, dt,
                     nsteps, T);
         std::printf("      t [s]   misfit quenched   envelope   misfit plain\n");
      }
      for (int s = 1; s <= nsteps; s++)
      {
         real_t dts = dt;
         odeq.Step(Uq, Vq, tq, dts);
         if (s * dt <= T0 + 0.5 * dt)
         {
            dts = dt;
            ode0.Step(U0, V0, t0, dts);
         }
         const bool report = (s % every == 0 || s == nsteps);
         const bool near5 = (m5 < 0.0 && tq >= 5.0 * tau1);
         if (!report && !near5) { continue; }
         const real_t mq = plain.Misfit(Uq);
         const real_t env = Envelope(d, tq);
         if (tq >= 2.0 * tau1) { worst_env = std::max(worst_env, mq / env); }
         if (near5) { m5 = mq; }
         if (s * dt <= T0 + 0.5 * dt)
         {
            const real_t m0 = plain.Misfit(U0);
            if (tq >= 500.0) { last0 = std::max(last0, m0); }
            if (report && g_rank == 0)
            {
               std::printf("   %8.1f   %.4e        %.3e   %.4e\n", tq, mq, env, m0);
            }
            if (report)
            {
               m600q = mq;   // the last report with both runs: about 600 s
               m600_0 = m0;
            }
         }
         else
         {
            if (tq >= T - 100.0) { lastq = std::max(lastq, mq); }
            if (report && g_rank == 0) { std::printf("   %8.1f   %.4e        %.3e\n", tq, mq, env); }
         }
      }
      if (Want("D"))
      {
         Check(worst_env <= 10.0, "D: misfit <= 10 x the modal envelope from 2 tau_1 on", worst_env);
         Check(m5 >= 0.0 && m5 < 1e-2, "D: misfit below 1 % at 5 tau_1", m5);
         Check(last0 > 1e-2, "D: the unquenched misfit over 500-600 s stays above 1 %", last0);
         Check(m600q >= 0.0 && m600_0 > 0.0 && m600q < 1e-3 * m600_0,
               "D: quenched misfit at about 600 s below 1e-3 of the unquenched", m600q / m600_0);
      }
      if (Want("S"))
      {
         const real_t mend = plain.Misfit(Uq);
         Check(mend < 1e-8, "S: the quenched field at 1200 s = the static solution (1e-8)", mend);
         Check(lastq < 1e-8, "S: ... and over 1100-1200 s", lastq);
      }
   }

   // ---- S with the band: the quenched field settles to the K0 static solution
   if (Want("S"))
   {
      Settle band;
      band.Build(pm, p, L, fa, Hin);
      const real_t T = 1200.0;
      real_t dt = 0.0;
      auto opq = band.Dyn(kv_deck, dt);
      opq->SetQuench(d);
      const int N = opq->Height();
      Vector Uq(N), Vq(N);
      Uq = 0.0;
      Vq = 0.0;
      CentralDifferenceSolver odeq;
      odeq.Init(*opq);
      real_t tq = 0.0;
      const int nsteps = static_cast<int>(std::lround(T / dt));
      const int every = std::max(1, static_cast<int>(std::lround(100.0 / dt)));
      if (g_rank == 0)
      {
         std::printf("  S (band): CFS band L %.0f km, f_alpha %.3f Hz, KV %.2f, dt %.4e s, %d "
                     "steps\n", L / 1e3, fa, kv_deck, dt, nsteps);
      }
      real_t m600 = -1.0;
      for (int s = 1; s <= nsteps; s++)
      {
         real_t dts = dt;
         odeq.Step(Uq, Vq, tq, dts);
         if (s % every == 0 || s == nsteps)
         {
            const real_t mq = band.Misfit(Uq);
            if (m600 < 0.0 && tq >= 600.0 - 0.5 * dt) { m600 = mq; }
            if (g_rank == 0) { std::printf("   %8.1f   misfit vs K0 static %.4e\n", tq, mq); }
         }
      }
      const real_t mend = band.Misfit(Uq);
      Check(std::isfinite(mend), "S (band): finite", mend);
      Check(mend < 1e-6, "S (band): the quenched field at 1200 s = the K0 static solution (1e-6)",
            mend);
   }

   // ---- R: the round trip
   if (Want("R"))
   {
      real_t dta = 0.0, dtb = 0.0;
      auto opa = plain.Dyn(kv_deck, dta);
      auto opb = plain.Dyn(kv_deck, dtb);
      const int N = opa->Height();
      Vector Ua(N), Va(N), Ub(N), Vb(N);
      Ua = 0.0; Va = 0.0; Ub = 0.0; Vb = 0.0;
      CentralDifferenceSolver odea, odeb;
      odea.Init(*opa);
      odeb.Init(*opb);
      real_t ta = 0.0, tb = 0.0;
      for (int s = 1; s <= 60; s++)
      {
         real_t dts = dta;
         odea.Step(Ua, Va, ta, dts);
         if (s == 21)
         {
            opb->SetQuench(d);
            opb->ClearQuench();
         }
         dts = dtb;
         odeb.Step(Ub, Vb, tb, dts);
      }
      Check(!opb->Quenched() && !opb->HasDamping() && opb->DampedMass() == nullptr,
            "R: no damping block and no damped mass after ClearQuench (none before)", 0.0);
      Check(GSame(Ua, Ub, c) && GSame(Va, Vb, c), "R: the run after the round trip = the plain run "
            "bit for bit", 0.0);
   }

   if (g_rank == 0)
   {
      std::printf("%s: %d of %d checks failed\n", num_failed ? "FAILED" : "ALL PASS", num_failed,
                  num_tests);
   }
   return num_failed ? 1 : 0;
}
