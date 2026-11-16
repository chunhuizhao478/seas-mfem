// tests/pml/test_fw_damping.cpp - the damper on the 2.5D term
// ([cycle.finite_width] damping_ratio, option A of the ring-train note;
// Chunhui Zhao 2026-11-15: "let's try option A first"), at the operator
// level: NitscheWaveOperator::SetMassDamping with the element rates
// d = 2 zeta omega_c, omega_c = sqrt(beta / rho), next to SetFiniteWidth:
//
//   rho u'' = div sigma - beta (u - u_pl) - 2 zeta sqrt(beta rho) u'
//
// W = 15 km, c = (1 - nu) pi / 2 ("wa2020"), zeta = 0.2 as the 1-D study of
// the note; 2-D, p = 1, h = 10 km.
//
// Parts:
//   U  the uniform mode (free walls, u_pl = 0, start: a uniform displacement
//      at rest): K u = 0 on it, so every DOF is the damped oscillator
//        u(t) = u0 exp(-zeta omega_c t) (cos w t + (zeta omega_c / w) sin w t),
//        w = omega_c sqrt(1 - zeta^2),
//      over five cutoff periods to 2e-3 u0 (the step is 0.05 x the stable
//      one); the field stays uniform; without the damper the same start
//      keeps its amplitude (the ring of the note)
//   E  the discrete energy identity (held walls, a body load from rest):
//        E^{n+1/2} - E^{n-1/2} = -dt vbar.C vbar,
//      with K + M_beta in the energy, to 1e-10 at every step; the damper
//      never adds energy
//   D  decay and the static state (held walls, the body load from rest):
//      the residual b - (K + M_beta) U in the M^{-1} norm is below
//      1.1 exp(-zeta omega_c t) of its start at every report (every mode is
//      underdamped with the rate d / 2: its frequency is at least omega_c),
//      and below 1e-6 at the end: the static state is the one of the
//      undamped operator; without the damper the residual stays above 1e-2
//   R  the quench on top (SetQuench + ClearQuench after 20 steps): the run
//      equals the damped run bit for bit (the quench puts the damper's
//      blocks back)
//
// Usage (from miniapps/seas): mpirun -np N ./seas_test_fw_damping [--only U,E,D,R]

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

const real_t kLam = 3.2e10, kMu = 3.2e10, kRho = 2670.0, kW = 15e3, kZeta = 0.2;
real_t BetaOverMu()
{
   const real_t nu = 0.5 * kLam / (kLam + kMu), c = 0.5 * M_PI * (1.0 - nu);
   return 4.0 * c * c / (kW * kW);
}
real_t OmegaC() { return std::sqrt(BetaOverMu() * kMu / kRho); }

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

/// held = true: every wall Dirichlet, g = 0; false: every wall free.
BoundaryConfig BC(bool held)
{
   BoundaryConfig bc;
   bc.fault_attr = 3;   // absent
   if (held) { bc.dirichlet_attrs = {5}; }
   bc.default_dirichlet_func = [](const Vector &, real_t, Vector &u)
   {
      u.SetSize(2);
      u = 0.0;
   };
   return bc;
}

/// omega_max^2 of M^{-1} (K + M_beta) (power iteration).
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

/// The operator with the 2.5D term (u_pl = 0), the damper of ratio zeta
/// (0: none) and the step frac x the stable one.
std::unique_ptr<NitscheWaveOperator> Op(ParMesh &pm, bool held, real_t zeta, real_t frac,
                                        real_t &dt)
{
   auto op = std::make_unique<NitscheWaveOperator>(
      pm, 1, MaterialField::MakeConstant(kLam, kMu, kRho), BC(held));
   Vector Upl(op->Height());
   Upl = 0.0;
   op->SetFiniteWidth(BetaOverMu(), Upl, nullptr);
   if (zeta > 0.0)
   {
      Vector rate(pm.GetNE());
      rate = 2.0 * zeta * OmegaC();
      op->SetMassDamping(rate);
   }
   dt = frac * 2.0 / std::sqrt(OmegaMax2(*op));
   op->SetTimeStep(dt);
   op->SetTime(0.0);
   return op;
}

/// The body load of the settle problem (x direction, Gaussian patch).
void BodyLoad(NitscheWaveOperator &op, Vector &b)
{
   const ParFiniteElementSpace &fes = op.Elasticity().GetFESpace();
   ParLinearForm lf(const_cast<ParFiniteElementSpace *>(&fes));
   VectorFunctionCoefficient fc(2, [](const Vector &x, Vector &f)
   {
      const real_t r2 = x(0) * x(0) + x(1) * x(1);
      f = 0.0;
      f(0) = 1.0 * std::exp(-r2 / (20e3 * 20e3));   // N/m^3
   });
   lf.AddDomainIntegrator(new VectorDomainLFIntegrator(fc));
   lf.Assemble();
   b.SetSize(op.Height());
   b = lf;
}

} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   const char *only = "";
   const char *dump = "";
   OptionsParser args(argc, argv);
   args.AddOption(&dump, "-dump", "--dump", "Prefix for the time series of U and D (text "
                  "files <prefix>_U<run>.txt, <prefix>_D<run>.txt; rank 0, 1 rank); empty: none.");
   args.AddOption(&only, "-o", "--only", "Comma-separated parts to run (U, E, D, R).");
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

   const real_t H = 120e3, h = 10e3;
   const int n = static_cast<int>(std::lround(2 * H / h));
   const real_t wc = OmegaC(), d = 2.0 * kZeta * wc;
   Mesh sm = Box(n, H);
   ParMesh pm(MPI_COMM_WORLD, sm);
   MPI_Comm c = pm.GetComm();
   if (g_rank == 0)
   {
      std::printf("test_fw_damping: box +-%.0f km, h %.0f km, p 1, %d ranks; W %.0f km, "
                  "omega_c = %.4f rad/s (period %.2f s), zeta = %.2f, d = 2 zeta omega_c = "
                  "%.4f 1/s\n", H / 1e3, h / 1e3, Mpi::WorldSize(), kW / 1e3, wc,
                  2.0 * M_PI / wc, kZeta, d);
   }

   // ---- U: the uniform mode against the damped oscillator
   if (Want("U"))
   {
      const real_t T = 5.0 * 2.0 * M_PI / wc, u0 = 1.0;
      real_t worst[2] = {0.0, 0.0}, spread[2] = {0.0, 0.0}, last_amp[2] = {0.0, 0.0};
      for (int run = 0; run < 2; run++)
      {
         const real_t zeta = (run == 0) ? kZeta : 0.0;
         real_t dt = 0.0;
         auto op = Op(pm, false, zeta, 0.05, dt);
         const ParFiniteElementSpace &fes = op->Elasticity().GetFESpace();
         ParGridFunction gf(const_cast<ParFiniteElementSpace *>(&fes));
         Vector e1(2);
         e1(0) = u0;
         e1(1) = 0.0;
         VectorConstantCoefficient cf(e1);
         gf.ProjectCoefficient(cf);
         const int N = op->Height();
         Vector U(N), V(N), U0(N);
         U = gf;
         U0 = U;
         V = 0.0;
         CentralDifferenceSolver ode;
         ode.Init(*op);
         real_t t = 0.0;
         const int nsteps = static_cast<int>(std::lround(T / dt));
         const real_t wd = wc * std::sqrt(1.0 - zeta * zeta);
         const real_t t_last = T - 2.0 * M_PI / wc;
         FILE *fd = (dump[0] && g_rank == 0)
                    ? std::fopen((std::string(dump) + "_U" + std::to_string(run) + ".txt").c_str(), "w")
                    : nullptr;
         for (int s = 0; s < nsteps; s++)
         {
            real_t dts = dt;
            ode.Step(U, V, t, dts);
            const real_t f = std::exp(-zeta * wc * t)
                             * (std::cos(wd * t) + zeta * wc / wd * std::sin(wd * t));
            real_t err = 0.0, sp = 0.0, amp = 0.0;
            for (int i = 0; i < N; i++)
            {
               err = std::max(err, std::abs(U(i) - f * U0(i)));
               amp = std::max(amp, std::abs(U(i)));
            }
            // uniform: the x DOFs all equal (the y DOFs zero is in err)
            real_t lo = std::numeric_limits<real_t>::infinity(), hi = -lo;
            for (int i = 0; i < N; i++)
            {
               if (U0(i) == 0.0) { continue; }
               lo = std::min(lo, U(i));
               hi = std::max(hi, U(i));
            }
            if (N > 0) { sp = hi - lo; }
            if (fd && N > 0) { std::fprintf(fd, "%.9e %.9e %.9e\n", t, hi, f * u0); }
            worst[run] = std::max(worst[run], GMax(err, c));
            spread[run] = std::max(spread[run], GMax(sp, c));
            if (t >= t_last) { last_amp[run] = std::max(last_amp[run], GMax(amp, c)); }
         }
         if (fd) { std::fclose(fd); }
         if (g_rank == 0)
         {
            std::printf("  U (zeta %.2f): %d steps of %.4e s to %.1f s; max |U - exact| %.3e, "
                        "spread over the DOFs %.3e, amplitude in the last period %.4e\n", zeta,
                        nsteps, dt, t, worst[run], spread[run], last_amp[run]);
         }
      }
      const real_t want_last = std::exp(-kZeta * wc * 4.0 * 2.0 * M_PI / wc);
      Check(worst[0] <= 2e-3 * u0, "U: the damped uniform mode is the oscillator "
            "exp(-zeta omega_c t)(cos w t + ...) to 2e-3 u0", worst[0]);
      Check(spread[0] <= 1e-9 * u0, "U: the field stays uniform (1e-9 u0)", spread[0]);
      Check(last_amp[0] <= 1.05 * want_last * u0, "U: amplitude after four periods <= 1.05 "
            "exp(-8 pi zeta) u0", last_amp[0] / (want_last * u0));
      Check(worst[1] <= 2e-3 * u0 && last_amp[1] >= 0.99 * u0, "U control: without the damper "
            "the mode keeps its amplitude (the ring)", last_amp[1]);
   }

   // ---- E: the discrete energy identity (held walls, the body load from rest)
   if (Want("E"))
   {
      real_t dt = 0.0;
      auto op = Op(pm, true, kZeta, 0.5, dt);
      Vector b;
      BodyLoad(*op, b);
      op->SetConstantLoad(b);
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
      const int nsteps = static_cast<int>(std::lround(120.0 / dt));
      auto energy = [&](const Vector &Un, const Vector &Un1)
      {
         vh = Un1;
         vh -= Un;
         vh /= dt;
         M.Mult(vh, Mv);
         op->Elasticity().ApplyStiffness(Un1, KU);   // K + M_beta
         const real_t kin = 0.5 * Dot(vh, Mv, c), pot = 0.5 * Dot(Un, KU, c);
         const real_t work = 0.5 * (Dot(b, Un, c) + Dot(b, Un1, c));
         return kin + pot - work;
      };
      std::vector<real_t> res, dis;
      real_t e_prev = 0.0, escale = 0.0;
      Um = U;
      for (int s = 1; s <= nsteps; s++)
      {
         Up = U;
         real_t dts = dt;
         ode.Step(U, V, t, dts);
         if (s == 1)
         {
            e_prev = energy(Up, U);
            Um = Up;
            escale = std::abs(e_prev);
            continue;
         }
         const real_t e_now = energy(Up, U);
         vbar = U;
         vbar -= Um;
         vbar /= (2.0 * dt);
         op->ApplyDamping(vbar, Cv);
         const real_t diss = dt * Dot(vbar, Cv, c);
         res.push_back(std::abs(e_now - e_prev + diss));
         dis.push_back(diss);
         escale = std::max(escale, std::abs(e_now));
         e_prev = e_now;
         Um = Up;
      }
      real_t worst = 0.0, diss_min = 0.0, diss_sum = 0.0;
      for (std::size_t i = 0; i < res.size(); i++)
      {
         worst = std::max(worst, res[i] / escale);
         diss_min = std::min(diss_min, dis[i] / escale);
         diss_sum += dis[i] / escale;
      }
      if (g_rank == 0)
      {
         std::printf("  E: %d steps of %.4e s; identity residual max %.3e of the energy scale, "
                     "smallest dissipation %.3e, total %.3e\n", nsteps, dt, worst, diss_min,
                     diss_sum);
      }
      Check(worst <= 1e-10, "E: E^{n+1/2} - E^{n-1/2} = -dt vbar.C vbar at every step (1e-10)",
            worst);
      Check(diss_min >= 0.0 && diss_sum > 0.1, "E: the damper never adds energy, and removes "
            "some", diss_min);
   }

   // ---- D: decay to the static state of the undamped operator
   if (Want("D"))
   {
      const real_t T = 16.5 / (kZeta * wc);   // exp(-zeta omega_c T) = 7e-8
      real_t end[2] = {0.0, 0.0}, worst_ratio = 0.0;
      for (int run = 0; run < 2; run++)
      {
         const real_t zeta = (run == 0) ? kZeta : 0.0;
         real_t dt = 0.0;
         auto op = Op(pm, true, zeta, 0.5, dt);
         Vector b;
         BodyLoad(*op, b);
         op->SetConstantLoad(b);
         const int N = op->Height();
         Vector U(N), V(N), KU(N), r(N), Mr(N);
         U = 0.0;
         V = 0.0;
         CentralDifferenceSolver ode;
         ode.Init(*op);
         real_t t = 0.0;
         auto resid = [&]()
         {
            op->Elasticity().ApplyStiffness(U, KU);
            r = b;
            r -= KU;
            op->Mass().Apply(r, Mr);
            return std::sqrt(Dot(r, Mr, c));
         };
         const real_t r0 = resid();
         FILE *fd = (dump[0] && g_rank == 0)
                    ? std::fopen((std::string(dump) + "_D" + std::to_string(run) + ".txt").c_str(), "w")
                    : nullptr;
         const int nsteps = static_cast<int>(std::lround(T / dt));
         const int every = std::max(1, static_cast<int>(std::lround(1.0 / dt)));
         for (int s = 1; s <= nsteps; s++)
         {
            real_t dts = dt;
            ode.Step(U, V, t, dts);
            if (s % every && s != nsteps) { continue; }
            const real_t rr = resid() / r0;
            if (fd) { std::fprintf(fd, "%.9e %.9e\n", t, rr); }
            if (run == 0)
            {
               // the last step of the residual is rounding: compare above 1e-9
               const real_t env = 1.1 * std::exp(-zeta * wc * t);
               if (env > 1e-9) { worst_ratio = std::max(worst_ratio, rr / env); }
            }
            end[run] = rr;
         }
         if (fd) { std::fclose(fd); }
         if (g_rank == 0)
         {
            std::printf("  D (zeta %.2f): %d steps of %.4e s to %.1f s; residual at the end "
                        "%.3e of its start%s\n", zeta, nsteps, dt, t, end[run],
                        run == 0 ? "" : " (no damper)");
         }
      }
      if (g_rank == 0)
      {
         std::printf("  D: largest residual / (1.1 exp(-zeta omega_c t)) over the reports %.4f\n",
                     worst_ratio);
      }
      Check(worst_ratio <= 1.0, "D: residual <= 1.1 exp(-zeta omega_c t) x its start at every "
            "report", worst_ratio);
      Check(end[0] <= 1e-6, "D: the damped run ends at the static state of K + M_beta (1e-6)",
            end[0]);
      Check(end[1] >= 1e-2, "D control: without the damper the residual stays above 1e-2",
            end[1]);
   }

   // ---- R: the quench on top of the damper, and back
   if (Want("R"))
   {
      real_t dt = 0.0;
      Vector UA, UB, UC;
      bool flags = true;
      for (int run = 0; run < 3; run++)
      {
         auto op = Op(pm, true, kZeta, 0.5, dt);
         Vector b;
         BodyLoad(*op, b);
         op->SetConstantLoad(b);
         const int N = op->Height();
         Vector U(N), V(N);
         U = 0.0;
         V = 0.0;
         CentralDifferenceSolver ode;
         ode.Init(*op);
         real_t t = 0.0;
         for (int s = 0; s < 60; s++)
         {
            if (s == 20 && run >= 1)
            {
               op->SetQuench(0.5);
               flags = flags && op->Quenched() && op->HasMassDamping();
               if (run == 1) { op->ClearQuench(); flags = flags && !op->Quenched(); }
            }
            real_t dts = dt;
            ode.Step(U, V, t, dts);
         }
         (run == 0 ? UA : run == 1 ? UB : UC) = U;
      }
      Check(flags && GSame(UA, UB, c), "R: SetQuench + ClearQuench on the damped run: the damped "
            "run bit for bit", 0.0);
      Check(!GSame(UA, UC, c), "R: the quench left on changes the damped run", 0.0);
   }

   if (g_rank == 0)
   {
      std::printf("%s: %d of %d checks failed\n", num_failed ? "FAILED" : "ALL PASS", num_failed,
                  num_tests);
   }
   return num_failed ? 1 : 0;
}
