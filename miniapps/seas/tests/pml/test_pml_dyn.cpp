// tests/pml/test_pml_dyn.cpp - PML plan rev2, P2: the dynamic PML on the
// Nitsche wave operator (Sections 3.1-3.7, tests T2-T4 of Section 7.3).
//
// Box [-H, H]^D, every wall Dirichlet (attribute 5), the band of thickness L
// on every wall.  Operators:
//   K0   the P1 band operator (C0 in the band, penalty x kappa): the QD side
//   K_p  the dynamic operator's stiffness (C, penalty x kappa)
//   B    PmlDynamicTerms: band force from the memory, b_band at entry
// Parts:
//   S0  static identity: K_p u - F_p(g) + StaticBand(u) = K0 u - F0(g) to
//       rounding, random u, nonzero wall data g (2-D and 3-D, p = 1, 2)
//   S1  the filters' zero-frequency limit (Proposition 3, discrete): u_s = 0,
//       a constant w held for alpha t >> 1: B(w, Phi) + E -> StaticBand(w)
//   T2  rest stays at rest: u_s the K0 static solution under a loaded wall
//       (plate shear g = gamma y on y = +-H), Enter, 1e4 steps at the run's
//       dt with the walls held: max |V| / (gamma c_s) <= 1e-12, or within
//       2x of the same run without the band (its rounding floor) when that
//       floor is above 1e-12 (2-D; 2e3 steps in 3-D)
//   T3  a pulse into the band (2-D, p = 2): the inner-box solution against a
//       3x larger box with plain walls (no reflection back by T), and the
//       same small box without the band: the band's reflection measure is
//       printed and must be < 5% of the no-band one
//   T4  no growth: the band box (h = 1 km) run 10 S-wave transits; the
//       inner energy over the second half below 1e-3 of its peak, finite
//   CFL the stable step with K_p (penalty x kappa) against plain K by power
//       iteration: at least 0.95 / sqrt(kappa_max) of plain (the bound the
//       dynamic phase applies); the runs here step at 0.5 x the K_p bound
//
// Usage (from miniapps/seas): ./seas_test_pml_dyn [--only S0,T3,...]

#include "mfem.hpp"
#include "../../domain/nitsche_wave_operator.hpp"
#include "../../solver/nitsche_cfl.hpp"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>

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

/// [-H, H]^D, n cells per side, every wall attribute 5.
Mesh Box(int D, int n, real_t H)
{
   Mesh m = (D == 2) ? Mesh::MakeCartesian2D(n, n, Element::TRIANGLE, true, 2 * H, 2 * H)
                     : Mesh::MakeCartesian3D(n, n, n, Element::TETRAHEDRON, 2 * H, 2 * H, 2 * H);
   for (int v = 0; v < m.GetNV(); v++)
   {
      real_t *x = m.GetVertex(v);
      for (int i = 0; i < D; i++) { x[i] -= H; }
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

using DFn = std::function<void(const Vector &, real_t, Vector &)>;

BoundaryConfig BC(int D, DFn g)
{
   BoundaryConfig bc;
   bc.fault_attr = 3;   // absent
   bc.dirichlet_attrs = {5};
   bc.default_dirichlet_func = [g, D](const Vector &x, real_t t, Vector &u)
   {
      u.SetSize(D);
      g(x, t, u);
   };
   return bc;
}

std::unique_ptr<NitscheWaveOperator> Op(ParMesh &pm, int p, DFn g)
{
   return std::make_unique<NitscheWaveOperator>(
      pm, p, MaterialField::MakeConstant(kLam, kMu, kRho), BC(pm.Dimension(), g));
}

std::shared_ptr<const PmlBand> MakeBand(ParMesh &pm, real_t L, real_t f_alpha)
{
   return std::make_shared<PmlBand>(pm, std::vector<int>{5}, Params(L, f_alpha), pm.GetComm());
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

// ------------------------------------------------------------------ S0
void TestS0(int D, int p)
{
   const real_t H = 10e3, L = 4e3;
   const int n = (D == 2) ? 10 : 5;
   Mesh sm = Box(D, n, H);
   ParMesh pm(MPI_COMM_WORLD, sm);
   auto band = MakeBand(pm, L, 0.3 * Cp() / L);
   DFn g = [D](const Vector &x, real_t t, Vector &u)
   {
      for (int i = 0; i < D; i++)
      { u(i) = 1e-3 * (std::sin(1e-4 * x(0) + 0.3 * i) + 0.5 * x((i + 1) % D) * 1e-4) * (1 + t); }
   };
   auto op0 = Op(pm, p, g);
   op0->Elasticity().SetPmlBand(band, false);   // K0
   auto op1 = Op(pm, p, g);
   op1->SetPml(band);                             // K_p + B
   op1->SetTimeStep(1e-3);
   const int N = op0->Height();
   Vector u(N), r0(N), r1(N), Ku(N), F(N), sb;
   u.Randomize(17 + g_rank);
   u -= 0.5;
   u *= 1e-3;
   const real_t t = 0.25;
   op0->Elasticity().ApplyStiffness(u, Ku);
   op0->Elasticity().AssembleDirichletOnlyRHS(F, t);
   r0 = Ku;
   r0 -= F;
   op1->Elasticity().ApplyStiffness(u, Ku);
   op1->Elasticity().AssembleDirichletOnlyRHS(F, t);
   op1->Pml()->StaticBand(u, t, sb);
   r1 = Ku;
   r1 -= F;
   r1 += sb;
   Vector d(r1);
   d -= r0;
   MPI_Comm c = pm.GetComm();
   const real_t rel = Norm(d, c) / Norm(r0, c);
   int nbe = op1->Pml()->NumBandElements(), nbf = op1->Pml()->NumBandFaces();
   MPI_Allreduce(MPI_IN_PLACE, &nbe, 1, MPI_INT, MPI_SUM, c);
   MPI_Allreduce(MPI_IN_PLACE, &nbf, 1, MPI_INT, MPI_SUM, c);
   std::ostringstream s;
   s << "S0 " << D << "-D p=" << p << ": K_p u - F_p + StaticBand = K0 u - F0 ("
     << nbe << " band elements, " << nbf << " band faces)";
   Check(rel < 1e-12, s.str(), rel);
   // and without the band force they differ (the test sees the band)
   Vector d2(Ku);
   d2 -= F;
   d2 -= r0;
   const real_t rel2 = Norm(d2, c) / Norm(r0, c);
   Check(rel2 > 1e-3, "S0   ... and K_p u - F_p alone differs from K0 u - F0", rel2);
}

// ------------------------------------------------------------------ S1
void TestS1(int D)
{
   const real_t H = 10e3, L = 4e3;
   const int n = (D == 2) ? 10 : 4;
   const int p = 2;
   Mesh sm = Box(D, n, H);
   ParMesh pm(MPI_COMM_WORLD, sm);
   const real_t fa = 0.3 * Cp() / L;
   auto band = MakeBand(pm, L, fa);
   DFn g0 = [](const Vector &, real_t, Vector &u) { u = 0.0; };
   auto op = Op(pm, p, g0);
   op->SetPml(band);
   const real_t alpha = band->Alpha();
   const real_t dt = 0.05 / alpha;
   op->SetTimeStep(dt);
   const int N = op->Height();
   Vector us(N), w(N), F(N), sb;
   us = 0.0;
   op->SetTime(0.0);
   op->PmlEnter(us);
   w.Randomize(5 + g_rank);
   w -= 0.5;
   const int nsteps = static_cast<int>(40.0 / (alpha * dt));
   for (int s = 1; s <= nsteps; s++) { op->Pml()->Advance(w, s * dt); }
   F = 0.0;
   op->Pml()->AddForce(F);
   op->Pml()->StaticBand(w, nsteps * dt, sb);
   Vector d(F);
   d -= sb;
   MPI_Comm c = pm.GetComm();
   const real_t rel = Norm(d, c) / Norm(sb, c);
   std::ostringstream s;
   s << "S1 " << D << "-D: constant w, alpha t = 40 (" << nsteps
     << " steps): B + E -> StaticBand(w) (zero-frequency limit)";
   Check(rel < 1e-10, s.str(), rel);
}

// ------------------------------------------------------------------ T2
// K0 u_s = F0(g) by CG on the assembled K0 (symmetric: C0 keeps the major
// symmetry), then the dynamic run from (u_s, 0).
// One rest run: the static solution of the operator the run starts from
// (K0 with the band, plain K without), then nsteps from (u_s, 0); returns
// max |V| / (gamma c_s).
real_t RestRun(int D, int nsteps, bool with_band)
{
   const real_t H = 10e3, L = 4e3, gam = 1e-4;
   const int n = (D == 2) ? 10 : 4;
   const int p = (D == 2) ? 2 : 1;
   Mesh sm = Box(D, n, H);
   ParMesh pm(MPI_COMM_WORLD, sm);
   auto band = MakeBand(pm, L, 0.3 * Cp() / L);
   DFn g = [gam](const Vector &x, real_t, Vector &u) { u = 0.0; u(0) = gam * x(1); };
   auto op0 = Op(pm, p, g);
   if (with_band) { op0->Elasticity().SetPmlBand(band, false); }
   const int N = op0->Height();
   Vector us(N), F(N);
   op0->Elasticity().AssembleDirichletOnlyRHS(F, 0.0);
   us = 0.0;
   {
      HypreParMatrix &K = op0->Elasticity().ParStiffness();
      HypreBoomerAMG amg(K);
      amg.SetPrintLevel(0);
      amg.SetSystemsOptions(D, true);
      CGSolver cg(pm.GetComm());
      cg.SetOperator(K);
      cg.SetPreconditioner(amg);
      cg.SetRelTol(1e-15);
      cg.SetAbsTol(0.0);
      cg.SetMaxIter(5000);
      cg.SetPrintLevel(0);
      cg.Mult(F, us);
      Vector r(N);
      K.Mult(us, r);
      r -= F;
      if (g_rank == 0 && std::getenv("PML_T2_TRACE"))
      {
         std::printf("    CG %d its, |K u_s - F| / |F| = %.3e\n", cg.GetNumIterations(),
                     Norm(r, pm.GetComm()) / Norm(F, pm.GetComm()));
      }
   }
   // the same step in both runs: 0.5 x the K_p bound
   auto op = Op(pm, p, g);
   op->SetPml(band);
   const real_t dt = 0.5 * 2.0 / std::sqrt(OmegaMax2(*op));
   if (!with_band) { op = Op(pm, p, g); }
   op->SetTimeStep(dt);
   op->SetTime(0.0);
   if (with_band) { op->PmlEnter(us); }
   Vector U(us), V(N);
   V = 0.0;
   real_t t = 0.0;
   CentralDifferenceSolver ode;
   ode.Init(*op);
   MPI_Comm c = pm.GetComm();
   real_t vmax = 0.0;
   const real_t vref = gam * Cs();
   for (int s = 0; s < nsteps; s++)
   {
      real_t dts = dt;
      ode.Step(U, V, t, dts);
      if (s % 100 == 99 || s == nsteps - 1) { vmax = std::max(vmax, GMax(V.Normlinf(), c)); }
      if (g_rank == 0 && std::getenv("PML_T2_TRACE") && (s + 1) % (nsteps / 10) == 0)
      { std::printf("    step %d  max|V|/(gamma c_s) %.3e\n", s + 1, vmax / vref); }
   }
   return vmax / vref;
}

// The plan's bound is 1e-12; the rounding floor of the same run without the
// band is measured too, and the band run must stay within 2x of it when it
// is above 1e-12 (rounding, not the band).
void TestT2(int D, int nsteps)
{
   const real_t vb = RestRun(D, nsteps, true);
   const real_t vp = RestRun(D, nsteps, false);
   const int p = (D == 2) ? 2 : 1;
   if (g_rank == 0)
   {
      std::printf("  T2 %d-D p=%d, %d steps: max |V| / (gamma c_s) band %.3e, no band "
                  "(rounding floor) %.3e\n", D, p, nsteps, vb, vp);
   }
   std::ostringstream s;
   s << "T2 " << D << "-D: the static state stays at rest, max |V| / (gamma c_s) <= "
        "max(1e-12, 2 x the no-band floor)";
   Check(vb <= std::max(1e-12, 2.0 * vp), s.str(), vb);
}

// ------------------------------------------------------------------ T3 / T4
// Values at element centroids inside |x_i| < Hin.
void Centroids(ParMesh &pm, const ParGridFunction &u, real_t Hin,
               std::map<std::pair<long, long>, std::array<real_t, 2>> &out)
{
   out.clear();
   Vector val(2), x(2);
   IntegrationPoint ip;
   ip.Set2(1.0 / 3.0, 1.0 / 3.0);
   for (int e = 0; e < pm.GetNE(); e++)
   {
      ElementTransformation *T = pm.GetElementTransformation(e);
      T->Transform(ip, x);
      if (std::abs(x(0)) >= Hin || std::abs(x(1)) >= Hin) { continue; }
      u.GetVectorValue(e, ip, val);
      out[{std::lround(x(0) * 10), std::lround(x(1) * 10)}] = {val(0), val(1)};
   }
}

struct PulseRun
{
   std::vector<std::map<std::pair<long, long>, std::array<real_t, 2>>> snaps;
   std::vector<real_t> energy;   // inner sum |u|^2 (centroids) per snapshot
   bool finite = true;
};

PulseRun RunPulse(real_t Hbox, real_t h, int p, const std::shared_ptr<const PmlBand> *band_in,
                  real_t L, real_t fa, real_t Hin, real_t dt, real_t T, real_t snap_dt)
{
   const int n = static_cast<int>(std::lround(2 * Hbox / h));
   Mesh sm = Box(2, n, Hbox);
   ParMesh pm(MPI_COMM_WORLD, sm);
   DFn g0 = [](const Vector &, real_t, Vector &u) { u = 0.0; };
   auto op = Op(pm, p, g0);
   std::shared_ptr<const PmlBand> band;
   if (band_in) { band = MakeBand(pm, L, fa); op->SetPml(band); }
   op->SetTimeStep(dt);
   op->SetTime(0.0);
   auto &fes = op->Elasticity().GetFESpace();
   ParGridFunction gu(&fes);
   const real_t x0 = 1.5e3, y0 = 0.8e3, wd = 1.2e3;
   VectorFunctionCoefficient cu(2, [=](const Vector &x, Vector &u)
   {
      const real_t r2 = (x(0) - x0) * (x(0) - x0) + (x(1) - y0) * (x(1) - y0);
      const real_t a = 1e-3 * std::exp(-r2 / (wd * wd));
      u(0) = a;
      u(1) = 0.6 * a;
   });
   gu.ProjectCoefficient(cu);
   Vector U(gu), V(U.Size()), us(U.Size());
   V = 0.0;
   us = 0.0;
   if (band) { op->PmlEnter(us); }
   CentralDifferenceSolver ode;
   ode.Init(*op);
   real_t t = 0.0;
   PulseRun r;
   const int nsteps = static_cast<int>(std::lround(T / dt));
   const int every = std::max(1, static_cast<int>(std::lround(snap_dt / dt)));
   for (int s = 1; s <= nsteps; s++)
   {
      real_t dts = dt;
      ode.Step(U, V, t, dts);
      if (s % every == 0)
      {
         gu = U;
         r.snaps.emplace_back();
         Centroids(pm, gu, Hin, r.snaps.back());
         real_t e = 0.0;
         for (auto &kv : r.snaps.back()) { e += kv.second[0] * kv.second[0] + kv.second[1] * kv.second[1]; }
         r.energy.push_back(e);
         if (!std::isfinite(e)) { r.finite = false; break; }
      }
   }
   return r;
}

// max over snapshots of ||a - ref|| / max over snapshots of ||ref||
real_t Misfit(const PulseRun &a, const PulseRun &ref)
{
   real_t dmax = 0.0, nmax = 0.0;
   const size_t ns = std::min(a.snaps.size(), ref.snaps.size());
   for (size_t k = 0; k < ns; k++)
   {
      real_t d = 0.0, nn = 0.0;
      for (auto &kv : ref.snaps[k])
      {
         auto it = a.snaps[k].find(kv.first);
         if (it == a.snaps[k].end()) { continue; }
         for (int c = 0; c < 2; c++)
         {
            const real_t e = it->second[c] - kv.second[c];
            d += e * e;
            nn += kv.second[c] * kv.second[c];
         }
      }
      dmax = std::max(dmax, std::sqrt(d));
      nmax = std::max(nmax, std::sqrt(nn));
   }
   return dmax / nmax;
}

// 0.5 x the stable step of the band box at mesh size h (power iteration on
// K_p), printed against plain K and the 1 / sqrt(kappa_max) estimate.
real_t BandDt(real_t Hb, real_t h, int p, real_t L, real_t fa)
{
   const int n = static_cast<int>(std::lround(2 * Hb / h));
   Mesh sm = Box(2, n, Hb);
   ParMesh pm(MPI_COMM_WORLD, sm);
   DFn g0 = [](const Vector &, real_t, Vector &u) { u = 0.0; };
   auto opb = Op(pm, p, g0);
   auto band = MakeBand(pm, L, fa);
   opb->SetPml(band);
   auto opp = Op(pm, p, g0);
   const real_t wb = OmegaMax2(*opb), wp = OmegaMax2(*opp);
   const real_t km = opb->Pml()->KappaMax();
   const real_t ratio = std::sqrt(wp / wb);
   if (g_rank == 0)
   {
      std::printf("  CFL h = %.0f m: stable dt with K_p (s0 wall %.2f, kappa_max %.2f) %.4e s, "
                  "plain K %.4e s, ratio %.3f (1 / sqrt(kappa_max) = %.3f)\n", h,
                  band->S0Wall(), km, 2.0 / std::sqrt(wb), 2.0 / std::sqrt(wp), ratio,
                  1.0 / std::sqrt(km));
   }
   Check(ratio >= 0.95 / std::sqrt(km), "CFL: K_p's stable step >= 0.95 / sqrt(kappa_max) "
         "of plain K's (the bound DynPhase applies)", ratio * std::sqrt(km));
   return 0.5 * 2.0 / std::sqrt(wb);
}

void TestT3T4(bool t3, bool t4)
{
   const real_t Hin = 10e3, L = 5e3;
   const int p = 2;
   const real_t fa = 0.2;   // Hz: alpha well below the pulse band (~2 Hz)
   const real_t Hb = Hin + L;
   if (t3)
   {
      const real_t h = 500.0, T = 7.0, sdt = 0.1;
      const real_t dt = BandDt(Hb, h, p, L, fa);
      std::shared_ptr<const PmlBand> dummy;
      PulseRun ref = RunPulse(3 * Hin, h, p, nullptr, 0, 0, Hin, dt, T, sdt);
      PulseRun bnd = RunPulse(Hb, h, p, &dummy, L, fa, Hin, dt, T, sdt);
      PulseRun plain = RunPulse(Hin, h, p, nullptr, 0, 0, Hin, dt, T, sdt);
      const real_t mb = Misfit(bnd, ref), mp = Misfit(plain, ref);
      if (g_rank == 0)
      {
         std::printf("  T3 inner-box misfit vs the 3x box (max over t <= %.0f s): band %.4e, "
                     "no band (Dirichlet walls) %.4e, ratio %.4f\n", T, mb, mp, mb / mp);
      }
      Check(bnd.finite && mb < 0.05 * mp, "T3 the band's reflection < 5% of the bare walls'", mb / mp);
   }
   if (t4)
   {
      // stability, not accuracy: the coarser mesh, 10 S-wave transits of the box
      const real_t h = 1000.0, T = 10.0 * 2 * Hb / Cs(), sdt = 0.5;
      const real_t dt = BandDt(Hb, h, p, L, fa);
      std::shared_ptr<const PmlBand> dummy;
      PulseRun bnd = RunPulse(Hb, h, p, &dummy, L, fa, Hin, dt, T, sdt);
      real_t peak = 0.0;
      for (real_t e : bnd.energy) { peak = std::max(peak, e); }
      const real_t last = bnd.energy.empty() ? 0.0 : bnd.energy.back();
      real_t late_max = 0.0;
      for (size_t k = bnd.energy.size() / 2; k < bnd.energy.size(); k++)
      { late_max = std::max(late_max, bnd.energy[k]); }
      if (g_rank == 0)
      {
         std::printf("  T4 %.0f s, %ld steps (%zu snapshots): inner sum|u|^2 peak %.4e, last "
                     "%.4e, max over the second half %.4e\n", T, std::lround(T / dt),
                     bnd.energy.size(), peak, last, late_max);
      }
      Check(bnd.finite && late_max < 1e-3 * peak,
            "T4 long run: finite, second-half inner energy < 1e-3 of the peak",
            late_max / peak);
   }
}

}  // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   const char *only = "";
   OptionsParser args(argc, argv);
   args.AddOption(&only, "-o", "--only", "Comma-separated parts (S0,S1,T2,T3,T4).");
   args.Parse();
   if (!args.Good()) { if (g_rank == 0) { args.PrintUsage(std::cout); } return 1; }
   {
      std::stringstream ss(only);
      std::string it;
      while (std::getline(ss, it, ',')) { if (!it.empty()) { g_only.insert(it); } }
   }
   if (Want("S0"))
   {
      if (g_rank == 0) { std::cout << "S0: static identity\n"; }
      TestS0(2, 1); TestS0(2, 2); TestS0(3, 1); TestS0(3, 2);
   }
   if (Want("S1"))
   {
      if (g_rank == 0) { std::cout << "S1: zero-frequency limit of the filters\n"; }
      TestS1(2); TestS1(3);
   }
   if (Want("T2"))
   {
      if (g_rank == 0) { std::cout << "T2: rest stays at rest\n"; }
      TestT2(2, 10000); TestT2(3, 2000);
   }
   if (Want("T3") || Want("T4"))
   {
      if (g_rank == 0) { std::cout << "T3/T4: pulse into the band, long run\n"; }
      TestT3T4(Want("T3"), Want("T4"));
   }
   if (g_rank == 0)
   {
      if (num_failed) { std::cout << "FAILED: " << num_failed << " of " << num_tests << "\n"; }
      else { std::cout << "ALL PASS (" << num_tests << " checks)\n"; }
   }
   return num_failed ? 1 : 0;
}
