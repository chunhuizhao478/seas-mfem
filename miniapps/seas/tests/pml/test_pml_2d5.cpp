// tests/pml/test_pml_2d5.cpp - the 2.5D finite-width term inside the PML
// band (2.5D with the PML, piece F-a; plan
// /mnt/project-files/fully-dynamic-fault-opening/finite-width-plan/
// finite-width-pml-lts-plan-rev1.pdf, Section 4.1).
//
// Box [-H, H]^2, every wall Dirichlet (attribute 5), the band on every wall,
// the 2.5D term beta (u - u_pl), beta = beta_over_mu mu, u_pl(t) = (t_off +
// t) U_pl1.  Operators:
//   K0   the QD band operator with the S0-weighted M_beta and substrate load
//   K_p  the dynamic stiffness with the plain M_beta and substrate load
//   B    PmlDynamicTerms with SetFiniteWidth: the band part of the term
// Parts:
//   F0  static identity: K_p u - F_p(t) + [band force] = K0 u - F0(t) to
//       rounding, random u, wall data and substrate nonzero (p = 1, 2); and
//       the band part is seen (without it the two differ)
//   F1  zero-frequency limit: u_s = 0, a constant w held for alpha t >> 1:
//       AddForce -> StaticBand(w) - (t_off + t) [substrate band part]
//   F2  rest stays at rest: u_s the QD solution (K0, substrate load), Enter,
//       1e4 steps: max |V| / (gamma c_s) <= max(1e-12, 2 x the no-band floor)
//   F3  a pulse into the band, W = 15 km (cutoff 0.087 Hz) and W = 3 km
//       (cutoff 0.43 Hz, inside the pulse's band): inner-box misfit against
//       a 3x larger plain box; the band must be < 5 % of the bare walls',
//       and the band with the term left unstretched is printed beside it
//   F4  no growth: 10 S-wave transits; W = 15 km: finite, second-half inner
//       energy < 1e-3 of the peak; W = 3 km (waves near the cutoff leave
//       slowly): finite and the energy keeps falling
//
// Usage (from miniapps/seas): ./seas_test_pml_2d5 [--only F0,F1,...]
// Debugging: PML_F3_TRACE=1 prints the F3 runs' sizes; PML_F3_T=<s> shortens F3.

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


const real_t kWa = 1.178097245;   // (1 - nu) pi / 2, nu = 0.25 (wa2020)
real_t BetaOverMu(real_t W) { return 4.0 * kWa * kWa / (W * W); }

/// A smooth substrate shape on the operator's vdofs (a sheared block).
void Substrate(NitscheWaveOperator &op, real_t amp, Vector &U)
{
   ParGridFunction g(&op.Elasticity().GetFESpace());
   VectorFunctionCoefficient c(2, [amp](const Vector &x, Vector &u)
   {
      u(0) = amp * std::tanh(x(1) / 3e3);
      u(1) = 0.2 * amp * std::sin(2e-4 * x(0));
   });
   g.ProjectCoefficient(c);
   U = g;
}

// ------------------------------------------------------------------ F0
void TestF0(int p)
{
   const real_t H = 10e3, L = 4e3, W = 15e3;
   Mesh sm = Box(2, 10, H);
   ParMesh pm(MPI_COMM_WORLD, sm);
   auto band = MakeBand(pm, L, 0.3 * Cp() / L);
   DFn g = [](const Vector &x, real_t t, Vector &u)
   {
      for (int i = 0; i < 2; i++)
      { u(i) = 1e-3 * (std::sin(1e-4 * x(0) + 0.3 * i) + 0.5 * x((i + 1) % 2) * 1e-4) * (1 + t); }
   };
   auto toff = std::make_shared<const real_t>(3.0);
   auto op0 = Op(pm, p, g);
   Vector Upl;
   Substrate(*op0, 2e-4, Upl);
   // both orders of the two calls on the QD side
   auto op0b = Op(pm, p, g);
   op0->Elasticity().SetPmlBand(band, false);
   op0->Elasticity().SetFiniteWidth(BetaOverMu(W), Upl, toff);
   op0b->Elasticity().SetFiniteWidth(BetaOverMu(W), Upl, toff);
   op0b->Elasticity().SetPmlBand(band, false);
   auto op1 = Op(pm, p, g);
   op1->SetPml(band);
   op1->SetFiniteWidth(BetaOverMu(W), Upl, toff);
   op1->SetTimeStep(1e-3);
   const int N = op0->Height();
   Vector u(N), r0(N), r0b(N), r1(N), Ku(N), F(N), sb, zero(N);
   u.Randomize(17 + g_rank);
   u -= 0.5;
   u *= 1e-3;
   const real_t t = 0.25;
   op0->Elasticity().ApplyStiffness(u, Ku);
   op0->Elasticity().AssembleDirichletOnlyRHS(F, t);
   r0 = Ku;
   r0 -= F;
   op0b->Elasticity().ApplyStiffness(u, Ku);
   op0b->Elasticity().AssembleDirichletOnlyRHS(F, t);
   r0b = Ku;
   r0b -= F;
   r0b -= r0;
   MPI_Comm c = pm.GetComm();
   Check(Norm(r0b, c) == 0.0, "F0   SetPmlBand and SetFiniteWidth in either order: the same "
         "K0 u - F0", Norm(r0b, c));
   // the dynamic side at the switch: Enter(u) at time t, memory zero
   op1->SetTime(t);
   op1->PmlEnter(u);
   op1->Elasticity().ApplyStiffness(u, Ku);
   op1->Elasticity().AssembleDirichletOnlyRHS(F, t);
   Vector fb(N);
   fb = 0.0;
   op1->Pml()->AddForce(fb);
   r1 = Ku;
   r1 -= F;
   r1 += fb;
   Vector d(r1);
   d -= r0;
   const real_t rel = Norm(d, c) / Norm(r0, c);
   std::ostringstream s;
   s << "F0 p=" << p << ": at entry K_p u - F_p + band force = K0 u - F0 (the QD residual)";
   Check(rel < 1e-12, s.str(), rel);
   // the band part of the 2.5D term is seen: the plain band force (no
   // SetFiniteWidth on the band) misses it
   auto op2 = Op(pm, p, g);
   op2->SetPml(band);
   op2->Elasticity().SetFiniteWidth(BetaOverMu(W), Upl, toff);
   op2->SetTimeStep(1e-3);
   op2->SetTime(t);
   op2->PmlEnter(u);
   Vector fb2(N);
   fb2 = 0.0;
   op2->Pml()->AddForce(fb2);
   fb2 -= fb;
   const real_t rel2 = Norm(fb2, c) / Norm(r0, c);
   Check(rel2 > 1e-6, "F0   ... and the band force without the band's 2.5D part differs", rel2);
}

// ------------------------------------------------------------------ F1
void TestF1()
{
   const real_t H = 10e3, L = 4e3, W = 15e3;
   const int p = 2;
   Mesh sm = Box(2, 10, H);
   ParMesh pm(MPI_COMM_WORLD, sm);
   auto band = MakeBand(pm, L, 0.3 * Cp() / L);
   DFn g0 = [](const Vector &, real_t, Vector &u) { u = 0.0; };
   auto toff = std::make_shared<const real_t>(0.0);
   auto op = Op(pm, p, g0);
   Vector Upl, zero;
   Substrate(*op, 0.0, Upl);   // zero substrate: the memory part alone
   op->SetPml(band);
   op->SetFiniteWidth(BetaOverMu(W), Upl, toff);
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
   op->Pml()->StaticBand(w, nsteps * dt, sb);   // includes beta (S0 - 1) w
   Vector d(F);
   d -= sb;
   MPI_Comm c = pm.GetComm();
   const real_t rel = Norm(d, c) / Norm(sb, c);
   Check(rel < 1e-10, "F1: constant w, alpha t = 40: band force -> its zero-frequency value "
         "(elastic part + beta (S0 - 1) w)", rel);
   // the 2.5D share of that static force is not negligible
   auto opn = Op(pm, p, g0);
   opn->SetPml(band);
   opn->SetTimeStep(dt);
   Vector sbn;
   opn->Pml()->StaticBand(w, nsteps * dt, sbn);
   sbn -= sb;
   Check(Norm(sbn, c) > 0.0, "F1   ... the 2.5D part of the static band force is there "
         "(norm relative to the whole)", Norm(sbn, c) / Norm(sb, c));
}

// ------------------------------------------------------------------ F2
real_t RestRun(int nsteps, bool with_band)
{
   const real_t H = 10e3, L = 4e3, gam = 1e-4, W = 15e3;
   const int p = 2;
   Mesh sm = Box(2, 10, H);
   ParMesh pm(MPI_COMM_WORLD, sm);
   auto band = MakeBand(pm, L, 0.3 * Cp() / L);
   DFn g = [gam](const Vector &x, real_t, Vector &u) { u = 0.0; u(0) = gam * x(1); };
   // u_pl = (t_off + t) U_pl1 = about 1 m, and its change over the run is
   // below rounding (rest means rest)
   const real_t T0 = 1e15;
   auto toff = std::make_shared<const real_t>(T0);
   auto op0 = Op(pm, p, g);
   Vector Upl;
   Substrate(*op0, 1.0 / T0, Upl);
   if (with_band) { op0->Elasticity().SetPmlBand(band, false); }
   op0->Elasticity().SetFiniteWidth(BetaOverMu(W), Upl, toff);
   const int N = op0->Height();
   Vector us(N), F(N);
   op0->Elasticity().AssembleDirichletOnlyRHS(F, 0.0);
   us = 0.0;
   {
      HypreParMatrix &K = op0->Elasticity().ParStiffness();
      HypreBoomerAMG amg(K);
      amg.SetPrintLevel(0);
      amg.SetSystemsOptions(2, true);
      CGSolver cg(pm.GetComm());
      cg.SetOperator(K);
      cg.SetPreconditioner(amg);
      cg.SetRelTol(1e-15);
      cg.SetAbsTol(0.0);
      cg.SetMaxIter(5000);
      cg.SetPrintLevel(0);
      cg.Mult(F, us);
   }
   auto opd = Op(pm, p, g);
   opd->SetPml(band);
   const real_t dt = 0.5 * 2.0 / std::sqrt(OmegaMax2(*opd));
   auto op = Op(pm, p, g);
   if (with_band) { op->SetPml(band); }
   op->SetFiniteWidth(BetaOverMu(W), Upl, toff);
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
   for (int s = 0; s < nsteps; s++)
   {
      real_t dts = dt;
      ode.Step(U, V, t, dts);
      if (s % 100 == 99 || s == nsteps - 1) { vmax = std::max(vmax, GMax(V.Normlinf(), c)); }
   }
   return vmax / (gam * Cs());
}

void TestF2(int nsteps)
{
   const real_t vb = RestRun(nsteps, true);
   const real_t vp = RestRun(nsteps, false);
   if (g_rank == 0)
   {
      std::printf("  F2 p=2, %d steps: max |V| / (gamma c_s) band %.3e, no band (rounding "
                  "floor) %.3e\n", nsteps, vb, vp);
   }
   Check(vb <= std::max(1e-12, 2.0 * vp), "F2: the QD state with the 2.5D term stays at rest "
         "in the band run, <= max(1e-12, 2 x the no-band floor)", vb);
}

// ------------------------------------------------------------------ F3 / F4
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

// mode: 0 no band, 1 the band with the stretched 2.5D term, 2 the band with
// the term left plain (the elastic operator's M_beta alone)
PulseRun RunPulse(real_t Hbox, real_t h, int p, int mode, real_t W,
                  real_t L, real_t fa, real_t Hin, real_t dt, real_t T, real_t snap_dt)
{
   const int n = static_cast<int>(std::lround(2 * Hbox / h));
   Mesh sm = Box(2, n, Hbox);
   ParMesh pm(MPI_COMM_WORLD, sm);
   DFn g0 = [](const Vector &, real_t, Vector &u) { u = 0.0; };
   auto op = Op(pm, p, g0);
   Vector Upl;
   Substrate(*op, 0.0, Upl);
   std::shared_ptr<const PmlBand> band;
   if (mode != 0) { band = MakeBand(pm, L, fa); op->SetPml(band); }
   if (mode == 2) { op->Elasticity().SetFiniteWidth(BetaOverMu(W), Upl, nullptr); }
   else { op->SetFiniteWidth(BetaOverMu(W), Upl, nullptr); }
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

real_t BandDt(real_t Hb, real_t h, int p, real_t L, real_t fa)
{
   const int n = static_cast<int>(std::lround(2 * Hb / h));
   Mesh sm = Box(2, n, Hb);
   ParMesh pm(MPI_COMM_WORLD, sm);
   DFn g0 = [](const Vector &, real_t, Vector &u) { u = 0.0; };
   auto opb = Op(pm, p, g0);
   opb->SetPml(MakeBand(pm, L, fa));
   return 0.5 * 2.0 / std::sqrt(OmegaMax2(*opb));
}

void TestF3(real_t W, real_t T)
{
   const real_t Hin = 10e3, L = 5e3, fa = 0.2, Hb = Hin + L, h = 500.0, sdt = 0.1;
   const int p = 2;
   const real_t dt = BandDt(Hb, h, p, L, fa);
   // no return by T; a whole number of cells beyond the inner box (the same
   // centroids in every run)
   const real_t Href = Hin + h * std::ceil((0.5 * Cp() * T + 2e3) / h);
   PulseRun ref = RunPulse(Href, h, p, 0, W, 0, 0, Hin, dt, T, sdt);
   PulseRun bnd = RunPulse(Hb, h, p, 1, W, L, fa, Hin, dt, T, sdt);
   PulseRun plain = RunPulse(Hin, h, p, 0, W, 0, 0, Hin, dt, T, sdt);
   PulseRun uns = RunPulse(Hb, h, p, 2, W, L, fa, Hin, dt, T, sdt);
   const real_t mb = Misfit(bnd, ref), mp = Misfit(plain, ref), mu = Misfit(uns, ref);
   if (g_rank == 0 && std::getenv("PML_F3_TRACE"))
   {
      std::printf("    dt %.4e Href %.0f finite ref %d bnd %d plain %d uns %d; snaps %zu %zu; "
                  "points %zu %zu; e0 %.3e %.3e\n", dt, Href, ref.finite, bnd.finite,
                  plain.finite, uns.finite, ref.snaps.size(), bnd.snaps.size(),
                  ref.snaps.empty() ? 0 : ref.snaps[0].size(),
                  bnd.snaps.empty() ? 0 : bnd.snaps[0].size(),
                  ref.energy.empty() ? 0.0 : ref.energy[0], bnd.energy.empty() ? 0.0 : bnd.energy[0]);
   }
   const real_t fc = std::sqrt(BetaOverMu(W) * kMu / kRho) / (2 * M_PI);
   if (g_rank == 0)
   {
      std::printf("  F3 W = %.0f km (cutoff %.3f Hz), t <= %.0f s: inner-box misfit vs the "
                  "large box: band %.4e, bare walls %.4e (ratio %.4f); band with the 2.5D "
                  "term left unstretched %.4e\n", W / 1e3, fc, T, mb, mp, mb / mp, mu);
   }
   std::ostringstream s;
   s << "F3 W = " << W / 1e3 << " km: the band's reflection < 5% of the bare walls'";
   Check(bnd.finite && mb < 0.05 * mp, s.str(), mb / mp);
}

void TestF4(real_t W, bool decay_only)
{
   const real_t Hin = 10e3, L = 5e3, fa = 0.2, Hb = Hin + L, h = 1000.0, sdt = 0.5;
   const int p = 2;
   const real_t T = 10.0 * 2 * Hb / Cs();
   const real_t dt = BandDt(Hb, h, p, L, fa);
   PulseRun bnd = RunPulse(Hb, h, p, 1, W, L, fa, Hin, dt, T, sdt);
   const size_t n = bnd.energy.size();
   real_t peak = 0.0, q2 = 0.0, q4 = 0.0;
   for (real_t e : bnd.energy) { peak = std::max(peak, e); }
   for (size_t k = n / 4; k < n / 2; k++) { q2 = std::max(q2, bnd.energy[k]); }
   for (size_t k = n / 2; k < n; k++) { q4 = std::max(q4, bnd.energy[k]); }
   if (g_rank == 0)
   {
      std::printf("  F4 W = %.0f km, %.0f s, %ld steps: inner sum|u|^2 peak %.4e, max over "
                  "the second quarter %.4e, over the second half %.4e\n", W / 1e3, T,
                  std::lround(T / dt), peak, q2, q4);
   }
   std::ostringstream s;
   if (decay_only)
   {
      // waves near the cutoff leave slowly (group speed -> 0), so the bound
      // of the plain test does not apply: the energy must keep falling
      s << "F4 long run at W = " << W / 1e3 << " km (cutoff inside the pulse's band): "
           "finite, the second half below the second quarter (no growth)";
      Check(bnd.finite && q4 < q2, s.str(), q4 / q2);
   }
   else
   {
      s << "F4 long run at W = " << W / 1e3 << " km: finite, second-half inner energy < "
           "1e-3 of the peak";
      Check(bnd.finite && q4 < 1e-3 * peak, s.str(), q4 / peak);
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
   args.AddOption(&only, "-o", "--only", "Comma-separated parts (F0,F1,F2,F3,F4).");
   args.Parse();
   if (!args.Good()) { if (g_rank == 0) { args.PrintUsage(std::cout); } return 1; }
   {
      std::stringstream ss(only);
      std::string it;
      while (std::getline(ss, it, ',')) { if (!it.empty()) { g_only.insert(it); } }
   }
   if (Want("F0"))
   {
      if (g_rank == 0) { std::cout << "F0: static identity with the 2.5D term\n"; }
      TestF0(1); TestF0(2);
   }
   if (Want("F1"))
   {
      if (g_rank == 0) { std::cout << "F1: zero-frequency limit\n"; }
      TestF1();
   }
   if (Want("F2"))
   {
      if (g_rank == 0) { std::cout << "F2: rest stays at rest\n"; }
      TestF2(10000);
   }
   if (Want("F3"))
   {
      if (g_rank == 0) { std::cout << "F3: pulse into the band\n"; }
      const char *tenv = std::getenv("PML_F3_T");   // a shorter run (debugging)
      const real_t T = tenv ? std::atof(tenv) : 7.0;
      TestF3(15e3, T); TestF3(3e3, T);
   }
   if (Want("F4"))
   {
      if (g_rank == 0) { std::cout << "F4: long run\n"; }
      TestF4(15e3, false); TestF4(3e3, true);
   }
   if (g_rank == 0)
   {
      if (num_failed) { std::cout << "FAILED: " << num_failed << " of " << num_tests << "\n"; }
      else { std::cout << "ALL PASS (" << num_tests << " checks)\n"; }
   }
   return num_failed ? 1 : 0;
}
