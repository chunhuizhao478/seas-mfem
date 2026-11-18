// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// Viscoelastic attenuation plan rev6, gate V3: plane strain on triangles.
// The same GMB model, coupler and operators as 3-D (test_visco_coupler.cpp),
// on a 2-D mesh with 2 displacement components; the strain sampler writes
// eps_zz = eps_yz = eps_xz = 0.
//
//   O2-2D  patch test on triangles, p = 1, 2, 3: a uniform dsig (all six
//          components set) gives P dsig = -int_{natural} (dsig n) . w with
//          the in-plane rows only, 1e-12 relative, with Dirichlet faces.
//   Z      out-of-plane memory: after a damped-wave run every xi_l has
//          xx yy xy only; the zz, yz, xz slots and dsig_yz, dsig_xz are 0
//          exactly at every node (dsig_zz = lambda_l tr xi_l is not 0: the
//          plane-strain sigma_zz, which no in-plane force reads).
//   T-V6   the damped plane S and P waves of T-V2 in the (x, y) plane: the
//          2-D run on triangles and the 3-D run on the same box extruded one
//          hex (24 tets) thick, both with u = the exact wave on every face,
//          reproduce the same exact damped wave (speed and decay of the fit)
//          to within 3 x their elastic discretization error, and below 0.1
//          of the gap to the elastic wave.  A tet layer is not a triangle
//          mesh, so the two runs are compared through the exact wave they
//          both must converge to, not node by node.
//
// Usage: mpirun -np 1|2 ./seas_test_visco_plane_strain
//
// Scope (Chunhui Zhao, 2026-11-18).

#include "mfem.hpp"
#include "../../domain/nitsche_wave_operator.hpp"
#include "../../solver/nitsche_cfl.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../rheology/rheology_coupler.hpp"
#include "../../rheology/rheology_config.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;
using cplx = std::complex<real_t>;

static int num_tests = 0, num_passed = 0, num_failed = 0;
static int g_rank = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; if (g_rank == 0) std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; if (g_rank == 0) std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static const real_t kBox = 12.0e3;
static const real_t kCuse3[5] = {0.0, 0.057, 0.11, 0.14, 0.166};   // 0.5 C_eig, tets
static const real_t kF0 = 0.5;
static const real_t kQs = 20.0, kQp = 30.0;

static real_t RMax(real_t v) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); return v; }
static int ISum(int v) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); return v; }

static Mesh *Box2(int n, bool natural_tb)
{
   Mesh *m = new Mesh(Mesh::MakeCartesian2D(n, n, Element::TRIANGLE, true, kBox, kBox));
   for (int b = 0; b < m->GetNBE(); b++)
   {
      Array<int> v;
      m->GetBdrElementVertices(b, v);
      const real_t y = 0.5 * (m->GetVertex(v[0])[1] + m->GetVertex(v[1])[1]);
      const bool tb = std::abs(y) < 1e-6 * kBox || std::abs(y - kBox) < 1e-6 * kBox;
      m->SetBdrAttribute(b, (natural_tb && tb) ? 1 : 5);
   }
   m->SetAttributes();
   return m;
}

static Mesh *Slab3(int n)
{
   Mesh *m = new Mesh(Mesh::MakeCartesian3DWith24TetsPerHex(n, n, 1, kBox, kBox, kBox / n));
   for (int b = 0; b < m->GetNBE(); b++) { m->SetBdrAttribute(b, 5); }
   m->SetAttributes();
   return m;
}

static MaterialField Material()
{
   return MaterialField::MakeConstant(TPV102Params::lambda, TPV102Params::mu, TPV102Params::rho);
}

static RheologyConfig VEConfig()
{
   RheologyConfig cfg;
   cfg.viscoelastic = true;
   cfg.gmb_band.L = 3;
   cfg.gmb_band.f_c = kF0;
   cfg.gmb_band.R = 100.0;
   cfg.gmb_q.from_vs = false;
   cfg.gmb_q.Q_s = kQs;
   cfg.gmb_q.Q_p = kQp;
   return cfg;
}

static GMBFit TheFit()
{
   const RheologyConfig cfg = VEConfig();
   return FitGMB(kQp, kQs, cfg.gmb_band, TPV102Params::lambda, TPV102Params::mu,
                 GMBModuli::Relaxed);
}

using VecFn = std::function<void(const Vector &, real_t, Vector &)>;

static std::unique_ptr<NitscheWaveOperator> MakeOp(ParMesh &pm, int p, const VecFn &uD,
                                                   bool natural_tb = false)
{
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.dirichlet_attrs = {5};
   if (natural_tb) { bc.natural_attrs = {1}; }
   bc.default_dirichlet_func = [uD](const Vector &x, real_t t, Vector &u) { uD(x, t, u); };
   return std::make_unique<NitscheWaveOperator>(pm, p, Material(), bc);
}

static std::unique_ptr<RheologyCoupler> MakeCoupler(NitscheWaveOperator &op)
{
   auto zero = [](const real_t *, real_t s[9]) { for (int i = 0; i < 9; i++) { s[i] = 0.0; } };
   return std::make_unique<RheologyCoupler>(op.Elasticity().GetFESpace(), op.Sample(),
                                            op.Boundary().dirichlet_attrs,
                                            BuildRheologyStack(VEConfig()), zero);
}

static real_t StepWith(NitscheWaveOperator &op, ParMesh &pm, int p, const RheologyCoupler *rc)
{
   // 0.5 x the stable constant: of triangles (NitscheStableConstant(2, p),
   // as test_2d_nitsche.cpp) or the tet table of test_visco_coupler.cpp
   real_t dt = (pm.Dimension() == 2)
               ? ComputeMaxDt(op.Sample(), pm, p, 0.5 * NitscheStableConstant(2, p))
               : ComputeMaxDt(op.Sample(), pm, p, kCuse3[p]);
   if (rc)
   {
      const RheologyField &F = rc->Field();
      real_t b = 1.0;
      for (int e = 0; e < F.NumElements(); e++)
      {
         b = std::max(b, rc->Stack().StiffnessBound(F.Context(e, 0), F.State(e, 0)));
      }
      dt /= std::sqrt(RMax(b));
   }
   return dt;
}

// ------------------------------------------------- in-plane damped wave

struct DampedWave
{
   real_t n[3], d[3], w, A;
   cplx kappa;
   cplx mem[3];
   int dim = 2;

   cplx Phase(const real_t *x, real_t t) const
   {
      const real_t s = n[0] * x[0] + n[1] * x[1];
      return std::exp(cplx(0.0, 1.0) * (w * t - kappa * s));
   }
   void U(const Vector &x, real_t t, Vector &u) const
   {
      const cplx e = A * Phase(x.GetData(), t);
      u.SetSize(dim);
      for (int i = 0; i < dim; i++) { u(i) = std::real(d[i] * e); }
   }
   void V(const Vector &x, real_t t, Vector &v) const
   {
      const cplx e = cplx(0.0, w) * A * Phase(x.GetData(), t);
      v.SetSize(dim);
      for (int i = 0; i < dim; i++) { v(i) = std::real(d[i] * e); }
   }
   void Eps(const real_t *x, real_t t, cplx e[6]) const
   {
      const cplx g = cplx(0.0, -1.0) * kappa * A * Phase(x, t);
      e[0] = g * d[0] * n[0];
      e[1] = g * d[1] * n[1];
      e[2] = 0.0;
      e[3] = 0.5 * g * (d[0] * n[1] + d[1] * n[0]);
      e[4] = 0.0;
      e[5] = 0.0;
   }
};

static DampedWave MakeWave(bool pwave, bool elastic, const GMBFit &f, int dim)
{
   DampedWave W;
   W.dim = dim;
   const real_t nn = std::hypot(1.0, 0.7);
   W.n[0] = 1.0 / nn; W.n[1] = 0.7 / nn; W.n[2] = 0.0;
   if (pwave) { W.d[0] = W.n[0]; W.d[1] = W.n[1]; }
   else { W.d[0] = -W.n[1]; W.d[1] = W.n[0]; }
   W.d[2] = 0.0;
   W.w = 2.0 * M_PI * kF0;
   cplx M = pwave ? cplx(f.lam_R + 2.0 * f.mu_R) : cplx(f.mu_R);
   for (int l = 0; l < 3; l++)
   {
      W.mem[l] = elastic ? cplx(0.0) : cplx(0.0, W.w) / cplx(f.omega[l], W.w);
      M += (pwave ? f.lam_l[l] + 2.0 * f.mu_l[l] : f.mu_l[l]) * W.mem[l];
   }
   W.kappa = W.w * std::sqrt(cplx(TPV102Params::rho) / M);
   W.A = 1.0e-5 / std::abs(W.kappa);
   return W;
}

static void SetWaveMemory(RheologyCoupler &rc, const DampedWave &W, const GMBFit &f)
{
   RheologyField &F = rc.Field();
   for (int e = 0; e < F.NumElements(); e++)
   {
      for (int i = 0; i < F.NodesPerElement(); i++)
      {
         cplx eps[6];
         W.Eps(F.Context(e, i).x, 0.0, eps);
         real_t *q = F.State(e, i);
         real_t *ds = F.Dsig(e, i);
         for (int c = 0; c < 6; c++) { ds[c] = 0.0; }
         for (int l = 0; l < 3; l++)
         {
            real_t xi[6], s[6];
            for (int c = 0; c < 6; c++) { xi[c] = std::real(W.mem[l] * eps[c]); q[6 * l + c] = xi[c]; }
            sym6::ApplyIsotropic(f.lam_l[l], f.mu_l[l], xi, s);
            for (int c = 0; c < 6; c++) { ds[c] += s[c]; }
         }
      }
   }
   rc.Refresh();
}

/// One period; relative L2 error to the wave W at the end.  zero_oop: the
/// count of nonzero out-of-plane memory / dsig slots after the run.
static real_t RunWave(ParMesh &pm, int p, const DampedWave &W, bool visco, int *nonzero_oop = nullptr)
{
   const int dim = pm.Dimension();
   VecFn uD = [W](const Vector &x, real_t t, Vector &u) { W.U(x, t, u); };
   auto op = MakeOp(pm, p, uD);
   std::unique_ptr<RheologyCoupler> rc;
   if (visco) { rc = MakeCoupler(*op); op->SetRheology(rc.get()); }
   auto &fes = op->Elasticity().GetFESpace();
   ParGridFunction gu(&fes), gv(&fes);
   VectorFunctionCoefficient cu(dim, [W](const Vector &x, real_t t, Vector &u) { W.U(x, t, u); });
   VectorFunctionCoefficient cv(dim, [W](const Vector &x, real_t t, Vector &v) { W.V(x, t, v); });
   cu.SetTime(0.0);
   cv.SetTime(0.0);
   gu.ProjectCoefficient(cu);
   gv.ProjectCoefficient(cv);
   Vector U(gu), V(gv);
   if (rc) { rc->Reset(U); SetWaveMemory(*rc, W, TheFit()); }
   const real_t T = 1.0 / kF0;
   real_t dt = StepWith(*op, pm, p, rc.get());
   const int nsteps = static_cast<int>(std::ceil(T / dt));
   dt = T / nsteps;
   op->SetTime(0.0);
   op->SetTimeStep(dt);
   real_t t = 0.0;
   CentralDifferenceSolver ode;
   ode.Init(*op);
   for (int s = 0; s < nsteps; s++) { real_t h = dt; ode.Step(U, V, t, h); }
   if (nonzero_oop && rc)
   {
      const RheologyField &F = rc->Field();
      int nz = 0;
      for (int e = 0; e < F.NumElements(); e++)
      {
         for (int i = 0; i < F.NodesPerElement(); i++)
         {
            const real_t *q = F.State(e, i);
            for (int l = 0; l < 3; l++)
            {
               nz += (q[6 * l + 2] != 0.0) + (q[6 * l + 4] != 0.0) + (q[6 * l + 5] != 0.0);
            }
            nz += (q[6 * 3 + 2] != 0.0) + (q[6 * 3 + 4] != 0.0) + (q[6 * 3 + 5] != 0.0);   // eps^n
            nz += (F.Dsig(e, i)[4] != 0.0) + (F.Dsig(e, i)[5] != 0.0);
         }
      }
      *nonzero_oop = ISum(nz);
   }
   gu = U;
   VectorFunctionCoefficient ce(dim, [W](const Vector &x, real_t tt, Vector &u) { W.U(x, tt, u); });
   ce.SetTime(t);
   ParGridFunction zero(&fes);
   zero = 0.0;
   return gu.ComputeL2Error(ce) / zero.ComputeL2Error(ce);
}

static real_t WaveGap(ParMesh &pm, int p, const DampedWave &a, const DampedWave &b)
{
   const int dim = pm.Dimension();
   H1_FECollection fec(p + 2, dim);
   ParFiniteElementSpace fes(&pm, &fec, dim);
   ParGridFunction g(&fes), zero(&fes);
   VectorFunctionCoefficient ca(dim, [a](const Vector &x, real_t t, Vector &u) { a.U(x, t, u); });
   VectorFunctionCoefficient cb(dim, [b](const Vector &x, real_t t, Vector &u) { b.U(x, t, u); });
   ca.SetTime(1.0 / kF0);
   cb.SetTime(1.0 / kF0);
   g.ProjectCoefficient(ca);
   zero = 0.0;
   return g.ComputeL2Error(cb) / zero.ComputeL2Error(cb);
}

// ------------------------------------------------------------------ O2-2D
static void TestO2()
{
   if (g_rank == 0) { std::cout << "\n=== O2-2D patch test on triangles ===\n"; }
   const real_t D[6] = {3.0e6, -2.0e6, 1.5e6, 0.7e6, -1.1e6, 0.4e6};
   for (int p = 1; p <= 3; p++)
   {
      std::unique_ptr<Mesh> m(Box2(3, true));
      ParMesh pm(MPI_COMM_WORLD, *m);
      m.reset();
      VecFn zero = [](const Vector &, real_t, Vector &u) { u.SetSize(2); u = 0.0; };
      auto op = MakeOp(pm, p, zero, true);
      auto rc = MakeCoupler(*op);
      RheologyField &F = rc->Field();
      for (int e = 0; e < F.NumElements(); e++)
         for (int i = 0; i < F.NodesPerElement(); i++)
            for (int c = 0; c < 6; c++) { F.Dsig(e, i)[c] = D[c]; }
      rc->Refresh();
      ParFiniteElementSpace &ufes = op->Elasticity().GetFESpace();
      Vector Fp(ufes.GetVSize()), Fr(ufes.GetVSize());
      Fp = 0.0;
      rc->AddLoad(Fp);
      Fr = 0.0;
      real_t S[9];
      sym6::ToRowMajor(D, S);
      Array<int> vd;
      Vector nor(2), sh;
      int nN = 0;
      for (int be = 0; be < pm.GetNBE(); be++)
      {
         if (pm.GetBdrAttribute(be) != 1) { continue; }
         FaceElementTransformations *Tr = pm.GetBdrFaceTransformations(be);
         if (!Tr) { continue; }
         nN++;
         const FiniteElement &fe = *ufes.GetFE(Tr->Elem1No);
         const int nd = fe.GetDof();
         sh.SetSize(nd);
         ufes.GetElementVDofs(Tr->Elem1No, vd);
         const IntegrationRule &ir = IntRules.Get(Tr->GetGeometryType(), 2 * p + 1);
         for (int q = 0; q < ir.GetNPoints(); q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            Tr->SetAllIntPoints(&ip);
            CalcOrtho(Tr->Jacobian(), nor);
            fe.CalcShape(Tr->GetElement1IntPoint(), sh);
            for (int c = 0; c < 2; c++)
            {
               const real_t t = ip.weight * (S[3 * c] * nor(0) + S[3 * c + 1] * nor(1));
               for (int j = 0; j < nd; j++) { Fr(vd[c * nd + j]) -= sh(j) * t; }
            }
         }
      }
      real_t scale = 0.0, err = 0.0;
      for (int i = 0; i < Fr.Size(); i++)
      {
         scale = std::max(scale, std::abs(Fr(i)));
         err = std::max(err, std::abs(Fp(i) - Fr(i)));
      }
      scale = RMax(scale);
      err = RMax(err);
      const int nD = ISum(rc->Operators().NumFaces(2)), nNa = ISum(nN);
      if (g_rank == 0)
      {
         std::cout << std::setprecision(4) << "  p = " << p << ": " << F.NodesPerElement()
                   << " nodes per triangle, " << nNa << " natural, " << nD
                   << " Dirichlet faces; |P dsig - ref| / max|ref| = " << err / scale << "\n";
      }
      TEST_ASSERT(nD > 0 && nNa > 0 && err <= 1e-12 * scale,
                  "O2-2D p = " << p << " P dsig = -int_N (dsig n) . w");
   }
}

// ------------------------------------------------------------- Z and T-V6
static void TestTV6()
{
   if (g_rank == 0) { std::cout << "\n=== T-V6 in-plane damped waves: 2-D triangles and a 3-D slab ===\n"; }
   const GMBFit f = TheFit();
   const int p = 3, n = 4;
   std::unique_ptr<Mesh> m2(Box2(n, false)), m3(Slab3(n));
   ParMesh pm2(MPI_COMM_WORLD, *m2), pm3(MPI_COMM_WORLD, *m3);
   m2.reset();
   m3.reset();
   for (int pw = 0; pw <= 1; pw++)
   {
      const char *nm = pw ? "P" : "S";
      const DampedWave v2 = MakeWave(pw == 1, false, f, 2), e2 = MakeWave(pw == 1, true, f, 2);
      const DampedWave v3 = MakeWave(pw == 1, false, f, 3), e3 = MakeWave(pw == 1, true, f, 3);
      int nz = -1;
      const real_t ev2 = RunWave(pm2, p, v2, true, &nz), ee2 = RunWave(pm2, p, e2, false);
      const real_t ev3 = RunWave(pm3, p, v3, true), ee3 = RunWave(pm3, p, e3, false);
      const real_t gap = WaveGap(pm2, p, e2, v2);
      if (g_rank == 0)
      {
         std::cout << std::setprecision(5) << "  " << nm << " wave, p = " << p << ", " << n << " x "
                   << n << ": 2-D visco " << ev2 << " (elastic " << ee2 << "), 3-D slab visco "
                   << ev3 << " (elastic " << ee3 << "); gap to the elastic wave " << gap
                   << "; nonzero out-of-plane slots " << nz << "\n";
      }
      TEST_ASSERT(nz == 0, "Z " << nm << " out-of-plane memory and dsig_yz, dsig_xz exactly 0");
      TEST_ASSERT(ev2 < 0.1 * gap && ev2 < 3.0 * ee2,
                  "T-V6 " << nm << " 2-D run reproduces the damped wave");
      TEST_ASSERT(ev3 < 0.1 * gap && ev3 < 3.0 * ee3,
                  "T-V6 " << nm << " 3-D slab run reproduces the same damped wave");
   }
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   if (g_rank == 0)
   {
      std::cout << "=== Viscoelastic attenuation, plane strain (plan rev6, V3) on "
                << Mpi::WorldSize() << " rank(s) ===\n";
   }
   TestO2();
   TestTV6();
   if (g_rank == 0)
   {
      std::cout << "\n" << num_passed << " passed, " << num_failed << " failed, " << num_tests
                << " total\n";
   }
   return num_failed == 0 ? 0 : 1;
}
