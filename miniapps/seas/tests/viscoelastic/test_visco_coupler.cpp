// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// Viscoelastic attenuation plan rev6, gates V2 (3-D) and V3 (plane strain):
// the GMB model of rheology/viscoelastic.hpp in the rheology coupler and the
// Nitsche wave operator, global stepping, no fault.
//
// Mesh: a box of side 12 km (MakeCartesian3DWith24TetsPerHex, or
// MakeCartesian2D triangles), every boundary face Dirichlet (attr 5).
// Material as TPV102 (V_S 3464 m/s, V_P 6000 m/s, rho 2670 kg/m^3); band
// f_c = 0.5 Hz, R = 100, L = 3; fixed Q_S = 20, Q_P = 30 (strong, so that
// the attenuation is far above the discretization error).
//
//   T-V2 damped plane S and P waves, p = 3: u = Re[A d exp(i(w t - kappa n.x))]
//        with real w = 2 pi 0.5 Hz and the complex wavenumber
//        kappa = w sqrt(rho / M(w)) of (2.18a), M(w) = M_R + sum_l M_l i w /
//        (i w + w_l) (the GMB modulus of the fit, relaxed convention D1).
//        Started from the exact u, du/dt and the exact memory xi_l =
//        Re[i w / (i w + w_l) eps], the Dirichlet faces carry the exact u.
//        After one period the L2 error to it must be below 0.1 of the gap
//        between it and the elastic (relaxed-moduli) wave, and below 3 times
//        the error of the elastic run against its own exact wave: the speed
//        (2.18b) and the decay (2.18c) are those of the fit.  On 2 x 2 x 2
//        and 4 x 4 x 4 hexes, so the error also falls with h.
//   T-V3 energy: a Gaussian bump in a closed box (u = 0 on every face).
//        E_h = 1/2 v.M v + 1/2 u.K u - dt^2/8 (K u).M^-1 (K u) + 1/2 int
//        sum_l xi_l : C_l : xi_l, the energy (2.34) with the O(dt^2) term
//        that makes it the exact invariant of the central difference scheme
//        without attenuation (checked: the elastic run keeps it to 1e-10).  With
//        attenuation it never rises by more than 1e-10 E_h(0), it falls,
//        and the drop is int D dt to 5 %.
//   T-V4 rest round trip: Reset(U) and two steps of the coupler with the same
//        U: Delta sigma = 0 exactly at every node and P Delta sigma = 0.
//   T-V5 two dynamic phases in a row: phase 2 from Reset(U_1) of the coupler
//        that ran phase 1 is bit for bit phase 2 with a new coupler.
//   T-V6 (plane strain, gate V3): the 2-D S wave on triangles equals, to
//        1e-12, the 3-D run on the same box extruded one element thick ...
//        see test_visco_plane_strain.cpp.
//
// Usage: mpirun -np 1|2 ./seas_test_visco_coupler
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
#include <cstring>
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
static const real_t kCuse[5] = {0.0, 0.057, 0.11, 0.14, 0.166};   // 0.5 C_eig
static const real_t kF0 = 0.5;                                    // [Hz]
static const real_t kQs = 20.0, kQp = 30.0;

static real_t RSum(real_t v) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); return v; }
static real_t RMax(real_t v) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); return v; }
static int ISum(int v) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); return v; }

// ---------------------------------------------------------------- fixture

static Mesh *Box3(int n)
{
   Mesh *m = new Mesh(Mesh::MakeCartesian3DWith24TetsPerHex(n, n, n, kBox, kBox, kBox));
   for (int b = 0; b < m->GetNBE(); b++) { m->SetBdrAttribute(b, 5); }
   m->SetAttributes();
   return m;
}

static MaterialField Material()
{
   return MaterialField::MakeConstant(TPV102Params::lambda, TPV102Params::mu, TPV102Params::rho);
}

static RheologyConfig VEConfig(real_t Qs, real_t Qp)
{
   RheologyConfig cfg;
   cfg.viscoelastic = true;
   cfg.gmb_band.L = 3;
   cfg.gmb_band.f_c = kF0;
   cfg.gmb_band.R = 100.0;
   cfg.gmb_q.from_vs = false;
   cfg.gmb_q.Q_s = Qs;
   cfg.gmb_q.Q_p = Qp;
   return cfg;
}

static GMBFit TheFit()
{
   const RheologyConfig cfg = VEConfig(kQs, kQp);
   return FitGMB(kQp, kQs, cfg.gmb_band, TPV102Params::lambda, TPV102Params::mu,
                 GMBModuli::Relaxed);
}

using VecFn = std::function<void(const Vector &, real_t, Vector &)>;

static std::unique_ptr<NitscheWaveOperator> MakeOp(ParMesh &pm, int p, const VecFn &uD)
{
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.dirichlet_attrs = {5};
   bc.default_dirichlet_func = [uD](const Vector &x, real_t t, Vector &u) { uD(x, t, u); };
   return std::make_unique<NitscheWaveOperator>(pm, p, Material(), bc);
}

static std::unique_ptr<RheologyCoupler> MakeCoupler(NitscheWaveOperator &op, real_t Qs, real_t Qp)
{
   auto zero = [](const real_t *, real_t s[9]) { for (int i = 0; i < 9; i++) { s[i] = 0.0; } };
   return std::make_unique<RheologyCoupler>(op.Elasticity().GetFESpace(), op.Sample(),
                                            op.Boundary().dirichlet_attrs,
                                            BuildRheologyStack(VEConfig(Qs, Qp)), zero);
}

/// The global time step: 0.5 C_eig times 1 / sqrt(StiffnessBound) (plan 5.9).
static real_t StepWith(NitscheWaveOperator &op, ParMesh &pm, int p, const RheologyCoupler *rc)
{
   real_t dt = ComputeMaxDt(op.Sample(), pm, p, kCuse[p]);
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

// ------------------------------------------------------- the damped plane wave

struct DampedWave
{
   real_t n[3], d[3], w, A;
   cplx kappa;                     ///< complex wavenumber, Im < 0 (decay along n)
   cplx mem[3];                    ///< i w / (i w + w_l)
   bool elastic = false;

   cplx Phase(const real_t *x, real_t t) const
   {
      const real_t s = n[0] * x[0] + n[1] * x[1] + n[2] * x[2];
      return std::exp(cplx(0.0, 1.0) * (w * t - kappa * s));
   }
   void U(const Vector &x, real_t t, Vector &u) const
   {
      const cplx e = A * Phase(x.GetData(), t);
      u.SetSize(3);
      for (int i = 0; i < 3; i++) { u(i) = std::real(d[i] * e); }
   }
   void V(const Vector &x, real_t t, Vector &v) const
   {
      const cplx e = cplx(0.0, w) * A * Phase(x.GetData(), t);
      v.SetSize(3);
      for (int i = 0; i < 3; i++) { v(i) = std::real(d[i] * e); }
   }
   /// The complex strain (xx yy zz xy yz xz, tensor shear) at (x, t).
   void Eps(const real_t *x, real_t t, cplx e[6]) const
   {
      const cplx g = cplx(0.0, -1.0) * kappa * A * Phase(x, t);   // d/dx_j = -i kappa n_j
      cplx H[3][3];
      for (int i = 0; i < 3; i++)
         for (int j = 0; j < 3; j++) { H[i][j] = g * d[i] * n[j]; }
      e[0] = H[0][0]; e[1] = H[1][1]; e[2] = H[2][2];
      e[3] = 0.5 * (H[0][1] + H[1][0]);
      e[4] = 0.5 * (H[1][2] + H[2][1]);
      e[5] = 0.5 * (H[0][2] + H[2][0]);
   }
};

static DampedWave MakeWave(bool pwave, bool elastic, const GMBFit &f)
{
   DampedWave W;
   real_t n[3] = {1.0, 0.7, 0.4};
   const real_t nn = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
   for (int i = 0; i < 3; i++) { W.n[i] = n[i] / nn; }
   if (pwave) { for (int i = 0; i < 3; i++) { W.d[i] = W.n[i]; } }
   else
   {
      // d = n x e_z / |n x e_z|
      real_t c[3] = {W.n[1], -W.n[0], 0.0};
      const real_t cn = std::hypot(c[0], c[1]);
      for (int i = 0; i < 3; i++) { W.d[i] = c[i] / cn; }
   }
   W.w = 2.0 * M_PI * kF0;
   W.elastic = elastic;
   cplx M = pwave ? cplx(f.lam_R + 2.0 * f.mu_R) : cplx(f.mu_R);
   for (int l = 0; l < 3; l++)
   {
      W.mem[l] = elastic ? cplx(0.0) : cplx(0.0, W.w) / cplx(f.omega[l], W.w);
      const real_t Ml = pwave ? f.lam_l[l] + 2.0 * f.mu_l[l] : f.mu_l[l];
      M += Ml * W.mem[l];
   }
   W.kappa = W.w * std::sqrt(cplx(TPV102Params::rho) / M);
   W.A = 1.0e-5 / std::abs(W.kappa);          // strain amplitude 1e-5
   return W;
}

/// xi_l and dsig at every node from the exact wave at t = 0 (after Reset).
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

/// One period of the wave W on the mesh pm; with a coupler when visco.
/// Returns the relative L2 error to `exact` at the end.
static real_t RunWave(ParMesh &pm, int p, const DampedWave &W, bool visco,
                      const DampedWave &exact, real_t *dt_used = nullptr)
{
   VecFn uD = [W](const Vector &x, real_t t, Vector &u) { W.U(x, t, u); };
   auto op = MakeOp(pm, p, uD);
   std::unique_ptr<RheologyCoupler> rc;
   if (visco) { rc = MakeCoupler(*op, kQs, kQp); op->SetRheology(rc.get()); }
   auto &fes = op->Elasticity().GetFESpace();
   ParGridFunction gu(&fes), gv(&fes);
   VectorFunctionCoefficient cu(3, [W](const Vector &x, real_t t, Vector &u) { W.U(x, t, u); });
   VectorFunctionCoefficient cv(3, [W](const Vector &x, real_t t, Vector &v) { W.V(x, t, v); });
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
   if (dt_used) { *dt_used = dt; }
   op->SetTime(0.0);
   op->SetTimeStep(dt);
   real_t t = 0.0;
   CentralDifferenceSolver ode;
   ode.Init(*op);
   for (int s = 0; s < nsteps; s++) { real_t h = dt; ode.Step(U, V, t, h); }
   gu = U;
   VectorFunctionCoefficient ce(3, [exact](const Vector &x, real_t tt, Vector &u) { exact.U(x, tt, u); });
   ce.SetTime(t);
   ParGridFunction zero(&fes);
   zero = 0.0;
   return gu.ComputeL2Error(ce) / zero.ComputeL2Error(ce);
}

/// || u_a(T) - u_b(T) || / || u_b(T) || of two exact waves (on mesh pm).
static real_t WaveGap(ParMesh &pm, int p, const DampedWave &a, const DampedWave &b)
{
   H1_FECollection fec(p + 2, 3);
   ParFiniteElementSpace fes(&pm, &fec, 3);
   ParGridFunction g(&fes), zero(&fes);
   VectorFunctionCoefficient ca(3, [a](const Vector &x, real_t t, Vector &u) { a.U(x, t, u); });
   VectorFunctionCoefficient cb(3, [b](const Vector &x, real_t t, Vector &u) { b.U(x, t, u); });
   ca.SetTime(1.0 / kF0);
   cb.SetTime(1.0 / kF0);
   g.ProjectCoefficient(ca);
   zero = 0.0;
   return g.ComputeL2Error(cb) / zero.ComputeL2Error(cb);
}

static void TestTV2()
{
   if (g_rank == 0) { std::cout << "\n=== T-V2 damped plane waves (speed and decay of the fit) ===\n"; }
   const GMBFit f = TheFit();
   TEST_ASSERT(f.error.empty(), "T-V2 fit Q_S = 20, Q_P = 30 accepted (" << f.error << ")");
   const int p = 3;
   for (int pw = 0; pw <= 1; pw++)
   {
      const DampedWave Wv = MakeWave(pw == 1, false, f), We = MakeWave(pw == 1, true, f);
      const real_t c_ph = Wv.w / std::real(Wv.kappa), c_R = We.w / std::real(We.kappa);
      const real_t decay = -std::imag(Wv.kappa);
      real_t ev[2], ee[2];
      real_t gap = 0.0;
      for (int m = 0; m < 2; m++)
      {
         std::unique_ptr<Mesh> sm(Box3(m == 0 ? 2 : 4));
         ParMesh pm(MPI_COMM_WORLD, *sm);
         sm.reset();
         ev[m] = RunWave(pm, p, Wv, true, Wv);
         ee[m] = RunWave(pm, p, We, false, We);
         if (m == 1) { gap = WaveGap(pm, p, We, Wv); }
      }
      if (g_rank == 0)
      {
         std::cout << std::setprecision(5) << "  " << (pw ? "P" : "S") << " wave: phase speed "
                   << c_ph << " m/s (relaxed " << c_R << ", ratio " << c_ph / c_R
                   << "), amplitude decay " << decay * 1e3 << " /km (2 Q c / w = "
                   << 1.0 / decay / 1e3 << " km)\n"
                   << "    visco run vs visco wave: " << ev[0] << " (2^3 hexes), " << ev[1]
                   << " (4^3)\n    elastic run vs elastic wave: " << ee[0] << ", " << ee[1]
                   << "\n    gap between the two exact waves after one period: " << gap << "\n";
      }
      TEST_ASSERT(ev[1] < 0.1 * gap, "T-V2 " << (pw ? "P" : "S")
                  << " error to the damped wave < 0.1 x its gap to the elastic wave");
      TEST_ASSERT(ev[1] < 3.0 * ee[1] + 1e-12, "T-V2 " << (pw ? "P" : "S")
                  << " error within 3 x the elastic discretization error");
      TEST_ASSERT(ev[1] < 0.5 * ev[0], "T-V2 " << (pw ? "P" : "S") << " error falls with h");
   }
}

// ------------------------------------------------------------------ T-V3

/// 1/2 int sum_l xi_l : C_l : xi_l (exact integral of the nodal interpolant)
/// and int D = int sum_l w_l xi_l : C_l : xi_l, all ranks.
static void MemoryEnergy(const RheologyCoupler &rc, ParMesh &pm, const GMBFit &f,
                         real_t &Em, real_t &Dm)
{
   const RheologyField &F = rc.Field();
   const FiniteElement &fe = F.NodeElement();
   MassIntegrator mi;
   DenseMatrix Me;
   const int nr = F.NodesPerElement();
   Em = Dm = 0.0;
   for (int e = 0; e < F.NumElements(); e++)
   {
      ElementTransformation *T = pm.GetElementTransformation(e);
      mi.AssembleElementMatrix(fe, *T, Me);
      for (int l = 0; l < 3; l++)
      {
         real_t s = 0.0;
         for (int i = 0; i < nr; i++)
         {
            real_t ci[6];
            sym6::ApplyIsotropic(f.lam_l[l], f.mu_l[l], F.State(e, i) + 6 * l, ci);
            for (int j = 0; j < nr; j++) { s += Me(i, j) * sym6::DoubleDot(ci, F.State(e, j) + 6 * l); }
         }
         Em += 0.5 * s;
         Dm += f.omega[l] * s;
      }
   }
   Em = RSum(Em);
   Dm = RSum(Dm);
}

struct EnergyRun
{
   std::vector<real_t> E, D;
   real_t dt = 0.0;
};

static EnergyRun RunEnergy(int p, bool visco, int nsteps)
{
   EnergyRun R;
   std::unique_ptr<Mesh> sm(Box3(2));
   ParMesh pm(MPI_COMM_WORLD, *sm);
   sm.reset();
   VecFn zero = [](const Vector &, real_t, Vector &u) { u.SetSize(3); u = 0.0; };
   auto op = MakeOp(pm, p, zero);
   std::unique_ptr<RheologyCoupler> rc;
   if (visco) { rc = MakeCoupler(*op, kQs, kQp); op->SetRheology(rc.get()); }
   auto &fes = op->Elasticity().GetFESpace();
   ParGridFunction gu(&fes);
   VectorFunctionCoefficient cu(3, [](const Vector &x, Vector &u) {
      const real_t c = 0.5 * kBox, s2 = std::pow(0.15 * kBox, 2);
      const real_t r2 = std::pow(x(0) - c, 2) + std::pow(x(1) - 0.45 * kBox, 2) + std::pow(x(2) - 0.55 * kBox, 2);
      const real_t g = 1.0e-2 * std::exp(-r2 / s2);
      u.SetSize(3);
      u(0) = g; u(1) = -0.5 * g; u(2) = 0.3 * g;
   });
   gu.ProjectCoefficient(cu);
   Vector U(gu), V(U.Size());
   V = 0.0;
   if (rc) { rc->Reset(U); }
   // a fixed dt for both runs (the visco one's), so the two are comparable
   std::unique_ptr<RheologyCoupler> probe = MakeCoupler(*op, kQs, kQp);
   const real_t dt = StepWith(*op, pm, p, probe.get());
   probe.reset();
   R.dt = dt;
   op->SetTime(0.0);
   op->SetTimeStep(dt);
   ElasticityDomainOperator<ParMesh> &E = op->Elasticity();
   ParBilinearForm mform(&fes);
   ConstantCoefficient rho(TPV102Params::rho);
   mform.AddDomainIntegrator(new VectorMassIntegrator(rho));
   mform.Assemble();
   mform.Finalize();
   const SparseMatrix &M = mform.SpMat();
   const GMBFit f = TheFit();
   Vector Ku(U.Size()), Kv(U.Size()), Mv(U.Size());
   auto energy = [&](real_t &En, real_t &Dn) {
      E.ApplyStiffness(U, Ku);
      E.ApplyStiffness(V, Kv);
      M.Mult(V, Mv);
      op->Mass().Apply(Ku, Kv);                      // Kv := M^-1 K u
      real_t e = 0.5 * (V * Mv) + 0.5 * (U * Ku) - 0.125 * dt * dt * (Ku * Kv);
      e = RSum(e);
      real_t Em = 0.0, Dm = 0.0;
      if (rc) { MemoryEnergy(*rc, pm, f, Em, Dm); }
      En = e + Em;
      Dn = Dm;
   };
   real_t t = 0.0;
   CentralDifferenceSolver ode;
   ode.Init(*op);
   real_t En, Dn;
   energy(En, Dn);
   R.E.push_back(En);
   R.D.push_back(Dn);
   for (int s = 0; s < nsteps; s++)
   {
      real_t h = dt;
      ode.Step(U, V, t, h);
      energy(En, Dn);
      R.E.push_back(En);
      R.D.push_back(Dn);
   }
   return R;
}

static void TestTV3()
{
   if (g_rank == 0) { std::cout << "\n=== T-V3 energy (2.34) never rises ===\n"; }
   const int p = 2, nsteps = 400;
   EnergyRun el = RunEnergy(p, false, nsteps);
   EnergyRun ve = RunEnergy(p, true, nsteps);
   real_t el_dev = 0.0, ve_rise = 0.0, intD = 0.0;
   for (int s = 1; s <= nsteps; s++)
   {
      el_dev = std::max(el_dev, std::abs(el.E[s] - el.E[0]));
      ve_rise = std::max(ve_rise, ve.E[s] - ve.E[s - 1]);
      intD += 0.5 * ve.dt * (ve.D[s] + ve.D[s - 1]);
   }
   const real_t drop = ve.E[0] - ve.E[nsteps];
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6) << "  p = " << p << ", " << nsteps << " steps of "
                << ve.dt << " s: elastic |E - E0| / E0 max " << el_dev / el.E[0]
                << "; visco: largest one-step rise / E0 " << ve_rise / ve.E[0]
                << ", drop / E0 " << drop / ve.E[0] << ", int D dt / E0 " << intD / ve.E[0]
                << ", drop / int D " << drop / intD << "\n";
   }
   TEST_ASSERT(el_dev < 1e-10 * el.E[0], "T-V3 the elastic run keeps E_h (1e-10)");
   TEST_ASSERT(ve_rise < 1e-10 * ve.E[0], "T-V3 E_h with attenuation never rises (1e-10 E_0)");
   TEST_ASSERT(drop > 1e-3 * ve.E[0], "T-V3 E_h with attenuation falls");
   TEST_ASSERT(std::abs(drop / intD - 1.0) < 0.05, "T-V3 the drop is int D dt to 5 %");
}

// ------------------------------------------------------------------ T-V4

static void TestTV4()
{
   if (g_rank == 0) { std::cout << "\n=== T-V4 rest round trip ===\n"; }
   for (int p = 1; p <= 3; p++)
   {
      std::unique_ptr<Mesh> sm(Box3(2));
      ParMesh pm(MPI_COMM_WORLD, *sm);
      sm.reset();
      VecFn zero = [](const Vector &, real_t, Vector &u) { u.SetSize(3); u = 0.0; };
      auto op = MakeOp(pm, p, zero);
      auto rc = MakeCoupler(*op, kQs, kQp);
      auto &fes = op->Elasticity().GetFESpace();
      ParGridFunction gu(&fes);
      VectorFunctionCoefficient cu(3, [](const Vector &x, Vector &u) {
         u.SetSize(3);
         u(0) = 1e-4 * std::sin(x(1) / 3e3) + 2e-5 * x(2) / 1e3;
         u(1) = -3e-5 * std::cos(x(0) / 4e3);
         u(2) = 1e-5 * x(0) * x(1) / 1e7;
      });
      gu.ProjectCoefficient(cu);
      Vector U(gu);
      rc->Reset(U);
      rc->Advance(U, 0.01, 0.01);
      rc->Advance(U, 0.02, 0.01);
      const RheologyField &F = rc->Field();
      int nz = 0;
      for (int e = 0; e < F.NumElements(); e++)
         for (int i = 0; i < F.NodesPerElement(); i++)
            for (int c = 0; c < 6; c++) { nz += F.Dsig(e, i)[c] != 0.0; }
      nz = ISum(nz);
      Vector Fp(U.Size());
      Fp = 0.0;
      rc->AddLoad(Fp);
      const real_t fmax = RMax(Fp.Normlinf());
      if (g_rank == 0)
      {
         std::cout << "  p = " << p << ": nonzero dsig components " << nz << ", max |P dsig| "
                   << fmax << "\n";
      }
      TEST_ASSERT(nz == 0 && fmax == 0.0, "T-V4 p = " << p << " Delta sigma = 0 and P Delta sigma = 0 exactly");
   }
}

// ------------------------------------------------------------------ T-V5

static bool BitEqual(const Vector &a, const Vector &b)
{
   if (a.Size() != b.Size()) { return false; }
   return std::memcmp(a.GetData(), b.GetData(), sizeof(real_t) * a.Size()) == 0;
}

static void TestTV5()
{
   if (g_rank == 0) { std::cout << "\n=== T-V5 two dynamic phases in a row (the reset) ===\n"; }
   const int p = 2, n1 = 60, n2 = 60;
   std::unique_ptr<Mesh> sm(Box3(2));
   ParMesh pm(MPI_COMM_WORLD, *sm);
   sm.reset();
   VecFn zero = [](const Vector &, real_t, Vector &u) { u.SetSize(3); u = 0.0; };
   VectorFunctionCoefficient cu(3, [](const Vector &x, Vector &u) {
      const real_t r2 = std::pow(x(0) - 5e3, 2) + std::pow(x(1) - 6e3, 2) + std::pow(x(2) - 7e3, 2);
      u.SetSize(3);
      u(0) = 1e-2 * std::exp(-r2 / 4e6); u(1) = 0.0; u(2) = -5e-3 * std::exp(-r2 / 4e6);
   });
   auto run_phase = [&](NitscheWaveOperator &op, Vector &U, Vector &V, int n, real_t dt) {
      op.SetTime(0.0);
      op.SetTimeStep(dt);
      real_t t = 0.0;
      CentralDifferenceSolver ode;
      ode.Init(op);
      for (int s = 0; s < n; s++) { real_t h = dt; ode.Step(U, V, t, h); }
   };
   // a: one coupler through both phases, Reset between them
   auto opa = MakeOp(pm, p, zero);
   auto rca = MakeCoupler(*opa, kQs, kQp);
   opa->SetRheology(rca.get());
   auto &fes = opa->Elasticity().GetFESpace();
   ParGridFunction gu(&fes);
   gu.ProjectCoefficient(cu);
   Vector Ua(gu), Va(Ua.Size());
   Va = 0.0;
   const real_t dt = StepWith(*opa, pm, p, rca.get());
   rca->Reset(Ua);
   run_phase(*opa, Ua, Va, n1, dt);
   const Vector U1(Ua), V1(Va);
   rca->Reset(Ua);
   run_phase(*opa, Ua, Va, n2, dt);
   // b: phase 2 alone with a new operator and coupler from (U1, V1)
   auto opb = MakeOp(pm, p, zero);
   auto rcb = MakeCoupler(*opb, kQs, kQp);
   opb->SetRheology(rcb.get());
   Vector Ub(U1), Vb(V1);
   rcb->Reset(Ub);
   run_phase(*opb, Ub, Vb, n2, dt);
   int same = (BitEqual(Ua, Ub) && BitEqual(Va, Vb)) ? 0 : 1;
   const RheologyField &Fa = rca->Field(), &Fb = rcb->Field();
   for (int e = 0; e < Fa.NumElements() && same == 0; e++)
   {
      for (int i = 0; i < Fa.NodesPerElement(); i++)
      {
         if (std::memcmp(Fa.State(e, i), Fb.State(e, i), sizeof(real_t) * Fa.StateVarsPerNode()) != 0
             || std::memcmp(Fa.Dsig(e, i), Fb.Dsig(e, i), sizeof(real_t) * 6) != 0) { same = 1; }
      }
   }
   same = ISum(same);
   const real_t vmax = RMax(Va.Normlinf());
   if (g_rank == 0) { std::cout << "  p = " << p << ", " << n1 << " + " << n2 << " steps, max |v| " << vmax << "\n"; }
   TEST_ASSERT(vmax > 0.0 && same == 0, "T-V5 phase 2 after Reset = a fresh phase 2, bit for bit (U, Udot, memory, dsig)");
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   if (g_rank == 0)
   {
      std::cout << "=== Viscoelastic attenuation in the Nitsche operator (plan rev6, V2) on "
                << Mpi::WorldSize() << " rank(s) ===\n";
   }
   std::string only = (argc > 2 && std::string(argv[1]) == "--only") ? argv[2] : "";
   auto run = [&](const char *n) { return only.empty() || only.find(n) != std::string::npos; };
   if (run("TV4")) { TestTV4(); }
   if (run("TV5")) { TestTV5(); }
   if (run("TV3")) { TestTV3(); }
   if (run("TV2")) { TestTV2(); }
   if (g_rank == 0)
   {
      std::cout << "\n" << num_passed << " passed, " << num_failed << " failed, " << num_tests
                << " total\n";
   }
   return num_failed == 0 ? 0 : 1;
}
