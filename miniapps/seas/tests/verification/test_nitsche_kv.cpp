// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_nitsche_kv.cpp — Kelvin-Voigt damping (NitscheWaveOperator::
// SetKelvinVoigt, Day et al. 2005 eq. A8), fault-free boxes, damping in every
// element (width < 0).  beta_e = eta_bar (1/3) (h_e/p) / c_s,e (KVScale::Node,
// h_e = mean edge length or 6V/A) or eta_bar dt (KVScale::Dt).
//
//   (a) eta_bar = 0: bit-identical to the undamped CentralDifferenceSolver,
//       p = 1, 2, 3 (200 steps, random U, V).
//   (b) K_eta symmetric positive semidefinite (random vectors), K_eta r = 0
//       for the six rigid motions (traction-free box), and the energy of the
//       random-field run: P_kv = -W.K_eta W <= 0 and the staggered energy
//       E^{n+1/2} = 1/2 W.M W + 1/2 U^{n+1}.K U^n (W = (U^{n+1}-U^n)/dt, the
//       half-step velocity the damping uses) with its exact balance
//       E^{n+1/2} - E^{n-1/2} = -dt/2 (W^{n+1/2} + W^{n-1/2}).K_eta W^{n-1/2}.
//   (c) largest stable C in dt = C h_min/(p^2 c_p) by bisection (blow-up over
//       1000 steps), damped / undamped, for eta_bar = 0.05, 0.1, 0.2, both
//       h_e definitions, p = 1, 2, 3; boxes "A" (3x3x2 hexes of 20x20x15 km),
//       "B" (3x3x3 hexes of 1 km) and "S" (box A with one tet cap of height
//       0.1 x the half cube: a sliver-like element).
//   (d) mode test (setup of the 10-14 velocity-jump test): 8^3 cubes of 1 km
//       split into tets, free boundaries, u0 = 0, initial x velocity
//       (i) node sawtooth, (ii) element flip, (iii) smooth resolved mode
//       sin(2 pi z / lambda), lambda = 5 h/p; E(t)/E(0) at 0.5, 1, 2 h/c_s,
//       damped over undamped, effective damping ratio, p = 1, 2.
//   (e) plane S wave (120 km), L2 error and order on 2x2x1 / 4x4x2 / 8x8x4
//       boxes, eta_bar = 0.1, p = 1, 2 (Dirichlet velocity data passed to the
//       damping term).
//   (f) cost per step damped vs undamped (box 4x4x2, p = 1, 2).
//
// Usage: mpirun -np N ./seas_test_nitsche_kv [--gate abcdef] [--p 1,2,3]
//        [--mesh ABS] [--etas 0.05,0.1,0.2] [--quick]
// Scope (Chunhui Zhao, 2026-10-15).

#include "mfem.hpp"
#include "../../domain/nitsche_wave_operator.hpp"
#include "../../solver/nitsche_cfl.hpp"
#include "../../config/tpv102_params.hpp"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; if (g_rank == 0) std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; if (g_rank == 0) std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)

static int num_tests = 0, num_passed = 0, num_failed = 0;
static int g_rank = 0;
static const real_t kL = 30.0e3;
static const real_t kCeig[4] = {0.0, 0.113, 0.212, 0.282};   // P3 measured C(p)

using KVScale = NitscheWaveOperator::KVScale;
using KVLength = NitscheWaveOperator::KVLength;
using ExactFn = std::function<void(const Vector &, real_t, Vector &)>;

struct KVSpec
{
   real_t eta = 0.0;
   KVScale scale = KVScale::Node;
   KVLength length = KVLength::Inscribed;
   std::string Name() const
   {
      if (eta == 0.0) { return "undamped"; }
      std::ostringstream o;
      o << (scale == KVScale::Dt ? "dt" : length == KVLength::MeanEdge ? "edge" : "6V/A") << " " << eta;
      return o.str();
   }
};

static real_t Dot(const Vector &a, const Vector &b, MPI_Comm comm)
{
   real_t loc = a * b, glob = 0.0;
   MPI_Allreduce(&loc, &glob, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, comm);
   return glob;
}
static real_t MaxAbs(const Vector &a, MPI_Comm comm)
{
   real_t loc = a.Size() ? a.Normlinf() : 0.0, glob = 0.0;
   MPI_Allreduce(&loc, &glob, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
   return glob;
}

// ------------------------------------------------------------------ meshes
static Mesh *Box24(int nx, int ny, int nz, real_t Lx, real_t Ly, real_t Lz, int attr)
{
   Mesh *mesh = new Mesh(Mesh::MakeCartesian3DWith24TetsPerHex(nx, ny, nz, Lx, Ly, Lz));
   for (int b = 0; b < mesh->GetNBE(); b++) { mesh->SetBdrAttribute(b, attr); }
   mesh->SetAttributes();
   return mesh;
}
/// Box A with the centre vertex of hex (1, 1, 0) moved to 0.1 x the half
/// height below its top face: the four tets on that face become flat caps.
static Mesh *SliverBox()
{
   const real_t hx = 2 * kL / 3, hz = kL / 2;
   Mesh *mesh = Box24(3, 3, 2, 2 * kL, 2 * kL, kL, 5);
   const real_t c[3] = {1.5 * hx, 1.5 * hx, 0.5 * hz};
   int best = -1;
   real_t bd = 1e300;
   for (int v = 0; v < mesh->GetNV(); v++)
   {
      const real_t *x = mesh->GetVertex(v);
      const real_t d = std::hypot(std::hypot(x[0] - c[0], x[1] - c[1]), x[2] - c[2]);
      if (d < bd) { bd = d; best = v; }
   }
   MFEM_VERIFY(bd < 1.0, "SliverBox: no vertex at the hex centre");
   mesh->GetVertex(best)[2] = hz - 0.1 * 0.5 * hz;
   return mesh;
}
static Mesh *MakeMesh(char which)
{
   if (which == 'A') { return Box24(3, 3, 2, 2 * kL, 2 * kL, kL, 5); }
   if (which == 'B') { return Box24(3, 3, 3, 3e3, 3e3, 3e3, 5); }
   return SliverBox();
}

static MaterialField TPV102Material()
{
   return MaterialField::MakeConstant(TPV102Params::lambda, TPV102Params::mu, TPV102Params::rho);
}
static const ExactFn kZero = [](const Vector &, real_t, Vector &u) { u.SetSize(3); u = 0.0; };

/// Dirichlet (attr 5, data u_D) box, or traction-free (attr 1) when free.
static std::unique_ptr<NitscheWaveOperator> MakeOp(ParMesh &pm, int p, bool free = false,
                                                   const ExactFn &uD = kZero)
{
   BoundaryConfig bc;
   bc.fault_attr = 3;   // absent
   if (free) { bc.natural_attrs = {1}; }
   else
   {
      bc.dirichlet_attrs = {5};
      bc.default_dirichlet_func = [uD](const Vector &x, real_t t, Vector &u) { uD(x, t, u); };
   }
   return std::make_unique<NitscheWaveOperator>(pm, p, TPV102Material(), bc);
}
static void SetKV(NitscheWaveOperator &op, const KVSpec &s, const BoundaryConfig *bcv = nullptr)
{
   op.SetKelvinVoigt(s.eta, s.scale, -1.0, s.length, bcv);
}

struct MassMatrix
{
   std::unique_ptr<ParBilinearForm> form;
   void Build(NitscheWaveOperator &op)
   {
      form = std::make_unique<ParBilinearForm>(&op.Elasticity().GetFESpace());
      form->AddDomainIntegrator(new VectorMassIntegrator(
         const_cast<NitscheMaterialSample &>(op.Sample()).Rho()));
      form->Assemble();
      form->Finalize();
   }
   void Mult(const Vector &x, Vector &y) const { y.SetSize(x.Size()); form->SpMat().Mult(x, y); }
};

static real_t DtUnit(NitscheWaveOperator &op, ParMesh &pm, int p)
{
   return ComputeMaxDt(op.Sample(), pm, p, 1.0);
}

// ------------------------------------------------------------------ (a)
static void GateA(int p)
{
   std::unique_ptr<Mesh> sm(MakeMesh('A'));
   ParMesh pm(MPI_COMM_WORLD, *sm);
   auto op1 = MakeOp(pm, p), op2 = MakeOp(pm, p);
   const real_t dt = 0.5 * kCeig[p] * DtUnit(*op1, pm, p);
   op1->SetTimeStep(dt);
   op2->SetTimeStep(dt);
   SetKV(*op2, KVSpec{0.0, KVScale::Node});
   const int n = op1->Height();
   Vector U1(n), V1(n);
   U1.Randomize(11 + g_rank); U1 -= 0.5;
   V1.Randomize(12 + g_rank); V1 -= 0.5;
   Vector U2(U1), V2(V1);
   real_t t1 = 0.0, t2 = 0.0;
   RunCentralDifference(*op1, U1, V1, t1, dt, 200);
   RunCentralDifference(*op2, U2, V2, t2, dt, 200);
   int diff = (std::memcmp(U1.GetData(), U2.GetData(), sizeof(real_t) * n) != 0
               || std::memcmp(V1.GetData(), V2.GetData(), sizeof(real_t) * n) != 0) ? 1 : 0;
   MPI_Allreduce(MPI_IN_PLACE, &diff, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
   TEST_ASSERT(diff == 0 && !op2->HasKelvinVoigt(),
               "(a) p=" << p << " eta_bar=0: bit-identical to the undamped CentralDifferenceSolver (200 steps)");
}

static real_t BisectC(ParMesh &pm, int p, const KVSpec &kv, real_t dt_unit, real_t lo,
                      real_t hi);

// ------------------------------------------------------------------ (b)
static void GateB(int p, bool quick)
{
   MPI_Comm comm = MPI_COMM_WORLD;
   const KVSpec kv{0.1, KVScale::Node};
   // operator checks on a traction-free 2x2x2 km box (rigid motions in the kernel)
   {
      std::unique_ptr<Mesh> sm(Box24(2, 2, 2, 2e3, 2e3, 2e3, 1));
      ParMesh pm(comm, *sm);
      auto op = MakeOp(pm, p, true);
      SetKV(*op, kv);
      const int n = op->Height();
      Vector x(n), y(n), Kx, Ky;
      real_t sym = 0.0, rq_min = 1e300, rq_max = 0.0;
      for (int k = 0; k < 20; k++)
      {
         x.Randomize(100 + 2 * k + 50 * g_rank); x -= 0.5;
         y.Randomize(101 + 2 * k + 50 * g_rank); y -= 0.5;
         op->ApplyKVStiffness(x, Kx);
         op->ApplyKVStiffness(y, Ky);
         const real_t xKx = Dot(x, Kx, comm);
         sym = std::max(sym, std::abs(Dot(y, Kx, comm) - Dot(x, Ky, comm)) / xKx);
         const real_t rq = xKx / Dot(x, x, comm);
         rq_min = std::min(rq_min, rq);
         rq_max = std::max(rq_max, rq);
      }
      // extreme eigenvalues of K_eta: Lanczos with full reorthogonalisation
      const int nl = quick ? 60 : 300;
      std::vector<Vector> Q;
      std::vector<real_t> al, be;
      Vector q(n), w, Kq;
      q.Randomize(999 + g_rank); q -= 0.5;
      q /= std::sqrt(Dot(q, q, comm));
      for (int j = 0; j < nl; j++)
      {
         Q.push_back(q);
         op->ApplyKVStiffness(q, Kq);
         w = Kq;
         al.push_back(Dot(q, w, comm));
         for (int r = 0; r < 2; r++)
            for (const Vector &qq : Q) { w.Add(-Dot(qq, w, comm), qq); }
         const real_t bj = std::sqrt(Dot(w, w, comm));
         if (bj < 1e-14 * std::abs(al[0])) { break; }
         be.push_back(bj);
         q = w;
         q /= bj;
      }
      const int m = (int) al.size();
      DenseMatrix Tm(m);
      Tm = 0.0;
      for (int j = 0; j < m; j++)
      {
         Tm(j, j) = al[j];
         if (j + 1 < m) { Tm(j, j + 1) = Tm(j + 1, j) = be[j]; }
      }
      Vector ev;
      Tm.Eigenvalues(ev);
      const real_t lam_lo = ev.Min(), lam_hi = ev.Max();
      // rigid motions: three translations, three rotations (exact at p >= 1)
      ParGridFunction g(&op->Elasticity().GetFESpace());
      real_t rig = 0.0;
      for (int m = 0; m < 6; m++)
      {
         VectorFunctionCoefficient rc(3, [m](const Vector &X, Vector &u) {
            u.SetSize(3); u = 0.0;
            const real_t xs = X(0) / 1e3, ys = X(1) / 1e3, zs = X(2) / 1e3;
            if (m < 3) { u(m) = 1.0; }
            else if (m == 3) { u(0) = -ys; u(1) = xs; }
            else if (m == 4) { u(1) = -zs; u(2) = ys; }
            else { u(2) = -xs; u(0) = zs; } });
         g.ProjectCoefficient(rc);
         Vector Kg;
         op->ApplyKVStiffness(g, Kg);
         rig = std::max(rig, MaxAbs(Kg, comm) / (lam_hi * MaxAbs(g, comm)));
      }
      if (g_rank == 0)
      {
         std::cout << std::setprecision(4) << "  (b) p = " << p << " free 2 km box, eta_bar = 0.1 node: beta "
                   << op->KelvinVoigtBetaMin() << " .. " << op->KelvinVoigtBetaMax()
                   << " s; symmetry defect " << sym << "; random x.Kx/x.x in [" << rq_min << ", "
                   << rq_max << "]; Lanczos (" << m << " steps) spectrum [" << lam_lo << ", " << lam_hi
                   << "], lambda_min/lambda_max = " << lam_lo / lam_hi << "; rigid |K r| / (lambda_max |r|) "
                   << rig << "\n";
      }
      TEST_ASSERT(sym < 1e-12, "(b) p=" << p << " K_eta symmetric");
      TEST_ASSERT(rq_min > 0.0 && lam_lo > -1e-10 * lam_hi, "(b) p=" << p << " K_eta positive semidefinite");
      TEST_ASSERT(rig < 1e-12, "(b) p=" << p << " K_eta r = 0 for the six rigid motions");
   }
   // energy on a random initial field, Dirichlet box A, production dt
   {
      std::unique_ptr<Mesh> sm(MakeMesh('A'));
      ParMesh pm(comm, *sm);
      auto op = MakeOp(pm, p);
      const real_t dt_unit = DtUnit(*op, pm, p);
      const real_t C_kv = BisectC(pm, p, kv, dt_unit, 0.02 * kCeig[p], kCeig[p]);
      const real_t dt = 0.5 * C_kv * dt_unit;
      op->SetTimeStep(dt);
      SetKV(*op, kv);
      MassMatrix M;
      M.Build(*op);
      const int n = op->Height(), nsteps = quick ? 100 : 400;
      Vector U(n), V(n), Uold, W(n), Wold(n), KU, MW, KvW, KvWo;
      U.Randomize(777 + g_rank); U -= 0.5;
      V = 0.0;
      real_t t = 0.0;
      CentralDifferenceSolver ode;
      ode.Init(*op);
      real_t Eold = 0.0, Pkv_max = -1e300, Pkv_scale = 0.0, rho_max = 0.0, inc_max = 0.0, E0 = 0.0;
      int incs = 0;
      for (int s = 0; s < nsteps; s++)
      {
         Uold = U;
         real_t h = dt;
         ode.Step(U, V, t, h);
         W = U;
         W -= Uold;
         W /= dt;
         op->Elasticity().ApplyStiffness(Uold, KU);
         M.Mult(W, MW);
         const real_t E = 0.5 * Dot(W, MW, comm) + 0.5 * Dot(U, KU, comm);
         op->ApplyKVStiffness(W, KvW);
         const real_t Pkv = -Dot(W, KvW, comm);
         Pkv_max = std::max(Pkv_max, Pkv);
         Pkv_scale = std::max(Pkv_scale, std::abs(Pkv));
         if (s == 0) { E0 = E; }
         else
         {
            // exact balance: E+ - E- = -dt/2 (W+ + W-).K_eta W-
            op->ApplyKVStiffness(Wold, KvWo);
            const real_t bal = -0.5 * dt * (Dot(W, KvWo, comm) + Dot(Wold, KvWo, comm));
            rho_max = std::max(rho_max, std::abs(E - Eold - bal) / E0);
            if (E > Eold) { incs++; inc_max = std::max(inc_max, (E - Eold) / Eold); }
         }
         Eold = E;
         Wold = W;
      }
      if (g_rank == 0)
      {
         std::cout << std::setprecision(4) << "  (b) p = " << p << " random field, box A, damped stable C = " << C_kv
                   << ", run at dt = 0.5 x that (production 0.5 C(p) = " << 0.5 * kCeig[p] << "), "
                   << nsteps << " steps: E^{n+1/2} " << E0 << " -> " << Eold << " (ratio "
                   << Eold / E0 << "); increases " << incs << " (max rel " << inc_max
                   << "); max P_kv " << Pkv_max << " (scale " << Pkv_scale
                   << "); balance residual max |rho|/E0 " << rho_max << "\n";
      }
      TEST_ASSERT(Pkv_max <= 1e-12 * Pkv_scale, "(b) p=" << p << " P_kv <= 0 at every step");
      TEST_ASSERT(incs == 0, "(b) p=" << p << " staggered energy never increases (random field)");
      TEST_ASSERT(rho_max < 1e-9, "(b) p=" << p << " energy balance with the half-step velocity holds to round-off");
   }
}

// ------------------------------------------------------------------ (c)
static bool BlowsUp(ParMesh &pm, int p, real_t dt, const KVSpec &kv, int nsteps)
{
   auto op = MakeOp(pm, p);
   op->SetTimeStep(dt);
   SetKV(*op, kv);
   const int n = op->Height();
   Vector U(n), V(n);
   U.Randomize(777 + g_rank); U -= 0.5;
   V = 0.0;
   MPI_Comm comm = pm.GetComm();
   const real_t n0 = std::sqrt(Dot(U, U, comm));
   real_t t = 0.0;
   CentralDifferenceSolver ode;
   ode.Init(*op);
   for (int s = 0; s < nsteps; ++s)
   {
      real_t h = dt;
      ode.Step(U, V, t, h);
      if (s % 50 == 49 || s == nsteps - 1)
      {
         const real_t nu = std::sqrt(Dot(U, U, comm));
         if (!std::isfinite(nu) || nu > 1e3 * n0) { return true; }
      }
   }
   return false;
}
static real_t BisectC(ParMesh &pm, int p, const KVSpec &kv, real_t dt_unit, real_t lo, real_t hi)
{
   const int nsteps = 1000;
   int guard = 0;
   while (BlowsUp(pm, p, lo * dt_unit, kv, nsteps) && guard++ < 10) { lo *= 0.7; }
   guard = 0;
   while (!BlowsUp(pm, p, hi * dt_unit, kv, nsteps) && guard++ < 10) { hi *= 1.3; }
   while ((hi - lo) / hi > 0.01)
   {
      const real_t mid = 0.5 * (lo + hi);
      if (BlowsUp(pm, p, mid * dt_unit, kv, nsteps)) { hi = mid; } else { lo = mid; }
   }
   return lo;
}
static void GateC(char mesh_id, int p, const std::vector<real_t> &etas)
{
   std::unique_ptr<Mesh> sm(MakeMesh(mesh_id));
   ParMesh pm(MPI_COMM_WORLD, *sm);
   const long long nglob = pm.GetGlobalNE();
   auto op = MakeOp(pm, p);
   const real_t dt_unit = DtUnit(*op, pm, p);
   // element lengths: min 6V/A, min mean edge, max (mean edge)/(6V/A)
   real_t hs[3] = {1e300, 1e300, 0.0};
   for (int e = 0; e < pm.GetNE(); e++)
   {
      const real_t he = op->KVElementLength(e, KVLength::MeanEdge);
      const real_t hi = op->KVElementLength(e, KVLength::Inscribed);
      hs[0] = std::min(hs[0], hi);
      hs[1] = std::min(hs[1], he);
      hs[2] = std::max(hs[2], he / hi);
   }
   MPI_Allreduce(MPI_IN_PLACE, hs, 2, MPITypeMap<real_t>::mpi_type, MPI_MIN, MPI_COMM_WORLD);
   MPI_Allreduce(MPI_IN_PLACE, &hs[2], 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, MPI_COMM_WORLD);
   const real_t C0 = BisectC(pm, p, KVSpec{}, dt_unit, 0.5 * kCeig[p], 1.2 * kCeig[p]);
   if (g_rank == 0)
   {
      std::cout << std::setprecision(5) << "  (c) mesh " << mesh_id << " (" << nglob
                << " tets; min 6V/A " << hs[0] << " m, min mean edge " << hs[1]
                << " m, max edge/(6V/A) " << hs[2] << "), p = " << p << ": undamped C0 = " << C0
                << " (C(p) = " << kCeig[p] << ")\n" << std::flush;
   }
   for (int il = 0; il < 2; il++)
   {
      for (real_t eta : etas)
      {
         const KVSpec kv{eta, KVScale::Node, il == 0 ? KVLength::MeanEdge : KVLength::Inscribed};
         const real_t C = BisectC(pm, p, kv, dt_unit, 0.1 * C0, C0);
         if (g_rank == 0)
         {
            std::cout << std::setprecision(4) << "  (c) mesh " << mesh_id << " p = " << p << " eta_bar = "
                      << eta << " h_e = " << (il == 0 ? "mean edge" : "6V/A     ") << ": C = " << C
                      << "  factor " << C / C0 << "  vs production 0.5 C(p) = " << 0.5 * kCeig[p]
                      << ": margin " << C / (0.5 * kCeig[p]) << "x\n" << std::flush;
         }
      }
   }
}

// ------------------------------------------------------------------ (d)
static real_t PlainEnergy(NitscheWaveOperator &op, const MassMatrix &M, const Vector &U,
                          const Vector &V, MPI_Comm comm)
{
   Vector KU, MV;
   op.Elasticity().ApplyStiffness(U, KU);
   M.Mult(V, MV);
   return 0.5 * Dot(V, MV, comm) + 0.5 * Dot(U, KU, comm);
}
static void GateD(int p)
{
   MPI_Comm comm = MPI_COMM_WORLD;
   const real_t h = 1000.0, Lb = 8 * h;
   std::unique_ptr<Mesh> sm(new Mesh(Mesh::MakeCartesian3D(8, 8, 8, Element::TETRAHEDRON, Lb, Lb, Lb)));
   for (int b = 0; b < sm->GetNBE(); b++) { sm->SetBdrAttribute(b, 1); }
   sm->SetAttributes();
   ParMesh pm(comm, *sm);
   const long long nglob = pm.GetGlobalNE();
   const real_t cs = std::sqrt(TPV102Params::mu / TPV102Params::rho);
   const real_t tcross = h / cs;
   const char *names[3] = {"node sawtooth", "element flip", "smooth 5 n/wl"};
   const real_t ts[3] = {0.5, 1.0, 2.0};
   const int NS = 5;
   const KVSpec specs[NS] = {KVSpec{}, KVSpec{0.1, KVScale::Node, KVLength::Inscribed},
                             KVSpec{0.15, KVScale::Node, KVLength::Inscribed},
                             KVSpec{0.04, KVScale::Node, KVLength::MeanEdge},
                             KVSpec{0.1, KVScale::Dt}};
   real_t res[3][NS][3], om[3] = {0, 0, 0}, beta[NS][2] = {};
   real_t dt = 0.0;
   for (int is = 0; is < NS; is++)
   {
      auto op = MakeOp(pm, p, true);
      MassMatrix M;
      M.Build(*op);
      dt = 0.5 * kCeig[p] * DtUnit(*op, pm, p);
      op->SetTimeStep(dt);
      SetKV(*op, specs[is]);
      if (is > 0) { beta[is][0] = op->KelvinVoigtBetaMin(); beta[is][1] = op->KelvinVoigtBetaMax(); }
      auto &fes = op->Elasticity().GetFESpace();
      const int n = op->Height(), nsc = fes.GetNDofs();
      ParGridFunction X(&fes);
      VectorFunctionCoefficient xc(3, [](const Vector &x, Vector &u) { u = x; });
      X.ProjectCoefficient(xc);
      for (int m = 0; m < 3; m++)
      {
         Vector U(n), V(n);
         U = 0.0;
         V = 0.0;
         const real_t lam = 5.0 * h / p;
         for (int i = 0; i < nsc; i++)
         {
            const real_t x = X(i), y = X(nsc + i), z = X(2 * nsc + i);
            real_t v = 0.0;
            if (m == 0)
            {
               const long ix = std::lround(x * p / h), iy = std::lround(y * p / h),
                          iz = std::lround(z * p / h);
               v = ((ix + iy + iz) % 2 == 0) ? 1.0 : -1.0;
            }
            else if (m == 2) { v = std::sin(2.0 * M_PI * z / lam); }
            V(i) = v;
         }
         if (m == 1)
         {
            Array<int> vd;
            for (int e = 0; e < pm.GetNE(); e++)
            {
               Vector c(3);
               pm.GetElementCenter(e, c);
               const long cx = (long) std::floor(c(0) / h), cy = (long) std::floor(c(1) / h),
                          cz = (long) std::floor(c(2) / h);
               const real_t sgn = ((cx + cy + cz) % 2 == 0) ? 1.0 : -1.0;
               fes.GetElementDofs(e, vd);
               for (int k = 0; k < vd.Size(); k++) { V(vd[k]) = sgn; }
            }
         }
         if (is == 0)
         {
            Vector KV, MV;
            op->Elasticity().ApplyStiffness(V, KV);
            M.Mult(V, MV);
            om[m] = std::sqrt(Dot(V, KV, comm) / Dot(V, MV, comm));
         }
         const real_t E0 = PlainEnergy(*op, M, U, V, comm);
         real_t t = 0.0;
         int next = 0;
         const int nsteps = (int) std::ceil(2.0 * tcross / dt) + 1;
         CentralDifferenceSolver ode;
         ode.Init(*op);
         for (int s = 0; s < nsteps && next < 3; s++)
         {
            real_t hh = dt;
            ode.Step(U, V, t, hh);
            while (next < 3 && t >= ts[next] * tcross - 1e-12)
            {
               res[m][is][next] = PlainEnergy(*op, M, U, V, comm) / E0;
               next++;
            }
         }
      }
   }
   if (g_rank == 0)
   {
      std::cout << "\n=== (d) mode test, p = " << p << ": " << nglob
                << " tets (8^3 cubes of 1 km, free), h/cs = " << tcross << " s, dt = " << dt
                << " s ===\n  beta range (s):";
      for (int is = 1; is < NS; is++)
      {
         std::cout << " " << specs[is].Name() << " " << beta[is][0] << ".." << beta[is][1] << ";";
      }
      std::cout << "\n  E(t)/E(0) at t = 0.5, 1, 2 h/cs\n  mode            f (Hz)";
      for (int is = 0; is < NS; is++) { std::cout << "  " << std::left << std::setw(23) << specs[is].Name() << std::right; }
      std::cout << "\n";
      for (int m = 0; m < 3; m++)
      {
         std::cout << "  " << std::left << std::setw(15) << names[m] << std::right << std::setw(7)
                   << std::setprecision(3) << std::fixed << om[m] / (2 * M_PI);
         for (int is = 0; is < NS; is++)
         {
            std::cout << "  ";
            for (int k = 0; k < 3; k++) { std::cout << std::setw(7) << res[m][is][k]; }
         }
         std::cout << std::defaultfloat << "\n";
      }
      std::cout << "  damped/undamped energy at 2 h/cs -> decay rate r = -ln(E/E_u)/t (1/s), damping ratio "
                   "zeta = r/(2 omega) [KV theory beta_mean omega/2]\n";
      for (int m = 0; m < 3; m++)
      {
         std::cout << "  " << std::left << std::setw(15) << names[m] << std::right;
         for (int is = 1; is < NS; is++)
         {
            const real_t ratio = res[m][is][2] / res[m][0][2];
            const real_t r = -std::log(ratio) / (2.0 * tcross);
            const real_t bm = 0.5 * (beta[is][0] + beta[is][1]);
            std::cout << std::setprecision(3) << "  " << specs[is].Name().c_str() << ": E/E_u " << ratio
                      << ", r " << r << ", zeta " << r / (2 * om[m]) << " [" << 0.5 * bm * om[m] << "]";
         }
         std::cout << "\n";
      }
   }
   bool ok = true;
   for (int m = 0; m < 3; m++)
      for (int is = 1; is < NS; is++)
         for (int k = 0; k < 3; k++) { if (res[m][is][k] > res[m][0][k] * (1 + 1e-6)) { ok = false; } }
   TEST_ASSERT(ok, "(d) p=" << p << " damping never leaves more energy than undamped (all modes, times)");
}

// ------------------------------------------------------------------ (e)
static void GateE(int p, bool quick)
{
   MPI_Comm comm = MPI_COMM_WORLD;
   const real_t cs = std::sqrt(TPV102Params::mu / TPV102Params::rho);
   real_t nv[3] = {1.0, 0.7, 0.4};
   const real_t nn = std::sqrt(nv[0] * nv[0] + nv[1] * nv[1] + nv[2] * nv[2]);
   for (auto &c : nv) { c /= nn; }
   const real_t kmag = 2.0 * M_PI / 120.0e3, om = cs * kmag;
   real_t d[3] = {nv[1], -nv[0], 0.0};
   const real_t dn = std::sqrt(d[0] * d[0] + d[1] * d[1]);
   for (auto &c : d) { c /= dn; }
   ExactFn ue = [=](const Vector &x, real_t t, Vector &u) {
      const real_t ph = kmag * (nv[0] * x(0) + nv[1] * x(1) + nv[2] * x(2)) - om * t;
      u.SetSize(3); for (int i = 0; i < 3; i++) { u(i) = d[i] * std::sin(ph); } };
   ExactFn ve = [=](const Vector &x, real_t t, Vector &v) {
      const real_t ph = kmag * (nv[0] * x(0) + nv[1] * x(1) + nv[2] * x(2)) - om * t;
      v.SetSize(3); for (int i = 0; i < 3; i++) { v(i) = -om * d[i] * std::cos(ph); } };
   BoundaryConfig bcv;
   bcv.fault_attr = 3;
   bcv.dirichlet_attrs = {5};
   bcv.default_dirichlet_func = ve;
   const int nref[3][3] = {{2, 2, 1}, {4, 4, 2}, {8, 8, 4}};
   const int nm = quick ? 2 : 3;
   std::unique_ptr<ParMesh> pm[3];
   for (int m = 0; m < nm; m++)
   {
      std::unique_ptr<Mesh> sm(Box24(nref[m][0], nref[m][1], nref[m][2], 2 * kL, 2 * kL, kL, 5));
      for (int v = 0; v < sm->GetNV(); v++)
      {
         real_t *c = sm->GetVertex(v);
         c[0] -= kL; c[1] -= kL; c[2] -= kL;
      }
      pm[m] = std::make_unique<ParMesh>(comm, *sm);
   }
   const real_t T = 2.0;
   real_t dt = 0.0;
   {
      auto opf = MakeOp(*pm[nm - 1], p, false, ue);
      dt = 0.25 * kCeig[p] * DtUnit(*opf, *pm[nm - 1], p);
   }
   const int nsteps = (int) std::ceil(T / dt);
   dt = T / nsteps;
   const KVSpec specs[3] = {KVSpec{}, KVSpec{0.1, KVScale::Node, KVLength::Inscribed},
                            KVSpec{0.1, KVScale::Dt}};
   real_t err[3][3], beta[3][3];
   for (int is = 0; is < 3; is++)
   {
      for (int m = 0; m < nm; m++)
      {
         auto op = MakeOp(*pm[m], p, false, ue);
         op->SetTimeStep(dt);
         SetKV(*op, specs[is], &bcv);
         beta[is][m] = is ? op->KelvinVoigtBetaMax() : 0.0;
         auto &fes = op->Elasticity().GetFESpace();
         ParGridFunction gu(&fes), gv(&fes);
         VectorFunctionCoefficient cu(3, ue), cv(3, ve);
         cu.SetTime(0.0); cv.SetTime(0.0);
         gu.ProjectCoefficient(cu); gv.ProjectCoefficient(cv);
         Vector U(gu), V(gv);
         real_t t = 0.0;
         op->SetTime(0.0);
         RunCentralDifference(*op, U, V, t, dt, nsteps);
         gu = U;
         cu.SetTime(t);
         ParGridFunction zero(&fes);
         zero = 0.0;
         err[is][m] = gu.ComputeL2Error(cu) / zero.ComputeL2Error(cu);
      }
   }
   if (g_rank == 0)
   {
      std::cout << std::setprecision(4) << "  (e) p = " << p << ", plane S wave (120 km, omega = " << om
                << " rad/s), T = " << T << " s, dt = " << dt << " s:\n";
      for (int is = 0; is < 3; is++)
      {
         std::cout << "     " << std::left << std::setw(10) << specs[is].Name().c_str() << std::right << ": rel L2 error";
         for (int m = 0; m < nm; m++) { std::cout << " " << err[is][m]; }
         if (nm == 3)
         {
            std::cout << "; order " << std::log2(err[is][0] / err[is][1]) << ", "
                      << std::log2(err[is][1] / err[is][2]);
         }
         if (is)
         {
            std::cout << "; beta_max";
            for (int m = 0; m < nm; m++) { std::cout << " " << beta[is][m]; }
            std::cout << " s; KV amplitude loss beta omega^2 T/2";
            for (int m = 0; m < nm; m++) { std::cout << " " << 0.5 * beta[is][m] * om * om * T; }
         }
         std::cout << "\n";
      }
   }
   if (nm == 3)
   {
      const real_t r2 = std::log2(err[1][1] / err[1][2]);
      TEST_ASSERT(r2 > p + 0.5, "(e) p=" << p << " eta_bar=0.1 node/6V/A: L2 order on the finest pair > p + 0.5");
   }
}

// ------------------------------------------------------------------ (f)
static void GateF(int p)
{
   std::unique_ptr<Mesh> sm(Box24(4, 4, 2, 2 * kL, 2 * kL, kL, 5));
   ParMesh pm(MPI_COMM_WORLD, *sm);
   const long long nglob = pm.GetGlobalNE();
   real_t tps[2];
   const int nsteps = 300;
   for (int is = 0; is < 2; is++)
   {
      auto op = MakeOp(pm, p);
      const real_t dt = 0.5 * kCeig[p] * DtUnit(*op, pm, p);
      op->SetTimeStep(dt);
      if (is) { SetKV(*op, KVSpec{0.1, KVScale::Node}); }
      const int n = op->Height();
      Vector U(n), V(n);
      U.Randomize(5 + g_rank); U -= 0.5;
      V = 0.0;
      real_t t = 0.0;
      RunCentralDifference(*op, U, V, t, dt, 10);   // warm-up
      MPI_Barrier(MPI_COMM_WORLD);
      const double t0 = MPI_Wtime();
      RunCentralDifference(*op, U, V, t, dt, nsteps);
      MPI_Barrier(MPI_COMM_WORLD);
      tps[is] = (MPI_Wtime() - t0) / nsteps;
   }
   if (g_rank == 0)
   {
      std::cout << std::setprecision(4) << "  (f) p = " << p << ", box 4x4x2 (" << nglob
                << " tets), np " << Mpi::WorldSize() << ": " << 1e3 * tps[0] << " ms/step undamped, "
                << 1e3 * tps[1] << " ms/step damped (x" << tps[1] / tps[0] << ")\n";
   }
}


// ------------------------------------------------------------------ (g)
// Recheck of (d): spectral content.  Lanczos in the M inner product on
// M^{-1} K (self-adjoint there) from a start vector v gives Ritz values
// theta_j = omega_j^2 and weights w_j (squared first components) = the
// fraction of the M-norm of v at omega_j.  For Kelvin-Voigt a mode's energy
// decays as exp(-beta omega^2 t), so the predicted energy left of v is
// sum_j w_j exp(-beta theta_j t) (velocity data, u0 = 0, light damping).
struct Spectrum { std::vector<real_t> th, w; };
static Spectrum LanczosM(NitscheWaveOperator &op, const MassMatrix &M, const Vector &v0, int nl,
                         MPI_Comm comm)
{
   const int n = v0.Size();
   std::vector<Vector> Q;
   std::vector<real_t> al, be;
   Vector q(v0), Mq, Kq, w(n), Mw;
   M.Mult(q, Mq);
   q /= std::sqrt(Dot(q, Mq, comm));
   for (int j = 0; j < nl; j++)
   {
      Q.push_back(q);
      op.Elasticity().ApplyStiffness(q, Kq);
      op.Mass().Apply(Kq, w);
      al.push_back(Dot(q, Kq, comm));
      for (int r = 0; r < 2; r++)
         for (const Vector &qq : Q)
         {
            M.Mult(qq, Mq);
            w.Add(-Dot(Mq, w, comm), qq);
         }
      M.Mult(w, Mw);
      const real_t bj = std::sqrt(std::max(Dot(w, Mw, comm), 0.0));
      if (bj < 1e-12 * std::abs(al[0])) { break; }
      be.push_back(bj);
      q = w;
      q /= bj;
   }
   const int m = (int) al.size();
   DenseMatrix T(m), V(m);
   T = 0.0;
   for (int j = 0; j < m; j++)
   {
      T(j, j) = al[j];
      if (j + 1 < m) { T(j, j + 1) = T(j + 1, j) = be[j]; }
   }
   Vector ev;
   T.Eigensystem(ev, V);
   Spectrum S;
   for (int k = 0; k < m; k++) { S.th.push_back(std::max(ev(k), 0.0)); S.w.push_back(V(0, k) * V(0, k)); }
   return S;
}
static void GateG(int p)
{
   MPI_Comm comm = MPI_COMM_WORLD;
   const real_t h = 1000.0, Lb = 8 * h;
   std::unique_ptr<Mesh> sm(new Mesh(Mesh::MakeCartesian3D(8, 8, 8, Element::TETRAHEDRON, Lb, Lb, Lb)));
   for (int b = 0; b < sm->GetNBE(); b++) { sm->SetBdrAttribute(b, 1); }
   sm->SetAttributes();
   ParMesh pm(comm, *sm);
   const real_t cs = std::sqrt(TPV102Params::mu / TPV102Params::rho), tcross = h / cs, T2 = 2 * tcross;
   auto op0 = MakeOp(pm, p, true);
   MassMatrix M;
   M.Build(*op0);
   const real_t dt = 0.5 * kCeig[p] * DtUnit(*op0, pm, p);
   const int n = op0->Height();
   auto &fes = op0->Elasticity().GetFESpace();
   const int nsc = fes.GetNDofs();
   ParGridFunction X(&fes);
   VectorFunctionCoefficient xc(3, [](const Vector &x, Vector &u) { u = x; });
   X.ProjectCoefficient(xc);
   // start vectors (x velocity): top eigenmode, seeded sawtooth, cos(n pi z/L) at the noise frequencies
   std::vector<std::string> names;
   std::vector<Vector> vs;
   {
      Vector v(n), Kv, w(n);
      v.Randomize(4242 + g_rank); v -= 0.5;
      for (int it = 0; it < 3000; it++)
      {
         op0->Elasticity().ApplyStiffness(v, Kv);
         op0->Mass().Apply(Kv, w);
         v = w;
         v /= std::sqrt(Dot(v, v, comm));
      }
      names.push_back("top eigenmode (power it.)"); vs.push_back(v);
   }
   {
      Vector v(n);
      v = 0.0;
      for (int i = 0; i < nsc; i++)
      {
         const long ix = std::lround(X(i) * p / h), iy = std::lround(X(nsc + i) * p / h),
                    iz = std::lround(X(2 * nsc + i) * p / h);
         v(i) = ((ix + iy + iz) % 2 == 0) ? 1.0 : -1.0;
      }
      names.push_back("seeded node sawtooth"); vs.push_back(v);
   }
   // TPV102 200 m noise: 12 Hz (p=1), 23 Hz (p=2) -> omega h/cs at h = 200 m
   const real_t wn[2] = {2 * M_PI * 12.0 * 200.0 / cs, 2 * M_PI * 23.0 * 200.0 / cs};
   for (int k = 0; k < 2; k++)
   {
      const int nz = std::max(1, (int) std::lround(wn[k] * Lb / (h * M_PI)));   // k_z = n pi / L
      Vector v(n);
      v = 0.0;
      for (int i = 0; i < nsc; i++) { v(i) = std::cos(nz * M_PI * X(2 * nsc + i) / Lb); }
      std::ostringstream o;
      o << "shear cos(" << nz << " pi z/L) [TPV " << (k == 0 ? "12 Hz p1" : "23 Hz p2") << ", omega h/cs "
        << std::setprecision(3) << nz * M_PI * h / Lb << " vs " << wn[k] << "]";
      names.push_back(o.str()); vs.push_back(v);
   }
   const KVSpec kv{0.1, KVScale::Node, KVLength::Inscribed};
   auto opd = MakeOp(pm, p, true);
   opd->SetTimeStep(dt);
   SetKV(*opd, kv);
   const real_t beta = 0.5 * (opd->KelvinVoigtBetaMin() + opd->KelvinVoigtBetaMax());
   op0->SetTimeStep(dt);
   auto run = [&](NitscheWaveOperator &op, const Vector &v0) {
      Vector U(n), V(v0);
      U = 0.0;
      const real_t E0 = PlainEnergy(op, M, U, V, comm);
      real_t t = 0.0;
      CentralDifferenceSolver ode;
      ode.Init(op);
      const int ns = (int) std::ceil(T2 / dt);
      for (int s = 0; s < ns; s++) { real_t hh = dt; ode.Step(U, V, t, hh); }
      return PlainEnergy(op, M, U, V, comm) / E0;
   };
   if (g_rank == 0)
   {
      std::cout << "\n=== (g) spectral recheck, p = " << p << ", 8^3 1 km box (free), dt = " << dt
                << " s, eta_bar 0.1 6V/A beta = " << beta << " s, t = 2 h/cs = " << T2 << " s ===\n";
   }
   const real_t bins[] = {0, 2, 4, 6, 8, 10, 14, 20, 30, 50, 1e9};
   const int nb = 10;
   for (size_t m = 0; m < vs.size(); m++)
   {
      Spectrum S = LanczosM(*op0, M, vs[m], 200, comm);
      real_t wsum = 0, wmean = 0, pred = 0, wtop = 0;
      real_t hist[10] = {0};
      for (size_t j = 0; j < S.th.size(); j++)
      {
         const real_t om = std::sqrt(S.th[j]), x = om * h / cs;
         wsum += S.w[j];
         wmean += S.w[j] * S.th[j];
         pred += S.w[j] * std::exp(-beta * S.th[j] * T2);
         wtop = std::max(wtop, x);
         for (int b = 0; b < nb; b++) { if (x >= bins[b] && x < bins[b + 1]) { hist[b] += S.w[j]; } }
      }
      const real_t Eu = run(*op0, vs[m]), Ed = run(*opd, vs[m]);
      if (g_rank == 0)
      {
         std::cout << "  " << names[m] << ":\n     Rayleigh omega h/cs = " << std::setprecision(4)
                   << std::sqrt(wmean / wsum) * h / cs << " (" << std::sqrt(wmean / wsum) / (2 * M_PI)
                   << " Hz); top Ritz omega h/cs " << wtop << "\n     spectral weight by omega h/cs bin:";
         for (int b = 0; b < nb; b++)
         {
            if (hist[b] / wsum > 0.005) { std::cout << " [" << bins[b] << "," << (b + 1 < nb ? bins[b + 1] : 999) << "): " << std::setprecision(3) << hist[b] / wsum; }
         }
         std::cout << "\n     energy left at 2 h/cs: undamped " << std::setprecision(4) << Eu << ", damped " << Ed
                   << ", ratio " << Ed / Eu << "; predicted from the spectrum sum w exp(-beta omega^2 t) "
                   << pred / wsum << "\n";
         if (m >= 2)
         {
            const real_t om = std::sqrt(wmean / wsum);
            std::cout << "     damping ratio zeta = -ln(ratio)/(2 omega t) " << -std::log(Ed / Eu) / (2 * om * T2)
                      << " [beta omega/2 = " << 0.5 * beta * om << "]\n";
         }
      }
   }
}

static std::vector<real_t> ParseList(const std::string &a)
{
   std::vector<real_t> r;
   size_t q = 0;
   while (q <= a.size())
   {
      size_t e = a.find(',', q);
      if (e == std::string::npos) { e = a.size(); }
      r.push_back(std::atof(a.substr(q, e - q).c_str()));
      q = e + 1;
   }
   return r;
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   bool quick = false;
   std::string gates = "abcdef", meshes = "ABS";
   std::vector<real_t> ps = {1, 2, 3}, etas = {0.05, 0.1, 0.2};
   for (int i = 1; i < argc; i++)
   {
      const std::string a = argv[i];
      if (a == "--quick") { quick = true; }
      else if (a == "--gate" && i + 1 < argc) { gates = argv[++i]; }
      else if (a == "--p" && i + 1 < argc) { ps = ParseList(argv[++i]); }
      else if (a == "--mesh" && i + 1 < argc) { meshes = argv[++i]; }
      else if (a == "--etas" && i + 1 < argc) { etas = ParseList(argv[++i]); }
   }
   auto has = [&](char g) { return gates.find(g) != std::string::npos; };
   for (real_t pr : ps)
   {
      const int p = (int) pr;
      if (g_rank == 0) { std::cout << "\n=== p = " << p << " ===\n" << std::flush; }
      if (has('a')) { GateA(p); }
      if (has('b')) { GateB(p, quick); }
      if (has('c')) { for (char m : meshes) { GateC(m, p, etas); } }
      if (has('d') && p <= 2) { GateD(p); }
      if (has('e') && p <= 2) { GateE(p, quick); }
      if (has('f') && p <= 2) { GateF(p); }
      if (has('g') && p <= 2) { GateG(p); }
   }
   if (g_rank == 0)
   {
      std::cout << "\n=== Summary: " << num_passed << "/" << num_tests << " passed, "
                << num_failed << " failed ===\n";
   }
   return num_failed == 0 ? 0 : 1;
}
