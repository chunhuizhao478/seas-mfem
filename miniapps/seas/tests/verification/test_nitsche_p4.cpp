// Nitsche dynamic rupture, phase P4 (absorbing boundary: first-order
// Lysmer-Kuhlemeyer dashpot).  Plan items C-03, N-06, C-02 (damping block),
// tests T-03 and T-12(a)-(d).
//
//   T-03:   an absorbing tag lands in neither the Dirichlet nor the fault
//           set; absorbing faces are collected; the sets are disjoint.
//   C-02:   (M + dt/2 C) * Apply(r) = r, M and C formed independently.
//
// Scope (Chunhui Zhao, 2026-10-14): tetrahedral meshes only; PML is future
// work.

#include "mfem.hpp"
#include "../../domain/nitsche_wave_operator.hpp"
#include "../../solver/nitsche_cfl.hpp"
#include "../../config/tpv102_params.hpp"

#include <array>
#include <functional>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>
#include <map>
#include <set>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;
static int g_rank = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; if (g_rank == 0) std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; if (g_rank == 0) std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static const real_t kL = 30.0e3;   // box [-L, L] x [-L, L] x [-L, 0]

static int GlobalSum(int v, MPI_Comm comm)
{
   int g = 0;
   MPI_Allreduce(&v, &g, 1, MPI_INT, MPI_SUM, comm);
   return g;
}

// 24-tet TPV102 box with the fault on Y=0 (attr 3), free surface Z=0
// (attr 1), absorbing bottom Z=-L (attr 6), Dirichlet sides (attr 5).
static Mesh *FaultBoxAbsorbingBottom(int nx, int ny, int nz)
{
   Mesh *mesh = new Mesh(Mesh::MakeCartesian3DWith24TetsPerHex(
      nx, ny, nz, 2*kL, 2*kL, kL));
   for (int v = 0; v < mesh->GetNV(); v++)
   {
      real_t *c = mesh->GetVertex(v);
      c[0] -= kL; c[1] -= kL; c[2] -= kL;
   }
   const real_t tol = kL / nz * 0.1;
   for (int b = 0; b < mesh->GetNBE(); b++)
   {
      Array<int> verts;
      mesh->GetBdrElementVertices(b, verts);
      real_t zc = 0.0;
      for (int i = 0; i < verts.Size(); i++) { zc += mesh->GetVertex(verts[i])[2]; }
      zc /= verts.Size();
      int attr = 5;
      if (std::abs(zc) < tol) { attr = 1; }
      else if (std::abs(zc + kL) < tol) { attr = 6; }
      mesh->SetBdrAttribute(b, attr);
   }
   const real_t ytol = 1e-6 * kL;
   std::vector<std::array<int, 3>> tris;
   for (int f = 0; f < mesh->GetNumFaces(); f++)
   {
      if (mesh->GetFaceInformation(f).IsBoundary()) { continue; }
      Array<int> fv;
      mesh->GetFaceVertices(f, fv);
      bool on_plane = true;
      for (int i = 0; i < fv.Size(); i++)
      {
         if (std::abs(mesh->GetVertex(fv[i])[1]) > ytol) { on_plane = false; }
      }
      if (on_plane) { tris.push_back({fv[0], fv[1], fv[2]}); }
   }
   for (const auto &t : tris) { mesh->AddBdrTriangle(t[0], t[1], t[2], 3); }
   mesh->FinalizeTopology();
   mesh->Finalize();
   mesh->SetAttributes();
   return mesh;
}


// ---------------------------------------------------------------------------
// Fault-free box [-Lx/2, Lx/2] x [-Ly/2, Ly/2] x [-Lz, 0], 24 tets per hex;
// boundary attrs from the face centroid via attr_of.
// ---------------------------------------------------------------------------
static Mesh *BoxMeshAttr(int nx, int ny, int nz, real_t Lx, real_t Ly,
                         real_t Lz, const std::function<int(const Vector &)> &attr_of)
{
   Mesh *mesh = new Mesh(Mesh::MakeCartesian3DWith24TetsPerHex(
      nx, ny, nz, Lx, Ly, Lz));
   for (int v = 0; v < mesh->GetNV(); v++)
   {
      real_t *c = mesh->GetVertex(v);
      c[0] -= 0.5 * Lx; c[1] -= 0.5 * Ly; c[2] -= Lz;
   }
   for (int b = 0; b < mesh->GetNBE(); b++)
   {
      Array<int> verts;
      mesh->GetBdrElementVertices(b, verts);
      Vector c(3);
      c = 0.0;
      for (int i = 0; i < verts.Size(); i++)
         for (int d = 0; d < 3; d++) { c(d) += mesh->GetVertex(verts[i])[d]; }
      c /= verts.Size();
      mesh->SetBdrAttribute(b, attr_of(c));
   }
   mesh->SetAttributes();
   return mesh;
}

// Pseudo-random density in [1000, 10000] (contrast up to 10x).
static real_t RhoField(const Vector &x)
{
   const real_t s = std::sin(1.3e-3 * x(0) + 2.1e-3 * x(1) + 0.7e-3 * x(2)
                             + 1.7e-7 * x(0) * x(2));
   return 1000.0 + 9000.0 * s * s;
}

static real_t Dot(const Vector &a, const Vector &b, MPI_Comm comm)
{
   real_t loc = a * b, glob = 0.0;
   MPI_Allreduce(&loc, &glob, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, comm);
   return glob;
}

static real_t Norm(const Vector &a, MPI_Comm comm) { return std::sqrt(Dot(a, a, comm)); }

// Independent dashpot assembly (test-local): face transformation from
// GetFaceElementTransformations(face), quadrature 2p+2, and the block as
// N^T K N with N = (phi (x) I_3) and K = Zp n n^T + Zs (I - n n^T).
static void IndependentDamping(NitscheWaveOperator &op, const Vector &v, Vector &Cv)
{
   ParMesh &mesh = op.GetMesh();
   auto &fes = op.Elasticity().GetFESpace();
   const NitscheMaterialSample &smp = op.Sample();
   const ParFiniteElementSpace *mfes = smp.RhoGF().ParFESpace();
   const int n = fes.GetNDofs();
   Cv.SetSize(v.Size());
   Cv = 0.0;
   Array<int> mdofs, dofs;
   for (int f : op.Elasticity().GetAbsorbingBdrFaces())
   {
      FaceElementTransformations *Tr = mesh.GetFaceElementTransformations(f);
      const int e = Tr->Elem1No;
      const FiniteElement &fe = *fes.GetFE(e);
      const int nd = fe.GetDof();
      mfes->GetElementDofs(e, mdofs);
      const real_t lam = smp.LambdaGF()(mdofs[0]), mu = smp.MuGF()(mdofs[0]);
      const real_t rho = smp.RhoGF()(mdofs[0]);
      const real_t Zp = rho * std::sqrt((lam + 2 * mu) / rho);
      const real_t Zs = rho * std::sqrt(mu / rho);
      const IntegrationRule &ir = IntRules.Get(Tr->GetGeometryType(),
                                               2 * fe.GetOrder() + 2);
      DenseMatrix Ce(3 * nd), N(3, 3 * nd), KN(3, 3 * nd), K(3);
      Ce = 0.0;
      Vector shape(nd), nor(3);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         Tr->SetAllIntPoints(&ip);
         fe.CalcShape(Tr->GetElement1IntPoint(), shape);
         CalcOrtho(Tr->Jacobian(), nor);
         const real_t nl = nor.Norml2();
         nor /= nl;
         N = 0.0;
         for (int c = 0; c < 3; c++)
            for (int k = 0; k < nd; k++) { N(c, c * nd + k) = shape(k); }
         for (int c = 0; c < 3; c++)
            for (int d = 0; d < 3; d++)
            {
               K(c, d) = Zp * nor(c) * nor(d) + Zs * ((c == d) - nor(c) * nor(d));
            }
         Mult(K, N, KN);
         AddMult_a_AtB(ip.weight * nl, N, KN, Ce);
      }
      fes.GetElementDofs(e, dofs);
      Vector x(3 * nd), y(3 * nd);
      for (int c = 0; c < 3; c++)
         for (int k = 0; k < nd; k++) { x(c * nd + k) = v(c * n + dofs[k]); }
      Ce.Mult(x, y);
      for (int c = 0; c < 3; c++)
         for (int k = 0; k < nd; k++) { Cv(c * n + dofs[k]) += y(c * nd + k); }
   }
}

static int AbsorbingBottomAndSides(const Vector &c, real_t Lz)
{
   return (std::abs(c(2)) < 1e-6 * Lz) ? 1 : 6;     // top natural, rest absorbing
}

// ------------------------------------------------------- C-02 damped inverse
static void TestDampedMass(int p)
{
   if (g_rank == 0) { std::cout << "\n=== C-02 (M + dt/2 C) inverse, p = " << p << " ===\n"; }
   std::unique_ptr<Mesh> m(BoxMeshAttr(3, 3, 2, 2 * kL, 2 * kL, kL,
      [](const Vector &c) { return AbsorbingBottomAndSides(c, kL); }));
   ParMesh pm(MPI_COMM_WORLD, *m);
   MPI_Comm comm = pm.GetComm();
   ConstantCoefficient lam_c(TPV102Params::lambda), mu_c(TPV102Params::mu);
   FunctionCoefficient rho_c(RhoField);
   MaterialField mat = MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c);
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.absorbing_attrs = {6};
   NitscheWaveOperator op(pm, p, mat, bc);
   const real_t dt = 0.05 * ComputeMaxDt(op.Sample(), pm, p, 1.0);
   op.SetTimeStep(dt);
   const NitscheMass *Md = op.DampedMass();
   int ndamped = Md ? Md->NumDampedElements() : 0;
   const int g_damped = GlobalSum(ndamped, comm);
   const int g_faces = GlobalSum(op.Elasticity().GetAbsorbingBdrFaces().Size(), comm);
   TEST_ASSERT(Md != nullptr && g_damped > 0, "damped mass built for absorbing elements");

   // Independent M (MFEM VectorMassIntegrator over the sampled rho).
   ParBilinearForm mform(&op.Elasticity().GetFESpace());
   mform.AddDomainIntegrator(new VectorMassIntegrator(
      const_cast<NitscheMaterialSample &>(op.Sample()).Rho()));
   mform.Assemble();
   mform.Finalize();

   Vector r(op.Height()), a, Ma(op.Height()), Ca, Cn;
   r.Randomize(91 + p + 100 * g_rank);
   r -= 0.5;
   Md->Apply(r, a);
   mform.SpMat().Mult(a, Ma);
   IndependentDamping(op, a, Ca);
   op.ApplyDamping(a, Cn);
   Vector res(Ma);
   res.Add(0.5 * dt, Ca);
   res -= r;
   const real_t rel = Norm(res, comm) / Norm(r, comm);
   Vector dC(Ca);
   dC -= Cn;
   const real_t relC = Norm(dC, comm) / Norm(Ca, comm);
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6) << "  absorbing faces = " << g_faces
                << ", damped elements = " << g_damped << ", dt = " << dt
                << " s\n  ||(M + dt/2 C) Apply(r) - r|| / ||r|| = " << rel
                << "\n  N-06 C vs independent C: ||dC a|| / ||C a|| = " << relC
                << "\n";
   }
   TEST_ASSERT(rel < 1e-12, "(M + dt/2 C) * Apply(r) == r to round-off");
   TEST_ASSERT(relC < 1e-12, "N-06 dashpot block == independent assembly");
}


// ===========================================================================
// T-12: dynamics with the dashpot.
// ===========================================================================

static const real_t kCuse[4] = {0.0, 0.057, 0.11, 0.14};   // P3: 0.5 C_eig

static MaterialField TPV102Material()
{
   return MaterialField::MakeConstant(TPV102Params::lambda, TPV102Params::mu,
                                      TPV102Params::rho);
}

// Assembled density-weighted mass (independent of NitscheMass).
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

// Smooth displacement bump centred in the box, zero velocity.
static void BumpU(const Vector &x, Vector &u)
{
   const real_t r2 = x(0) * x(0) + x(1) * x(1) + (x(2) + 0.5 * kL) * (x(2) + 0.5 * kL);
   const real_t g = std::exp(-r2 / (8.0e3 * 8.0e3));
   u.SetSize(3);
   u(0) = 0.6 * g; u(1) = -0.3 * g; u(2) = 1.0 * g;
}

static std::unique_ptr<NitscheWaveOperator> MakeAbsorbingBoxOp(ParMesh &pm, int p)
{
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.absorbing_attrs = {6};
   return std::make_unique<NitscheWaveOperator>(pm, p, TPV102Material(), bc);
}

static void InitBump(NitscheWaveOperator &op, Vector &U, Vector &V)
{
   ParGridFunction g(&op.Elasticity().GetFESpace());
   VectorFunctionCoefficient c(3, BumpU);
   g.ProjectCoefficient(c);
   U = g;
   V.SetSize(U.Size());
   V = 0.0;
}

// E_h = 1/2 V.MV + 1/2 U.KU - dt^2/8 A.MA with A = M^{-1}(-K U - C V).
struct EnergyEval
{
   NitscheWaveOperator &op; const MassMatrix &M; real_t dt; MPI_Comm comm;
   real_t Eh(const Vector &U, const Vector &V) const
   {
      Vector KU, CV, MV, F, A, MA;
      op.Elasticity().ApplyStiffness(U, KU);
      op.ApplyDamping(V, CV);
      F = KU; F += CV; F.Neg();
      op.Mass().Apply(F, A);
      M.Mult(V, MV);
      M.Mult(A, MA);
      return 0.5 * Dot(V, MV, comm) + 0.5 * Dot(U, KU, comm)
             - dt * dt / 8.0 * Dot(A, MA, comm);
   }
   real_t VCV(const Vector &V) const
   {
      Vector CV;
      op.ApplyDamping(V, CV);
      return Dot(V, CV, comm);
   }
};

// Run n steps with MFEM's CentralDifferenceSolver, tracking E_h and the
// trapezoidal dissipation sum.  Returns final U, V.
struct DampedRun
{
   real_t Eh0 = 0, EhN = 0, diss = 0, max_up = 0;   // max_up: max E_h increase
};

static DampedRun RunDamped(ParMesh &pm, int p, real_t dt, int n, Vector &U,
                           Vector &V, bool track)
{
   auto op = MakeAbsorbingBoxOp(pm, p);
   op->SetTimeStep(dt);
   MassMatrix M;
   if (track) { M.Build(*op); }
   InitBump(*op, U, V);
   EnergyEval ev{*op, M, dt, pm.GetComm()};
   DampedRun r;
   real_t Eh_prev = 0, vcv_prev = 0;
   if (track) { r.Eh0 = Eh_prev = ev.Eh(U, V); vcv_prev = ev.VCV(V); }
   real_t t = 0.0;
   CentralDifferenceSolver ode;
   ode.Init(*op);
   for (int s = 0; s < n; ++s)
   {
      real_t dts = dt;
      ode.Step(U, V, t, dts);
      if (track)
      {
         const real_t Eh = ev.Eh(U, V), vcv = ev.VCV(V);
         r.diss += dt * 0.5 * (vcv_prev + vcv);
         r.max_up = std::max(r.max_up, Eh - Eh_prev);
         Eh_prev = Eh; vcv_prev = vcv;
      }
   }
   r.EhN = Eh_prev;
   return r;
}

// --------------------------------------------------------- T-12(b) and (d)
static void TestDampedEnergy(int p)
{
   if (g_rank == 0) { std::cout << "\n=== T-12(b,d) energy with dashpot, 3x3x2 box, p = " << p << " ===\n"; }
   std::unique_ptr<Mesh> m(BoxMeshAttr(3, 3, 2, 2 * kL, 2 * kL, kL,
      [](const Vector &c) { return AbsorbingBottomAndSides(c, kL); }));
   ParMesh pm(MPI_COMM_WORLD, *m);
   auto op0 = MakeAbsorbingBoxOp(pm, p);
   const real_t dt_unit = ComputeMaxDt(op0->Sample(), pm, p, 1.0);
   const real_t T = 8.0;                       // waves reach and leave walls
   // Base step C_use/4 (decision, Chunhui Zhao 2026-10-14): C_use itself is
   // pre-asymptotic for the dt-scaling of the balance residual.
   int n = (int) std::ceil(T / (0.25 * kCuse[p] * dt_unit));
   const real_t dt = T / n;
   Vector U, V;
   DampedRun r1 = RunDamped(pm, p, dt, n, U, V, true);
   DampedRun r2 = RunDamped(pm, p, 0.5 * dt, 2 * n, U, V, true);
   const real_t res1 = (r1.EhN - r1.Eh0) + r1.diss;
   const real_t res2 = (r2.EhN - r2.Eh0) + r2.diss;
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6) << "  dt = " << dt << " s, T = " << T
                << " s, " << n << " steps\n"
                << "  E_h: " << r1.Eh0 << " -> " << r1.EhN << " (E_h,N/E_h,0 = "
                << r1.EhN / r1.Eh0 << "), max per-step increase / E_h,0 = "
                << r1.max_up / r1.Eh0 << "\n"
                << "  dt/2 run: max per-step increase / E_h,0 = "
                << r2.max_up / r2.Eh0 << "\n"
                << "  balance residual (E_h,N - E_h,0 + sum dt/2 [VCV_n + VCV_n+1]) / E_h,0:"
                << " dt " << res1 / r1.Eh0 << ", dt/2 " << res2 / r2.Eh0
                << ", ratio " << res1 / res2 << "\n";
   }
   TEST_ASSERT(r1.max_up <= 1e-12 * r1.Eh0 && r2.max_up <= 1e-12 * r2.Eh0,
               "T-12(b) E_h decreases monotonically (step increase <= 1e-12 E_h,0)");
   TEST_ASSERT(r1.EhN < 0.5 * r1.Eh0, "T-12(b) the dashpot removes energy");
   TEST_ASSERT(std::abs(res1 / res2) > 3.5 && std::abs(res1 / res2) < 4.5,
               "T-12(d) balance residual falls by ~4 when dt halves");
}

// ------------------------------------------------------------------ T-12(c)
// Negative control: explicit old-velocity damping with the plain M^{-1}
// (test-local wrapper, not an operator option).
class OldVelocityOp : public SecondOrderTimeDependentOperator
{
public:
   explicit OldVelocityOp(NitscheWaveOperator &op)
      : SecondOrderTimeDependentOperator(op.Height(), 0.0), op_(op) {}
   void Mult(const Vector &, Vector &) const override { MFEM_ABORT("second order"); }
   void Mult(const Vector &x, const Vector &v, Vector &a) const override
   {
      op_.SetTime(GetTime());
      op_.Mult(x, v, a);
   }
   void ImplicitSolve(const real_t, const real_t, const Vector &x,
                      const Vector &v, Vector &k) override
   {
      Vector KU, CV, F;
      op_.Elasticity().ApplyStiffness(x, KU);
      op_.ApplyDamping(v, CV);
      F = KU; F += CV; F.Neg();
      op_.Mass().Apply(F, k);
   }
private:
   NitscheWaveOperator &op_;
};

static void TestDampedOrder(int p)
{
   if (g_rank == 0) { std::cout << "\n=== T-12(c) temporal order with dashpot, 3x3x2 box, p = " << p << " ===\n"; }
   std::unique_ptr<Mesh> m(BoxMeshAttr(3, 3, 2, 2 * kL, 2 * kL, kL,
      [](const Vector &c) { return AbsorbingBottomAndSides(c, kL); }));
   ParMesh pm(MPI_COMM_WORLD, *m);
   MPI_Comm comm = pm.GetComm();
   auto op0 = MakeAbsorbingBoxOp(pm, p);
   const real_t dt_unit = ComputeMaxDt(op0->Sample(), pm, p, 1.0);
   const real_t T = 4.0;
   // Base step C_use/4 (decision, Chunhui Zhao 2026-10-14).
   int n0 = (int) std::ceil(T / (0.25 * kCuse[p] * dt_unit));
   const real_t dt0 = T / n0;
   auto run = [&](int refine, bool old_velocity, Vector &U, Vector &V)
   {
      auto op = MakeAbsorbingBoxOp(pm, p);
      const real_t dt = dt0 / refine;
      op->SetTimeStep(dt);
      InitBump(*op, U, V);
      real_t t = 0.0;
      if (old_velocity)
      {
         OldVelocityOp ov(*op);
         RunCentralDifference(ov, U, V, t, dt, n0 * refine);
      }
      else
      {
         RunCentralDifference(*op, U, V, t, dt, n0 * refine);
      }
   };
   Vector Ur, Vr;
   run(16, false, Ur, Vr);
   real_t eU[2][3], eV[2][3];
   const int ref[3] = {1, 2, 4};
   for (int v = 0; v < 2; v++)
      for (int k = 0; k < 3; k++)
      {
         Vector U, V;
         run(ref[k], v == 1, U, V);
         U -= Ur; V -= Vr;
         eU[v][k] = Norm(U, comm);
         eV[v][k] = Norm(V, comm);
      }
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6) << "  dt0 = " << dt0 << " s, T = " << T << " s\n";
      for (int v = 0; v < 2; v++)
      {
         std::cout << "  " << (v == 0 ? "implicit-damping (M + dt/2 C)" : "old-velocity control   ")
                   << " U ratios " << eU[v][0] / eU[v][1] << ", " << eU[v][1] / eU[v][2]
                   << "; Udot ratios " << eV[v][0] / eV[v][1] << ", " << eV[v][1] / eV[v][2]
                   << "\n";
      }
   }
   auto in = [](real_t r, real_t lo, real_t hi) { return r >= lo && r <= hi; };
   TEST_ASSERT(in(eU[0][0] / eU[0][1], 3.5, 4.5) && in(eU[0][1] / eU[0][2], 3.5, 4.5)
               && in(eV[0][0] / eV[0][1], 3.5, 4.5) && in(eV[0][1] / eV[0][2], 3.5, 4.5),
               "T-12(c) second order with the (M + dt/2 C) step");
   // The control's O(dt) error comes from the absorbing faces only; it shows
   // in U (the interior O(dt^2) error dominates the Udot norm).
   TEST_ASSERT(in(eU[1][0] / eU[1][1], 1.6, 2.5) && in(eU[1][1] / eU[1][2], 1.6, 2.5),
               "T-12(c) old-velocity control is first order (U ratios ~2)");
}


// ------------------------------------------------------------------ T-12(a)
// Reflection off the absorbing bottom.  Method: run A has the absorbing
// bottom at z = -30 km and exact incident-only Dirichlet data on every other
// wall; reference run B is the same column three times deeper (z = -90 km)
// with exact incident data on every wall, so nothing reflects back into
// z > -30 km within the run.  D = U_A - U_B on A's mesh (elements matched by
// centroid, DOFs by nodal coordinates) is the reflected field.  D is fitted
// by least squares (L2, over a region near the bottom) to the analytic
// unit-amplitude reflected P and S pulses: the fit coefficients are the
// measured R_PP and R_PS.
struct Incident
{
   real_t th = 0;                // incidence angle from the vertical (x-z plane)
   bool pwave = true;
   real_t c = 0, sigma = 0;      // speed, Gaussian width [s]
   real_t x0[3] = {0, 0, 0};     // pulse centre at t = 0
   real_t xb[3] = {0, 0, -kL};   // where the centre hits the bottom
   real_t t_hit = 0;
   real_t s[3], d[3];            // slowness, polarisation
   real_t g(real_t tau) const { const real_t a = tau / sigma; return std::exp(-a * a); }
   real_t dg(real_t tau) const { const real_t a = tau / sigma; return -2 * a / sigma * std::exp(-a * a); }
   void U(const Vector &x, real_t t, Vector &u, bool vel) const
   {
      real_t sx = 0;
      for (int i = 0; i < 3; i++) { sx += s[i] * (x(i) - x0[i]); }
      const real_t a = vel ? dg(t - sx) : g(t - sx);
      u.SetSize(3);
      for (int i = 0; i < 3; i++) { u(i) = d[i] * a; }
   }
};

// Reflected unit pulse through xb at t_hit with slowness sr, polarisation dr.
static void ReflectedShape(const Incident &I, const real_t sr[3], const real_t dr[3],
                           const Vector &x, real_t t, Vector &u)
{
   real_t sx = 0;
   for (int i = 0; i < 3; i++) { sx += sr[i] * (x(i) - I.xb[i]); }
   const real_t a = I.g(t - I.t_hit - sx);
   u.SetSize(3);
   for (int i = 0; i < 3; i++) { u(i) = dr[i] * a; }
}

// Analytic first-order (Lysmer-Kuhlemeyer) reflection of a plane P wave at
// angle th off the dashpot boundary: solves sigma_xz = Zs v_x, sigma_zz =
// Zp v_z for (R_PP, R_PS).  Boundary normal -z, medium above.
static void AnalyticRP(real_t th, real_t &Rpp, real_t &Rps)
{
   const real_t lam = TPV102Params::lambda, mu = TPV102Params::mu, rho = TPV102Params::rho;
   const real_t cp = std::sqrt((lam + 2 * mu) / rho), cs = std::sqrt(mu / rho);
   const real_t Zp = rho * cp, Zs = rho * cs;
   const real_t sphi = cs / cp * std::sin(th), cphi = std::sqrt(1 - sphi * sphi);
   struct W { real_t dx, dz, sx, sz; };
   const W wi{std::sin(th), -std::cos(th), std::sin(th) / cp, -std::cos(th) / cp};
   const W wp{std::sin(th), std::cos(th), std::sin(th) / cp, std::cos(th) / cp};
   const W ws{cphi, -sphi, sphi / cs, cphi / cs};
   auto e1 = [&](const W &w) { return -mu * (w.dx * w.sz + w.dz * w.sx) - Zs * w.dx; };
   auto e2 = [&](const W &w) { return -(lam * (w.dx * w.sx + w.dz * w.sz) + 2 * mu * w.dz * w.sz) - Zp * w.dz; };
   // [e1(wp) e1(ws); e2(wp) e2(ws)] [Rpp; Rps] = -[e1(wi); e2(wi)]
   const real_t a11 = e1(wp), a12 = e1(ws), a21 = e2(wp), a22 = e2(ws);
   const real_t b1 = -e1(wi), b2 = -e2(wi);
   const real_t det = a11 * a22 - a12 * a21;
   Rpp = (b1 * a22 - a12 * b2) / det;
   Rps = (a11 * b2 - a21 * b1) / det;
}

static Mesh *ReflectionMesh(int nx, int ny, int nz, real_t Lx, real_t Ly,
                            real_t Lz, bool absorbing_bottom, real_t xregion)
{
   Mesh *mesh = BoxMeshAttr(nx, ny, nz, Lx, Ly, Lz,
      [=](const Vector &c) { return (absorbing_bottom && std::abs(c(2) + Lz) < 1e-6 * Lz) ? 6 : 5; });
   for (int e = 0; e < mesh->GetNE(); e++)
   {
      Array<int> ev;
      mesh->GetElementVertices(e, ev);
      real_t c[3] = {0, 0, 0};
      for (int i = 0; i < ev.Size(); i++)
         for (int d = 0; d < 3; d++) { c[d] += mesh->GetVertex(ev[i])[d]; }
      for (int d = 0; d < 3; d++) { c[d] /= ev.Size(); }
      const bool region = (c[2] > -kL && c[2] < -10.0e3 && std::abs(c[0]) < xregion);
      mesh->SetAttribute(e, region ? 2 : 1);
   }
   mesh->SetAttributes();
   return mesh;
}

// U_B (on mesh B) mapped onto mesh A's vector DOFs; serial.
static void MapBtoA(NitscheWaveOperator &opA, NitscheWaveOperator &opB,
                    const Vector &UB, Vector &UBA)
{
   ParMesh &mA = opA.GetMesh(), &mB = opB.GetMesh();
   auto &fA = opA.Elasticity().GetFESpace();
   auto &fB = opB.Elasticity().GetFESpace();
   const int nA = fA.GetNDofs(), nB = fB.GetNDofs();
   auto key = [](const real_t c[3]) {
      return std::array<long long, 3>{std::llround(c[0]), std::llround(c[1]), std::llround(c[2])};
   };
   auto centroid = [](Mesh &m, int e, real_t c[3]) {
      Array<int> ev; m.GetElementVertices(e, ev);
      c[0] = c[1] = c[2] = 0;
      for (int i = 0; i < ev.Size(); i++)
         for (int d = 0; d < 3; d++) { c[d] += m.GetVertex(ev[i])[d]; }
      for (int d = 0; d < 3; d++) { c[d] /= ev.Size(); }
   };
   std::map<std::array<long long, 3>, int> bmap;
   for (int e = 0; e < mB.GetNE(); e++) { real_t c[3]; centroid(mB, e, c); bmap[key(c)] = e; }
   UBA.SetSize(3 * nA);
   UBA = 0.0;
   Array<int> dA, dB;
   for (int e = 0; e < mA.GetNE(); e++)
   {
      real_t c[3]; centroid(mA, e, c);
      auto it = bmap.find(key(c));
      MFEM_VERIFY(it != bmap.end(), "MapBtoA: element " << e << " not in mesh B");
      const int eb = it->second;
      const FiniteElement &fe = *fA.GetFE(e);
      const IntegrationRule &nodes = fe.GetNodes();
      ElementTransformation *TA = mA.GetElementTransformation(e);
      std::vector<Vector> pa(nodes.GetNPoints()), pb(nodes.GetNPoints());
      for (int k = 0; k < nodes.GetNPoints(); k++) { TA->Transform(nodes.IntPoint(k), pa[k]); }
      ElementTransformation *TB = mB.GetElementTransformation(eb);
      for (int k = 0; k < nodes.GetNPoints(); k++) { TB->Transform(nodes.IntPoint(k), pb[k]); }
      fA.GetElementDofs(e, dA);
      fB.GetElementDofs(eb, dB);
      for (int k = 0; k < nodes.GetNPoints(); k++)
      {
         int match = -1;
         for (int l = 0; l < nodes.GetNPoints(); l++)
         {
            Vector dd(pa[k]); dd -= pb[l];
            if (dd.Norml2() < 1e-3) { match = l; break; }
         }
         MFEM_VERIFY(match >= 0, "MapBtoA: node " << k << " of element " << e << " unmatched");
         for (int cc = 0; cc < 3; cc++) { UBA(cc * nA + dA[k]) = UB(cc * nB + dB[match]); }
      }
   }
}

static real_t RegionDot(ParFiniteElementSpace &fes, const Vector &a, const Vector &b)
{
   ParMesh *pm = fes.GetParMesh();
   Array<int> marker(pm->attributes.Max());
   marker = 0;
   marker[1] = 1;      // attribute 2 = measurement region
   ConstantCoefficient one(1.0);
   ParBilinearForm m(&fes);
   m.AddDomainIntegrator(new VectorMassIntegrator(one), marker);
   m.Assemble();
   m.Finalize();
   Vector Mb(b.Size());
   m.SpMat().Mult(b, Mb);
   return Dot(a, Mb, pm->GetComm());
}

struct ReflResult { real_t Rpp = 0, Rps = 0, dnorm = 0; };

static ReflResult MeasureReflection(const Incident &I, int p, int nx, int ny,
                                    int nz, real_t Lx, real_t Ly, real_t xregion)
{
   std::unique_ptr<Mesh> mA(ReflectionMesh(nx, ny, nz, Lx, Ly, kL, true, xregion));
   std::unique_ptr<Mesh> mB(ReflectionMesh(nx, ny, 3 * nz, Lx, Ly, 3 * kL, false, xregion));
   ParMesh pA(MPI_COMM_WORLD, *mA), pB(MPI_COMM_WORLD, *mB);
   auto uex = [I](const Vector &x, real_t t, Vector &u) { I.U(x, t, u, false); };
   auto vex = [I](const Vector &x, real_t t, Vector &v) { I.U(x, t, v, true); };
   auto make = [&](ParMesh &pm, bool absorbing) {
      BoundaryConfig bc;
      bc.fault_attr = 3;
      bc.dirichlet_attrs = {5};
      if (absorbing) { bc.absorbing_attrs = {6}; }
      bc.default_dirichlet_func = uex;
      return std::make_unique<NitscheWaveOperator>(pm, p, TPV102Material(), bc);
   };
   auto opA = make(pA, true), opB = make(pB, false);
   const real_t dtA = ComputeMaxDt(opA->Sample(), pA, p, kCuse[p]);
   const real_t dtB = ComputeMaxDt(opB->Sample(), pB, p, kCuse[p]);
   const real_t Tend = I.t_hit + 10.0e3 / I.c;
   const int n = (int) std::ceil(Tend / std::min(dtA, dtB));
   const real_t dt = Tend / n;
   auto run = [&](NitscheWaveOperator &op, Vector &U) {
      auto &fes = op.Elasticity().GetFESpace();
      ParGridFunction gu(&fes), gv(&fes);
      VectorFunctionCoefficient cu(3, uex), cv(3, vex);
      cu.SetTime(0.0); cv.SetTime(0.0);
      gu.ProjectCoefficient(cu); gv.ProjectCoefficient(cv);
      U = gu; Vector V(gv);
      op.SetTimeStep(dt);
      real_t t = 0.0;
      RunCentralDifference(op, U, V, t, dt, n);
   };
   Vector UA, UB, UBA;
   run(*opA, UA);
   run(*opB, UB);
   MapBtoA(*opA, *opB, UB, UBA);
   Vector D(UA);
   D -= UBA;

   // Analytic reflected shapes at Tend on A.
   const real_t lam = TPV102Params::lambda, mu = TPV102Params::mu, rho = TPV102Params::rho;
   const real_t cp = std::sqrt((lam + 2 * mu) / rho), cs = std::sqrt(mu / rho);
   real_t sP[3], dP[3], sS[3], dS[3];
   if (I.pwave)
   {
      const real_t sphi = cs / cp * std::sin(I.th), cphi = std::sqrt(1 - sphi * sphi);
      sP[0] = std::sin(I.th) / cp; sP[1] = 0; sP[2] = std::cos(I.th) / cp;
      dP[0] = std::sin(I.th); dP[1] = 0; dP[2] = std::cos(I.th);
      sS[0] = sphi / cs; sS[1] = 0; sS[2] = cphi / cs;
      dS[0] = cphi; dS[1] = 0; dS[2] = -sphi;
   }
   else   // normal-incidence S (x-polarised): reflected S only
   {
      sP[0] = 0; sP[1] = 0; sP[2] = 1.0 / cs; dP[0] = 1; dP[1] = 0; dP[2] = 0;
      sS[0] = 0; sS[1] = 0; sS[2] = 1.0 / cp; dS[0] = 0; dS[1] = 0; dS[2] = 1;
   }
   auto &fes = opA->Elasticity().GetFESpace();
   auto shape = [&](const real_t sr[3], const real_t dr[3], Vector &out) {
      ParGridFunction g(&fes);
      VectorFunctionCoefficient cf(3, [&I, sr, dr](const Vector &x, real_t t, Vector &u)
                                   { ReflectedShape(I, sr, dr, x, t, u); });
      cf.SetTime(dt * n);
      g.ProjectCoefficient(cf);
      out = g;
   };
   Vector PhiP, PhiS;
   shape(sP, dP, PhiP);
   shape(sS, dS, PhiS);
   const real_t gpp = RegionDot(fes, PhiP, PhiP), gss = RegionDot(fes, PhiS, PhiS);
   const real_t gps = RegionDot(fes, PhiP, PhiS);
   const real_t bp = RegionDot(fes, D, PhiP), bs = RegionDot(fes, D, PhiS);
   const real_t det = gpp * gss - gps * gps;
   ReflResult r;
   r.Rpp = (bp * gss - gps * bs) / det;
   r.Rps = (gpp * bs - gps * bp) / det;
   r.dnorm = std::sqrt(RegionDot(fes, D, D) / gpp);
   return r;
}

static Incident MakeIncident(real_t th_deg, bool pwave)
{
   const real_t lam = TPV102Params::lambda, mu = TPV102Params::mu, rho = TPV102Params::rho;
   Incident I;
   I.th = th_deg * M_PI / 180.0;
   I.pwave = pwave;
   I.c = pwave ? std::sqrt((lam + 2 * mu) / rho) : std::sqrt(mu / rho);
   I.sigma = pwave ? 1.2 : 2.0;                 // P3's widened pulses
   const real_t path = 22.0e3 / std::cos(I.th); // from z = -8 km to the bottom
   I.x0[0] = -22.0e3 * std::tan(I.th); I.x0[1] = 0; I.x0[2] = -8.0e3;
   I.t_hit = path / I.c;
   if (pwave)
   {
      I.s[0] = std::sin(I.th) / I.c; I.s[1] = 0; I.s[2] = -std::cos(I.th) / I.c;
      I.d[0] = std::sin(I.th); I.d[1] = 0; I.d[2] = -std::cos(I.th);
   }
   else
   {
      I.s[0] = 0; I.s[1] = 0; I.s[2] = -1.0 / I.c;
      I.d[0] = 1; I.d[1] = 0; I.d[2] = 0;
   }
   return I;
}

static void TestReflection()
{
   if (Mpi::WorldSize() > 1)
   {
      if (g_rank == 0) { std::cout << "\n  T-12(a) runs serially only (cross-mesh DOF map); skipped\n"; }
      return;
   }
   std::cout << "\n=== T-12(a) reflection off the dashpot (p = 2) ===\n"
             << "  method: A (absorbing bottom at -30 km) minus deep reference B\n"
                "  (-90 km, exact incident data on all walls), fitted to analytic\n"
                "  reflected P and S pulses over z in (-30, -10) km\n";
   // Normal incidence, column 10 x 10 km, two resolutions.
   for (int w = 0; w < 2; w++)
   {
      const bool pw = (w == 0);
      const Incident I = MakeIncident(0.0, pw);
      ReflResult r10 = MeasureReflection(I, 2, 1, 1, 3, 10.0e3, 10.0e3, 1e9);
      ReflResult r5 = MeasureReflection(I, 2, 2, 2, 6, 10.0e3, 10.0e3, 1e9);
      std::cout << std::setprecision(5) << "  normal " << (pw ? "P" : "S")
                << ": |D|/|pulse| = " << r10.dnorm << " (10 km hexes), "
                << r5.dnorm << " (5 km hexes); fitted R_same = " << r10.Rpp
                << ", " << r5.Rpp << " (analytic 0)\n";
      TEST_ASSERT(std::abs(r5.Rpp) < 0.05 && r5.dnorm < r10.dnorm,
                  "T-12(a) normal-incidence reflection small and falling with resolution");
   }
   // Oblique P, 60 x 10 x 30 km, 5 km hexes: unasserted diagnostics
   // (decision, Chunhui Zhao 2026-10-14).  Not cleanly measurable in this
   // geometry: the incident-only Dirichlet data on the lateral walls turns the
   // bottom corners into diffraction sources whose waves reach the
   // measurement region before the measurement time.
   std::cout << "  oblique P (NOT cleanly measurable in this geometry: incident-only\n"
                "  lateral-wall data makes the bottom corners diffraction sources;\n"
                "  printed for the record, not asserted):\n";
   for (real_t th : {30.0, 45.0, 60.0})
   {
      const Incident I = MakeIncident(th, true);
      ReflResult r = MeasureReflection(I, 2, 12, 2, 6, 60.0e3, 10.0e3, 15.0e3);
      real_t Rpp, Rps;
      AnalyticRP(I.th, Rpp, Rps);
      std::cout << std::setprecision(5) << "    P at " << th << " deg: measured R_PP = "
                << r.Rpp << ", R_PS = " << r.Rps << "; analytic first-order R_PP = "
                << Rpp << ", R_PS = " << Rps << "\n";
   }
}

// ---------------------------------------------------------------------- T-03
static void TestT03(ParMesh &pmesh)
{
   if (g_rank == 0) { std::cout << "\n=== T-03 absorbing face classification ===\n"; }
   MPI_Comm comm = pmesh.GetComm();
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.dirichlet_attrs = {5};
   bc.absorbing_attrs = {6};
   bc.default_dirichlet_func = [](const Vector &, real_t, Vector &u)
   { u.SetSize(3); u = 0.0; };
   LinearElastic model(TPV102Params::lambda, TPV102Params::mu);
   ElasticityDomainOperator<ParMesh> op(pmesh, 2, model, 0.0, kL, 2.0 * kL,
                                        bc, DGMethod::IP, SolverType::CG_AMG);

   // Ground truth: boundary elements tagged 6.
   std::set<int> tagged;
   for (int be = 0; be < pmesh.GetNBE(); be++)
   {
      if (pmesh.GetBdrAttribute(be) == 6)
      {
         tagged.insert(pmesh.GetBdrElementFaceIndex(be));
      }
   }
   const Array<int> &abs_faces = op.GetAbsorbingBdrFaces();
   std::set<int> abs_set(abs_faces.begin(), abs_faces.end());
   const Array<int> &fault = op.GetFaultInteriorFaces();
   int in_fault = 0, not_boundary = 0, bad_attr = 0;
   for (int f : abs_faces)
   {
      for (int g : fault) { if (g == f) { in_fault++; } }
      if (!pmesh.GetFaceInformation(f).IsBoundary()) { not_boundary++; }
   }
   for (int be : op.GetAbsorbingBdrElements())
   {
      const int a = pmesh.GetBdrAttribute(be);
      if (a != 6 || bc.dirichlet_attrs.count(a) || a == bc.fault_attr) { bad_attr++; }
   }
   const int g_tagged = GlobalSum((int) tagged.size(), comm);
   const int g_abs = GlobalSum(abs_faces.Size(), comm);
   const int g_same = GlobalSum(abs_set == tagged ? 0 : 1, comm);
   const int g_fault = GlobalSum(op.GetFaultInteriorFaces().Size()
                                 + op.GetFaultSharedFaces().Size(), comm);
   in_fault = GlobalSum(in_fault, comm);
   not_boundary = GlobalSum(not_boundary, comm);
   bad_attr = GlobalSum(bad_attr, comm);
   if (g_rank == 0)
   {
      std::cout << "  attr-6 boundary faces = " << g_tagged
                << ", absorbing faces collected = " << g_abs
                << ", fault faces = " << g_fault << "\n";
   }
   TEST_ASSERT(g_tagged > 0 && g_abs == g_tagged && g_same == 0,
               "absorbing face list == attr-6 boundary faces");
   TEST_ASSERT(in_fault == 0, "no absorbing face is a fault face");
   TEST_ASSERT(not_boundary == 0, "every absorbing face is an exterior boundary face");
   TEST_ASSERT(bad_attr == 0, "absorbing attr in neither the Dirichlet nor the fault set");
   TEST_ASSERT(g_fault > 0, "fault faces still classified alongside absorbing faces");
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   if (g_rank == 0)
   {
      std::cout << "=== Nitsche P4 (dashpot absorbing boundary) on "
                << Mpi::WorldSize() << " rank(s) ===\n";
   }
   {
      std::unique_ptr<Mesh> m(FaultBoxAbsorbingBottom(4, 4, 2));
      ParMesh pm(MPI_COMM_WORLD, *m);
      TestT03(pm);
   }
   for (int p = 1; p <= 3; p++) { TestDampedMass(p); }
   const bool t12 = (argc > 1 && std::string(argv[1]) == "--t12");
   if (t12)
   {
      for (int p = 1; p <= 3; p++) { TestDampedEnergy(p); }
      for (int p = 1; p <= 3; p++) { TestDampedOrder(p); }
   }
   const bool t12a = (argc > 1 && std::string(argv[1]) == "--t12a");
   if (t12a) { TestReflection(); }
   if (g_rank == 0)
   {
      std::cout << "\n=== Summary: " << num_passed << "/" << num_tests
                << " passed, " << num_failed << " failed ===\n";
   }
   return (num_failed == 0) ? 0 : 1;
}
