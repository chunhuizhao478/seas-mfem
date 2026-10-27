// Nitsche dynamic rupture, phase P3 (N-01 fault-free NitscheWaveOperator,
// N-02 ComputeMaxDt + CentralDifferenceSolver loop).  Fault-free tet boxes,
// every outer boundary Dirichlet with zero data (F = 0), TPV102 material.
//
//   T-09(c): CFL.  Largest stable C in dt = C h_min/(p^2 c_p), by bisection
//            on blow-up over 1000 steps (1% resolution), beside
//            2/omega_max (power iteration, omega_max^2 = max eig M^{-1}K).
//            3x3x2-hex box (432 tets), p = 1..4.  Power iteration also on
//            the 6x6x3 box (5184 tets), p = 1..3.
//   T-09(d): Verlet pin.  Hand loop mirroring NewmarkSolver::Step
//            (beta = 0, gamma = 1/2) vs CentralDifferenceSolver, 1000 steps,
//            bitwise.  6x6x3 box, p = 1..3.
//   T-10:    modified energy E_h = E - dt^2/8 A^T M A conserved to 1e-10
//            over 10^4 steps.  3x3x2 box, p = 1..3.
//   T-09(b): temporal order: dt, dt/2, dt/4 vs dt/16; error ratios in
//            [3.5, 4.5].  3x3x2 box, p = 1..3.
//   ComputeMaxDt: equals a serial-mesh reference bitwise (np-independent).
//
// Usage: mpirun -np 1 ./seas_test_nitsche_p3            (base suite)
//        mpirun -np 1 ./seas_test_nitsche_p3 --t18      (T-18 only)
//        mpirun -np 1 ./seas_test_nitsche_p3 --waves    (T-09(a), T-18:
//                                                       time-dependent
//                                                       Dirichlet data)
//        mpirun -np 2 ./seas_test_nitsche_p3 --np2      (T-09(d), T-10,
//                                                       ComputeMaxDt)
// The time step for T-09(d)/T-10 is C_use * h_min/(p^2 c_p) with C_use =
// 0.5 * C_eig rounded to 2 significant digits (np-independent).
//
// Scope (Chunhui Zhao, 2026-10-13): tetrahedral meshes only.

#include "mfem.hpp"
#include "../../domain/nitsche_wave_operator.hpp"
#include "../../solver/nitsche_cfl.hpp"
#include "../../config/tpv102_params.hpp"

#include <chrono>
#include <functional>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>
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

static real_t Dot(const Vector &a, const Vector &b, MPI_Comm comm)
{
   real_t loc = a * b, glob = 0.0;
   MPI_Allreduce(&loc, &glob, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, comm);
   return glob;
}

static real_t Norm(const Vector &a, MPI_Comm comm)
{
   return std::sqrt(Dot(a, a, comm));
}

static double Now()
{
   return std::chrono::duration<double>(
      std::chrono::steady_clock::now().time_since_epoch()).count();
}

// Fault-free 24-tet box: every boundary element attr 5 (Dirichlet).
static Mesh *BoxMesh(int nx, int ny, int nz)
{
   Mesh *mesh = new Mesh(Mesh::MakeCartesian3DWith24TetsPerHex(
      nx, ny, nz, 2*kL, 2*kL, kL));
   for (int v = 0; v < mesh->GetNV(); v++)
   {
      real_t *c = mesh->GetVertex(v);
      c[0] -= kL; c[1] -= kL; c[2] -= kL;
   }
   for (int b = 0; b < mesh->GetNBE(); b++) { mesh->SetBdrAttribute(b, 5); }
   mesh->SetAttributes();
   return mesh;
}

static BoundaryConfig MakeBC()
{
   BoundaryConfig bc;
   bc.fault_attr = 3;           // absent from the mesh
   bc.dirichlet_attrs = {5};
   bc.default_dirichlet_func = [](const Vector &, real_t, Vector &u)
   { u.SetSize(3); u = 0.0; };
   return bc;
}

static MaterialField TPV102Material()
{
   return MaterialField::MakeConstant(TPV102Params::lambda, TPV102Params::mu,
                                      TPV102Params::rho);
}

static std::unique_ptr<NitscheWaveOperator> MakeOp(ParMesh &pmesh, int p,
                                                    double *t_setup = nullptr,
                                                    double *t_first = nullptr,
                                                    double *t_matvec = nullptr)
{
   MPI_Barrier(pmesh.GetComm());
   double t0 = Now();
   auto op = std::make_unique<NitscheWaveOperator>(pmesh, p, TPV102Material(),
                                                   MakeBC());
   MPI_Barrier(pmesh.GetComm());
   double t1 = Now();
   Vector x(op->Height()), y;
   x = 1.0;
   op->Elasticity().ApplyStiffness(x, y);   // assembly + solver objects
   MPI_Barrier(pmesh.GetComm());
   double t2 = Now();
   op->Elasticity().ApplyStiffness(x, y);   // pure mat-vec
   MPI_Barrier(pmesh.GetComm());
   double t3 = Now();
   if (t_setup) { *t_setup = t1 - t0; }
   if (t_first) { *t_first = t2 - t1; }
   if (t_matvec) { *t_matvec = t3 - t2; }
   return op;
}

// Assembled density-weighted mass (independent of NitscheMass), for energies.
struct MassMatrix
{
   std::unique_ptr<ParBilinearForm> form;
   void Build(NitscheWaveOperator &op)
   {
      auto &fes = op.Elasticity().GetFESpace();
      form = std::make_unique<ParBilinearForm>(&fes);
      form->AddDomainIntegrator(new VectorMassIntegrator(
         const_cast<NitscheMaterialSample &>(op.Sample()).Rho()));
      form->Assemble();
      form->Finalize();
   }
   void Mult(const Vector &x, Vector &y) const
   {
      y.SetSize(x.Size());
      form->SpMat().Mult(x, y);
   }
};

// Smooth initial data vanishing on the box boundary.
static void SmoothU(const Vector &x, Vector &u)
{
   const real_t s = std::sin(M_PI * (x(0) + kL) / (2 * kL))
                    * std::sin(M_PI * (x(1) + kL) / (2 * kL))
                    * std::sin(M_PI * (x(2) + kL) / kL);
   u.SetSize(3);
   u(0) = 1.0 * s; u(1) = 0.5 * s; u(2) = -0.3 * s;
}

static void SmoothV(const Vector &x, Vector &v)
{
   const real_t s = std::sin(2 * M_PI * (x(0) + kL) / (2 * kL))
                    * std::sin(M_PI * (x(1) + kL) / (2 * kL))
                    * std::sin(M_PI * (x(2) + kL) / kL);
   v.SetSize(3);
   v(0) = -0.2 * s; v(1) = 0.1 * s; v(2) = 0.05 * s;
}

static void SmoothData(NitscheWaveOperator &op, Vector &U, Vector &V)
{
   ParGridFunction gu(&op.Elasticity().GetFESpace());
   ParGridFunction gv(&op.Elasticity().GetFESpace());
   VectorFunctionCoefficient cu(3, SmoothU), cv(3, SmoothV);
   gu.ProjectCoefficient(cu);
   gv.ProjectCoefficient(cv);
   U = gu;
   V = gv;
}

// omega_max^2 = max eig of M^{-1} K by power iteration (Rayleigh quotient
// v.Kv / v.Mv), converged to rel. change < tol.
static real_t PowerIteration(NitscheWaveOperator &op, const MassMatrix &M,
                             real_t tol, int max_it, int &iters)
{
   MPI_Comm comm = op.GetMesh().GetComm();
   const int n = op.Height();
   Vector v(n), Kv(n), Mv(n), w(n);
   v.Randomize(4242 + g_rank);
   v -= 0.5;
   v /= Norm(v, comm);
   real_t lam_old = 0.0, lam = 0.0;
   for (iters = 1; iters <= max_it; ++iters)
   {
      op.Elasticity().ApplyStiffness(v, Kv);
      M.Mult(v, Mv);
      lam = Dot(v, Kv, comm) / Dot(v, Mv, comm);
      if (iters > 1 && std::abs(lam - lam_old) < tol * std::abs(lam)) { break; }
      lam_old = lam;
      op.Mass().Apply(Kv, w);
      v = w;
      v /= Norm(v, comm);
   }
   return lam;
}

static int Sig2Round(real_t x, real_t &out)
{
   const real_t e = std::floor(std::log10(x));
   const real_t s = std::pow(10.0, e - 1);
   out = std::round(x / s) * s;
   return 0;
}

// Serial-mesh reference for ComputeMaxDt with constant material.
static real_t SerialDtRef(Mesh &smesh, int p, real_t cfl)
{
   const real_t cp = std::sqrt((TPV102Params::lambda + 2.0 * TPV102Params::mu)
                               / TPV102Params::rho);
   const real_t p2 = real_t(p) * real_t(p);
   real_t dmin = std::numeric_limits<real_t>::infinity();
   for (int e = 0; e < smesh.GetNE(); ++e)
   {
      const real_t vol = smesh.GetElementVolume(e);
      Array<int> vert;
      smesh.GetElementVertices(e, vert);
      Vector v0(smesh.GetVertex(vert[0]), 3), v1(smesh.GetVertex(vert[1]), 3);
      Vector v2(smesh.GetVertex(vert[2]), 3), v3(smesh.GetVertex(vert[3]), 3);
      auto tri_area = [](const Vector &a, const Vector &b,
                         const Vector &c) -> real_t {
         real_t e1[3] = {b(0)-a(0), b(1)-a(1), b(2)-a(2)};
         real_t e2[3] = {c(0)-a(0), c(1)-a(1), c(2)-a(2)};
         real_t cx = e1[1]*e2[2] - e1[2]*e2[1];
         real_t cy = e1[2]*e2[0] - e1[0]*e2[2];
         real_t cz = e1[0]*e2[1] - e1[1]*e2[0];
         return 0.5 * std::sqrt(cx*cx + cy*cy + cz*cz);
      };
      const real_t A = tri_area(v0, v1, v2) + tri_area(v0, v2, v3)
                       + tri_area(v0, v3, v1) + tri_area(v1, v2, v3);
      const real_t h = (A > 0) ? 6.0 * vol / A : std::pow(vol, 1.0 / 3.0);
      dmin = std::min(dmin, cfl * h / (p2 * cp));
   }
   return dmin;
}

// Run nsteps with dt from random U, V = 0; true if it blew up.
static bool BlowsUp(ParMesh &pmesh, int p, real_t dt, int nsteps)
{
   auto op = std::make_unique<NitscheWaveOperator>(pmesh, p, TPV102Material(),
                                                   MakeBC());
   op->SetTimeStep(dt);
   const int n = op->Height();
   Vector U(n), V(n);
   U.Randomize(777 + g_rank);
   U -= 0.5;
   V = 0.0;
   MPI_Comm comm = pmesh.GetComm();
   const real_t n0 = Norm(U, comm);
   real_t t = 0.0;
   CentralDifferenceSolver ode;
   ode.Init(*op);
   for (int s = 0; s < nsteps; ++s)
   {
      real_t dts = dt;
      ode.Step(U, V, t, dts);
      if (s % 50 == 49 || s == nsteps - 1)
      {
         const real_t nu = Norm(U, comm);
         if (!std::isfinite(nu) || nu > 1e3 * n0) { return true; }
      }
   }
   return false;
}

struct CflResult { real_t C_bis = 0, C_eig = 0, dt_unit = 0; int it = 0; };

static CflResult MeasureCFL(ParMesh &pmesh, Mesh &smesh, int p,
                            const char *mesh_name)
{
   CflResult r;
   double ts, tf, tm;
   auto op = MakeOp(pmesh, p, &ts, &tf, &tm);
   r.dt_unit = ComputeMaxDt(op->Sample(), pmesh, p, 1.0);
   const real_t dref = SerialDtRef(smesh, p, 1.0);
   MassMatrix M;
   M.Build(*op);
   const real_t w2 = PowerIteration(*op, M, 1e-6, 3000, r.it);
   r.C_eig = (2.0 / std::sqrt(w2)) / r.dt_unit;
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6) << "  [" << mesh_name << ", p=" << p
                << "] vsize=" << op->Height() << ", construct " << ts
                << " s, first K*x (assembly + AMG objects) " << tf
                << " s, K*x " << tm << " s\n"
                << std::setprecision(17) << "    dt_unit = h_min/(p^2 c_p) = "
                << r.dt_unit << " (serial ref " << dref << ")\n"
                << std::setprecision(8) << "    omega_max^2 = " << w2
                << " after " << r.it << " power iterations; C_eig = "
                << r.C_eig << "\n";
   }
   TEST_ASSERT(r.dt_unit == dref, "ComputeMaxDt == serial-mesh reference (bitwise)");
   return r;
}

// ------------------------------------------------------------------- T-09(c)
static void TestCFL(ParMesh &pmesh, Mesh &smesh, int p, real_t &C_bis_out,
                    real_t &C_eig_out)
{
   if (g_rank == 0) { std::cout << "\n=== T-09(c) CFL, 3x3x2 box, p = " << p << " ===\n"; }
   CflResult r = MeasureCFL(pmesh, smesh, p, "3x3x2");
   real_t lo = 0.8 * r.C_eig, hi = 1.2 * r.C_eig;
   const int nsteps = 1000;
   int guard = 0;
   while (BlowsUp(pmesh, p, lo * r.dt_unit, nsteps) && guard++ < 10) { lo *= 0.7; }
   guard = 0;
   while (!BlowsUp(pmesh, p, hi * r.dt_unit, nsteps) && guard++ < 10) { hi *= 1.3; }
   while ((hi - lo) / hi > 0.01)
   {
      const real_t mid = 0.5 * (lo + hi);
      if (BlowsUp(pmesh, p, mid * r.dt_unit, nsteps)) { hi = mid; } else { lo = mid; }
   }
   r.C_bis = lo;
   const real_t rel = std::abs(r.C_bis - r.C_eig) / r.C_eig;
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6) << "    bisection: largest stable C = "
                << lo << " (first unstable " << hi << "), C_eig = " << r.C_eig
                << ", |C_bis - C_eig|/C_eig = " << rel << "\n";
   }
   TEST_ASSERT(rel <= 0.0101, "C_CFL bisection agrees with 2/omega_max to 1%");
   C_bis_out = r.C_bis;
   C_eig_out = r.C_eig;
}

// ------------------------------------------------------------------- T-09(d)
static void TestVerletPin(ParMesh &pmesh, int p, real_t C_use)
{
   if (g_rank == 0) { std::cout << "\n=== T-09(d) Verlet pin, 6x6x3 box, p = " << p << " ===\n"; }
   auto op_a = MakeOp(pmesh, p);
   auto op_b = MakeOp(pmesh, p);
   const real_t dt = C_use * ComputeMaxDt(op_a->Sample(), pmesh, p, 1.0);
   op_a->SetTimeStep(dt);
   op_b->SetTimeStep(dt);
   const int n = op_a->Height(), nsteps = 1000;
   Vector U0(n), V0(n);
   SmoothData(*op_a, U0, V0);

   // (1) MFEM's CentralDifferenceSolver via the N-02 loop helper.
   Vector Ua(U0), Va(V0);
   real_t ta = 0.0;
   RunCentralDifference(*op_a, Ua, Va, ta, dt, nsteps);

   // (2) Hand loop mirroring NewmarkSolver::Step (linalg/ode.cpp) with
   //     beta = 0, gamma = 1/2, operation for operation.
   const real_t beta = 0.0, gamma = 0.5;
   const real_t fac0 = 0.5 - beta, fac2 = 1.0 - gamma;
   const real_t fac3 = beta, fac4 = gamma;
   Vector x(U0), v(V0), a(n);
   real_t t = 0.0;
   op_b->Mult(x, v, a);                       // A^0, once
   for (int s = 0; s < nsteps; ++s)
   {
      op_b->SetTime(t + dt);
      x.Add(dt, v);
      x.Add(fac0*dt*dt, a);
      v.Add(fac2*dt, a);
      op_b->SetTime(t + dt);
      op_b->ImplicitSolve(fac3*dt*dt, fac4*dt, x, v, a);
      x.Add(fac3*dt*dt, a);
      v.Add(fac4*dt, a);
      t += dt;
   }

   int du = 0, dv = 0;
   for (int i = 0; i < n; i++)
   {
      if (std::memcmp(&Ua(i), &x(i), sizeof(real_t)) != 0) { du++; }
      if (std::memcmp(&Va(i), &v(i), sizeof(real_t)) != 0) { dv++; }
   }
   int d[2] = {du, dv};
   MPI_Allreduce(MPI_IN_PLACE, d, 2, MPI_INT, MPI_SUM, pmesh.GetComm());
   const real_t nU = Norm(Ua, pmesh.GetComm());
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6) << "  dt = " << dt << " s (C = "
                << C_use << "), ||U(T)|| = " << nU << ", t_end equal: "
                << (ta == t) << "\n  differing entries: U " << d[0]
                << ", Udot " << d[1] << "\n";
   }
   TEST_ASSERT(d[0] == 0 && d[1] == 0 && ta == t,
               "hand Newmark loop == CentralDifferenceSolver bitwise (1000 steps)");
}

// ---------------------------------------------------------------------- T-10
static void TestEnergy(ParMesh &pmesh, int p, real_t C_use)
{
   if (g_rank == 0) { std::cout << "\n=== T-10 modified energy, 3x3x2 box, p = " << p << " ===\n"; }
   MPI_Comm comm = pmesh.GetComm();
   auto op = MakeOp(pmesh, p);
   const real_t dt = C_use * ComputeMaxDt(op->Sample(), pmesh, p, 1.0);
   op->SetTimeStep(dt);
   MassMatrix M;
   M.Build(*op);
   const int n = op->Height(), nsteps = 10000;
   Vector U(n), V(n), KU(n), MV(n), A(n);
   SmoothData(*op, U, V);

   auto energies = [&](real_t &E, real_t &Eh, real_t &AtMA)
   {
      op->Elasticity().ApplyStiffness(U, KU);
      op->Mass().Apply(KU, A);           // A = M^{-1} K U = -acceleration
      M.Mult(V, MV);
      const real_t kin = 0.5 * Dot(V, MV, comm);
      const real_t pot = 0.5 * Dot(U, KU, comm);
      AtMA = Dot(KU, A, comm);           // (KU)^T M^{-1} (KU) = A^T M A
      E = kin + pot;
      Eh = E - dt * dt / 8.0 * AtMA;
   };

   real_t E0, Eh0, a0;
   energies(E0, Eh0, a0);
   real_t max_dEh = 0, Emin = E0, Emax = E0, max_AtMA = a0, max_gap = 0;
   // Linear fit of E_n vs n.
   double sn = 0, sE = 0, snn = 0, snE = 0;
   real_t t = 0.0;
   CentralDifferenceSolver ode;
   ode.Init(*op);
   for (int s = 1; s <= nsteps; ++s)
   {
      real_t dts = dt;
      ode.Step(U, V, t, dts);
      real_t E, Eh, AtMA;
      energies(E, Eh, AtMA);
      max_dEh = std::max(max_dEh, std::abs(Eh - Eh0) / Eh0);
      Emin = std::min(Emin, E);
      Emax = std::max(Emax, E);
      max_AtMA = std::max(max_AtMA, AtMA);
      max_gap = std::max(max_gap, std::abs(E - Eh0));
      sn += s; sE += E; snn += double(s) * s; snE += double(s) * E;
   }
   const double N = nsteps;
   const double slope = (N * snE - sn * sE) / (N * snn - sn * sn);
   const real_t bound = dt * dt / 8.0 * max_AtMA;
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6) << "  dt = " << dt << " s (C = "
                << C_use << "), E_0 = " << E0 << ", E_h,0 = " << Eh0 << "\n"
                << "  max |E_h - E_h,0|/E_h,0 = " << max_dEh << "\n"
                << "  E oscillation (max - min)/E_0 = " << (Emax - Emin) / E0
                << ", dt^2/8 max A^T M A / E_0 = " << bound / E0
                << ", max |E - E_h,0| / E_0 = " << max_gap / E0 << "\n"
                << "  linear-fit drift of E: " << slope / E0
                << " E_0 per step\n";
   }
   TEST_ASSERT(max_dEh < 1e-10, "modified energy E_h conserved to 1e-10");
   TEST_ASSERT(max_gap <= bound * (1.0 + 1e-6) + 1e-10 * Eh0,
               "physical E within dt^2/8 max A^T M A of E_h,0");
   TEST_ASSERT(std::abs(slope) < 1e-10 * E0, "E drift slope < 1e-10 E_0 per step");
}

// ------------------------------------------------------------------- T-09(b)
static void TestTemporalOrder(ParMesh &pmesh, int p, real_t C_use)
{
   if (g_rank == 0) { std::cout << "\n=== T-09(b) temporal order, 3x3x2 box, p = " << p << " ===\n"; }
   MPI_Comm comm = pmesh.GetComm();
   auto op0 = MakeOp(pmesh, p);
   const real_t dt0 = C_use * ComputeMaxDt(op0->Sample(), pmesh, p, 1.0);
   const int n0 = 64;                      // steps at dt0; T = 64 dt0
   Vector U0(op0->Height()), V0(op0->Height());
   SmoothData(*op0, U0, V0);

   auto run = [&](int refine, Vector &U, Vector &V)
   {
      auto op = MakeOp(pmesh, p);
      const real_t dt = dt0 / refine;
      op->SetTimeStep(dt);
      U = U0; V = V0;
      real_t t = 0.0;
      RunCentralDifference(*op, U, V, t, dt, n0 * refine);
   };
   Vector Ur, Vr;
   run(16, Ur, Vr);
   real_t eU[3], eV[3];
   const int ref[3] = {1, 2, 4};
   for (int k = 0; k < 3; k++)
   {
      Vector U, V;
      run(ref[k], U, V);
      U -= Ur; V -= Vr;
      eU[k] = Norm(U, comm);
      eV[k] = Norm(V, comm);
   }
   const real_t rU1 = eU[0] / eU[1], rU2 = eU[1] / eU[2];
   const real_t rV1 = eV[0] / eV[1], rV2 = eV[1] / eV[2];
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6) << "  dt0 = " << dt0 << " s (C = "
                << C_use << "), T = " << n0 * dt0 << " s\n"
                << "  ||U err|| dt, dt/2, dt/4: " << eU[0] << ", " << eU[1]
                << ", " << eU[2] << "  ratios " << rU1 << ", " << rU2 << "\n"
                << "  ||V err|| dt, dt/2, dt/4: " << eV[0] << ", " << eV[1]
                << ", " << eV[2] << "  ratios " << rV1 << ", " << rV2 << "\n";
   }
   auto in = [](real_t r) { return r >= 3.5 && r <= 4.5; };
   TEST_ASSERT(in(rU1) && in(rU2), "U error ratios in [3.5, 4.5]");
   TEST_ASSERT(in(rV1) && in(rV2), "Udot error ratios in [3.5, 4.5]");
}


// ===========================================================================
// Plane-wave tests with time-dependent Dirichlet data (T-09(a), T-18).
// The exact solution is supplied as BoundaryConfig::default_dirichlet_func.
// ===========================================================================

// Box [-Lx/2, Lx/2] x [-Ly/2, Ly/2] x [-Lz, 0], 24 tets per hex, every
// boundary element attr 5.  If z_iface is finite, element attribute 1 below
// it and 2 above (it must be a plane of element faces).
static Mesh *BoxMeshExt(int nx, int ny, int nz, real_t Lx, real_t Ly,
                        real_t Lz, real_t z_iface)
{
   Mesh *mesh = new Mesh(Mesh::MakeCartesian3DWith24TetsPerHex(
      nx, ny, nz, Lx, Ly, Lz));
   for (int v = 0; v < mesh->GetNV(); v++)
   {
      real_t *c = mesh->GetVertex(v);
      c[0] -= 0.5 * Lx; c[1] -= 0.5 * Ly; c[2] -= Lz;
   }
   for (int b = 0; b < mesh->GetNBE(); b++) { mesh->SetBdrAttribute(b, 5); }
   if (std::isfinite(z_iface))
   {
      for (int e = 0; e < mesh->GetNE(); e++)
      {
         Array<int> ev;
         mesh->GetElementVertices(e, ev);
         real_t zc = 0.0;
         for (int i = 0; i < ev.Size(); i++) { zc += mesh->GetVertex(ev[i])[2]; }
         zc /= ev.Size();
         mesh->SetAttribute(e, (zc < z_iface) ? 1 : 2);
      }
   }
   mesh->SetAttributes();
   return mesh;
}

using ExactFn = std::function<void(const Vector &, real_t, Vector &)>;

static std::unique_ptr<NitscheWaveOperator> MakeWaveOp(
   ParMesh &pmesh, int p, const MaterialField &mat, const ExactFn &u_exact)
{
   BoundaryConfig bc;
   bc.fault_attr = 3;           // absent from the mesh
   bc.dirichlet_attrs = {5};
   bc.default_dirichlet_func = [u_exact](const Vector &x, real_t t, Vector &u)
   { u_exact(x, t, u); };
   return std::make_unique<NitscheWaveOperator>(pmesh, p, mat, bc);
}

// Run from the exact (u, du/dt) at t = 0 to T; return ||u_h - u||/||u|| (L2).
static real_t RunExact(NitscheWaveOperator &op, const ExactFn &u_exact,
                       const ExactFn &v_exact, real_t dt, int nsteps,
                       Vector *U_out = nullptr, real_t *E_ratio = nullptr,
                       const MassMatrix *M = nullptr)
{
   auto &fes = op.Elasticity().GetFESpace();
   ParGridFunction gu(&fes), gv(&fes);
   VectorFunctionCoefficient cu(3, u_exact), cv(3, v_exact);
   cu.SetTime(0.0); cv.SetTime(0.0);
   gu.ProjectCoefficient(cu);
   gv.ProjectCoefficient(cv);
   Vector U(gu), V(gv);
   op.SetTime(0.0);
   op.SetTimeStep(dt);
   MPI_Comm comm = op.GetMesh().GetComm();
   auto energy = [&]()
   {
      Vector KU, MV;
      op.Elasticity().ApplyStiffness(U, KU);
      M->Mult(V, MV);
      return 0.5 * Dot(V, MV, comm) + 0.5 * Dot(U, KU, comm);
   };
   real_t E0 = (M && E_ratio) ? energy() : 0.0, Emax = E0;
   real_t t = 0.0;
   CentralDifferenceSolver ode;
   ode.Init(op);
   for (int s = 0; s < nsteps; ++s)
   {
      real_t dts = dt;
      ode.Step(U, V, t, dts);
      if (M && E_ratio && (s % 20 == 19)) { Emax = std::max(Emax, energy()); }
   }
   if (M && E_ratio) { *E_ratio = Emax / E0; }
   gu = U;
   cu.SetTime(t);
   ParGridFunction zero(&fes);
   zero = 0.0;
   const real_t err = gu.ComputeL2Error(cu);
   const real_t nrm = zero.ComputeL2Error(cu);
   if (U_out) { *U_out = U; }
   return err / nrm;
}

// ------------------------------------------------------------------- T-09(a)
struct PlaneWave
{
   real_t k[3], d[3], omega, A;
   void U(const Vector &x, real_t t, Vector &u) const
   {
      const real_t ph = k[0]*x(0) + k[1]*x(1) + k[2]*x(2) - omega * t;
      u.SetSize(3);
      for (int i = 0; i < 3; i++) { u(i) = A * d[i] * std::sin(ph); }
   }
   void V(const Vector &x, real_t t, Vector &v) const
   {
      const real_t ph = k[0]*x(0) + k[1]*x(1) + k[2]*x(2) - omega * t;
      v.SetSize(3);
      for (int i = 0; i < 3; i++) { v(i) = -A * omega * d[i] * std::cos(ph); }
   }
};

static PlaneWave MakePlaneWave(bool pwave)
{
   const real_t cp = std::sqrt((TPV102Params::lambda + 2.0 * TPV102Params::mu)
                               / TPV102Params::rho);
   const real_t cs = std::sqrt(TPV102Params::mu / TPV102Params::rho);
   real_t n[3] = {1.0, 0.7, 0.4};
   const real_t nn = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
   for (auto &c : n) { c /= nn; }
   const real_t kmag = 2.0 * M_PI / 120.0e3;     // wavelength 120 km
   PlaneWave w;
   for (int i = 0; i < 3; i++) { w.k[i] = kmag * n[i]; }
   if (pwave)
   {
      for (int i = 0; i < 3; i++) { w.d[i] = n[i]; }
      w.omega = cp * kmag;
   }
   else
   {
      // d = normalize(n x e_z), perpendicular to n.
      real_t d[3] = {n[1], -n[0], 0.0};
      const real_t dn = std::sqrt(d[0]*d[0] + d[1]*d[1]);
      for (int i = 0; i < 3; i++) { w.d[i] = d[i] / dn; }
      w.omega = cs * kmag;
   }
   w.A = 1.0;
   return w;
}

static void TestPlaneWaveConvergence(int p, real_t C_use, bool pwave)
{
   const char *name = pwave ? "P" : "S";
   if (g_rank == 0)
   {
      std::cout << "\n=== T-09(a) plane " << name << " wave convergence, p = "
                << p << " ===\n";
   }
   const PlaneWave w = MakePlaneWave(pwave);
   ExactFn ue = [w](const Vector &x, real_t t, Vector &u) { w.U(x, t, u); };
   ExactFn ve = [w](const Vector &x, real_t t, Vector &v) { w.V(x, t, v); };
   const int nref[3][3] = {{2, 2, 1}, {4, 4, 2}, {8, 8, 4}};
   std::unique_ptr<ParMesh> pm[3];
   for (int m = 0; m < 3; m++)
   {
      std::unique_ptr<Mesh> sm(BoxMeshExt(nref[m][0], nref[m][1], nref[m][2],
                                          2 * kL, 2 * kL, kL, NAN));
      pm[m] = std::make_unique<ParMesh>(MPI_COMM_WORLD, *sm);
   }
   const real_t T = 2.0;
   // dt from the finest mesh, halved until the finest-mesh error changes by
   // < 1% between dt and dt/2 (time error negligible).
   auto op_f = MakeWaveOp(*pm[2], p, TPV102Material(), ue);
   real_t dt = C_use * ComputeMaxDt(op_f->Sample(), *pm[2], p, 1.0);
   int nsteps = (int) std::ceil(T / dt);
   dt = T / nsteps;
   real_t e_dt = RunExact(*op_f, ue, ve, dt, nsteps);
   real_t change = 1.0;
   int halvings = 0;
   for (; halvings < 4; ++halvings)
   {
      auto op_h = MakeWaveOp(*pm[2], p, TPV102Material(), ue);
      const real_t e_half = RunExact(*op_h, ue, ve, 0.5 * dt, 2 * nsteps);
      change = std::abs(e_dt - e_half) / e_half;
      if (g_rank == 0)
      {
         std::cout << std::setprecision(6) << "  finest mesh: dt = " << dt
                   << " s err " << e_dt << ", dt/2 err " << e_half
                   << ", relative change " << change << "\n";
      }
      if (change < 0.01) { break; }
      dt *= 0.5; nsteps *= 2; e_dt = e_half;
   }
   TEST_ASSERT(change < 0.01, "time error negligible at finest mesh (<1%)");
   real_t err[3];
   err[2] = e_dt;
   for (int m = 0; m < 2; m++)
   {
      auto op = MakeWaveOp(*pm[m], p, TPV102Material(), ue);
      err[m] = RunExact(*op, ue, ve, dt, nsteps);
   }
   const real_t r1 = std::log2(err[0] / err[1]), r2 = std::log2(err[1] / err[2]);
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6) << "  dt = " << dt << " s, T = " << T
                << " s, rel L2 error on 2x2x1 / 4x4x2 / 8x8x4 hex boxes: "
                << err[0] << ", " << err[1] << ", " << err[2]
                << "\n  observed L2 order: " << r1 << ", " << r2
                << " (expected ~" << p + 1 << ")\n";
   }
   TEST_ASSERT(r2 > p + 0.5, "plane-wave L2 order on the finest pair > p + 0.5");
}

// ---------------------------------------------------------------------- T-18
struct TwoLayer
{
   real_t z0;                 // interface
   real_t c1, c2, Z1, Z2;     // wave speeds and impedances of the mode
   real_t sigma, t0, A;       // Gaussian pulse g(tau) = exp(-((tau - t0)/sigma)^2)
   bool pwave;                // P: displacement along z; S: along x
   real_t R() const { return (Z1 - Z2) / (Z1 + Z2); }
   real_t Tc() const { return 2.0 * Z1 / (Z1 + Z2); }
   real_t g(real_t tau) const { const real_t a = (tau - t0) / sigma; return std::exp(-a * a); }
   real_t dg(real_t tau) const { const real_t a = (tau - t0) / sigma; return -2.0 * a / sigma * std::exp(-a * a); }
   // part: 0 = full exact, 1 = incident only, 2 = reflected shape (unit), 3 = transmitted shape (unit)
   real_t Scalar(real_t z, real_t t, int part, bool vel) const
   {
      auto G = [&](real_t tau) { return vel ? dg(tau) : g(tau); };
      const real_t inc = G(t - (z - z0) / c1), ref = G(t + (z - z0) / c1);
      const real_t tr = G(t - (z - z0) / c2);
      if (part == 1) { return (z < z0) ? A * inc : 0.0; }
      if (part == 2) { return (z < z0) ? A * ref : 0.0; }
      if (part == 3) { return (z < z0) ? 0.0 : A * tr; }
      return (z < z0) ? A * (inc + R() * ref) : A * Tc() * tr;
   }
   void Vec(const Vector &x, real_t t, Vector &u, int part, bool vel) const
   {
      u.SetSize(3); u = 0.0;
      u(pwave ? 2 : 0) = Scalar(x(2), t, part, vel);
   }
};

static MaterialField LayeredMaterial(real_t z0, real_t s,
                                     std::unique_ptr<FunctionCoefficient> c[3])
{
   // Layer 1 (z < z0): TPV102.  Layer 2: rho, cp, cs each x sqrt(s), so the
   // P and S impedance contrast Z2/Z1 is s.
   const real_t f = std::sqrt(s);
   const real_t rho1 = TPV102Params::rho, lam1 = TPV102Params::lambda,
                mu1 = TPV102Params::mu;
   const real_t rho2 = rho1 * f, mu2 = mu1 * f * f * f, lam2 = lam1 * f * f * f;
   c[0] = std::make_unique<FunctionCoefficient>(
      [=](const Vector &x) { return x(2) < z0 ? lam1 : lam2; });
   c[1] = std::make_unique<FunctionCoefficient>(
      [=](const Vector &x) { return x(2) < z0 ? mu1 : mu2; });
   c[2] = std::make_unique<FunctionCoefficient>(
      [=](const Vector &x) { return x(2) < z0 ? rho1 : rho2; });
   return MaterialField::MakeCoefficient(c[0].get(), c[1].get(), c[2].get());
}

static TwoLayer MakeTwoLayer(real_t z0, real_t s, bool pwave)
{
   const real_t f = std::sqrt(s);
   const real_t rho1 = TPV102Params::rho;
   const real_t cp1 = std::sqrt((TPV102Params::lambda + 2.0 * TPV102Params::mu) / rho1);
   const real_t cs1 = std::sqrt(TPV102Params::mu / rho1);
   TwoLayer L;
   L.z0 = z0;
   L.pwave = pwave;
   L.c1 = pwave ? cp1 : cs1;
   L.c2 = L.c1 * f;
   L.Z1 = rho1 * L.c1;
   L.Z2 = rho1 * f * L.c2;
   // Pulse widths doubled (P 1.2 s, S 2.0 s; about 7 km in space) so p = 1
   // is in its asymptotic regime on the 1x1x2 / 2x2x4 / 4x4x8 meshes.
   L.sigma = pwave ? 1.2 : 2.0;
   L.t0 = pwave ? 1.5 : 2.5;
   L.A = 1.0;
   return L;
}

// Inner product <a, b> restricted to elements of attribute `attr`.
static real_t LayerDot(ParFiniteElementSpace &fes, int attr, const Vector &a,
                       const Vector &b)
{
   ParMesh *pm = fes.GetParMesh();
   Array<int> marker(pm->attributes.Max());
   marker = 0;
   marker[attr - 1] = 1;
   ConstantCoefficient one(1.0);
   ParBilinearForm m(&fes);
   m.AddDomainIntegrator(new VectorMassIntegrator(one), marker);
   m.Assemble();
   m.Finalize();
   Vector Mb(b.Size());
   m.SpMat().Mult(b, Mb);
   return Dot(a, Mb, pm->GetComm());
}

struct T18Result { real_t err, R, T, E_ratio; };

static T18Result RunTwoLayer(ParMesh &pmesh, int p, real_t s, bool pwave,
                             real_t C_use, Vector *U_out = nullptr,
                             bool homogeneous_material = false,
                             int nsteps_override = -1)
{
   const real_t z0 = -0.5 * kL;
   const TwoLayer L = MakeTwoLayer(z0, s, pwave);
   std::unique_ptr<FunctionCoefficient> c[3];
   MaterialField mat = homogeneous_material ? TPV102Material()
                                            : LayeredMaterial(z0, s, c);
   ExactFn ue = [L](const Vector &x, real_t t, Vector &u) { L.Vec(x, t, u, 0, false); };
   ExactFn ve = [L](const Vector &x, real_t t, Vector &v) { L.Vec(x, t, v, 0, true); };
   auto op = MakeWaveOp(pmesh, p, mat, ue);
   MassMatrix M;
   M.Build(*op);
   real_t dt = C_use * ComputeMaxDt(op->Sample(), pmesh, p, 1.0);
   const real_t Tend = pwave ? 2.5 : 4.0;
   int nsteps = (int) std::ceil(Tend / dt);
   dt = Tend / nsteps;
   if (nsteps_override > 0) { nsteps = nsteps_override; }
   T18Result r;
   Vector U;
   r.err = RunExact(*op, ue, ve, dt, nsteps, &U, &r.E_ratio, &M);
   if (U_out) { *U_out = U; }
   const real_t tf = dt * nsteps;
   // Measured R and T: fit the numerical field to the analytic shapes.
   auto &fes = op->Elasticity().GetFESpace();
   auto proj = [&](int part, Vector &out)
   {
      ParGridFunction g(&fes);
      VectorFunctionCoefficient cf(3, [L, part](const Vector &x, real_t t, Vector &u)
                                   { L.Vec(x, t, u, part, false); });
      cf.SetTime(tf);
      g.ProjectCoefficient(cf);
      out = g;
   };
   Vector inc, ref, tr;
   proj(1, inc); proj(2, ref); proj(3, tr);
   Vector Ur(U);
   Ur -= inc;
   r.R = LayerDot(fes, 1, Ur, ref) / LayerDot(fes, 1, ref, ref);
   r.T = LayerDot(fes, 2, U, tr) / LayerDot(fes, 2, tr, tr);
   return r;
}

static void TestTwoLayer(real_t C_use[5])
{
   const real_t z0 = -0.5 * kL;
   const real_t Lh = 15.0e3;   // lateral box size (1-D problem)
   const int nref[3][3] = {{1, 1, 2}, {2, 2, 4}, {4, 4, 8}};
   std::unique_ptr<ParMesh> pm[3];
   for (int m = 0; m < 3; m++)
   {
      std::unique_ptr<Mesh> sm(BoxMeshExt(nref[m][0], nref[m][1], nref[m][2],
                                          Lh, Lh, kL, z0));
      pm[m] = std::make_unique<ParMesh>(MPI_COMM_WORLD, *sm);
   }
   if (g_rank == 0)
   {
      std::cout << "\n=== T-18 two-layer reflection/transmission ===\n"
                << "  box 15 x 15 x 30 km, interface z = -15 km; meshes "
                   "1x1x2 / 2x2x4 / 4x4x8 hexes (24 tets each)\n";
   }

   // (1) Error vs analytic and its order, contrast 2, P and S, p = 1..3.
   for (int w = 0; w < 2; w++)
   {
      const bool pwave = (w == 0);
      for (int p = 1; p <= 3; p++)
      {
         real_t e[3];
         T18Result rf{};
         for (int m = 0; m < 3; m++)
         {
            T18Result r = RunTwoLayer(*pm[m], p, 2.0, pwave, C_use[p]);
            e[m] = r.err;
            if (m == 2) { rf = r; }
         }
         const real_t r1 = std::log2(e[0] / e[1]), r2 = std::log2(e[1] / e[2]);
         const TwoLayer L = MakeTwoLayer(z0, 2.0, pwave);
         if (g_rank == 0)
         {
            std::cout << std::setprecision(6) << "  [" << (pwave ? "P" : "S")
                      << ", contrast 2, p=" << p << "] rel L2 error "
                      << e[0] << ", " << e[1] << ", " << e[2] << "  order "
                      << r1 << ", " << r2 << "\n    finest: R = " << rf.R
                      << " (exact " << L.R() << "), T = " << rf.T << " (exact "
                      << L.Tc() << ")\n";
         }
         TEST_ASSERT(r2 > p + 0.5, "two-layer L2 order on the finest pair > p + 0.5");
      }
   }

   // (2) R, T, stability at contrasts 1, 2, 4 (penalty_factor unchanged),
   //     middle mesh, p = 2.
   const real_t contrasts[3] = {1.0, 2.0, 4.0};
   for (int w = 0; w < 2; w++)
   {
      const bool pwave = (w == 0);
      for (real_t s : contrasts)
      {
         const TwoLayer L = MakeTwoLayer(z0, s, pwave);
         T18Result r = RunTwoLayer(*pm[1], 2, s, pwave, C_use[2]);
         const bool finite = std::isfinite(r.err) && std::isfinite(r.E_ratio);
         if (g_rank == 0)
         {
            std::cout << std::setprecision(6) << "  [" << (pwave ? "P" : "S")
                      << ", contrast " << s << ", 2x2x4, p=2] R = " << r.R
                      << " (exact " << L.R() << "), T = " << r.T << " (exact "
                      << L.Tc() << "), rel err " << r.err
                      << ", max E / E_0 = " << r.E_ratio << "\n";
         }
         TEST_ASSERT(finite && r.E_ratio < 10.0, "two-layer run stable, energy bounded");
         TEST_ASSERT(std::abs(r.R - L.R()) < 0.02 && std::abs(r.T - L.Tc()) < 0.02,
                     "measured R and T within 0.02 of analytic");
      }
   }

   // (3) Equal materials (layered field, contrast 1) == homogeneous run bitwise.
   {
      Vector Ul, Uh;
      RunTwoLayer(*pm[0], 2, 1.0, true, C_use[2], &Ul, false, 200);
      RunTwoLayer(*pm[0], 2, 1.0, true, C_use[2], &Uh, true, 200);
      int nd = 0;
      for (int i = 0; i < Ul.Size(); i++)
      {
         if (std::memcmp(&Ul(i), &Uh(i), sizeof(real_t)) != 0) { nd++; }
      }
      MPI_Allreduce(MPI_IN_PLACE, &nd, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      if (g_rank == 0)
      {
         std::cout << "  equal-material layered vs homogeneous (1x1x2, p=2, "
                      "200 steps): differing entries " << nd << "\n";
      }
      TEST_ASSERT(nd == 0, "equal materials == homogeneous run bitwise");
   }
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   const bool np2 = (argc > 1 && std::string(argv[1]) == "--np2");
   const bool waves = (argc > 1 && std::string(argv[1]) == "--waves");
   const bool t18 = (argc > 1 && std::string(argv[1]) == "--t18");
   if (g_rank == 0)
   {
      std::cout << "=== Nitsche P3 (fault-free explicit elastodynamics) on "
                << Mpi::WorldSize() << " rank(s)"
                << (np2 ? " [--np2 subset]" : "") << " ===\n";
   }

   std::unique_ptr<Mesh> small(BoxMesh(3, 3, 2));
   std::unique_ptr<Mesh> big(BoxMesh(6, 6, 3));
   ParMesh psmall(MPI_COMM_WORLD, *small);
   ParMesh pbig(MPI_COMM_WORLD, *big);

   // C_use(p) = 0.5 * C_eig on the small box, rounded to 2 sig. digits.
   real_t C_use[5] = {0, 0, 0, 0, 0};
   for (int p = 1; p <= 3; p++)
   {
      if (g_rank == 0) { std::cout << "\n=== ComputeMaxDt / omega_max, p = " << p << " ===\n"; }
      CflResult rs = MeasureCFL(psmall, *small, p, "3x3x2");
      CflResult rb = MeasureCFL(pbig, *big, p, "6x6x3");
      Sig2Round(0.5 * rs.C_eig, C_use[p]);
      if (g_rank == 0)
      {
         std::cout << std::setprecision(6) << "  C_eig 3x3x2 = " << rs.C_eig
                   << ", 6x6x3 = " << rb.C_eig << " (ratio "
                   << rb.C_eig / rs.C_eig << "); C_use = " << C_use[p] << "\n";
      }
   }

   if (t18)
   {
      TestTwoLayer(C_use);
      if (g_rank == 0)
      {
         std::cout << "\n=== Summary: " << num_passed << "/" << num_tests
                   << " passed, " << num_failed << " failed ===\n";
      }
      return (num_failed == 0) ? 0 : 1;
   }

   if (waves)
   {
      for (int p = 1; p <= 3; p++)
      {
         TestPlaneWaveConvergence(p, C_use[p], true);
         TestPlaneWaveConvergence(p, C_use[p], false);
      }
      TestTwoLayer(C_use);
      if (g_rank == 0)
      {
         std::cout << "\n=== Summary: " << num_passed << "/" << num_tests
                   << " passed, " << num_failed << " failed ===\n";
      }
      return (num_failed == 0) ? 0 : 1;
   }

   if (!np2)
   {
      real_t Cb[5] = {0}, Ce[5] = {0};
      for (int p = 1; p <= 4; p++) { TestCFL(psmall, *small, p, Cb[p], Ce[p]); }
      if (g_rank == 0)
      {
         std::cout << "\n  C_CFL table (3x3x2 box, dt = C h_min/(p^2 c_p)):\n"
                   << "    p   C_bisection   C_eig=2/omega_max\n";
         for (int p = 1; p <= 4; p++)
         {
            std::cout << "    " << p << "   " << std::setprecision(5) << Cb[p]
                      << "        " << Ce[p] << "\n";
         }
      }
   }

   for (int p = 1; p <= 3; p++) { TestVerletPin(pbig, p, C_use[p]); }
   for (int p = 1; p <= 3; p++) { TestEnergy(psmall, p, C_use[p]); }
   if (!np2)
   {
      for (int p = 1; p <= 3; p++) { TestTemporalOrder(psmall, p, C_use[p]); }
   }

   if (g_rank == 0)
   {
      std::cout << "\n=== Summary: " << num_passed << "/" << num_tests
                << " passed, " << num_failed << " failed ===\n";
   }
   return (num_failed == 0) ? 0 : 1;
}
