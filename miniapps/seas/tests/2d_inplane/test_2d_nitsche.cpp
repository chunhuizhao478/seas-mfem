// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// tests/2d_inplane/test_2d_nitsche.cpp — 2D in-plane design rev 3, phase C
// (part 1, C1): the Nitsche dynamic operator on triangle meshes (design
// section 5.5).  The 3D tests it follows are tests/verification/
// test_nitsche_p3.cpp (N1 to N3), test_nitsche_p4.cpp (N4) and
// test_nitsche_lsw_inline.cpp (N5).  Material TPV102 (= TPV205) throughout;
// the meshes are built here, each quadrilateral of a grid split into 4
// triangles through its centre (mirror symmetric about every grid line).
//
//   N1  plane P and S waves along (1, 0.7) off the grid axes, the exact
//       solution as time-dependent Dirichlet data on every wall, on the
//       4 x 4, 8 x 8 and 16 x 16 grids of [-30, 30]^2 km, p = 1, 2, 3, to
//       T = 2 s: the finest mesh's time error below 1% of its error, and the
//       observed L2 order on the finest pair above p + 0.5.
//   N2  the CFL constant on triangles (design section 5.5, derivation
//       10.3): C_eig = 2 / (omega_max dt_unit), omega_max^2 = max eig
//       M^{-1} K by power iteration and dt_unit = h_min / (p^2 c_p) with
//       h = 4 A / P, and the largest stable C by bisection on blow-up over
//       1000 steps (1% resolution); they agree to 1%, p = 1 to 4, on the
//       4 x 4 grid, and the power iteration on the 8 x 8 grid gives the same
//       C_eig.  The 2D column of NitscheStableConstant (solver/
//       nitsche_cfl.hpp) is the measured min(C_eig, bisection) rounded down;
//       printed against the 3D column.  Every run here steps at half of it.
//       N2M (opt-in): C_eig on the committed TPV205-2D and TPV105-2D meshes.
//   N3  the step: the modified energy E_h = E - dt^2/8 A^T M A conserved to
//       1e-10 over 10^4 steps and the physical E within its dt^2 bound
//       (4 x 4 grid); the hand Newmark loop equals CentralDifferenceSolver
//       bit for bit over 1000 steps (8 x 8 grid); temporal order from half
//       the step, error ratios in [3.5, 4.5]; ComputeMaxDt equals a serial
//       4 A / P reference bit for bit.  p = 1, 2, 3.
//   N4  the absorbing boundary: (M + dt/2 C) Apply(r) = r with M and C
//       formed independently (p = 1, 2, 3, variable density); the
//       reflection of plane pulses off the dashpot, measured as A (absorbing
//       bottom) minus B (three times as deep, exact incident data on every
//       wall) and fitted to the analytic reflected P and S pulses: normal
//       incidence P and S (fitted R_same small and |D| falling with
//       resolution), and P at 45 degrees against the analytic first-order
//       R_PP and R_PS (within 0.002 at both of two resolutions).  The 2D
//       box is wide enough that the incident-only data of the side walls
//       cannot reach the fit region before the fit time (the reason the 3D
//       test could only print its oblique values).
//   N5  the fault: an in-plane slip-weakening rupture (TPV205 friction,
//       tau_0 = 70 MPa, sigma_0 = 120 MPa, the nucleation dtau F(r) G(t) of
//       test_nitsche_lsw_inline.cpp at x = 0) on the 24 x 24 grid of
//       [-30, 30]^2 km with the fault on y = 0, |x| < 15 km (attr 3), the
//       walls absorbing, up = [0, 0, 1] (the slip in the strike slot):
//       Proposition 1 of the design (the out-of-plane component of d, V, the
//       trial traction, the prestress and the nucleation traction is
//       exactly 0 at every DOF and step, and so is the dip slot of X); the
//       mirror y -> -y (the response to tangential slip leaves the fault
//       normal stress alone: |sigma_eff - 120 MPa| < 1e-8 x 120 MPa at every
//       DOF and step) and the half turn (x, y) -> (-x, -y), which leaves the
//       stress, the mesh and the nucleation unchanged (slip, slip rate and
//       slip path at x and -x to 1e-9 of their largest values); the first DOFs to
//       slip inside the nucleation patch and inside its rise time; slip
//       above 1 m by 3 s.  The same problem rotated by 90 degrees with
//       up = [0, 1] (the slip in the dip slot) equals it to 1e-9; the fast
//       fault kernels, dense blocks and fused sweep equal the CSR products
//       and ComputeTractionDiagnostics to 1e-11; a run on more ranks equals
//       the 1-rank run to 1e-11 (on 2 and 4 ranks the partition puts the
//       two sides of the fault on different ranks, so every fault face is a
//       rank-boundary face).  Kelvin-Voigt (whole domain and a layer)
//       and Prakash-Clifton run finite; Prakash-Clifton with L = 1e-200 is
//       bit for bit the run without it.  p = 1, 2.
//
// Usage: mpirun -np 1 ./seas_test_2d_nitsche              (N1 to N5)
//        mpirun -np 2|4 ./seas_test_2d_nitsche --np       (N3, N5: the
//                                                          rank-invariant parts)
//        ./seas_test_2d_nitsche --only N2[,N4,...]        (some parts; N2M only
//                                                          this way)
//
// Scope (Chunhui Zhao, 2026-10-25).

#include "mfem.hpp"
#include "../../domain/nitsche_wave_operator.hpp"
#include "../../solver/nitsche_cfl.hpp"
#include "../../friction/slip_weakening_closure.hpp"
#include "../../config/tpv102_params.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_failed = 0;
static int g_rank = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { if (g_rank == 0) std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; if (g_rank == 0) std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static const real_t kL = 30.0e3;   // boxes [-L, L]^2 unless stated
static const real_t kLam = TPV102Params::lambda, kMu = TPV102Params::mu, kRho = TPV102Params::rho;
static real_t Cp() { return std::sqrt((kLam + 2.0 * kMu) / kRho); }
static real_t Cs() { return std::sqrt(kMu / kRho); }

/// The step of every run here: half the stable constant of N2.
static real_t CUse(int p) { return 0.5 * NitscheStableConstant(2, p); }

static real_t Dot(const Vector &a, const Vector &b, MPI_Comm comm)
{
   real_t loc = a * b, glob = 0.0;
   MPI_Allreduce(&loc, &glob, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, comm);
   return glob;
}
static real_t Norm(const Vector &a, MPI_Comm comm) { return std::sqrt(Dot(a, a, comm)); }
static int ISum(int v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_INT, MPI_SUM, c); return v; }
static real_t RMax(real_t v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, c); return v; }
static real_t RMin(real_t v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPITypeMap<real_t>::mpi_type, MPI_MIN, c); return v; }

static MaterialField Material() { return MaterialField::MakeConstant(kLam, kMu, kRho); }

// ------------------------------------------------------------------ meshes
// [x0, x0 + nx hx] x [y0, y0 + ny hy], every cell split into 4 triangles
// through its centre.  Outer edges get bdr_attr(edge midpoint); an interior
// grid edge whose midpoint satisfies fault(xm, ym) becomes a boundary
// segment of attribute 3.  With rotate the vertices are mapped to (-y, x)
// afterwards (exact in floating point; the triangles keep their order).
using AttrFn = std::function<int(real_t, real_t)>;
using FaultFn = std::function<bool(real_t, real_t)>;

static Mesh *CrossGrid(int nx, int ny, real_t x0, real_t y0, real_t hx, real_t hy,
                       const AttrFn &bdr_attr, const FaultFn &fault = nullptr,
                       bool rotate = false)
{
   const int ngv = (nx + 1) * (ny + 1);
   Mesh *m = new Mesh(2, ngv + nx * ny, 4 * nx * ny, 2 * (nx + ny), 2);
   auto gv = [&](int i, int j) { return i + (nx + 1) * j; };
   auto cv = [&](int i, int j) { return ngv + i + nx * j; };
   auto put = [&](real_t x, real_t y) { if (rotate) { m->AddVertex(-y, x); } else { m->AddVertex(x, y); } };
   for (int j = 0; j <= ny; j++)
      for (int i = 0; i <= nx; i++) { put(x0 + i * hx, y0 + j * hy); }
   for (int j = 0; j < ny; j++)
      for (int i = 0; i < nx; i++) { put(x0 + (i + 0.5) * hx, y0 + (j + 0.5) * hy); }
   for (int j = 0; j < ny; j++)
      for (int i = 0; i < nx; i++)
      {
         const int c = cv(i, j);
         m->AddTriangle(gv(i, j), gv(i + 1, j), c, 1);
         m->AddTriangle(gv(i + 1, j), gv(i + 1, j + 1), c, 1);
         m->AddTriangle(gv(i + 1, j + 1), gv(i, j + 1), c, 1);
         m->AddTriangle(gv(i, j + 1), gv(i, j), c, 1);
      }
   for (int i = 0; i < nx; i++)
   {
      const real_t xm = x0 + (i + 0.5) * hx;
      m->AddBdrSegment(gv(i, 0), gv(i + 1, 0), bdr_attr(xm, y0));
      m->AddBdrSegment(gv(i + 1, ny), gv(i, ny), bdr_attr(xm, y0 + ny * hy));
   }
   for (int j = 0; j < ny; j++)
   {
      const real_t ym = y0 + (j + 0.5) * hy;
      m->AddBdrSegment(gv(nx, j), gv(nx, j + 1), bdr_attr(x0 + nx * hx, ym));
      m->AddBdrSegment(gv(0, j + 1), gv(0, j), bdr_attr(x0, ym));
   }
   if (fault)
   {
      for (int j = 1; j < ny; j++)
         for (int i = 0; i < nx; i++)
         {
            if (fault(x0 + (i + 0.5) * hx, y0 + j * hy)) { m->AddBdrSegment(gv(i, j), gv(i + 1, j), 3); }
         }
   }
   m->FinalizeTopology();
   m->Finalize(false, true);
   m->SetAttributes();
   return m;
}

/// [-L, L]^2 in n x n cells, every wall attr 5.
static Mesh *DirichletBox(int n)
{
   return CrossGrid(n, n, -kL, -kL, 2 * kL / n, 2 * kL / n, [](real_t, real_t) { return 5; });
}

static BoundaryConfig ZeroDirichletBC()
{
   BoundaryConfig bc;
   bc.fault_attr = 3;           // absent from the mesh
   bc.dirichlet_attrs = {5};
   bc.default_dirichlet_func = [](const Vector &, real_t, Vector &u) { u.SetSize(2); u = 0.0; };
   return bc;
}

static std::unique_ptr<NitscheWaveOperator> MakeOp(ParMesh &pm, int p)
{
   return std::make_unique<NitscheWaveOperator>(pm, p, Material(), ZeroDirichletBC());
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

// Smooth data vanishing on the walls of [-L, L]^2.
static void SmoothU(const Vector &x, Vector &u)
{
   const real_t s = std::sin(M_PI * (x(0) + kL) / (2 * kL)) * std::sin(M_PI * (x(1) + kL) / (2 * kL));
   u.SetSize(2);
   u(0) = 1.0 * s;
   u(1) = 0.5 * s;
}
static void SmoothV(const Vector &x, Vector &v)
{
   const real_t s = std::sin(2 * M_PI * (x(0) + kL) / (2 * kL)) * std::sin(M_PI * (x(1) + kL) / (2 * kL));
   v.SetSize(2);
   v(0) = -0.2 * s;
   v(1) = 0.1 * s;
}
static void SmoothData(NitscheWaveOperator &op, Vector &U, Vector &V)
{
   ParGridFunction gu(&op.Elasticity().GetFESpace()), gv(&op.Elasticity().GetFESpace());
   VectorFunctionCoefficient cu(2, SmoothU), cv(2, SmoothV);
   gu.ProjectCoefficient(cu);
   gv.ProjectCoefficient(cv);
   U = gu;
   V = gv;
}

// ====================================================================== N2
// omega_max^2 = max eig M^{-1} K (Rayleigh quotient v.Kv / v.Mv), converged
// to a relative change below tol.
static real_t PowerIteration(NitscheWaveOperator &op, const MassMatrix &M, real_t tol,
                             int max_it, int &iters)
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

// Run nsteps from random U, V = 0; true if it blew up.
static bool BlowsUp(ParMesh &pm, int p, real_t dt, int nsteps)
{
   auto op = MakeOp(pm, p);
   op->SetTimeStep(dt);
   const int n = op->Height();
   Vector U(n), V(n);
   U.Randomize(777 + g_rank);
   U -= 0.5;
   V = 0.0;
   MPI_Comm comm = pm.GetComm();
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

static real_t CEig(ParMesh &pm, int p, int &iters)
{
   auto op = MakeOp(pm, p);
   const real_t dt_unit = ComputeMaxDt(op->Sample(), pm, p, 1.0);
   MassMatrix M;
   M.Build(*op);
   const real_t w2 = PowerIteration(*op, M, 1e-7, 20000, iters);
   return (2.0 / std::sqrt(w2)) / dt_unit;
}

static void TestCFL()
{
   if (g_rank == 0) { std::cout << "\n=== N2 the CFL constant on triangles ===\n"; }
   std::unique_ptr<Mesh> s4(DirichletBox(4)), s8(DirichletBox(8));
   ParMesh p4(MPI_COMM_WORLD, *s4), p8(MPI_COMM_WORLD, *s8);
   for (int p = 1; p <= 4; p++)
   {
      int it4 = 0, it8 = 0;
      const real_t c4 = CEig(p4, p, it4), c8 = CEig(p8, p, it8);
      auto op = MakeOp(p4, p);
      const real_t dt_unit = ComputeMaxDt(op->Sample(), p4, p, 1.0);
      real_t lo = 0.8 * c4, hi = 1.2 * c4;
      int guard = 0;
      while (BlowsUp(p4, p, lo * dt_unit, 1000) && guard++ < 10) { lo *= 0.7; }
      guard = 0;
      while (!BlowsUp(p4, p, hi * dt_unit, 1000) && guard++ < 10) { hi *= 1.3; }
      while ((hi - lo) / hi > 0.01)
      {
         const real_t mid = 0.5 * (lo + hi);
         if (BlowsUp(p4, p, mid * dt_unit, 1000)) { hi = mid; } else { lo = mid; }
      }
      const real_t rel = std::abs(lo - c4) / c4, rel8 = std::abs(c8 - c4) / c4;
      const real_t cmin = std::min(lo, c4), c2 = NitscheStableConstant(2, p);
      const real_t c3 = NitscheStableConstant(3, p);
      if (g_rank == 0)
      {
         std::cout << std::setprecision(6) << "  p = " << p << ": C_eig = " << c4 << " (4 x 4, "
                   << it4 << " iterations), " << c8 << " (8 x 8, " << it8 << "); bisection: "
                   << "largest stable C = " << lo << ", first unstable " << hi
                   << "\n    C(p): 2D " << c2 << ", 3D " << c3 << "; measured / 3D = " << cmin / c3
                   << "\n";
      }
      TEST_ASSERT(rel <= 0.0101, "N2 p = " << p << ": bisection agrees with 2/omega_max to 1%");
      TEST_ASSERT(rel8 <= 0.01, "N2 p = " << p << ": C_eig on the 8 x 8 grid within 1% of 4 x 4");
      // The bisection stops at a 1% bracket and its blow-up test near the
      // limit depends on how round-off seeds the top mode, so its result
      // differs by platform within that 1% (p = 3: 0.461572 in the cloud
      // container, 0.464472 on the Mac); the column holds the smallest
      // measurement, rounded down.  Check that it is safe here and no
      // further below this platform's measurement than the bracket.
      TEST_ASSERT(c2 <= cmin && c2 >= 0.99 * cmin, "N2 p = " << p << ": the 2D C(p) is at most "
                  "the measured min(C_eig, bisection) and within the bisection's 1% of it");
   }
}

// N2 on the committed unstructured meshes (100 m, every wall Dirichlet, the
// fault curve left as interior faces): C_eig by power iteration; the decks'
// cfl = 0.5 must stay below C_eig / C(p).  Opt-in (--only N2M; minutes).
static void TestCFLMeshes()
{
   if (g_rank == 0) { std::cout << "\n=== N2M the CFL constant on the committed 2D meshes ===\n"; }
   struct M { const char *file; int wall; };
   const M ms[2] = {{"tpv205/mesh/tpv205_2d_100m.msh", 105}, {"tpv105/mesh/tpv105_2d_100m.msh", 5}};
   for (const M &mm : ms)
   {
      Mesh smesh(mm.file, 1, 1);
      ParMesh pm(MPI_COMM_WORLD, smesh);
      for (int p = 1; p <= 3; p++)
      {
         BoundaryConfig bc;
         bc.fault_attr = 999;         // absent: the fault curve is interior
         bc.dirichlet_attrs = {mm.wall};
         bc.default_dirichlet_func = [](const Vector &, real_t, Vector &u) { u.SetSize(2); u = 0.0; };
         NitscheWaveOperator op(pm, p, Material(), bc);
         const real_t dt_unit = ComputeMaxDt(op.Sample(), pm, p, 1.0);
         MassMatrix M;
         M.Build(op);
         int it = 0;
         const real_t w2 = PowerIteration(op, M, 1e-6, 3000, it);
         const real_t ce = (2.0 / std::sqrt(w2)) / dt_unit, c2 = NitscheStableConstant(2, p);
         if (g_rank == 0)
         {
            std::cout << std::setprecision(6) << "  " << mm.file << " (" << smesh.GetNE()
                      << " triangles), p = " << p << ": C_eig = " << ce << " (" << it
                      << " iterations), C_eig / C(p) = " << ce / c2 << "\n";
         }
         TEST_ASSERT(ce > 0.5 * c2, "N2M " << mm.file << " p = " << p << ": cfl 0.5 x C(p) "
                     "below C_eig");
      }
   }
}

// ====================================================================== N3
// Serial-mesh reference for ComputeMaxDt with constant material: 4 A / P.
static real_t SerialDtRef(Mesh &smesh, int p, real_t cfl)
{
   const real_t p2 = real_t(p) * real_t(p);
   real_t dmin = std::numeric_limits<real_t>::infinity();
   for (int e = 0; e < smesh.GetNE(); ++e)
   {
      const real_t vol = smesh.GetElementVolume(e);
      Array<int> v;
      smesh.GetElementVertices(e, v);
      real_t P = 0.0;
      for (int i = 0; i < 3; i++)
      {
         const real_t *a = smesh.GetVertex(v[i]), *b = smesh.GetVertex(v[(i + 1) % 3]);
         P += std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]));
      }
      dmin = std::min(dmin, cfl * ((P > 0) ? 4.0 * vol / P : std::sqrt(vol)) / (p2 * Cp()));
   }
   return dmin;
}

static void TestStep(int p)
{
   if (g_rank == 0) { std::cout << "\n=== N3 the step, p = " << p << " ===\n"; }
   MPI_Comm comm = MPI_COMM_WORLD;
   std::unique_ptr<Mesh> s4(DirichletBox(4)), s8(DirichletBox(8));
   ParMesh p4(comm, *s4), p8(comm, *s8);
   const real_t C = CUse(p);

   // ComputeMaxDt against the serial 4A/P reference
   {
      auto op = MakeOp(p4, p);
      const real_t dt = ComputeMaxDt(op->Sample(), p4, p, 1.0), ref = SerialDtRef(*s4, p, 1.0);
      if (g_rank == 0)
      {
         std::cout << std::setprecision(17) << "  ComputeMaxDt(C = 1) = " << dt << ", serial 4A/P "
                   << "reference " << ref << "\n";
      }
      TEST_ASSERT(dt == ref, "N3 ComputeMaxDt == serial 4A/P reference (bitwise)");
   }

   // Verlet pin
   {
      auto op_a = MakeOp(p8, p), op_b = MakeOp(p8, p);
      const real_t dt = C * ComputeMaxDt(op_a->Sample(), p8, p, 1.0);
      op_a->SetTimeStep(dt);
      op_b->SetTimeStep(dt);
      const int n = op_a->Height(), nsteps = 1000;
      Vector U0(n), V0(n);
      SmoothData(*op_a, U0, V0);
      Vector Ua(U0), Va(V0);
      real_t ta = 0.0;
      RunCentralDifference(*op_a, Ua, Va, ta, dt, nsteps);
      const real_t fac0 = 0.5, fac2 = 0.5, fac3 = 0.0, fac4 = 0.5;
      Vector x(U0), v(V0), a(n);
      real_t t = 0.0;
      op_b->Mult(x, v, a);
      for (int s = 0; s < nsteps; ++s)
      {
         op_b->SetTime(t + dt);
         x.Add(dt, v);
         x.Add(fac0 * dt * dt, a);
         v.Add(fac2 * dt, a);
         op_b->SetTime(t + dt);
         op_b->ImplicitSolve(fac3 * dt * dt, fac4 * dt, x, v, a);
         x.Add(fac3 * dt * dt, a);
         v.Add(fac4 * dt, a);
         t += dt;
      }
      int d[2] = {0, 0};
      for (int i = 0; i < n; i++)
      {
         if (std::memcmp(&Ua(i), &x(i), sizeof(real_t)) != 0) { d[0]++; }
         if (std::memcmp(&Va(i), &v(i), sizeof(real_t)) != 0) { d[1]++; }
      }
      MPI_Allreduce(MPI_IN_PLACE, d, 2, MPI_INT, MPI_SUM, comm);
      const real_t nU = Norm(Ua, comm);
      if (g_rank == 0)
      {
         std::cout << std::setprecision(6) << "  Verlet pin (8 x 8): dt = " << dt << " s, ||U(T)|| = "
                   << nU << ", differing entries U " << d[0] << ", Udot " << d[1] << "\n";
      }
      TEST_ASSERT(d[0] == 0 && d[1] == 0 && ta == t,
                  "N3 hand Newmark loop == CentralDifferenceSolver bit for bit (1000 steps)");
   }

   // modified energy
   {
      auto op = MakeOp(p4, p);
      const real_t dt = C * ComputeMaxDt(op->Sample(), p4, p, 1.0);
      op->SetTimeStep(dt);
      MassMatrix M;
      M.Build(*op);
      const int n = op->Height(), nsteps = 10000;
      Vector U(n), V(n), KU(n), MV(n), A(n);
      SmoothData(*op, U, V);
      auto energies = [&](real_t &E, real_t &Eh, real_t &AtMA)
      {
         op->Elasticity().ApplyStiffness(U, KU);
         op->Mass().Apply(KU, A);
         M.Mult(V, MV);
         AtMA = Dot(KU, A, comm);
         E = 0.5 * Dot(V, MV, comm) + 0.5 * Dot(U, KU, comm);
         Eh = E - dt * dt / 8.0 * AtMA;
      };
      real_t E0, Eh0, a0;
      energies(E0, Eh0, a0);
      real_t max_dEh = 0, max_AtMA = a0, max_gap = 0;
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
         max_AtMA = std::max(max_AtMA, AtMA);
         max_gap = std::max(max_gap, std::abs(E - Eh0));
         sn += s; sE += E; snn += double(s) * s; snE += double(s) * E;
      }
      const double N = nsteps, slope = (N * snE - sn * sE) / (N * snn - sn * sn);
      const real_t bound = dt * dt / 8.0 * max_AtMA;
      if (g_rank == 0)
      {
         std::cout << std::setprecision(6) << "  energy (4 x 4): dt = " << dt << " s, max |E_h - "
                   << "E_h,0| / E_h,0 = " << max_dEh << ", max |E - E_h,0| / E_0 = " << max_gap / E0
                   << " (bound " << bound / E0 << "), drift " << slope / E0 << " E_0 per step, "
                   << "over the run " << std::abs(slope) * nsteps / bound << " of the bound\n";
      }
      TEST_ASSERT(max_dEh < 1e-10, "N3 modified energy E_h conserved to 1e-10 (10^4 steps)");
      TEST_ASSERT(max_gap <= bound * (1.0 + 1e-6) + 1e-10 * Eh0,
                  "N3 physical E within dt^2/8 max A^T M A of E_h,0");
      // The 3D test's slope bound (1e-10 E_0 per step) is not a drift test
      // at the 2D step: E oscillates by up to dt^2/8 A^T M A (1.2e-3 E_0 at
      // p = 1, the step 1.9 times P3's), so the fitted slope of a drift-free
      // E reaches 1e-10.  The drift over the run must stay below 1% of that
      // oscillation instead.
      TEST_ASSERT(std::abs(slope) * nsteps < 0.01 * bound,
                  "N3 E drift over the run < 1% of dt^2/8 max A^T M A");
   }

   // temporal order, from half the production step: at the production step
   // (half the stable one) the coarsest pair is not yet asymptotic on
   // triangles (U ratio 3.05 at p = 1); the finer pairs are 3.9 to 4.1
   {
      auto op0 = MakeOp(p4, p);
      const real_t dt0 = 0.5 * C * ComputeMaxDt(op0->Sample(), p4, p, 1.0);
      const int n0 = 64;
      Vector U0(op0->Height()), V0(op0->Height());
      SmoothData(*op0, U0, V0);
      auto run = [&](int refine, Vector &U, Vector &V)
      {
         auto op = MakeOp(p4, p);
         const real_t dt = dt0 / refine;
         op->SetTimeStep(dt);
         U = U0;
         V = V0;
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
         U -= Ur;
         V -= Vr;
         eU[k] = Norm(U, comm);
         eV[k] = Norm(V, comm);
      }
      const real_t rU1 = eU[0] / eU[1], rU2 = eU[1] / eU[2], rV1 = eV[0] / eV[1], rV2 = eV[1] / eV[2];
      if (g_rank == 0)
      {
         std::cout << std::setprecision(6) << "  temporal order (4 x 4, dt0 = " << dt0 << " s): U error ratios " << rU1
                   << ", " << rU2 << "; Udot " << rV1 << ", " << rV2 << "\n";
      }
      auto in = [](real_t r) { return r >= 3.5 && r <= 4.5; };
      TEST_ASSERT(in(rU1) && in(rU2) && in(rV1) && in(rV2), "N3 U and Udot error ratios in [3.5, 4.5]");
   }
}

// ====================================================================== N1
using ExactFn = std::function<void(const Vector &, real_t, Vector &)>;

static std::unique_ptr<NitscheWaveOperator> MakeWaveOp(ParMesh &pm, int p, const ExactFn &u_exact)
{
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.dirichlet_attrs = {5};
   bc.default_dirichlet_func = [u_exact](const Vector &x, real_t t, Vector &u) { u_exact(x, t, u); };
   return std::make_unique<NitscheWaveOperator>(pm, p, Material(), bc);
}

// From the exact (u, du/dt) at t = 0 to nsteps dt; ||u_h - u|| / ||u|| (L2).
static real_t RunExact(NitscheWaveOperator &op, const ExactFn &u_exact, const ExactFn &v_exact,
                       real_t dt, int nsteps)
{
   auto &fes = op.Elasticity().GetFESpace();
   ParGridFunction gu(&fes), gv(&fes);
   VectorFunctionCoefficient cu(2, u_exact), cv(2, v_exact);
   cu.SetTime(0.0);
   cv.SetTime(0.0);
   gu.ProjectCoefficient(cu);
   gv.ProjectCoefficient(cv);
   Vector U(gu), V(gv);
   op.SetTime(0.0);
   op.SetTimeStep(dt);
   real_t t = 0.0;
   CentralDifferenceSolver ode;
   ode.Init(op);
   for (int s = 0; s < nsteps; ++s)
   {
      real_t dts = dt;
      ode.Step(U, V, t, dts);
   }
   gu = U;
   cu.SetTime(t);
   ParGridFunction zero(&fes);
   zero = 0.0;
   return gu.ComputeL2Error(cu) / zero.ComputeL2Error(cu);
}

struct PlaneWave
{
   real_t k[2], d[2], omega;
   void U(const Vector &x, real_t t, Vector &u) const
   {
      const real_t ph = k[0] * x(0) + k[1] * x(1) - omega * t;
      u.SetSize(2);
      for (int i = 0; i < 2; i++) { u(i) = d[i] * std::sin(ph); }
   }
   void V(const Vector &x, real_t t, Vector &v) const
   {
      const real_t ph = k[0] * x(0) + k[1] * x(1) - omega * t;
      v.SetSize(2);
      for (int i = 0; i < 2; i++) { v(i) = -omega * d[i] * std::cos(ph); }
   }
};

static PlaneWave MakePlaneWave(bool pwave)
{
   real_t n[2] = {1.0, 0.7};
   const real_t nn = std::hypot(n[0], n[1]);
   n[0] /= nn;
   n[1] /= nn;
   const real_t kmag = 2.0 * M_PI / 120.0e3;     // wavelength 120 km
   PlaneWave w;
   w.k[0] = kmag * n[0];
   w.k[1] = kmag * n[1];
   if (pwave) { w.d[0] = n[0]; w.d[1] = n[1]; w.omega = Cp() * kmag; }
   else { w.d[0] = -n[1]; w.d[1] = n[0]; w.omega = Cs() * kmag; }
   return w;
}

static void TestPlaneWave(int p, bool pwave)
{
   if (g_rank == 0)
   {
      std::cout << "\n=== N1 plane " << (pwave ? "P" : "S") << " wave convergence, p = " << p << " ===\n";
   }
   const PlaneWave w = MakePlaneWave(pwave);
   ExactFn ue = [w](const Vector &x, real_t t, Vector &u) { w.U(x, t, u); };
   ExactFn ve = [w](const Vector &x, real_t t, Vector &v) { w.V(x, t, v); };
   const int nref[3] = {4, 8, 16};
   std::unique_ptr<ParMesh> pm[3];
   for (int m = 0; m < 3; m++)
   {
      std::unique_ptr<Mesh> sm(DirichletBox(nref[m]));
      pm[m] = std::make_unique<ParMesh>(MPI_COMM_WORLD, *sm);
   }
   const real_t T = 2.0;
   auto op_f = MakeWaveOp(*pm[2], p, ue);
   real_t dt = CUse(p) * ComputeMaxDt(op_f->Sample(), *pm[2], p, 1.0);
   int nsteps = (int) std::ceil(T / dt);
   dt = T / nsteps;
   real_t e_dt = RunExact(*op_f, ue, ve, dt, nsteps), change = 1.0;
   for (int h = 0; h < 4; ++h)
   {
      auto op_h = MakeWaveOp(*pm[2], p, ue);
      const real_t e_half = RunExact(*op_h, ue, ve, 0.5 * dt, 2 * nsteps);
      change = std::abs(e_dt - e_half) / e_half;
      if (g_rank == 0)
      {
         std::cout << std::setprecision(6) << "  16 x 16: dt = " << dt << " s err " << e_dt
                   << ", dt/2 err " << e_half << ", relative change " << change << "\n";
      }
      if (change < 0.01) { break; }
      dt *= 0.5;
      nsteps *= 2;
      e_dt = e_half;
   }
   TEST_ASSERT(change < 0.01, "N1 time error negligible on the finest mesh (< 1%)");
   real_t err[3];
   err[2] = e_dt;
   for (int m = 0; m < 2; m++)
   {
      auto op = MakeWaveOp(*pm[m], p, ue);
      err[m] = RunExact(*op, ue, ve, dt, nsteps);
   }
   const real_t r1 = std::log2(err[0] / err[1]), r2 = std::log2(err[1] / err[2]);
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6) << "  rel L2 error on 4 x 4 / 8 x 8 / 16 x 16: " << err[0]
                << ", " << err[1] << ", " << err[2] << "; observed order " << r1 << ", " << r2
                << " (expected ~" << p + 1 << ")\n";
   }
   TEST_ASSERT(r2 > p + 0.5, "N1 plane-wave L2 order on the finest pair > p + 0.5");
}

// ====================================================================== N4
// Pseudo-random density in [1000, 10000].
static real_t RhoField(const Vector &x)
{
   const real_t s = std::sin(1.3e-3 * x(0) + 2.1e-3 * x(1) + 1.7e-7 * x(0) * x(1));
   return 1000.0 + 9000.0 * s * s;
}

// Independent dashpot (test-local): quadrature 2p + 2 and C_e = N^T K N with
// N = phi (x) I_2, K = Zp n n^T + Zs (I - n n^T).
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
      const real_t Zp = rho * std::sqrt((lam + 2 * mu) / rho), Zs = rho * std::sqrt(mu / rho);
      const IntegrationRule &ir = IntRules.Get(Tr->GetGeometryType(), 2 * fe.GetOrder() + 2);
      DenseMatrix Ce(2 * nd), N(2, 2 * nd), KN(2, 2 * nd), K(2);
      Ce = 0.0;
      Vector shape(nd), nor(2);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         Tr->SetAllIntPoints(&ip);
         fe.CalcShape(Tr->GetElement1IntPoint(), shape);
         CalcOrtho(Tr->Jacobian(), nor);
         const real_t nl = nor.Norml2();
         nor /= nl;
         N = 0.0;
         for (int c = 0; c < 2; c++)
            for (int k = 0; k < nd; k++) { N(c, c * nd + k) = shape(k); }
         for (int c = 0; c < 2; c++)
            for (int d = 0; d < 2; d++) { K(c, d) = Zp * nor(c) * nor(d) + Zs * ((c == d) - nor(c) * nor(d)); }
         Mult(K, N, KN);
         AddMult_a_AtB(ip.weight * nl, N, KN, Ce);
      }
      fes.GetElementDofs(e, dofs);
      Vector x(2 * nd), y(2 * nd);
      for (int c = 0; c < 2; c++)
         for (int k = 0; k < nd; k++) { x(c * nd + k) = v(c * n + dofs[k]); }
      Ce.Mult(x, y);
      for (int c = 0; c < 2; c++)
         for (int k = 0; k < nd; k++) { Cv(c * n + dofs[k]) += y(c * nd + k); }
   }
}

static void TestDampedMass(int p)
{
   if (g_rank == 0) { std::cout << "\n=== N4 (M + dt/2 C) inverse, p = " << p << " ===\n"; }
   // top natural, the other walls absorbing
   std::unique_ptr<Mesh> m(CrossGrid(4, 4, -kL, -kL, 2 * kL / 4, 2 * kL / 4,
      [](real_t, real_t y) { return (std::abs(y - kL) < 1.0) ? 1 : 6; }));
   ParMesh pm(MPI_COMM_WORLD, *m);
   MPI_Comm comm = pm.GetComm();
   ConstantCoefficient lam_c(kLam), mu_c(kMu);
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
   const int g_damped = ISum(Md ? Md->NumDampedElements() : 0, comm);
   const int g_faces = ISum(op.Elasticity().GetAbsorbingBdrFaces().Size(), comm);
   TEST_ASSERT(Md != nullptr && g_damped == 12 && g_faces == 12,
               "N4 damped mass built: 12 absorbing edges, 12 damped triangles");
   if (!Md) { return; }
   ParBilinearForm mform(&op.Elasticity().GetFESpace());
   mform.AddDomainIntegrator(new VectorMassIntegrator(const_cast<NitscheMaterialSample &>(op.Sample()).Rho()));
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
      std::cout << std::setprecision(6) << "  ||(M + dt/2 C) Apply(r) - r|| / ||r|| = " << rel
                << ", ||C a - C_indep a|| / ||C a|| = " << relC << "\n";
   }
   TEST_ASSERT(rel < 1e-12, "N4 (M + dt/2 C) Apply(r) == r to round-off");
   TEST_ASSERT(relC < 1e-12, "N4 dashpot block == independent assembly");
}

// The reflection measurement (test_nitsche_p4.cpp T-12(a) in 2D; y up, the
// absorbing bottom y = -H, the medium above).  A plane pulse u = d g(t -
// s.(x - x0)), g(tau) = exp(-(tau/sigma)^2), whose centre hits the bottom at
// (0, -H) at t_hit.
struct Incident
{
   real_t th = 0;
   bool pwave = true;
   real_t c = 0, sigma = 0;
   real_t x0[2] = {0, 0}, xb[2] = {0, 0};
   real_t t_hit = 0;
   real_t s[2] = {0, 0}, d[2] = {0, 0};
   real_t g(real_t tau) const { const real_t a = tau / sigma; return std::exp(-a * a); }
   real_t dg(real_t tau) const { const real_t a = tau / sigma; return -2 * a / sigma * std::exp(-a * a); }
   void U(const Vector &x, real_t t, Vector &u, bool vel) const
   {
      const real_t sx = s[0] * (x(0) - x0[0]) + s[1] * (x(1) - x0[1]);
      const real_t a = vel ? dg(t - sx) : g(t - sx);
      u.SetSize(2);
      u(0) = d[0] * a;
      u(1) = d[1] * a;
   }
};

// Analytic first-order reflection of a plane P wave at angle th off the
// dashpot (as test_nitsche_p4.cpp, z -> y).
static void AnalyticRP(real_t th, real_t &Rpp, real_t &Rps)
{
   const real_t cp = Cp(), cs = Cs(), Zp = kRho * cp, Zs = kRho * cs;
   const real_t sphi = cs / cp * std::sin(th), cphi = std::sqrt(1 - sphi * sphi);
   struct W { real_t dx, dz, sx, sz; };
   const W wi{std::sin(th), -std::cos(th), std::sin(th) / cp, -std::cos(th) / cp};
   const W wp{std::sin(th), std::cos(th), std::sin(th) / cp, std::cos(th) / cp};
   const W ws{cphi, -sphi, sphi / cs, cphi / cs};
   auto e1 = [&](const W &w) { return -kMu * (w.dx * w.sz + w.dz * w.sx) - Zs * w.dx; };
   auto e2 = [&](const W &w) { return -(kLam * (w.dx * w.sx + w.dz * w.sz) + 2 * kMu * w.dz * w.sz) - Zp * w.dz; };
   const real_t a11 = e1(wp), a12 = e1(ws), a21 = e2(wp), a22 = e2(ws), b1 = -e1(wi), b2 = -e2(wi);
   const real_t det = a11 * a22 - a12 * a21;
   Rpp = (b1 * a22 - a12 * b2) / det;
   Rps = (a11 * b2 - a21 * b1) / det;
}

// [xl, xr] x [-H, 0] in cells of size h; the bottom attr 6 when absorbing,
// else 5; the other walls 5.  Element attr 2: the fit region (centroid
// -H < y < y_top, |x - xc| < xw), else 1.
static Mesh *ReflectionMesh(real_t xl, real_t xr, real_t H, real_t h, bool absorbing,
                            real_t y_top, real_t xc, real_t xw)
{
   const int nx = (int) std::lround((xr - xl) / h), ny = (int) std::lround(H / h);
   Mesh *m = CrossGrid(nx, ny, xl, -H, h, h, [=](real_t, real_t y)
                       { return (absorbing && std::abs(y + H) < 1e-6 * H) ? 6 : 5; });
   for (int e = 0; e < m->GetNE(); e++)
   {
      Array<int> ev;
      m->GetElementVertices(e, ev);
      real_t c[2] = {0, 0};
      for (int i = 0; i < ev.Size(); i++)
         for (int d = 0; d < 2; d++) { c[d] += m->GetVertex(ev[i])[d] / ev.Size(); }
      m->SetAttribute(e, (c[1] < y_top && std::abs(c[0] - xc) < xw) ? 2 : 1);
   }
   m->SetAttributes();
   return m;
}

// U_B on mesh B mapped onto mesh A's DOFs (elements matched by centroid,
// DOFs by node position); serial.
static void MapBtoA(NitscheWaveOperator &opA, NitscheWaveOperator &opB, const Vector &UB, Vector &UBA)
{
   ParMesh &mA = opA.GetMesh(), &mB = opB.GetMesh();
   auto &fA = opA.Elasticity().GetFESpace();
   auto &fB = opB.Elasticity().GetFESpace();
   const int nA = fA.GetNDofs(), nB = fB.GetNDofs();
   auto centroid = [](Mesh &m, int e) {
      Array<int> ev;
      m.GetElementVertices(e, ev);
      real_t c[2] = {0, 0};
      for (int i = 0; i < ev.Size(); i++)
         for (int d = 0; d < 2; d++) { c[d] += m.GetVertex(ev[i])[d]; }
      return std::array<long long, 2>{std::llround(c[0] / ev.Size()), std::llround(c[1] / ev.Size())};
   };
   std::map<std::array<long long, 2>, int> bmap;
   for (int e = 0; e < mB.GetNE(); e++) { bmap[centroid(mB, e)] = e; }
   UBA.SetSize(2 * nA);
   UBA = 0.0;
   Array<int> dA, dB;
   for (int e = 0; e < mA.GetNE(); e++)
   {
      auto it = bmap.find(centroid(mA, e));
      MFEM_VERIFY(it != bmap.end(), "MapBtoA: element " << e << " not in mesh B");
      const int eb = it->second;
      const IntegrationRule &nodes = fA.GetFE(e)->GetNodes();
      const int nn = nodes.GetNPoints();
      std::vector<Vector> pa(nn), pb(nn);
      for (int k = 0; k < nn; k++)
      {
         mA.GetElementTransformation(e)->Transform(nodes.IntPoint(k), pa[k]);
         mB.GetElementTransformation(eb)->Transform(nodes.IntPoint(k), pb[k]);
      }
      fA.GetElementDofs(e, dA);
      fB.GetElementDofs(eb, dB);
      for (int k = 0; k < nn; k++)
      {
         int match = -1;
         for (int l = 0; l < nn && match < 0; l++)
         {
            Vector dd(pa[k]);
            dd -= pb[l];
            if (dd.Norml2() < 1e-3) { match = l; }
         }
         MFEM_VERIFY(match >= 0, "MapBtoA: node " << k << " of element " << e << " unmatched");
         for (int cc = 0; cc < 2; cc++) { UBA(cc * nA + dA[k]) = UB(cc * nB + dB[match]); }
      }
   }
}

static real_t RegionDot(ParFiniteElementSpace &fes, const Vector &a, const Vector &b)
{
   ParMesh *pm = fes.GetParMesh();
   Array<int> marker(pm->attributes.Max());
   marker = 0;
   marker[1] = 1;      // attribute 2: the fit region
   ConstantCoefficient one(1.0);
   ParBilinearForm m(&fes);
   m.AddDomainIntegrator(new VectorMassIntegrator(one), marker);
   m.Assemble();
   m.Finalize();
   Vector Mb(b.Size());
   m.SpMat().Mult(b, Mb);
   return Dot(a, Mb, pm->GetComm());
}

struct Geometry2D { real_t xl, xr, H, y_top, xw, t_after; };
struct ReflResult { real_t Rpp = 0, Rps = 0, dnorm = 0; };

static ReflResult MeasureReflection(const Incident &I, int p, real_t h, const Geometry2D &G)
{
   std::unique_ptr<Mesh> mA(ReflectionMesh(G.xl, G.xr, G.H, h, true, G.y_top, 0.0, G.xw));
   std::unique_ptr<Mesh> mB(ReflectionMesh(G.xl, G.xr, 3 * G.H, h, false, G.y_top, 0.0, G.xw));
   ParMesh pA(MPI_COMM_WORLD, *mA), pB(MPI_COMM_WORLD, *mB);
   ExactFn uex = [I](const Vector &x, real_t t, Vector &u) { I.U(x, t, u, false); };
   ExactFn vex = [I](const Vector &x, real_t t, Vector &v) { I.U(x, t, v, true); };
   auto make = [&](ParMesh &pm, bool absorbing) {
      BoundaryConfig bc;
      bc.fault_attr = 3;
      bc.dirichlet_attrs = {5};
      if (absorbing) { bc.absorbing_attrs = {6}; }
      bc.default_dirichlet_func = uex;
      return std::make_unique<NitscheWaveOperator>(pm, p, Material(), bc);
   };
   auto opA = make(pA, true), opB = make(pB, false);
   const real_t dtA = ComputeMaxDt(opA->Sample(), pA, p, CUse(p));
   const real_t dtB = ComputeMaxDt(opB->Sample(), pB, p, CUse(p));
   const real_t Tend = I.t_hit + G.t_after;
   const int n = (int) std::ceil(Tend / std::min(dtA, dtB));
   const real_t dt = Tend / n;
   auto run = [&](NitscheWaveOperator &op, Vector &U) {
      auto &fes = op.Elasticity().GetFESpace();
      ParGridFunction gu(&fes), gv(&fes);
      VectorFunctionCoefficient cu(2, uex), cv(2, vex);
      cu.SetTime(0.0);
      cv.SetTime(0.0);
      gu.ProjectCoefficient(cu);
      gv.ProjectCoefficient(cv);
      U = gu;
      Vector V(gv);
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
   const real_t cp = Cp(), cs = Cs();
   real_t sP[2], dP[2], sS[2], dS[2];
   if (I.pwave)
   {
      const real_t sphi = cs / cp * std::sin(I.th), cphi = std::sqrt(1 - sphi * sphi);
      sP[0] = std::sin(I.th) / cp; sP[1] = std::cos(I.th) / cp;
      dP[0] = std::sin(I.th);      dP[1] = std::cos(I.th);
      sS[0] = sphi / cs;           sS[1] = cphi / cs;
      dS[0] = cphi;                dS[1] = -sphi;
   }
   else   // normal-incidence S (x-polarised): R_same is the reflected S
   {
      sP[0] = 0; sP[1] = 1.0 / cs; dP[0] = 1; dP[1] = 0;
      sS[0] = 0; sS[1] = 1.0 / cp; dS[0] = 0; dS[1] = 1;
   }
   auto &fes = opA->Elasticity().GetFESpace();
   auto shape = [&](const real_t sr[2], const real_t dr[2], Vector &out) {
      ParGridFunction g(&fes);
      VectorFunctionCoefficient cf(2, [&I, sr, dr](const Vector &x, real_t t, Vector &u) {
         const real_t sx = sr[0] * (x(0) - I.xb[0]) + sr[1] * (x(1) - I.xb[1]);
         const real_t a = I.g(t - I.t_hit - sx);
         u.SetSize(2);
         u(0) = dr[0] * a;
         u(1) = dr[1] * a;
      });
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

// A pulse whose centre passes through (x0) at t = 0 and hits (0, -H) at t_hit.
static Incident MakeIncident(real_t th_deg, bool pwave, real_t H, real_t y0)
{
   Incident I;
   I.th = th_deg * M_PI / 180.0;
   I.pwave = pwave;
   I.c = pwave ? Cp() : Cs();
   I.sigma = pwave ? 1.2 : 2.0;
   const real_t dy = H + y0;                  // from y0 to the bottom
   I.x0[0] = -dy * std::tan(I.th);
   I.x0[1] = y0;
   I.xb[0] = 0.0;
   I.xb[1] = -H;
   I.t_hit = dy / std::cos(I.th) / I.c;
   if (pwave)
   {
      I.s[0] = std::sin(I.th) / I.c; I.s[1] = -std::cos(I.th) / I.c;
      I.d[0] = std::sin(I.th);       I.d[1] = -std::cos(I.th);
   }
   else
   {
      I.s[0] = 0; I.s[1] = -1.0 / I.c;
      I.d[0] = 1; I.d[1] = 0;
   }
   return I;
}

static void TestReflection()
{
   if (Mpi::WorldSize() > 1) { return; }   // MapBtoA is serial
   std::cout << "\n=== N4 reflection off the dashpot (p = 2) ===\n"
                "  A (absorbing bottom) minus deep reference B (3 H, exact incident data on\n"
                "  every wall), fitted to the analytic reflected P and S pulses\n";
   // normal incidence: a 10 km column, H = 30 km, pulse from y = -8 km,
   // fit over -30 < y < -10 km at 10 km / c after the hit
   for (int w = 0; w < 2; w++)
   {
      const bool pw = (w == 0);
      const Incident I = MakeIncident(0.0, pw, 30.0e3, -8.0e3);
      const Geometry2D G{-5.0e3, 5.0e3, 30.0e3, -10.0e3, 1e9, 10.0e3 / I.c};
      const ReflResult r5 = MeasureReflection(I, 2, 5.0e3, G), r25 = MeasureReflection(I, 2, 2.5e3, G);
      std::cout << std::setprecision(5) << "  normal " << (pw ? "P" : "S") << ": |D|/|pulse| = "
                << r5.dnorm << " (5 km cells), " << r25.dnorm << " (2.5 km); fitted R_same = "
                << r5.Rpp << ", " << r25.Rpp << " (analytic 0)\n";
      TEST_ASSERT(std::abs(r25.Rpp) < 0.05 && r25.dnorm < r5.dnorm,
                  "N4 normal-incidence " << (pw ? "P" : "S") << " reflection small and falling "
                  "with resolution");
   }
   // P at 45 degrees: H = 50 km, the pulse centre on the top wall at t = 0
   // (it reaches the bottom first at x = -2 H = -100 km, outside the box),
   // the box x in [-76, 30] km: the side walls' incident-only data reaches
   // the fit region (|x| < 10 km, -50 < y < -35 km) about 2 s after the hit,
   // the fit is 1 s after it.
   {
      const Incident I = MakeIncident(45.0, true, 50.0e3, 0.0);
      const Geometry2D G{-76.0e3, 30.0e3, 50.0e3, -35.0e3, 10.0e3, 1.0};
      real_t Rpp, Rps;
      AnalyticRP(I.th, Rpp, Rps);
      const ReflResult r2 = MeasureReflection(I, 2, 2.0e3, G), r1 = MeasureReflection(I, 2, 1.0e3, G);
      std::cout << std::setprecision(5) << "  P at 45 deg: measured R_PP = " << r2.Rpp << ", "
                << r1.Rpp << " (2 km, 1 km cells), R_PS = " << r2.Rps << ", " << r1.Rps
                << "; analytic first-order R_PP = " << Rpp << ", R_PS = " << Rps << "\n";
      // first run (10-25): within 5e-5 (R_PP) and 1.9e-4 (R_PS) at both
      // resolutions; the bound is ten times that
      auto near = [&](const ReflResult &r) {
         return std::abs(r.Rpp - Rpp) < 2e-3 && std::abs(r.Rps - Rps) < 2e-3;
      };
      TEST_ASSERT(near(r2) && near(r1), "N4 P at 45 deg: R_PP and R_PS within 0.002 of the "
                  "analytic first-order values on 2 km and 1 km cells");
   }
}

// ====================================================================== N5
static const real_t kMuS = 0.677, kMuD = 0.525, kDc = 0.40;
static const real_t kTau0 = 70.0e6, kSigma0 = 120.0e6;
static const real_t kNucR = 9.0e3, kNucDtau = 15.0e6;
static const real_t kLf = 15.0e3;     // the fault |x| < 15 km
static const int kNf = 24;            // 2.5 km cells

static real_t NucF(real_t r)
{
   if (r >= kNucR) { return 0.0; }
   const real_t r2 = r * r, R2 = kNucR * kNucR;
   return std::exp(r2 / (r2 - R2));
}

/// rotate = false: the fault on y = 0, up = [0, 0, 1] (strike slot).
/// rotate = true: the same problem turned by +90 degrees, (x, y) -> (-y, x):
/// the fault on x = 0, up = [0, 1] (dip slot).
static Mesh *FaultBox(bool rotate)
{
   const real_t h = 2 * kL / kNf;
   return CrossGrid(kNf, kNf, -kL, -kL, h, h, [](real_t, real_t) { return 5; },
                    [](real_t x, real_t y) { return y == 0.0 && std::abs(x) < kLf; }, rotate);
}

// sigma in the rotated frame: R s R^T with R = [[0, -1], [1, 0]].
static void RotateStress(bool rotate, real_t s[9])
{
   if (!rotate) { return; }
   const real_t xx = s[0], xy = s[1], yy = s[4];
   s[0] = yy;
   s[4] = xx;
   s[1] = s[3] = -xy;
}

static NitscheFaultSetup FaultSetup(const IFaultClosure &cl, bool rotate)
{
   NitscheFaultSetup fs;
   fs.closure = &cl;
   const real_t ref[3] = {rotate ? 1.0 : 0.0, rotate ? 0.0 : -1.0, 0.0};   // R (0, -1)
   const real_t up[3] = {0.0, rotate ? 1.0 : 0.0, rotate ? 0.0 : 1.0};
   for (int d = 0; d < 3; d++)
   {
      fs.orientation.ref_normal[d] = ref[d];
      fs.orientation.up[d] = up[d];
      fs.up[d] = up[d];
   }
   fs.orientation.fault_attrs = {3};
   fs.prestress = [rotate](const real_t[3], real_t s[9]) {
      for (int i = 0; i < 9; i++) { s[i] = 0.0; }
      s[1] = s[3] = kTau0;
      s[4] = -kSigma0;
      RotateStress(rotate, s);
   };
   fs.nucleation_spatial = [rotate](const real_t x[3], real_t s[9]) {
      for (int i = 0; i < 9; i++) { s[i] = 0.0; }
      const real_t along = rotate ? x[1] : x[0];
      s[1] = s[3] = kNucDtau * NucF(std::abs(along));
      RotateStress(rotate, s);
   };
   fs.nucleation_temporal = [](real_t t) { return NucleationTemporal(t); };
   fs.parameters = [](const real_t[3], FaultPointState &st) {
      st.mu_s = kMuS;
      st.mu_d = kMuD;
      st.Dc = kDc;
      st.cohesion = 0.0;
   };
   fs.sigma_min = 0.0;
   return fs;
}

struct FaultRunConfig
{
   int p = 1;
   bool rotate = false;
   real_t tfinal = 3.0;
   bool plain = false;        ///< fast kernels, dense blocks and fused sweep off
   real_t kv_eta = 0.0;       ///< > 0: Kelvin-Voigt
   real_t kv_width = -1.0;    ///< its layer (< 0: every element)
   real_t pc_length = 0.0;    ///< > 0: Prakash-Clifton (v* = 1 m/s)
   bool self = false;         ///< on MPI_COMM_SELF (rank 0 only)
};

/// Per owned DOF: (x, y) of the unrotated problem, the face midpoint, and
/// the values d (2, unrotated frame), V (2), slip path, sigma_eff.
struct FaultRow { std::array<long long, 4> key; std::array<real_t, 6> val; };

struct FaultRunResult
{
   bool finite = true;
   long steps = 0;
   real_t dt = 0;
   int fault_dofs = 0;
   int shared_faces = 0;         ///< fault faces on a rank boundary
   int prop1_bad = 0;            ///< nonzero out-of-plane entries seen
   real_t sig_dev = 0;           ///< max |sigma_eff - sigma_0| over DOFs and steps
   real_t first_t = -1, first_dist = -1;
   real_t dmax = 0, Vmax = 0;
   long long kv_layer = 0;
   std::vector<FaultRow> rows;   ///< gathered on rank 0 of the run's comm, sorted
};

static FaultRunResult RunFault(const FaultRunConfig &rc)
{
   FaultRunResult R;
   MPI_Comm comm = rc.self ? MPI_COMM_SELF : MPI_COMM_WORLD;
   std::unique_ptr<Mesh> m(FaultBox(rc.rotate));
   int crank = 0, csize = 1;
   MPI_Comm_rank(comm, &crank);
   MPI_Comm_size(comm, &csize);
   // On 2 and 4 ranks a partition that cuts along the fault (the two sides
   // on different ranks; on 4 also across it at x = 0), so that every fault
   // face is a rank-boundary face; METIS otherwise.
   std::vector<int> part;
   if (csize == 2 || csize == 4)
   {
      part.resize(m->GetNE());
      for (int e = 0; e < m->GetNE(); e++)
      {
         Vector c;
         m->GetElementCenter(e, c);
         const real_t along = rc.rotate ? c(1) : c(0), across = rc.rotate ? -c(0) : c(1);
         part[e] = (across < 0.0 ? 1 : 0) + ((csize == 4 && along < 0.0) ? 2 : 0);
      }
   }
   ParMesh pm(comm, *m, part.empty() ? nullptr : part.data());
   m.reset();
   SlipWeakeningClosure cl;
   NitscheFaultSetup fs = FaultSetup(cl, rc.rotate);
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.absorbing_attrs = {5};
   NitscheWaveOperator op(pm, rc.p, Material(), bc, {}, &fs);
   if (rc.plain)
   {
      op.SetFastFaultKernels(false);
      op.SetDenseBlocks(false);
      op.SetFusedStiffness(false);
   }
   const real_t dt = ComputeMaxDt(op.Sample(), pm, rc.p, CUse(rc.p));
   op.SetTimeStep(dt);
   if (rc.kv_eta > 0.0)
   {
      op.SetKelvinVoigt(rc.kv_eta, NitscheWaveOperator::KVScale::Node, rc.kv_width);
      R.kv_layer = op.KelvinVoigtLayerElements();
   }
   const int no = op.NumOwnedFaultDOFs();
   {
      Vector d(3 * no), V(3 * no), psi(no);
      d = 0.0;
      V = 0.0;
      for (int j = 0; j < no; j++)
      {
         FaultPointState st = op.FaultParams(j);
         op.Closure().SeedEquilibrium(kTau0, 0.0, kSigma0, st);
         psi(j) = st.psi;
      }
      op.SetFaultState(d, V, psi);
   }
   if (rc.pc_length > 0.0)
   {
      Vector S, sig0(no);
      op.BackgroundTraction(0.0, S);
      for (int j = 0; j < no; j++)
      {
         const real_t *f = &op.FaultFrames()(9 * j);
         real_t X0, X1;
         fault_frame::ToLocal(&S(3 * j), 0.0, f, f + 3, f + 6, X0, X1, sig0(j));
      }
      op.SetPrakashClifton(1.0, rc.pc_length);
      op.SetPrakashCliftonState(sig0);
   }
   R.dt = dt;
   R.fault_dofs = ISum(no, comm);
   R.shared_faces = op.NumSharedFaultFacesGlobal();
   const int slot_out = rc.rotate ? 1 : 0;   // the out-of-plane slot of X: strike (rotated) or dip
   // positions in the unrotated frame: (x, y) = (y', -x') for the rotated run
   const Vector &Xp = op.FaultPositions();
   auto ux = [&](int j) { return rc.rotate ? Xp(3 * j + 1) : Xp(3 * j); };
   auto uy = [&](int j) { return rc.rotate ? -Xp(3 * j) : Xp(3 * j + 1); };
   int bad = 0;
   auto prop1 = [&](const Vector &v3) {
      for (int j = 0; j < no; j++) { if (v3(3 * j + 2) != 0.0) { bad++; } }
   };
   prop1(op.PrestressTraction());
   prop1(op.NucleationTraction());
   for (int j = 0; j < no; j++) { if (Xp(3 * j + 2) != 0.0) { bad++; } }
   std::vector<long> first_step(no, -1);
   const int n = op.Height();
   Vector U(n), Ud(n);
   U = 0.0;
   Ud = 0.0;
   real_t t = 0.0;
   op.SetTime(t);
   CentralDifferenceSolver ode;
   ode.Init(op);
   const long nsteps = std::lround(std::ceil(rc.tfinal / dt - 1e-9));
   real_t sig_dev = 0.0;
   for (long s = 0; s < nsteps; s++)
   {
      real_t h_s = dt;
      ode.Step(U, Ud, t, h_s);
      const NitscheStepRecord &rec = op.LastStep();
      if (!std::isfinite(rec.max_V) || !std::isfinite(rec.min_sigma_eff)) { R.finite = false; break; }
      prop1(op.Slip());
      prop1(op.SlipRate());
      prop1(op.LastTrialTraction());
      const Vector &X = op.LastShear(), &sg = op.LastSigmaEff();
      for (int j = 0; j < no; j++)
      {
         if (X(2 * j + slot_out) != 0.0) { bad++; }
         sig_dev = std::max(sig_dev, std::abs(sg(j) - kSigma0));
         const real_t *V = &op.SlipRate()(3 * j);
         if (first_step[j] < 0 && (V[0] != 0.0 || V[1] != 0.0)) { first_step[j] = s + 1; }
      }
   }
   R.steps = op.LastStep().step;
   R.finite = (ISum(R.finite ? 0 : 1, comm) == 0) && std::isfinite(RMax(U.Normlinf(), comm));
   R.prop1_bad = ISum(bad, comm);
   R.sig_dev = RMax(sig_dev, comm);
   {
      long smin = std::numeric_limits<long>::max();
      for (int j = 0; j < no; j++) { if (first_step[j] > 0) { smin = std::min(smin, first_step[j]); } }
      MPI_Allreduce(MPI_IN_PLACE, &smin, 1, MPI_LONG, MPI_MIN, comm);
      real_t far = -1.0;
      for (int j = 0; j < no; j++) { if (first_step[j] == smin) { far = std::max(far, std::abs(ux(j))); } }
      far = RMax(far, comm);
      if (smin != std::numeric_limits<long>::max()) { R.first_t = smin * dt; R.first_dist = far; }
   }
   // rows: the unrotated frame's values, d and V as (x, y) components
   ElasticityDomainOperator<ParMesh> &E = op.Elasticity();
   const int nd = E.GetNumFaultDOFs(), nbf = E.GetNbfPerFace();
   Vector c_loc(3 * nd), c_own;
   for (int k = 0; k < nd; k++)
   {
      long long key[3];
      real_t xv[3][3];
      op.Orientation().FaceVertexData(k / nbf, key, xv);
      for (int d = 0; d < 3; d++) { c_loc(3 * k + d) = 0.5 * (xv[0][d] + xv[1][d]); }
   }
   E.RestrictToOwnedFault(c_loc, c_own, 3);
   std::vector<double> buf(10 * no);
   real_t dmax = 0.0, Vmax = 0.0;
   for (int j = 0; j < no; j++)
   {
      double *b = &buf[10 * j];
      const real_t *d = &op.Slip()(3 * j), *V = &op.SlipRate()(3 * j);
      b[0] = ux(j);
      b[1] = uy(j);
      b[2] = rc.rotate ? c_own(3 * j + 1) : c_own(3 * j);
      b[3] = rc.rotate ? -c_own(3 * j) : c_own(3 * j + 1);
      b[4] = rc.rotate ? d[1] : d[0];
      b[5] = rc.rotate ? -d[0] : d[1];
      b[6] = rc.rotate ? V[1] : V[0];
      b[7] = rc.rotate ? -V[0] : V[1];
      b[8] = op.State()(j);
      b[9] = op.LastSigmaEff()(j);
      dmax = std::max(dmax, std::hypot(d[0], d[1]));
      Vmax = std::max(Vmax, std::hypot(V[0], V[1]));
   }
   R.dmax = RMax(dmax, comm);
   R.Vmax = RMax(Vmax, comm);
   int nloc = 10 * no;
   std::vector<int> counts(csize), displs(csize, 0);
   MPI_Gather(&nloc, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm);
   std::vector<double> all;
   if (crank == 0)
   {
      for (int r = 1; r < csize; r++) { displs[r] = displs[r - 1] + counts[r - 1]; }
      all.resize(displs[csize - 1] + counts[csize - 1]);
   }
   MPI_Gatherv(buf.data(), nloc, MPI_DOUBLE, all.data(), counts.data(), displs.data(), MPI_DOUBLE, 0, comm);
   if (crank == 0)
   {
      const size_t nr = all.size() / 10;
      R.rows.resize(nr);
      for (size_t k = 0; k < nr; k++)
      {
         const double *b = &all[10 * k];
         for (int c = 0; c < 4; c++) { R.rows[k].key[c] = std::llround(b[c] * 1e3); }
         R.rows[k].val = {b[4], b[5], b[6], b[7], b[8], b[9]};
      }
      std::sort(R.rows.begin(), R.rows.end(),
                [](const FaultRow &a, const FaultRow &b) { return a.key < b.key; });
   }
   return R;
}

/// max over matched rows of |a - b| per value column, relative to the
/// largest |a| of the column group (slip, slip rate, slip path, sigma_eff);
/// -1 when the keys differ.
static real_t RowsDiff(const FaultRunResult &A, const FaultRunResult &B)
{
   if (A.rows.size() != B.rows.size()) { return -1.0; }
   real_t scale[4] = {0, 0, 0, 0}, diff[4] = {0, 0, 0, 0};
   const int grp[6] = {0, 0, 1, 1, 2, 3};
   for (size_t k = 0; k < A.rows.size(); k++)
   {
      if (A.rows[k].key != B.rows[k].key) { return -1.0; }
      for (int c = 0; c < 6; c++)
      {
         scale[grp[c]] = std::max(scale[grp[c]], std::abs(A.rows[k].val[c]));
         diff[grp[c]] = std::max(diff[grp[c]], std::abs(A.rows[k].val[c] - B.rows[k].val[c]));
      }
   }
   real_t r = 0.0;
   for (int g = 0; g < 4; g++) { if (scale[g] > 0) { r = std::max(r, diff[g] / scale[g]); } }
   return r;
}

/// The half turn (x, y) -> (-x, -y): it maps the problem onto itself and
/// the + side of the fault onto the - side, so d(x) = d(-x) and V, psi
/// likewise; d_x, V_x and the slip path at x against -x, relative to the
/// largest values.
static real_t MirrorX(const FaultRunResult &A)
{
   std::map<std::array<long long, 4>, const FaultRow *> by;
   for (const FaultRow &r : A.rows) { by[r.key] = &r; }
   real_t scale[3] = {0, 0, 0}, diff[3] = {0, 0, 0};
   for (const FaultRow &r : A.rows)
   {
      auto it = by.find({-r.key[0], r.key[1], -r.key[2], r.key[3]});
      if (it == by.end()) { return -1.0; }
      const FaultRow &q = *it->second;
      const real_t a[3] = {r.val[0], r.val[2], r.val[4]}, b[3] = {q.val[0], q.val[2], q.val[4]};
      for (int g = 0; g < 3; g++)
      {
         scale[g] = std::max(scale[g], std::abs(a[g]));
         diff[g] = std::max(diff[g], std::abs(a[g] - b[g]));
      }
   }
   real_t rr = 0.0;
   for (int g = 0; g < 3; g++) { if (scale[g] > 0) { rr = std::max(rr, diff[g] / scale[g]); } }
   return rr;
}

static void TestFault(int p)
{
   if (g_rank == 0) { std::cout << "\n=== N5 the fault, p = " << p << " ===\n"; }
   const int np = Mpi::WorldSize();
   FaultRunConfig rc;
   rc.p = p;
   const FaultRunResult A = RunFault(rc);
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6) << "  strike slot: " << A.fault_dofs << " fault DOFs, dt = "
                << A.dt << " s, " << A.steps << " steps; max |d| = " << A.dmax << " m, max |V| = "
                << A.Vmax << " m/s; first slip at t = " << A.first_t << " s, |x| <= " << A.first_dist
                << " m; out-of-plane nonzeros " << A.prop1_bad << "; max |sigma_eff - sigma_0| = "
                << A.sig_dev << " Pa\n";
   }
   TEST_ASSERT(A.finite && A.steps > 0, "N5 p = " << p << ": the 2D rupture runs finite");
   TEST_ASSERT(A.prop1_bad == 0, "N5 p = " << p << ": Proposition 1, the out-of-plane slot is "
               "exactly 0 (d, V, trial traction, prestress, nucleation, X)");
   TEST_ASSERT(A.sig_dev < 1e-8 * kSigma0, "N5 p = " << p << ": mirror symmetry about the fault, "
               "|sigma_eff - 120 MPa| < 1e-8 x 120 MPa at every DOF and step");
   TEST_ASSERT(A.first_t > 0 && A.first_t <= TPV102Params::nuc_T && A.first_dist < kNucR,
               "N5 p = " << p << ": the first DOFs to slip lie in the nucleation patch, within its rise time");
   TEST_ASSERT(A.dmax > 1.0, "N5 p = " << p << ": slip above 1 m by 3 s");
   real_t mx = 0.0;
   if (g_rank == 0) { mx = MirrorX(A); }
   MPI_Bcast(&mx, 1, MPITypeMap<real_t>::mpi_type, 0, MPI_COMM_WORLD);
   if (g_rank == 0) { std::cout << "  half turn, x against -x: " << mx << "\n"; }
   TEST_ASSERT(mx >= 0 && mx < 1e-9, "N5 p = " << p << ": half-turn symmetry, slip, slip rate "
               "and slip path at x and -x to 1e-9");

   FaultRunConfig rr = rc;
   rr.rotate = true;
   const FaultRunResult B = RunFault(rr);
   FaultRunConfig rp = rc;
   rp.plain = true;
   const FaultRunResult C = RunFault(rp);
   real_t dB = 0.0, dC = 0.0;
   if (g_rank == 0) { dB = RowsDiff(A, B); dC = RowsDiff(A, C); }
   real_t dd[2] = {dB, dC};
   MPI_Bcast(dd, 2, MPITypeMap<real_t>::mpi_type, 0, MPI_COMM_WORLD);
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6) << "  dip slot (rotated 90 deg, up = [0, 1]): out-of-plane "
                << "nonzeros " << B.prop1_bad << ", max |sigma_eff - sigma_0| = " << B.sig_dev
                << " Pa; against the strike-slot run " << dd[0] << "\n  CSR products and "
                << "ComputeTractionDiagnostics against the fast path: " << dd[1] << "\n";
   }
   TEST_ASSERT(B.finite && B.prop1_bad == 0 && B.sig_dev < 1e-8 * kSigma0,
               "N5 p = " << p << ": dip slot: finite, Proposition 1, mirror symmetry");
   TEST_ASSERT(dd[0] >= 0 && dd[0] < 1e-9, "N5 p = " << p << ": the rotated problem (dip slot) "
               "equals the strike-slot run to 1e-9");
   TEST_ASSERT(dd[1] >= 0 && dd[1] < 1e-11, "N5 p = " << p << ": fast fault kernels, dense blocks "
               "and fused sweep equal the CSR products and ComputeTractionDiagnostics to 1e-11");

   if (np > 1)
   {
      FaultRunResult S;
      FaultRunConfig rs = rc;
      rs.self = true;
      if (g_rank == 0) { S = RunFault(rs); }
      real_t dS = 0.0;
      if (g_rank == 0) { dS = RowsDiff(A, S); }
      MPI_Bcast(&dS, 1, MPITypeMap<real_t>::mpi_type, 0, MPI_COMM_WORLD);
      if (g_rank == 0)
      {
         std::cout << "  " << np << " ranks (" << A.shared_faces << " of " << A.fault_dofs / (p + 1)
                   << " fault faces on a rank boundary) against 1 rank: " << dS << "\n";
      }
      TEST_ASSERT(dS >= 0 && dS < 1e-11, "N5 p = " << p << ": " << np << " ranks equal 1 rank to 1e-11");
      TEST_ASSERT(np > 4 || A.shared_faces == A.fault_dofs / (p + 1), "N5 p = " << p
                  << ": on 2 and 4 ranks every fault face is a rank-boundary face");
   }

   if (p == 1)
   {
      FaultRunConfig rk = rc;
      rk.kv_eta = 0.1;
      const FaultRunResult K = RunFault(rk);
      rk.kv_width = 2.0 * (2 * kL / kNf);
      const FaultRunResult Kl = RunFault(rk);
      FaultRunConfig rq = rc;
      rq.pc_length = 0.01;
      const FaultRunResult Q = RunFault(rq);
      rq.pc_length = 1e-200;
      const FaultRunResult Q0 = RunFault(rq);
      real_t dq = 0.0;
      if (g_rank == 0) { dq = RowsDiff(A, Q0); }
      MPI_Bcast(&dq, 1, MPITypeMap<real_t>::mpi_type, 0, MPI_COMM_WORLD);
      if (g_rank == 0)
      {
         std::cout << std::setprecision(6) << "  Kelvin-Voigt eta 0.1: whole domain max |d| = " << K.dmax
                   << " m; layer of " << rk.kv_width << " m (" << Kl.kv_layer << " triangles) max |d| = "
                   << Kl.dmax << " m; Prakash-Clifton L = 0.01 m max |d| = " << Q.dmax
                   << " m, L = 1e-200 against none " << dq << "\n";
      }
      TEST_ASSERT(K.finite && Kl.finite && K.prop1_bad == 0 && Kl.prop1_bad == 0 && Kl.kv_layer > 0
                  && Kl.kv_layer < 4 * kNf * kNf,
                  "N5 Kelvin-Voigt (whole domain, layer) runs finite with Proposition 1");
      TEST_ASSERT(Q.finite && Q.prop1_bad == 0 && dq == 0.0,
                  "N5 Prakash-Clifton runs finite; L = 1e-200 is bit for bit the run without it");
   }
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   std::set<std::string> only;
   bool np_mode = false;
   for (int i = 1; i < argc; i++)
   {
      const std::string a = argv[i];
      if (a == "--np") { np_mode = true; }
      else if (a == "--only" && i + 1 < argc)
      {
         std::string s = argv[++i];
         size_t b = 0;
         while (b <= s.size())
         {
            const size_t e = s.find(',', b);
            only.insert(s.substr(b, e == std::string::npos ? std::string::npos : e - b));
            if (e == std::string::npos) { break; }
            b = e + 1;
         }
      }
   }
   auto run = [&](const char *part) {
      if (!only.empty()) { return only.count(part) > 0; }
      if (np_mode) { return std::string(part) == "N3" || std::string(part) == "N5"; }
      return true;
   };
   if (g_rank == 0)
   {
      std::cout << "=== 2D Nitsche operator (2D in-plane design rev 3, phase C) on "
                << Mpi::WorldSize() << " rank(s) ===\n";
   }
   if (run("N2")) { TestCFL(); }
   if (only.count("N2M")) { TestCFLMeshes(); }
   if (run("N3")) { for (int p = 1; p <= 3; p++) { TestStep(p); } }
   if (run("N1"))
   {
      for (int p = 1; p <= 3; p++) { TestPlaneWave(p, true); TestPlaneWave(p, false); }
   }
   if (run("N4"))
   {
      for (int p = 1; p <= 3; p++) { TestDampedMass(p); }
      TestReflection();
   }
   if (run("N5")) { for (int p = 1; p <= 2; p++) { TestFault(p); } }
   if (g_rank == 0)
   {
      std::cout << "\n=== Summary: " << num_tests - num_failed << "/" << num_tests << " passed ===\n";
   }
   return num_failed == 0 ? 0 : 1;
}
