// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_elasticity_heterogeneous_ctor.cpp — Phase 2b Stage 1b of
// document/spatial_seas_dev/PLAN_spatial_seas_quasidynamic_driver_2026-05-31.md
//
// Proves the new heterogeneous `ElasticityDomainOperator(Coefficient& lambda,
// Coefficient& mu, ...)` ctor (IP method) reproduces the existing constant
// `LinearElastic(lambda, mu)` ctor BIT-FOR-BIT when the external coefficients
// are constant.  This validates the Stage-1a `MaterialCoefficient` delegation
// path end-to-end: external Coefficient -> MaterialCoefficient -> the DG
// integrators' per-qp `Eval(T, ip)` (ElasticityIntegrator for the volume K;
// DGElasticityIPCombinedIntegrator for the slip/Dirichlet RHS and the traction
// recovery).  A genuinely spatially-varying coefficient is exercised by the
// Phase-2b Stage-3 2-layer slab test.
//
// `MaterialCoefficient` in external mode delegating to `ConstantCoefficient(c)`
// returns exactly `c` per qp — the same double the constant ctor's constant
// mode returns — so every assembled entry, the CG solve, and the recovered
// traction must be bit-identical (max abs diff == 0).
//
// Serial (Mesh) fixture, no MPI — mirrors test_elasticity_operator.cpp.

#include "mfem.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../config/bp5_params.hpp"
#include "fault_mesh_fixture.hpp"   // test::CreateTestMesh3D (R-503)

#include <cmath>
#include <iostream>
#include <string>

using namespace mfem;
using namespace mfem::seas;

namespace
{
int g_fail = 0;
int g_checks = 0;

void Check(bool ok, const std::string &msg)
{
   g_checks++;
   if (!ok) { g_fail++; std::cerr << "FAILED: " << msg << "\n"; }
}

// Max |a-b| over a Vector pair; 1e300 if sizes differ.
real_t MaxAbsDiff(const Vector &a, const Vector &b)
{
   if (a.Size() != b.Size()) { return 1e300; }
   real_t m = 0.0;
   for (int i = 0; i < a.Size(); i++) { m = std::max(m, std::abs(a(i) - b(i))); }
   return m;
}

using test::CreateTestMesh3D;   // shared fault-mesh fixture (R-503)

BoundaryConfig MakeBC(real_t Vp)
{
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.dirichlet_attrs = {5};
   bc.default_dirichlet_func = MakeBP5DirichletFunc(Vp);
   return bc;
}
}  // namespace

int main()
{
   const real_t Lx = 2.0, Ly = 2.0, Lz = 2.0;
   BP5Params bp5;
   const real_t lambda = bp5.lambda();
   const real_t mu     = bp5.mu();
   const real_t Vp     = bp5.Vp;
   const real_t Wf = Lz, lf = 2.0 * Lx;
   const real_t time = 1.0e7;   // > 0 so Dirichlet plate loading is active

   // Two meshes (each operator owns/curves its mesh independently, but the
   // structure — and thus DOF numbering — is identical).
   Mesh mesh_c = CreateTestMesh3D(2, 2, 2, Lx, Ly, Lz);
   Mesh mesh_h = CreateTestMesh3D(2, 2, 2, Lx, Ly, Lz);

   BoundaryConfig bc_c = MakeBC(Vp);
   BoundaryConfig bc_h = MakeBC(Vp);
   DomainConfig dcfg;   // defaults (IP, residual_check etc.)

   // (A) constant ctor (LinearElastic) — the reference path.
   LinearElastic le(lambda, mu);
   ElasticityDomainOperator<Mesh> op_const(
      mesh_c, /*order=*/1, le, Vp, Wf, lf, bc_c,
      DGMethod::IP, SolverType::CG_AMG, dcfg);

   // (B) heterogeneous ctor with CONSTANT external coefficients — exercises the
   //     MaterialCoefficient -> external delegation through the IP integrators.
   ConstantCoefficient lam_c(lambda), mu_c(mu);
   ElasticityDomainOperator<Mesh> op_het(
      mesh_h, /*order=*/1, lam_c, mu_c, Vp, Wf, lf, bc_h,
      DGMethod::IP, SolverType::CG_AMG, dcfg);

   const int nfd_c = op_const.GetNumFaultDOFs();
   const int nfd_h = op_het.GetNumFaultDOFs();
   Check(nfd_c > 0, "constant op has fault DOFs");
   Check(nfd_c == nfd_h, "both ctors detect the same fault-DOF count");
   // R-504: the bit-for-bit Vector comparisons below are index-by-index and
   // therefore rely on identical FE-space SIZE *and* DOF ordering between the
   // two operators.  That holds because both meshes are built by the same
   // deterministic CreateTestMesh3D call; assert the size equality explicitly
   // so a future ordering/numbering change surfaces here rather than as a
   // spurious value mismatch.
   Check(op_const.GetFESpace().GetVSize() == op_het.GetFESpace().GetVSize(),
         "both operators have identical FE-space size (identical meshes)");

   if (nfd_c == nfd_h && nfd_c > 0)
   {
      // Non-trivial slip (dip + strike) to exercise both components.
      Vector slip_bc(2 * nfd_c);
      for (int i = 0; i < nfd_c; i++)
      {
         slip_bc(2 * i)     = 5.0e-4;   // dip
         slip_bc(2 * i + 1) = 1.0e-3;   // strike
      }

      // (1) slip-only RHS (DGElasticityIPCombinedIntegrator::AssembleSlipFaceRHS)
      Vector rhs_slip_c, rhs_slip_h;
      op_const.AssembleSlipOnlyRHS(rhs_slip_c, slip_bc);
      op_het.AssembleSlipOnlyRHS(rhs_slip_h, slip_bc);
      Check(MaxAbsDiff(rhs_slip_c, rhs_slip_h) == 0.0,
            "slip-only RHS bit-for-bit (slip integrator, coeff path)");

      // (2) Dirichlet-loading RHS (boundary face integrator)
      Vector rhs_dir_c, rhs_dir_h;
      op_const.AssembleDirichletOnlyRHS(rhs_dir_c, time);
      op_het.AssembleDirichletOnlyRHS(rhs_dir_h, time);
      Check(MaxAbsDiff(rhs_dir_c, rhs_dir_h) == 0.0,
            "Dirichlet-loading RHS bit-for-bit (boundary integrator, coeff path)");

      // (3) full Solve K u = b(slip, Dirichlet@time) — exercises the volume
      //     stiffness K (ElasticityIntegrator, coeff path) + the CG solve.
      GridFunction u_c(&op_const.GetFESpace());
      GridFunction u_h(&op_het.GetFESpace());
      u_c = 0.0; u_h = 0.0;
      op_const.Solve(time, slip_bc, u_c);
      op_het.Solve(time, slip_bc, u_h);
      Check(u_c.Size() == u_h.Size(), "solved displacement same size");
      Check(MaxAbsDiff(u_c, u_h) == 0.0,
            "solved displacement bit-for-bit (volume K + RHS + CG, coeff path)");

      // (4) traction recovery (ComputeTractionImpl IP branch -> trac integrator)
      Vector trac_c, trac_h;
      op_const.ComputeTraction(u_c, slip_bc, trac_c);
      op_het.ComputeTraction(u_h, slip_bc, trac_h);
      Check(trac_c.Size() == 2 * nfd_c, "traction size == 2*nfd");
      Check(MaxAbsDiff(trac_c, trac_h) == 0.0,
            "recovered fault traction bit-for-bit (traction integrator, coeff path)");
   }

   std::cout << "=== test_elasticity_heterogeneous_ctor ===\n"
             << "  fault DOFs: " << nfd_c << "\n"
             << "  checks: " << g_checks << ", failures: " << g_fail << "\n"
             << (g_fail == 0 ? "PASS\n" : "FAIL\n");
   return (g_fail == 0) ? 0 : 1;
}
