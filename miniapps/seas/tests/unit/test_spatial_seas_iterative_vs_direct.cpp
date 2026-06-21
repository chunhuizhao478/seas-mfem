// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_spatial_seas_iterative_vs_direct.cpp — Phase 4 of
// document/spatial_seas_dev/PLAN_spatial_seas_quasidynamic_driver_2026-05-31.md
//
// Solver-equivalence acceptance test for the Phase-4 AMG-preconditioned Krylov
// dispatch.  Two ElasticityDomainOperators are built on the SAME ParMesh with
// the SAME slip + Dirichlet RHS but DIFFERENT solver setups, and their solved
// displacement + recovered fault traction must agree:
//
//   ||u_CGAMG - u_ref|| / ||u_ref||           < 1e-7   and
//   ||trac_CGAMG - trac_ref|| / ||trac_ref||  < 1e-7   (global L2).
//
// LOCAL reference = SolverType::GMRES_AMG (a DIFFERENT preconditioner setup:
// CG_AMG uses SetSystemsOptions(3,byNODES)+SetElasticityOptions, GMRES_AMG uses
// SetSystemsOptions only).  Both Krylov methods converge to the UNIQUE solution
// u* = K^{-1} b (K is SPD), so agreement to a tight tolerance validates that the
// CG_AMG path — including the R-006 SetElasticityOptions-on-DG decision — solves
// correctly (not just that two identical paths agree).
//
// DEVIATION from the plan's literal "CG_AMG vs MUMPS direct": MUMPS bus-errors
// inside dmumps_scatter_dist_rhs_ on small/degenerate matrices (verified on the
// 2x2x2 Cartesian fixture, np=1 and np=4) — it is stable only on the production
// BP5 mesh, which is Frontera-only (project memory: no local production-mesh
// runs).  The CG_AMG-vs-MUMPS comparison on the BP5 mesh is therefore part of
// the Phase-4 Frontera perfgraph/scaling job (Stage C), not this local unit
// test.  The MUMPS path is untouched by Phase 4, so this is purely a
// test-harness reference choice.
//
// DG spaces have no rank-shared DOFs (GetVSize == GetTrueVSize), so a local
// sum-of-squares + MPI_Allreduce(SUM) gives the correct global L2 norm.

#include "mfem.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../config/bp5_params.hpp"
#include "fault_mesh_fixture.hpp"   // test::CreateTestMesh3D

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

namespace
{
int g_rank = 0;
int g_fail = 0;
int g_checks = 0;

void Check(bool ok, const std::string &msg)
{
   g_checks++;
   if (!ok)
   {
      g_fail++;
      std::cerr << "[rank " << g_rank << "] FAILED: " << msg << "\n";
   }
}

// Global relative L2 of (a - b) vs b, over a no-overlap (DG) vector pair.
real_t GlobalRelL2(const Vector &a, const Vector &b, MPI_Comm comm)
{
   MFEM_VERIFY(a.Size() == b.Size(),
               "GlobalRelL2: size mismatch " << a.Size() << " vs " << b.Size());
   real_t ld = 0.0, lb = 0.0;
   const int n = a.Size();
   for (int i = 0; i < n; i++)
   {
      const real_t d = a(i) - b(i);
      ld += d * d;
      lb += b(i) * b(i);
   }
   real_t gd = 0.0, gb = 0.0;
   MPI_Allreduce(&ld, &gd, 1, MPI_DOUBLE, MPI_SUM, comm);
   MPI_Allreduce(&lb, &gb, 1, MPI_DOUBLE, MPI_SUM, comm);
   return (gb > 0.0) ? std::sqrt(gd) / std::sqrt(gb) : std::sqrt(gd);
}
}  // namespace

int main(int argc, char *argv[])
{
#ifdef MFEM_USE_MPI
   MPI_Init(&argc, &argv);
   MPI_Comm comm = MPI_COMM_WORLD;
   MPI_Comm_rank(comm, &g_rank);
#else
#  error "test_spatial_seas_iterative_vs_direct requires MFEM_USE_MPI=YES."
#endif

   int global_fail = 0;
   {
   const real_t Lx = 2.0, Ly = 2.0, Lz = 2.0;
   BP5Params bp5;
   const real_t lambda = bp5.lambda(), mu = bp5.mu(), Vp = bp5.Vp;
   const real_t Wf = Lz, lf = 2.0 * Lx;
   const real_t time = 1.0e7;   // > 0 ⇒ Dirichlet plate loading active

   Mesh smesh = test::CreateTestMesh3D(2, 2, 2, Lx, Ly, Lz);
   ParMesh pmesh(comm, smesh);
   smesh.Clear();

   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.dirichlet_attrs = {5};
   bc.default_dirichlet_func = MakeBP5DirichletFunc(Vp);

   LinearElastic le(lambda, mu);

   auto make_cfg = [](real_t rtol)
   {
      DomainConfig c;
      c.ksp_rtol = rtol;
      c.ksp_atol = 0.0;
      c.ksp_maxit = 5000;
      c.amg_print_level = 0;
      c.check_residual = true;   // mandatory residual check on the AMG paths
      return c;
   };

   // (A) CG + BoomerAMG (SetSystemsOptions + SetElasticityOptions, R-006).
   ElasticityDomainOperator<ParMesh> op_cg(
      pmesh, /*order=*/1, le, Vp, Wf, lf, bc,
      DGMethod::IP, SolverType::CG_AMG, make_cfg(1.0e-12));

   // (B) GMRES + BoomerAMG (SetSystemsOptions only) — independent reference.
   ElasticityDomainOperator<ParMesh> op_ref(
      pmesh, /*order=*/1, le, Vp, Wf, lf, bc,
      DGMethod::IP, SolverType::GMRES_AMG, make_cfg(1.0e-12));

   const int nfd = op_cg.GetNumFaultDOFs();
   Check(op_cg.GetFESpace().GetVSize() == op_ref.GetFESpace().GetVSize(),
         "CG and GMRES operators share FE-space size (same mesh)");

   Vector slip_bc(2 * nfd);
   for (int i = 0; i < nfd; i++)
   {
      slip_bc(2 * i)     = 5.0e-4;   // dip
      slip_bc(2 * i + 1) = 1.0e-3;   // strike
   }

   ParGridFunction u_cg(&op_cg.GetFESpace());
   ParGridFunction u_ref(&op_ref.GetFESpace());
   u_cg = 0.0; u_ref = 0.0;
   op_cg.Solve(time, slip_bc, u_cg);
   op_ref.Solve(time, slip_bc, u_ref);

   const real_t u_rel = GlobalRelL2(u_cg, u_ref, comm);
   Check(u_rel < 1.0e-7,
         "||u_CGAMG - u_GMRESAMG|| / ||u_GMRESAMG|| < 1e-7 (solver equivalence)");

   Vector trac_cg, trac_ref;
   op_cg.ComputeTraction(u_cg, slip_bc, trac_cg);
   op_ref.ComputeTraction(u_ref, slip_bc, trac_ref);
   Check(trac_cg.Size() == trac_ref.Size(), "traction sizes match");
   const real_t t_rel = GlobalRelL2(trac_cg, trac_ref, comm);
   Check(t_rel < 1.0e-7,
         "||trac_CGAMG - trac_GMRESAMG|| / ||trac_GMRESAMG|| < 1e-7 (fault traction)");

   // R-402: preconditioner reuse — a 2nd solve (different loading time) must
   // NOT rebuild K or the AMG setup (guarded by stiffness_assembled_).
   ParGridFunction u_cg2(&op_cg.GetFESpace());
   u_cg2 = 0.0;
   op_cg.Solve(2.0 * time, slip_bc, u_cg2);
   Check(op_cg.NumStiffnessAssemblies() == 1,
         "AMG/solver setup runs exactly once over 2 solves (preconditioner reuse)");

   int gf = g_fail, gc = g_checks, gnfd = nfd;
   MPI_Allreduce(MPI_IN_PLACE, &gf, 1, MPI_INT, MPI_SUM, comm);
   MPI_Allreduce(MPI_IN_PLACE, &gc, 1, MPI_INT, MPI_SUM, comm);
   MPI_Allreduce(MPI_IN_PLACE, &gnfd, 1, MPI_INT, MPI_SUM, comm);
   global_fail = gf;
   if (g_rank == 0)
   {
      std::cout << "=== test_spatial_seas_iterative_vs_direct ===\n"
                << "  global fault DOFs: " << gnfd << "\n"
                << "  ||du||/||u|| = " << u_rel << ", ||dT||/||T|| = " << t_rel << "\n"
                << "  (local ref = GMRES_AMG; CG_AMG-vs-MUMPS on BP5 mesh = Frontera)\n"
                << "  checks (all ranks): " << gc << ", failures: " << global_fail << "\n"
                << (global_fail == 0 ? "PASS\n" : "FAIL\n");
   }
   }

   MPI_Finalize();
   return (global_fail == 0) ? 0 : 1;
}
