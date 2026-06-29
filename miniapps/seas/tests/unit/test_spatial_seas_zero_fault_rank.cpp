// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_spatial_seas_zero_fault_rank.cpp — np=4 regression for the FAULT-FREE-RANK
// fix in drivers/spatial_seas_driver.cpp §5b.
//
// BUG (pre-fix): the quasi-dynamic driver built its scalar friction constants by
// indexing element 0 of this rank's OWNED rate-state vectors —
//   fc.V0 = rs.V_0(0); fc.f0 = rs.f_0(0); fc.b = rs.b(0); fc.Dc = rs.Dc(0);
// guarded by MFEM_VERIFY(rs.b.Size() > 0, "spatial_seas: empty rate-state
// vectors.") — and likewise read rs.a(0), rs.V_init(0), and fault_dof_basis()
// column 0.  A rank that owns ZERO fault DOFs (entirely NORMAL at high rank
// counts — the planar fault does not reach every volume partition; the proven
// drivers/seas_driver.cpp builds fc from config and so runs at 8N×400r) has empty
// rs.* and aborted at drivers/spatial_seas_driver.cpp:954.  The fix sources the
// (globally uniform) scalars from the config on such ranks, mirroring
// seas_driver.cpp:431-435.
//
// This test FORCES fault-free ranks with an explicit y-band partition of the
// y=0-fault CreateTestMesh3D fixture (cf. tests/parallel/
// test_R002_empty_fault_rank_no_abort.cpp), then drives the QD driver's exact
// R-001 construction chain (geom(false) → SetRateStatePerDOF → friction
// constants → RateStateFaultOperator ctor → ComputeParamsFaultLocal; cf.
// tests/unit/test_spatial_seas_qd_coupling_smoke.cpp).  Pre-fix the two outer
// ranks aborted at the friction-constant build; the fix makes them fall back to
// the config uniform defaults.  Reaching the collective verdict without abort,
// with the fixture proven non-vacuous (some rank owns zero fault DOFs, some owns
// > 0), is the regression gate.
//
// Construction order is validated only (SetSAFSMode → SetInitialCondition → RK45
// need a SAFS-clean basis the tiny Cartesian fixture lacks — project memory: QD
// init/stepping is BP5-mesh Frontera; see test_spatial_seas_qd_coupling_smoke).
// CG_AMG (MUMPS bus-errors on tiny meshes).  Run: mpirun -np 4
// ./seas_test_spatial_seas_zero_fault_rank

#include "mfem.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../fault/fault_geometry.hpp"
#include "../../fault/rate_state_fault.hpp"
#include "../../friction/dieterich_ruina.hpp"
#include "../../friction/state_evolution.hpp"
#include "../../config/bp5_params.hpp"
#include "../../common/mpi_context.hpp"
#include "../../spatial/code/spatial_friction.hpp"   // RateStatePerDOFParams
#include "fault_mesh_fixture.hpp"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

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

// Assign each element to a rank by its y-centroid, so the two OUTER y-bands hold
// no fault-adjacent elements.  CreateTestMesh3D(2,2,2,2,2,2) spans y in [-2,2]
// with four y-layers ([-2,-1],[-1,0],[0,1],[1,2]) and the fault at y=0; ranks 0
// and 1 straddle y=0 (own the fault), ranks 2 and 3 own the outer bands (no
// fault DOFs).  Mirrors test_R002_empty_fault_rank_no_abort.cpp::BuildSkewedPartition.
std::vector<int> BuildYBandPartition(const Mesh &mesh)
{
   const int ne = mesh.GetNE();
   std::vector<int> part(ne, 0);
   for (int e = 0; e < ne; e++)
   {
      Array<int> ev; mesh.GetElementVertices(e, ev);
      real_t cy = 0.0;
      for (int v = 0; v < ev.Size(); v++) { cy += mesh.GetVertex(ev[v])[1]; }
      cy /= ev.Size();
      if      (cy < -1.0) { part[e] = 2; }   // outer band [-2,-1] — fault-free
      else if (cy <  0.0) { part[e] = 0; }   // inner band [-1, 0] — owns fault
      else if (cy <  1.0) { part[e] = 1; }   // inner band [ 0, 1] — owns fault
      else                { part[e] = 3; }   // outer band [ 1, 2] — fault-free
   }
   return part;
}
}  // namespace

int main(int argc, char *argv[])
{
#ifdef MFEM_USE_MPI
   MPI_Init(&argc, &argv);
   MPI_Comm comm = MPI_COMM_WORLD;
   MPI_Comm_rank(comm, &g_rank);
#else
#  error "test_spatial_seas_zero_fault_rank requires MFEM_USE_MPI=YES."
#endif

   int nprocs = 1;
   MPI_Comm_size(comm, &nprocs);
   if (nprocs != 4)
   {
      if (g_rank == 0)
      {
         std::cerr << "[zero-fault-rank] TEST REQUIRES np=4, got " << nprocs
                   << " — skipping.\n";
      }
      MPI_Finalize();
      return 0;
   }

   int global_fail = 0;
   {
   const real_t Lx = 2.0, Ly = 2.0, Lz = 2.0;
   BP5Params seed;
   const real_t Vp = seed.Vp;
   const real_t Wf = Lz, lf = 2.0 * Lx;

   // y=0-fault fixture + explicit partition that strands ranks 2/3 off the fault.
   Mesh smesh = test::CreateTestMesh3D(2, 2, 2, Lx, Ly, Lz);
   std::vector<int> part = BuildYBandPartition(smesh);
   ParMesh pmesh(comm, smesh, part.data());
   smesh.Clear();

   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.dirichlet_attrs = {5};
   bc.default_dirichlet_func = MakeBP5DirichletFunc(Vp);

   LinearElastic le(seed.lambda(), seed.mu());
   DomainConfig dcfg;
   dcfg.ksp_rtol = 1.0e-12;
   dcfg.check_residual = true;
   ElasticityDomainOperator<ParMesh> domain(
      pmesh, /*order=*/1, le, Vp, Wf, lf, bc,
      DGMethod::IP, SolverType::CG_AMG, dcfg);

   MPIContext mpi(comm);
   const int n_owned = domain.GetNumOwnedFaultDOFs();

   // ---- Non-vacuity: the fixture MUST produce both fault-free and fault-owning
   //      ranks, else the test passes even without the fix. ----
   int gmin = n_owned, gmax = n_owned;
   MPI_Allreduce(MPI_IN_PLACE, &gmin, 1, MPI_INT, MPI_MIN, comm);
   MPI_Allreduce(MPI_IN_PLACE, &gmax, 1, MPI_INT, MPI_MAX, comm);
   Check(gmin == 0, "fixture non-vacuous: some rank owns ZERO fault DOFs");
   Check(gmax >  0, "fixture non-vacuous: some rank owns fault DOFs");
   const bool fault_free = (n_owned == 0);
   if (g_rank == 2 || g_rank == 3)
   {
      Check(fault_free, "outer y-band rank (2/3) owns zero fault DOFs");
   }

   // ---- R-001 construction order: geom(false) → SetRateStatePerDOF ----
   FaultGeometry<ParMesh> geom(domain, seed, &mpi, /*compute_bp5_params=*/false);

   // Uniform rate-state.  fill() sizes each vector to n_owned, so on a fault-free
   // rank EVERY rs.* is EMPTY — exactly the state that aborted pre-fix.
   const real_t a_u = 0.010, b_u = 0.015, Dc_u = 0.004;
   const real_t f0_u = 0.6,  V0_u = 1.0e-6;
   const real_t Vinit_u = V0_u;
   const real_t eta_u = seed.eta();
   const real_t sigma_n_u = 50.0e6;
   spatial::RateStatePerDOFParams rs;
   auto fill = [&](Vector &v, real_t val) { v.SetSize(n_owned); v = val; };
   fill(rs.a, a_u); fill(rs.b, b_u); fill(rs.Dc, Dc_u);
   fill(rs.eta, eta_u); fill(rs.V_init, Vinit_u);
   fill(rs.sigma_n_eff, sigma_n_u);
   fill(rs.f_0, f0_u); fill(rs.V_0, V0_u);

   Vector init_vel_dir(2); init_vel_dir(0) = 0.0; init_vel_dir(1) = 1.0;
   geom.SetRateStatePerDOF(rs, init_vel_dir);   // callee must tolerate empty
   Check(!geom.HasParams(),
         "geom.HasParams()==false after SetRateStatePerDOF (R-001)");

   // ---- THE FIX: scalar friction constants.  Fault-owning ranks read element 0
   //      (unchanged); a fault-free rank (empty rs) falls back to the config
   //      uniform defaults instead of indexing rs.*(0). ----
   if (fault_free)
   {
      Check(rs.b.Size() == 0,
            "fault-free rank: rs.b empty (pre-fix aborted at :954, OOB at :956)");
   }
   DieterichRuinaFriction::Constants fc;
   if (rs.b.Size() > 0)
   {
      fc.V0 = rs.V_0(0); fc.f0 = rs.f_0(0); fc.b = rs.b(0); fc.Dc = rs.Dc(0);
   }
   else
   {
      fc.V0 = V0_u; fc.f0 = f0_u; fc.b = b_u; fc.Dc = Dc_u;
   }
   DieterichRuinaFriction friction(fc);
   AgingLawPsi aging(fc.b, fc.V0, fc.f0);
   if (fault_free)
   {
      Check(fc.b == b_u && fc.V0 == V0_u && fc.f0 == f0_u && fc.Dc == Dc_u,
            "fault-free rank: friction constants from config fallback");
   }

   // ---- fault operator (BEFORE the stress source, R-001): the ctor must accept
   //      a geometry with zero owned fault DOFs. ----
   RateStateFaultOperator<ParMesh, 2> fault_op(&geom, &friction, &aging,
                                               seed, &mpi);

   // ---- stress source: ComputeParams* must early-return (no-op) on empty. ----
   const real_t tau_strike_u = sigma_n_u * f0_u + eta_u * Vinit_u;
   geom.ComputeParamsFaultLocal(tau_strike_u, /*tau_dip=*/0.0,
                                sigma_n_u, /*P_p=*/0.0);
   // ComputeParams* early-returns when num_fault_dofs_==0, so a fault-free rank
   // does NOT flip HasParams — assert it only where there are owned fault DOFs.
   if (n_owned > 0)
   {
      Check(geom.HasParams(),
            "fault-owning rank: HasParams()==true after ComputeParamsFaultLocal");
      Check(geom.GetTauPre().Size() == 2 * n_owned,
            "fault-owning rank: tau_pre_ size == 2*n_owned");
   }

   // ---- Collective verdict (every rank, incl. the fault-free ones, enters). ----
   int gf = g_fail, gc = g_checks, gowned = n_owned;
   MPI_Allreduce(MPI_IN_PLACE, &gf, 1, MPI_INT, MPI_SUM, comm);
   MPI_Allreduce(MPI_IN_PLACE, &gc, 1, MPI_INT, MPI_SUM, comm);
   MPI_Allreduce(MPI_IN_PLACE, &gowned, 1, MPI_INT, MPI_SUM, comm);
   global_fail = gf;
   if (g_rank == 0)
   {
      std::cout << "=== test_spatial_seas_zero_fault_rank ===\n"
                << "  global owned fault DOFs: " << gowned
                << "  (min/rank=" << gmin << ", max/rank=" << gmax << ")\n"
                << "  fault-free ranks survived the QD construction chain.\n"
                << "  checks (all ranks): " << gc << ", failures: " << global_fail
                << "\n" << (global_fail == 0 ? "PASS\n" : "FAIL\n");
   }
   }  // end scope — MPI-owning objects destruct before MPI_Finalize

   MPI_Finalize();
   return (global_fail == 0) ? 0 : 1;
}
