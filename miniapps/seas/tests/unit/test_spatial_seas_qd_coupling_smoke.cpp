// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_spatial_seas_qd_coupling_smoke.cpp — Phase 5 of
// document/spatial_seas_dev/PLAN_spatial_seas_quasidynamic_driver_2026-05-31.md
//
// Tiny-fixture local smoke for the quasi-dynamic coupling CONSTRUCTION ORDER
// the QD driver wires (R-001): geom(false) → SetRateStatePerDOF →
// RateStateFaultOperator ctor (asserts !geom.HasParams()) → ComputeParams*
// (the stress source; here ComputeParamsFaultLocal).  Asserts the HasParams
// transitions and that the stress source populated per-DOF tau_pre_/σ_n.
//
// SetSAFSMode → SetInitialCondition → RK45 are NOT exercised here: SetSAFSMode
// requires a SAFS-CLEAN fault-DOF basis (no zero-normal / t1-degenerate DOFs),
// and the small Cartesian elasticity-test fixture has degenerate fault-plane-
// edge DOFs that the guard correctly refuses.  Only the BP5 / SAFS production
// meshes are clean, so the equilibrium<1e-6 4-phase init + the 50-step run are
// the Phase-5 Frontera acceptance (project memory: no local production-mesh
// runs).  Uses CG_AMG (MUMPS bus-errors on tiny meshes — Phase 4).  np=1, np=4.

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
}  // namespace

int main(int argc, char *argv[])
{
#ifdef MFEM_USE_MPI
   MPI_Init(&argc, &argv);
   MPI_Comm comm = MPI_COMM_WORLD;
   MPI_Comm_rank(comm, &g_rank);
#else
#  error "test_spatial_seas_qd_coupling_smoke requires MFEM_USE_MPI=YES."
#endif

   int global_fail = 0;
   {
   const real_t Lx = 2.0, Ly = 2.0, Lz = 2.0;
   BP5Params seed;
   const real_t Vp = seed.Vp;
   const real_t Wf = Lz, lf = 2.0 * Lx;

   Mesh smesh = test::CreateTestMesh3D(2, 2, 2, Lx, Ly, Lz);
   ParMesh pmesh(comm, smesh);
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

   // ---- R-001 construction order: geom(false) → SetRateStatePerDOF ----
   FaultGeometry<ParMesh> geom(domain, seed, &mpi, /*compute_bp5_params=*/false);

   // Uniform rate-state (built directly; no resolver).  a-b < 0 (VW core).
   const real_t a_u = 0.010, b_u = 0.015, Dc_u = 0.004;
   const real_t f0_u = 0.6,   V0_u = 1.0e-6;
   const real_t Vinit_u = V0_u;            // reference state ⇒ f_ss = f0
   const real_t eta_u = seed.eta();        // radiation damping mu/(2 cs)
   const real_t sigma_n_u = 50.0e6;        // effective normal stress
   spatial::RateStatePerDOFParams rs;
   auto fill = [&](Vector &v, real_t val) { v.SetSize(n_owned); v = val; };
   fill(rs.a, a_u); fill(rs.b, b_u); fill(rs.Dc, Dc_u);
   fill(rs.eta, eta_u); fill(rs.V_init, Vinit_u);
   fill(rs.sigma_n_eff, sigma_n_u);
   fill(rs.f_0, f0_u); fill(rs.V_0, V0_u);

   Vector init_vel_dir(2); init_vel_dir(0) = 0.0; init_vel_dir(1) = 1.0;
   geom.SetRateStatePerDOF(rs, init_vel_dir);
   Check(!geom.HasParams(),
         "geom.HasParams()==false after SetRateStatePerDOF (R-001)");

   // ---- friction + fault op (BEFORE the stress source, R-001) ----
   DieterichRuinaFriction::Constants fc;
   fc.V0 = V0_u; fc.f0 = f0_u; fc.b = b_u; fc.Dc = Dc_u;
   DieterichRuinaFriction friction(fc);
   AgingLawPsi aging(fc.b, fc.V0, fc.f0);
   RateStateFaultOperator<ParMesh, 2> fault_op(&geom, &friction, &aging,
                                               seed, &mpi);

   // ---- stress source: reference-state-consistent uniform prestress ----
   // tau_strike ≈ sigma_n*f0 + eta*V_init (steady state at V=V0) so the
   // 4-phase psi solve equilibrates.
   const real_t tau_strike_u = sigma_n_u * f0_u + eta_u * Vinit_u;
   geom.ComputeParamsFaultLocal(tau_strike_u, /*tau_dip=*/0.0,
                                sigma_n_u, /*P_p=*/0.0);
   Check(geom.HasParams(), "geom.HasParams()==true after ComputeParamsFaultLocal");
   Check(geom.GetTauPre().Size() == 2 * n_owned,
         "tau_pre_ size == 2*n_owned after ComputeParamsFaultLocal");
   Check(geom.sigma_n_per_dof().Size() == n_owned,
         "sigma_n_per_dof_ size == n_owned");
   // The stress source overwrote σ_n with the FaultLocalPrestress value.
   {
      const Vector &sn = geom.sigma_n_per_dof();
      real_t sn_err = 0.0;
      for (int i = 0; i < sn.Size(); i++)
      { sn_err = std::max(sn_err, std::abs(sn(i) - sigma_n_u)); }
      Check(sn_err <= 1e-6 * sigma_n_u, "per-DOF effective σ_n == sigma_n_u");
   }

   // NOTE: SetSAFSMode → SetInitialCondition → RK45 are intentionally NOT
   // exercised here.  SetSAFSMode requires a SAFS-CLEAN fault-DOF basis
   // (NumZeroNormalFallbacks()==0 && NumT1Fallbacks()==0), and the small
   // Cartesian elasticity-test fixture has degenerate fault-plane-edge DOFs
   // (1 zero-normal here) that the SAFS guard correctly refuses.  Only the
   // BP5 / SAFS production meshes provide a clean basis, so the equilibrium<1e-6
   // 4-phase init + the 50-step run are the Phase-5 Frontera acceptance (the
   // QD driver wires that exact path; see the Frontera Stage-D sbatch).
   const int n_zero_normal = geom.NumZeroNormalFallbacks();
   if (g_rank == 0 && n_zero_normal > 0)
   {
      std::cout << "[note] tiny Cartesian fixture has " << n_zero_normal
                << " zero-normal fault DOF(s) — SAFS-mode init/stepping is the "
                << "BP5-mesh Frontera test (construction order validated here).\n";
   }

   int gf = g_fail, gc = g_checks;
   MPI_Allreduce(MPI_IN_PLACE, &gf, 1, MPI_INT, MPI_SUM, comm);
   MPI_Allreduce(MPI_IN_PLACE, &gc, 1, MPI_INT, MPI_SUM, comm);
   global_fail = gf;
   if (g_rank == 0)
   {
      std::cout << "=== test_spatial_seas_qd_coupling_smoke ===\n"
                << "  R-001 construction order validated (geom→rate-state→"
                << "fault-op→ComputeParamsFaultLocal).\n"
                << "  checks (all ranks): " << gc << ", failures: "
                << global_fail << "\n"
                << (global_fail == 0 ? "PASS\n" : "FAIL\n");
   }
   }

   MPI_Finalize();
   return (global_fail == 0) ? 0 : 1;
}
