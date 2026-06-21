// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_spatial_seas_faultgeom_parity.cpp — Phase 3 of
// document/spatial_seas_dev/PLAN_spatial_seas_quasidynamic_driver_2026-05-31.md
//
// Proves the spatial -> quasi-dynamic fault-geometry wiring + owned-DOF
// ordering on the BP5 mesh, two ways:
//
//   PART A (1e-10 parity, R-007 split).  Build FaultGeometry two ways:
//     (a) the analytic BP5 ctor  (domain, params, &mpi)  -> ComputeBP5Params;
//     (b) the new ctor (domain, seed, &mpi, false) + SetRateStatePerDOF, with
//         the per-DOF rate-state values evaluated from bp5_params at the
//         OWNED coordinates produced by the NEW owned-order getters
//         (GetFaultDOFCoords3D + RestrictToOwnedFault).
//   Assert a_values_ / dc_values_ / eta_values_ / V_init_vec_ agree to 1e-10.
//   Because path (b) reads the GETTER-derived owned coords while path (a)
//   reads FaultGeometry's internal owned coords, agreement proves the getter
//   walk and the geometry walk produce the SAME owned-DOF ordering (R-002) and
//   that SetRateStatePerDOF maps rs -> geom arrays correctly (incl. the R-010
//   V_init dip/strike decomposition).  tau_pre_/sigma_n_per_dof_ are NOT
//   compared element-wise (ComputeBP5Params does not populate sigma_n_per_dof_;
//   its tau_pre_ is heterogeneous — R-007).
//
//   PART B (resolver wiring).  Run the REAL SpatialFrictionResolver::
//   ResolveRateState on the owned-order getter tables with a UNIFORM
//   [friction.rate_state] block, assert the per-DOF output is uniform with the
//   expected values (incl. auto-eta = 0.5*sqrt(mu*rho)), then SetRateStatePerDOF
//   it into a fresh FaultGeometry and assert the geom arrays match.  (The
//   resolver cannot reproduce BP5's 2-D a(x2,x3) boxcar, so the 1e-10 BP5
//   parity uses the direct-bp5 rs above; PART B exercises the resolver path.)
//
// Runs at np=1 and np=4 (acceptance: owned/shared layout consistent).  Each
// rank checks its OWNED fault DOFs; a global MPI_Allreduce(MAX) of the local
// failure count gives the collective verdict.

#include "mfem.hpp"

#include "../../domain/elasticity_operator.hpp"
#include "../../fault/fault_geometry.hpp"
#include "../../config/bp5_params.hpp"
#include "../../common/mpi_context.hpp"
#include "../../spatial/code/spatial_friction.hpp"

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

// Max |a-b| over a Vector pair (returns 0 for empty / size-0 ranks).
real_t MaxAbsDiff(const Vector &a, const Vector &b)
{
   if (a.Size() != b.Size()) { return 1e300; }
   real_t m = 0.0;
   for (int i = 0; i < a.Size(); i++)
   {
      m = std::max(m, std::abs(a(i) - b(i)));
   }
   return m;
}
}  // namespace

int main(int argc, char *argv[])
{
#ifdef MFEM_USE_MPI
   MPI_Init(&argc, &argv);
   MPI_Comm comm = MPI_COMM_WORLD;
   MPI_Comm_rank(comm, &g_rank);
#else
#  error "test_spatial_seas_faultgeom_parity requires MFEM_USE_MPI=YES."
#endif

   const char *mesh_path = "bp5/mesh/reference/bp5_tandem_exact.msh";
   const real_t tol = 1e-10;

   int global_fail = 0;
   // Scope so every MFEM/MPI-owning object (ParMesh, operator, FaultGeometry,
   // MPIContext) is destroyed BEFORE MPI_Finalize — otherwise their dtors
   // free MPI communicators after finalize (illegal per the MPI standard).
   {
   // ---- Mesh + elasticity operator (BP5 moduli + bc) ----
   BP5Params bp5;     // analytic BP5 parameters (also the "seed")
   Mesh smesh(mesh_path, 1, 1);
   MFEM_VERIFY(smesh.Dimension() == 3, "parity test: need a 3D BP5 mesh");
   ParMesh pmesh(comm, smesh);
   smesh.Clear();
   pmesh.SetCurvature(1);

   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.dirichlet_attrs = {5};
   bc.default_dirichlet_func = MakeBP5DirichletFunc(bp5.Vp);

   LinearElastic le(bp5.lambda(), bp5.mu());
   DomainConfig dcfg;
   ElasticityDomainOperator<ParMesh> domain(
      pmesh, /*order=*/1, le, bp5.Vp, bp5.Wf, bp5.lf,
      bc, DGMethod::IP, SolverType::CG_AMG, dcfg);

   MPIContext mpi(comm);
   const int n_owned = domain.GetNumOwnedFaultDOFs();

   // ---- Owned-order getter tables (the production wiring inputs) ----
   Array<int> dof_to_elem, dof_to_attr;
   domain.GetFaultDOFToElem(dof_to_elem);
   domain.GetFaultDOFToAttr(dof_to_attr);
   std::vector<IntegrationPoint> dof_ips;
   domain.GetFaultDOFIntegrationPoints(dof_ips);
   Vector local_coords_3d, dof_coords_3d;
   domain.GetFaultDOFCoords3D(local_coords_3d);
   domain.RestrictToOwnedFault(local_coords_3d, dof_coords_3d, /*comps=*/3);

   Check(dof_to_elem.Size() == n_owned, "GetFaultDOFToElem size == owned");
   Check(dof_to_attr.Size() == n_owned, "GetFaultDOFToAttr size == owned");
   Check(static_cast<int>(dof_ips.size()) == n_owned,
         "GetFaultDOFIntegrationPoints size == owned");
   Check(dof_coords_3d.Size() == 3 * n_owned,
         "owned dof_coords_3d size == 3*owned");
   for (int i = 0; i < dof_to_attr.Size(); i++)
   {
      Check(dof_to_attr[i] == bc.fault_attr, "dof_to_attr == fault_attr");
   }

   // =================================================================
   // PART A — 1e-10 parity (path a: ComputeBP5Params; path b: getter
   //          coords + direct-bp5 rs + SetRateStatePerDOF).
   // =================================================================
   FaultGeometry<ParMesh> geom_a(domain, bp5, &mpi);   // ComputeBP5Params
   FaultGeometry<ParMesh> geom_b(domain, bp5, &mpi, /*compute_bp5_params=*/false);

   Check(geom_a.NumFaultDOFs() == n_owned, "geom_a NumFaultDOFs == owned");
   Check(geom_b.NumFaultDOFs() == n_owned, "geom_b NumFaultDOFs == owned");
   Check(!geom_b.HasParams(),
         "geom_b HasParams()==false after compute=false ctor (R-001)");

   // Build rs from bp5_params evaluated at the GETTER-derived owned coords.
   spatial::RateStatePerDOFParams rs;
   rs.a.SetSize(n_owned);  rs.Dc.SetSize(n_owned);
   rs.eta.SetSize(n_owned); rs.V_init.SetSize(n_owned);
   rs.sigma_n_eff.SetSize(n_owned);
   rs.b.SetSize(n_owned); rs.f_0.SetSize(n_owned); rs.V_0.SetSize(n_owned);
   const real_t bp5_eta = bp5.eta();
   for (int i = 0; i < n_owned; i++)
   {
      const real_t x2 = dof_coords_3d(3 * i + 0);    // along-strike = X
      const real_t x3 = -dof_coords_3d(3 * i + 2);   // depth = -Z (Tandem)
      rs.a(i)   = bp5.a_of_x2_x3(x2, x3);
      rs.Dc(i)  = bp5.Dc_of_x2_x3(x2, x3);
      rs.eta(i) = bp5_eta;
      real_t Vi[2];
      bp5.V_init_vec(x2, x3, Vi);
      // BP5's dip velocity is a negligible floor (V_zero = 1e-20), far below
      // the 1e-10 parity tolerance, so the pure-strike (0,1) decomposition is
      // valid: the omitted dip introduces at most |V_zero| error.
      Check(std::abs(Vi[0]) <= tol,
            "bp5 V_init dip component negligible (<= parity tol)");
      rs.V_init(i)     = Vi[1];     // strike magnitude
      rs.sigma_n_eff(i) = bp5.sigma_n;
      rs.b(i) = bp5.b; rs.f_0(i) = bp5.f0; rs.V_0(i) = bp5.V0;
   }

   Vector init_vel_dir(2); init_vel_dir(0) = 0.0; init_vel_dir(1) = 1.0;  // BP5: strike
   geom_b.SetRateStatePerDOF(rs, init_vel_dir);

   Check(!geom_b.HasParams(),
         "geom_b HasParams()==false after SetRateStatePerDOF (R-001)");
   Check(MaxAbsDiff(geom_b.GetAValues(),  geom_a.GetAValues())  <= tol,
         "a_values_ parity (path a == path b)");
   Check(MaxAbsDiff(geom_b.GetDcValues(), geom_a.GetDcValues()) <= tol,
         "dc_values_ parity");
   Check(MaxAbsDiff(geom_b.GetEtaValues(), geom_a.GetEtaValues()) <= tol,
         "eta_values_ parity");
   Check(MaxAbsDiff(geom_b.GetVInit(),    geom_a.GetVInit())     <= tol,
         "V_init_vec_ parity (incl. dip/strike decomposition)");

   // =================================================================
   // PART B — real ResolveRateState (uniform config) wiring.
   // =================================================================
   spatial::RateStateBlock blk;     // uniform defaults (a=0.010,b=0.015,Dc=0.004,...)
   blk.eta_auto = true;
   MaterialField material =
      MaterialField::MakeConstant(bp5.lambda(), bp5.mu(), bp5.rho);
   spatial::SpatialFrictionResolver resolver;
   spatial::PorePressureSpec pp;        // zero (sigma_n_default is effective)
   Vector empty_sigma;                  // size 0 -> resolver uses sigma_n_default
   spatial::RateStatePerDOFParams rsu = resolver.ResolveRateState(
      blk, dof_coords_3d, dof_to_elem, dof_to_attr, material, pmesh, pp,
      empty_sigma);

   const real_t eta_auto_expected = 0.5 * std::sqrt(bp5.mu() * bp5.rho);
   for (int i = 0; i < n_owned; i++)
   {
      Check(std::abs(rsu.a(i) - blk.a_default) <= 1e-14, "resolver a uniform");
      Check(std::abs(rsu.Dc(i) - blk.Dc_default) <= 1e-14, "resolver Dc uniform");
      Check(std::abs(rsu.V_init(i) - blk.V_init_default) <= 1e-14,
            "resolver V_init uniform");
      Check(std::abs(rsu.sigma_n_eff(i) - blk.sigma_n_default) <= 1e-6,
            "resolver sigma_n_eff == sigma_n_default (empty total + zero P_p)");
      Check(std::abs(rsu.eta(i) - eta_auto_expected) <= 1e-3,
            "resolver auto-eta == 0.5*sqrt(mu*rho)");
   }

   FaultGeometry<ParMesh> geom_c(domain, bp5, &mpi, /*compute_bp5_params=*/false);
   geom_c.SetRateStatePerDOF(rsu, init_vel_dir);
   for (int i = 0; i < n_owned; i++)
   {
      Check(std::abs(geom_c.GetAValues()(i) - blk.a_default) <= 1e-14,
            "geom_c a == a_default (resolver->SetRateStatePerDOF)");
      Check(std::abs(geom_c.GetVInit()(2 * i) - 0.0) <= 1e-30,
            "geom_c V_init dip == 0");
      Check(std::abs(geom_c.GetVInit()(2 * i + 1) - blk.V_init_default) <= 1e-14,
            "geom_c V_init strike == V_init_default");
   }

   // ---- Collective verdict ----
   int gf = g_fail, global_checks = g_checks, global_owned = n_owned;
   MPI_Allreduce(MPI_IN_PLACE, &gf, 1, MPI_INT, MPI_SUM, comm);
   MPI_Allreduce(MPI_IN_PLACE, &global_checks, 1, MPI_INT, MPI_SUM, comm);
   MPI_Allreduce(MPI_IN_PLACE, &global_owned, 1, MPI_INT, MPI_SUM, comm);
   global_fail = gf;

   if (g_rank == 0)
   {
      std::cout << "=== test_spatial_seas_faultgeom_parity ===\n"
                << "  global owned fault DOFs: " << global_owned << "\n"
                << "  checks (all ranks): " << global_checks
                << ", failures: " << global_fail << "\n"
                << (global_fail == 0 ? "PASS\n" : "FAIL\n");
   }
   }  // end scope — all MPI-owning objects destruct here, before MPI_Finalize

   MPI_Finalize();
   return (global_fail == 0) ? 0 : 1;
}
