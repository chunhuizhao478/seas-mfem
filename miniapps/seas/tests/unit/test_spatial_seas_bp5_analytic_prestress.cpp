// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_spatial_seas_bp5_analytic_prestress.cpp — Phase 3c of
// document/spatial_seas_dev/PLAN_spatial_seas_quasidynamic_driver_2026-05-31.md
//
// Closes the tau_pre_ half of the Phase-3 parity gap (R-007).  The plan's
// Phase-3 parity test deliberately does NOT compare tau_pre_ / sigma_n_per_dof_
// element-wise, because a UNIFORM FaultLocalPrestress cannot match BP5's
// a(x2,x3)-heterogeneous, nucleation-carrying analytic tau0 (R-007).  This test
// supplies the missing per-DOF source: spatial::Bp5AnalyticStressSource fed the
// SAME BP5Params that ComputeBP5Params uses, run through the generic
// FaultGeometry::ComputeParams<StressSource> Cauchy projection, must reproduce
// the analytic ComputeBP5Params tau_pre_ to round-off, and produce a per-DOF
// sigma_n_per_dof_ equal to the BP5 effective sigma_n.
//
// Construction (per fault DOF, owned view):
//   path (a) — analytic: FaultGeometry(domain, bp5, &mpi) -> ComputeBP5Params
//              writes tau_pre_(2i)=tau0_vec[dip], tau_pre_(2i+1)=tau0_vec[strike].
//   path (b) — projection: build Bp5AnalyticStressSource from the geometry's
//              CONSTANT (planar BP5) fault basis (n,t1,t2) taken from
//              fault_dof_basis() col 0 (constness asserted), then
//              geom.ComputeParams(source) projects S onto (n,t1,t2):
//              n.S.n=sigma_n, t1.S.n=tau_dip, t2.S.n=tau_strike.
//   Assert tau_pre_(b) == tau_pre_(a) to 1e-10 (relative), sigma_n_per_dof_(b)
//   == bp5.sigma_n, HasParams()==true, and (frame check) |dip| << |strike|
//   (BP5 is pure strike-slip; a dip/strike basis swap would fail this).
//
// Runs at np=1 and np=4.  Each rank checks its OWNED fault DOFs; ranks with no
// owned fault DOFs skip the per-DOF block and contribute only to the collective
// MPI_Allreduce(SUM) verdict.

#include "mfem.hpp"

#include "../../domain/elasticity_operator.hpp"
#include "../../fault/fault_geometry.hpp"
#include "../../fault/fault_geometry_safs_templated.inl"  // ComputeParams<StressSource>
#include "../../config/bp5_params.hpp"
#include "../../common/mpi_context.hpp"
#include "../../spatial/code/spatial_stress.hpp"           // Bp5AnalyticStressSource

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

real_t MaxAbs(const Vector &v)
{
   real_t m = 0.0;
   for (int i = 0; i < v.Size(); i++) { m = std::max(m, std::abs(v(i))); }
   return m;
}

real_t MaxAbsDiff(const Vector &a, const Vector &b)
{
   if (a.Size() != b.Size()) { return 1e300; }
   real_t m = 0.0;
   for (int i = 0; i < a.Size(); i++) { m = std::max(m, std::abs(a(i) - b(i))); }
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
#  error "test_spatial_seas_bp5_analytic_prestress requires MFEM_USE_MPI=YES."
#endif

   const char *mesh_path = "bp5/mesh/reference/bp5_tandem_exact.msh";
   const real_t rel_tol  = 1e-10;

   int global_fail = 0;
   // Scope so every MFEM/MPI-owning object destructs BEFORE MPI_Finalize.
   {
   // ---- Mesh + elasticity operator (BP5 moduli + bc), same fixture as the
   //      Phase-3 parity test. ----
   BP5Params bp5;
   Mesh smesh(mesh_path, 1, 1);
   MFEM_VERIFY(smesh.Dimension() == 3, "bp5_analytic test: need a 3D BP5 mesh");
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

   // path (a) — analytic reference tau_pre_ (ComputeBP5Params).
   FaultGeometry<ParMesh> geom_ref(domain, bp5, &mpi);
   // path (b) — projection target (same analytic ctor; ComputeParams overwrites
   //            tau_pre_/sigma_n_per_dof_ but keeps a/eta/Dc/V_init).
   FaultGeometry<ParMesh> geom_proj(domain, bp5, &mpi);

   Check(geom_ref.NumFaultDOFs() == n_owned, "geom_ref NumFaultDOFs == owned");
   Check(geom_proj.NumFaultDOFs() == n_owned, "geom_proj NumFaultDOFs == owned");

   if (n_owned > 0)
   {
      // Snapshot the analytic reference tau_pre_ before the projection overwrites
      // the SECOND geometry (the two share the same domain -> same owned ordering,
      // coords, and basis).
      Vector tau_ref(geom_ref.GetTauPre());   // deep copy (size 2N)
      Check(tau_ref.Size() == 2 * n_owned, "tau_ref size == 2*owned");

      // ---- Extract the CONSTANT planar-BP5 fault basis from the geometry. ----
      const DenseMatrix &B = geom_ref.fault_dof_basis();
      Check(B.Height() == 9 && B.Width() == n_owned,
            "fault_dof_basis shape == (9, owned)");
      const real_t col0[9] = { B(0, 0), B(1, 0), B(2, 0),
                               B(3, 0), B(4, 0), B(5, 0),
                               B(6, 0), B(7, 0), B(8, 0) };
      real_t n[3]  = { col0[0], col0[1], col0[2] };
      real_t t1[3] = { col0[3], col0[4], col0[5] };
      real_t t2[3] = { col0[6], col0[7], col0[8] };

      // The BP5 fault is planar, so every owned DOF's (n,t1,t2) triple equals
      // col 0 UP TO A GLOBAL SIGN FLIP of the whole triple: the CalcOrtho /
      // FaultBasis orientation can flip the normal (and, to stay right-handed,
      // t1,t2) on some faces/partitions (fault_geometry.hpp:1452-1467).  A full
      // ±flip is projection-INVARIANT (t1.S.n == (-t1).S.(-n)), so the single
      // col-0 basis Bp5AnalyticStressSource embeds still recovers each DOF's
      // tau0_vec exactly (verified by the tau_pre_ parity below).  A PARTIAL
      // flip or in-plane rotation would break BOTH this ±check and the parity,
      // so this guards the single-basis assumption without false positives.
      real_t basis_dev = 0.0;
      for (int i = 0; i < n_owned; i++)
      {
         real_t dev_plus = 0.0, dev_minus = 0.0;
         dev_plus  = std::max(dev_plus,  std::abs(B(0, i) - col0[0]));
         dev_plus  = std::max(dev_plus,  std::abs(B(1, i) - col0[1]));
         dev_plus  = std::max(dev_plus,  std::abs(B(2, i) - col0[2]));
         dev_plus  = std::max(dev_plus,  std::abs(B(3, i) - col0[3]));
         dev_plus  = std::max(dev_plus,  std::abs(B(4, i) - col0[4]));
         dev_plus  = std::max(dev_plus,  std::abs(B(5, i) - col0[5]));
         dev_plus  = std::max(dev_plus,  std::abs(B(6, i) - col0[6]));
         dev_plus  = std::max(dev_plus,  std::abs(B(7, i) - col0[7]));
         dev_plus  = std::max(dev_plus,  std::abs(B(8, i) - col0[8]));
         dev_minus = std::max(dev_minus, std::abs(B(0, i) + col0[0]));
         dev_minus = std::max(dev_minus, std::abs(B(1, i) + col0[1]));
         dev_minus = std::max(dev_minus, std::abs(B(2, i) + col0[2]));
         dev_minus = std::max(dev_minus, std::abs(B(3, i) + col0[3]));
         dev_minus = std::max(dev_minus, std::abs(B(4, i) + col0[4]));
         dev_minus = std::max(dev_minus, std::abs(B(5, i) + col0[5]));
         dev_minus = std::max(dev_minus, std::abs(B(6, i) + col0[6]));
         dev_minus = std::max(dev_minus, std::abs(B(7, i) + col0[7]));
         dev_minus = std::max(dev_minus, std::abs(B(8, i) + col0[8]));
         basis_dev = std::max(basis_dev, std::min(dev_plus, dev_minus));
      }
      Check(basis_dev <= 1e-12,
            "BP5 fault basis is +/-constant across owned DOFs (planar, full-flip ok)");

      // ---- path (b): Bp5AnalyticStressSource through the Cauchy projection. ----
      spatial::Bp5AnalyticStressSource source(bp5, n, t1, t2);
      geom_proj.ComputeParams(source);   // P_p = 0 (bp5.sigma_n is effective)

      Check(geom_proj.HasParams(), "geom_proj HasParams() after ComputeParams");

      const Vector &tau_proj = geom_proj.GetTauPre();
      const Vector &sn_proj  = geom_proj.sigma_n_per_dof();
      Check(tau_proj.Size() == 2 * n_owned, "projected tau_pre_ size == 2*owned");
      Check(sn_proj.Size() == n_owned, "projected sigma_n_per_dof_ size == owned");

      // (1) tau_pre_ parity to round-off (relative to the reference magnitude).
      const real_t tau_scale = MaxAbs(tau_ref) + 1.0;
      Check(MaxAbsDiff(tau_proj, tau_ref) <= rel_tol * tau_scale,
            "tau_pre_ projection == analytic ComputeBP5Params (R-007)");

      // (2) per-DOF effective sigma_n == BP5 sigma_n.
      real_t sn_err = 0.0;
      for (int i = 0; i < n_owned; i++)
      {
         sn_err = std::max(sn_err, std::abs(sn_proj(i) - bp5.sigma_n));
      }
      Check(sn_err <= rel_tol * bp5.sigma_n,
            "sigma_n_per_dof_ == bp5.sigma_n (effective, compression +)");

      // (3) frame check: BP5 is pure strike-slip, so the dip slot (2i) must be
      //     negligible vs the strike slot (2i+1).  A dip/strike basis swap or a
      //     wrong (x2,x3) mapping would surface here.
      real_t max_dip = 0.0, max_strike = 0.0;
      for (int i = 0; i < n_owned; i++)
      {
         max_dip    = std::max(max_dip,    std::abs(tau_proj(2 * i)));
         max_strike = std::max(max_strike, std::abs(tau_proj(2 * i + 1)));
      }
      Check(max_strike > 0.0, "strike prestress is non-zero");
      Check(max_dip <= 1e-6 * (max_strike + 1.0),
            "dip prestress negligible vs strike (pure strike-slip frame)");
   }

   // ---- Collective verdict ----
   int gf = g_fail, global_checks = g_checks, global_owned = n_owned;
   MPI_Allreduce(MPI_IN_PLACE, &gf, 1, MPI_INT, MPI_SUM, comm);
   MPI_Allreduce(MPI_IN_PLACE, &global_checks, 1, MPI_INT, MPI_SUM, comm);
   MPI_Allreduce(MPI_IN_PLACE, &global_owned, 1, MPI_INT, MPI_SUM, comm);
   global_fail = gf;

   if (g_rank == 0)
   {
      std::cout << "=== test_spatial_seas_bp5_analytic_prestress ===\n"
                << "  global owned fault DOFs: " << global_owned << "\n"
                << "  checks (all ranks): " << global_checks
                << ", failures: " << global_fail << "\n"
                << (global_fail == 0 ? "PASS\n" : "FAIL\n");
   }
   }  // end scope — MPI-owning objects destruct before MPI_Finalize

   MPI_Finalize();
   return (global_fail == 0) ? 0 : 1;
}
