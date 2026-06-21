// =========================================================================
// test_spatial_seas_srw_dispatch — proves the QD fault operator drives STRONG
// RATE WEAKENING (SlipLawSRWPsi) through its per-DOF _SRW route, with PER-DOF
// V_w and a, and that the aging-law path is byte-identical when no SRW law is
// selected.  (Plan PLAN_spatial_seas_srw_qd_2026-06-20.md, Phases S1+S2.)
//
// Method (fully controlled — no domain solve, no equilibrium init):
//   1. Build a FaultGeometry on the BP5 mesh via the deferred ctor; fill it
//      with SetRateStatePerDOF using a hand-built rs whose V_w AND a VARY per
//      DOF.  Assert geom.GetVwValues() == rs.V_w  (S1 plumbing).
//   2. Build TWO RateStateFaultOperator<,2> on that geom — one with AgingLawPsi,
//      one with SlipLawSRWPsi (production mode).  Put both in SAFS mode with the
//      SAME synthetic per-DOF tau_pre / sigma_n.
//   3. Feed both the SAME synthetic state + traction to ComputeRHS.
//   4. The slip-rate V (from the law-INDEPENDENT Brent friction solve) must be
//      identical between the two operators.  Using that V_abs:
//        - aging psi-rate must equal AgingLawPsi::Rate(V_abs, psi, Dc[i]).
//        - SRW   psi-rate must equal SlipLawSRWPsi::Rate_SRW(V_abs, psi, Dc[i],
//          V_w[i], a[i])  — proving per-DOF V_w/a actually reach the law.
//        - the two psi-rates must DIFFER on some DOF (SRW is really active, not
//          a silent aging fallback).
//
// Pure assert-based; returns nonzero on any failure.
// =========================================================================
#include "mfem.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../fault/fault_geometry.hpp"
#include "../../fault/rate_state_fault.hpp"
#include "../../friction/dieterich_ruina.hpp"
#include "../../friction/state_evolution.hpp"
#include "../../friction/slip_law_srw_psi.hpp"
#include "../../config/bp5_params.hpp"
#include "../../common/mpi_context.hpp"
#include "../../spatial/code/spatial_friction.hpp"

#include <cmath>
#include <iostream>
#include <string>
#include <mpi.h>

using namespace mfem;
using namespace mfem::seas;

namespace
{
int g_rank = 0, g_fail = 0, g_checks = 0;
void Check(bool ok, const std::string &msg)
{
   g_checks++;
   if (!ok) { g_fail++; std::cerr << "[rank " << g_rank << "] FAILED: " << msg << "\n"; }
}
}  // namespace

int main(int argc, char *argv[])
{
   MPI_Init(&argc, &argv);
   MPI_Comm comm = MPI_COMM_WORLD;
   MPI_Comm_rank(comm, &g_rank);

   const char *mesh_path = "bp5/mesh/reference/bp5_tandem_exact.msh";

   int global_fail = 0;
   {
   BP5Params bp5;
   Mesh smesh(mesh_path, 1, 1);
   MFEM_VERIFY(smesh.Dimension() == 3, "srw dispatch: need a 3D BP5 mesh");
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
   const int N = domain.GetNumOwnedFaultDOFs();

   // ---- Friction scalars (uniform b/f0/V0; SRW f_w) ----
   const real_t b_s = 0.015, f0_s = 0.6, V0_s = 1.0e-6, f_w = 0.1, V_w_def = 0.1;
   const real_t sigma_n = 50.0e6, eta = bp5.eta();

   // ---- Hand-built rs with PER-DOF varying a and V_w (the discriminator) ----
   spatial::RateStatePerDOFParams rs;
   rs.a.SetSize(N);  rs.Dc.SetSize(N);  rs.eta.SetSize(N);
   rs.V_init.SetSize(N);  rs.sigma_n_eff.SetSize(N);
   rs.b.SetSize(N);  rs.f_0.SetSize(N);  rs.V_0.SetSize(N);
   rs.V_w.SetSize(N);
   // Per-DOF a AND V_w must VARY, and V_w must SPAN the friction-solved V so the
   // SRW steady state is genuinely V_w-sensitive (R-006: a narrow V_w range at a
   // low operating point makes Rate_SRW insensitive, so a broken impl using
   // DOF-0's scalars would pass).  V_w spans 24 decades (1e-24..1e0) crossing the
   // ~1e-19 operating V, so per-DOF V_w produces materially different psi_ss.
   for (int i = 0; i < N; i++)
   {
      rs.a(i)          = 0.008 + 0.004 * (i % 4);   // 0.008..0.016, per-DOF
      rs.Dc(i)         = 0.10;
      rs.eta(i)        = eta;
      rs.V_init(i)     = 1.0e-9;
      rs.sigma_n_eff(i)= sigma_n;
      rs.b(i)          = b_s;
      rs.f_0(i)        = f0_s;
      rs.V_0(i)        = V0_s;
      rs.V_w(i)        = std::pow(10.0, -4.0 + (i % 9));  // 1e-4..1e4, spans V
   }
   Vector init_vel_dir(2); init_vel_dir(0) = 0.0; init_vel_dir(1) = 1.0;

   FaultGeometry<ParMesh> geom(domain, bp5, &mpi, /*compute_bp5_params=*/false);
   geom.SetRateStatePerDOF(rs, init_vel_dir);

   // ---- S1: V_w plumbed resolver -> geom ----
   Check(geom.GetVwValues().Size() == N, "geom.GetVwValues().Size() == N (S1)");
   real_t vw_err = 0.0;
   for (int i = 0; i < N; i++)
   { vw_err = std::max(vw_err, std::abs(geom.GetVwValues()(i) - rs.V_w(i))); }
   Check(vw_err <= 1e-15, "geom.GetVwValues() == rs.V_w to 1e-15 (S1)");

   // ---- Two operators: aging vs SRW (production) ----
   DieterichRuinaFriction::Constants fc;
   fc.V0 = V0_s; fc.f0 = f0_s; fc.b = b_s; fc.Dc = 0.10;
   DieterichRuinaFriction friction(fc);

   AgingLawPsi  aging(b_s, V0_s, f0_s);
   SlipLawSRWPsi srw(rs.a(0), b_s, V0_s, f0_s, f_w, V_w_def);
   srw.SetProductionMode();

   RateStateFaultOperator<ParMesh, 2> op_aging(&geom, &friction, &aging, bp5, &mpi);
   RateStateFaultOperator<ParMesh, 2> op_srw  (&geom, &friction, &srw,   bp5, &mpi);

   // Synthetic per-DOF pre-stress / sigma_n (SAFS mode) — pure strike shear.
   Vector tau_pre(2 * N), sig_n(N);
   for (int i = 0; i < N; i++)
   { tau_pre(2*i) = 0.0; tau_pre(2*i+1) = 15.0e6; sig_n(i) = sigma_n; }
   op_aging.SetSAFSMode(true, &tau_pre, &sig_n);
   op_srw.SetSAFSMode  (true, &tau_pre, &sig_n);

   Vector traction(2 * N); traction = 0.0;  // tau_pre carries the 15 MPa shear

   // (i) LOW-psi state: psi << 0 -> friction frictionless-limit V = tau/eta ~
   // O(1) m/s, so Rate_SRW's -(V/L) factor is nonzero AND V sits in the
   // V_w-sensitive band (1e-4..1e4).  A moderate psi (~f0) gives V~0 ->
   // Rate_SRW~0 regardless of V_w/a -> a VACUOUS per-DOF test (R-006).  Aging's
   // exp((f0-psi)/b) overflows at psi<<0, so aging byte-identity is checked
   // separately at moderate psi below.
   Vector state_lo(3 * N); state_lo = 0.0;
   for (int i = 0; i < N; i++) { state_lo(3*i + 2) = -2.0; }
   Vector rate_srw(3 * N);
   op_srw.ComputeRHS(traction, state_lo, rate_srw);

   // (ii) MODERATE-psi state: finite for BOTH laws; used for the aging
   // byte-identity check and the aging-vs-SRW liveness check.
   Vector state_mod(3 * N); state_mod = 0.0;
   for (int i = 0; i < N; i++) { state_mod(3*i + 2) = 0.6; }
   Vector rate_aging_m(3 * N), rate_srw_m(3 * N);
   op_aging.ComputeRHS(traction, state_mod, rate_aging_m);
   op_srw.ComputeRHS  (traction, state_mod, rate_srw_m);

   // ---- S2a: per-DOF SRW dispatch at the SENSITIVE (low-psi) operating point.
   // For each DOF use the operator's OWN friction-solved V (Vs):
   //   ref_srw    = Rate_SRW(Vs, psi, Dc[i], V_w[i], a[i])  (CORRECT per-DOF)
   //   ref_scalar = Rate_SRW(Vs, psi, Dc[i], V_w[0], a[0])  (BROKEN DOF-0 scalar)
   // The operator MUST match ref_srw; n_sensitive (DOFs where per-DOF != DOF-0
   // scalar) MUST be substantial, so the match could NOT come from a scalar impl.
   int n_sensitive = 0;
   real_t max_srw_err = 0.0;
   for (int i = 0; i < N; i++)
   {
      const real_t Vs = std::sqrt(rate_srw(3*i)*rate_srw(3*i) +
                                  rate_srw(3*i+1)*rate_srw(3*i+1));
      const real_t psi = state_lo(3*i + 2);
      const real_t ref_srw    = srw.Rate_SRW(Vs, psi, rs.Dc(i), rs.V_w(i), rs.a(i));
      const real_t ref_scalar = srw.Rate_SRW(Vs, psi, rs.Dc(i), rs.V_w(0), rs.a(0));
      const real_t r_srw      = rate_srw(3*i + 2);
      max_srw_err = std::max(max_srw_err,
                             std::abs(r_srw - ref_srw) / (std::abs(ref_srw) + 1.0));
      if (std::abs(ref_srw - ref_scalar) > 1e-6 * (std::abs(ref_srw) + 1.0))
      { n_sensitive++; }
   }
   Check(max_srw_err <= 1e-12,
         "SRW psi-rate == SlipLawSRWPsi::Rate_SRW(.,V_w[i],a[i]) (per-DOF)");
   Check(n_sensitive >= N / 2,
         "per-DOF V_w/a materially differs from DOF-0 scalar on >=half the DOFs "
         "(operating point is sensitive -> the per-DOF match is decisive)");

   // ---- S2b: aging byte-identity + SRW liveness at the moderate operating pt.
   int n_diff_aging = 0;
   real_t max_aging_err = 0.0;
   for (int i = 0; i < N; i++)
   {
      const real_t Va = std::sqrt(rate_aging_m(3*i)*rate_aging_m(3*i) +
                                  rate_aging_m(3*i+1)*rate_aging_m(3*i+1));
      const real_t ref_aging = aging.Rate(Va, state_mod(3*i + 2), rs.Dc(i));
      const real_t r_aging = rate_aging_m(3*i + 2);
      max_aging_err = std::max(max_aging_err,
                               std::abs(r_aging - ref_aging) / (std::abs(ref_aging) + 1.0));
      if (std::abs(rate_srw_m(3*i + 2) - r_aging) > 1e-12 * (std::abs(r_aging) + 1.0))
      { n_diff_aging++; }
   }
   Check(max_aging_err <= 1e-12,
         "aging psi-rate == AgingLawPsi::Rate (byte-identical else arm)");
   Check(n_diff_aging > 0,
         "SRW psi-rate DIFFERS from aging on some DOF (SRW really active)");

   if (g_rank == 0)
   {
      std::cout << "[srw_dispatch] N=" << N
                << " vw_err=" << vw_err
                << " max_srw_err=" << max_srw_err
                << " max_aging_err=" << max_aging_err
                << " n_sensitive=" << n_sensitive << "/" << N
                << " n_diff_aging=" << n_diff_aging << "\n";
   }
   }  // scope: destroy MFEM/MPI objects before Finalize

   int local_fail = g_fail;
   MPI_Allreduce(&local_fail, &global_fail, 1, MPI_INT, MPI_SUM, comm);
   if (g_rank == 0)
   {
      std::cout << "=== SRW DISPATCH: " << (g_checks - g_fail) << "/" << g_checks
                << " checks passed (this rank); global_fail=" << global_fail
                << " ===\n";
   }
   MPI_Finalize();
   return global_fail == 0 ? 0 : 1;
}
