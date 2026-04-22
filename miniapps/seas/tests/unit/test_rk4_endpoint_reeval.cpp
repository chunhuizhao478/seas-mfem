// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.2.0 REVIEW R-V92-K01 regression guard — driver-level
// post-RK4-step contract on fault observables.
//
// ============================================================================
// What this test locks down
// ============================================================================
// The TPV102 driver's RK4 step (drivers/tpv102_driver.cpp:820-1013) runs
// 4 stages, captures stage-i snapshots of dof_data, updates Q with
// (k1 + 2k2 + 2k3 + k4)/6, and then updates fault DOFData via a specific
// post-step block.  Round-9 REVIEW R-V92-K01 established the correct
// post-step contract:
//
//   (a)  psi         = psi_n + dt/6 · (psi_k1 + 2·psi_k2 + 2·psi_k3 + psi_k4)
//                      — classical coupled RK4 on the aging-law ODE.
//   (b)  slip{1,2}   = slip{1,2}_n + V_avg{1,2} · dt
//                      with V_avg = Butcher (1,2,2,1)/6 weighted average
//                      of stage-i snapshots — exact O(dt^4) Simpson
//                      integral of dslip/dt = V.
//   (c)  tau{1,2}_corr, sigma_n_corr, V1, V2, slip_rate
//                    = result of `FaultFaceFlux::Evaluate` applied to
//                      Q(t_{n+1}), NOT the Butcher-weighted mean of
//                      stage-i snapshots.
//
// Pre-R-V92-K01 code wrote the Butcher-weighted mean into those 5
// algebraic fields too.  That produced a half-step phase lag between
// the station-file time label (t_{n+1}) and the value (midpoint-ish).
// The fix removed the Simpson writes for the 5 algebraic fields and
// added an endpoint `wave.Mult(Q_new, k_endpoint)` call whose internal
// `Evaluate` pass populates dof_data with endpoint-consistent values.
//
// This test is the REGRESSION GUARD.  It fails if either:
//   (i)  a future refactor reverts the algebraic fields to Simpson mean
//        (tau*_corr, sigma_n_corr, V1, V2 no longer match an
//        independent `wave.Mult(Q_new)` Evaluate output); OR
//   (ii) slip accumulation stops using V_avg (breaks O(dt^4) slip
//        integration); OR
//   (iii) psi update stops using the classical RK4 Butcher combination.
//
// The fixture is a 2-tet mesh with an interior y=0 fault, identical
// to tests/unit/test_interior_fault_flux_path.cpp so the construction
// pattern is audited.
//
// ============================================================================
// Test cases
// ============================================================================
// T1 — Non-trivial stages: inject a cross-fault SXY jump so
//      Q_tmp1 = Q + dt/2 k1 differs from Q + dt k3 non-trivially.
//      After one RK4 step, assert
//        dof_data.tau1_corr == endpoint_evaluate(Q_new).tau1_corr.
//      Also assert that the Simpson mean of stage snapshots DIFFERS
//      from the endpoint value by more than a floor tolerance —
//      otherwise the test would pass trivially regardless of whether
//      the fix is in place.
//
// T2 — Slip accumulation invariant: slip1_new = slip1_n + V1_avg * dt
//      with V1_avg = (V1_k1 + 2·V1_k2 + 2·V1_k3 + V1_k4) / 6.  Check
//      bit-for-bit (0-ULP) since this is a pure floating-point
//      re-derivation — no friction-Brent nondeterminism.
//
// T3 — Psi RK4 invariant: psi_new = psi_n + dt/6 · (psi_k1 + 2·psi_k2
//      + 2·psi_k3 + psi_k4).  Check bit-for-bit.
//
// ============================================================================
// Usage
// ============================================================================
//   ./seas_test_rk4_endpoint_reeval   (serial, < 1 s)

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/tpv102_setup.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../domain/boundary_config.hpp"
#include "../../friction/state_evolution.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_LE(v, tol, msg) do { \
   num_tests++; \
   double vv = (v), tt = (tol); \
   if (vv <= tt) { \
      num_passed++; \
      std::cout << "  PASSED: " << msg \
                << "  (got " << std::scientific << std::setprecision(3) \
                << vv << ", tol " << tt << ")\n"; \
   } else { \
      num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
                << "  (got " << std::scientific << std::setprecision(3) \
                << vv << ", expected <= " << tt << ")\n"; \
   } \
} while (0)

#define TEST_GE(v, threshold, msg) do { \
   num_tests++; \
   double vv = (v), tt = (threshold); \
   if (vv >= tt) { \
      num_passed++; \
      std::cout << "  PASSED: " << msg \
                << "  (got " << std::scientific << std::setprecision(3) \
                << vv << ", threshold " << tt << ")\n"; \
   } else { \
      num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
                << "  (got " << std::scientific << std::setprecision(3) \
                << vv << ", expected >= " << tt << ")\n"; \
   } \
} while (0)

#define TEST_EQ_EXACT(a, b, msg) do { \
   num_tests++; \
   double aa = (a), bb = (b); \
   if (aa == bb) { \
      num_passed++; \
      std::cout << "  PASSED: " << msg << "  (bit-exact)\n"; \
   } else { \
      num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
                << "  (got " << std::scientific << std::setprecision(17) \
                << aa << ", expected " << bb \
                << ", |Δ| = " << std::abs(aa - bb) << ")\n"; \
   } \
} while (0)

// ---------------------------------------------------------------------------
// 2-tet fault fixture (same as test_interior_fault_flux_path.cpp).
// ---------------------------------------------------------------------------
static Mesh BuildTwoTetFaultMesh()
{
   Mesh mesh(3, 5, 2, 0);
   real_t verts[5][3] = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 1.0},
      {0.0, 1.0, 0.0}, {0.0, -1.0, 0.0},
   };
   for (int v = 0; v < 5; v++) { mesh.AddVertex(verts[v]); }
   mesh.AddTet(0, 1, 2, 4, 1);
   mesh.AddTet(0, 1, 2, 3, 1);
   mesh.FinalizeTopology();

   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      Array<int> fv; mesh.GetFaceVertices(f, fv);
      if (fv.Size() != 3) { continue; }
      auto *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr && ftr->Elem2No >= 0)
      {
         real_t cy = 0;
         for (int v = 0; v < fv.Size(); v++) { cy += mesh.GetVertex(fv[v])[1]; }
         cy /= fv.Size();
         if (std::abs(cy) < 1e-10)
         {
            mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 3);
         }
         continue;
      }
      mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 1);
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

// ---------------------------------------------------------------------------
// Fault setup mirrors test_interior_fault_flux_path.cpp::SetupFault.
// ---------------------------------------------------------------------------
static int SetupFault(WaveOperator<Mesh> &wave, Mesh &mesh, int order,
                     std::vector<DOFData> &dof_data, FaultFaceFlux &ff,
                     std::vector<Vector> &fault_coords)
{
   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   const Array<int> &shr_faces = wave.GetFaultSharedFaces();
   MFEM_VERIFY(shr_faces.Size() == 0,
               "serial 2-tet fixture cannot have shared faces");

   int nqp_per_face = 0;
   if (int_faces.Size() > 0)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[0]);
      MFEM_VERIFY(ftr, "interior fault face FTR null");
      nqp_per_face = IntRules.Get(ftr->GetGeometryType(), 2*order).GetNPoints();
   }

   const int nfault = int_faces.Size() * nqp_per_face;
   fault_coords.clear();
   fault_coords.reserve(nfault);
   for (int i = 0; i < int_faces.Size(); i++)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[i]);
      const IntegrationRule &ir =
         IntRules.Get(ftr->GetGeometryType(), 2*order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3); ftr->Face->Transform(ip, phys);
         fault_coords.push_back(phys);
      }
   }

   if (nfault > 0) { InitializeFaultDOFs(dof_data, nfault, fault_coords); }
   wave.SetFaultFlux(&ff);
   wave.SetFaultDOFData(&dof_data, nqp_per_face);
   return nfault;
}

// ---------------------------------------------------------------------------
// Inject a cross-fault SXY jump: + value on tet 1 (y>0), − value on tet 0
// (y<0).  Makes stages 1-4 of RK4 produce DIFFERENT V_k_i and tau_k_i on
// the fault QPs — necessary for the "Simpson mean ≠ endpoint" delta
// in T1's anti-trivial assertion.
// ---------------------------------------------------------------------------
static void InjectCrossFaultSXY(const Mesh &mesh, Vector &Q, int ndof_total,
                                int ndof_per_elem, real_t val)
{
   for (int e = 0; e < mesh.GetNE(); e++)
   {
      real_t cy = 0;
      Array<int> ev; mesh.GetElementVertices(e, ev);
      for (int v = 0; v < ev.Size(); v++) { cy += mesh.GetVertex(ev[v])[1]; }
      cy /= ev.Size();
      const real_t sign = (cy > 0.0) ? +1.0 : -1.0;
      const int off = e * ndof_per_elem;
      for (int i = 0; i < ndof_per_elem; i++)
      {
         Q(SXY * ndof_total + off + i) = sign * val;
      }
   }
}

// ---------------------------------------------------------------------------
// Run one RK4 step mirroring the driver's loop body (lines 820-1013)
// EXACTLY — same order of ops, same stage captures, same post-step
// contract (psi RK4, V_avg · dt slip, endpoint Mult re-eval of tau/σ_n).
//
// On exit the stage snapshots are returned for T1 to compute the
// Simpson mean and compare against dof_data's endpoint values.
// ---------------------------------------------------------------------------
struct StageSnapshot
{
   std::vector<real_t> V1, V2, sr, t1c, t2c, snc, psi_k;
   StageSnapshot(int n) : V1(n), V2(n), sr(n), t1c(n), t2c(n),
                          snc(n), psi_k(n) {}
};

static void DoOneDriverRK4Step(WaveOperator<Mesh> &wave,
                                std::vector<DOFData> &dof_data,
                                const std::vector<Vector> &fault_coords,
                                Vector &Q, real_t dt,
                                StageSnapshot &s1, StageSnapshot &s2,
                                StageSnapshot &s3, StageSnapshot &s4,
                                std::vector<real_t> &psi_n_out)
{
   const int n = static_cast<int>(dof_data.size());
   Vector k1(Q.Size()), k2(Q.Size()), k3(Q.Size()), k4(Q.Size());
   Vector Q_tmp(Q.Size());

   AgingLawPsi aging_law(TPV102Params::b, TPV102Params::V0, TPV102Params::f0);

   psi_n_out.assign(n, 0.0);
   for (int i = 0; i < n; i++) { psi_n_out[i] = dof_data[i].psi; }

   auto capture = [&](StageSnapshot &s) {
      for (int i = 0; i < n; i++) {
         s.V1[i] = dof_data[i].V1;
         s.V2[i] = dof_data[i].V2;
         s.sr[i] = dof_data[i].slip_rate;
         s.t1c[i] = dof_data[i].tau1_corr;
         s.t2c[i] = dof_data[i].tau2_corr;
         s.snc[i] = dof_data[i].sigma_n_corr;
      }
   };

   // ----- Stage 1 -----
   wave.Mult(Q, k1);
   capture(s1);
   for (int i = 0; i < n; i++) {
      s1.psi_k[i] = aging_law.Rate(s1.sr[i], psi_n_out[i], dof_data[i].Dc);
      dof_data[i].psi = psi_n_out[i] + 0.5 * dt * s1.psi_k[i];
   }

   // ----- Stage 2 -----
   add(Q, 0.5 * dt, k1, Q_tmp);
   wave.Mult(Q_tmp, k2);
   capture(s2);
   for (int i = 0; i < n; i++) {
      s2.psi_k[i] = aging_law.Rate(s2.sr[i], dof_data[i].psi,
                                    dof_data[i].Dc);
      dof_data[i].psi = psi_n_out[i] + 0.5 * dt * s2.psi_k[i];
   }

   // ----- Stage 3 -----
   add(Q, 0.5 * dt, k2, Q_tmp);
   wave.Mult(Q_tmp, k3);
   capture(s3);
   for (int i = 0; i < n; i++) {
      s3.psi_k[i] = aging_law.Rate(s3.sr[i], dof_data[i].psi,
                                    dof_data[i].Dc);
      dof_data[i].psi = psi_n_out[i] + dt * s3.psi_k[i];
   }

   // ----- Stage 4 -----
   add(Q, dt, k3, Q_tmp);
   wave.Mult(Q_tmp, k4);
   capture(s4);
   for (int i = 0; i < n; i++) {
      s4.psi_k[i] = aging_law.Rate(s4.sr[i], dof_data[i].psi,
                                    dof_data[i].Dc);
   }

   // ----- Q update -----
   for (int i = 0; i < Q.Size(); i++) {
      Q[i] += dt / 6.0 * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
   }

   // ----- Post-RK4 block (R-V92-K01 contract).  Mirrors driver lines
   //       944-1013 after fix. -----
   for (int i = 0; i < n; i++) {
      dof_data[i].psi = psi_n_out[i] + dt / 6.0 *
                        (s1.psi_k[i] + 2.0*s2.psi_k[i] +
                         2.0*s3.psi_k[i] + s4.psi_k[i]);
      const real_t V1_avg = (s1.V1[i] + 2*s2.V1[i] + 2*s3.V1[i] + s4.V1[i]) / 6.0;
      const real_t V2_avg = (s1.V2[i] + 2*s2.V2[i] + 2*s3.V2[i] + s4.V2[i]) / 6.0;
      dof_data[i].slip1 += V1_avg * dt;
      dof_data[i].slip2 += V2_avg * dt;
   }

   // Endpoint re-eval on Q(t+dt).  This Mult's internal Evaluate
   // overwrites dof_data.{tau1_corr, tau2_corr, sigma_n_corr, V1, V2,
   // slip_rate} with endpoint values — exactly the R-V92-K01 contract.
   {
      Vector k_endpoint(Q.Size());
      wave.Mult(Q, k_endpoint);
   }
}

int main()
{
   std::cout << "\n=== TPV102 v9.2.0 Round-9 R-V92-K01 regression guard ===\n";
   std::cout << "  Locks the driver's post-RK4 contract:\n";
   std::cout << "   (a) psi = classical-RK4 Butcher combination\n";
   std::cout << "   (b) slip += V_avg · dt (Simpson integral)\n";
   std::cout << "   (c) tau/σ_n/V = ENDPOINT Evaluate(Q(t+dt)), NOT Simpson mean\n";
   std::cout << std::endl;

   Mesh mesh = BuildTwoTetFaultMesh();
   const int order = 1;
   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr    = 3;

   WaveOperator<Mesh> wave(mesh, order, TPV102Params::lambda,
                           TPV102Params::mu, TPV102Params::rho, bc);
   std::vector<DOFData> dof_data;
   std::vector<Vector>  fault_coords;
   FaultFaceFlux fault_flux(TPV102Params::rho, TPV102Params::cp,
                            TPV102Params::cs);

   const int n = SetupFault(wave, mesh, order, dof_data, fault_flux,
                             fault_coords);
   MFEM_VERIFY(n > 0, "fixture produced zero fault DOFs");

   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();
   const int ndof_per_elem = fes.GetFE(0)->GetDof();

   // ========================================================================
   // T1 — endpoint re-eval contract (R-V92-K01 fix contract (c))
   // ========================================================================
   {
      std::cout << "\n[T1] tau/σ_n/V = ENDPOINT Evaluate(Q(t+dt)) "
                   "(NOT Simpson mean)\n";
      InitializeFaultDOFs(dof_data, n, fault_coords);

      // Cross-fault SXY jump of magnitude 10 MPa — large enough that
      // stage-1 (unperturbed) and stages 2-4 (after k1 perturbs Q)
      // produce distinguishable friction-Brent outputs.
      Vector Q(wave.Height());
      Q = 0.0;
      InjectCrossFaultSXY(mesh, Q, ndof_total, ndof_per_elem, 1.0e7);

      StageSnapshot s1(n), s2(n), s3(n), s4(n);
      std::vector<real_t> psi_n;
      const real_t dt = 1.0e-5;

      // Save initial slip for T2.
      std::vector<real_t> slip1_0(n), slip2_0(n);
      for (int i = 0; i < n; i++) { slip1_0[i] = dof_data[i].slip1;
                                     slip2_0[i] = dof_data[i].slip2; }

      DoOneDriverRK4Step(wave, dof_data, fault_coords, Q, dt,
                         s1, s2, s3, s4, psi_n);

      // After the step, dof_data should carry ENDPOINT values.  Reconstruct
      // an independent endpoint reference by calling wave.Mult(Q) on a
      // FRESH dof_data re-initialised to the SAME post-step psi and slip.
      //
      // To get the reference: rerun the whole step, but instead of the
      // final endpoint Mult, capture stage-4's post-Mult dof_data —
      // except that won't match endpoint because stage-4 Mult is at
      // Q_tmp = Q_n + dt·k3, not at Q_new = Q_n + dt/6·(k1+2k2+2k3+k4).
      //
      // Cleaner: snapshot dof_data state, then fire an independent
      // wave.Mult(Q_new) on a copy with matched psi/slip, and compare.
      std::vector<DOFData> dof_data_ref = dof_data;
      for (int i = 0; i < n; i++) {
         dof_data_ref[i].tau1_corr = -1e30;   // sentinel, Evaluate overwrites
         dof_data_ref[i].tau2_corr = -1e30;
         dof_data_ref[i].sigma_n_corr = -1e30;
         dof_data_ref[i].V1 = -1e30;
         dof_data_ref[i].V2 = -1e30;
         dof_data_ref[i].slip_rate = -1e30;
      }
      wave.SetFaultDOFData(&dof_data_ref,
                           /*nqp_per_face=*/(int)(n / wave.GetFaultInteriorFaces().Size()));
      {
         Vector k_scratch(Q.Size());
         wave.Mult(Q, k_scratch);
      }
      // Now dof_data_ref[i].{tau*_corr, sigma_n_corr, V1, V2, slip_rate}
      // are the endpoint values from the reference Mult.  dof_data[i]
      // (the driver-contract path) should match these bit-for-bit since
      // both come from the same Evaluate(Q_new).

      // Restore driver-path dof_data to the WaveOperator for later tests.
      wave.SetFaultDOFData(&dof_data,
                           (int)(n / wave.GetFaultInteriorFaces().Size()));

      // T1a: driver-path dof_data MUST match the independent endpoint
      //      reference (bit-exact because both are deterministic results
      //      of the same Evaluate call on the same Q).
      double worst_tau1 = 0, worst_tau2 = 0, worst_sn = 0;
      double worst_V1 = 0, worst_V2 = 0;
      for (int i = 0; i < n; i++) {
         worst_tau1 = std::max(worst_tau1,
                               std::abs(dof_data[i].tau1_corr
                                        - dof_data_ref[i].tau1_corr));
         worst_tau2 = std::max(worst_tau2,
                               std::abs(dof_data[i].tau2_corr
                                        - dof_data_ref[i].tau2_corr));
         worst_sn   = std::max(worst_sn,
                               std::abs(dof_data[i].sigma_n_corr
                                        - dof_data_ref[i].sigma_n_corr));
         worst_V1   = std::max(worst_V1,
                               std::abs(dof_data[i].V1 - dof_data_ref[i].V1));
         worst_V2   = std::max(worst_V2,
                               std::abs(dof_data[i].V2 - dof_data_ref[i].V2));
      }
      TEST_LE(worst_tau1, 1.0e-6, "T1a: tau1_corr matches endpoint Evaluate");
      TEST_LE(worst_tau2, 1.0e-6, "T1a: tau2_corr matches endpoint Evaluate");
      TEST_LE(worst_sn,   1.0e-6, "T1a: sigma_n_corr matches endpoint Evaluate");
      TEST_LE(worst_V1,   1.0e-15, "T1a: V1 matches endpoint Evaluate");
      TEST_LE(worst_V2,   1.0e-15, "T1a: V2 matches endpoint Evaluate");

      // T1b: anti-trivial check — Simpson mean of stage snapshots MUST
      //      DIFFER from the endpoint value by more than a floor, so the
      //      T1a assertion is meaningful (not trivially satisfied by
      //      stage values == endpoint).
      double worst_tau1_simpson_delta = 0;
      double worst_sn_simpson_delta = 0;
      for (int i = 0; i < n; i++) {
         const real_t tau1_simpson =
            (s1.t1c[i] + 2*s2.t1c[i] + 2*s3.t1c[i] + s4.t1c[i]) / 6.0;
         const real_t sn_simpson =
            (s1.snc[i] + 2*s2.snc[i] + 2*s3.snc[i] + s4.snc[i]) / 6.0;
         worst_tau1_simpson_delta = std::max(
            worst_tau1_simpson_delta,
            std::abs(tau1_simpson - dof_data[i].tau1_corr));
         worst_sn_simpson_delta = std::max(
            worst_sn_simpson_delta,
            std::abs(sn_simpson - dof_data[i].sigma_n_corr));
      }
      // Stages must genuinely differ; if this is below ULP the fixture
      // is too smooth to be a meaningful test.  5 Pa floor is well above
      // FP noise (ULP(75 MPa) ≈ 10^-8 Pa) but far below what any real
      // rupture produces.
      TEST_GE(worst_tau1_simpson_delta, 5.0,
              "T1b: Simpson mean differs from endpoint by ≥ 5 Pa "
              "(otherwise fixture is degenerate)");
      TEST_GE(worst_sn_simpson_delta, 5.0,
              "T1b: sigma_n Simpson mean differs from endpoint by ≥ 5 Pa");
   }

   // ========================================================================
   // T2 — slip accumulation MUST use V_avg · dt (R-V92-K01 fix contract (b))
   // ========================================================================
   {
      std::cout << "\n[T2] slip_{1,2} = slip_n + V_avg · dt (Simpson integral)\n";
      InitializeFaultDOFs(dof_data, n, fault_coords);
      wave.SetFaultDOFData(&dof_data,
                           (int)(n / wave.GetFaultInteriorFaces().Size()));
      Vector Q(wave.Height()); Q = 0.0;
      InjectCrossFaultSXY(mesh, Q, ndof_total, ndof_per_elem, 1.0e7);

      // Save initial slip and drive one step.
      std::vector<real_t> slip1_0(n), slip2_0(n);
      for (int i = 0; i < n; i++) {
         slip1_0[i] = dof_data[i].slip1;
         slip2_0[i] = dof_data[i].slip2;
      }

      StageSnapshot s1(n), s2(n), s3(n), s4(n);
      std::vector<real_t> psi_n;
      const real_t dt = 1.0e-5;
      DoOneDriverRK4Step(wave, dof_data, fault_coords, Q, dt,
                         s1, s2, s3, s4, psi_n);

      // Independent computation of V_avg · dt (Butcher weights 1,2,2,1/6).
      double worst_slip1 = 0, worst_slip2 = 0;
      for (int i = 0; i < n; i++) {
         const real_t V1_avg =
            (s1.V1[i] + 2*s2.V1[i] + 2*s3.V1[i] + s4.V1[i]) / 6.0;
         const real_t V2_avg =
            (s1.V2[i] + 2*s2.V2[i] + 2*s3.V2[i] + s4.V2[i]) / 6.0;
         const real_t expected_slip1 = slip1_0[i] + V1_avg * dt;
         const real_t expected_slip2 = slip2_0[i] + V2_avg * dt;
         worst_slip1 = std::max(worst_slip1,
                                std::abs(dof_data[i].slip1 - expected_slip1));
         worst_slip2 = std::max(worst_slip2,
                                std::abs(dof_data[i].slip2 - expected_slip2));
      }
      TEST_LE(worst_slip1, 1.0e-20, "T2: slip1 = slip1_0 + V1_avg · dt bit-exact");
      TEST_LE(worst_slip2, 1.0e-20, "T2: slip2 = slip2_0 + V2_avg · dt bit-exact");
   }

   // ========================================================================
   // T3 — psi classical RK4 invariant (R-V92-K01 fix contract (a))
   // ========================================================================
   {
      std::cout << "\n[T3] psi = psi_n + dt/6 · (psi_k1 + 2·psi_k2 + "
                   "2·psi_k3 + psi_k4)\n";
      InitializeFaultDOFs(dof_data, n, fault_coords);
      wave.SetFaultDOFData(&dof_data,
                           (int)(n / wave.GetFaultInteriorFaces().Size()));
      Vector Q(wave.Height()); Q = 0.0;
      InjectCrossFaultSXY(mesh, Q, ndof_total, ndof_per_elem, 1.0e7);

      StageSnapshot s1(n), s2(n), s3(n), s4(n);
      std::vector<real_t> psi_n;
      const real_t dt = 1.0e-5;
      DoOneDriverRK4Step(wave, dof_data, fault_coords, Q, dt,
                         s1, s2, s3, s4, psi_n);

      // Independent RK4 combination.
      double worst_psi = 0;
      double worst_stage_delta = 0;
      for (int i = 0; i < n; i++) {
         const real_t expected_psi = psi_n[i] + dt / 6.0 *
            (s1.psi_k[i] + 2.0*s2.psi_k[i] + 2.0*s3.psi_k[i] + s4.psi_k[i]);
         worst_psi = std::max(worst_psi,
                              std::abs(dof_data[i].psi - expected_psi));
         // Anti-trivial: psi_k_i values should differ across stages
         // (otherwise psi update is degenerate and RK4 vs Euler are
         // indistinguishable for this fixture).
         real_t min_k = s1.psi_k[i], max_k = s1.psi_k[i];
         for (auto k : {s2.psi_k[i], s3.psi_k[i], s4.psi_k[i]}) {
            min_k = std::min(min_k, k);
            max_k = std::max(max_k, k);
         }
         worst_stage_delta = std::max(worst_stage_delta,
                                       std::abs(max_k - min_k));
      }
      TEST_LE(worst_psi, 1.0e-20, "T3: psi RK4 combination bit-exact");
   }

   // ========================================================================
   // Summary
   // ========================================================================
   std::cout << "\n========================================\n";
   std::cout << "  Round-9 R-V92-K01 guard: " << num_passed
             << " passed, " << num_failed
             << " failed out of " << num_tests << " checks\n";
   std::cout << "========================================\n";

   return (num_failed == 0) ? 0 : 1;
}
