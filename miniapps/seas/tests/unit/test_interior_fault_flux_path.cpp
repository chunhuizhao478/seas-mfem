// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.2.0 plan Step 6 — interior-fault-path Pelties eq. (7) probe.
// Hypothesis under test: H-V92-U — interior-fault-face branch of
// WaveOperator::ComputeFaceFluxRHS (dynamic/wave_operator.inl:746-950)
// does NOT produce the Pelties 2012 eq. (7) trial traction at a fault QP
// when the bulk Q carries a simple velocity perturbation.
//
// ============================================================================
// Why this test is needed
// ============================================================================
// §17.6.5 of the plan notes that the INTERIOR-fault branch handles the
// majority of fault faces on a production TPV102 mesh, yet it had no
// dedicated unit test that drives `wave.Mult(Q, k)` with a nontrivial Q
// and then verifies the fault QP DOFData against the closed-form Pelties
// eq. (7) result.  Every previous test either:
//   - bypassed `wave.Mult` and called `FaultFaceFlux::Evaluate` directly
//     (test_fault_face_flux*.cpp — passes trivially because the input
//     Q_plus, Q_minus are set by the test, not by the wave operator);
//   - drove `wave.Mult` on Q=0 (test_no_penalty_dynamic_rupture.cpp T1);
//   - ran the shared-fault branch (test_shared_fault_dof_data_consistency).
// This file closes the coverage gap directly.
//
// ============================================================================
// Fixture
// ============================================================================
// 2-tet mesh with ONE interior face at y=0 tagged as fault (attr=3):
//
//       vertex 3 = (0, +1, 0)                            [+ side, y>0]
//       vertex 4 = (0, -1, 0)                            [- side, y<0]
//       common face vertices 0,1,2 = (0,0,0),(1,0,0),(0,0,1)  [fault @ y=0]
//
//   tet 0 = {0,1,2,4}   — centroid y<0  → − side (ref_normal = (0,-1,0))
//   tet 1 = {0,1,2,3}   — centroid y>0  → + side
//
// Boundary: all other external triangles get attr=1 (free-surface).
// Background stress lives ONLY in DOFData (sigma_n0=120 MPa, tau2_0=75 MPa);
// Q carries only perturbations (TPV102 driver convention).
//
// ============================================================================
// Analytic expectations (Pelties 2012 eq. 7, homogeneous material)
// ============================================================================
// In the canonical fault-local frame (BP5 / Tandem):
//   can_n  = (0, -1,  0)   (fault normal, + → -)
//   can_t1 = (0,  0, -1)   (dip)
//   can_t2 = (+1, 0,  0)   (strike)
// So for a stress tensor σ in global coords, the rotated components satisfy:
//   local SXX (σ_nn)       = σ_yy
//   local SXY (τ_1 = σ_nt1) = -σ_yz  (sign from can_n⊗can_t1 = -e_y⊗-e_z = +e_y⊗e_z)
//   local SXZ (τ_2 = σ_nt2) = -σ_yx = -σ_xy
// For a velocity v:
//   local VX  (v_n)  = -v_y
//   local VY  (v_t1) = -v_z
//   local VZ  (v_t2) = +v_x
//
// Test patterns (each with a clean expected tau*_trial / sigma_n_trial):
//
//   T0 — Q = 0 everywhere.  Trial = 0.  data.{tau1_corr,tau2_corr}
//        reduce to eta_s·V_ini (O(1e-9) Pa) + pre-stress.
//        Check |data.sigma_n_corr - sigma_n0| <= 1 Pa and
//              |data.tau2_corr - tau2_0| <= 1 Pa.
//
//   T1 — anti-symmetric v_y: Q[VY] = +V_test on + side, -V_test on − side.
//        Local v_n⁺ = -V_test, v_n⁻ = +V_test, so v_n⁻ - v_n⁺ = +2V_test.
//        Pelties 7a: sigma_n_trial = eta_p · 2V_test = Zp · V_test.
//        Expected: data.sigma_n_corr - sigma_n0 ≈ +Zp·V_test.
//
//   T2 — anti-symmetric v_x: Q[VX] = +V_test on + side, -V_test on − side.
//        Local v_t2⁺ = +V_test, v_t2⁻ = -V_test, so v_t2⁻ - v_t2⁺ = -2V_test.
//        Pelties 7c: tau_2_trial = eta_s · (-2V_test) = -Zs · V_test.
//        Perturbation 30 Pa on 75 MPa pre-stress — friction solver stays at
//        equilibrium (V ≈ V_ini), so V2 term in Eq.10 is O(1e-9) Pa.
//        Expected: data.tau2_corr - tau2_0 ≈ -Zs·V_test.
//
//   T3 — anti-symmetric v_z: Q[VZ] = +V_test on + side, -V_test on − side.
//        Local v_t1⁺ = -V_test, v_t1⁻ = +V_test, so v_t1⁻ - v_t1⁺ = +2V_test.
//        Pelties 7b: tau_1_trial = eta_s · 2V_test = +Zs·V_test.
//        Expected: data.tau1_corr - tau1_0 ≈ +Zs·V_test.
//
// ============================================================================
// What a FAIL tells you
// ============================================================================
// T0 FAIL  → friction solver diverged from equilibrium; check ComputeInitialPsi.
// T1 FAIL  → Pelties 7a (normal mode) broken — check the η_p path and the
//            rotation of VY ↔ local VX.  Directly confirms H-V92-U if the
//            shared-fault companion test (test_shared_fault_dof_data_consistency
//            Phase D) is green.
// T2 FAIL  → Pelties 7c (strike shear) broken — most likely culprit for the
//            observed TPV102 symptom since strike pre-stress = 75 MPa is the
//            dominant loaded channel.
// T3 FAIL  → Pelties 7b (dip shear) broken — if T2 passes but T3 fails, the
//            swap between tangent1 and tangent2 in the interior-fault rotation
//            (wave_operator.inl:822-832) is suspected.
//
// Usage:
//   ./seas_test_interior_fault_flux_path   (serial, < 1 s)

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/tpv102_setup.hpp"
// v9.3.0 Phase 4 (I-06 migration): this test originally used the
// fluctuation-Q path (Q = 0, DOFData carrying pre-stress).  Under the
// v9.3.0 migration, the TPV102 fault dispatch calls EvaluateTotal on
// total-stress Q with DOFData pre-stress zeroed.  The Pelties eq. (7)
// predictions are identical in magnitude under either representation;
// only the bulk-Q / DOFData partition changes.  tpv102_setup_total.hpp
// supplies InitializeStateTotal and ZeroDOFDataPreStressTotal.
#include "../../dynamic/tpv102_setup_total.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
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

// ---------------------------------------------------------------------------
// 2-tet mesh sharing one interior face at y=0, tagged as FAULT (attr=3).
// Other external faces are tagged as attr=1 (free surface / natural BC).
// Pattern mirrors test_shared_fault_dof_data_consistency's BuildFixtureMesh.
// ---------------------------------------------------------------------------
static Mesh BuildTwoTetFaultMesh()
{
   Mesh mesh(3, 5, 2, 0);
   real_t verts[5][3] = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 1.0},
      {0.0, 1.0, 0.0},    // vertex 3 on +y side
      {0.0, -1.0, 0.0},   // vertex 4 on -y side
   };
   for (int v = 0; v < 5; v++) { mesh.AddVertex(verts[v]); }
   mesh.AddTet(0, 1, 2, 4, 1);   // tet 0: − side
   mesh.AddTet(0, 1, 2, 3, 1);   // tet 1: + side
   mesh.FinalizeTopology();

   // External triangles (only Elem1, no Elem2) get free-surface attr=1.
   // The interior triangle (0,1,2) at y=0 gets fault attr=3 via a bdr
   // triangle overlay (same trick as test_shared_fault_dof_data_consistency.cpp
   // line 168-174).
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      Array<int> fv;
      mesh.GetFaceVertices(f, fv);
      if (fv.Size() != 3) { continue; }
      auto *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr && ftr->Elem2No >= 0)
      {
         // Interior face.  Centroid-y test: if at y=0 this is the fault.
         real_t cy = 0;
         for (int v = 0; v < fv.Size(); v++) { cy += mesh.GetVertex(fv[v])[1]; }
         cy /= fv.Size();
         if (std::abs(cy) < 1e-10)
         {
            mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 3);   // FAULT
         }
         continue;
      }
      mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 1);          // FREE
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

// ---------------------------------------------------------------------------
// Fill fault DOFData exactly the way the TPV102 driver does: background
// stress in sigma_n0/tau2_0, friction params from TPV102Params.  For this
// test we use the driver's InitializeFaultDOFs helper so the DOFData
// layout is byte-identical to production.
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

   if (nfault > 0)
   {
      InitializeFaultDOFs(dof_data, nfault, fault_coords);
   }
   wave.SetFaultFlux(&ff);
   wave.SetFaultDOFData(&dof_data, nqp_per_face);
   return nfault;
}

// Project a uniform-per-element value into a specific component of Q.
// tet 0 (centroid y<0) gets val_minus; tet 1 (centroid y>0) gets val_plus.
static void SetAntiSymQComponent(const Mesh &mesh, Vector &Q, int comp,
                                  int ndof_total, int ndof_per_elem,
                                  real_t val_plus, real_t val_minus)
{
   for (int e = 0; e < mesh.GetNE(); e++)
   {
      real_t cy = 0;
      Array<int> ev; mesh.GetElementVertices(e, ev);
      for (int v = 0; v < ev.Size(); v++) { cy += mesh.GetVertex(ev[v])[1]; }
      cy /= ev.Size();
      const real_t val = (cy > 0.0) ? val_plus : val_minus;
      const int off = e * ndof_per_elem;
      for (int i = 0; i < ndof_per_elem; i++)
      {
         Q(comp * ndof_total + off + i) = val;
      }
   }
}

// Re-initialise DOFData so the friction solver starts from equilibrium.
// Needed between test cases because `wave.Mult` advances psi via the
// friction solver and we want each T_i to start from clean state.
// Under the v9.3.0 Phase 4 total-stress migration, also zero the DOFData
// pre-stress fields (sigma_n0, tau1_0, tau2_0) because the bulk Q now
// carries those values — EvaluateTotal would double-count otherwise.
static void ResetDOFData(std::vector<DOFData> &dof_data,
                          const std::vector<Vector> &fault_coords)
{
   const int ndof = static_cast<int>(dof_data.size());
   InitializeFaultDOFs(dof_data, ndof, fault_coords);
   ZeroDOFDataPreStressTotal(dof_data, ndof);
}

// Initialize bulk Q for the total-stress migration.  Replaces `Q = 0.0`
// from the pre-migration fluctuation convention.  Callers that want to
// add a velocity-jump perturbation do so ON TOP of this initial state.
static void InitTotalStateQ(Vector &Q, int ndof_total)
{
   InitializeStateTotal(Q, ndof_total,
                        TPV102Params::sigma_n, TPV102Params::tau_ini);
}

int main()
{
   std::cout << "\n=== TPV102 v9.2.0 Step 6 — Interior-Fault-Path Pelties Eq. (7) ===\n";

   Mesh mesh = BuildTwoTetFaultMesh();

   int n_fault_bdr = 0;
   for (int b = 0; b < mesh.GetNBE(); b++)
   {
      if (mesh.GetBdrAttribute(b) == 3) { n_fault_bdr++; }
   }
   std::cout << "  Mesh       : 2 tets, " << mesh.GetNumFaces()
             << " faces, " << mesh.GetNBE() << " bdr elems, "
             << n_fault_bdr << " fault bdr\n";

   BoundaryConfig bc;
   bc.natural_attrs   = {1};
   bc.fault_attr      = 3;
   bc.absorbing_attrs = {};

   const int order = 1;
   WaveOperator<Mesh> wave(mesh, order,
                           TPV102Params::lambda,
                           TPV102Params::mu,
                           TPV102Params::rho, bc);
   const int size = wave.Height();
   const auto &fes = wave.GetFESpace();
   const int ndof_total    = fes.GetNDofs();
   const int ndof_per_elem = fes.GetFE(0)->GetDof();

   FaultFaceFlux ff(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   std::vector<DOFData> dof_data;
   std::vector<Vector> fault_coords;
   const int nfault = SetupFault(wave, mesh, order, dof_data, ff, fault_coords);
   std::cout << "  Fault QPs  : " << nfault << "\n";

   // Post-R-001: the RK4 fault dispatch now selects EvaluateTotal vs
   // Evaluate via has_bulk_bg_.  This test drives a total-Q bulk state
   // (InitTotalStateQ + ZeroDOFDataPreStressTotal), so wire the
   // background so the dispatch keeps calling EvaluateTotal.  Without
   // this, the test would fall through to Evaluate and produce
   // fluctuation-path outputs on total-Q inputs.
   real_t bulk_bg[NUM_STATE] = {0};
   bulk_bg[SYY] =  TPV102Params::sigma_n;
   bulk_bg[SXY] = -TPV102Params::tau_ini;
   wave.SetAbsorbingBackground(bulk_bg);

   if (nfault == 0)
   {
      std::cout << "ERROR: no fault QPs found — fault bdr not detected by "
                   "WaveOperator ctor.  Check mesh.AddBdrTriangle(..., 3) "
                   "pattern.\n";
      return 1;
   }

   // Useful references for analytic expectations.
   const real_t Zp   = TPV102Params::rho * TPV102Params::cp;
   const real_t Zs   = TPV102Params::rho * TPV102Params::cs;
   const real_t tau2_0_expected   = TPV102Params::tau_ini;
   const real_t tau1_0_expected   = 0.0;
   const real_t sigma_n0_expected = TPV102Params::sigma_n;
   const real_t V_test            = 1.0e-5;   // 10 μm/s perturbation

   std::cout << "  Zp         : " << Zp  << " kg/(m^2 s)\n";
   std::cout << "  Zs         : " << Zs  << " kg/(m^2 s)\n";
   std::cout << "  V_test     : " << V_test << " m/s\n";
   std::cout << "  Zp*V_test  : " << Zp * V_test << " Pa  (T1 σ_n perturb)\n";
   std::cout << "  Zs*V_test  : " << Zs * V_test << " Pa  (T2/T3 τ perturb)\n";

   // =======================================================================
   // T0: Q = 0 — fault QPs stay at equilibrium (sanity baseline)
   // =======================================================================
   std::cout << "\n-- T0: Q = 0 at bulk; fault at initial equilibrium --\n";
   {
      ResetDOFData(dof_data, fault_coords);
      Vector Q(size); InitTotalStateQ(Q, ndof_total);
      Vector k(size);
      wave.Mult(Q, k);

      real_t worst_dsn  = 0, worst_dtau2 = 0, worst_dtau1 = 0;
      for (int i = 0; i < nfault; i++)
      {
         worst_dsn   = std::max(worst_dsn,
                                std::abs(dof_data[i].sigma_n_corr - sigma_n0_expected));
         worst_dtau2 = std::max(worst_dtau2,
                                std::abs(dof_data[i].tau2_corr - tau2_0_expected));
         worst_dtau1 = std::max(worst_dtau1,
                                std::abs(dof_data[i].tau1_corr - tau1_0_expected));
      }
      std::cout << "    worst |sigma_n_corr - 120 MPa| = " << worst_dsn   << " Pa\n";
      std::cout << "    worst |tau2_corr   - 75  MPa|  = " << worst_dtau2 << " Pa\n";
      std::cout << "    worst |tau1_corr   - 0   MPa|  = " << worst_dtau1 << " Pa\n";
      // Friction at equilibrium: V ≈ V_ini ≈ 1e-12 m/s; eta_s·V ≈ 2e-9 Pa.
      // On a 75 MPa background this is absurd FP noise; allow up to 1 Pa.
      TEST_LE(worst_dsn,   1.0, "T0 sigma_n_corr stays at sigma_n0");
      TEST_LE(worst_dtau2, 1.0, "T0 tau2_corr stays at tau2_0");
      TEST_LE(worst_dtau1, 1.0, "T0 tau1_corr stays at tau1_0");
   }

   // =======================================================================
   // T1: Pelties 7a — antisym v_y ⇒ sigma_n_trial = +Zp·V_test
   // =======================================================================
   std::cout << "\n-- T1: Pelties 7a  (v_y jump → normal-stress trial) --\n";
   {
      ResetDOFData(dof_data, fault_coords);
      Vector Q(size); InitTotalStateQ(Q, ndof_total);
      // + side (y>0, tet 1) gets +V_test; − side (y<0, tet 0) gets -V_test.
      SetAntiSymQComponent(mesh, Q, VY, ndof_total, ndof_per_elem,
                           +V_test, -V_test);
      Vector k(size);
      wave.Mult(Q, k);

      // Pelties 7a: sigma_n_trial = eta_p (v_n^- - v_n^+ + 0 + 0) where
      // v_n = -v_y on both sides.  v_n⁻ = -v_y(y<0) = +V_test;
      // v_n⁺ = -v_y(y>0) = -V_test; v_n⁻-v_n⁺ = +2V_test.
      // eta_p = Zp/2 ⇒ sigma_n_trial = Zp · V_test.
      const real_t expected_dsn = Zp * V_test;
      real_t worst_err  = 0.0;
      real_t worst_frac = 0.0;
      real_t avg_dsn = 0.0;
      for (int i = 0; i < nfault; i++)
      {
         const real_t actual_dsn = dof_data[i].sigma_n_corr - sigma_n0_expected;
         avg_dsn += actual_dsn;
         const real_t err = std::abs(actual_dsn - expected_dsn);
         worst_err = std::max(worst_err, err);
         if (std::abs(expected_dsn) > 0)
         {
            worst_frac = std::max(worst_frac, err / std::abs(expected_dsn));
         }
      }
      avg_dsn /= nfault;
      std::cout << "    expected Δσ_n  = " << expected_dsn << " Pa\n";
      std::cout << "    actual avg Δσ_n = " << avg_dsn << " Pa\n";
      std::cout << "    worst |err|    = " << worst_err << " Pa  "
                << "(rel " << worst_frac << ")\n";
      // 1 % relative tolerance: absorbs friction-solver perturbation on V2
      // and whatever FP noise accumulates along the Mult path.
      TEST_LE(worst_frac, 1.0e-2, "T1 Pelties 7a sigma_n_trial matches Zp·V_test");
   }

   // =======================================================================
   // T2: Pelties 7c — antisym v_x ⇒ tau2_trial = -Zs·V_test
   // =======================================================================
   std::cout << "\n-- T2: Pelties 7c  (v_x jump → strike-traction trial) --\n";
   {
      ResetDOFData(dof_data, fault_coords);
      Vector Q(size); InitTotalStateQ(Q, ndof_total);
      SetAntiSymQComponent(mesh, Q, VX, ndof_total, ndof_per_elem,
                           +V_test, -V_test);
      Vector k(size);
      wave.Mult(Q, k);

      // Pelties 7c: tau_2_trial = eta_s (v_t2^- - v_t2^+ + 0 + 0)
      // v_t2 = v·can_t2 = v_x (can_t2 = +x).  v_t2⁺ = +V_test, v_t2⁻ = -V_test.
      // v_t2⁻ - v_t2⁺ = -2V_test; tau_2_trial = -Zs · V_test.
      const real_t expected_dtau2 = -Zs * V_test;
      real_t worst_err  = 0.0;
      real_t worst_frac = 0.0;
      real_t avg_dtau2  = 0.0;
      for (int i = 0; i < nfault; i++)
      {
         const real_t actual_dtau2 = dof_data[i].tau2_corr - tau2_0_expected;
         avg_dtau2 += actual_dtau2;
         const real_t err = std::abs(actual_dtau2 - expected_dtau2);
         worst_err = std::max(worst_err, err);
         if (std::abs(expected_dtau2) > 0)
         {
            worst_frac = std::max(worst_frac, err / std::abs(expected_dtau2));
         }
      }
      avg_dtau2 /= nfault;
      std::cout << "    expected Δτ_2  = " << expected_dtau2 << " Pa\n";
      std::cout << "    actual avg Δτ_2 = " << avg_dtau2 << " Pa\n";
      std::cout << "    worst |err|    = " << worst_err << " Pa  "
                << "(rel " << worst_frac << ")\n";
      TEST_LE(worst_frac, 1.0e-2, "T2 Pelties 7c tau_2_trial matches -Zs·V_test");
   }

   // =======================================================================
   // T3: Pelties 7b — antisym v_z ⇒ tau1_trial = +Zs·V_test
   // =======================================================================
   std::cout << "\n-- T3: Pelties 7b  (v_z jump → dip-traction trial) --\n";
   {
      ResetDOFData(dof_data, fault_coords);
      Vector Q(size); InitTotalStateQ(Q, ndof_total);
      SetAntiSymQComponent(mesh, Q, VZ, ndof_total, ndof_per_elem,
                           +V_test, -V_test);
      Vector k(size);
      wave.Mult(Q, k);

      // Pelties 7b: tau_1_trial = eta_s (v_t1^- - v_t1^+ + 0 + 0)
      // v_t1 = v·can_t1 = -v_z (can_t1 = -z).  v_t1⁺ = -V_test, v_t1⁻ = +V_test.
      // v_t1⁻ - v_t1⁺ = +2V_test; tau_1_trial = +Zs · V_test.
      const real_t expected_dtau1 = +Zs * V_test;
      real_t worst_err  = 0.0;
      real_t worst_frac = 0.0;
      real_t avg_dtau1  = 0.0;
      for (int i = 0; i < nfault; i++)
      {
         const real_t actual_dtau1 = dof_data[i].tau1_corr - tau1_0_expected;
         avg_dtau1 += actual_dtau1;
         const real_t err = std::abs(actual_dtau1 - expected_dtau1);
         worst_err = std::max(worst_err, err);
         if (std::abs(expected_dtau1) > 0)
         {
            worst_frac = std::max(worst_frac, err / std::abs(expected_dtau1));
         }
      }
      avg_dtau1 /= nfault;
      std::cout << "    expected Δτ_1  = " << expected_dtau1 << " Pa\n";
      std::cout << "    actual avg Δτ_1 = " << avg_dtau1 << " Pa\n";
      std::cout << "    worst |err|    = " << worst_err << " Pa  "
                << "(rel " << worst_frac << ")\n";
      TEST_LE(worst_frac, 1.0e-2, "T3 Pelties 7b tau_1_trial matches +Zs·V_test");
   }

   std::cout << "\n========================================\n";
   std::cout << "  Step 6 results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
