// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 pepper-bug unit-test plan (2026-04-22) — Phase 3 P3-T3 (NEW):
// Interior-fault probe AT RUPTURE-REGIME AMPLITUDES.
//
// ============================================================================
// Motivation — why this test is new
// ============================================================================
// test_interior_fault_flux_path.cpp and its asymmetric sibling
// P3-T1 (test_interior_fault_flux_path_asymmetric.cpp) both drive the
// bulk with V_test = 1e-5 m/s — seven orders of magnitude below the
// slip-rate regime during TPV102 rupture (V ≈ 1 m/s at rupture front).
// The friction solver at V_test = 1e-5 stays near the pre-nucleation
// equilibrium (V ≈ V_ini = 1e-12).  All existing interior-fault tests
// therefore probe the LINEAR-TRIAL-TRACTION path but NOT the
// strength-saturated regime the pepper-report plots show (MFEM vs
// GROG3D overlay plot tpv102_flt_0_3.png: τ_strike drops from 75 MPa
// to ~65 MPa during rupture, and τ_dip shows spurious MPa-level
// oscillations on a pure strike-slip problem).
//
// ============================================================================
// What this test does
// ============================================================================
// Identical to P3-T1 (asymmetric one-sided v_x perturbation on the
// − side of a 2-tet fault fixture), but with V_test = 1.0 m/s instead
// of 1e-5.  This drives the trial strike traction tau_2_trial into the
// few-MPa regime — a fraction of the 75 MPa pre-stress — which puts
// the friction solver in the transient rupture regime where the pepper
// pattern is observed.
//
// The test checks TWO invariants:
//
//   I1. tau1_corr (DIP traction) stays ~0 Pa (to within eta_s·V ≈ 3e5 Pa).
//       TPV102 is PURE STRIKE-SLIP: tau1_0 = 0 and all bulk Q deviations
//       are in the VX (strike velocity) component.  Any MPa-level
//       spillover into tau1_corr is a canonical-frame bug — either a
//       tangent swap (t1/t2 swapped) or a sign inconsistency.  The
//       pepper plot's spurious τ_dip oscillation matches this signature.
//
//   I2. tau2_corr (STRIKE traction) sign matches the Pelties 7c
//       prediction.  With Q[VX] = +V_test on − side only (tet 0),
//       v_t2⁻ = +V_test, v_t2⁺ = 0, so v_t2⁻ - v_t2⁺ = +V_test, giving
//       tau_2_trial = +eta_s · V_test = +0.5·Zs·V_test ≈ +4.6 MPa.
//       After friction correction, tau_2_corr must still be POSITIVE
//       (or at minimum not flip to the wrong sign of the pre-stress
//       tau2_0 = 75 MPa).  A sign-flip here is a direct pepper bug.
//
// ============================================================================
// Failure interpretation
// ============================================================================
// I1 FAIL  ⇒ MPa-level tau_dip from a pure strike-slip input.  This is
//            a canonical-frame bug — t1 (dip) and t2 (strike) are
//            swapped in the interior branch (hotspot
//            wave_operator.inl:1107-1110 or FaultFaceFlux::ComputeTrialTraction
//            line 39-64).  Directly matches the pepper-plot signature.
// I2 FAIL  ⇒ Sign inconsistency in Pelties 7c strike coupling.  Also
//            a direct pepper-plot match (tau_strike anomalies).
//
// ============================================================================
// Usage
// ============================================================================
//   ./seas_test_interior_fault_flux_path_rupture_regime   (serial, < 1 s)

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/tpv102_setup.hpp"
#include "../../dynamic/tpv102_setup_total.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../domain/boundary_config.hpp"

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

#define TEST_GE(v, bound, msg) do { \
   num_tests++; \
   double vv = (v), bb = (bound); \
   if (vv >= bb) { \
      num_passed++; \
      std::cout << "  PASSED: " << msg \
                << "  (got " << std::scientific << std::setprecision(3) \
                << vv << ", lower bound " << bb << ")\n"; \
   } else { \
      num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
                << "  (got " << std::scientific << std::setprecision(3) \
                << vv << ", expected >= " << bb << ")\n"; \
   } \
} while (0)

static Mesh BuildTwoTetFaultMesh()
{
   Mesh mesh(3, 5, 2, 0);
   real_t verts[5][3] = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 1.0},
      {0.0,  1.0, 0.0},
      {0.0, -1.0, 0.0},
   };
   for (int v = 0; v < 5; v++) { mesh.AddVertex(verts[v]); }
   mesh.AddTet(0, 1, 2, 4, 1);
   mesh.AddTet(0, 1, 2, 3, 1);
   mesh.FinalizeTopology();

   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      Array<int> fv;
      mesh.GetFaceVertices(f, fv);
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

static void SetOneSidedQComponent(const Mesh &mesh, Vector &Q, int comp,
                                   int ndof_total, int ndof_per_elem,
                                   real_t val, int side_y_sign)
{
   for (int e = 0; e < mesh.GetNE(); e++)
   {
      real_t cy = 0;
      Array<int> ev; mesh.GetElementVertices(e, ev);
      for (int v = 0; v < ev.Size(); v++) { cy += mesh.GetVertex(ev[v])[1]; }
      cy /= ev.Size();
      const bool this_side = (side_y_sign < 0) ? (cy < 0.0) : (cy > 0.0);
      if (!this_side) { continue; }
      const int off = e * ndof_per_elem;
      for (int i = 0; i < ndof_per_elem; i++)
      {
         Q(comp * ndof_total + off + i) += val;
      }
   }
}

static void ResetDOFData(std::vector<DOFData> &dof_data,
                          const std::vector<Vector> &fault_coords)
{
   const int ndof = static_cast<int>(dof_data.size());
   InitializeFaultDOFs(dof_data, ndof, fault_coords);
   ZeroDOFDataPreStressTotal(dof_data, ndof);
}

static void InitTotalStateQ(Vector &Q, int ndof_total)
{
   InitializeStateTotal(Q, ndof_total,
                        TPV102Params::sigma_n, TPV102Params::tau_ini);
}

int main()
{
   std::cout << "\n=== TPV102 Pepper-Bug Phase 3 P3-T3: "
             << "Interior-Fault AT RUPTURE-REGIME AMPLITUDES ===\n";

   Mesh mesh = BuildTwoTetFaultMesh();

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
   if (nfault == 0) { std::cout << "ERROR: no fault QPs\n"; return 1; }

   real_t bulk_bg[NUM_STATE] = {0};
   bulk_bg[SYY] =  TPV102Params::sigma_n;
   bulk_bg[SXY] = -TPV102Params::tau_ini;
   wave.SetAbsorbingBackground(bulk_bg);

   const real_t Zs   = TPV102Params::rho * TPV102Params::cs;
   const real_t tau1_0_expected   = 0.0;

   // Rupture-regime amplitude: V_test = 1 m/s (matches TPV102 rupture slip rate).
   // This drives tau2_trial = 0.5·Zs·V_test ≈ +4.6 MPa — a few percent of the
   // 75 MPa pre-stress, putting the friction solver in the strength-saturated
   // regime where the pepper pattern is observed.
   const real_t V_rupture = 1.0;     // m/s
   const real_t tau_trial_target = 0.5 * Zs * V_rupture;   // ~4.6 MPa
   std::cout << "  V_rupture  : " << V_rupture << " m/s\n";
   std::cout << "  tau2_trial : " << tau_trial_target << " Pa\n";
   std::cout << "  eta_s·V    : " << 0.5 * Zs * V_rupture << " Pa\n";

   // =========================================================================
   // T1-RUP: pure strike-slip perturbation at rupture amplitude
   //         Q[VX] = +V_rupture on − side (tet 0), 0 on + side.
   //
   // Invariants:
   //   I1: |tau1_corr - tau1_0| < eta_s · V_rupture
   //       (dip traction should be ~0; any MPa-level dip from pure strike
   //        input is a canonical-frame bug — pepper-plot signature).
   //   I2: (tau2_corr - tau2_0) > 0
   //       (strike traction increases under +v_t2⁻ jump per Pelties 7c).
   // =========================================================================
   std::cout << "\n-- T1-RUP: Q[VX] = 1 m/s on − side; tau2 drives strike only --\n";
   {
      ResetDOFData(dof_data, fault_coords);
      Vector Q(size); InitTotalStateQ(Q, ndof_total);
      SetOneSidedQComponent(mesh, Q, VX, ndof_total, ndof_per_elem,
                            +V_rupture, /*side=*/-1);
      Vector k(size);
      wave.Mult(Q, k);

      // I1: dip traction bleed-through
      real_t worst_dtau1 = 0.0;
      real_t avg_dtau1 = 0.0;
      for (int i = 0; i < nfault; i++)
      {
         const real_t d = dof_data[i].tau1_corr - tau1_0_expected;
         avg_dtau1 += d;
         worst_dtau1 = std::max(worst_dtau1, std::abs(d));
      }
      avg_dtau1 /= nfault;

      // eta_s · V_rupture is the friction-correction scale.  Allow 10% of
      // that as tolerance — a correct implementation should give 0 exactly
      // for pure strike-slip input but we absorb numerical noise.
      const real_t tol_dip = 0.1 * (0.5 * Zs * V_rupture);
      std::cout << "    avg tau1_corr    = " << avg_dtau1 << " Pa\n";
      std::cout << "    worst |Δτ_1|     = " << worst_dtau1 << " Pa\n";
      std::cout << "    tol (10% η_s·V)  = " << tol_dip << " Pa\n";
      TEST_LE(worst_dtau1, tol_dip,
              "I1: tau1_corr (dip) stays ~0 under pure strike-slip input "
              "— pepper-plot signature test");

      // I2: strike traction sign
      real_t worst_dtau2 = 0.0;
      real_t worst_dsn = 0.0;
      real_t avg_dtau2 = 0.0;
      real_t avg_dsn = 0.0;
      for (int i = 0; i < nfault; i++)
      {
         avg_dtau2 += dof_data[i].tau2_corr - TPV102Params::tau_ini;
         avg_dsn   += dof_data[i].sigma_n_corr - TPV102Params::sigma_n;
         worst_dtau2 = std::max(worst_dtau2,
                                 std::abs(dof_data[i].tau2_corr - TPV102Params::tau_ini));
         worst_dsn   = std::max(worst_dsn,
                                 std::abs(dof_data[i].sigma_n_corr - TPV102Params::sigma_n));
      }
      avg_dtau2 /= nfault;
      avg_dsn   /= nfault;
      std::cout << "    avg Δτ_2         = " << avg_dtau2 << " Pa\n";
      std::cout << "    avg Δσ_n         = " << avg_dsn << " Pa\n";
      std::cout << "    worst |Δτ_2|     = " << worst_dtau2 << " Pa\n";

      // We don't assert the exact value of tau2_corr because the friction
      // solver returns a corrected (not trial) traction: tau_corr =
      // tau_trial - eta_s · V.  V can be as large as V_rupture when the
      // friction reaches breakaway — we only check that the SIGN matches.
      // For Q[VX] = +V_rupture on − side the trial tau_2 is POSITIVE, so
      // we assert avg_dtau2 > -0.1·eta_s·V_rupture (allowing for correction
      // down to ~0 but not a sign flip to negative).
      const real_t lower = -0.1 * (0.5 * Zs * V_rupture);
      TEST_GE(avg_dtau2, lower,
              "I2: tau2_corr sign matches Pelties 7c prediction "
              "(positive for +v_t2⁻ input)");

      // I3 (bonus): normal stress should be ~unchanged for a PURE strike
      // perturbation.  Q[VX] doesn't enter Pelties 7a (v_n = -v_y, unchanged).
      // Any MPa-level Δσ_n here is a tangent-swap bug.
      std::cout << "    worst |Δσ_n|     = " << worst_dsn << " Pa\n";
      TEST_LE(worst_dsn, tol_dip,
              "I3: sigma_n_corr stays ~pre-stress under pure strike-slip input");
   }

   std::cout << "\n========================================\n";
   std::cout << "  P3-T3 rupture-regime results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
