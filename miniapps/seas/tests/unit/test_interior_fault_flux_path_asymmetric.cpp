// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 pepper-bug unit-test plan (2026-04-22) — Phase 3 P3-T1:
// asymmetric-Q probe for the interior-fault branch.
//
// ============================================================================
// Motivation
// ============================================================================
// The existing test_interior_fault_flux_path.cpp uses ANTI-SYMMETRIC
// velocity perturbations (Q[comp] = +V on + side, -V on − side).  In that
// regime the bulk state satisfies Q_plus = -Q_minus componentwise, so the
// interior-branch (+, −) routing at wave_operator.inl:1131-1134 reads
//
//    Q_plus_local  = elem1_on_plus ? Q_self_can : Q_nbr_can
//    Q_minus_local = elem1_on_plus ? Q_nbr_can  : Q_self_can
//
// but the subsequent Pelties eq. (7) depends only on the difference
// Q_minus - Q_plus which is the SAME regardless of which entry is in
// which slot.  The anti-symmetric test therefore cannot detect a bug
// where the (+, −) labels are swapped — a bug that would directly
// produce the observed triangle-local "pepper" on TPV102 fault
// outputs because the pre-stress + friction solve would be evaluated
// on the wrong side.
//
// This test uses ONE-SIDED perturbations (Q[comp] = +V on − side, 0 on
// + side) so Q_plus ≠ Q_minus in magnitude.  The analytic Pelties 7a-c
// predictions now have magnitude |V|/2 · Z instead of |V|/2 · 2Z, but
// critically the SIGN of the deviation from pre-stress flips if (+, −)
// are swapped.
//
// ============================================================================
// Expected values — BP5 canonical frame for TPV102
// ============================================================================
// can_n  = (0, -1, 0)  (fault normal, + side at y>0, − side at y<0)
// can_t1 = (0,  0, -1) (dip)
// can_t2 = (+1, 0, 0)  (strike)
//
// Rotation to fault-local (see test_canonical_rotation_pure_strikeslip):
//   local v_n  = -v_y
//   local v_t1 = -v_z
//   local v_t2 = +v_x
//
// T1a: Q[VY] = +V_test on − side (tet 0), 0 on + side (tet 1).
//   v_n⁻ = -V_test, v_n⁺ = 0;  v_n⁻ - v_n⁺ = -V_test
//   Pelties 7a: σ_n_trial = η_p · (-V_test) = -(Zp/2)·V_test
//   Expected: data.sigma_n_corr - σ_n0 ≈ -0.5 · Zp · V_test
//
// T2a: Q[VX] = +V_test on − side, 0 on + side.
//   v_t2⁻ = +V_test, v_t2⁺ = 0;  v_t2⁻ - v_t2⁺ = +V_test
//   Pelties 7c: τ_2_trial = η_s · (+V_test) = +(Zs/2)·V_test
//   Expected: data.tau2_corr - τ2_0 ≈ +0.5 · Zs · V_test
//
// T3a: Q[VZ] = +V_test on − side, 0 on + side.
//   v_t1⁻ = -V_test, v_t1⁺ = 0;  v_t1⁻ - v_t1⁺ = -V_test
//   Pelties 7b: τ_1_trial = η_s · (-V_test) = -(Zs/2)·V_test
//   Expected: data.tau1_corr - τ1_0 ≈ -0.5 · Zs · V_test
//
// T4a: Q[VY] = +V_test on + side, 0 on − side (opposite asymmetry).
//   v_n⁻ = 0, v_n⁺ = -V_test;  v_n⁻ - v_n⁺ = +V_test
//   Expected: data.sigma_n_corr - σ_n0 ≈ +0.5 · Zp · V_test
//   Sign FLIPS from T1a.  If the (+, −) swap is buggy, T1a and T4a will
//   produce identical or identically-wrong results.
//
// ============================================================================
// Failure interpretation
// ============================================================================
// T1a FAIL + T1 (symmetric) PASS  ⇒ (+, −) routing bug in interior branch.
//                                     Specific hotspot: wave_operator.inl:1117,
//                                     1131-1134 — the elem1_on_plus / swap logic.
// T4a produces same sign as T1a   ⇒ (+, −) routing is not sign-sensitive
//                                     at all; swap bug is confirmed.
// T2a FAIL  ⇒ strike-direction (dominant TPV102 load) is broken.
// T3a FAIL  ⇒ dip-direction is broken.
//
// ============================================================================
// Usage:  ./seas_test_interior_fault_flux_path_asymmetric   (serial, < 1 s)

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

// 2-tet fault mesh shared with test_interior_fault_flux_path.cpp
// Tet 0 = {V0, V1, V2, V4}  centroid y<0  − side
// Tet 1 = {V0, V1, V2, V3}  centroid y>0  + side
static Mesh BuildTwoTetFaultMesh()
{
   Mesh mesh(3, 5, 2, 0);
   real_t verts[5][3] = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 1.0},
      {0.0,  1.0, 0.0},
      {0.0, -1.0, 0.0},
   };
   for (int v = 0; v < 5; v++) { mesh.AddVertex(verts[v]); }
   mesh.AddTet(0, 1, 2, 4, 1);   // tet 0: − side
   mesh.AddTet(0, 1, 2, 3, 1);   // tet 1: + side
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

// Asymmetric: set Q[comp] on ONE element only (selected by `side_y_sign`).
// side_y_sign = -1 → set on tet 0 (y<0, − side); +1 → set on tet 1 (y>0, + side).
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
   std::cout << "\n=== TPV102 Pepper-Bug Phase 3 P3-T1: "
             << "Interior-Fault Asymmetric-Q Probe ===\n";

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
   if (nfault == 0)
   {
      std::cout << "ERROR: no fault QPs found\n";
      return 1;
   }

   real_t bulk_bg[NUM_STATE] = {0};
   bulk_bg[SYY] =  TPV102Params::sigma_n;
   bulk_bg[SXY] = -TPV102Params::tau_ini;
   wave.SetAbsorbingBackground(bulk_bg);

   const real_t Zp   = TPV102Params::rho * TPV102Params::cp;
   const real_t Zs   = TPV102Params::rho * TPV102Params::cs;
   const real_t tau2_0_expected   = TPV102Params::tau_ini;
   const real_t tau1_0_expected   = 0.0;
   const real_t sigma_n0_expected = TPV102Params::sigma_n;
   const real_t V_test            = 1.0e-5;   // 10 μm/s perturbation

   std::cout << "  V_test     : " << V_test << " m/s\n";
   std::cout << "  0.5·Zp·V   : " << 0.5 * Zp * V_test << " Pa\n";
   std::cout << "  0.5·Zs·V   : " << 0.5 * Zs * V_test << " Pa\n";

   // =========================================================================
   // T1a: Q[VY] on − side only — Pelties 7a with one-sided v_y
   // Expected:  Δσ_n = -0.5·Zp·V_test (SIGN: negative)
   // =========================================================================
   std::cout << "\n-- T1a: Q[VY] on − side only (asymmetric v_y) --\n";
   real_t avg_dsn_T1a = 0.0;
   {
      ResetDOFData(dof_data, fault_coords);
      Vector Q(size); InitTotalStateQ(Q, ndof_total);
      SetOneSidedQComponent(mesh, Q, VY, ndof_total, ndof_per_elem,
                            +V_test, /*side=*/-1);
      Vector k(size);
      wave.Mult(Q, k);

      const real_t expected_dsn = -0.5 * Zp * V_test;
      real_t worst_err = 0.0, worst_frac = 0.0;
      for (int i = 0; i < nfault; i++)
      {
         const real_t actual = dof_data[i].sigma_n_corr - sigma_n0_expected;
         avg_dsn_T1a += actual;
         const real_t err = std::abs(actual - expected_dsn);
         worst_err = std::max(worst_err, err);
         if (std::abs(expected_dsn) > 0)
         {
            worst_frac = std::max(worst_frac, err / std::abs(expected_dsn));
         }
      }
      avg_dsn_T1a /= nfault;
      std::cout << "    expected Δσ_n  = " << expected_dsn << " Pa\n";
      std::cout << "    actual avg Δσ_n = " << avg_dsn_T1a << " Pa\n";
      std::cout << "    worst rel err   = " << worst_frac << "\n";
      TEST_LE(worst_frac, 1.0e-2, "T1a Pelties 7a σ_n = -0.5·Zp·V (one-sided v_y on − side)");
   }

   // =========================================================================
   // T2a: Q[VX] on − side only — Pelties 7c with one-sided v_x (strike)
   // Expected:  Δτ_2 = +0.5·Zs·V_test
   // =========================================================================
   std::cout << "\n-- T2a: Q[VX] on − side only (asymmetric v_x / strike) --\n";
   {
      ResetDOFData(dof_data, fault_coords);
      Vector Q(size); InitTotalStateQ(Q, ndof_total);
      SetOneSidedQComponent(mesh, Q, VX, ndof_total, ndof_per_elem,
                            +V_test, /*side=*/-1);
      Vector k(size);
      wave.Mult(Q, k);

      const real_t expected_dtau2 = +0.5 * Zs * V_test;
      real_t worst_err = 0.0, worst_frac = 0.0;
      real_t avg = 0.0;
      for (int i = 0; i < nfault; i++)
      {
         const real_t actual = dof_data[i].tau2_corr - tau2_0_expected;
         avg += actual;
         const real_t err = std::abs(actual - expected_dtau2);
         worst_err = std::max(worst_err, err);
         if (std::abs(expected_dtau2) > 0)
         {
            worst_frac = std::max(worst_frac, err / std::abs(expected_dtau2));
         }
      }
      avg /= nfault;
      std::cout << "    expected Δτ_2  = " << expected_dtau2 << " Pa\n";
      std::cout << "    actual avg Δτ_2 = " << avg << " Pa\n";
      std::cout << "    worst rel err   = " << worst_frac << "\n";
      TEST_LE(worst_frac, 1.0e-2, "T2a Pelties 7c τ_2 = +0.5·Zs·V (one-sided v_x on − side)");
   }

   // =========================================================================
   // T3a: Q[VZ] on − side only — Pelties 7b with one-sided v_z (dip)
   // Expected:  Δτ_1 = -0.5·Zs·V_test
   // =========================================================================
   std::cout << "\n-- T3a: Q[VZ] on − side only (asymmetric v_z / dip) --\n";
   {
      ResetDOFData(dof_data, fault_coords);
      Vector Q(size); InitTotalStateQ(Q, ndof_total);
      SetOneSidedQComponent(mesh, Q, VZ, ndof_total, ndof_per_elem,
                            +V_test, /*side=*/-1);
      Vector k(size);
      wave.Mult(Q, k);

      const real_t expected_dtau1 = -0.5 * Zs * V_test;
      real_t worst_err = 0.0, worst_frac = 0.0;
      real_t avg = 0.0;
      for (int i = 0; i < nfault; i++)
      {
         const real_t actual = dof_data[i].tau1_corr - tau1_0_expected;
         avg += actual;
         const real_t err = std::abs(actual - expected_dtau1);
         worst_err = std::max(worst_err, err);
         if (std::abs(expected_dtau1) > 0)
         {
            worst_frac = std::max(worst_frac, err / std::abs(expected_dtau1));
         }
      }
      avg /= nfault;
      std::cout << "    expected Δτ_1  = " << expected_dtau1 << " Pa\n";
      std::cout << "    actual avg Δτ_1 = " << avg << " Pa\n";
      std::cout << "    worst rel err   = " << worst_frac << "\n";
      TEST_LE(worst_frac, 1.0e-2, "T3a Pelties 7b τ_1 = -0.5·Zs·V (one-sided v_z on − side)");
   }

   // =========================================================================
   // T4a: Q[VY] on + side only — opposite asymmetry from T1a
   // Expected:  Δσ_n = +0.5·Zp·V_test (SIGN: positive — FLIPPED from T1a)
   //
   // This is the critical (+, −) swap detector: if the interior branch
   // mislabels sides, T1a and T4a will produce identical sigma_n_corr
   // deviations instead of equal-and-opposite.
   // =========================================================================
   std::cout << "\n-- T4a: Q[VY] on + side only (opposite asymmetry from T1a) --\n";
   real_t avg_dsn_T4a = 0.0;
   {
      ResetDOFData(dof_data, fault_coords);
      Vector Q(size); InitTotalStateQ(Q, ndof_total);
      SetOneSidedQComponent(mesh, Q, VY, ndof_total, ndof_per_elem,
                            +V_test, /*side=*/+1);
      Vector k(size);
      wave.Mult(Q, k);

      const real_t expected_dsn = +0.5 * Zp * V_test;
      real_t worst_err = 0.0, worst_frac = 0.0;
      for (int i = 0; i < nfault; i++)
      {
         const real_t actual = dof_data[i].sigma_n_corr - sigma_n0_expected;
         avg_dsn_T4a += actual;
         const real_t err = std::abs(actual - expected_dsn);
         worst_err = std::max(worst_err, err);
         if (std::abs(expected_dsn) > 0)
         {
            worst_frac = std::max(worst_frac, err / std::abs(expected_dsn));
         }
      }
      avg_dsn_T4a /= nfault;
      std::cout << "    expected Δσ_n  = " << expected_dsn << " Pa\n";
      std::cout << "    actual avg Δσ_n = " << avg_dsn_T4a << " Pa\n";
      std::cout << "    worst rel err   = " << worst_frac << "\n";
      TEST_LE(worst_frac, 1.0e-2, "T4a Pelties 7a σ_n = +0.5·Zp·V (one-sided v_y on + side)");
   }

   // =========================================================================
   // T5: SIGN ANTI-SYMMETRY between T1a and T4a — the (+, −) swap detector.
   //
   // Under a correct (+, −) labelling, flipping which side of the fault
   // carries the perturbation must flip the sign of the trial normal
   // traction.  avg_dsn_T1a should be ≈ -avg_dsn_T4a.  If they match
   // sign OR are grossly unequal in magnitude, the (+, −) routing has a
   // bug — a direct pepper-bug signature.
   // =========================================================================
   std::cout << "\n-- T5: Sign anti-symmetry T1a vs T4a (pepper-bug detector) --\n";
   {
      const real_t sum = avg_dsn_T1a + avg_dsn_T4a;
      const real_t mag = 0.5 * (std::abs(avg_dsn_T1a) + std::abs(avg_dsn_T4a));
      const real_t rel_sum = (mag > 0.0) ? std::abs(sum) / mag : 0.0;
      std::cout << "    avg Δσ_n T1a = " << avg_dsn_T1a << " Pa\n";
      std::cout << "    avg Δσ_n T4a = " << avg_dsn_T4a << " Pa\n";
      std::cout << "    T1a + T4a     = " << sum << " Pa  (relative " << rel_sum << ")\n";
      TEST_LE(rel_sum, 1.0e-2,
              "T5 avg(T1a) + avg(T4a) ≈ 0 — (+, −) routing is sign-consistent");
   }

   std::cout << "\n========================================\n";
   std::cout << "  Phase 3 P3-T1 results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
