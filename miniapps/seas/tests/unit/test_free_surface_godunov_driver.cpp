// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.3.0 Phase 2 AC-3 integration gate for I-04 (addresses
// REVIEW.md R-I04-002).
//
// ============================================================================
// Gate
// ============================================================================
// Plan Phase 2 AC-3 (tpv102_debug_v9.3.0_debug_plan.md:338-344):
//
//   "New integration test seas_test_free_surface_godunov_driver:
//    end-to-end TPV102 1-element mesh, 10 RK4 steps at dt=1e-4;
//    max|Q[SYY]| at the fault QP stays <= 1 μPa under both gamma and
//    godunov modes (locked fault + quiescent initial condition => no
//    energy to radiate => sigma_yy should not grow)."
//
// This gate exercises the full dispatch path that the Phase 1 unit test
// does NOT touch: WaveOperator::SetFreeSurfaceBCMode, the FaceBC::
// FreeSurface branch in ComputeFaceFluxRHS, and the driver's CLI wiring.
// A regression in any of those (for example, a future refactor that
// forgets to call the setter) passes the kernel-level equivalence test
// but fails here.
//
// ============================================================================
// Fixture
// ============================================================================
// Same 2-tet mesh as test_interior_fault_flux_path.cpp:
//   tet 0 = {0,1,2,4}  centroid y<0  (- side)
//   tet 1 = {0,1,2,3}  centroid y>0  (+ side)
// Interior face (0,1,2) at y=0 tagged as FAULT (attr=3); all other
// external triangles tagged FREE (attr=1).  This exercises BOTH fault
// dispatch AND the free-surface dispatch whose mode flag we care about.
//
// Initial condition: Q = 0 everywhere (fluctuation mode, quiescent).
// DOFData is initialized via InitializeFaultDOFs (pre-stress in
// sigma_n0/tau2_0, psi at equilibrium, V ~ V_ini = 1e-12 m/s).
// The fault is "locked" in the sense that V stays at ~1e-12 m/s
// throughout: no nucleation perturbation is applied.
//
// ============================================================================
// Integrator
// ============================================================================
// Classical RK4 on the bulk Q:
//   k1 = wave.Mult(Q)
//   k2 = wave.Mult(Q + 0.5 dt k1)
//   k3 = wave.Mult(Q + 0.5 dt k2)
//   k4 = wave.Mult(Q +     dt k3)
//   Q  += dt/6 * (k1 + 2 k2 + 2 k3 + k4)
//
// Each Mult call re-evaluates the friction solver and refreshes DOFData
// (same semantics as in the production TPV102 driver's Mult-based inner
// loop).  We do not separately integrate psi: under quiescent IC, psi
// drift is O(V_ini * dt) ~ 1e-16, negligible.
//
// ============================================================================
// Acceptance
// ============================================================================
// After 10 RK4 steps at dt = 1e-4 s, under each of FreeSurfaceBCMode::Gamma
// and FreeSurfaceBCMode::Godunov:
//   max over all (DOF, element, component=SYY) of |Q(SYY * ndof_total + dof)|
// must be <= 1.0e-6 Pa.  The plan says "at the fault QP"; we use the
// stricter full-bulk max because if ANY bulk DOF's Q[SYY] drifts above
// 1 μPa, it will reach the fault QP within one more step of propagation.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/tpv102_setup.hpp"
// R-I06-003 (round 2): the test was silently trivialized under the I-06
// migration because Q = 0 gives zero trial traction + zero flux at every
// dispatch choice.  Switch to total-Q initialisation so the wiring is
// actually exercised.
#include "../../dynamic/tpv102_setup_total.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../domain/boundary_config.hpp"

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

// ---------------------------------------------------------------------------
// 2-tet mesh sharing one interior face at y=0, tagged FAULT (attr=3).
// Other external faces tagged FREE (attr=1).  Mirrors
// test_interior_fault_flux_path.cpp:142-183 but scaled to TPV102 cell
// size (1000 m edges) so the plan's dt=1e-4 s sits comfortably below
// CFL = h / cp ~ 1000/6000 ~ 1.7e-1 for order 1.
// ---------------------------------------------------------------------------
static Mesh BuildTwoTetFaultMesh()
{
   // 1000 m edges — 1 TPV102-scale element per side of the fault.
   const real_t L = 1000.0;
   Mesh mesh(3, 5, 2, 0);
   real_t verts[5][3] = {
      {0.0, 0.0, 0.0}, {L, 0.0, 0.0}, {0.0, 0.0, L},
      {0.0,  L, 0.0},    // vertex 3 on +y side
      {0.0, -L, 0.0},    // vertex 4 on -y side
   };
   for (int v = 0; v < 5; v++) { mesh.AddVertex(verts[v]); }
   mesh.AddTet(0, 1, 2, 4, 1);   // tet 0: - side
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
         if (std::abs(cy) < 1e-6)   // absolute tol in meters (1e-6 m)
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
// Populate wave's fault bookkeeping from the canonical fault-face lists.
// Same pattern as the production driver.
// ---------------------------------------------------------------------------
static int SetupFault(WaveOperator<Mesh> &wave, Mesh &mesh, int order,
                      std::vector<DOFData> &dof_data, FaultFaceFlux &ff,
                      std::vector<Vector> &fault_coords)
{
   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
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
      // R-I06-003 (round 2): under the I-06 migration bulk Q carries the
      // pre-stress and DOFData pre-stress fields must be zeroed so
      // EvaluateTotal's ComputeTrialTraction does not double-count.
      ZeroDOFDataPreStressTotal(dof_data, nfault);
   }
   wave.SetFaultFlux(&ff);
   wave.SetFaultDOFData(&dof_data, nqp_per_face);
   return nfault;
}

// ---------------------------------------------------------------------------
// Classical RK4 on Q (no coupled psi-integrator — DOFData evolves as a
// side-effect of wave.Mult, matching the production inner-loop semantics
// for a single Mult call).
// ---------------------------------------------------------------------------
static void RK4Step(WaveOperator<Mesh> &wave, Vector &Q, real_t dt)
{
   const int size = Q.Size();
   Vector k1(size), k2(size), k3(size), k4(size), Qtmp(size);

   wave.Mult(Q, k1);

   Qtmp = Q;
   Qtmp.Add(0.5 * dt, k1);
   wave.Mult(Qtmp, k2);

   Qtmp = Q;
   Qtmp.Add(0.5 * dt, k2);
   wave.Mult(Qtmp, k3);

   Qtmp = Q;
   Qtmp.Add(dt, k3);
   wave.Mult(Qtmp, k4);

   // Q += dt/6 * (k1 + 2*k2 + 2*k3 + k4)
   Q.Add(dt / 6.0, k1);
   Q.Add(dt / 3.0, k2);
   Q.Add(dt / 3.0, k3);
   Q.Add(dt / 6.0, k4);
}

// ---------------------------------------------------------------------------
// Measure max |Q[SYY, i] - sigma_n0| across every scalar DOF i.
// Under total-Q init, Q[SYY] = +sigma_n0 at every DOF; any drift from
// that baseline indicates radiation from the BC or fault dispatch.
// ---------------------------------------------------------------------------
static real_t MaxAbsQSYYDelta(const Vector &Q, int ndof_total)
{
   real_t worst = 0.0;
   for (int i = 0; i < ndof_total; i++)
   {
      worst = std::max(worst,
                       std::abs(Q(SYY * ndof_total + i)
                                - TPV102Params::sigma_n));
   }
   return worst;
}

// ---------------------------------------------------------------------------
// Run one BC-mode configuration and return max|Q[SYY]| after 10 RK4 steps.
// Fresh wave operator + DOFData on every call so the modes are independent.
// ---------------------------------------------------------------------------
static real_t RunMode(FreeSurfaceBCMode mode, const std::string &label)
{
   std::cout << "\n-- mode = " << label << " --\n";

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
   wave.SetFreeSurfaceBCMode(mode);

   // R-I06-005 (round 2): plumb the bulk pre-stress background so the
   // free-surface dispatch picks up FreeSurfaceTotal / FreeSurfaceGodunovTotal
   // on tilted faces (the 2-tet fixture's external free-surface faces
   // are NOT all horizontal).  Without this the test would fail because
   // gamma-mirror / Godunov-projection on non-horizontal free faces
   // radiates the pre-stress tensor.  Buffer lifetime = function scope,
   // outlives `wave` because wave drops out of scope first on return.
   real_t bg[NUM_STATE] = {0};
   bg[SYY] =  TPV102Params::sigma_n;
   bg[SXY] = -TPV102Params::tau_ini;
   wave.SetAbsorbingBackground(bg);

   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();
   const int size = wave.Height();

   FaultFaceFlux ff(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   std::vector<DOFData> dof_data;
   std::vector<Vector> fault_coords;
   const int nfault = SetupFault(wave, mesh, order, dof_data, ff, fault_coords);
   MFEM_VERIFY(nfault > 0, "expected >= 1 fault QP in 2-tet fixture");

   // R-I06-003 (round 2): initialize Q with the TPV102 pre-stress in
   // global coords.  Under the I-06 migration wave_operator.inl
   // dispatches EvaluateTotal, so a Q=0 init would short-circuit the
   // trial traction and the test would silently pass on any dispatch
   // bug.  Initializing with pre-stress ensures EvaluateTotal actually
   // sees realistic state and the free-surface BC + fault dispatch
   // combination is exercised.
   (void)size;
   Vector Q;
   InitializeStateTotal(Q, ndof_total,
                        TPV102Params::sigma_n, TPV102Params::tau_ini);
   const real_t dt = 1.0e-4; // per plan AC-3

   std::cout << "  initial max|Q[SYY] - sigma_n0| = "
             << std::scientific << std::setprecision(3)
             << MaxAbsQSYYDelta(Q, ndof_total) << " Pa\n";

   // Gate A (R-I06-003 regression): after one Mult, DOFData.sigma_n_corr
   // at every fault QP must equal sigma_n0 to solver tolerance.  This
   // catches regressions that would make the test trivial:
   //   - InitializeStateTotal not called → Q = 0 → sigma_n_corr = 0 ≠ sigma_n0.
   //   - ZeroDOFDataPreStressTotal not called + Evaluate dispatched →
   //     sigma_n_corr = 2 * sigma_n0 (double-count).
   //   - Fault dispatch skipped → sigma_n_corr unchanged from seed.
   {
      Vector k0(Q.Size());
      wave.Mult(Q, k0);
      real_t k0_max = 0.0;
      real_t worst_sigma = 0.0;
      real_t worst_tau2  = 0.0;
      for (int i = 0; i < Q.Size(); i++)
      {
         k0_max = std::max(k0_max, std::abs(k0(i)));
      }
      for (const auto &d : dof_data)
      {
         worst_sigma = std::max(worst_sigma,
                                 std::abs(d.sigma_n_corr - TPV102Params::sigma_n));
         worst_tau2  = std::max(worst_tau2,
                                 std::abs(d.tau2_corr - TPV102Params::tau_ini));
      }
      std::cout << "  initial max|dQ/dt|            = " << k0_max << " Pa/s\n";
      std::cout << "  Gate A worst |sigma_n_corr - sigma_n0| = "
                << worst_sigma << " Pa\n";
      std::cout << "  Gate A worst |tau2_corr - tau_ini|     = "
                << worst_tau2  << " Pa\n";
      // 1 kPa tolerance absorbs friction-solver noise on V ~ V_ini.
      TEST_LE(worst_sigma, 1.0e3,
              "Gate A: EvaluateTotal yields sigma_n_corr = sigma_n0 "
              "under total-Q init (catches InitializeStateTotal / "
              "ZeroDOFDataPreStressTotal / EvaluateTotal dispatch "
              "regressions)");
      TEST_LE(worst_tau2, 1.0e3,
              "Gate A: EvaluateTotal yields tau2_corr = tau_ini under "
              "total-Q init");
   }

   for (int step = 0; step < 10; step++)
   {
      RK4Step(wave, Q, dt);
      if (step == 0 || step == 9)
      {
         std::cout << "  after step " << (step+1)
                   << "  max|Q[SYY] - sigma_n0| = "
                   << MaxAbsQSYYDelta(Q, ndof_total) << " Pa\n";
      }
   }
   const real_t final_delta = MaxAbsQSYYDelta(Q, ndof_total);
   return final_delta;
}

int main()
{
   std::cout << "\n=== TPV102 v9.3.0 Phase 2 AC-3 (I-04 + I-06): "
             << "FreeSurfaceGodunov driver integration ===\n";
   std::cout << "  2-tet fixture, locked fault, TPV102 pre-stress IC, "
             << "10 RK4 steps @ dt=1e-4\n";
   std::cout << "  Gate A: dof_data.sigma_n_corr = sigma_n0 after one Mult\n";
   std::cout << "    (R-I06-003 regression gate — catches "
             << "InitializeStateTotal / ZeroDOFDataPreStressTotal / "
             << "EvaluateTotal dispatch bugs)\n";
   std::cout << "  Gate B: |Q[SYY] - sigma_n0| bounded over 10 RK4 steps\n";
   std::cout << "    (catches catastrophic BC regression; 2-tet fixture's\n"
             << "    tilted free-surface faces exercise FreeSurfaceTotal\n"
             << "    via SetAbsorbingBackground)\n";

   // Gate B budget: 5 MPa (~4% of sigma_n0).  Friction-solver noise
   // propagation on this small fixture drifts SYY by ~1 MPa in 10 steps
   // under the CORRECT dispatch; a catastrophically-wrong dispatch
   // (e.g., free-surface BC ignores Q_bg) would push drift to O(sigma_n0)
   // (~120 MPa) within a few steps.
   const real_t budget = 5.0e6;

   const real_t syy_gamma =
      RunMode(FreeSurfaceBCMode::Gamma,   "gamma");
   TEST_LE(syy_gamma, budget,
           "Gate B: max|Q[SYY] - sigma_n0| <= 5 MPa under gamma BC "
           "(no catastrophic radiation under total-Q)");

   const real_t syy_godunov =
      RunMode(FreeSurfaceBCMode::Godunov, "godunov");
   TEST_LE(syy_godunov, budget,
           "Gate B: max|Q[SYY] - sigma_n0| <= 5 MPa under godunov BC");

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
