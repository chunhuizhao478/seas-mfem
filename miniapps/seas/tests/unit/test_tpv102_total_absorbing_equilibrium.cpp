// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.3.0 R-I06-001 + R-I06-006 regression gate: under the total-
// stress migration, bulk Q carries the pre-stress tensor at every DOF,
// including DOFs adjacent to absorbing (lateral / bottom / far-Y)
// boundaries.  The pre-migration GodunovFlux::Absorbing uses Q_ghost=0
// which is correct ONLY for fluctuation drivers; under total-Q it
// radiates the ~120 MPa pre-stress as a constant-in-time source of
// outgoing waves.
//
// The fix: GodunovFlux::AbsorbingTotal(nor, Q_self, Q_bg, F_h) takes
// the caller-supplied background and reduces to the interior identity
// when Q_self = Q_bg — preserving uniform pre-stress as a true
// equilibrium.  WaveOperator::SetAbsorbingBackground plumbs the
// background through; when non-null, the absorbing-BC dispatch in
// ComputeFaceFluxRHS calls AbsorbingTotal.
//
// This test builds a 2-tet fixture with:
//   - one interior face tagged as FAULT (attr=3)
//   - every other external face tagged as ABSORBING (attr=5)
// initializes Q with the TPV102 pre-stress, and verifies that
// wave.Mult(Q, k) produces k ~ 0 at every DOF — the equilibrium is
// preserved by the total-stress absorbing flux.
//
// Regression-direction:
//   - With the fix, k should be bounded by FP noise (~1e-6 * Zp *
//     sigma_n ~ 2 Pa/s).
//   - Without the fix (before R-I06-001 lands), k at absorbing-adjacent
//     DOFs contains the A^- Q_pre term, magnitude ~O(cp * sigma_n / h)
//     ~O(1e9 Pa/s) — 9 orders of magnitude above the budget.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/tpv102_setup.hpp"
#include "../../dynamic/tpv102_setup_total.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_LE(v, tol, msg) do { \
   num_tests++; \
   double vv = (v), tt = (tol); \
   if (vv <= tt) { num_passed++; \
      std::cout << "  PASSED: " << msg << "  (got " << std::scientific \
                << std::setprecision(3) << vv << ", tol " << tt << ")\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
                << "  (got " << std::scientific << std::setprecision(3) \
                << vv << ", expected <= " << tt << ")\n"; } \
} while (0)

// ---------------------------------------------------------------------------
// 2-tet mesh sharing one interior face at y=0 (tagged FAULT attr=3).
// All OTHER external triangles tagged ABSORBING (attr=5) — this is the
// key difference from BuildTwoTetFaultMesh in the other tests which use
// free-surface (attr=1).  The 1000 m edge scale keeps dt well below CFL.
// ---------------------------------------------------------------------------
static Mesh BuildTwoTetAbsorbingMesh(real_t L = 1000.0)
{
   Mesh mesh(3, 5, 2, 0);
   real_t verts[5][3] = {
      {0.0, 0.0, 0.0}, {L, 0.0, 0.0}, {0.0, 0.0, L},
      {0.0,  L, 0.0},
      {0.0, -L, 0.0},
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
         if (std::abs(cy) < 1e-6)
         {
            mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 3);   // FAULT
         }
         continue;
      }
      mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 5);          // ABSORBING
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

int main()
{
   std::cout << "\n=== TPV102 v9.3.0 R-I06-001 + R-I06-006 (I-06): "
             << "Total-Q absorbing BC preserves pre-stress equilibrium ===\n";

   Mesh mesh = BuildTwoTetAbsorbingMesh();
   BoundaryConfig bc;
   bc.natural_attrs   = {};
   bc.fault_attr      = 3;
   bc.absorbing_attrs = {5};

   const int order = 1;
   WaveOperator<Mesh> wave(mesh, order,
                           TPV102Params::lambda,
                           TPV102Params::mu,
                           TPV102Params::rho, bc);

   // Plumb the absorbing background (the fix).
   real_t absorbing_bg[NUM_STATE] = {0};
   absorbing_bg[SYY] =  TPV102Params::sigma_n;
   absorbing_bg[SXY] = -TPV102Params::tau_ini;
   wave.SetAbsorbingBackground(absorbing_bg);

   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();

   // Fault bookkeeping so the fault-face dispatch has valid DOFData.
   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   int nqp_per_face = 0;
   if (int_faces.Size() > 0)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[0]);
      nqp_per_face = IntRules.Get(ftr->GetGeometryType(), 2*order).GetNPoints();
   }
   const int nfault = int_faces.Size() * nqp_per_face;

   std::vector<Vector> fault_coords;
   fault_coords.reserve(nfault);
   for (int i = 0; i < int_faces.Size(); i++)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[i]);
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2*order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3); ftr->Face->Transform(ip, phys);
         fault_coords.push_back(phys);
      }
   }

   std::vector<DOFData> dof_data;
   InitializeFaultDOFs(dof_data, nfault, fault_coords);
   ZeroDOFDataPreStressTotal(dof_data, nfault);

   FaultFaceFlux ff(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   wave.SetFaultFlux(&ff);
   wave.SetFaultDOFData(&dof_data, nqp_per_face);

   // Initialize Q with the pre-stress (global coords).
   Vector Q;
   InitializeStateTotal(Q, ndof_total,
                        TPV102Params::sigma_n, TPV102Params::tau_ini);

   // One wave.Mult — the fix makes this equilibrium.
   Vector k(Q.Size());
   wave.Mult(Q, k);

   // Equilibrium gate.  Budget: 1e-6 * Zp * sigma_n (~2 Pa/s).  Under
   // the fix, k is bounded by FP noise from the rotation + Interior
   // chain + friction-solver perturbation on a V~V_ini locked fault.
   real_t worst_k = 0.0;
   int worst_i = -1;
   for (int i = 0; i < k.Size(); i++)
   {
      real_t a = std::abs(k(i));
      if (a > worst_k) { worst_k = a; worst_i = i; }
   }
   const real_t Zp = TPV102Params::Zp;
   const real_t budget = 1.0e-6 * Zp * TPV102Params::sigma_n;
   std::cout << "  worst |k| = " << std::scientific << std::setprecision(3)
             << worst_k << " Pa/s (at DOF " << worst_i << ")\n"
             << "  budget    = " << budget << " Pa/s (1e-6 * Zp * sigma_n)\n";

   TEST_LE(worst_k, budget,
           "R-I06-001: AbsorbingTotal preserves uniform pre-stress "
           "equilibrium (no constant-in-time radiation from absorbing "
           "faces under total-Q)");

   // Round-6 R-002: fluctuation-mode regression check removed — the
   // wave operator no longer supports the fluctuation dispatch path
   // (SetAbsorbingBackground(nullptr) would cause wave.Mult to abort
   // on the "Total-Q only" MFEM_VERIFY).  Equivalent coverage for
   // Q_bg=0 (bit-identical with the old Absorbing(Q_self) kernel) is
   // provided by the direct flux-kernel tests in
   // test_free_surface_total_direct.cpp (Gate 1/2) and the round-5
   // test_R001_rk4_fluctuation_fault_dispatch standalone flux kernel
   // tests.

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
