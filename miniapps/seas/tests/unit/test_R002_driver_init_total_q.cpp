// Round-6 R-002 regression: the driver's bulk-Q initial condition must
// carry the TPV102 pre-stress tensor (total-Q formulation), not zero
// (fluctuation-Q).
//
// The driver itself is not easily launched as a subprocess from a unit
// test (needs a mesh file + full mpirun plumbing).  This test reproduces
// the driver's init pattern standalone:
//   - Call InitializeFaultDOFs + ZeroDOFDataPreStressTotal (purpose #2).
//   - Call InitializeStateTotal (purpose #2).
//   - Call wave.SetAbsorbingBackground (purpose #2).
//
// Then verify:
//   G1: Q[SYY] ≈ +sigma_n at every DOF (was 0 under the pre-fix
//       `InitializeState`).
//   G2: Q[SXY] ≈ -tau_ini at every DOF.
//   G3: dof_data pre-stress fields all zero (ZeroDOFDataPreStressTotal
//       contract).
//   G4: wave.GetAbsorbingBackground() != nullptr (has_bulk_bg_ set).
//
// All gates fail under the pre-round-6 driver setup (which used
// `InitializeState` + no `SetAbsorbingBackground` + no
// `ZeroDOFDataPreStressTotal`).

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

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)

// 2-tet fixture with one fault face (standard pattern).
static Mesh BuildTwoTetFaultMesh(real_t L = 1000.0)
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

int main()
{
   std::cout << "\n=== Round-6 R-002 regression: driver init uses total-Q ===\n";

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
   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();

   // Fault bookkeeping.
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

   // ===== Reproduce the Round-6 driver init pattern =====
   std::vector<DOFData> dof_data;
   InitializeFaultDOFs(dof_data, nfault, fault_coords);
   ZeroDOFDataPreStressTotal(dof_data, nfault);  // Round-6 R-002

   FaultFaceFlux ff(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   wave.SetFaultFlux(&ff);
   wave.SetFaultDOFData(&dof_data, nqp_per_face);

   // Round-6 R-002: supply the bulk background so total-Q-aware BCs
   // reach EvaluateTotal / AbsorbingTotal / FreeSurface*Total paths.
   real_t Q_bg[NUM_STATE] = {0};
   Q_bg[SYY] =  TPV102Params::sigma_n;
   Q_bg[SXY] = -TPV102Params::tau_ini;
   wave.SetAbsorbingBackground(Q_bg);

   // Round-6 R-002: bulk Q = pre-stress.
   Vector Q;
   InitializeStateTotal(Q, ndof_total,
                        TPV102Params::sigma_n, TPV102Params::tau_ini);

   // ===== Gates =====
   // G1: Q[SYY] == +sigma_n at every DOF.
   real_t worst_syy = 0.0;
   for (int i = 0; i < ndof_total; i++)
   {
      worst_syy = std::max(worst_syy,
                           std::abs(Q[SYY * ndof_total + i]
                                    - TPV102Params::sigma_n));
   }
   TEST_LE(worst_syy, 1.0,
           "G1 (R-002): Q[SYY] = +sigma_n at every DOF "
           "(pre-fix: Q = 0 everywhere → fails)");

   // G2: Q[SXY] == -tau_ini at every DOF.
   real_t worst_sxy = 0.0;
   for (int i = 0; i < ndof_total; i++)
   {
      worst_sxy = std::max(worst_sxy,
                           std::abs(Q[SXY * ndof_total + i]
                                    - (-TPV102Params::tau_ini)));
   }
   TEST_LE(worst_sxy, 1.0,
           "G2 (R-002): Q[SXY] = -tau_ini at every DOF");

   // G3: dof_data pre-stress fields zeroed.
   real_t worst_pre = 0.0;
   for (int i = 0; i < nfault; i++)
   {
      worst_pre = std::max({worst_pre,
                            std::abs(dof_data[i].sigma_n0),
                            std::abs(dof_data[i].tau1_0),
                            std::abs(dof_data[i].tau2_0)});
   }
   TEST_LE(worst_pre, 0.0,
           "G3 (R-002): DOFData pre-stress fields (sigma_n0, tau*_0) "
           "all zeroed by ZeroDOFDataPreStressTotal");

   // G4: has_bulk_bg_ set.
   TEST_ASSERT(wave.GetAbsorbingBackground() != nullptr,
               "G4 (R-002): wave.GetAbsorbingBackground() != nullptr "
               "(SetAbsorbingBackground was called)");

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
