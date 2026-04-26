// Round-11 R-1103: SetMixedFluxMode(None) bit-identity contract.
//
// Contract gate (from MIXED_FLUX_PLAN.md Phase 5 R-1103):
//   - Path A: wave_a.AdvanceADER(Q, dt, 2, Q_new_a) on a fresh mesh
//             without ever calling SetMixedFluxMode (mode at default None).
//   - Path B: wave_b.AdvanceADER(Q, dt, 2, Q_new_b) on a fresh mesh
//             AFTER calling wave_b.SetMixedFluxMode(None) (idempotent).
//   - Assert max |Q_new_a − Q_new_b| == 0 to bit precision.
//
// This guards against an implementation that DOES set state when
// SetMixedFluxMode(None) is called (e.g., wrongly populating the central
// set or mutating mixed_flux_mode_ in a way that affects dispatch).

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

namespace
{

Mesh BuildTwoTetFaultMesh(real_t L = 1000.0)
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

real_t MaxAbsDiff(const Vector &a, const Vector &b)
{
   real_t m = 0.0;
   for (int i = 0; i < a.Size(); i++)
   {
      real_t d = std::abs(a(i) - b(i));
      if (d > m) { m = d; }
   }
   return m;
}

// Run one ADER-2 macro-step from a TPV102-pre-stress IC.  Caller chooses
// whether to call SetMixedFluxMode beforehand.
void RunOneStep(WaveOperator<Mesh> &wave, std::vector<DOFData> &dof_data,
                const Mesh &mesh, int order, int ader_order,
                Vector &Q_init_out, Vector &Q_new_out)
{
   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();
   const int n_local_qps = wave.GetNumLocalFaultQPs();
   const int nbf = wave.GetNbfPerFace();

   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   std::vector<Vector> fault_coords;
   fault_coords.reserve(n_local_qps);
   for (int i = 0; i < int_faces.Size(); i++)
   {
      auto *ftr =
         const_cast<Mesh &>(mesh).GetInteriorFaceTransformations(int_faces[i]);
      const IntegrationRule &ir = IntRules.Get(
         ftr->GetGeometryType(), 2 * order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3); ftr->Face->Transform(ip, phys);
         fault_coords.push_back(phys);
      }
   }

   InitializeFaultDOFs(dof_data, n_local_qps, fault_coords);
   ZeroDOFDataPreStressTotal(dof_data, n_local_qps);

   wave.SetFaultDOFData(&dof_data, nbf);

   InitializeStateTotal(Q_init_out, ndof_total,
                        TPV102Params::sigma_n, TPV102Params::tau_ini);
   real_t bulk_bg[NUM_STATE] = {0};
   bulk_bg[SYY] =  TPV102Params::sigma_n;
   bulk_bg[SXY] = -TPV102Params::tau_ini;
   wave.SetAbsorbingBackground(bulk_bg);

   const real_t dt = 1e-4;
   wave.AdvanceADER(Q_init_out, dt, ader_order, Q_new_out);
}

} // anonymous

int main()
{
   std::cout << "\n=== Round-11 R-1103: SetMixedFluxMode(None) bit-identity ===\n";

   const int order = 1;
   const int ader_order = 2;

   BoundaryConfig bc;
   bc.natural_attrs   = {1};
   bc.fault_attr      = 3;
   bc.absorbing_attrs = {};

   // Path A: never call SetMixedFluxMode (mode at default None).
   Mesh mesh_a = BuildTwoTetFaultMesh();
   WaveOperator<Mesh> wave_a(mesh_a, order,
                             TPV102Params::lambda,
                             TPV102Params::mu,
                             TPV102Params::rho, bc);
   FaultFaceFlux ff_a(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   wave_a.SetFaultFlux(&ff_a);
   std::vector<DOFData> dof_a;
   Vector Q_init_a, Q_new_a;
   RunOneStep(wave_a, dof_a, mesh_a, order, ader_order, Q_init_a, Q_new_a);

   // Path B: explicitly call SetMixedFluxMode(None) before AdvanceADER.
   Mesh mesh_b = BuildTwoTetFaultMesh();
   WaveOperator<Mesh> wave_b(mesh_b, order,
                             TPV102Params::lambda,
                             TPV102Params::mu,
                             TPV102Params::rho, bc);
   FaultFaceFlux ff_b(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   wave_b.SetFaultFlux(&ff_b);
   std::vector<DOFData> dof_b;
   wave_b.SetMixedFluxMode(MixedFluxMode::None);
   Vector Q_init_b, Q_new_b;
   RunOneStep(wave_b, dof_b, mesh_b, order, ader_order, Q_init_b, Q_new_b);

   // Gate 1: Q_init bit-identical (sanity).
   const real_t err_init = MaxAbsDiff(Q_init_a, Q_init_b);
   TEST_LE(err_init, 0.0,
           "Q_init bit-identical between paths");

   // Gate 2: Q_new bit-identical — the main contract.
   const real_t err_new = MaxAbsDiff(Q_new_a, Q_new_b);
   std::cout << "  max |Q_new_a − Q_new_b| = " << std::scientific
             << err_new << "\n";
   TEST_LE(err_new, 0.0,
           "Q_new bit-identical between paths after AdvanceADER "
           "(SetMixedFluxMode(None) is a true no-op)");

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
