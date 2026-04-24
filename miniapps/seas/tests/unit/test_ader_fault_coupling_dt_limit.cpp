// ADER bulk-fault coupling regression: on a faulted TPV102 fixture,
// AdvanceADER must reduce to wave.Mult in the dt -> 0 limit.
//
// Why this matters:
//   Existing ADER tests cover the CK predictor on bulk waves and the
//   standalone FaultFaceFlux::EvaluateADERTotal wrapper, but they do not
//   cover the live WaveOperator path that:
//     1. reconstructs fault-local traces from bulk Q,
//     2. rotates them into the canonical frame,
//     3. dispatches EvaluateADERTotal,
//     4. rotates imposed states back, and
//     5. assembles the per-side flux into the bulk RHS.
//
// On a smooth problem, a single ADER step satisfies
//   (Q_new - Q) / dt = wave.Mult(Q) + O(dt)
// even on a faulted mesh. If that defect does NOT decay with dt on a
// 2-tet fixture that isolates the fault coupling, the ADER fault path is
// a concrete suspect for the production pepper pattern.

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

constexpr real_t kFixtureSize = 1000.0;
constexpr real_t kVelAmp      = 2.0e-5;
constexpr real_t kStressAmp   = 5.0e4;

const char *StateName(int c)
{
   static const char *kNames[NUM_STATE] =
   {
      "SXX", "SYY", "SZZ", "SXY", "SYZ", "SXZ", "VX", "VY", "VZ"
   };
   return (c >= 0 && c < NUM_STATE) ? kNames[c] : "?";
}

Mesh BuildTwoTetFaultMesh(real_t L = kFixtureSize)
{
   Mesh mesh(3, 5, 2, 0);
   real_t verts[5][3] = {
      {0.0, 0.0, 0.0}, {L, 0.0, 0.0}, {0.0, 0.0, L},
      {0.0,  L, 0.0}, {0.0, -L, 0.0},
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
         real_t cy = 0.0;
         for (int v = 0; v < fv.Size(); v++) { cy += mesh.GetVertex(fv[v])[1]; }
         cy /= fv.Size();
         if (std::abs(cy) < 1e-8)
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

int SetupFault(WaveOperator<Mesh> &wave, Mesh &mesh, int order,
               std::vector<DOFData> &dof_data, FaultFaceFlux &ff,
               std::vector<Vector> &fault_coords)
{
   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   MFEM_VERIFY(int_faces.Size() > 0, "2-tet fixture must expose a fault face");

   auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[0]);
   MFEM_VERIFY(ftr, "interior fault face FTR null");
   const int nqp_per_face =
      IntRules.Get(ftr->GetGeometryType(), 2 * order).GetNPoints();

   fault_coords.clear();
   for (int i = 0; i < int_faces.Size(); i++)
   {
      auto *face_tr = mesh.GetInteriorFaceTransformations(int_faces[i]);
      const IntegrationRule &ir =
         IntRules.Get(face_tr->GetGeometryType(), 2 * order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         face_tr->SetAllIntPoints(&ip);
         Vector phys(3);
         face_tr->Face->Transform(ip, phys);
         fault_coords.push_back(phys);
      }
   }

   const int nfault = int_faces.Size() * nqp_per_face;
   InitializeFaultDOFs(dof_data, nfault, fault_coords);
   ZeroDOFDataPreStressTotal(dof_data, nfault);
   wave.SetFaultFlux(&ff);
   wave.SetFaultDOFData(&dof_data, nqp_per_face);
   return nfault;
}

void AddSideConstant(const Mesh &mesh, const FiniteElementSpace &fes,
                     Vector &Q, int comp, real_t val_plus, real_t val_minus)
{
   const int ndof_total = fes.GetNDofs();
   for (int e = 0; e < mesh.GetNE(); e++)
   {
      Array<int> ev;
      mesh.GetElementVertices(e, ev);
      real_t cy = 0.0;
      for (int v = 0; v < ev.Size(); v++) { cy += mesh.GetVertex(ev[v])[1]; }
      cy /= ev.Size();

      const real_t val = (cy > 0.0) ? val_plus : val_minus;
      Array<int> edofs;
      fes.GetElementDofs(e, edofs);
      for (int j = 0; j < edofs.Size(); j++)
      {
         Q[comp * ndof_total + edofs[j]] += val;
      }
   }
}

void BuildExcitedState(const Mesh &mesh, const WaveOperator<Mesh> &wave, Vector &Q)
{
   const FiniteElementSpace &fes = wave.GetFESpace();
   InitializeStateTotal(Q, fes.GetNDofs(),
                        TPV102Params::sigma_n, TPV102Params::tau_ini);

   // Mixed perturbation: excite normal, strike, and dip channels without
   // pushing the 2-tet fixture into a rupture-scale regime.
   AddSideConstant(mesh, fes, Q, VX,  +kVelAmp,        -0.5 * kVelAmp);
   AddSideConstant(mesh, fes, Q, VY,  -0.7 * kVelAmp,  +0.4 * kVelAmp);
   AddSideConstant(mesh, fes, Q, VZ,  +0.3 * kVelAmp,  -0.2 * kVelAmp);
   AddSideConstant(mesh, fes, Q, SXY, +kStressAmp,     -0.6 * kStressAmp);
   AddSideConstant(mesh, fes, Q, SXZ, -0.5 * kStressAmp, +0.4 * kStressAmp);
}

struct DtResult
{
   real_t defect = 0.0;
   real_t rhs_scale = 0.0;
   int max_comp = -1;
   int max_dof = -1;
   real_t max_abs_diff = 0.0;
};

DtResult RunOneDt(real_t dt)
{
   const int order = 1;

   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr = 3;
   bc.absorbing_attrs = {};

   real_t bulk_bg[NUM_STATE] = {0.0};
   bulk_bg[SYY] =  TPV102Params::sigma_n;
   bulk_bg[SXY] = -TPV102Params::tau_ini;

   Mesh mesh_mult = BuildTwoTetFaultMesh();
   WaveOperator<Mesh> wave_mult(mesh_mult, order,
                                TPV102Params::lambda,
                                TPV102Params::mu,
                                TPV102Params::rho, bc);
   wave_mult.SetAbsorbingBackground(bulk_bg);

   FaultFaceFlux ff_mult(TPV102Params::rho,
                         TPV102Params::cp,
                         TPV102Params::cs);
   std::vector<DOFData> dof_mult;
   std::vector<Vector> fault_coords_mult;
   SetupFault(wave_mult, mesh_mult, order, dof_mult, ff_mult, fault_coords_mult);

   Vector Q_mult;
   BuildExcitedState(mesh_mult, wave_mult, Q_mult);
   Vector rhs_ref(Q_mult.Size());
   wave_mult.Mult(Q_mult, rhs_ref);

   Mesh mesh_ader = BuildTwoTetFaultMesh();
   WaveOperator<Mesh> wave_ader(mesh_ader, order,
                                TPV102Params::lambda,
                                TPV102Params::mu,
                                TPV102Params::rho, bc);
   wave_ader.SetAbsorbingBackground(bulk_bg);

   FaultFaceFlux ff_ader(TPV102Params::rho,
                         TPV102Params::cp,
                         TPV102Params::cs);
   std::vector<DOFData> dof_ader;
   std::vector<Vector> fault_coords_ader;
   SetupFault(wave_ader, mesh_ader, order, dof_ader, ff_ader, fault_coords_ader);

   Vector Q_ader;
   BuildExcitedState(mesh_ader, wave_ader, Q_ader);
   Vector Q_new(Q_ader.Size());
   wave_ader.AdvanceADER(Q_ader, dt, /*order=*/2, Q_new);

   Vector rhs_ader(Q_ader.Size());
   subtract(Q_new, Q_ader, rhs_ader);
   rhs_ader /= dt;

   DtResult out;
   const int ndof_total = wave_mult.GetFESpace().GetNDofs();
   for (int i = 0; i < rhs_ref.Size(); i++)
   {
      const real_t diff = std::abs(rhs_ader(i) - rhs_ref(i));
      if (diff > out.max_abs_diff)
      {
         out.max_abs_diff = diff;
         out.max_comp = i / ndof_total;
         out.max_dof  = i % ndof_total;
      }
      out.rhs_scale = std::max(out.rhs_scale, std::abs(rhs_ref(i)));
   }
   out.defect = out.max_abs_diff / std::max(out.rhs_scale, real_t(1.0));
   return out;
}

} // namespace

int main()
{
   std::cout << "\n=== ADER bulk-fault coupling: dt -> 0 live consistency ===\n";

   const real_t dt_list[3] =
   {
      1.0e-4,
      5.0e-5,
      2.5e-5
   };

   real_t defect[3] = {0.0, 0.0, 0.0};
   for (int i = 0; i < 3; i++)
   {
      DtResult r = RunOneDt(dt_list[i]);
      defect[i] = r.defect;
      std::cout << "  dt = " << std::scientific << std::setprecision(6)
                << dt_list[i]
                << "  rel defect = " << defect[i]
                << "  max_abs_diff = " << r.max_abs_diff
                << " at " << StateName(r.max_comp)
                << "[dof " << r.max_dof << "]"
                << "  rhs_scale = " << r.rhs_scale << "\n";
   }

   TEST_LE(defect[2], 5.0e-2,
           "Smallest-dt ADER/live defect <= 5e-2 on faulted fixture");

   const bool at_roundoff = (defect[2] < 1.0e-8);
   if (at_roundoff)
   {
      std::cout << "  defect is already at machine-precision scale; "
                << "skipping convergence-rate gate\n";
   }
   else
   {
      for (int i = 0; i < 2; i++)
      {
         const real_t rate = std::log(defect[i] / defect[i + 1]) / std::log(2.0);
         std::cout << "  observed rate(dt" << i << "->dt" << (i + 1)
                   << ") = " << std::fixed << std::setprecision(2)
                   << rate << "\n";
         TEST_LE(0.8 - rate, 0.4,
                 "ADER/live fault-coupling defect decays with dt");
      }
   }

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
