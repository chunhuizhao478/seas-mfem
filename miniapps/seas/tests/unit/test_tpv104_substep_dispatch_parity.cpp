// Round-10 R-1001: substep DRIVER-LEVEL dispatch parity at O=1.
//
// The T_TPV104_SSI_3 contract (per tpv104_substep_iterator.hpp:79-81):
//   "The single-sub-step case `deltaT = {dt_macro}, time_weights = {1.0}`
//    is the TPV102 one-shot limit (T_TPV104_SSI_3)."
//
// At the iterator level, this is enforced by feeding Q_pointwise[0] = Q̄
// (the time-AVERAGE).  At the driver level, the substep orchestrator must
// pick tau_nodes such that ComputeADERSubStepStates returns Q̄, NOT Q(τ=dt).
//
// Pre-R-1001 fix: tau_nodes[0] = dt_step (cumulative end-of-sub-step) →
// substep dispatch differs from one-shot at O=1 by (dt/2) · L(Q).
//
// Post-R-1001 fix: tau_nodes[0] = dt_step / 2 (sub-step midpoint) →
// for ADER-2 (Q linear in τ), Q(τ=dt/2) = Q̄ exactly, restoring SSI_3.
//
// This test runs the SAME 2-tet TPV102 locked-fault fixture used by
// `test_ader_tpv102_smoke`, advances ONE macro-step under both dispatch
// paths, and asserts Q_new bit-equal to FP precision.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/tpv102_setup.hpp"
#include "../../dynamic/tpv102_setup_total.hpp"
#include "../../dynamic/tpv104_substep_iterator.hpp"
#include "../../friction/slip_law_srw_psi.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../config/tpv104_params.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <numeric>
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

// Replicates the driver's AdvanceADERWithSubStep helper inline so this
// test does not need to link against the driver TU.  Logic is identical
// to drivers/tpv104_driver.cpp::AdvanceADERWithSubStep (post-R-1001/R-1002):
//   - midpoint sub-step nodes,
//   - per-call deltaT rescaling.
template <typename WaveOpT>
void RunSubStepDispatch(WaveOpT &wave,
                        Tpv104SubStepIterator &iterator,
                        std::vector<DOFData> &dof_data,
                        const std::vector<Vector> &fault_coords,
                        const std::vector<real_t> &V_w,
                        const Vector &Q,
                        real_t dt_step,
                        int ader_order,
                        real_t t_step_start,
                        Vector &Q_new)
{
   const std::vector<real_t> configured_deltaT  = iterator.GetDeltaT();
   const std::vector<real_t> configured_weights = iterator.GetTimeWeights();
   const int O = static_cast<int>(configured_deltaT.size());
   const real_t configured_sum =
      std::accumulate(configured_deltaT.begin(), configured_deltaT.end(),
                      static_cast<real_t>(0));
   const real_t dt_scale = dt_step / configured_sum;
   std::vector<real_t> deltaT_scaled(O);
   for (int o = 0; o < O; o++)
   {
      deltaT_scaled[o] = configured_deltaT[o] * dt_scale;
   }
   iterator.SetSubSteps(deltaT_scaled, configured_weights);
   const std::vector<real_t> &deltaT = iterator.GetDeltaT();

   // R-1001 fix: midpoint nodes, NOT cumulative end.
   std::vector<real_t> tau_nodes(O);
   real_t acc = 0.0;
   for (int o = 0; o < O; o++)
   {
      tau_nodes[o] = acc + 0.5 * deltaT[o];
      acc += deltaT[o];
   }

   std::vector<Vector> Q_per_node;
   wave.ComputeADERSubStepStates(Q, dt_step, ader_order, tau_nodes,
                                 Q_per_node);
   const int n_local = wave.GetNumLocalFaultQPs();
   std::vector<std::vector<real_t>> Q_plus(O), Q_minus(O);
   for (int o = 0; o < O; o++)
   {
      wave.EvaluateBulkAtFaultQPsCanonical(Q_per_node[o], Q_plus[o], Q_minus[o]);
   }
   const size_t n_words = static_cast<size_t>(NUM_STATE) * n_local;
   std::vector<real_t> I_imp_p(n_words, 0.0), I_imp_m(n_words, 0.0);
   if (n_local > 0)
   {
      iterator.AdvanceWithSubStepStates(dof_data, fault_coords, V_w,
                                        Q_plus, Q_minus,
                                        dt_step, t_step_start,
                                        I_imp_p.data(), I_imp_m.data(),
                                        FrictionSolver::Method::NewtonRaphsonStable);
   }
   wave.SetSubStepFaultImposedStates(
      n_local > 0 ? I_imp_p.data() : nullptr,
      n_local > 0 ? I_imp_m.data() : nullptr,
      n_local);
   wave.AdvanceADER(Q, dt_step, ader_order, Q_new);
   wave.ResetSubStepFaultImposedStates();
}

} // anonymous

int main()
{
   std::cout << "\n=== Round-10 R-1001: substep dispatch parity at O=1 ===\n";

   Mesh mesh = BuildTwoTetFaultMesh();

   BoundaryConfig bc;
   bc.natural_attrs   = {1};
   bc.fault_attr      = 3;
   bc.absorbing_attrs = {};

   const int order = 1;
   const int ader_order = 2;
   const real_t dt = 1e-4;

   // Two parallel WaveOperator instances + two DOFData copies — one for
   // each dispatch path.  They MUST be independent so neither path's
   // state leaks into the other.
   WaveOperator<Mesh> wave_one(mesh, order,
                               TPV102Params::lambda,
                               TPV102Params::mu,
                               TPV102Params::rho, bc);
   WaveOperator<Mesh> wave_sub(mesh, order,
                               TPV102Params::lambda,
                               TPV102Params::mu,
                               TPV102Params::rho, bc);

   const auto &fes = wave_one.GetFESpace();
   const int ndof_total = fes.GetNDofs();
   const int n_local_qps = wave_one.GetNumLocalFaultQPs();
   const int nbf = wave_one.GetNbfPerFace();

   const Array<int> &int_faces = wave_one.GetFaultInteriorFaces();
   std::vector<Vector> fault_coords;
   fault_coords.reserve(n_local_qps);
   for (int i = 0; i < int_faces.Size(); i++)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[i]);
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

   std::vector<DOFData> dof_one, dof_sub;
   InitializeFaultDOFs(dof_one, n_local_qps, fault_coords);
   ZeroDOFDataPreStressTotal(dof_one, n_local_qps);
   dof_sub = dof_one;

   FaultFaceFlux ff(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   wave_one.SetFaultFlux(&ff);
   wave_one.SetFaultDOFData(&dof_one, nbf);
   wave_sub.SetFaultFlux(&ff);
   wave_sub.SetFaultDOFData(&dof_sub, nbf);

   Vector Q_init;
   InitializeStateTotal(Q_init, ndof_total,
                        TPV102Params::sigma_n, TPV102Params::tau_ini);
   real_t bulk_bg[NUM_STATE] = {0};
   bulk_bg[SYY] =  TPV102Params::sigma_n;
   bulk_bg[SXY] = -TPV102Params::tau_ini;
   wave_one.SetAbsorbingBackground(bulk_bg);
   wave_sub.SetAbsorbingBackground(bulk_bg);

   // SubStepIterator setup at O=1 with single sub-step covering the macro-step.
   SlipLawSRWPsi state_evo(
      TPV104Params::a_in, TPV104Params::b, TPV104Params::V0,
      TPV104Params::f0, TPV104Params::muW, TPV104Params::V_w_out);
   state_evo.SetProductionMode();
   Tpv104SubStepIterator iterator(ff, state_evo);
   iterator.SetSubSteps({dt}, {1.0});

   std::vector<real_t> V_w(n_local_qps, TPV104Params::V_w_in);

   // ------------------------------------------------------------------
   // Path A: legacy one-shot wave.AdvanceADER.
   // ------------------------------------------------------------------
   Vector Q_new_one(Q_init.Size());
   wave_one.AdvanceADER(Q_init, dt, ader_order, Q_new_one);

   // ------------------------------------------------------------------
   // Path B: substep dispatch with O=1 deltaT={dt}, weights={1}.
   //         Post-R-1001 fix: tau_nodes[0] = dt/2 (midpoint).
   //         Path B's iterator must NOT mutate Path A's dof_data.
   // ------------------------------------------------------------------
   Vector Q_new_sub(Q_init.Size());
   RunSubStepDispatch(wave_sub, iterator, dof_sub, fault_coords, V_w,
                      Q_init, dt, ader_order, /*t_step_start=*/0.0,
                      Q_new_sub);

   // ------------------------------------------------------------------
   // Gate 1: Q_new bit-equal between paths.  Tolerance: small multiple
   // of FP precision × wave amplitude.  TPV102 pre-stress is ~1e8 Pa;
   // 1e-8 absolute matches ~1e-16 relative.
   // ------------------------------------------------------------------
   std::cout << "\n-- Gate 1: Q_new bit-equal between paths at O=1 --\n";
   const real_t err = MaxAbsDiff(Q_new_one, Q_new_sub);
   std::cout << "  max |Q_new_one - Q_new_sub| = " << std::scientific
             << err << "\n";
   TEST_LE(err, 1.0e-8,
           "R-1001: substep dispatch matches one-shot at O=1 within FP precision");

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
