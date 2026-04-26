// Round-7 R-602/R-603: SubStep dispatch guard contract.
//
// This test validates two contracts of the substep integration in
// `wave_operator.inl::ComputeADERFaceFluxRHS`:
//
//   1. NO-OP default: when SetSubStepFaultImposedStates has NEVER been
//      called (default state), AdvanceADER produces results identical
//      to the pre-R-602 baseline.  We verify by running once via the
//      default path, then once after Set/Reset cycle (which leaves the
//      pointers null again), and asserting bit-identical Q_new.
//
//   2. GUARD-ACTIVE substitution: when the iterator side-channel
//      pointers are set, the fault branch substitutes the supplied
//      imposed states for the inline EvaluateADER output.  We verify
//      by setting the pointers to a CHOSEN non-default value and
//      asserting Q_new differs from the default.
//
// Mesh: same 2-tet locked-fault fixture as test_ader_tpv102_smoke.
// The fault flux is wired so the inline path runs, but with a fixed
// total-Q pre-stress so the resulting Q_new is non-trivial enough to
// detect a substitution.

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

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)

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

// 2-tet fault mesh, same as test_ader_tpv102_smoke.
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

// Run one ADER step with the WaveOperator already configured.  Returns
// Q_new.  This wraps the boilerplate to make the gates readable.
Vector RunOneStep(WaveOperator<Mesh> &wave, const Vector &Q,
                  real_t dt, int order)
{
   Vector Q_new(Q.Size());
   wave.AdvanceADER(Q, dt, order, Q_new);
   return Q_new;
}

} // anonymous

int main()
{
   std::cout << "\n=== Round-7 R-602/R-603: SubStep dispatch guard ===\n";

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

   const int n_local_qps = wave.GetNumLocalFaultQPs();
   std::cout << "  ndof_total = " << ndof_total
             << ", n_local_fault_qps = " << n_local_qps << "\n";
   TEST_ASSERT(n_local_qps > 0,
               "2-tet fixture has at least one local fault QP");

   // Per-QP fault coords (driver pattern).
   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   const int nbf = wave.GetNbfPerFace();
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

   std::vector<DOFData> dof_data;
   InitializeFaultDOFs(dof_data, n_local_qps, fault_coords);
   ZeroDOFDataPreStressTotal(dof_data, n_local_qps);

   FaultFaceFlux ff(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   wave.SetFaultFlux(&ff);
   wave.SetFaultDOFData(&dof_data, nbf);

   Vector Q;
   InitializeStateTotal(Q, ndof_total,
                        TPV102Params::sigma_n, TPV102Params::tau_ini);
   real_t bulk_bg[NUM_STATE] = {0};
   bulk_bg[SYY] =  TPV102Params::sigma_n;
   bulk_bg[SXY] = -TPV102Params::tau_ini;
   wave.SetAbsorbingBackground(bulk_bg);

   const real_t dt = 1e-4;
   const int ader_order = 2;

   // ------------------------------------------------------------------
   // Gate 1: NO-OP DEFAULT — first AdvanceADER call (no Set/Reset cycle).
   // ------------------------------------------------------------------
   std::cout << "\n-- Gate 1: NO-OP default — baseline AdvanceADER --\n";
   std::vector<DOFData> dof_data_baseline = dof_data;  // snapshot
   wave.SetFaultDOFData(&dof_data_baseline, nbf);
   Vector Q_baseline = RunOneStep(wave, Q, dt, ader_order);

   // ------------------------------------------------------------------
   // Gate 2: SET / RESET CYCLE WITH NULLPTR (ResetSubStepFaultImposedStates
   //         leaves pointers null) — bit-identical to baseline.
   // ------------------------------------------------------------------
   std::cout << "\n-- Gate 2: Reset() then AdvanceADER bit-identical --\n";
   std::vector<DOFData> dof_data_reset = dof_data;
   wave.SetFaultDOFData(&dof_data_reset, nbf);
   wave.ResetSubStepFaultImposedStates();
   Vector Q_reset = RunOneStep(wave, Q, dt, ader_order);

   const real_t err_reset = MaxAbsDiff(Q_baseline, Q_reset);
   TEST_LE(err_reset, 0.0,
           "Q_reset bit-identical to Q_baseline (Reset leaves pointers "
           "null → inline EvaluateADER path)");

   // ------------------------------------------------------------------
   // Gate 3: GUARD-ACTIVE — install non-zero per-QP imposed states.
   //         Verify Q_new differs from baseline (substitution observed).
   //         The injected I_imp is scaled to match the time-integrated
   //         shape (dt × Q_imp) so the fault branch's downstream
   //         T_can rotation produces a sane bulk update.
   // ------------------------------------------------------------------
   std::cout << "\n-- Gate 3: Guard active — injection differs from baseline --\n";
   std::vector<DOFData> dof_data_substep = dof_data;
   wave.SetFaultDOFData(&dof_data_substep, nbf);

   const size_t n_words = static_cast<size_t>(NUM_STATE) * n_local_qps;
   std::vector<real_t> I_imp_plus(n_words, 0.0);
   std::vector<real_t> I_imp_minus(n_words, 0.0);
   // Inject CONSTANT non-zero canonical-frame imposed state on every QP.
   // For TPV102's locked fixture, Q_imp ≈ (sigma_n, tau_ini, 0, ...)
   // in canonical frame — pre-stress only.  Multiply by dt to put into
   // I-units, then scale by 0.5 so the substitution is detectably
   // different from the inline EvaluateADER value (which produces the
   // physically-correct pre-stress Q_imp).
   const real_t inject_scale = 0.5 * dt;
   for (int i = 0; i < n_local_qps; i++)
   {
      I_imp_plus [i * NUM_STATE + SXX] =  inject_scale * TPV102Params::sigma_n;
      I_imp_plus [i * NUM_STATE + SXY] = -inject_scale * TPV102Params::tau_ini;
      I_imp_minus[i * NUM_STATE + SXX] =  inject_scale * TPV102Params::sigma_n;
      I_imp_minus[i * NUM_STATE + SXY] = -inject_scale * TPV102Params::tau_ini;
   }
   wave.SetSubStepFaultImposedStates(I_imp_plus.data(),
                                     I_imp_minus.data(),
                                     n_local_qps);
   Vector Q_substep = RunOneStep(wave, Q, dt, ader_order);
   wave.ResetSubStepFaultImposedStates();

   const real_t err_substep = MaxAbsDiff(Q_baseline, Q_substep);
   std::cout << "    max |Q_substep − Q_baseline| = " << std::scientific
             << err_substep << "\n";
   TEST_ASSERT(err_substep > 0.0,
               "Q_substep differs from Q_baseline when guard is active "
               "(substitution observed)");

   // ------------------------------------------------------------------
   // Gate 4: AFTER RESET, NEXT AdvanceADER reverts to baseline.
   // ------------------------------------------------------------------
   std::cout << "\n-- Gate 4: After Reset, AdvanceADER returns to baseline --\n";
   // Note: dof_data_substep was already mutated by Gate 3.  Use a fresh
   // copy to isolate the wave-operator side-channel state.
   std::vector<DOFData> dof_data_post = dof_data;
   wave.SetFaultDOFData(&dof_data_post, nbf);
   // Reset is idempotent and called twice (once after Gate 3, once here).
   wave.ResetSubStepFaultImposedStates();
   Vector Q_post = RunOneStep(wave, Q, dt, ader_order);
   const real_t err_post = MaxAbsDiff(Q_baseline, Q_post);
   TEST_LE(err_post, 0.0,
           "Q_post bit-identical to Q_baseline after Reset (no leftover "
           "side-channel state)");

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
