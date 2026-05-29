// ===========================================================================
// Phase 13 — C-6 fault-bearing parity test.
//
// The single tripwire that catches the whole `flux_`-placeholder leak class
// (REVIEW rounds 3-4 R-001).  On a mesh WITH a fault, a `BimaterialWaveOperator`
// built from a HOMOGENEOUS Coefficient material must match the scalar
// `WaveOperator` `Mult` (RK4) and `AdvanceADER` (predictor-corrector) to LU
// rounding (1e-9), because every material access on the bimaterial object
// routes through the per-element flux pool (which holds the EXACT same
// constants as the scalar `flux_`).
//
// The bimaterial object's inherited `flux_` is the deliberately-unphysical
// (1,1,1) sentinel (set by BimaterialWaveOperator's ctor delegation).  So if
// ANY skeleton site — the fault imposed-state flux (the 8 R-001 sites), the
// boundary flux, the bulk volume RHS, or the ADER CK recursion — bare-used
// `flux_` instead of the overridden `FluxForElem_(e)`, its contribution would
// be computed from (1,1,1) and DIFFER from the real-material scalar reference
// → this test fails at the offending DOFs.  A missed material site therefore
// becomes a hard parity failure, never silent wrong physics.  The test also
// asserts the bimaterial output is finite (no NaN) — the (1,1,1) sentinel is
// valid-but-wrong, never NaN, so a leak shows as a finite parity miss.
//
// What this covers that C-5 (test_phaseh_wave_operator_constant_parity, a
// fault-FREE box) does NOT: the fault imposed-state flux.  C-5 covers the
// interior non-fault faces, boundary, volume, and CK on the matrix path; this
// adds the fault.
//
// SCOPE / COVERAGE NOTES (REVIEW R-004, R-005):
//  - HOMOGENEOUS material: with equal material on both fault-adjacent elements
//    the per-side elem_plus/elem_minus selection is a no-op, so this test does
//    NOT exercise the per-side / bi-material fault-flux asymmetry — that
//    correctness rests on faithfulness to hrs-ref (git diff hrs-ref HEAD), not
//    on C-6.  A true bi-material fault regression case would be needed.
//  - SERIAL Mult + AdvanceADER only: the substep dispatch fault sites
//    (ComputeADERFaceFluxRHS per-QP `elem_plus_qq`, and the parallel-only
//    ComputeADERSharedFaceFluxRHS `qa.local_elem`) are NOT exercised here; the
//    latter needs np>1.  Their R-001 routing is verified by hrs-ref
//    faithfulness + inspection only.
// ===========================================================================

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/bimaterial_wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/heterogeneous_material.hpp"
#include "../../dynamic/tpv102_setup.hpp"
#include "../../dynamic/tpv102_setup_total.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do {                                         \
   num_tests++;                                                             \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; }    \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: "      \
                                  << msg << "\n"; }                         \
} while (0)

// ---------------------------------------------------------------------------
// 2-tet mesh sharing one interior face at y=0, tagged FAULT (attr=3); the
// external faces are FREE (attr=1).  Identical to
// test_interior_fault_flux_path.cpp::BuildTwoTetFaultMesh.
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

// Wire the fault into a WaveOperator (works for both the scalar base and the
// bimaterial subclass through a base reference).  Returns the fault QP count.
static int SetupFault(WaveOperator<Mesh> &wave, Mesh &mesh, int order,
                      std::vector<DOFData> &dof_data, FaultFaceFlux &ff,
                      std::vector<Vector> &fault_coords)
{
   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   MFEM_VERIFY(wave.GetFaultSharedFaces().Size() == 0,
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
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2*order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3); ftr->Face->Transform(ip, phys);
         fault_coords.push_back(phys);
      }
   }

   if (nfault > 0) { InitializeFaultDOFs(dof_data, nfault, fault_coords); }
   wave.SetFaultFlux(&ff);
   wave.SetFaultDOFData(&dof_data, nqp_per_face);
   return nfault;
}

// Reset DOFData to clean total-stress equilibrium (zero pre-stress fields —
// the bulk Q carries them under the total-stress migration).
static void ResetDOFData(std::vector<DOFData> &dof_data,
                         const std::vector<Vector> &fault_coords)
{
   const int ndof = static_cast<int>(dof_data.size());
   InitializeFaultDOFs(dof_data, ndof, fault_coords);
   ZeroDOFDataPreStressTotal(dof_data, ndof);
}

// tet 1 (centroid y>0) gets val_plus; tet 0 (centroid y<0) gets val_minus.
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

// FIELD-SCALE relative difference: |a-b| / max_i(|a_i|, |b_i|).
//
// Why not the per-component metric C-5 uses (|a_i-b_i|/max(|a_i|,|b_i|,1)):
// C-5 drives a fully non-equilibrium state (every dQ/dt component O(1e8)), so
// the per-component metric is well-conditioned there.  C-6's fault fixture has
// near-equilibrium DOFs where dQ/dt ≈ 0 is computed as a near-cancellation of
// O(1e8) imposed-state terms; the matrix path's per-element GodunovFlux carries
// a 1-ULP-different material (the GodunovFluxPool stores round_sig(.,6) values),
// so those cancelling DOFs land ~6e-4 apart in ABSOLUTE terms — a huge per-
// component "relative" error on a ~0 quantity, but ~1e-15 of the field's
// dominant magnitude.  The field-scale metric is the correct relative-error
// measure for a mixed-magnitude vector field and is what makes "matches to
// 1e-9" meaningful here.  A genuine (1,1,1) sentinel leak changes the fault/
// boundary/volume flux by O(1) of the field magnitude, i.e. field-scale rel
// ~ O(1) >> 1e-9 — so this metric still catches the leak class decisively.
static real_t MaxRelDiff(const Vector &a, const Vector &b, int &n_over_tol,
                         real_t tol, bool &all_finite)
{
   real_t scale = 1.0;
   all_finite = true;
   for (int i = 0; i < a.Size(); i++)
   {
      if (!std::isfinite(b(i)) || !std::isfinite(a(i))) { all_finite = false; }
      scale = std::max(scale, std::max(std::abs(a(i)), std::abs(b(i))));
   }
   real_t max_rel = 0.0;
   n_over_tol = 0;
   for (int i = 0; i < a.Size(); i++)
   {
      const real_t rel = std::abs(a(i) - b(i)) / scale;
      max_rel = std::max(max_rel, rel);
      if (rel > tol) { ++n_over_tol; }
   }
   return max_rel;
}

int main()
{
   std::cout << "\n=== Phase 13 C-6: BimaterialWaveOperator fault-bearing "
                "parity (R-001 tripwire) ===\n";

   Mesh mesh = BuildTwoTetFaultMesh();

   BoundaryConfig bc;
   bc.natural_attrs   = {1};
   bc.fault_attr      = 3;
   bc.absorbing_attrs = {};

   const int order = 1;
   // Bulk material MUST be exactly representable in 6 significant figures so
   // the GodunovFluxPool dedup (round_sig(.,6)) is a no-op and pool->At(e)
   // equals the scalar flux_ bit-for-bit — otherwise the homogeneous parity
   // floor is the ~1e-6 dedup rounding (amplified in the RHS), not 1e-9.  TPV102's
   // mu = rho*cs^2 / lambda = rho*cp^2-2mu are NOT 6-sig-fig-exact, so we use
   // round, TPV102-flavoured constants (matching C-5's 32.0e9 choice) and derive
   // the fault impedance cp/cs from them.  The fault friction physics is
   // identical for both operators, so the parity is driven purely by the bulk
   // material — physical self-consistency of the fixture is not required.
   const real_t lam = 32.0e9;     // 6-sig-fig-exact
   const real_t mu  = 32.0e9;     // 6-sig-fig-exact
   const real_t rho = 2670.0;     // 6-sig-fig-exact
   const real_t cp  = std::sqrt((lam + 2.0 * mu) / rho);
   const real_t cs  = std::sqrt(mu / rho);

   // Scalar reference operator (real material via flux_).
   WaveOperator<Mesh> wave_scalar(mesh, order, lam, mu, rho, bc);

   // Bimaterial operator from a HOMOGENEOUS Coefficient material: per-element
   // pool holds the EXACT (lam, mu, rho); the inherited flux_ is the (1,1,1)
   // sentinel.  Coefficients must outlive the operator.
   ConstantCoefficient lam_c(lam), mu_c(mu), rho_c(rho);
   BimaterialWaveOperator<Mesh> wave_bimat(
      mesh, order, MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c), bc);

   TEST_ASSERT(wave_scalar.UsesGodunovFluxPool() == false,
               "scalar operator: no flux pool");
   TEST_ASSERT(wave_bimat.UsesGodunovFluxPool() == true,
               "bimaterial operator: per-element flux pool present");

   const int size = wave_scalar.Height();
   const auto &fes = wave_scalar.GetFESpace();
   const int ndof_total    = fes.GetNDofs();
   const int ndof_per_elem = fes.GetFE(0)->GetDof();

   // Separate fault state per operator (each Mult advances its own DOFData via
   // the friction solve), identical FaultFaceFlux impedance for both.
   FaultFaceFlux ff_s(rho, cp, cs);
   FaultFaceFlux ff_b(rho, cp, cs);
   std::vector<DOFData> dof_s, dof_b;
   std::vector<Vector>  coords_s, coords_b;
   const int nfault_s = SetupFault(wave_scalar, mesh, order, dof_s, ff_s, coords_s);
   const int nfault_b = SetupFault(wave_bimat,  mesh, order, dof_b, ff_b, coords_b);
   TEST_ASSERT(nfault_s > 0 && nfault_s == nfault_b,
               "both operators see the same nonzero fault QP count");
   std::cout << "  Fault QPs  : " << nfault_s << "\n";

   // Total-Q background so the fault dispatch uses EvaluateTotal (matches the
   // production driver path that the matrix operator runs under).
   real_t bulk_bg[NUM_STATE] = {0};
   bulk_bg[SYY] =  TPV102Params::sigma_n;
   bulk_bg[SXY] = -TPV102Params::tau_ini;
   wave_scalar.SetAbsorbingBackground(bulk_bg);
   wave_bimat.SetAbsorbingBackground(bulk_bg);

   const real_t k_tol = 1e-9;
   const real_t V_test = 1.0e-5;   // strike-direction slip-rate perturbation

   // -----------------------------------------------------------------------
   // C-6a: Mult (RK4 path) parity — exercises ComputeFaceFluxRHS's interior-
   // fault imposed-state flux (R-001 site 1), boundary flux, volume RHS.
   // -----------------------------------------------------------------------
   std::cout << "\n-- C-6a: Mult (RK4) parity --\n";
   {
      ResetDOFData(dof_s, coords_s);
      ResetDOFData(dof_b, coords_b);
      Vector Q(size); InitializeStateTotal(Q, ndof_total,
                                           TPV102Params::sigma_n,
                                           TPV102Params::tau_ini);
      // Strike-direction (VX) velocity jump across the y=0 fault.
      SetAntiSymQComponent(mesh, Q, VX, ndof_total, ndof_per_elem,
                           +0.5 * V_test, -0.5 * V_test);

      Vector k_s(size), k_b(size);
      wave_scalar.Mult(Q, k_s);
      wave_bimat.Mult(Q, k_b);

      // Determinism control: a SECOND scalar operator with its own fault state
      // must reproduce wave_scalar bit-for-bit, so any scalar-vs-bimat diff is
      // genuinely the matrix path (not harness nondeterminism).
      {
         WaveOperator<Mesh> wave_s2(mesh, order, lam, mu, rho, bc);
         FaultFaceFlux ff_s2(rho, cp, cs);
         std::vector<DOFData> dof_s2; std::vector<Vector> coords_s2;
         SetupFault(wave_s2, mesh, order, dof_s2, ff_s2, coords_s2);
         wave_s2.SetAbsorbingBackground(bulk_bg);
         ResetDOFData(dof_s2, coords_s2);
         Vector k_s2(size); wave_s2.Mult(Q, k_s2);
         real_t worst2 = 0;
         for (int i = 0; i < size; i++)
         { worst2 = std::max(worst2, std::abs(k_s(i) - k_s2(i))); }
         TEST_ASSERT(worst2 == 0.0,
                     "C-6a: two scalar operators are bit-identical (determinism)");
      }

      int n_over = 0; bool finite = true;
      const real_t max_rel = MaxRelDiff(k_s, k_b, n_over, k_tol, finite);
      TEST_ASSERT(finite, "C-6a: bimaterial Mult output is finite (no NaN "
                          "from the (1,1,1) sentinel)");
      if (n_over != 0)
      {
         std::cerr << "  DIFF: " << n_over << " / " << size
                   << " over tol; max rel = " << max_rel << "\n";
      }
      TEST_ASSERT(n_over == 0,
                  "C-6a: bimaterial (homogeneous Coefficient) Mult matches "
                  "scalar to 1e-9 (fault imposed-state flux R-001 + boundary "
                  "+ volume)");
   }

   // -----------------------------------------------------------------------
   // C-6b: AdvanceADER (predictor-corrector) parity — exercises
   // ComputeADERFaceFluxRHS's interior-fault imposed-state flux (R-001 site 3),
   // the ADER CK recursion (ApplyElementJacobian_), and the ADER volume update.
   // -----------------------------------------------------------------------
   std::cout << "\n-- C-6b: AdvanceADER (order 2) parity --\n";
   {
      ResetDOFData(dof_s, coords_s);
      ResetDOFData(dof_b, coords_b);
      Vector Q(size); InitializeStateTotal(Q, ndof_total,
                                           TPV102Params::sigma_n,
                                           TPV102Params::tau_ini);
      SetAntiSymQComponent(mesh, Q, VX, ndof_total, ndof_per_elem,
                           +0.5 * V_test, -0.5 * V_test);

      const real_t dt = 1.0e-4;
      Vector Qn_s(size), Qn_b(size);
      wave_scalar.AdvanceADER(Q, dt, 2, Qn_s);
      wave_bimat.AdvanceADER(Q, dt, 2, Qn_b);

      int n_over = 0; bool finite = true;
      const real_t max_rel = MaxRelDiff(Qn_s, Qn_b, n_over, k_tol, finite);
      TEST_ASSERT(finite, "C-6b: bimaterial AdvanceADER output is finite "
                          "(no NaN from the (1,1,1) sentinel)");
      if (n_over != 0)
      {
         std::cerr << "  DIFF: " << n_over << " / " << size
                   << " over tol; max rel = " << max_rel << "\n";
      }
      TEST_ASSERT(n_over == 0,
                  "C-6b: bimaterial (homogeneous Coefficient) AdvanceADER "
                  "matches scalar to 1e-9 (fault imposed-state flux R-001 + "
                  "CK recursion + ADER volume)");
   }

   // -----------------------------------------------------------------------
   // C-6c: ComputeMaxDt parity — the per-element walk over the homogeneous
   // pool must equal the scalar h_min_/cp formula.
   // -----------------------------------------------------------------------
   std::cout << "\n-- C-6c: ComputeMaxDt parity --\n";
   {
      const real_t cfl = 0.4;
      const real_t dt_s = wave_scalar.ComputeMaxDt(cfl);
      const real_t dt_b = wave_bimat.ComputeMaxDt(cfl);
      const real_t rel = std::abs(dt_s - dt_b)
         / std::max(std::abs(dt_s), std::max(std::abs(dt_b), real_t(1.0)));
      TEST_ASSERT(std::isfinite(dt_b) && rel <= 1e-12,
                  "C-6c: bimaterial ComputeMaxDt matches scalar (homogeneous)");
   }

   std::cout << "\n==============================================\n";
   std::cout << "Summary: " << num_passed << " / " << num_tests
             << " passed (" << num_failed << " failed)\n";
   std::cout << "==============================================\n";
   return (num_failed > 0) ? 1 : 0;
}
