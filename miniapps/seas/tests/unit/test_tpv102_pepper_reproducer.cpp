// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 pepper reproducer — R-001 gate for the v9.4.0 split-prestress refactor
// (companion to debug_document/tpv102_debug_document/
//  tpv102_debug_v9.4.0_debug_plan.md).
//
// CLAIM UNDER TEST (v9.4.0 plan §1, "Root Cause Hypothesis"):
//
//   For a uniform bulk state Q[c] = C (with C a large scalar like
//   sigma_n0 ~ 1.2e8 Pa under total-Q), the fault-QP trace
//   reconstruction
//
//       Q_self(q) = sum_i shape_i(q) * Q[dof_i]
//
//   should produce identical values at every face quadrature point
//   because shape functions form a partition of unity.  In floating
//   point, however, the shape-weight permutation across QPs (i.e., at
//   different QPs, the element DOFs carry different shape values in
//   possibly different orders in the accumulation) leads to ULP-scale
//   per-QP dispersion that scales linearly with C.
//
//   The plan argues this dispersion is the SEED for the pepper
//   artifact: it is re-created every time step (forced response), and
//   the nonlinear friction solve amplifies it into visible per-element
//   noise under rupture drive.
//
// WHY GATE THE REFACTOR ON THIS TEST:
//   Under total-Q, bulk Q carries C ~ 1.2e8 Pa.
//   Under split-prestress (fluctuation-Q), bulk Q carries dynamic
//   fluctuations of O(1e6 Pa).  If dispersion scales linearly with C,
//   split-prestress cuts the seed by ~100x, which is the plan's
//   justification for a 5-commit architectural refactor.
//
//   If dispersion does NOT exist (a), is NOT ULP-bounded (b), or does
//   NOT scale linearly with C (c), the ULP hypothesis is wrong — the
//   pepper seed lives elsewhere (friction conditioning, rotation
//   rebuild, MPI drift, …) and the v9.4.0 refactor is misdirected.
//
// THREE SUBTESTS (R-001 acceptance):
//   (a) Dispersion is nonzero for at least one C in the sweep (the ULP
//       mechanism manifests on this platform/compiler).
//   (b) Dispersion < 100 * ULP(C) for every C (bounded at the expected
//       ULP scale; a larger dispersion would indicate a genuine
//       non-partition-of-unity bug, not an ULP artifact).
//   (c) Relative dispersion is roughly constant across C (linear
//       scaling; ratio of max rel-dispersion to min rel-dispersion
//       over the sweep < 10).
//
// FIXTURE:
//   Matches TPV102 production as closely as possible:
//     - Cartesian 3D mesh split into tetrahedra (MFEM's default tet
//       generation gives the same element topology as the TPV102 mesh
//       at a fault QP — 3 shape values active per face from a 4-vertex
//       tet).
//     - Order-1 L2 finite element, BasisType::GaussLobatto (matches
//       wave_operator.inl:34 construction).
//     - Dunavant face rule via IntRules.Get(TRIANGLE, 2*order).
//
// Runs in <1 s on the local laptop, no MPI needed.

#include "mfem.hpp"

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

using namespace mfem;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)

namespace {

// Per-C sweep result returned from ProbeOneMagnitude.
struct SweepResult
{
   real_t C;
   real_t min_Q;
   real_t max_Q;
   real_t abs_dispersion;   // max - min
   real_t rel_dispersion;   // abs_dispersion / C
};

// Probe one magnitude C: set every element DOF to C, interpolate at each
// face QP, return min/max/dispersion of the interpolated value.
SweepResult ProbeOneMagnitude(const FiniteElement &fe,
                              FaceElementTransformations &ftr,
                              const IntegrationRule &ir,
                              real_t C)
{
   const int ndof = fe.GetDof();
   Vector Q_dofs(ndof);
   Q_dofs = C;

   const int nqp = ir.GetNPoints();
   std::vector<real_t> Q_at_qps(nqp);

   real_t max_val = std::numeric_limits<real_t>::lowest();
   real_t min_val = std::numeric_limits<real_t>::max();

   Vector shape(ndof);
   for (int q = 0; q < nqp; q++)
   {
      const IntegrationPoint &ip = ir.IntPoint(q);
      ftr.SetAllIntPoints(&ip);

      IntegrationPoint ip_elem;
      ftr.Loc1.Transform(ip, ip_elem);

      fe.CalcShape(ip_elem, shape);

      // Per-QP interpolation: Q_self(q) = sum_i shape_i(q) * Q_dofs[i].
      // Deliberately written as a naive sequential accumulator so the
      // round-off behaviour matches wave_operator.inl's fault-QP loop
      // (which also does sequential accumulation across element DOFs,
      // see fault_face_flux.cpp:39-65 ComputeTrialTraction).
      real_t sum = 0.0;
      for (int i = 0; i < ndof; i++)
      {
         sum += shape(i) * Q_dofs(i);
      }
      Q_at_qps[q] = sum;
      if (sum > max_val) { max_val = sum; }
      if (sum < min_val) { min_val = sum; }
   }

   SweepResult r;
   r.C = C;
   r.min_Q = min_val;
   r.max_Q = max_val;
   r.abs_dispersion = max_val - min_val;
   r.rel_dispersion = (C != 0.0) ? (r.abs_dispersion / std::abs(C)) : 0.0;

   std::cout << "    C = " << std::scientific << std::setprecision(3) << C
             << "  min = " << std::setprecision(15) << min_val
             << "  max = " << max_val
             << std::scientific << std::setprecision(3)
             << "  |Δ| = " << r.abs_dispersion
             << "  rel = " << r.rel_dispersion << "\n";
   return r;
}

} // anonymous namespace

int main()
{
   std::cout << "\n=== TPV102 pepper reproducer (R-001 v9.4.0 gate) ===\n";

   // ---- Fixture: single hex-split-into-tets, order-1 L2 GaussLobatto ----
   // Use Cartesian3D with Element::TETRAHEDRON: the hex is split into
   // tets internally, producing tet elements that match TPV102 mesh
   // topology at the element level (face = triangle, 4 vertex DOFs per
   // tet, face QPs at Dunavant points).
   Mesh mesh = Mesh::MakeCartesian3D(1, 1, 1, Element::TETRAHEDRON,
                                     1.0, 1.0, 1.0);
   const int order = 1;
   L2_FECollection fec(order, 3, BasisType::GaussLobatto);
   FiniteElementSpace fes(&mesh, &fec);

   const FiniteElement *fe = fes.GetFE(0);
   const int ndof_per_el = fe->GetDof();
   const Geometry::Type face_geom = Geometry::TRIANGLE;
   const IntegrationRule &ir = IntRules.Get(face_geom, 2 * order);
   const int nqp = ir.GetNPoints();

   std::cout << "Fixture:\n"
             << "  element type   : TETRAHEDRON\n"
             << "  order          : " << order << "\n"
             << "  ndof_per_el    : " << ndof_per_el << "\n"
             << "  face QP count  : " << nqp
             << " (Dunavant " << nqp << "-point on triangle)\n"
             << "  basis          : L2 + GaussLobatto\n\n";

   // MFEM's hex→tets gives an element whose face 0 is a triangle face.
   // Pick any face — all faces expose the same per-QP interpolation
   // arithmetic.  We need an interior face (shared between two tets)
   // to get a FaceElementTransformations with Loc1.Transform populated.
   int chosen_face = -1;
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      FaceElementTransformations *ftr_probe =
         mesh.GetFaceElementTransformations(f);
      if (ftr_probe && ftr_probe->Elem2No >= 0 &&
          ftr_probe->GetGeometryType() == face_geom)
      {
         chosen_face = f;
         break;
      }
   }
   MFEM_VERIFY(chosen_face >= 0,
               "pepper reproducer: no interior triangle face found in the "
               "Cartesian3D-tet fixture");

   FaceElementTransformations *ftr =
      mesh.GetFaceElementTransformations(chosen_face);
   MFEM_VERIFY(ftr, "null FaceElementTransformations on chosen face");

   std::cout << "Sweep per-QP reconstruction Q_self(q) = Σ_i shape_i(q)·C:\n";

   // ---- C-sweep: measure dispersion at 4 magnitudes ----
   const std::vector<real_t> C_values = {1.0e6, 1.0e7, 1.0e8, 1.0e9};
   std::vector<SweepResult> results;
   results.reserve(C_values.size());
   for (real_t C : C_values)
   {
      results.push_back(ProbeOneMagnitude(*fe, *ftr, ir, C));
   }

   // ---- (a) non-zero dispersion exists ----
   std::cout << "\nSubtest (a): per-QP dispersion exists on at least one C\n";
   bool any_nonzero = false;
   for (const auto &r : results)
   {
      if (r.abs_dispersion > 0.0) { any_nonzero = true; break; }
   }
   TEST_ASSERT(any_nonzero,
               "(a) at least one C in the sweep shows nonzero per-QP "
               "dispersion (the ULP mechanism is present on this platform)");

   // If (a) fails, (b) and (c) are not meaningful — emit a diagnostic
   // hint so the reviewer knows the refactor is not justified by this
   // reproducer.
   if (!any_nonzero)
   {
      std::cout
         << "\n  [R-001 GATE DIAGNOSTIC] All four magnitudes produced zero\n"
         << "  per-QP dispersion.  The v9.4.0 ULP-seed hypothesis is NOT\n"
         << "  supported by this platform/compiler.  Either:\n"
         << "    (i)  shape-function partition-of-unity is exact at the\n"
         << "         Dunavant face QPs on this tet topology (GaussLobatto\n"
         << "         may round to 1.0 exactly at low order), OR\n"
         << "    (ii) the compiler is folding the sum_i shape(i)*C loop\n"
         << "         before rounding can differ across QPs.\n"
         << "  Do NOT proceed with the v9.4.0 split-prestress refactor\n"
         << "  without a different reproducer that captures the pepper\n"
         << "  seed through the ACTUAL code path exercised on Frontera.\n";
   }

   // ---- (b) dispersion is ULP-bounded for every C (only meaningful if (a)) ----
   std::cout << "\nSubtest (b): |Δ| < 100 * ULP(C) for every C\n";
   for (const auto &r : results)
   {
      const real_t ulp_bound =
         100.0 * std::numeric_limits<real_t>::epsilon() * std::abs(r.C);
      char buf[256];
      std::snprintf(buf, sizeof(buf),
         "(b) dispersion at C=%.0e bounded by 100·ULP(C)=%.3e (got %.3e)",
         r.C, ulp_bound, r.abs_dispersion);
      TEST_ASSERT(r.abs_dispersion <= ulp_bound, buf);
   }

   // ---- (c) dispersion scales linearly with C ----
   std::cout << "\nSubtest (c): relative dispersion roughly constant across C\n";
   if (any_nonzero)
   {
      real_t max_rel = 0.0;
      real_t min_rel = std::numeric_limits<real_t>::max();
      for (const auto &r : results)
      {
         if (r.abs_dispersion > 0.0)
         {
            if (r.rel_dispersion > max_rel) { max_rel = r.rel_dispersion; }
            if (r.rel_dispersion < min_rel) { min_rel = r.rel_dispersion; }
         }
      }
      const real_t ratio =
         (min_rel > 0.0) ? (max_rel / min_rel)
                         : std::numeric_limits<real_t>::infinity();
      std::cout << "    max_rel = " << std::scientific << std::setprecision(3)
                << max_rel
                << "  min_rel = " << min_rel
                << "  ratio   = " << ratio << "\n";
      // Allow up to 10× spread — strictly linear scaling would give
      // ratio = 1, but small C values may have better cancellation if
      // some floating-point accidents zero out the last-place error.
      TEST_ASSERT(ratio < 10.0,
                  "(c) max/min relative dispersion ratio < 10 "
                  "(dispersion scales linearly with C)");
   }
   else
   {
      std::cout << "    skipped — (a) failed, (c) not meaningful\n";
   }

   // ---- Quantitative seed-amplitude comparison TotalQ vs FluctuationQ ----
   // The plan's whole argument hinges on: dispersion(C=1.2e8 Pa) >>
   // dispersion(C=1e6 Pa).  Print that ratio explicitly so the reviewer
   // can read it off.
   std::cout
      << "\n  [R-001 QUANTITATIVE TAKEAWAY]\n"
      << "    Total-Q  bulk-Q scale  C ~ 1.2e8 Pa  →  |Δ| ≈ "
      << std::scientific << std::setprecision(3) << results.back().abs_dispersion
      << " Pa\n"
      << "    Fluct.-Q bulk-Q scale  C ~ 1.0e6 Pa  →  |Δ| ≈ "
      << results.front().abs_dispersion << " Pa\n";
   if (results.front().abs_dispersion > 0.0)
   {
      std::cout << "    Seed-amplitude reduction under split-prestress: ~"
                << (results.back().abs_dispersion /
                    results.front().abs_dispersion) << "×\n";
   }

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
