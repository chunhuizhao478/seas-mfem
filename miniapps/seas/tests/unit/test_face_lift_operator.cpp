// MFEM face-lift operator diagnostic for the TPV102 pepper bug.
//
// The wave operator's per-QP face flux deposition implements the
// face-lift operator L: face-QP-space → bulk-element-DOF-space via
//
//   rhs[dof_i] -= sum_q w_q * shape1_q(i) * F_h_q[c]
//
// L is a (ndof × nqp) matrix per element (per face), with entries
// L[i, q] = w_q * shape1_q(i).
//
// This test injects ONE-HOT perturbations (F_h = δ_qq0 * unit) per QP
// and records the per-DOF response.  That gives one column of L per QP.
// We then probe two pepper-relevant invariants:
//
//   (A) Uniform F_h (= 1 at every QP for a single component) lifts to
//       what per-DOF distribution?  By divergence theorem, uniform F_h
//       on a planar face SHOULD give a per-DOF distribution that
//       integrates exactly to F_h * face_area.  Per-DOF entries are
//       w_q * shape1_q(i) summed over q == ∫_face shape1(i) * F_h dS.
//
//   (B) Per-QP one-hot lifts: how do they differ for adjacent QPs of
//       the same face?  If the per-QP lift columns are mutually nearly
//       orthogonal, the lift faithfully separates QPs.  If they're
//       highly coupled (large off-diagonal entries in L^T L), per-QP
//       perturbations bleed into the same DOFs.
//
// All read-only — no code modifications outside tests/unit/.

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

namespace mfem { namespace seas { int g_seas_my_rank = 0; } }

namespace {
constexpr real_t kL = 1000.0;
constexpr int    kOrder = 1;

Mesh BuildTwoTetFaultMesh()
{
   Mesh mesh(3, 5, 2, 0);
   real_t verts[5][3] = {
      {0.0, 0.0, 0.0}, {kL, 0.0, 0.0}, {0.0, 0.0, kL},
      {0.0,  kL, 0.0}, {0.0, -kL, 0.0},
   };
   for (int v = 0; v < 5; v++) { mesh.AddVertex(verts[v]); }
   mesh.AddTet(0, 1, 2, 4, 1);
   mesh.AddTet(0, 1, 2, 3, 1);
   mesh.FinalizeTopology();
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      Array<int> fv; mesh.GetFaceVertices(f, fv);
      if (fv.Size() != 3) { continue; }
      auto *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr && ftr->Elem2No >= 0)
      {
         real_t cy = 0;
         for (int v = 0; v < fv.Size(); v++) { cy += mesh.GetVertex(fv[v])[1]; }
         cy /= fv.Size();
         if (std::abs(cy) < 1e-8)
         { mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 3); }
         continue;
      }
      mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 1);
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}
} // anonymous

int main()
{
   std::cout << "\n=== MFEM face-lift operator diagnostic ===\n";
   std::cout << std::scientific << std::setprecision(15);

   Mesh mesh = BuildTwoTetFaultMesh();
   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr = 3;
   bc.absorbing_attrs = {};
   WaveOperator<Mesh> wave(mesh, kOrder, TPV102Params::lambda,
                            TPV102Params::mu, TPV102Params::rho, bc);
   const auto &fes = wave.GetFESpace();
   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   if (int_faces.Size() == 0) { return 1; }

   const int face_idx = int_faces[0];
   auto *ftr = mesh.GetInteriorFaceTransformations(face_idx);
   const int e1 = ftr->Elem1No;
   const int e2 = ftr->Elem2No;
   const FiniteElement *fe1 = fes.GetFE(e1);
   const FiniteElement *fe2 = fes.GetFE(e2);
   const int ndof = fe1->GetDof();
   const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2*kOrder);
   const int nqp = ir.GetNPoints();

   std::cout << "  Fixture: 2-tet, fault face " << face_idx
             << ", e1=" << e1 << " e2=" << e2
             << ", ndof_per_el=" << ndof << ", nqp=" << nqp << "\n";

   // Build the lift matrix L_e1[i, q] = w_q * shape1_q(i)
   // and L_e2[i, q] = w_q * shape2_q(i).
   DenseMatrix L_e1(ndof, nqp), L_e2(ndof, nqp);
   L_e1 = 0.0; L_e2 = 0.0;
   std::vector<real_t> w_arr(nqp);
   std::vector<DenseMatrix> shape_e1(1, DenseMatrix(ndof, nqp));
   std::vector<DenseMatrix> shape_e2(1, DenseMatrix(ndof, nqp));
   shape_e1[0] = 0.0; shape_e2[0] = 0.0;

   for (int q = 0; q < nqp; q++)
   {
      const IntegrationPoint &ip = ir.IntPoint(q);
      ftr->SetAllIntPoints(&ip);
      Vector n_raw(3); CalcOrtho(ftr->Face->Jacobian(), n_raw);
      const real_t nl = n_raw.Norml2();
      const real_t w = ip.weight * nl;
      w_arr[q] = w;
      IntegrationPoint ip1, ip2;
      ftr->Loc1.Transform(ip, ip1);
      ftr->Loc2.Transform(ip, ip2);
      Vector s1(ndof), s2(ndof);
      fe1->CalcShape(ip1, s1);
      fe2->CalcShape(ip2, s2);
      for (int i = 0; i < ndof; i++)
      {
         L_e1(i, q) = w * s1(i);
         L_e2(i, q) = w * s2(i);
         shape_e1[0](i, q) = s1(i);
         shape_e2[0](i, q) = s2(i);
      }
   }

   std::cout << "\n-- shape1 evaluated at QPs (rows = DOF, cols = QP) --\n";
   for (int i = 0; i < ndof; i++)
   {
      std::cout << "    DOF " << i << ":";
      for (int q = 0; q < nqp; q++)
      { std::cout << "  " << shape_e1[0](i, q); }
      std::cout << "\n";
   }
   std::cout << "  per-QP weights w_q:";
   for (int q = 0; q < nqp; q++) { std::cout << "  " << w_arr[q]; }
   std::cout << "\n";

   // ============================================================
   // (A) Uniform F_h check: lift unit F_h at every QP
   // ============================================================
   std::cout << "\n-- (A) Uniform F_h = 1 at every QP, lifted to DOFs --\n";
   {
      Vector F_h_uniform(nqp); F_h_uniform = 1.0;
      Vector dof_lift_e1(ndof);
      L_e1.Mult(F_h_uniform, dof_lift_e1);
      std::cout << "  per-DOF lift on e1 (sum over q of w_q*shape1_q(i)):\n";
      real_t sum = 0;
      for (int i = 0; i < ndof; i++)
      {
         std::cout << "    DOF " << i << ": " << dof_lift_e1(i) << "\n";
         sum += dof_lift_e1(i);
      }
      // Theoretical: sum_i ∫ shape(i) dS = ∫ 1 dS = face_area
      std::cout << "  sum(per-DOF) = " << sum << "  ; face_area expected = "
                << w_arr[0]*nqp << " (= sum_q w_q for uniform shape weights)\n";
      // For TET face quadrature with nqp=3 (Dunavant order 2):
      // each w_q = (1/3) * face_area, so sum_q w_q = face_area.
   }

   // ============================================================
   // (B) Per-QP one-hot lifts → columns of L
   // ============================================================
   std::cout << "\n-- (B) One-hot QP perturbations → per-DOF response --\n";
   for (int q0 = 0; q0 < nqp; q0++)
   {
      Vector F_h_onehot(nqp); F_h_onehot = 0.0; F_h_onehot(q0) = 1.0;
      Vector lift(ndof); L_e1.Mult(F_h_onehot, lift);
      std::cout << "  inject F_h[q=" << q0 << "] = 1 → DOF lift:";
      for (int i = 0; i < ndof; i++) { std::cout << "  " << lift(i); }
      std::cout << "\n";
   }

   // ============================================================
   // (C) L^T L (Gram matrix of lift columns).  Off-diagonal entries
   //     measure how strongly different QPs couple to the same DOFs.
   // ============================================================
   std::cout << "\n-- (C) Gram matrix L^T L (qp × qp) --\n";
   DenseMatrix LTL(nqp, nqp);
   MultAtB(L_e1, L_e1, LTL);
   for (int q1 = 0; q1 < nqp; q1++)
   {
      std::cout << "    [";
      for (int q2 = 0; q2 < nqp; q2++)
      {
         std::cout << " " << LTL(q1, q2);
      }
      std::cout << " ]\n";
   }
   // Off-diagonal / diagonal ratio quantifies QP coupling.
   real_t max_off = 0, min_diag = LTL(0,0), max_diag = LTL(0,0);
   for (int q1 = 0; q1 < nqp; q1++)
   {
      for (int q2 = 0; q2 < nqp; q2++)
      {
         if (q1 == q2)
         {
            min_diag = std::min(min_diag, LTL(q1,q2));
            max_diag = std::max(max_diag, LTL(q1,q2));
         }
         else
         {
            max_off = std::max(max_off, std::abs(LTL(q1,q2)));
         }
      }
   }
   std::cout << "    max|off-diag| = " << max_off
             << ",  diag in [" << min_diag << " ... " << max_diag << "]\n";
   std::cout << "    coupling ratio max_off/min_diag = "
             << (max_off / std::max(min_diag, 1e-30)) << "\n";

   // ============================================================
   // (D) L L^T (DOF × DOF Gram).  Off-diagonal measures DOF coupling
   //     induced by the face-lift: per-QP perturbations at the same
   //     QP can affect multiple DOFs.
   // ============================================================
   std::cout << "\n-- (D) Gram matrix L L^T (dof × dof) --\n";
   DenseMatrix LLT(ndof, ndof);
   MultABt(L_e1, L_e1, LLT);
   for (int i = 0; i < ndof; i++)
   {
      std::cout << "    [";
      for (int j = 0; j < ndof; j++)
      { std::cout << " " << LLT(i, j); }
      std::cout << " ]\n";
   }

   // ============================================================
   // (E) FP-rounding test: lift uniform F_h via different orderings
   //     of accumulation (sum-then-multiply vs multiply-then-sum).
   //     If results differ at the ULP, the lift's accumulation order
   //     introduces ULP-scale per-DOF noise — the SEED of pepper.
   // ============================================================
   std::cout << "\n-- (E) FP-order sensitivity of lift sum --\n";
   {
      // Order 1: rhs[i] = sum_q w_q * shape1_q(i) * F_h_q for fixed F_h_q=1
      // Order 2: rhs[i] = (sum_q w_q * shape1_q(i)) * 1.0 (mathematically equiv)
      Vector lift_o1(ndof); lift_o1 = 0.0;
      Vector lift_o2(ndof); lift_o2 = 0.0;
      // Order 1: pairwise w_q * shape1_q(i) * 1.0
      for (int i = 0; i < ndof; i++)
      {
         real_t s = 0;
         for (int q = 0; q < nqp; q++)
         { s += w_arr[q] * shape_e1[0](i, q) * 1.0; }
         lift_o1(i) = s;
      }
      // Order 2: pairwise w_q * shape1_q(i), then multiply by 1.0
      for (int i = 0; i < ndof; i++)
      {
         real_t s = 0;
         for (int q = 0; q < nqp; q++)
         { s += w_arr[q] * shape_e1[0](i, q); }
         lift_o2(i) = s * 1.0;
      }
      real_t max_diff = 0;
      for (int i = 0; i < ndof; i++)
      { max_diff = std::max(max_diff, std::abs(lift_o1(i) - lift_o2(i))); }
      std::cout << "    max|lift_order1 - lift_order2| = " << max_diff
                << "\n";

      // Order 3: REVERSE the QP iteration order
      Vector lift_rev(ndof); lift_rev = 0.0;
      for (int i = 0; i < ndof; i++)
      {
         real_t s = 0;
         for (int q = nqp - 1; q >= 0; q--)
         { s += w_arr[q] * shape_e1[0](i, q); }
         lift_rev(i) = s;
      }
      real_t max_diff_rev = 0;
      for (int i = 0; i < ndof; i++)
      {
         max_diff_rev = std::max(max_diff_rev,
                                 std::abs(lift_o1(i) - lift_rev(i)));
      }
      std::cout << "    max|lift_forward - lift_reverse| = " << max_diff_rev
                << "  (FP order matters at this scale)\n";
   }

   // ============================================================
   // (F) THE SMOKING-GUN PROBE: do TWO different shape evaluations
   //     of the SAME unit-Q field at different QPs give different
   //     Q_self values?  This was the T2 finding (1.49e-08 drift
   //     at q=0 only).
   // ============================================================
   std::cout << "\n-- (F) shape-weighted sum of UNIFORM input per QP --\n";
   {
      const real_t Q_unit = 1.20e8;  // sigma_n
      for (int q = 0; q < nqp; q++)
      {
         real_t s_forward = 0, s_reverse = 0;
         for (int i = 0; i < ndof; i++)
         { s_forward += shape_e1[0](i, q) * Q_unit; }
         for (int i = ndof - 1; i >= 0; i--)
         { s_reverse += shape_e1[0](i, q) * Q_unit; }
         std::cout << "    q=" << q
                   << ": forward sum = " << s_forward
                   << "  reverse sum = " << s_reverse
                   << "  diff = " << (s_forward - s_reverse) << "\n";
         std::cout << "         expected = " << Q_unit
                   << "  drift_fwd = " << (s_forward - Q_unit) << "\n";
      }
   }

   std::cout << "\n========================================\n"
             << "  Lift diagnostic complete\n"
             << "========================================\n";
   return 0;
}
