// Diagnostic-only test for the TPV102 per-tet pepper investigation.
//
// Per the user's correction (2026-04-22): we should NOT modify code to
// fix the bug; we need to LOCALIZE the bug via diagnostics first.
// Per-DOF amplification of FP-roundoff into visible per-tet pepper has
// been confirmed phenomenologically (test_adjacent_triangle_*).  This
// file isolates WHERE in the wave operator the amplification happens.
//
// Tests (all read-only — no code modifications):
//   T1: Per-QP T_can identity check — for a planar fault face, the
//       canonical-frame rotation matrices computed at every QP must be
//       BIT-IDENTICAL.  If they differ even at 1 ULP, the runtime
//       rebuild is FP-sensitive.
//   T2: shape1·Q identity check — Q_self at every QP must equal the
//       uniform Q_bg value when Q is uniform.  If they differ at FP
//       scale, the shape-weighted interpolation introduces per-QP
//       roundoff that compounds.
//   T3: Per-DOF amplification factor — set bulk Q to Q_bg PLUS a
//       hand-chosen 1e-10 Pa perturbation at ONE DOF.  Run a SINGLE
//       wave.Mult.  Measure the max per-DOF rhs response.  The ratio
//       (max |rhs|) / (1e-10) is the linear amplification factor of
//       the wave operator's spatial discretization.  If this is much
//       larger than O(c/h * dt^-1), the operator is non-physical.
//   T4: Per-QP F_h variation under per-DOF Q variation — hand-set
//       different per-DOF Q values and call flux_.Interior(can_n,
//       Q_self_q, Q_self_q, F_h) for each QP.  Measure how much F_h
//       differs across QPs.  Quantifies the per-QP F_h variability
//       per unit per-DOF Q variation.

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

const char *kCompName[NUM_STATE] = {"SXX","SYY","SZZ","SXY","SYZ","SXZ","VX","VY","VZ"};

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
   std::cout << "\n=== TPV102 pepper diagnostic chain ===\n";
   std::cout << std::scientific << std::setprecision(15);

   Mesh mesh = BuildTwoTetFaultMesh();
   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr = 3;
   bc.absorbing_attrs = {};
   real_t bulk_bg[NUM_STATE] = {0};
   bulk_bg[SYY] =  TPV102Params::sigma_n;
   bulk_bg[SXY] = -TPV102Params::tau_ini;

   WaveOperator<Mesh> wave(mesh, kOrder, TPV102Params::lambda,
                            TPV102Params::mu, TPV102Params::rho, bc);
   wave.SetAbsorbingBackground(bulk_bg);

   FaultFaceFlux ff(TPV102Params::lambda, TPV102Params::mu, TPV102Params::rho);
   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();
   const int size = wave.Height();

   // Setup minimal fault DOFs (no nucleation, no friction drive)
   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   if (int_faces.Size() == 0)
   {
      std::cout << "FAIL: no interior fault faces\n"; return 1;
   }
   const int face_idx = int_faces[0];
   auto *ftr = mesh.GetInteriorFaceTransformations(face_idx);
   const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2*kOrder);
   const int nqp = ir.GetNPoints();
   std::cout << "  Fixture: 2-tet fault face, nqp = " << nqp << "\n";

   std::vector<Vector> fault_coords(nqp, Vector(3));
   for (int q = 0; q < nqp; q++)
   {
      const IntegrationPoint &ip = ir.IntPoint(q);
      ftr->SetAllIntPoints(&ip);
      ftr->Face->Transform(ip, fault_coords[q]);
   }
   std::vector<DOFData> dof_data;
   InitializeFaultDOFs(dof_data, nqp, fault_coords);
   ZeroDOFDataPreStressTotal(dof_data, nqp);
   wave.SetFaultFlux(&ff);
   wave.SetFaultDOFData(&dof_data, nqp);

   // ============================================================
   // T1: Per-QP T_can identity check on this planar face
   // ============================================================
   std::cout << "\n-- T1: per-QP T_can bit-identity (planar fault) --\n";
   {
      DenseMatrix T_q[3];   // up to 3 QPs
      for (int q = 0; q < nqp; q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector n_raw(3); CalcOrtho(ftr->Face->Jacobian(), n_raw);
         const real_t nl = n_raw.Norml2();
         n_raw /= nl;
         // Use BuildFrame (the per-QP fault path uses ComputeOrientedFrame
         // which has the same geometry inputs).
         real_t n_arr[3] = {n_raw(0), n_raw(1), n_raw(2)};
         real_t t1[3], t2[3];
         GodunovFlux::BuildFrame(n_arr, t1, t2);
         T_q[q].SetSize(NUM_STATE, NUM_STATE);
         GodunovFlux::BuildRotation(n_arr, t1, t2, T_q[q]);
      }
      // Compare q=0 vs q=1 vs q=2 entry-wise
      real_t max_diff_01 = 0, max_diff_02 = 0, max_diff_12 = 0;
      int worst_i01=0,worst_j01=0,worst_i02=0,worst_j02=0,worst_i12=0,worst_j12=0;
      for (int i = 0; i < NUM_STATE; i++)
      {
         for (int j = 0; j < NUM_STATE; j++)
         {
            const real_t d01 = std::abs(T_q[0](i,j) - T_q[1](i,j));
            const real_t d02 = std::abs(T_q[0](i,j) - T_q[2](i,j));
            const real_t d12 = std::abs(T_q[1](i,j) - T_q[2](i,j));
            if (d01 > max_diff_01) { max_diff_01 = d01; worst_i01=i; worst_j01=j; }
            if (d02 > max_diff_02) { max_diff_02 = d02; worst_i02=i; worst_j02=j; }
            if (d12 > max_diff_12) { max_diff_12 = d12; worst_i12=i; worst_j12=j; }
         }
      }
      std::cout << "    max|T(q=0)-T(q=1)| = " << max_diff_01
                << " at (" << worst_i01 << "," << worst_j01 << ")\n"
                << "    max|T(q=0)-T(q=2)| = " << max_diff_02
                << " at (" << worst_i02 << "," << worst_j02 << ")\n"
                << "    max|T(q=1)-T(q=2)| = " << max_diff_12
                << " at (" << worst_i12 << "," << worst_j12 << ")\n";
      std::cout << "    Verdict: "
                << ((max_diff_01==0 && max_diff_02==0 && max_diff_12==0)
                    ? "BIT-IDENTICAL (rotation NOT a pepper source)"
                    : "DIFFER — rotation IS a per-QP source") << "\n";
   }

   // ============================================================
   // T2: shape1·Q identity check at uniform Q
   // ============================================================
   std::cout << "\n-- T2: shape1·Q on uniform Q_bg (per-QP Q_self) --\n";
   {
      Vector Q(size);
      InitializeStateTotal(Q, ndof_total, TPV102Params::sigma_n,
                           TPV102Params::tau_ini);
      const real_t *Q_data = Q.GetData();
      const FiniteElement *fe1 = fes.GetFE(ftr->Elem1No);
      const int ndof = fe1->GetDof();
      const int dof_offset1 = ftr->Elem1No * fe1->GetDof();
      for (int q = 0; q < nqp; q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         IntegrationPoint ip1; ftr->Loc1.Transform(ip, ip1);
         Vector s1(ndof); fe1->CalcShape(ip1, s1);
         real_t Q_self_SXY = 0, Q_self_SYY = 0;
         for (int i = 0; i < ndof; i++)
         {
            Q_self_SXY += s1(i) * Q_data[SXY*ndof_total + dof_offset1 + i];
            Q_self_SYY += s1(i) * Q_data[SYY*ndof_total + dof_offset1 + i];
         }
         std::cout << "    q=" << q << "  shape1=[";
         for (int i = 0; i < ndof; i++)
         { std::cout << s1(i) << (i<ndof-1?", ":""); }
         std::cout << "]  Q_self[SXY] = " << Q_self_SXY
                   << "  Q_self[SYY] = " << Q_self_SYY << "\n";
         std::cout << "    diff from Q_bg[SXY] = " << (Q_self_SXY - bulk_bg[SXY])
                   << "  diff from Q_bg[SYY] = " << (Q_self_SYY - bulk_bg[SYY])
                   << "\n";
      }
   }

   // ============================================================
   // T3: Per-DOF amplification factor — perturbation at ONE DOF
   // ============================================================
   std::cout << "\n-- T3: amplification factor (1e-10 Pa perturb at 1 DOF) --\n";
   {
      const real_t kPerturb = 1.0e-10;
      Vector Q(size);
      InitializeStateTotal(Q, ndof_total, TPV102Params::sigma_n,
                           TPV102Params::tau_ini);
      // Perturb Q[SXY, DOF 0 of element 0] by +kPerturb.
      const int target_dof = SXY * ndof_total + 0;
      Q[target_dof] += kPerturb;

      Vector k(size); wave.Mult(Q, k);
      // Subtract baseline (rhs from Q_bg without perturbation)
      Vector Q0(size);
      InitializeStateTotal(Q0, ndof_total, TPV102Params::sigma_n,
                           TPV102Params::tau_ini);
      Vector k0(size); wave.Mult(Q0, k0);
      // delta_rhs = k - k0
      real_t max_delta = 0; int max_c = -1, max_d = -1;
      for (int i = 0; i < k.Size(); i++)
      {
         const real_t d = std::abs(k(i) - k0(i));
         if (d > max_delta)
         { max_delta = d; max_c = i / ndof_total; max_d = i % ndof_total; }
      }
      std::cout << "    perturb input: Q[SXY, dof=0] += "
                << kPerturb << " Pa\n"
                << "    max |delta rhs| = " << max_delta
                << " Pa/s at " << kCompName[max_c]
                << "[dof " << max_d << "]\n"
                << "    amplification ratio = " << (max_delta / kPerturb)
                << " (per second)\n";
      std::cout << "    Per Newton's law of motion + DG, expected ~ "
                << "cs/h_face ~ " << (TPV102Params::cs / kL)
                << " (1/s) for SXY perturbation propagating at cs.\n";
   }

   // ============================================================
   // T4: Per-QP F_h variation when Q_self differs per QP
   // ============================================================
   std::cout << "\n-- T4: F_h variation across QPs for hand-varied Q_self --\n";
   {
      // Build a representative can_n (use q=0 of the fault face).
      const IntegrationPoint &ip0 = ir.IntPoint(0);
      ftr->SetAllIntPoints(&ip0);
      Vector n_raw(3); CalcOrtho(ftr->Face->Jacobian(), n_raw);
      n_raw /= n_raw.Norml2();
      real_t can_n[3] = {n_raw(0), n_raw(1), n_raw(2)};

      // Hand-vary Q across "QPs": one with bulk_bg, others with bulk_bg
      // perturbed by 1e-10 Pa in SXY, SXZ.
      real_t Qb[NUM_STATE]; for (int c=0;c<NUM_STATE;c++) Qb[c]=bulk_bg[c];
      real_t Qa[NUM_STATE]; for (int c=0;c<NUM_STATE;c++) Qa[c]=bulk_bg[c];
      Qa[SXY] += 1e-10; Qa[SXZ] += 1e-10;

      real_t F_h_b[NUM_STATE], F_h_a[NUM_STATE];
      ff.GetSolver();  // ensure not unused
      // Use wave's flux (FaultFaceFlux has no public Interior; use wave)
      const GodunovFlux &gf = wave.GetFlux();
      gf.Interior(can_n, Qb, Qb, F_h_b);
      gf.Interior(can_n, Qa, Qa, F_h_a);

      std::cout << "    Q_self_a - Q_self_b = +1e-10 Pa at SXY,SXZ\n";
      for (int c = 0; c < NUM_STATE; c++)
      {
         const real_t d = F_h_a[c] - F_h_b[c];
         if (std::abs(d) > 1e-15)
         {
            std::cout << "    delta F_h[" << kCompName[c]
                      << "] = " << d << "\n";
         }
      }
   }

   std::cout << "\n========================================\n"
             << "  Diagnostic dump complete\n"
             << "========================================\n";
   return 0;
}
