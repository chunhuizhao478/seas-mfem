// TPV102 pepper diagnostic — TRUE per-QP perturbation response.
//
// Per user direction (2026-04-22): use a SOLVED face interpolation
// system to inject a one-QP perturbation in the RECONSTRUCTED fault
// trace, then trace the response through three operator levels.
//
// Setup:
//   3 face-active DOFs (the tet vertices on the fault face) carry
//   shape weights shape_i(q) for q ∈ [0,nqp).  The face interpolation
//   matrix S(q, i) = shape_i(q) is square (3×3 for nqp=3).  Solve
//
//     S · delta_dof = eps * e_{q_t}
//
//   to obtain the volume-DOF perturbation that produces eps at QP q_t
//   on the fault face's reconstructed Q_self in ONE channel.
//
// Three diagnostic levels:
//   L1 — Interpolation: confirm Q_self(q) = eps * delta_{q,q_t} after
//        applying the prescribed delta_dof.
//   L2 — Fault-face operator: rotate Q_self/Q_nbr, run EvaluateTotal,
//        rotate Q_imp back, lift via shape1.  Measure per-DOF rhs
//        contribution from THIS face only.  Compare to baseline.
//   L3 — Full wave.Mult: end-to-end rhs response.  L3 − L2 isolates
//        the volume operator's contribution.
//
// Channel order: SYY (normal stress, cleanest), SXY (strike shear),
// SXZ (dip shear).

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
const char *kCompName[NUM_STATE] = {
   "SXX","SYY","SZZ","SXY","SYZ","SXZ","VX","VY","VZ"};

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
   std::cout << "\n=== TPV102 pepper: prescribed-QP perturbation response ===\n";
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

   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   if (int_faces.Size() == 0) { return 1; }
   const int face_idx = int_faces[0];
   auto *ftr = mesh.GetInteriorFaceTransformations(face_idx);
   const int e1 = ftr->Elem1No, e2 = ftr->Elem2No;
   const FiniteElement *fe1 = fes.GetFE(e1);
   const FiniteElement *fe2 = fes.GetFE(e2);
   const int ndof = fe1->GetDof();
   const int dof_offset1 = e1 * ndof;
   const int dof_offset2 = e2 * ndof;
   const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2*kOrder);
   const int nqp = ir.GetNPoints();

   // Set up fault DOFData (uniform background, no nucleation).
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

   // Build face interpolation matrix S(q, i) = shape1_q(i) for face-active
   // DOFs only.  In our 2-tet fixture, DOFs 0,1,2 are face-active (the
   // 3 vertices of the fault triangle); DOF 3 has shape = 0 on the face.
   // Cache shape values per QP.
   std::vector<int> face_active_dofs;
   {
      // Determine face-active DOFs from shape1 at q=0 (any q works).
      const IntegrationPoint &ip0 = ir.IntPoint(0);
      ftr->SetAllIntPoints(&ip0);
      IntegrationPoint ip1; ftr->Loc1.Transform(ip0, ip1);
      Vector s1(ndof); fe1->CalcShape(ip1, s1);
      for (int i = 0; i < ndof; i++)
      { if (std::abs(s1(i)) > 1e-12) face_active_dofs.push_back(i); }
   }
   const int nactive = static_cast<int>(face_active_dofs.size());
   std::cout << "  Face-active DOFs (e1): [";
   for (int i = 0; i < nactive; i++)
   { std::cout << face_active_dofs[i] << (i<nactive-1?", ":""); }
   std::cout << "]  nqp=" << nqp << "\n";
   if (nactive != nqp)
   {
      std::cout << "  WARN: nactive != nqp; face S not square — solve needs LSQ\n";
   }

   // Build S(q, i_active).
   DenseMatrix S(nqp, nactive);
   for (int q = 0; q < nqp; q++)
   {
      const IntegrationPoint &ip = ir.IntPoint(q);
      ftr->SetAllIntPoints(&ip);
      IntegrationPoint ip1; ftr->Loc1.Transform(ip, ip1);
      Vector s1(ndof); fe1->CalcShape(ip1, s1);
      for (int j = 0; j < nactive; j++)
      { S(q, j) = s1(face_active_dofs[j]); }
   }
   std::cout << "  Face interpolation matrix S(q, dof_active):\n";
   for (int q = 0; q < nqp; q++)
   {
      std::cout << "    q=" << q << ": [";
      for (int j = 0; j < nactive; j++)
      { std::cout << " " << S(q, j); }
      std::cout << " ]\n";
   }

   // Pre-solve S^{-1} (3x3 inverse).
   DenseMatrix Sinv(nactive, nqp);
   {
      DenseMatrix Stmp(S);
      DenseMatrixInverse Sinv_solver(Stmp);
      Sinv_solver.GetInverseMatrix(Sinv);
   }

   // Helper: build a perturbation for QP q_t in channel comp with magnitude
   // eps, then return the full-size Vector (size = NUM_STATE * ndof_total)
   // containing only the perturbation (to be ADDED to baseline Q).
   auto build_perturbation =
      [&](int comp, int q_t, real_t eps, Vector &delta_Q) {
      delta_Q.SetSize(size); delta_Q = 0.0;
      Vector rhs(nqp); rhs = 0.0; rhs(q_t) = eps;
      Vector dof_pert(nactive);
      Sinv.Mult(rhs, dof_pert);
      for (int j = 0; j < nactive; j++)
      {
         const int i_local = face_active_dofs[j];
         delta_Q[comp * ndof_total + dof_offset1 + i_local] = dof_pert(j);
      }
   };

   // Verify L1 (interpolation) for SYY perturbation at q=0.
   auto verify_interpolation = [&](int comp, int q_t, real_t eps,
                                    const Vector &Q_perturbed) {
      const real_t *Q_data = Q_perturbed.GetData();
      std::cout << "  L1 verify (channel " << kCompName[comp]
                << ", target q=" << q_t << ", eps=" << eps << "):\n";
      for (int q = 0; q < nqp; q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         IntegrationPoint ip1; ftr->Loc1.Transform(ip, ip1);
         Vector s1(ndof); fe1->CalcShape(ip1, s1);
         real_t Q_self_pert = 0;
         for (int i = 0; i < ndof; i++)
         {
            Q_self_pert += s1(i) *
               Q_data[comp*ndof_total + dof_offset1 + i];
         }
         // Subtract baseline (uniform Q_bg interpolated to q): bulk_bg[comp].
         const real_t delta_Q_self = Q_self_pert - bulk_bg[comp];
         const real_t expected = (q == q_t) ? eps : 0.0;
         std::cout << "    q=" << q << ": delta Q_self[" << kCompName[comp]
                   << "] = " << delta_Q_self
                   << " (expected " << expected
                   << ", err = " << (delta_Q_self - expected) << ")\n";
      }
   };

   // L3 driver: full wave.Mult delta-rhs.
   auto run_full_mult = [&](int comp, int q_t, real_t eps) {
      Vector Q0(size), Q1(size);
      InitializeStateTotal(Q0, ndof_total, TPV102Params::sigma_n,
                           TPV102Params::tau_ini);
      Vector delta_Q;
      build_perturbation(comp, q_t, eps, delta_Q);
      add(Q0, 1.0, delta_Q, Q1);

      Vector k0(size), k1(size);
      wave.Mult(Q0, k0);
      wave.Mult(Q1, k1);
      Vector dk(size); subtract(k1, k0, dk);

      // L1 verification
      verify_interpolation(comp, q_t, eps, Q1);

      // Find max |dk| per channel
      real_t max_per_c[NUM_STATE] = {0};
      int max_dof_per_c[NUM_STATE]; for (int c=0;c<NUM_STATE;c++) max_dof_per_c[c]=-1;
      for (int i = 0; i < dk.Size(); i++)
      {
         const int c = i / ndof_total;
         const real_t a = std::abs(dk(i));
         if (a > max_per_c[c])
         { max_per_c[c] = a; max_dof_per_c[c] = i % ndof_total; }
      }
      std::cout << "  L3 full wave.Mult delta-rhs (max |dk| per channel):\n";
      for (int c = 0; c < NUM_STATE; c++)
      {
         std::cout << "    " << kCompName[c] << "[" << max_dof_per_c[c]
                   << "] = " << max_per_c[c]
                   << "  ratio/eps = " << (max_per_c[c] / eps) << "\n";
      }
   };

   // L2 driver: ONLY the fault-face operator (rotate + EvaluateTotal +
   // lift), not the volume term or other faces.  Mimics the wave_operator
   // interior-fault dispatch but on a single face, and accumulates the
   // per-DOF flux contribution into a local rhs vector.
   auto run_fault_face_only = [&](int comp, int q_t, real_t eps,
                                   real_t out_dk[NUM_STATE][16]) {
      // Build perturbed Q.
      Vector Q0(size), Q1(size);
      InitializeStateTotal(Q0, ndof_total, TPV102Params::sigma_n,
                           TPV102Params::tau_ini);
      Vector delta_Q; build_perturbation(comp, q_t, eps, delta_Q);
      add(Q0, 1.0, delta_Q, Q1);

      // Helper to assemble fault-face rhs contribution into local arrays.
      auto fault_rhs = [&](const Vector &Q,
                           real_t rhs_e1[NUM_STATE][16],
                           real_t rhs_e2[NUM_STATE][16]) {
         for (int c = 0; c < NUM_STATE; c++)
            for (int i = 0; i < 16; i++) { rhs_e1[c][i] = 0; rhs_e2[c][i] = 0; }
         const real_t *Q_data = Q.GetData();

         // Reconstruct canonical fault frame from FaultBasis q=0.
         // (Same as wave_operator's per-QP reconstruction.)
         const Array<int> &intf = wave.GetFaultInteriorFaces();
         (void)intf;
         // Use BuildFrame for a simple, reproducible canonical frame.
         const IntegrationPoint &ip0 = ir.IntPoint(0);
         ftr->SetAllIntPoints(&ip0);
         Vector n_raw(3); CalcOrtho(ftr->Face->Jacobian(), n_raw);
         n_raw /= n_raw.Norml2();
         real_t can_n[3] = {n_raw(0), n_raw(1), n_raw(2)};
         real_t can_t1[3], can_t2[3];
         GodunovFlux::BuildFrame(can_n, can_t1, can_t2);
         DenseMatrix T_can(NUM_STATE), Tinv_can(NUM_STATE);
         GodunovFlux::BuildRotation(can_n, can_t1, can_t2, T_can);
         GodunovFlux::BuildRotationInverse(can_n, can_t1, can_t2, Tinv_can);

         for (int q = 0; q < nqp; q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            ftr->SetAllIntPoints(&ip);
            real_t nl = n_raw.Norml2();  // already normalized; use Jacobian
            Vector n_q(3); CalcOrtho(ftr->Face->Jacobian(), n_q);
            nl = n_q.Norml2();
            const real_t w = ip.weight * nl;
            IntegrationPoint ip1, ip2;
            ftr->Loc1.Transform(ip, ip1);
            ftr->Loc2.Transform(ip, ip2);
            Vector s1(ndof), s2(ndof);
            fe1->CalcShape(ip1, s1);
            fe2->CalcShape(ip2, s2);
            real_t Q_self[NUM_STATE], Q_nbr[NUM_STATE];
            for (int c = 0; c < NUM_STATE; c++)
            {
               Q_self[c] = 0; Q_nbr[c] = 0;
               for (int i = 0; i < ndof; i++)
               {
                  Q_self[c] += s1(i) * Q_data[c*ndof_total + dof_offset1 + i];
                  Q_nbr[c]  += s2(i) * Q_data[c*ndof_total + dof_offset2 + i];
               }
            }
            // Rotate to canonical frame.
            real_t Q_self_can[NUM_STATE], Q_nbr_can[NUM_STATE];
            for (int c = 0; c < NUM_STATE; c++)
            {
               Q_self_can[c] = 0; Q_nbr_can[c] = 0;
               for (int k = 0; k < NUM_STATE; k++)
               {
                  Q_self_can[c] += Tinv_can(c,k) * Q_self[k];
                  Q_nbr_can[c]  += Tinv_can(c,k) * Q_nbr[k];
               }
            }
            // Friction.
            DOFData fdata = dof_data[q];
            real_t Q_imp_plus[NUM_STATE], Q_imp_minus[NUM_STATE];
            ff.EvaluateTotal(fdata, Q_self_can, Q_nbr_can,
                             Q_imp_plus, Q_imp_minus);
            // Rotate Q_imp back to global.
            real_t Q_imp_plus_g[NUM_STATE], Q_imp_minus_g[NUM_STATE];
            for (int c = 0; c < NUM_STATE; c++)
            {
               Q_imp_plus_g[c] = 0; Q_imp_minus_g[c] = 0;
               for (int k = 0; k < NUM_STATE; k++)
               {
                  Q_imp_plus_g [c] += T_can(c,k) * Q_imp_plus[k];
                  Q_imp_minus_g[c] += T_can(c,k) * Q_imp_minus[k];
               }
            }
            // Compute per-side flux F_h.
            real_t F_h_plus[NUM_STATE], F_h_minus[NUM_STATE];
            wave.GetFlux().Interior(can_n, Q_imp_plus_g,
                                    Q_imp_plus_g, F_h_plus);
            wave.GetFlux().Interior(can_n, Q_imp_minus_g,
                                    Q_imp_minus_g, F_h_minus);
            // Lift (assuming e1 is on minus side, like the test fixture).
            // Just record both sides; caller can interpret.
            for (int c = 0; c < NUM_STATE; c++)
            {
               for (int i = 0; i < ndof; i++)
               {
                  rhs_e1[c][i] += w * s1(i) * F_h_minus[c];
                  rhs_e2[c][i] += w * s2(i) * F_h_plus[c];
               }
            }
         }
      };
      real_t r0_e1[NUM_STATE][16], r0_e2[NUM_STATE][16];
      real_t r1_e1[NUM_STATE][16], r1_e2[NUM_STATE][16];
      fault_rhs(Q0, r0_e1, r0_e2);
      fault_rhs(Q1, r1_e1, r1_e2);
      std::cout << "  L2 fault-face-only delta-rhs (max |dk| per channel, e1):\n";
      for (int c = 0; c < NUM_STATE; c++)
      {
         real_t mx = 0; int mi = -1;
         for (int i = 0; i < ndof; i++)
         {
            real_t d = std::abs(r1_e1[c][i] - r0_e1[c][i]);
            if (d > mx) { mx = d; mi = i; }
         }
         out_dk[c][0] = mx;
         std::cout << "    " << kCompName[c] << "[" << mi << "] = " << mx
                   << "  ratio/eps = " << (mx / eps) << "\n";
      }
   };

   // ---- Run the diagnostic for SYY → SXY → SXZ ----
   const real_t kEps = 1.0e-6;  // probe size; large enough to be visible
                                 // above FP roundoff but small enough to
                                 // stay in linear regime
   real_t dk_dummy[NUM_STATE][16];

   for (int probe : {SYY, SXY, SXZ})
   {
      std::cout << "\n=== Probe channel: " << kCompName[probe]
                << ", eps = " << kEps << " Pa, target q=0 ===\n";
      run_fault_face_only(probe, 0, kEps, dk_dummy);
      run_full_mult(probe, 0, kEps);
   }

   std::cout << "\n========================================\n"
             << "  Three-level QP perturbation diagnostic complete\n"
             << "========================================\n";
   return 0;
}
