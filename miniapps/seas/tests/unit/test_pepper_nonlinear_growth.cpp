// TPV102 pepper diagnostic — nonlinear growth across rupture regimes.
//
// Per user direction (2026-04-22):
//   (a) Sharpened L1 with delta-only formulation: compute delta_Q_self,
//       delta_F_h, delta_rhs WITHOUT subtracting from a 1.2e8 Pa
//       background — eliminates the cancellation ambiguity.
//   (b) eps ladder (1, 1e-3, 1e-6, 1e-9 Pa) to find the FP linearity
//       floor of the wave operator's response.
//   (c) Multistep growth-rate map: run baseline vs baseline+δ across
//       three regimes (nucleation OFF / WEAK / STRONG), measure
//       γ_n = ||δ_{n+1}|| / ||δ_n||.  γ > 1 in the strong regime only
//       isolates the rupture-friction nonlinear amplifier.

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

real_t Norm2(const Vector &v) { return v.Norml2(); }

} // anonymous

int main()
{
   std::cout << "\n=== TPV102 pepper nonlinear-growth diagnostic ===\n";
   std::cout << std::scientific << std::setprecision(6);

   Mesh mesh = BuildTwoTetFaultMesh();
   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr = 3;
   bc.absorbing_attrs = {};
   real_t bulk_bg[NUM_STATE] = {0};
   bulk_bg[SYY] =  TPV102Params::sigma_n;
   bulk_bg[SXY] = -TPV102Params::tau_ini;

   const auto build_wave = [&](real_t tau2_nuc_amp,
                                std::vector<DOFData> *dof_data_out,
                                FaultFaceFlux *ff_out,
                                int *nqp_out) -> WaveOperator<Mesh>* {
      WaveOperator<Mesh> *wave = new WaveOperator<Mesh>(
         mesh, kOrder, TPV102Params::lambda, TPV102Params::mu,
         TPV102Params::rho, bc);
      wave->SetAbsorbingBackground(bulk_bg);
      const Array<int> &intf = wave->GetFaultInteriorFaces();
      auto *ftr = mesh.GetInteriorFaceTransformations(intf[0]);
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(),
                                                2*kOrder);
      const int nqp = ir.GetNPoints();
      *nqp_out = nqp;
      std::vector<Vector> coords(nqp, Vector(3));
      for (int q = 0; q < nqp; q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         ftr->Face->Transform(ip, coords[q]);
      }
      InitializeFaultDOFs(*dof_data_out, nqp, coords);
      ZeroDOFDataPreStressTotal(*dof_data_out, nqp);
      for (int q = 0; q < nqp; q++)
      { (*dof_data_out)[q].tau2_nuc = tau2_nuc_amp; }
      wave->SetFaultFlux(ff_out);
      wave->SetFaultDOFData(dof_data_out, nqp);
      return wave;
   };

   // -----------------------------------------------------------------
   // (a) Sharpened L1: delta-only computation.  S * delta_dof = eps*e_{q_t}.
   //     Compute delta_Q_self via shape * delta_dof directly (no baseline).
   // -----------------------------------------------------------------
   std::cout << "\n--- (a) Delta-only L1 verification ---\n";
   {
      FaultFaceFlux ff(TPV102Params::lambda, TPV102Params::mu,
                        TPV102Params::rho);
      std::vector<DOFData> dof_data;
      int nqp;
      auto *wave = build_wave(0.0, &dof_data, &ff, &nqp);
      const auto &fes = wave->GetFESpace();
      const Array<int> &int_faces = wave->GetFaultInteriorFaces();
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[0]);
      const FiniteElement *fe1 = fes.GetFE(ftr->Elem1No);
      const int ndof = fe1->GetDof();
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(),
                                                2*kOrder);
      // Build S(q, j_active) and S^{-1}.
      std::vector<int> active;
      {
         const IntegrationPoint &ip0 = ir.IntPoint(0);
         ftr->SetAllIntPoints(&ip0);
         IntegrationPoint ip1; ftr->Loc1.Transform(ip0, ip1);
         Vector s1(ndof); fe1->CalcShape(ip1, s1);
         for (int i = 0; i < ndof; i++)
         { if (std::abs(s1(i)) > 1e-12) active.push_back(i); }
      }
      const int nact = active.size();
      DenseMatrix S(nqp, nact), Sinv(nact, nqp);
      for (int q = 0; q < nqp; q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         IntegrationPoint ip1; ftr->Loc1.Transform(ip, ip1);
         Vector s1(ndof); fe1->CalcShape(ip1, s1);
         for (int j = 0; j < nact; j++) { S(q, j) = s1(active[j]); }
      }
      DenseMatrixInverse(S).GetInverseMatrix(Sinv);

      const real_t eps_ladder[] = {1.0, 1e-3, 1e-6, 1e-9};
      for (real_t eps : eps_ladder)
      {
         Vector rhs_eps(nqp); rhs_eps = 0.0; rhs_eps(0) = eps;
         Vector delta_dof(nact); Sinv.Mult(rhs_eps, delta_dof);
         std::cout << "  eps=" << eps << "  delta_Q_self per QP (delta-only):\n";
         for (int q = 0; q < nqp; q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            ftr->SetAllIntPoints(&ip);
            IntegrationPoint ip1; ftr->Loc1.Transform(ip, ip1);
            Vector s1(ndof); fe1->CalcShape(ip1, s1);
            real_t dQ = 0;
            for (int j = 0; j < nact; j++) { dQ += s1(active[j]) * delta_dof(j); }
            const real_t expected = (q == 0) ? eps : 0.0;
            std::cout << "    q=" << q << "  delta_Q_self = " << dQ
                      << "  expected = " << expected
                      << "  err/eps = " << ((dQ - expected) / eps) << "\n";
         }
      }
      delete wave;
   }

   // -----------------------------------------------------------------
   // (a-cont) eps ladder: response in full wave.Mult, normalized.
   //          Confirms linearity range and FP floor.
   // -----------------------------------------------------------------
   std::cout << "\n--- (a) eps-ladder linearity (full Mult, channel SYY) ---\n";
   {
      FaultFaceFlux ff(TPV102Params::lambda, TPV102Params::mu,
                        TPV102Params::rho);
      std::vector<DOFData> dof_data;
      int nqp;
      auto *wave = build_wave(0.0, &dof_data, &ff, &nqp);
      const auto &fes = wave->GetFESpace();
      const Array<int> &int_faces = wave->GetFaultInteriorFaces();
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[0]);
      const int e1 = ftr->Elem1No;
      const FiniteElement *fe1 = fes.GetFE(e1);
      const int ndof = fe1->GetDof();
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(),
                                                2*kOrder);
      const int ndof_total = fes.GetNDofs();
      const int size = wave->Height();
      const int dof_offset1 = e1 * ndof;
      // Active DOFs and Sinv (same as above)
      std::vector<int> active;
      {
         const IntegrationPoint &ip0 = ir.IntPoint(0);
         ftr->SetAllIntPoints(&ip0);
         IntegrationPoint ip1; ftr->Loc1.Transform(ip0, ip1);
         Vector s1(ndof); fe1->CalcShape(ip1, s1);
         for (int i = 0; i < ndof; i++)
         { if (std::abs(s1(i)) > 1e-12) active.push_back(i); }
      }
      const int nact = active.size();
      DenseMatrix S(nqp, nact), Sinv(nact, nqp);
      for (int q = 0; q < nqp; q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         IntegrationPoint ip1; ftr->Loc1.Transform(ip, ip1);
         Vector s1(ndof); fe1->CalcShape(ip1, s1);
         for (int j = 0; j < nact; j++) { S(q, j) = s1(active[j]); }
      }
      DenseMatrixInverse(S).GetInverseMatrix(Sinv);

      Vector Q0(size), k0(size);
      InitializeStateTotal(Q0, ndof_total, TPV102Params::sigma_n,
                           TPV102Params::tau_ini);
      wave->Mult(Q0, k0);
      std::cout << "  channel " << kCompName[SYY]
                << ", target q=0, ladder ratios = ||dk||_inf/eps:\n";
      const real_t eps_ladder[] = {1.0, 1e-3, 1e-6, 1e-9, 1e-12};
      for (real_t eps : eps_ladder)
      {
         Vector rhs_eps(nqp); rhs_eps = 0.0; rhs_eps(0) = eps;
         Vector delta_dof(nact); Sinv.Mult(rhs_eps, delta_dof);
         Vector Q1 = Q0;
         for (int j = 0; j < nact; j++)
         {
            Q1[SYY * ndof_total + dof_offset1 + active[j]] += delta_dof(j);
         }
         Vector k1(size); wave->Mult(Q1, k1);
         Vector dk(size); subtract(k1, k0, dk);
         real_t mxx = 0; for (int i = 0; i < size; i++) { mxx = std::max(mxx, std::abs(dk(i))); }
         std::cout << "    eps=" << eps
                   << "  max|dk| = " << mxx
                   << "  ratio/eps = " << (mxx / eps)
                   << (mxx / eps > 1e6 ? "  <-- LINEAR" : "  <-- FP-DOMINATED")
                   << "\n";
      }
      delete wave;
   }

   // -----------------------------------------------------------------
   // (b) Multistep growth-rate map across rupture regimes.
   //     Baseline vs baseline+δ.  δ injected ONCE at step 0 in SYY
   //     channel via the prescribed-QP mechanism.  Track ||δ_n|| of
   //     several quantities per step.
   // -----------------------------------------------------------------
   std::cout << "\n--- (b) Multistep growth across rupture regimes ---\n";
   const real_t kPertEps = 1.0;  // 1 Pa SYY perturbation, well above FP floor
   const int kNStep = 30;
   const real_t kDt = 5.0e-5;
   struct Regime { const char *name; real_t tau2_nuc; };
   Regime regimes[] = {
      {"OFF (no nucleation)", 0.0},
      {"WEAK (1 MPa)",        1.0e6},
      {"STRONG (25 MPa)",     TPV102Params::nuc_dtau},
   };

   for (const Regime &r : regimes)
   {
      std::cout << "\n  --- regime: " << r.name
                << "  (tau2_nuc = " << r.tau2_nuc << " Pa) ---\n";
      // Baseline run
      FaultFaceFlux ff_b(TPV102Params::lambda, TPV102Params::mu,
                          TPV102Params::rho);
      std::vector<DOFData> dof_b;
      int nqp_b;
      auto *wave_b = build_wave(r.tau2_nuc, &dof_b, &ff_b, &nqp_b);
      const auto &fes_b = wave_b->GetFESpace();
      const int ndof_total = fes_b.GetNDofs();
      const int size = wave_b->Height();
      Vector Q_b(size), Q_b_new(size);
      InitializeStateTotal(Q_b, ndof_total, TPV102Params::sigma_n,
                           TPV102Params::tau_ini);

      // Perturbed run (separate wave + dof_data instances).
      FaultFaceFlux ff_p(TPV102Params::lambda, TPV102Params::mu,
                          TPV102Params::rho);
      std::vector<DOFData> dof_p;
      int nqp_p;
      auto *wave_p = build_wave(r.tau2_nuc, &dof_p, &ff_p, &nqp_p);
      Vector Q_p(size), Q_p_new(size);
      InitializeStateTotal(Q_p, ndof_total, TPV102Params::sigma_n,
                           TPV102Params::tau_ini);

      // Inject perturbation into Q_p[SYY] at face DOFs via S^{-1}.
      const Array<int> &int_faces = wave_p->GetFaultInteriorFaces();
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[0]);
      const int e1 = ftr->Elem1No;
      const FiniteElement *fe1 = fes_b.GetFE(e1);
      const int ndof = fe1->GetDof();
      const int dof_offset1 = e1 * ndof;
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(),
                                                2*kOrder);
      std::vector<int> active;
      {
         const IntegrationPoint &ip0 = ir.IntPoint(0);
         ftr->SetAllIntPoints(&ip0);
         IntegrationPoint ip1; ftr->Loc1.Transform(ip0, ip1);
         Vector s1(ndof); fe1->CalcShape(ip1, s1);
         for (int i = 0; i < ndof; i++)
         { if (std::abs(s1(i)) > 1e-12) active.push_back(i); }
      }
      const int nact = active.size();
      DenseMatrix S(nqp_b, nact), Sinv(nact, nqp_b);
      for (int q = 0; q < nqp_b; q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         IntegrationPoint ip1; ftr->Loc1.Transform(ip, ip1);
         Vector s1(ndof); fe1->CalcShape(ip1, s1);
         for (int j = 0; j < nact; j++) { S(q, j) = s1(active[j]); }
      }
      DenseMatrixInverse(S).GetInverseMatrix(Sinv);
      Vector pert_qps(nqp_b); pert_qps = 0.0; pert_qps(0) = kPertEps;
      Vector pert_dof(nact); Sinv.Mult(pert_qps, pert_dof);
      for (int j = 0; j < nact; j++)
      {
         Q_p[SYY * ndof_total + dof_offset1 + active[j]] += pert_dof(j);
      }

      // Measure norms
      auto qp_norm = [&](const std::vector<DOFData> &a,
                         const std::vector<DOFData> &b,
                         std::vector<real_t> &out) {
         out.resize(4);
         real_t dV = 0, dT1 = 0, dT2 = 0, dSn = 0;
         for (size_t q = 0; q < a.size(); q++)
         {
            dV  = std::max(dV,  std::abs(a[q].slip_rate    - b[q].slip_rate));
            dT1 = std::max(dT1, std::abs(a[q].tau1_corr    - b[q].tau1_corr));
            dT2 = std::max(dT2, std::abs(a[q].tau2_corr    - b[q].tau2_corr));
            dSn = std::max(dSn, std::abs(a[q].sigma_n_corr - b[q].sigma_n_corr));
         }
         out[0]=dV; out[1]=dT1; out[2]=dT2; out[3]=dSn;
      };

      Vector dQ(size); subtract(Q_p, Q_b, dQ);
      const real_t dQ0_norm = Norm2(dQ);
      std::cout << "    initial ||δQ||_2 = " << dQ0_norm << "\n";
      std::cout << "    step  ||δQ||         ||δV||        ||δT1||"
                << "        ||δT2||        ||δSn||       gamma_dQ\n";
      real_t prev_dQ_norm = dQ0_norm;
      for (int step = 0; step < kNStep; step++)
      {
         wave_b->AdvanceADER(Q_b, kDt, 2, Q_b_new); Q_b.Swap(Q_b_new);
         wave_p->AdvanceADER(Q_p, kDt, 2, Q_p_new); Q_p.Swap(Q_p_new);
         subtract(Q_p, Q_b, dQ);
         const real_t dQn = Norm2(dQ);
         std::vector<real_t> df; qp_norm(dof_p, dof_b, df);
         const real_t gamma = (prev_dQ_norm > 0) ? (dQn / prev_dQ_norm) : 0;
         std::cout << "    " << std::setw(4) << step
                   << "  " << dQn
                   << "  " << df[0] << "  " << df[1]
                   << "  " << df[2] << "  " << df[3]
                   << "  γ=" << gamma << "\n";
         prev_dQ_norm = dQn;
      }
      delete wave_b;
      delete wave_p;
   }

   std::cout << "\n========================================\n"
             << "  Nonlinear-growth diagnostic complete\n"
             << "========================================\n";
   return 0;
}
