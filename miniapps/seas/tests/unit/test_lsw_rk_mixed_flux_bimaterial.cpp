// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// tests/unit/test_lsw_rk_mixed_flux_bimaterial.cpp
//
// Phase 4 of PLAN_mixed_flux_seissol_port_2026-06-19.md — the numerics gate for
// the SAFS v3-ALT production path: matrix (BimaterialWaveOperator, heterogeneous
// CVM material) + mixed_flux=Adjacent + RK4 + LSW.
//
// The port verification (workflow wf_ad2d4057) REFUTED the earlier "matrix+LSW+RK
// needs new C++" claim: AdvanceRKCoupledLSW_Spatial is templated on MeshType and
// takes the BASE WaveOperator<MeshType>&; BimaterialWaveOperator does not override
// Mult; the driver binds the bimaterial operator through a base reference
// (spatial_dyn_driver.cpp:1288).  This test upgrades that STRUCTURAL confidence to
// a RUNTIME guarantee BEFORE any Frontera allocation: it runs the exact production
// combination on a tiny heterogeneous serial mesh and asserts no abort + finite +
// bounded energy + the LSW arm actually fired through the matrix operator.
//
// This is ADDITIVE — it touches no production header, no Mult, no stepper body, so
// the TPV/BP5 byte-exact contract is preserved by construction.
//
// Helpers (BuildSmallFaultTetMesh, MakeLSWPatchDOF, SetupLSWFault, AllFinite) are
// the proven idioms from the scalar sibling test_lsw_rk_mixed_flux.cpp; the only
// new ingredient is the BimaterialWaveOperator construction + base-reference call.
//
// Usage:  ./seas_test_lsw_rk_mixed_flux_bimaterial   (serial, < 5 s)

#include "mfem.hpp"

#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/bimaterial_wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/godunov_flux_bimaterial.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/tpv205_friction.hpp"
#include "../../dynamic/tpv102_setup.hpp"
#include "../../dynamic/rk_time_stepper.hpp"
#include "../../dynamic/nucleation_method.hpp"
#include "../../dynamic/heterogeneous_material.hpp"
#include "../../config/tpv205_params.hpp"
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

#define TEST_TRUE(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else      { num_failed++; std::cout << "  FAILED [" << __LINE__ \
                       << "]: " << msg << "\n"; } \
} while (0)

namespace
{

// TPV205 nucleation-patch LSW DOFData (verbatim idiom from the scalar sibling).
DOFData MakeLSWPatchDOF()
{
   DOFData d;
   d.Zp_plus  = TPV205Params::Zp;
   d.Zp_minus = TPV205Params::Zp;
   d.Zs_plus  = TPV205Params::Zs;
   d.Zs_minus = TPV205Params::Zs;
   d.eta_p    = TPV205Params::eta_p;
   d.eta_s    = TPV205Params::eta_s;

   d.sigma_n0    = TPV205Params::sigma_n;     // 120 MPa
   d.tau1_0      = 0.0;                        // pure strike-slip
   d.tau2_0      = TPV205Params::tau_nuc;      // 81.6 MPa
   d.sigma_n_nuc = 0.0;
   d.tau1_nuc    = 0.0;
   d.tau2_nuc    = 0.0;

   d.lsw_mu_s = TPV205Params::mu_s;
   d.lsw_mu_d = TPV205Params::mu_d;
   d.lsw_d_c  = TPV205Params::d_c;
   d.a = 0.0; d.psi = 0.0; d.Dc = 0.0;

   d.slip_rate = 0.0; d.V1 = 0.0; d.V2 = 0.0;
   d.slip1 = 0.0; d.slip2 = 0.0;

   d.tau1_corr    = 0.0;
   d.tau2_corr    = TPV205Params::tau_nuc;
   d.sigma_n_corr = TPV205Params::sigma_n;
   return d;
}

// 48-tet 3x3x3 vertex grid with an interior fault plane at y = L/2 (fault attr 3,
// free attr 1).  Verbatim from the scalar sibling / test_mixed_flux_dispatch_*.
Mesh BuildSmallFaultTetMesh(real_t L = 1000.0)
{
   const int Nv = 3 * 3 * 3;
   Mesh mesh(3, Nv, 0, 0);
   auto vid = [&](int i, int j, int k) { return i + 3*j + 9*k; };
   for (int k = 0; k < 3; k++)
      for (int j = 0; j < 3; j++)
         for (int i = 0; i < 3; i++)
         {
            real_t v[3] = {i*0.5*L, j*0.5*L, k*0.5*L};
            mesh.AddVertex(v);
         }
   auto vfunc = [&](int ci, int cj, int ck, int lv) {
      return vid(ci + (lv & 1), cj + ((lv >> 1) & 1), ck + ((lv >> 2) & 1));
   };
   const int pat[6][4] = {
      {0, 1, 3, 7}, {0, 3, 2, 7}, {0, 2, 6, 7},
      {0, 6, 4, 7}, {0, 4, 5, 7}, {0, 5, 1, 7}
   };
   for (int ck = 0; ck < 2; ck++)
      for (int cj = 0; cj < 2; cj++)
         for (int ci = 0; ci < 2; ci++)
            for (int t = 0; t < 6; t++)
            {
               mesh.AddTet(vfunc(ci, cj, ck, pat[t][0]),
                           vfunc(ci, cj, ck, pat[t][1]),
                           vfunc(ci, cj, ck, pat[t][2]),
                           vfunc(ci, cj, ck, pat[t][3]), 1);
            }
   mesh.FinalizeTopology();
   const real_t eps = 1e-6;
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      Array<int> fv; mesh.GetFaceVertices(f, fv);
      if (fv.Size() != 3) { continue; }
      auto *ftr = mesh.GetFaceElementTransformations(f);
      real_t cy = 0.0;
      for (int v = 0; v < fv.Size(); v++) { cy += mesh.GetVertex(fv[v])[1]; }
      cy /= fv.Size();
      if (ftr && ftr->Elem2No >= 0)
      {
         if (std::abs(cy - 0.5 * L) < eps)
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

// fault_coords + InitializeFaultDOFs, then overwrite every QP to a TPV205 LSW
// nucleation patch.  `ff` must OUTLIVE every Mult/stepper call (raw ptr held).
// Takes the BASE WaveOperator<Mesh>& so a BimaterialWaveOperator binds here via
// the production base reference (spatial_dyn_driver.cpp:1288).
void SetupLSWFault(WaveOperator<Mesh> &wave, FaultFaceFlux &ff,
                   std::vector<DOFData> &dof_data,
                   const Mesh &mesh, int order)
{
   const int n_local_qps = wave.GetNumLocalFaultQPs();
   const int nbf         = wave.GetNbfPerFace();
   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   std::vector<Vector> fault_coords;
   for (int i = 0; i < int_faces.Size(); i++)
   {
      auto *ftr =
         const_cast<Mesh &>(mesh).GetInteriorFaceTransformations(int_faces[i]);
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2 * order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3); ftr->Face->Transform(ip, phys);
         fault_coords.push_back(phys);
      }
   }
   InitializeFaultDOFs(dof_data, n_local_qps, fault_coords);
   for (auto &d : dof_data)
   {
      const DOFData lsw = MakeLSWPatchDOF();
      d.sigma_n0 = lsw.sigma_n0; d.tau1_0 = lsw.tau1_0; d.tau2_0 = lsw.tau2_0;
      d.sigma_n_nuc = 0.0; d.tau1_nuc = 0.0; d.tau2_nuc = 0.0;
      d.lsw_mu_s = lsw.lsw_mu_s; d.lsw_mu_d = lsw.lsw_mu_d; d.lsw_d_c = lsw.lsw_d_c;
      d.a = 0.0; d.psi = 0.0; d.Dc = 0.0;
      d.slip1 = 0.0; d.slip2 = 0.0; d.slip_rate = 0.0; d.V1 = 0.0; d.V2 = 0.0;
   }
   wave.SetFaultFlux(&ff);
   wave.SetFaultDOFData(&dof_data, nbf);
   real_t zero_bg[NUM_STATE] = {0};
   wave.SetAbsorbingBackground(zero_bg);
}

bool AllFinite(const Vector &v)
{
   for (int i = 0; i < v.Size(); i++) { if (!std::isfinite(v(i))) { return false; } }
   return true;
}

// ===========================================================================
// L8 — AdvanceRKCoupledLSW_Spatial on a BimaterialWaveOperator via base ref.
//      matrix (heterogeneous) + mixed_flux=Adjacent + RK4 + LSW.
// ===========================================================================
void L8_BimaterialMixedFluxRK()
{
   std::cout << "\n-- L8: matrix + Adjacent + RK4 + LSW (the production path) --\n";
   const int order = 1;
   const real_t L = 1000.0;
   Mesh mesh = BuildSmallFaultTetMesh(L);

   BoundaryConfig bc;
   bc.natural_attrs = {1}; bc.fault_attr = 3; bc.absorbing_attrs = {};

   // Heterogeneous (Mode::Coefficient) material with a contrast ACROSS the fault
   // (y-gradient): every element distinct => genuinely non-Constant, so the
   // matrix path is REQUIRED (the scalar path rejects non-Constant material).
   // rho constant.  This mirrors the production CVM-bulk + LSW-fault structure.
   auto lam_fn = [L](const Vector &x) -> real_t
                 { return 30.0e9 * (1.0 + 0.3 * (x(1) / L)); };
   auto mu_fn  = [L](const Vector &x) -> real_t
                 { return 30.0e9 * (1.0 + 0.3 * (x(1) / L)); };
   FunctionCoefficient lam_c(lam_fn), mu_c(mu_fn);
   ConstantCoefficient  rho_c(2670.0);

   BimaterialWaveOperator<Mesh> bimat(
      mesh, order, MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c), bc);
   bimat.SetFaultFrictionLaw(FaultFrictionLaw::LSW);
   // R-1408: set the mixed-flux mode AFTER the mesh is final.
   bimat.SetMixedFluxMode(MixedFluxMode::Adjacent);
   // Permit central flux under RK (the ADER path would abort in ComputeMaxDt).
   bimat.SetCflRkAware(true);

   TEST_TRUE(!bimat.GetCentralFluxFaceSet().empty(),
             "Adjacent mode populated a non-empty central_flux_face_set_ on the "
             "matrix operator (mixed flux is actually active)");

   // The PRODUCTION binding: drive the bimaterial operator through a base ref.
   WaveOperator<Mesh> &wave = bimat;
   TEST_TRUE(wave.GetFaultFrictionLaw() == FaultFrictionLaw::LSW,
             "base-ref GetFaultFrictionLaw() == LSW (matrix op, LSW arm selected)");

   FaultFaceFlux ff(TPV205Params::rho, TPV205Params::cp, TPV205Params::cs);
   std::vector<DOFData> dof_data;
   SetupLSWFault(wave, ff, dof_data, mesh, order);

   // Seed a modest Gaussian bulk stress pulse so the energy scale E0 is well
   // defined (a central-flux RK instability blows up >> E0; a stable scheme stays
   // O(E0) plus the modest fault-radiated growth).
   const int size = wave.Height();
   auto &fes = const_cast<FiniteElementSpace &>(wave.GetFESpace());
   const int ndof_total = fes.GetNDofs();
   const real_t sig = L / 6.0, amp = 1.0e6;
   FunctionCoefficient syy([sig, amp, L](const Vector &x) -> real_t {
      const real_t dx = x(0) - 0.5*L, dy = x(1) - 0.5*L, dz = x(2) - 0.5*L;
      return amp * std::exp(-(dx*dx + dy*dy + dz*dz) / (sig*sig)); });
   Vector Q(size); Q = 0.0;
   GridFunction gf(&fes); gf.ProjectCoefficient(syy);
   for (int d = 0; d < ndof_total; d++)
   { Q(SYY*ndof_total + d) = gf(d); Q(SXX*ndof_total + d) = 0.5*gf(d); }

   const real_t E0 = Q.Norml2();

   // March several RK4 steps via the base reference.  dt = 1e-4 s is the
   // CFL-safe step the scalar sibling uses on this mesh/material.
   const int n_steps = 20;
   const real_t dt = 1.0e-4;
   Vector Q_new(size);
   bool all_finite = true;
   real_t t = 0.0, E_max = E0;
   for (int s = 0; s < n_steps; s++)
   {
      AdvanceRKCoupledLSW_Spatial<Mesh>(wave, dof_data, Q, dt, t, Q_new,
                                        MakeRK4Tableau(), nullptr);
      if (!AllFinite(Q_new)) { all_finite = false; break; }
      Q = Q_new;
      t += dt;
      E_max = std::max(E_max, Q.Norml2());
   }

   const real_t E_final = Q.Norml2();
   std::cout << "    E0=" << std::scientific << std::setprecision(3) << E0
             << "  E_final=" << E_final << "  E_max=" << E_max << "\n";

   TEST_TRUE(all_finite,
             "20 RK4 steps via base ref: Q_new finite every step (no abort, no "
             "central-flux blow-up on the matrix operator)");
   TEST_TRUE(std::isfinite(E_final) && E_final > 0.0 && E_max < 1.0e4 * E0,
             "energy bounded: E stays O(E0) (< 1e4*E0), no exponential growth");

   bool any_slid = false, slip_finite = true;
   for (const auto &d : dof_data)
   {
      if (!std::isfinite(d.slip1) || !std::isfinite(d.slip2)) { slip_finite = false; }
      if (std::abs(d.slip2) > 1e-9) { any_slid = true; }
   }
   TEST_TRUE(slip_finite, "fault slip finite on every QP after the RK march");
   TEST_TRUE(any_slid,
             "fault slid (slip2 != 0) — EvaluateLSW fired through the matrix "
             "operator's virtual Mult under the base-ref RK stepper");
}

} // namespace

int main(int argc, char *argv[])
{
   (void)argc; (void)argv;
   std::cout << "\n=== matrix + RK + LSW + mixed-flux (bimaterial) — L8 ===\n";
   L8_BimaterialMixedFluxRK();
   std::cout << "\n========================================\n"
             << "  Results: " << num_passed << " passed, " << num_failed
             << " failed out of " << num_tests << " tests\n"
             << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
