// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// tests/unit/test_lsw_rk_mixed_flux.cpp — Phase 14 RK + linear slip-weakening
// (LSW) + scalar mixed flux.  Plan:
//   document/mixed_flux_dev/PLAN_rk45_lsw_mixed_flux_2026-05-29.md
//
// Serial unit coverage L1–L7 (the parallel shared-fault check L8 lives in
// test_lsw_rk_shared_fault_mpi.cpp):
//
//   L1  EvaluateLSW(Q±) == EvaluateADER_LSW(Q±·dt, dt) — bit-exact for
//       power-of-two dt (the I/dt round-trip is exact), rtol 1e-12 for a
//       general dt; Q_imp == I_imp/dt.
//   L2  EvaluateLSW is SLIP-STATELESS: data.slip1/slip2 unchanged by the call.
//   L3  Strength-barrier QP (lsw_mu_s = mu_s_barrier) ⇒ V = 0, τ_corr = τ_trial;
//       tensile σ_n (default floor 0) ⇒ free-slide V = |τ_total|/η_s.
//   L4  WaveOperator Mult dispatch on fault_friction_law_: LSW dof_data + LSW
//       flag ⇒ EvaluateLSW fires (finite, sliding); default RateAndState path
//       still runs (regression); GetFaultFrictionLaw round-trips.
//   L5  Bulk equivalence: AdvanceRKCoupledLSW_Spatial(RK4) == DoRK4Step on a
//       frictionless cube (the LSW stepper's bulk integration matches the proven
//       coupled-RK4 idiom; fault state absent).
//   L6  Coupled (Q, slip) RK4: AdvanceRKCoupledLSW_Spatial matches an independent
//       manual reference (Q AND slip) to round-off, and the ABSOLUTE slip combine
//       differs measurably from the buggy `+=` combine (guards the R-002 trap).
//   L7  FSAL predicate: classical RK4 is NOT FSAL (endpoint re-eval needed); DP45
//       IS FSAL (skipped) — the predicate AdvanceRKCoupledLSW_Spatial keys on.
//
// Usage:  ./seas_test_lsw_rk_mixed_flux   (serial, < 5 s)

#include "mfem.hpp"

#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/tpv205_friction.hpp"
#include "../../dynamic/tpv102_setup.hpp"
#include "../../dynamic/rk_time_stepper.hpp"
#include "../../dynamic/nucleation_method.hpp"
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

// Reference classical-RK4 bulk idiom (git:8461c67 / test_rk4_conservation.cpp).
// Forward-declared here so L5 (inside the anonymous namespace) can call it;
// defined at end of file.
void DoRK4Step_local(const WaveOperator<Mesh>& wave, Vector& Q, real_t dt);

namespace
{

// ---------------------------------------------------------------------------
// A TPV205 nucleation-patch LSW DOFData (rupture-area QP: tau_nuc = 81.6 MPa,
// sigma_n0 = 120 MPa, mu_s = 0.677, mu_d = 0.525, d_c = 0.40 m).  Mirrors
// MakeNucleationPatchDOF in test_tpv205_evaluate_ader_lsw_parity.cpp.  Slip
// defaults to 0 (μ = μ_s); the caller may set slip1/slip2 to exercise μ(δ).
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
   // Rate-and-state slots zeroed defensively (LSW must not read them).
   d.a = 0.0; d.psi = 0.0; d.Dc = 0.0;

   d.slip_rate = 0.0; d.V1 = 0.0; d.V2 = 0.0;
   d.slip1 = 0.0; d.slip2 = 0.0;

   d.tau1_corr    = 0.0;
   d.tau2_corr    = TPV205Params::tau_nuc;
   d.sigma_n_corr = TPV205Params::sigma_n;
   return d;
}

// A representative NONZERO fault-local bulk state (fluctuation), so
// ComputeTrialTraction produces a nonzero trial.  Values are arbitrary but
// fixed (deterministic).  Component order: SXX,SYY,SZZ,SXY,SYZ,SXZ,VX,VY,VZ.
void FillRepresentativeQ(real_t* Qp, real_t* Qm)
{
   for (int c = 0; c < NUM_STATE; ++c) { Qp[c] = 0.0; Qm[c] = 0.0; }
   Qp[SXX] =  1.0e6;  Qm[SXX] = -0.7e6;
   Qp[SXY] =  2.0e6;  Qm[SXY] = -1.3e6;
   Qp[SXZ] = -0.5e6;  Qm[SXZ] =  0.9e6;
   Qp[VX]  =  0.02;   Qm[VX]  = -0.015;
   Qp[VY]  = -0.03;   Qm[VY]  =  0.011;
   Qp[VZ]  =  0.007;  Qm[VZ]  = -0.004;
}

// 48-tet 3x3x3 vertex grid with an interior fault plane at y = L/2 (fault
// attr = 3, free = 1).  Verbatim from test_mixed_flux_dispatch_adjacent.cpp.
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

// Build fault_coords + initialize fault DOFData, then overwrite every QP to a
// TPV205 LSW nucleation-patch QP (lsw_* set; rate-state slots zeroed).  `ff`
// must OUTLIVE every wave.Mult / stepper call (the operator holds a raw ptr).
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
      // Keep the impedances InitializeFaultDOFs set (TPV102==TPV205 material),
      // overwrite the friction/prestress to the LSW nucleation patch.
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
// L1 — EvaluateLSW == EvaluateADER_LSW(·dt)/dt
// ===========================================================================
void L1_EvaluateLSW_vs_ADER()
{
   std::cout << "\n-- L1: EvaluateLSW(Q±) == EvaluateADER_LSW(Q±·dt, dt) --\n";
   FaultFaceFlux flux(TPV205Params::rho, TPV205Params::cp, TPV205Params::cs);

   real_t Qp[NUM_STATE], Qm[NUM_STATE];
   FillRepresentativeQ(Qp, Qm);

   // Exercise the slip-weakening regime: δ partway to d_c (μ between μ_s, μ_d).
   auto make = [&]() { DOFData d = MakeLSWPatchDOF();
                       d.slip1 = 0.05; d.slip2 = 0.12; return d; };

   const real_t pow2_dts[] = {1.0, 0.5, 0.25};
   bool all_pow2_exact = true;
   bool vimp_pow2_ok = true;  // (regression) EvaluateLSW must store v_imp (TPV6 stations)
   for (real_t dt : pow2_dts)
   {
      DOFData a = make(), b = make();
      real_t qip[NUM_STATE], qim[NUM_STATE];
      real_t Ip[NUM_STATE], Im[NUM_STATE], iip[NUM_STATE], iim[NUM_STATE];
      for (int c = 0; c < NUM_STATE; ++c) { Ip[c] = Qp[c] * dt; Im[c] = Qm[c] * dt; }
      flux.EvaluateLSW(a, Qp, Qm, qip, qim);
      flux.EvaluateADER_LSW(b, Ip, Im, dt, iip, iim);
      bool fields_exact =
         a.V1 == b.V1 && a.V2 == b.V2 && a.slip_rate == b.slip_rate &&
         a.tau1_corr == b.tau1_corr && a.tau2_corr == b.tau2_corr &&
         a.sigma_n_corr == b.sigma_n_corr;
      bool imp_exact = true;
      for (int c = 0; c < NUM_STATE; ++c)
      {
         if (qip[c] != iip[c] / dt || qim[c] != iim[c] / dt) { imp_exact = false; }
      }
      if (!(fields_exact && imp_exact)) { all_pow2_exact = false; }

      // (regression) EvaluateLSW must POPULATE DOFData::v_imp_{plus,minus} — the
      // TPV6/7 per-side station writer reads ONLY these.  They must be non-zero on
      // a sliding QP AND bit-exact equal to EvaluateADER_LSW's per-side store (the
      // ADER path was already correct).  Before the fix these stayed {0,0,0}, so the
      // station traces had zero velocity/displacement while stress was right.
      bool vimp_nonzero = false, vimp_match = true;
      for (int k = 0; k < 3; ++k)
      {
         if (a.v_imp_plus[k] != 0.0 || a.v_imp_minus[k] != 0.0) { vimp_nonzero = true; }
         if (a.v_imp_plus[k]  != b.v_imp_plus[k] ||
             a.v_imp_minus[k] != b.v_imp_minus[k]) { vimp_match = false; }
      }
      if (!(vimp_nonzero && vimp_match)) { vimp_pow2_ok = false; }
   }
   TEST_TRUE(all_pow2_exact,
             "bit-exact DOFData + Q_imp==I_imp/dt for dt in {1, 0.5, 0.25}");
   TEST_TRUE(vimp_pow2_ok,
             "EvaluateLSW populates DOFData::v_imp_{plus,minus} (non-zero, == "
             "EvaluateADER_LSW store) — TPV6/7 station zero-velocity regression");

   // General dt: rtol 1e-12 (the (Q·dt)·(1/dt) round-trip is ~1 ULP off).
   {
      const real_t dt = 0.1;
      DOFData a = make(), b = make();
      real_t qip[NUM_STATE], qim[NUM_STATE];
      real_t Ip[NUM_STATE], Im[NUM_STATE], iip[NUM_STATE], iim[NUM_STATE];
      for (int c = 0; c < NUM_STATE; ++c) { Ip[c] = Qp[c] * dt; Im[c] = Qm[c] * dt; }
      flux.EvaluateLSW(a, Qp, Qm, qip, qim);
      flux.EvaluateADER_LSW(b, Ip, Im, dt, iip, iim);
      auto rel = [](real_t x, real_t y) {
         return std::abs(x - y) / std::max({std::abs(x), std::abs(y), real_t(1)}); };
      bool ok = rel(a.V1, b.V1) < 1e-12 && rel(a.V2, b.V2) < 1e-12 &&
                rel(a.slip_rate, b.slip_rate) < 1e-12 &&
                rel(a.tau1_corr, b.tau1_corr) < 1e-12 &&
                rel(a.tau2_corr, b.tau2_corr) < 1e-12 &&
                rel(a.sigma_n_corr, b.sigma_n_corr) < 1e-12;
      TEST_TRUE(ok, "DOFData matches within rtol 1e-12 for dt = 0.1");
      TEST_TRUE(a.slip_rate > 0.0,
                "sliding QP exercised (slip_rate > 0, EvaluateLSW actually solved)");
   }
}

// ===========================================================================
// L2 — EvaluateLSW is slip-stateless
// ===========================================================================
void L2_SlipStateless()
{
   std::cout << "\n-- L2: EvaluateLSW leaves data.slip1/slip2 unchanged --\n";
   FaultFaceFlux flux(TPV205Params::rho, TPV205Params::cp, TPV205Params::cs);
   real_t Qp[NUM_STATE], Qm[NUM_STATE];
   FillRepresentativeQ(Qp, Qm);
   DOFData d = MakeLSWPatchDOF();
   d.slip1 = 0.037; d.slip2 = -0.091;
   const real_t s1 = d.slip1, s2 = d.slip2;
   real_t qip[NUM_STATE], qim[NUM_STATE];
   flux.EvaluateLSW(d, Qp, Qm, qip, qim);
   TEST_TRUE(d.slip1 == s1 && d.slip2 == s2,
             "data.slip1/slip2 bit-unchanged across EvaluateLSW");
   TEST_TRUE(d.slip_rate >= 0.0 && std::isfinite(d.slip_rate),
             "slip_rate written, finite, non-negative");
}

// ===========================================================================
// L3 — barrier lock + tensile free-slide
// ===========================================================================
void L3_BarrierAndTensile()
{
   std::cout << "\n-- L3: strength barrier (V=0) + tensile free-slide --\n";
   FaultFaceFlux flux(TPV205Params::rho, TPV205Params::cp, TPV205Params::cs);
   real_t Qp[NUM_STATE], Qm[NUM_STATE];
   FillRepresentativeQ(Qp, Qm);

   // Barrier QP: mu_s = mu_s_barrier ⇒ locked regardless of traction.
   {
      DOFData d = MakeLSWPatchDOF();
      d.lsw_mu_s = TPV205Params::mu_s_barrier;
      real_t qip[NUM_STATE], qim[NUM_STATE];
      flux.EvaluateLSW(d, Qp, Qm, qip, qim);
      TEST_TRUE(d.slip_rate == 0.0 && d.V1 == 0.0 && d.V2 == 0.0,
                "barrier QP locked: V = 0");
   }

   // Tensile sigma_n_total <= 0 with default floor (0) ⇒ strength 0 ⇒
   // free-slide V = |tau_total| / eta_s.  Use Q = 0 so total == prestress.
   {
      DOFData d = MakeLSWPatchDOF();
      d.sigma_n0 = -10.0e6;        // tensile background
      d.tau2_0   = 70.0e6;         // |tau_total| = 70 MPa (Q=0)
      real_t Qz[NUM_STATE] = {0};
      real_t qip[NUM_STATE], qim[NUM_STATE];
      flux.EvaluateLSW(d, Qz, Qz, qip, qim);
      const real_t expected = 70.0e6 / TPV205Params::eta_s;
      const real_t rel = std::abs(d.slip_rate - expected) / expected;
      TEST_TRUE(rel < 1e-12,
                "tensile sigma_n (floor 0) free-slides: V = |tau|/eta_s");
   }
}

// ===========================================================================
// L4 — Mult dispatch on fault_friction_law_
// ===========================================================================
void L4_MultDispatch()
{
   std::cout << "\n-- L4: WaveOperator::Mult dispatch (LSW vs RS default) --\n";
   const int order = 1;
   Mesh mesh = BuildSmallFaultTetMesh(1000.0);

   BoundaryConfig bc;
   bc.natural_attrs = {1}; bc.fault_attr = 3; bc.absorbing_attrs = {};

   // GetFaultFrictionLaw round-trip + default is RateAndState.
   {
      WaveOperator<Mesh> wave(mesh, order, TPV205Params::lambda,
                              TPV205Params::mu, TPV205Params::rho, bc);
      TEST_TRUE(wave.GetFaultFrictionLaw() == FaultFrictionLaw::RateAndState,
                "default fault_friction_law_ == RateAndState");
      wave.SetFaultFrictionLaw(FaultFrictionLaw::LSW);
      TEST_TRUE(wave.GetFaultFrictionLaw() == FaultFrictionLaw::LSW,
                "SetFaultFrictionLaw(LSW) round-trips via GetFaultFrictionLaw");
   }

   // LSW dispatch: LSW dof_data (a=psi=Dc=0) + LSW flag ⇒ EvaluateLSW fires.
   // If the RS Evaluate had (wrongly) fired on a=0 data, slip_rate would be
   // NaN; finite + sliding proves the LSW arm dispatched.
   {
      WaveOperator<Mesh> wave(mesh, order, TPV205Params::lambda,
                              TPV205Params::mu, TPV205Params::rho, bc);
      wave.SetFaultFrictionLaw(FaultFrictionLaw::LSW);
      FaultFaceFlux ff(TPV205Params::rho, TPV205Params::cp, TPV205Params::cs);
      std::vector<DOFData> dof_data;
      SetupLSWFault(wave, ff, dof_data, mesh, order);

      Vector Q(wave.Height()); Q = 0.0;     // at rest (fluctuation-Q)
      Vector dQdt(wave.Height());
      wave.Mult(Q, dQdt);

      bool slip_finite = true, any_slid = false;
      for (const auto &d : dof_data)
      {
         if (!std::isfinite(d.slip_rate) || d.slip_rate < 0.0) { slip_finite = false; }
         if (d.slip_rate > 1e-6) { any_slid = true; }
      }
      TEST_TRUE(AllFinite(dQdt), "LSW Mult: dQdt all finite (no NaN ⇒ not RS-on-a=0)");
      TEST_TRUE(slip_finite, "LSW Mult: every QP slip_rate finite & non-negative");
      TEST_TRUE(any_slid,
                "LSW Mult: nucleation patch slides (EvaluateLSW solved, V>0)");
   }

   // Default RateAndState path still runs (regression that Phase 2 left the
   // RS arm intact).  Standard TPV102 RS fault init.
   {
      WaveOperator<Mesh> wave(mesh, order, TPV102Params::lambda,
                              TPV102Params::mu, TPV102Params::rho, bc);
      // default law == RateAndState
      const int n_local_qps = wave.GetNumLocalFaultQPs();
      const int nbf = wave.GetNbfPerFace();
      const Array<int> &int_faces = wave.GetFaultInteriorFaces();
      std::vector<Vector> fault_coords;
      for (int i = 0; i < int_faces.Size(); i++)
      {
         auto *ftr =
            const_cast<Mesh &>(mesh).GetInteriorFaceTransformations(int_faces[i]);
         const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2*order);
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
      FaultFaceFlux ff(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
      wave.SetFaultFlux(&ff);
      wave.SetFaultDOFData(&dof_data, nbf);
      real_t zero_bg[NUM_STATE] = {0};
      wave.SetAbsorbingBackground(zero_bg);
      Vector Q(wave.Height()); Q = 0.0;
      Vector dQdt(wave.Height());
      wave.Mult(Q, dQdt);
      TEST_TRUE(AllFinite(dQdt),
                "RateAndState default Mult still runs (regression: dQdt finite)");
   }
}

// Reference coupled-RK4 on (Q, slip), independent of AdvanceRKCoupledLSW_Spatial:
// stage-local slip staging + ABSOLUTE combine.  Also returns the BUGGY (`+=`)
// slip the trap would produce.  dof_data is left in the post-step state.
struct LSWRef { Vector Q_new; std::vector<real_t> slip1_abs, slip2_abs,
                              slip1_bug, slip2_bug; };
LSWRef DoLSWRK4Reference(WaveOperator<Mesh> &wave, std::vector<DOFData> &dof,
                         const Vector &Q, real_t dt)
{
   const RKTableau tab = MakeRK4Tableau();
   const int s = tab.stages, h = Q.Size(), n = (int)dof.size();
   std::vector<Vector> k(s, Vector(h));
   std::vector<real_t> s1n(n), s2n(n);
   std::vector<std::vector<real_t>> V1(s, std::vector<real_t>(n, 0.0)),
                                    V2(s, std::vector<real_t>(n, 0.0));
   for (int m = 0; m < n; m++) { s1n[m] = dof[m].slip1; s2n[m] = dof[m].slip2; }
   Vector Qs(h);
   for (int i = 0; i < s; i++)
   {
      Qs = Q;
      for (int j = 0; j < i; j++)
         if (tab.a[i][j] != 0.0) { Qs.Add(dt * tab.a[i][j], k[j]); }
      for (int m = 0; m < n; m++)
      {
         real_t a1 = s1n[m], a2 = s2n[m];
         for (int j = 0; j < i; j++)
         { a1 += dt*tab.a[i][j]*V1[j][m]; a2 += dt*tab.a[i][j]*V2[j][m]; }
         dof[m].slip1 = a1; dof[m].slip2 = a2;
      }
      wave.Mult(Qs, k[i]);
      for (int m = 0; m < n; m++) { V1[i][m] = dof[m].V1; V2[i][m] = dof[m].V2; }
   }
   LSWRef r; r.Q_new.SetSize(h); r.Q_new = Q;
   for (int i = 0; i < s; i++)
      if (tab.b[i] != 0.0) { r.Q_new.Add(dt * tab.b[i], k[i]); }
   r.slip1_abs.resize(n); r.slip2_abs.resize(n);
   r.slip1_bug.resize(n); r.slip2_bug.resize(n);
   for (int m = 0; m < n; m++)
   {
      real_t inc1 = 0.0, inc2 = 0.0;
      for (int i = 0; i < s; i++) { inc1 += dt*tab.b[i]*V1[i][m];
                                    inc2 += dt*tab.b[i]*V2[i][m]; }
      r.slip1_abs[m] = s1n[m] + inc1;        r.slip2_abs[m] = s2n[m] + inc2;
      // `+=` bug: dof[m].slip after the stage loop holds the LAST stage's
      // staged slip, NOT s_n; a `+=` combine adds inc on top of that.
      r.slip1_bug[m] = dof[m].slip1 + inc1;  r.slip2_bug[m] = dof[m].slip2 + inc2;
   }
   return r;
}

// ===========================================================================
// L5 — bulk equivalence on a frictionless cube
// ===========================================================================
void L5_BulkEquivalence()
{
   std::cout << "\n-- L5: AdvanceRKCoupledLSW_Spatial(RK4) == DoRK4Step (bulk) --\n";
   const double L = 8.0e3; const int nx = 6;
   Mesh mesh = Mesh::MakeCartesian3D(nx, nx, nx, Element::HEXAHEDRON, L, L, L);
   for (int v = 0; v < mesh.GetNV(); v++)
   { real_t* x = mesh.GetVertex(v); x[0]-=0.5*L; x[1]-=0.5*L; x[2]-=0.5*L; }
   for (int be = 0; be < mesh.GetNBE(); be++) { mesh.SetBdrAttribute(be, 5); }
   mesh.SetAttributes();

   BoundaryConfig bc; bc.natural_attrs = {}; bc.fault_attr = 0;
   bc.absorbing_attrs = {5};
   const int order = 1;
   WaveOperator<Mesh> wave(mesh, order, TPV205Params::lambda,
                           TPV205Params::mu, TPV205Params::rho, bc);
   { real_t zb[NUM_STATE] = {0}; wave.SetAbsorbingBackground(zb); }
   // LSW stepper requires the LSW flag (no fault faces ⇒ bulk-only effect).
   wave.SetFaultFrictionLaw(FaultFrictionLaw::LSW);

   const int size = wave.Height();
   auto& fes = const_cast<FiniteElementSpace&>(wave.GetFESpace());
   const int ndof_total = fes.GetNDofs();
   const double sig = L/8.0, amp = 1.0e6;
   FunctionCoefficient syy([sig,amp](const Vector& x)->real_t {
      double r2 = x(0)*x(0)+x(1)*x(1)+x(2)*x(2);
      return amp*std::exp(-r2/(sig*sig)); });
   Vector Q(size); Q = 0.0;
   GridFunction gf(&fes); gf.ProjectCoefficient(syy);
   for (int d = 0; d < ndof_total; d++)
   { Q(SYY*ndof_total+d) = gf(d); Q(SXX*ndof_total+d) = 0.5*gf(d); }

   const real_t dt = 5.0e-5;
   Vector Q_ref(Q); DoRK4Step_local(wave, Q_ref, dt);

   std::vector<DOFData> dof_empty;
   Vector Q_rk(size);
   AdvanceRKCoupledLSW_Spatial<Mesh>(wave, dof_empty, Q, dt, 0.0, Q_rk,
                                     MakeRK4Tableau(), nullptr);
   real_t max_abs = 0.0, max_diff = 0.0;
   for (int i = 0; i < size; i++)
   { max_abs = std::max(max_abs, std::abs(Q_ref(i)));
     max_diff = std::max(max_diff, std::abs(Q_ref(i) - Q_rk(i))); }
   const real_t rel = max_diff / std::max(max_abs, real_t(1e-30));
   std::cout << "    rel = " << std::scientific << std::setprecision(4) << rel << "\n";
   TEST_TRUE(rel < 1e-13,
             "LSW RK4 bulk step reproduces DoRK4Step to round-off");
}

// ===========================================================================
// L6 — coupled (Q, slip) RK4 matches manual reference; absolute != buggy
// ===========================================================================
void L6_SlipCombine()
{
   std::cout << "\n-- L6: coupled (Q, slip) RK4 + absolute-combine (R-002) --\n";
   const int order = 1;
   Mesh mesh = BuildSmallFaultTetMesh(1000.0);
   BoundaryConfig bc; bc.natural_attrs = {1}; bc.fault_attr = 3;
   bc.absorbing_attrs = {};

   // Subject operator + dof_data.
   WaveOperator<Mesh> wave(mesh, order, TPV205Params::lambda,
                           TPV205Params::mu, TPV205Params::rho, bc);
   wave.SetFaultFrictionLaw(FaultFrictionLaw::LSW);
   FaultFaceFlux ff_sub(TPV205Params::rho, TPV205Params::cp, TPV205Params::cs);
   std::vector<DOFData> dof_sub;
   SetupLSWFault(wave, ff_sub, dof_sub, mesh, order);
   const int n = (int)dof_sub.size();

   Vector Q(wave.Height()); Q = 0.0;
   const real_t dt = 1.0e-4;

   // Reference operator + dof_data (fresh, identical init).
   WaveOperator<Mesh> wave_ref(mesh, order, TPV205Params::lambda,
                               TPV205Params::mu, TPV205Params::rho, bc);
   wave_ref.SetFaultFrictionLaw(FaultFrictionLaw::LSW);
   FaultFaceFlux ff_ref(TPV205Params::rho, TPV205Params::cp, TPV205Params::cs);
   std::vector<DOFData> dof_ref;
   SetupLSWFault(wave_ref, ff_ref, dof_ref, mesh, order);
   LSWRef ref = DoLSWRK4Reference(wave_ref, dof_ref, Q, dt);

   // Subject: the production stepper.
   Vector Q_rk(wave.Height());
   AdvanceRKCoupledLSW_Spatial<Mesh>(wave, dof_sub, Q, dt, 0.0, Q_rk,
                                     MakeRK4Tableau(), nullptr);

   // (a) Q matches the manual reference to round-off.
   real_t qa = 0.0, qd = 0.0;
   for (int i = 0; i < Q_rk.Size(); i++)
   { qa = std::max(qa, std::abs(ref.Q_new(i)));
     qd = std::max(qd, std::abs(ref.Q_new(i) - Q_rk(i))); }
   TEST_TRUE(qd / std::max(qa, real_t(1e-30)) < 1e-12,
             "stepper Q_new matches manual coupled-RK4 reference");

   // (b) slip matches the ABSOLUTE-combine reference to round-off.
   real_t sd = 0.0, slip_scale = 0.0;
   bool any_slid = false, discriminates = false;
   for (int m = 0; m < n; m++)
   {
      sd = std::max(sd, std::abs(dof_sub[m].slip1 - ref.slip1_abs[m]));
      sd = std::max(sd, std::abs(dof_sub[m].slip2 - ref.slip2_abs[m]));
      slip_scale = std::max({slip_scale, std::abs(ref.slip1_abs[m]),
                             std::abs(ref.slip2_abs[m])});
      if (std::abs(dof_sub[m].slip2) > 1e-9) { any_slid = true; }
      // (c) the absolute combine must differ measurably from the buggy `+=`.
      if (std::abs(ref.slip2_abs[m] - ref.slip2_bug[m]) > 1e-9) { discriminates = true; }
   }
   TEST_TRUE(sd <= 1e-12 * slip_scale + 1e-18,
             "stepper slip matches ABSOLUTE-combine reference (= not +=)");
   TEST_TRUE(any_slid, "fault slid (slip2 != 0) so the combine is exercised");
   TEST_TRUE(discriminates,
             "absolute-combine slip differs from buggy += slip (test discriminates)");

   // R-005: independent analytic anchor — the endpoint slip-rate must be the
   // LSW closed-form rupture-initiation value to ORDER OF MAGNITUDE.  At t=0
   // (Q=0) the analytic V = (tau_nuc - mu_s*sigma_n)/eta_s ≈ 0.0779 m/s; over
   // one dt the radiation-damping feedback perturbs it O(1) (the radiated trial
   // traction is comparable to the driving stress), so a TIGHT match is not
   // expected — but V staying within a factor of 2 of the analytic value is an
   // oracle that does NOT depend on the stepper's internal staging logic (it
   // comes from TPV205Params), so it catches a shared staging/physics error the
   // reference-match (which mirrors the stepper) could miss.
   const real_t V_analytic =
      (TPV205Params::tau_nuc - TPV205Params::mu_s * TPV205Params::sigma_n)
      / TPV205Params::eta_s;
   real_t Vmax = 0.0;
   for (int m = 0; m < n; m++) { Vmax = std::max(Vmax, std::abs(dof_sub[m].V2)); }
   TEST_TRUE(Vmax > 0.5 * V_analytic && Vmax < 2.0 * V_analytic,
             "endpoint slip-rate is the analytic LSW V to O(1) (independent oracle)");
}

// ===========================================================================
// L7 — FSAL predicate (the endpoint-re-eval gate)
// ===========================================================================
bool LastStageIsEndpoint(const RKTableau& tab)
{
   for (int j = 0; j < tab.stages; ++j)
   { if (tab.a[tab.stages - 1][j] != tab.b[j]) { return false; } }
   return true;
}
void L7_FSAL()
{
   std::cout << "\n-- L7: FSAL predicate (RK4 non-FSAL, DP45 FSAL) --\n";
   TEST_TRUE(!LastStageIsEndpoint(MakeRK4Tableau()),
             "RK4 is NOT FSAL ⇒ endpoint re-eval performed");
   TEST_TRUE(LastStageIsEndpoint(MakeDormandPrinceRK45Tableau()),
             "DP45 IS FSAL ⇒ endpoint re-eval skipped");
}

} // namespace

// Reference DoRK4Step (verbatim idiom) — declared after the anonymous namespace
// use above via a forward at file scope.
void DoRK4Step_local(const WaveOperator<Mesh>& wave, Vector& Q, real_t dt);

int main(int argc, char* argv[])
{
   (void)argc; (void)argv;
   std::cout << "\n=== Phase 14 — RK + LSW + mixed-flux unit tests (L1–L7) ===\n";
   L1_EvaluateLSW_vs_ADER();
   L2_SlipStateless();
   L3_BarrierAndTensile();
   L4_MultDispatch();
   L5_BulkEquivalence();
   L6_SlipCombine();
   L7_FSAL();
   std::cout << "\n========================================\n"
             << "  Results: " << num_passed << " passed, " << num_failed
             << " failed out of " << num_tests << " tests\n"
             << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}

void DoRK4Step_local(const WaveOperator<Mesh>& wave, Vector& Q, real_t dt)
{
   const int size = Q.Size();
   Vector k1(size), k2(size), k3(size), k4(size), Q_tmp(size);
   wave.Mult(Q, k1);
   add(Q, 0.5*dt, k1, Q_tmp); wave.Mult(Q_tmp, k2);
   add(Q, 0.5*dt, k2, Q_tmp); wave.Mult(Q_tmp, k3);
   add(Q, dt,     k3, Q_tmp); wave.Mult(Q_tmp, k4);
   Q.Add(dt/6.0, k1); Q.Add(dt/3.0, k2); Q.Add(dt/3.0, k3); Q.Add(dt/6.0, k4);
}
