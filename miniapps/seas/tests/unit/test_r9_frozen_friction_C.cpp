// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// Round-9 Step 1 (2026-04-23) + Round-11 extension (2026-04-23):
// FREEZE-A/B/C amplifier discriminator + dt-refinement baseline,
// cross-fixture (Kuhn + fixed D4).
//
// Round-9 spec (R-003): "frozen friction" is three distinct tests:
//   FREEZE-A: skip friction solve entirely (fault locked).  Requires
//             production hook in FaultFaceFlux::EvaluateADER guarded
//             by SEAS_TEST_FREEZE_A=1 (leaves DOFData unchanged and
//             passes Q_avg_plus/Q_avg_minus through as the imposed
//             state — no flux correction applied).
//   FREEZE-B: cache step-1 DOFData (psi, slip_rate, V1, V2,
//             tau1_corr, tau2_corr, sigma_n_corr, slip1, slip2),
//             restore it verbatim before every step ≥ 2.  IMPORTANT
//             (Round-12 Patch-0 correction): among these fields,
//             ONLY psi is a true causal input to the next
//             FaultFaceFlux::Evaluate() call.  The remaining
//             fields — slip_rate, V1/V2, tau*_corr, sigma_n_corr,
//             slip1, slip2 — are OUTPUTS of Evaluate() and play no
//             role in the next call (see fault_face_flux.cpp:70).
//             The actual per-step inputs to Evaluate() are bulk Q
//             (from AdvanceADER), data.psi, data.tau*_0,
//             data.tau*_nuc, data.sigma_n0/_nuc, and the static
//             impedances / material params.  FREEZE-B therefore
//             behaves identically to FREEZE-C in terms of friction
//             inputs — the only difference is that FREEZE-B also
//             restores the reported `slip_rate` / `tau*_corr`
//             values captured AFTER step-1's Collect(), which
//             freezes the reported WORST spread at the step-1
//             magnitude regardless of the live friction
//             trajectory.  The test's `Collect()` is what observes
//             these fields, so FREEZE-B's "ratio = 1.0000" was a
//             *measurement artifact*, not evidence that the full
//             state restore was neutralising the amplifier.
//   FREEZE-C: hold ψ constant after step 1 but solve normally.  Pure
//             test-level intervention — restores only data.psi.
//
// Round-11 decision question (from user plan): is the first
// nonuniformity created by the MFEM wave update even when friction
// is frozen (→ amplifier outside the fault), or only when live fault
// outputs are fed back (→ amplifier inside the fault)?
//
// dt-refinement (R-004): run BASELINE pepper guard at dt, dt/2, dt/4
// and record tau1_corr spread.  If spread ∝ dt^p with p > 0 → pepper
// is numerical error that vanishes as dt → 0 (sub-resolved physics,
// not a code bug).  If spread dt-insensitive → code bug.
//
// Cross-fixture (R-002): run all tests on BOTH Kuhn and fixed D4.
//
// Pre-committed Step-2 action matrix (R-005):
//   FREEZE-C CLOSES pepper AND dt-refinement gives spread ∝ dt^p >0 →
//       round-10 target is "sub-resolved rate-state physics".
//       Fix scope: dt refinement, ψ-integrator precision, or
//       explicit damping.  Line-level target:
//       friction/dieterich_ruina.hpp::UpdateStateAnalytic.
//   FREEZE-C CLOSES AND dt-refinement insensitive →
//       round-10 target is a CODE bug in the ψ-feedback coupling.
//       Same file but different line: the stage-to-stage ψ
//       interleave in tpv102_driver.cpp RK4 loop.
//   FREEZE-C does NOT close AND dt-refinement gives spread ∝ dt^p >0 →
//       amplifier is in wave kernel physics (sub-resolved), NOT in
//       rate-state.  Fix scope: mesh refinement / dt refinement /
//       SIPG tightening.
//   FREEZE-C does NOT close AND dt-refinement insensitive →
//       code bug OUTSIDE the friction solver AND outside ψ.
//       Requires flux-layer freeze unblock to disambiguate via
//       FREEZE-A / FREEZE-B.
//   NONE of FREEZE-C branches close → amplifier not in ψ.
//
// Gate: test exits with 0 always (diagnostic only); raw output is
// the verdict, not the exit code.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/tpv102_setup.hpp"
#include "../../dynamic/d4_tet_mesh.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace mfem { namespace seas { int g_seas_my_rank = 0; } }  // NOLINT

namespace
{

constexpr real_t kL = 1000.0;
constexpr int    kOrder = 1;
constexpr int    kAderOrder = 2;

Mesh BuildKuhnFaultMesh()
{
   Mesh mesh = Mesh::MakeCartesian3D(2, 2, 2, Element::TETRAHEDRON,
                                     kL, kL, kL, /*sfc_ordering=*/false);
   mesh.FinalizeTopology();
   mesh.Finalize();
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
         if (std::abs(cy - 0.5 * kL) < 1e-8)
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

Mesh BuildFaultMeshByEnv()
{
   const char *fx = std::getenv("SEAS_TEST_FIXTURE");
   if (fx && std::string(fx) == "d4")
   {
      std::cout << "  [SEAS_TEST_FIXTURE=d4] using BuildD4Mesh(true).\n";
      return BuildD4Mesh(/*add_fault=*/true, kL);
   }
   std::cout << "  [fixture=kuhn] (default)\n";
   return BuildKuhnFaultMesh();
}

struct Spreads
{
   real_t slip_rate = 0;
   real_t tau1_corr = 0;
   real_t tau2_corr = 0;
   real_t sigma_n_corr = 0;
};

// Round-11 Step 1: four-way amplifier discriminator.  FREEZE-A and
// FREEZE-B require authorization for a flux-layer hook; FREEZE-C is
// purely test-level.  See file header for semantics.
enum class FreezeMode { None, FreezeA, FreezeB, FreezeC };

const char *FreezeModeLabel(FreezeMode m)
{
   switch (m)
   {
      case FreezeMode::None:    return "none";
      case FreezeMode::FreezeA: return "FREEZE-A (skip EvaluateADER friction solve)";
      case FreezeMode::FreezeB: return "FREEZE-B (cache+restore full DOFData friction state)";
      case FreezeMode::FreezeC: return "FREEZE-C (ψ held after step 1)";
   }
   return "?";
}

// Snapshot of the friction-writeback fields held on DOFData.  Captured
// after step 1; restored before every subsequent step when FREEZE-B or
// FREEZE-C is active.  FREEZE-C uses only `psi`; FREEZE-B uses all.
struct FrictionSnapshot
{
   std::vector<real_t> psi;
   std::vector<real_t> slip_rate;
   std::vector<real_t> V1, V2;
   std::vector<real_t> tau1_corr, tau2_corr, sigma_n_corr;
   std::vector<real_t> slip1, slip2;
   bool captured = false;
};

void CaptureFrictionSnapshot(FrictionSnapshot &snap,
                             const std::vector<DOFData> &dof)
{
   const int n = static_cast<int>(dof.size());
   snap.psi.resize(n);
   snap.slip_rate.resize(n);
   snap.V1.resize(n); snap.V2.resize(n);
   snap.tau1_corr.resize(n); snap.tau2_corr.resize(n);
   snap.sigma_n_corr.resize(n);
   snap.slip1.resize(n); snap.slip2.resize(n);
   for (int i = 0; i < n; i++)
   {
      snap.psi[i]          = dof[i].psi;
      snap.slip_rate[i]    = dof[i].slip_rate;
      snap.V1[i]           = dof[i].V1;
      snap.V2[i]           = dof[i].V2;
      snap.tau1_corr[i]    = dof[i].tau1_corr;
      snap.tau2_corr[i]    = dof[i].tau2_corr;
      snap.sigma_n_corr[i] = dof[i].sigma_n_corr;
      snap.slip1[i]        = dof[i].slip1;
      snap.slip2[i]        = dof[i].slip2;
   }
   snap.captured = true;
}

void RestoreFrictionSnapshot(const FrictionSnapshot &snap,
                             std::vector<DOFData> &dof, FreezeMode mode)
{
   if (!snap.captured) { return; }
   const int n = static_cast<int>(dof.size());
   for (int i = 0; i < n; i++)
   {
      if (mode == FreezeMode::FreezeC)
      {
         dof[i].psi = snap.psi[i];
      }
      else if (mode == FreezeMode::FreezeB)
      {
         dof[i].psi          = snap.psi[i];
         dof[i].slip_rate    = snap.slip_rate[i];
         dof[i].V1           = snap.V1[i];
         dof[i].V2           = snap.V2[i];
         dof[i].tau1_corr    = snap.tau1_corr[i];
         dof[i].tau2_corr    = snap.tau2_corr[i];
         dof[i].sigma_n_corr = snap.sigma_n_corr[i];
         dof[i].slip1        = snap.slip1[i];
         dof[i].slip2        = snap.slip2[i];
      }
      // FreezeMode::None or FreezeA: nothing to restore at the test level.
   }
}

// Round-12 Patch 3: generic RAII env-var helper.  Replaces the
// previous FreezeAEnvScope; used both for SEAS_TEST_FREEZE_A and for
// SEAS_TEST_EVAL_FACE_AVG_STAGE.  Null/empty `var_name` → inactive,
// so callers can construct unconditionally and opt in by passing
// nullptr when the case doesn't need the env.
struct ScopedEnvVar
{
   std::string name;
   ScopedEnvVar() = default;
   ScopedEnvVar(const char *var_name, const char *value)
      : name(var_name ? var_name : "")
   {
      if (!name.empty())
      {
         setenv(name.c_str(), value ? value : "", 1);
      }
   }
   ~ScopedEnvVar()
   {
      if (!name.empty()) { unsetenv(name.c_str()); }
   }
   ScopedEnvVar(const ScopedEnvVar &) = delete;
   ScopedEnvVar &operator=(const ScopedEnvVar &) = delete;
};

// Round-12 Patch 3: per-face stage averaging modes consumed by the
// SEAS_TEST_EVAL_FACE_AVG_STAGE hook in wave_operator.inl.  Env
// strings must match the parser in ComputeADERFaceFluxRHS.
enum class AvgMode { None, Trial, Theta, Vabs, Tcorr };

const char *AvgModeEnvValue(AvgMode m)
{
   switch (m)
   {
      case AvgMode::None:  return "";
      case AvgMode::Trial: return "trial";
      case AvgMode::Theta: return "theta";
      case AvgMode::Vabs:  return "vabs";
      case AvgMode::Tcorr: return "tcorr";
   }
   return "";
}

const char *AvgModeLabel(AvgMode m)
{
   switch (m)
   {
      case AvgMode::None:  return "none";
      case AvgMode::Trial: return "AVG_TRIAL (face-mean trial traction)";
      case AvgMode::Theta: return "AVG_THETA (face-mean Θ magnitude)";
      case AvgMode::Vabs:  return "AVG_VABS  (face-mean |V|)";
      case AvgMode::Tcorr: return "AVG_TCORR (face-mean corrected traction)";
   }
   return "?";
}

Spreads Collect(const std::vector<DOFData> &dof)
{
   Spreads s;
   if (dof.empty()) { return s; }
   real_t sr_min = dof[0].slip_rate, sr_max = sr_min;
   real_t t1_min = dof[0].tau1_corr, t1_max = t1_min;
   real_t t2_min = dof[0].tau2_corr, t2_max = t2_min;
   real_t sn_min = dof[0].sigma_n_corr, sn_max = sn_min;
   for (const DOFData &d : dof)
   {
      sr_min = std::min(sr_min, d.slip_rate); sr_max = std::max(sr_max, d.slip_rate);
      t1_min = std::min(t1_min, d.tau1_corr); t1_max = std::max(t1_max, d.tau1_corr);
      t2_min = std::min(t2_min, d.tau2_corr); t2_max = std::max(t2_max, d.tau2_corr);
      sn_min = std::min(sn_min, d.sigma_n_corr); sn_max = std::max(sn_max, d.sigma_n_corr);
   }
   s.slip_rate = sr_max - sr_min;
   s.tau1_corr = t1_max - t1_min;
   s.tau2_corr = t2_max - t2_min;
   s.sigma_n_corr = sn_max - sn_min;
   return s;
}

// Run one pepper-guard simulation with specified dt, optional
// friction-freeze mode, optional face-averaging mode, and optional
// non-fault both-side symmetrization.  Freeze / averaging /
// nonfault-sym are orthogonal — a given call activates at most one
// of each via the RAII env scopes below.
Spreads RunPepperGuard(const char *label, real_t dt, int n_steps,
                       FreezeMode freeze_mode,
                       AvgMode avg_mode = AvgMode::None,
                       bool nonfault_both_sym = false)
{
   std::cout << "\n=== " << label << "  (dt=" << dt
             << ", n_steps=" << n_steps
             << ", freeze=" << FreezeModeLabel(freeze_mode)
             << ", avg=" << AvgModeLabel(avg_mode)
             << ", nf_both_sym=" << (nonfault_both_sym ? "yes" : "no")
             << ") ===\n";

   // FREEZE-A env flag: activates the hook in EvaluateADER for the
   // lifetime of this run only.  RAII guard unsets on scope exit so
   // later runs in the same process are unaffected.
   ScopedEnvVar freeze_a_scope(
      freeze_mode == FreezeMode::FreezeA ? "SEAS_TEST_FREEZE_A" : nullptr,
      "1");

   // Round-12 Patch 3: per-face stage averaging.  SEAS_TEST_EVAL_FACE_AVG_STAGE
   // triggers the face-averaging hook in ComputeADERFaceFluxRHS (wave_operator.inl).
   ScopedEnvVar avg_scope(
      avg_mode != AvgMode::None ? "SEAS_TEST_EVAL_FACE_AVG_STAGE" : nullptr,
      AvgModeEnvValue(avg_mode));

   // Round-13C Patch 1: non-fault branch n↔-n symmetrization hook.
   // Activates the `SEAS_TEST_NONFAULT_BOTH_SYM=1` path in
   // ComputeADERFaceFluxRHS (interior non-fault + boundary).
   ScopedEnvVar nonfault_sym_scope(
      nonfault_both_sym ? "SEAS_TEST_NONFAULT_BOTH_SYM" : nullptr, "1");

   Mesh mesh = BuildFaultMeshByEnv();
   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr    = 3;
   real_t bulk_bg[NUM_STATE] = {0};

   WaveOperator<Mesh> wave(mesh, kOrder, TPV102Params::lambda,
                            TPV102Params::mu, TPV102Params::rho, bc);
   wave.SetAbsorbingBackground(bulk_bg);

   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   int nqp = 0;
   if (int_faces.Size() > 0)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[0]);
      nqp = IntRules.Get(ftr->GetGeometryType(), 2*kOrder).GetNPoints();
   }
   std::vector<Vector> fault_coords;
   for (int i = 0; i < int_faces.Size(); i++)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[i]);
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2*kOrder);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3); ftr->Face->Transform(ip, phys);
         fault_coords.push_back(phys);
      }
   }
   const int n_fault = int_faces.Size() * nqp;
   std::vector<DOFData> dof_data;
   InitializeFaultDOFs(dof_data, n_fault, fault_coords);
   for (int i = 0; i < n_fault; i++)
   { dof_data[i].tau2_nuc = TPV102Params::nuc_dtau; }

   FaultFaceFlux ff(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   wave.SetFaultFlux(&ff);
   wave.SetFaultDOFData(&dof_data, nqp);

   Vector Q(wave.Height()), Q_new(wave.Height());
   Q = 0.0;

   FrictionSnapshot snap;
   const bool uses_test_level_restore =
      (freeze_mode == FreezeMode::FreezeB) ||
      (freeze_mode == FreezeMode::FreezeC);

   Spreads worst;
   for (int step = 0; step < n_steps; step++)
   {
      // FREEZE-B / FREEZE-C: restore cached friction state before
      // steps >= 2.  FREEZE-A operates inside EvaluateADER via the
      // env hook — DOFData is never touched so no restore is needed.
      if (uses_test_level_restore && step >= 2 && snap.captured)
      {
         RestoreFrictionSnapshot(snap, dof_data, freeze_mode);
      }

      wave.AdvanceADER(Q, dt, kAderOrder, Q_new);
      Q.Swap(Q_new);

      // Capture snapshot AFTER step 1 so step >= 2 starts from the
      // step-1 friction state.
      if (uses_test_level_restore && step == 1)
      {
         CaptureFrictionSnapshot(snap, dof_data);
      }

      Spreads s = Collect(dof_data);
      worst.slip_rate    = std::max(worst.slip_rate,    s.slip_rate);
      worst.tau1_corr    = std::max(worst.tau1_corr,    s.tau1_corr);
      worst.tau2_corr    = std::max(worst.tau2_corr,    s.tau2_corr);
      worst.sigma_n_corr = std::max(worst.sigma_n_corr, s.sigma_n_corr);
   }

   std::cout << "  worst spreads over " << n_steps << " steps:\n";
   std::cout << "    slip_rate    : " << std::scientific
             << std::setprecision(6) << worst.slip_rate << "\n";
   std::cout << "    tau1_corr    : " << worst.tau1_corr << "\n";
   std::cout << "    tau2_corr    : " << worst.tau2_corr << "\n";
   std::cout << "    sigma_n_corr : " << worst.sigma_n_corr << "\n";
   return worst;
}

} // anonymous

int main()
{
   std::cout << "\n========================================\n"
             << "  Round-12: FREEZE-A/B/C + AVG stage-averaging + dt\n"
             << "  Fixture: " << (std::getenv("SEAS_TEST_FIXTURE") ? std::getenv("SEAS_TEST_FIXTURE") : "kuhn") << "\n"
             << "========================================\n";

   // Baseline at the nominal dt used by the pepper guard.
   const real_t kDtNominal = 5.0e-5;
   const int    kNStepsNominal = 20;
   // For dt/2 and dt/4, run 2x and 4x the number of steps to reach
   // the same PHYSICAL time.  Worst-spread measurement captures the
   // max spread throughout the run; at dt/4 we expect the worst
   // spread to occur later in physical time but at the same
   // physical state.
   //
   // However, computing the full 80 steps at dt/4 might be slow and
   // also the eigenmode analysis wants spread-at-same-physical-time.
   // Simplification: run all three at the SAME number of steps (20),
   // so the physical time covered is different (shorter for dt/2,
   // dt/4).  The ratio of WORST-OVER-20-STEPS as a function of dt
   // still captures the dt-dependence of the numerical error per
   // physical time → useful for the verdict.
   //
   // Trade-off: shorter physical time may miss peak pepper.  But
   // for 20 dt/4 steps at dt/4 = 1.25e-5, physical time is 2.5e-4 s,
   // vs baseline 1.0e-3 s.  The linearly-growing wave only covers
   // 1/4 the distance, which is OK — the t^2 pepper growth means
   // pepper at step 5 in baseline corresponds to pepper at step 20
   // in dt/4 roughly.  Good enough for verdict.
   //
   // R-004 refinement: to make fair comparison, run all three to the
   // SAME PHYSICAL TIME: dt=5e-5 × 20 steps, dt=2.5e-5 × 40 steps,
   // dt=1.25e-5 × 80 steps.  Then the worst spread is measured at
   // the same physical time.
   const real_t t_final = kDtNominal * kNStepsNominal;  // 1e-3 s

   Spreads baseline = RunPepperGuard("BASELINE (dt=5e-5, 20 steps)",
                                       kDtNominal, kNStepsNominal,
                                       FreezeMode::None);

   Spreads freezeA = RunPepperGuard("FREEZE-A (dt=5e-5, 20 steps, skip fault friction solve)",
                                      kDtNominal, kNStepsNominal,
                                      FreezeMode::FreezeA);

   Spreads freezeB = RunPepperGuard("FREEZE-B (dt=5e-5, 20 steps, full-state restore)",
                                      kDtNominal, kNStepsNominal,
                                      FreezeMode::FreezeB);

   Spreads freezeC = RunPepperGuard("FREEZE-C (dt=5e-5, 20 steps, ψ held after step 1)",
                                      kDtNominal, kNStepsNominal,
                                      FreezeMode::FreezeC);

   // Round-12 AVG stage-averaging cases (all at nominal dt, no freeze).
   Spreads avgTrial = RunPepperGuard("AVG_TRIAL (face-mean trial traction)",
                                      kDtNominal, kNStepsNominal,
                                      FreezeMode::None, AvgMode::Trial);

   Spreads avgTheta = RunPepperGuard("AVG_THETA (face-mean Theta)",
                                      kDtNominal, kNStepsNominal,
                                      FreezeMode::None, AvgMode::Theta);

   Spreads avgVabs  = RunPepperGuard("AVG_VABS  (face-mean |V|)",
                                      kDtNominal, kNStepsNominal,
                                      FreezeMode::None, AvgMode::Vabs);

   Spreads avgTcorr = RunPepperGuard("AVG_TCORR (face-mean corrected traction)",
                                      kDtNominal, kNStepsNominal,
                                      FreezeMode::None, AvgMode::Tcorr);

   // Round-13C Patch 2: non-fault both-side symmetrization cases.
   Spreads nonfaultBothSym = RunPepperGuard(
      "NONFAULT_BOTH_SYM (dt=5e-5, 20 steps)",
      kDtNominal, kNStepsNominal,
      FreezeMode::None, AvgMode::None, /*nonfault_both_sym=*/true);

   Spreads freezeA_nfBothSym = RunPepperGuard(
      "FREEZE-A + NONFAULT_BOTH_SYM (dt=5e-5, 20 steps)",
      kDtNominal, kNStepsNominal,
      FreezeMode::FreezeA, AvgMode::None, /*nonfault_both_sym=*/true);

   const real_t kDtHalf = kDtNominal * 0.5;
   const int    kNStepsHalf = static_cast<int>(t_final / kDtHalf + 0.5);
   Spreads dt_half = RunPepperGuard("DT/2 baseline (dt=2.5e-5, 40 steps)",
                                      kDtHalf, kNStepsHalf,
                                      FreezeMode::None);

   const real_t kDtQuart = kDtNominal * 0.25;
   const int    kNStepsQuart = static_cast<int>(t_final / kDtQuart + 0.5);
   Spreads dt_quart = RunPepperGuard("DT/4 baseline (dt=1.25e-5, 80 steps)",
                                       kDtQuart, kNStepsQuart,
                                       FreezeMode::None);

   // Summary + verdict.
   std::cout << "\n========================================\n"
             << "  VERDICT TABLE\n"
             << "========================================\n";
   std::cout << std::scientific << std::setprecision(6);
   std::cout << "                 | slip_rate  | tau1_corr  | tau2_corr  | sigma_n_corr\n";
   std::cout << "  ---------------|------------|------------|------------|--------------\n";
   std::cout << "  BASELINE       | " << baseline.slip_rate
             << " | " << baseline.tau1_corr
             << " | " << baseline.tau2_corr
             << " | " << baseline.sigma_n_corr << "\n";
   std::cout << "  FREEZE-A       | " << freezeA.slip_rate
             << " | " << freezeA.tau1_corr
             << " | " << freezeA.tau2_corr
             << " | " << freezeA.sigma_n_corr << "\n";
   std::cout << "  FREEZE-B       | " << freezeB.slip_rate
             << " | " << freezeB.tau1_corr
             << " | " << freezeB.tau2_corr
             << " | " << freezeB.sigma_n_corr << "\n";
   std::cout << "  FREEZE-C       | " << freezeC.slip_rate
             << " | " << freezeC.tau1_corr
             << " | " << freezeC.tau2_corr
             << " | " << freezeC.sigma_n_corr << "\n";
   std::cout << "  AVG_TRIAL      | " << avgTrial.slip_rate
             << " | " << avgTrial.tau1_corr
             << " | " << avgTrial.tau2_corr
             << " | " << avgTrial.sigma_n_corr << "\n";
   std::cout << "  AVG_THETA      | " << avgTheta.slip_rate
             << " | " << avgTheta.tau1_corr
             << " | " << avgTheta.tau2_corr
             << " | " << avgTheta.sigma_n_corr << "\n";
   std::cout << "  AVG_VABS       | " << avgVabs.slip_rate
             << " | " << avgVabs.tau1_corr
             << " | " << avgVabs.tau2_corr
             << " | " << avgVabs.sigma_n_corr << "\n";
   std::cout << "  AVG_TCORR      | " << avgTcorr.slip_rate
             << " | " << avgTcorr.tau1_corr
             << " | " << avgTcorr.tau2_corr
             << " | " << avgTcorr.sigma_n_corr << "\n";
   std::cout << "  NF_BOTH_SYM    | " << nonfaultBothSym.slip_rate
             << " | " << nonfaultBothSym.tau1_corr
             << " | " << nonfaultBothSym.tau2_corr
             << " | " << nonfaultBothSym.sigma_n_corr << "\n";
   std::cout << "  FA + NF_BOTH_SYM| " << freezeA_nfBothSym.slip_rate
             << " | " << freezeA_nfBothSym.tau1_corr
             << " | " << freezeA_nfBothSym.tau2_corr
             << " | " << freezeA_nfBothSym.sigma_n_corr << "\n";
   std::cout << "  DT/2 baseline  | " << dt_half.slip_rate
             << " | " << dt_half.tau1_corr
             << " | " << dt_half.tau2_corr
             << " | " << dt_half.sigma_n_corr << "\n";
   std::cout << "  DT/4 baseline  | " << dt_quart.slip_rate
             << " | " << dt_quart.tau1_corr
             << " | " << dt_quart.tau2_corr
             << " | " << dt_quart.sigma_n_corr << "\n";

   auto ratio = [](real_t a, real_t b)
   { return (b > 0) ? a / b : 0.0; };
   const real_t r_freezeA = ratio(freezeA.tau1_corr, baseline.tau1_corr);
   const real_t r_freezeB = ratio(freezeB.tau1_corr, baseline.tau1_corr);
   const real_t r_freezeC = ratio(freezeC.tau1_corr, baseline.tau1_corr);
   const real_t r_avgTrial = ratio(avgTrial.tau1_corr, baseline.tau1_corr);
   const real_t r_avgTheta = ratio(avgTheta.tau1_corr, baseline.tau1_corr);
   const real_t r_avgVabs  = ratio(avgVabs.tau1_corr,  baseline.tau1_corr);
   const real_t r_avgTcorr = ratio(avgTcorr.tau1_corr, baseline.tau1_corr);
   const real_t r_nonfault_both_sym = ratio(nonfaultBothSym.tau1_corr, baseline.tau1_corr);
   const real_t r_freezeA_nfBothSym = ratio(freezeA_nfBothSym.tau1_corr, baseline.tau1_corr);
   const real_t r_dt2 = ratio(dt_half.tau1_corr, baseline.tau1_corr);
   const real_t r_dt4 = ratio(dt_quart.tau1_corr, baseline.tau1_corr);

   std::cout << "\n=== Key ratios (tau1_corr, relative to baseline) ===\n"
             << std::fixed << std::setprecision(4);
   std::cout << "  FREEZE-A / baseline = " << r_freezeA << "\n";
   std::cout << "  FREEZE-B / baseline = " << r_freezeB << "\n";
   std::cout << "  FREEZE-C / baseline = " << r_freezeC << "\n";
   std::cout << "  AVG_TRIAL/ baseline = " << r_avgTrial << "\n";
   std::cout << "  AVG_THETA/ baseline = " << r_avgTheta << "\n";
   std::cout << "  AVG_VABS / baseline = " << r_avgVabs  << "\n";
   std::cout << "  AVG_TCORR/ baseline = " << r_avgTcorr << "\n";
   std::cout << "  NONFAULT_BOTH_SYM / baseline        = "
             << r_nonfault_both_sym << "\n";
   std::cout << "  FREEZE-A + NONFAULT_BOTH_SYM / base = "
             << r_freezeA_nfBothSym << "\n";
   std::cout << "  DT/2     / baseline = " << r_dt2 << "\n";
   std::cout << "  DT/4     / baseline = " << r_dt4 << "\n";

   // Round-13C Patch 2 verdict interpretation:
   //  r_nonfault_both_sym < 0.2  → non-fault sign/orientation is a
   //                                dominant full-loop driver
   //  r_nonfault_both_sym ≈ 1.0  → SXZ cleanup is real but not the
   //                                main long-time pepper source
   //  r_freezeA_nfBothSym should still be 0.0 (FREEZE-A bypasses
   //                                the fault entirely).
   std::cout << "\n=== Round-13C NONFAULT_BOTH_SYM verdict ===\n";
   if (r_nonfault_both_sym < 0.2)
   {
      std::cout << "  → NONFAULT_BOTH_SYM reduces pepper by > 80%.\n"
                << "    Non-fault sign/orientation is a DOMINANT\n"
                << "    full-loop driver.  Extend the hook to shared\n"
                << "    faces next.\n";
   }
   else if (r_nonfault_both_sym < 0.9)
   {
      std::cout << "  → NONFAULT_BOTH_SYM reduces pepper by "
                << std::fixed << std::setprecision(1)
                << (1.0 - r_nonfault_both_sym) * 100.0 << "%.\n"
                << "    Non-fault sign symmetrization helps partially;\n"
                << "    remaining residual is elsewhere (likely the\n"
                << "    boundary SXY floor identified in Round 13B).\n";
   }
   else
   {
      std::cout << "  → NONFAULT_BOTH_SYM ratio ≈ 1.0.\n"
                << "    The SXZ-cleanup proved by Gate 4e does not\n"
                << "    translate to a visible long-time pepper drop.\n"
                << "    Target the boundary SXY floor directly.\n";
   }

   // Round-12 Patch 3 dedicated AVG-stage table (mode | tau1_corr | ratio).
   std::cout << "\n=== Round-12 AVG stage table ===\n"
             << std::scientific << std::setprecision(6);
   std::cout << "  mode      | tau1_corr    | ratio\n"
             << "  ----------|--------------|-------\n";
   std::cout << "  BASELINE  | " << baseline.tau1_corr << " | "
             << std::fixed << std::setprecision(4) << 1.0000 << "\n"
             << std::scientific << std::setprecision(6);
   std::cout << "  AVG_TRIAL | " << avgTrial.tau1_corr << " | "
             << std::fixed << std::setprecision(4) << r_avgTrial << "\n"
             << std::scientific << std::setprecision(6);
   std::cout << "  AVG_THETA | " << avgTheta.tau1_corr << " | "
             << std::fixed << std::setprecision(4) << r_avgTheta << "\n"
             << std::scientific << std::setprecision(6);
   std::cout << "  AVG_VABS  | " << avgVabs.tau1_corr  << " | "
             << std::fixed << std::setprecision(4) << r_avgVabs  << "\n"
             << std::scientific << std::setprecision(6);
   std::cout << "  AVG_TCORR | " << avgTcorr.tau1_corr << " | "
             << std::fixed << std::setprecision(4) << r_avgTcorr << "\n";

   std::cout << std::scientific << std::setprecision(3);
   std::cout << "\n=== dt-scaling fit (tau1_corr vs dt) ===\n";
   // If tau1_corr ∝ dt^p, then p = log(ratio) / log(dt_ratio).
   // Use dt and dt/4 to get a cleaner estimate.
   if (dt_quart.tau1_corr > 0 && baseline.tau1_corr > 0)
   {
      const real_t p = std::log(baseline.tau1_corr / dt_quart.tau1_corr)
                       / std::log(4.0);
      std::cout << "  tau1_corr ∝ dt^" << std::fixed << std::setprecision(3)
                << p << " (fit from baseline / dt/4 ratio)\n";
      std::cout << "  p > 0.5  → pepper vanishes as dt → 0 (sub-resolved physics)\n";
      std::cout << "  p ≈ 0    → dt-insensitive (code bug)\n";
   }

   std::cout << "\n=== Verdict per R-005 matrix (FREEZE-C + dt) ===\n";
   // Thresholds:
   //   FREEZE-* "closes" if ratio < 0.1 (90%+ reduction).
   //   DT "scales" if DT/4 ratio < 0.5 (spread drops with dt).
   const bool freezeA_closes = (r_freezeA < 0.1);
   const bool freezeB_closes = (r_freezeB < 0.1);
   const bool freezeC_closes = (r_freezeC < 0.1);
   const bool dt_scales = (r_dt4 < 0.5);

   std::cout << "  FREEZE-A closes pepper? "
             << (freezeA_closes ? "YES" : "NO") << "\n";
   std::cout << "  FREEZE-B closes pepper? "
             << (freezeB_closes ? "YES" : "NO") << "\n";
   std::cout << "  FREEZE-C closes pepper? "
             << (freezeC_closes ? "YES" : "NO") << "\n";
   std::cout << "  dt refinement scales pepper? "
             << (dt_scales ? "YES (spread decreases)" : "NO (dt-insensitive)")
             << "\n\n";

   // Round-11 decision tree: which of {fault friction, fault
   // feedback, MFEM wave update} generates the first asymmetry.
   //
   // FREEZE-A disables the fault friction solve AND the fault-to-
   // bulk flux correction (Q_imp = Q_avg means zero jump in imposed
   // state).  If pepper persists, the amplifier is entirely within
   // the DG wave update — interior non-fault faces, boundary faces,
   // or the ADER volume-term + mass-inverse chain.
   //
   // FREEZE-B freezes the full friction output (no feedback into
   // subsequent DOFData reads) but the Riemann flux correction still
   // runs using the cached tau*_corr / sigma_n_corr fields.  If
   // FREEZE-A closes but FREEZE-B does not, the amplifier is the
   // live re-coupling of friction outputs back to DOFData.  If both
   // close, the amplifier requires the friction solve to run.
   //
   // FREEZE-C holds only psi.  If it closes but FREEZE-B does not,
   // the ψ-feedback path dominates.
   std::cout << "=== Round-11 decision ===\n";
   if (freezeA_closes)
   {
      std::cout << "  → FREEZE-A closes pepper:\n"
                << "    Fault friction + flux correction together cause "
                << "the asymmetry.\n"
                << "    Disambiguate with FREEZE-B vs FREEZE-C:\n";
      if (freezeC_closes && !freezeB_closes)
      {
         std::cout << "      FREEZE-C closes but FREEZE-B does not → "
                   << "ψ feedback is the path.\n";
      }
      else if (freezeB_closes)
      {
         std::cout << "      FREEZE-B closes → full friction-state "
                   << "re-coupling is the path.\n";
      }
      else
      {
         std::cout << "      neither FREEZE-B nor FREEZE-C close → "
                   << "friction SOLVE itself is the seed; test "
                   << "atan2 / Brent / trial-traction.\n";
      }
   }
   else
   {
      std::cout << "  → FREEZE-A does NOT close pepper:\n"
                << "    The MFEM wave update generates the first "
                << "nonuniformity even with the fault locked.\n"
                << "    Amplifier is OUTSIDE the fault — in one of:\n"
                << "      (a) interior non-fault face flux,\n"
                << "      (b) boundary face flux,\n"
                << "      (c) ADER volume term / mass inverse.\n"
                << "    Use test_adjacent_triangle_fault_first_step_audit\n"
                << "    with SEAS_TEST_SKIP_BOUNDARY=1 and\n"
                << "    SEAS_TEST_SKIP_INTERIOR_NONFAULT=1 to isolate.\n";
   }

   std::cout << "\n=== Verdict per R-005 (ψ + dt) ===\n";
   if (freezeC_closes && dt_scales)
   {
      std::cout << "  → ψ-evolution IS amplifier AND dt-sensitive.\n"
                << "    Round-next target: friction/dieterich_ruina.hpp::\n"
                << "                       UpdateStateAnalytic precision.\n";
   }
   else if (freezeC_closes && !dt_scales)
   {
      std::cout << "  → ψ IS amplifier BUT dt-insensitive.\n"
                << "    Round-next target: CODE BUG in ψ-feedback coupling.\n";
   }
   else if (!freezeC_closes && dt_scales)
   {
      std::cout << "  → ψ NOT amplifier, but dt-sensitive — sub-resolved\n"
                << "    phenomenon in wave-kernel or ADER stage interleave.\n";
   }
   else
   {
      std::cout << "  → ψ NOT amplifier AND dt-insensitive.\n"
                << "    FREEZE-A/B verdict above dominates.\n";
   }

   // Round-12 AVG verdict: identify the FIRST closing averaging mode
   // along the pipeline (Trial → Theta → Vabs → Tcorr).
   const bool avgTrial_closes = (r_avgTrial < 0.1);
   const bool avgTheta_closes = (r_avgTheta < 0.1);
   const bool avgVabs_closes  = (r_avgVabs  < 0.1);
   const bool avgTcorr_closes = (r_avgTcorr < 0.1);

   std::cout << "\n=== Round-12 AVG verdict ===\n";
   std::cout << "  AVG_TRIAL closes pepper? "
             << (avgTrial_closes ? "YES" : "NO") << "\n";
   std::cout << "  AVG_THETA closes pepper? "
             << (avgTheta_closes ? "YES" : "NO") << "\n";
   std::cout << "  AVG_VABS  closes pepper? "
             << (avgVabs_closes  ? "YES" : "NO") << "\n";
   std::cout << "  AVG_TCORR closes pepper? "
             << (avgTcorr_closes ? "YES" : "NO") << "\n\n";

   if (avgTrial_closes)
   {
      std::cout << "  → AVG_TRIAL closes first:\n"
                << "    Seed enters through per-QP trial traction /\n"
                << "    trace sampling.  Round-next should target\n"
                << "    ComputeTrialTraction inputs: I_plus/I_minus\n"
                << "    reconstruction via shape · I, Loc1/Loc2 trace,\n"
                << "    and the rotation-to-canonical step.\n";
   }
   else if (avgTheta_closes)
   {
      std::cout << "  → AVG_THETA is first closing mode:\n"
                << "    The nonlinear magnitude map Θ = √(τ1² + τ2²)\n"
                << "    amplifies per-QP ULP asymmetry into O(1)-scale\n"
                << "    asymmetric Θ.  Round-next should test whether\n"
                << "    replacing sqrt with a rescaled form removes\n"
                << "    the amplification.\n";
   }
   else if (avgVabs_closes)
   {
      std::cout << "  → AVG_VABS is first closing mode:\n"
                << "    Root-solver (Brent on Eq. 8) is the first\n"
                << "    amplifier.  Round-next should instrument\n"
                << "    per-QP Brent iteration counts and residuals.\n";
   }
   else if (avgTcorr_closes)
   {
      std::cout << "  → AVG_TCORR is first (and only) closing mode:\n"
                << "    Imposed-state / Riemann injection is the first\n"
                << "    stage where the asymmetry becomes causally\n"
                << "    relevant.  Round-next should target Eq. 11/12\n"
                << "    construction in BuildImposedState and its use\n"
                << "    inside flux_.Interior.\n";
   }
   else
   {
      std::cout << "  → NONE of the four averaging modes close pepper.\n"
                << "    Asymmetry persists DOWNSTREAM of corrected\n"
                << "    traction — contradicts Round-11's FREEZE-A=0\n"
                << "    evidence.  Force re-check of the harness:\n"
                << "    verify the hook in wave_operator.inl fires\n"
                << "    (print a debug message per face at q==0) and\n"
                << "    confirm FREEZE-A still zeroes pepper on the\n"
                << "    same build.\n";
   }

   return 0;
}
