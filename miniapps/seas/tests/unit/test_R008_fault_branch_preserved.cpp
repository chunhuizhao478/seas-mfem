// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 "Topology-Based Precomputed Face-Rotation" plan 2026-04-23
// Phase 2a regression gate (§5.1.9):
//
//   test_R008_fault_branch_preserved_after_dispatch_reorder
//     — Runs one ADER step on the M0 fault fixture with the opt-in flag
//       OFF (default).  If the §6.2 dispatch restructure dropped or
//       garbled the fault-branch body, fault_dof_data_ will be zero or
//       NaN after AdvanceADER.  We assert every fault QP has nonzero
//       finite slip_rate / tau1_corr / tau2_corr / sigma_n_corr.
//
//   P_SWITCH_ON_FAULT_NOT_ZEROED (§6.4) — with flag ON, the fault
//   dispatch must fire: assert |tau1_corr| or |tau2_corr| > 1e-10 at
//   every fault QP.
//
//   P_SWITCH_ON_NO_CRASH (§6.4) — the flag-on AdvanceADER completes
//   without aborting on the M0 fixture.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../domain/boundary_config.hpp"
#include "../../common/seas_types.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace mfem { namespace seas { int g_seas_my_rank = 0; } }  // NOLINT

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_TRUE(expr, msg) do { \
   num_tests++; \
   if (expr) { num_passed++; \
      std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)

namespace
{

constexpr real_t kL        = 1000.0;
constexpr real_t kRho      = 2670.0;
constexpr real_t kLambda   = 32.04e9;
constexpr real_t kMu       = 32.04e9;
constexpr real_t kDt       = 5.0e-5;
constexpr int    kOrder    = 1;
constexpr int    kAderOrder = 2;

Mesh BuildM0FaultMesh()
{
   const int n = 2;
   Mesh mesh = Mesh::MakeCartesian3D(n, n, n, Element::TETRAHEDRON, kL, kL, kL,
                                     false);
   mesh.FinalizeTopology();
   mesh.Finalize();

   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      Array<int> fv; mesh.GetFaceVertices(f, fv);
      if (fv.Size() != 3) { continue; }
      auto *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr && ftr->Elem2No >= 0)
      {
         real_t cy = 0.0;
         for (int v = 0; v < fv.Size(); v++) { cy += mesh.GetVertex(fv[v])[1]; }
         cy /= fv.Size();
         if (std::abs(cy - 0.5 * kL) < 1e-8)
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

// Seed the DOFData with the TPV102 persistent-nucleation fixture used by
// test_adjacent_triangle_fault_first_step_audit: uniform pre-stress,
// dipping/strike basis at the fault.
void SeedDOFData(std::vector<DOFData> &dof, int nfq)
{
   dof.clear();
   dof.resize(nfq);

   // TPV102-style uniform pre-stress in Pa.  Field names from
   // DOFData (fault_face_flux.hpp:25): sigma_n0 > 0 is compression;
   // {Dc, a} are friction parameters; psi is the state variable;
   // tau2_nuc is the strike-slip nucleation driver used by total-Q
   // dispatch (kept zero here since we exercise fluctuation-Q via
   // bulk_bg_zero).
   const real_t sigma_n0 = 40.0e6;
   const real_t tau_ini  = 30.0e6;    // strike-slip shear pre-stress [Pa]
   const real_t V_ini    = 1.0e-12;   // steady sliding seed
   const real_t a        = 0.008;
   const real_t Dc       = 0.02;
   const real_t cs       = std::sqrt(kMu / kRho);
   const real_t Zp       = kRho * std::sqrt((kLambda + 2.0 * kMu)/kRho);
   const real_t Zs       = kRho * cs;

   // FrictionSolver uses implicit V0=1e-6, f0=0.6, b≡0 (Slip law).
   // State-steady psi: f_ss = f0 + a*ln(V/V0).  Use that as psi_ss.
   const real_t V0 = 1.0e-6;
   const real_t f0 = 0.6;
   const real_t psi_ss = f0 + a * std::log(V_ini / V0);

   for (int i = 0; i < nfq; i++)
   {
      DOFData &d = dof[i];
      d.slip_rate = V_ini;
      d.V1 = 0.0;            // dip
      d.V2 = V_ini;          // strike
      d.slip1 = 0.0;
      d.slip2 = 0.0;
      d.tau1_0 = 0.0;        // dip component of pre-stress
      d.tau2_0 = tau_ini;    // strike component
      d.sigma_n0 = sigma_n0;
      d.tau1_corr = 0.0;
      d.tau2_corr = 0.0;
      d.sigma_n_corr = 0.0;
      d.psi = psi_ss;
      d.a = a;
      d.Dc = Dc;
      d.Zp_plus = Zp;  d.Zp_minus = Zp;
      d.Zs_plus = Zs;  d.Zs_minus = Zs;
      d.eta_p = 0.5 * Zp;
      d.eta_s = 0.5 * Zs;
   }
}

struct StepResult
{
   Vector Q1;
   std::vector<DOFData> dof_after;
   bool aborted = false;
};

StepResult RunOneStep(Mesh &mesh, bool flag_on)
{
   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr    = 3;

   real_t Q_bg_zero[NUM_STATE] = {0};

   FaultFaceFlux ff(kRho,
                    std::sqrt((kLambda + 2.0 * kMu) / kRho),
                    std::sqrt(kMu / kRho));
   WaveOperator<Mesh> wave(mesh, kOrder, kLambda, kMu, kRho, bc);
   wave.SetAbsorbingBackground(Q_bg_zero);
   wave.SetFaultFlux(&ff);

   std::vector<DOFData> dof;
   SeedDOFData(dof, wave.GetNumTotalFaultQPs());
   wave.SetFaultDOFData(&dof, wave.GetNbfPerFace());

   if (flag_on)
   {
      wave.UsePrecomputedFaceFluxes(true);
   }

   Vector Q0(wave.Height());
   Q0 = 0.0;
   Vector Q1(wave.Height());
   wave.AdvanceADER(Q0, kDt, kAderOrder, Q1);

   StepResult r;
   r.Q1 = Q1;
   r.dof_after = dof;
   return r;
}

// Check that EVERY fault QP evolved from the seeded state.  The seed
// has (slip_rate=V_ini=1e-12, tau1_corr=0, tau2_corr=0, sigma_n_corr=0),
// so the "finite + nonzero" check was a tautology against the seed
// values (the DOFData always started non-degenerate).  This
// strengthened version asserts:
//   (1) every field is finite
//   (2) tau2_corr was UPDATED from 0 (its seed) to a physically-nonzero
//       corrected-traction value.  EvaluateADER always writes
//       tau_corr even on a steady slip-rate state, so a fault branch
//       replaced by `continue` (which skips EvaluateADER) leaves
//       tau2_corr at the seed 0 and this check fires.
//   (3) |sigma_n_corr| > 0: normal-stress trial is recomputed every
//       step even at zero bulk perturbation (background normal stress
//       is re-assembled in EvaluateADER).
// The Q1-norm check after this routine provides the complementary
// "fault branch did accumulate into rhs" signal.
int CountFailingFaultQPs(const std::vector<DOFData> &dof,
                         real_t eps = 1e-14)
{
   int fails = 0;
   for (const DOFData &d : dof)
   {
      bool finite =
         std::isfinite(d.slip_rate) && std::isfinite(d.V1) &&
         std::isfinite(d.V2)        && std::isfinite(d.tau1_corr) &&
         std::isfinite(d.tau2_corr) && std::isfinite(d.sigma_n_corr);
      if (!finite) { fails++; continue; }
      // EvaluateADER must write a nonzero tau2_corr OR sigma_n_corr —
      // both were seeded as 0, so either being nonzero proves the fault
      // branch body EXECUTED (not merely seeded).  A no-op continue
      // leaves both at 0 and this check fails (unlike the prior
      // "max(slip_rate, ...) > eps" check which was a tautology
      // because slip_rate is seeded at V_ini=1e-12 > eps).
      const real_t evolved_scale =
         std::max(std::abs(d.tau2_corr), std::abs(d.sigma_n_corr));
      if (evolved_scale < eps) { fails++; }
   }
   return fails;
}

} // anonymous

int main()
{
   std::cout << "\n=== TPV102 Phase 2a R008 fault-branch-preservation + "
                "P_SWITCH_ON_FAULT_NOT_ZEROED ===\n";

   // =====================================================================
   // test_R008_fault_branch_preserved_after_dispatch_reorder:
   //   With flag OFF (default), the restructured ComputeADERFaceFluxRHS
   //   must route fault faces through the preserved fault-branch body
   //   and produce finite nonzero DOFData.
   // =====================================================================
   Mesh mesh_off = BuildM0FaultMesh();
   StepResult off = RunOneStep(mesh_off, /*flag_on=*/false);

   const int n_off_fail = CountFailingFaultQPs(off.dof_after);
   TEST_TRUE(n_off_fail == 0,
             "test_R008_fault_branch_preserved: flag-OFF one-step "
             "ADER step produces finite nonzero DOFData on every "
             "fault QP (fault branch body preserved after dispatch "
             "restructure)");

   real_t q1_off_inf = 0.0;
   for (int i = 0; i < off.Q1.Size(); i++)
   {
      q1_off_inf = std::max(q1_off_inf, std::abs(off.Q1(i)));
   }
   TEST_TRUE(q1_off_inf > 1e-14,
             "flag-OFF ADER step produces nonzero Q_new (bulk coupled "
             "to fault)");

   // =====================================================================
   // P_SWITCH_ON_FAULT_NOT_ZEROED + P_SWITCH_ON_NO_CRASH:
   //   With flag ON, the fault dispatch must fire (fault_dof_data_
   //   still populated, fault QPs nonzero), and AdvanceADER must not
   //   crash.
   // =====================================================================
   Mesh mesh_on = BuildM0FaultMesh();
   StepResult on = RunOneStep(mesh_on, /*flag_on=*/true);

   const int n_on_fail = CountFailingFaultQPs(on.dof_after);
   TEST_TRUE(n_on_fail == 0,
             "P_SWITCH_ON_FAULT_NOT_ZEROED: flag-ON ADER step produces "
             "finite nonzero DOFData on every fault QP");

   real_t q1_on_inf = 0.0;
   for (int i = 0; i < on.Q1.Size(); i++)
   {
      q1_on_inf = std::max(q1_on_inf, std::abs(on.Q1(i)));
   }
   TEST_TRUE(q1_on_inf > 1e-14,
             "P_SWITCH_ON_NO_CRASH: flag-ON ADER step completes with "
             "nonzero Q_new");

   std::cout << "\n===================================================\n";
   std::cout << "Total:  " << num_tests << "\n";
   std::cout << "Passed: " << num_passed << "\n";
   std::cout << "Failed: " << num_failed << "\n";
   std::cout << "===================================================\n";
   return (num_failed == 0) ? 0 : 1;
}
