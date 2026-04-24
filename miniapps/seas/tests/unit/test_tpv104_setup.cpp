// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// Unit tests for dynamic/tpv104_setup.hpp (§4.10 Step 3 gates
// T_TPV104_SETUP_1..5).

#include "test_macros.hpp"
#include "../../dynamic/tpv104_setup.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../fault/fault_basis.hpp"
#include "../../config/tpv104_params.hpp"

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

// ---------------------------------------------------------------------------
// Helper: build a synthetic fault-coordinate list covering the TPV104
// VW core, transition, and strengthening regions.  Each QP is at
// (along_strike, y=0, -depth).
// ---------------------------------------------------------------------------
static void BuildSyntheticFaultCoords(std::vector<Vector> &fault_coords)
{
   // Grid of 5×5 QPs: x2 ∈ {-20, -12, 0, 12, 20} km,
   //                   z  ∈ {  0,  3,  7.5, 12, 20} km.
   const std::vector<real_t> xs = {-20e3, -12e3, 0.0, 12e3, 20e3};
   const std::vector<real_t> zs = {  0.0,  3e3, 7.5e3, 12e3, 20e3};

   fault_coords.clear();
   for (real_t x : xs)
   {
      for (real_t z : zs)
      {
         Vector c(3);
         c(0) = x;
         c(1) = 0.0;
         c(2) = -z;   // fault_coords uses z < 0 = depth
         fault_coords.push_back(c);
      }
   }
}

// ---------------------------------------------------------------------------
// T_TPV104_SETUP_1 — ψ init anchor at every fault QP.
// For every DOF, |ψ_init − ψ_anchor(a_i)| < 1e-8.  The anchor is
// ComputeInitialPsiTPV104(a_i), which at a_in = 0.01 returns 5.6359184e-1
// (T_TPV104_P_1).
// ---------------------------------------------------------------------------
void TestPsiInitPerQP()
{
   std::cout << "\n[T_TPV104_SETUP_1] ψ_init at every QP\n";

   std::vector<Vector> fault_coords;
   BuildSyntheticFaultCoords(fault_coords);
   const int ndof = static_cast<int>(fault_coords.size());

   std::vector<DOFData> dof_data;
   InitializeFaultDOFs_TPV104(dof_data, ndof, fault_coords);

   TEST_ASSERT(static_cast<int>(dof_data.size()) == ndof,
               "dof_data resized to ndof");

   int ok = 0;
   for (int i = 0; i < ndof; ++i)
   {
      const real_t a_i      = dof_data[i].a;
      const real_t psi_ref  = ComputeInitialPsiTPV104(a_i);
      const real_t psi_got  = dof_data[i].psi;
      if (std::abs(psi_got - psi_ref) < 1e-12) { ++ok; }
      else
      {
         std::cerr << "  FAIL i=" << i
                   << " a=" << a_i
                   << " psi_ref=" << psi_ref
                   << " psi_got=" << psi_got
                   << " diff=" << (psi_got - psi_ref) << "\n";
      }
   }
   TEST_ASSERT(ok == ndof,
               "every DOF has psi == ComputeInitialPsiTPV104(a_i)");

   // Additional anchor: at a hypocenter-depth QP in the VW core
   // (along_strike=0, down_dip=7.5 km), a = a_in = 0.01 and ψ ≈ 0.5636.
   int hypo = -1;
   for (int i = 0; i < ndof; ++i)
   {
      if (fault_coords[i](0) == 0.0 && std::abs(fault_coords[i](2)) == 7.5e3)
      {
         hypo = i;
         break;
      }
   }
   TEST_ASSERT(hypo >= 0, "synthetic grid contains the hypocenter QP");
   if (hypo >= 0)
   {
      TEST_NEAR(dof_data[hypo].a, TPV104Params::a_in, 1e-15,
                "hypo QP has a = a_in = 0.01");
      TEST_NEAR(dof_data[hypo].psi, 5.6359184e-01, 1e-8,
                "hypo QP has ψ_init ≈ 5.6359184e-1 (§4.1 anchor)");
   }
}

// ---------------------------------------------------------------------------
// T_TPV104_SETUP_2 — V_w[i] side-channel per-QP.
// Inside VW core (|x| < 15 km, 0 < z < 15 km): V_w = 0.1.
// Outside:                                     V_w = 1.0.
// Transition zone: strictly between the two.
// ---------------------------------------------------------------------------
void TestVwSideChannel()
{
   std::cout << "\n[T_TPV104_SETUP_2] V_w[i] side-channel population\n";

   std::vector<Vector> fault_coords;
   BuildSyntheticFaultCoords(fault_coords);
   std::vector<real_t> V_w;
   PopulateVwSideChannel_TPV104(V_w, fault_coords);

   TEST_ASSERT(V_w.size() == fault_coords.size(),
               "V_w.size() == fault_coords.size()");

   int vw_core_ok = 0;
   int vw_outside_ok = 0;
   int vw_core_count = 0;
   int vw_outside_count = 0;
   for (size_t i = 0; i < V_w.size(); ++i)
   {
      const real_t x = fault_coords[i](0);
      const real_t z = std::abs(fault_coords[i](2));

      // Deep VW core: |x| ≤ L_s = 15 km AND |z - hypo_down_dip| ≤ W/2 = 7.5 km.
      const bool in_vw_core = (std::abs(x) <= TPV104Params::Ls) &&
                              (std::abs(z - TPV104Params::hypo_down_dip)
                                 <= TPV104Params::W / 2.0);
      // Strengthening: outside the boxcar + transition on both axes.
      const bool in_strengthening =
         (std::abs(x) >= TPV104Params::Ls + TPV104Params::ws) ||
         (std::abs(z - TPV104Params::hypo_down_dip)
            >= TPV104Params::W / 2.0 + TPV104Params::w);

      if (in_vw_core)
      {
         ++vw_core_count;
         if (std::abs(V_w[i] - TPV104Params::V_w_in) < 1e-15) { ++vw_core_ok; }
         else
         {
            std::cerr << "  FAIL VW-core i=" << i << " x=" << x << " z=" << z
                      << " V_w=" << V_w[i] << " expected=" << TPV104Params::V_w_in
                      << "\n";
         }
      }
      else if (in_strengthening)
      {
         ++vw_outside_count;
         if (std::abs(V_w[i] - TPV104Params::V_w_out) < 1e-15) { ++vw_outside_ok; }
         else
         {
            std::cerr << "  FAIL outside i=" << i << " x=" << x << " z=" << z
                      << " V_w=" << V_w[i] << " expected=" << TPV104Params::V_w_out
                      << "\n";
         }
      }
      else
      {
         // Transition zone — must be strictly between V_w_in and V_w_out.
         TEST_ASSERT(V_w[i] >= TPV104Params::V_w_in
                     && V_w[i] <= TPV104Params::V_w_out,
                     "transition V_w ∈ [V_w_in, V_w_out]");
      }
   }
   TEST_ASSERT(vw_core_count > 0 && vw_core_ok == vw_core_count,
               "every VW-core QP has V_w = 0.1 m/s");
   TEST_ASSERT(vw_outside_count > 0 && vw_outside_ok == vw_outside_count,
               "every strengthening QP has V_w = 1.0 m/s");
}

// ---------------------------------------------------------------------------
// T_TPV104_SETUP_3 — coordinate-mapping invariant ([C1] + [C2]).
// For the TPV104 reference fault (n = (0, -1, 0), up = (0, 0, 1)), the
// FaultBasis oriented frame must produce:
//   tangent1 parallel to (0, 0, -1)  (dip — downward)
//   tangent2 parallel to (+1, 0, 0)  (strike — along-x)
//
// If this fails, BP5's `FaultBasis` has been modified — halt ([C2]).
// ---------------------------------------------------------------------------
void TestFaultBasisConvention()
{
   std::cout << "\n[T_TPV104_SETUP_3] FaultBasis oriented-frame convention\n";

   // Raw face normal consistent with the TPV104 ref (y = 0 fault plane,
   // CalcOrtho → (0, -1, 0) or (0, +1, 0) depending on orientation).
   // Use (0, -1, 0) directly — the sign flip inside ComputeOrientedFrame
   // handles any ref-normal disagreement.
   Vector n_raw(3);
   n_raw(0) = 0.0; n_raw(1) = -1.0; n_raw(2) = 0.0;
   Vector ref_normal(3);
   ref_normal(0) = 0.0; ref_normal(1) = -1.0; ref_normal(2) = 0.0;
   Vector up(3);
   up(0) = 0.0; up(1) = 0.0; up(2) = 1.0;

   real_t normal[3], tangent1[3], tangent2[3];
   bool sign_flipped = false;
   real_t nl = 0.0;

   FaultBasis::ComputeOrientedFrame(n_raw, /*dim=*/3, ref_normal, up,
                                    normal, tangent1, tangent2,
                                    sign_flipped, nl);

   TEST_NEAR(nl, 1.0, 1e-15, "raw normal length = 1");
   TEST_ASSERT(!sign_flipped, "ref-aligned input gives sign_flipped = false");

   // tangent1 (dip, downward) must be ∥ (0, 0, -1).
   TEST_NEAR(tangent1[0],  0.0, 1e-15, "tangent1.x = 0");
   TEST_NEAR(tangent1[1],  0.0, 1e-15, "tangent1.y = 0");
   TEST_NEAR(tangent1[2], -1.0, 1e-15, "tangent1.z = -1 (dip-down)");

   // tangent2 (strike) must be ∥ (+1, 0, 0).
   TEST_NEAR(tangent2[0],  1.0, 1e-15, "tangent2.x = +1 (along-strike)");
   TEST_NEAR(tangent2[1],  0.0, 1e-15, "tangent2.y = 0");
   TEST_NEAR(tangent2[2],  0.0, 1e-15, "tangent2.z = 0");

   // Normal stays aligned with the ref.
   TEST_NEAR(normal[0],  0.0, 1e-15, "normal.x = 0");
   TEST_NEAR(normal[1], -1.0, 1e-15, "normal.y = -1");
   TEST_NEAR(normal[2],  0.0, 1e-15, "normal.z = 0");
}

// ---------------------------------------------------------------------------
// T_TPV104_SETUP_4 — pre-stress layout (BP5 / Tandem canonical frame).
// At every DOF:
//   tau1_0    == 0       (no dip pre-stress)
//   tau2_0    == 4e7     (along-strike)
//   sigma_n0  == 1.2e8   (positive compression)
// A regression that wrote under the GodunovFlux t1=strike convention
// would swap tau1_0 ↔ tau2_0 and trip this test.
// ---------------------------------------------------------------------------
void TestPreStressLayout()
{
   std::cout << "\n[T_TPV104_SETUP_4] pre-stress layout (BP5 convention)\n";

   std::vector<Vector> fault_coords;
   BuildSyntheticFaultCoords(fault_coords);
   const int ndof = static_cast<int>(fault_coords.size());

   std::vector<DOFData> dof_data;
   InitializeFaultDOFs_TPV104(dof_data, ndof, fault_coords);

   int ok_t1 = 0, ok_t2 = 0, ok_sn = 0;
   int ok_vslip = 0, ok_Vinit = 0;
   int ok_L = 0;
   int ok_nuc = 0;
   for (int i = 0; i < ndof; ++i)
   {
      if (dof_data[i].tau1_0 == 0.0) { ++ok_t1; }
      if (dof_data[i].tau2_0 == TPV104Params::tau_ini) { ++ok_t2; }
      if (dof_data[i].sigma_n0 == TPV104Params::sigma_n) { ++ok_sn; }
      if (dof_data[i].V1 == 0.0 && dof_data[i].slip1 == 0.0) { ++ok_vslip; }
      if (dof_data[i].V2 == TPV104Params::V_ini) { ++ok_Vinit; }
      if (dof_data[i].Dc == TPV104Params::L) { ++ok_L; }
      if (dof_data[i].sigma_n_nuc == 0.0
          && dof_data[i].tau1_nuc == 0.0
          && dof_data[i].tau2_nuc == 0.0) { ++ok_nuc; }
   }
   TEST_ASSERT(ok_t1 == ndof,
               "every DOF has tau1_0 == 0 (no dip pre-stress)");
   TEST_ASSERT(ok_t2 == ndof,
               "every DOF has tau2_0 == 40 MPa (along-strike)");
   TEST_ASSERT(ok_sn == ndof,
               "every DOF has sigma_n0 == 120 MPa (positive compression)");
   TEST_ASSERT(ok_vslip == ndof,
               "every DOF has V1 == 0 and slip1 == 0 (no dip slip)");
   TEST_ASSERT(ok_Vinit == ndof,
               "every DOF has V2 == V_ini = 1e-16 m/s");
   TEST_ASSERT(ok_L == ndof,
               "every DOF has Dc == L = 0.4 m (TPV104 critical slip)");
   TEST_ASSERT(ok_nuc == ndof,
               "every DOF has nucleation channels zeroed at init");

   // Numerical anchors.
   TEST_NEAR(TPV104Params::tau_ini, 4.0e7, 1.0,  "τ_ini anchor");
   TEST_NEAR(TPV104Params::sigma_n, 1.2e8, 1.0,  "σ_n anchor");
   TEST_NEAR(TPV104Params::L,       0.4,   1e-15, "L anchor");
   TEST_NEAR(TPV104Params::V_ini,   1.0e-16, 1e-30, "V_ini anchor");

   // Initial corrected traction = background.
   int ok_corr = 0;
   for (int i = 0; i < ndof; ++i)
   {
      if (dof_data[i].tau1_corr == 0.0
          && dof_data[i].tau2_corr == TPV104Params::tau_ini
          && dof_data[i].sigma_n_corr == TPV104Params::sigma_n)
      {
         ++ok_corr;
      }
   }
   TEST_ASSERT(ok_corr == ndof,
               "every DOF has initial corrected traction = background");
}

// ---------------------------------------------------------------------------
// T_TPV104_SETUP_5 — no `InitializeStateTotal_TPV104` symbol.
// Static grep on `dynamic/tpv104_setup.hpp` must return zero hits for
// `InitializeStateTotal_TPV104` or `ApplyNucleationTotalPrestress_TPV104`.
// Enforces §3.10 directive (fluctuation-Q only).
// ---------------------------------------------------------------------------
void TestNoTotalQSymbols()
{
   std::cout << "\n[T_TPV104_SETUP_5] no total-Q symbols (§3.10)\n";

   auto locate = [](const std::string &tail) -> std::string
   {
      const std::vector<std::string> prefixes = {"", "../", "../../"};
      for (const auto &p : prefixes)
      {
         std::ifstream f(p + tail);
         if (f.is_open()) { return p + tail; }
      }
      return "";
   };
   const std::string path = locate("dynamic/tpv104_setup.hpp");
   TEST_ASSERT(!path.empty(), "locate dynamic/tpv104_setup.hpp");
   if (path.empty()) { return; }

   std::ifstream in(path);
   std::ostringstream oss;
   oss << in.rdbuf();
   const std::string body = oss.str();

   // Scan for the forbidden symbols at the top-level identifier scope,
   // ignoring their mention inside comments (the header has a prose
   // block explaining why we do NOT ship total-Q helpers).
   auto count_non_comment = [&body](const std::string &needle) -> int
   {
      int hits = 0;
      size_t pos = 0;
      while ((pos = body.find(needle, pos)) != std::string::npos)
      {
         const size_t line_start = body.rfind('\n', pos);
         const size_t line_begin = (line_start == std::string::npos)
                                   ? 0 : line_start + 1;
         const size_t line_end = body.find('\n', pos);
         const std::string line = body.substr(line_begin,
                                              line_end - line_begin);
         const size_t cs = line.find("//");
         const bool in_comment = cs != std::string::npos
                                 && cs <= (pos - line_begin);
         if (!in_comment) { ++hits; }
         pos += needle.size();
      }
      return hits;
   };

   TEST_ASSERT(count_non_comment("InitializeStateTotal_TPV104") == 0,
               "no InitializeStateTotal_TPV104 symbol at code scope");
   TEST_ASSERT(count_non_comment("ApplyNucleationTotalPrestress_TPV104") == 0,
               "no ApplyNucleationTotalPrestress_TPV104 symbol at code scope");
}

// ---------------------------------------------------------------------------
// Helper sanity tests.
// ---------------------------------------------------------------------------
void TestInitializeStateZerosQ()
{
   std::cout << "\n[T_TPV104_SETUP_6] InitializeState_TPV104 sets Q = 0\n";

   Vector Q;
   const int ndof_total = 100;
   InitializeState_TPV104(Q, ndof_total);
   TEST_ASSERT(Q.Size() == NUM_STATE * ndof_total,
               "Q resized to NUM_STATE * ndof_total");
   real_t max_abs = 0.0;
   for (int i = 0; i < Q.Size(); ++i)
   {
      max_abs = std::max(max_abs, std::abs(Q[i]));
   }
   TEST_NEAR(max_abs, 0.0, 1e-30, "Q = 0 after InitializeState_TPV104");
}

void TestDefaultStations()
{
   std::cout << "\n[T_TPV104_SETUP_7] DefaultStations_TPV104 builds 9 entries\n";

   const auto stations = DefaultStations_TPV104();
   TEST_ASSERT(stations.size() == 9,
               "DefaultStations_TPV104 returns 9 stations");

   bool found_hypo = false;
   for (const auto &s : stations)
   {
      if (s.along_strike == 0.0 && s.down_dip == TPV104Params::hypo_down_dip)
      {
         found_hypo = true;
         break;
      }
   }
   TEST_ASSERT(found_hypo,
               "hypocenter station (0, 7.5 km) is in the list");
}

int main(int argc, char *argv[])
{
   TestPsiInitPerQP();
   TestVwSideChannel();
   TestFaultBasisConvention();
   TestPreStressLayout();
   TestNoTotalQSymbols();
   TestInitializeStateZerosQ();
   TestDefaultStations();

   TEST_PRINT_RESULTS();
   return (num_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
