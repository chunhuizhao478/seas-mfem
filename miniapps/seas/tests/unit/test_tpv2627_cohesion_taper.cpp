// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_tpv2627_cohesion_taper.cpp — Phase 3 of
// PLAN_TPV26_27_spatial_dyn_driver_2026-06-28.md §5.
//
// TPV26/27 depth-dependent frictional cohesion (spec Part 4; PLAN §1.3):
//
//   C0(depth) = 0.40 MPa + 0.00072 MPa/m * (5000 - depth)   for depth <= 5000
//             = 0.40 MPa                                     for depth >  5000
//
// i.e. 4.0 MPa at the free surface, tapering to 0.40 MPa at/below 5 km (the
// spec adds this cohesion to "suppress free-surface effects").
//
// This is delivered CONFIG-ONLY: SpatialFrictionResolver::ResolveSlipWeakening
// already implements  C0 = max(floor, floor + grad * (ref_depth - depth))
// with depth = max(0, -z)  (spatial_friction.cpp:1949-1971).
//
// [CORR vs PLAN §5 Phase 3]: the plan says to put the four taper keys in the
// TOP-LEVEL [friction.slip_weakening] block.  They actually live on a
// SpatialRule and are parsed from a [[friction.slip_weakening.spatial]] rule
// (spatial_friction.hpp:418-431, spatial_friction.cpp:540-559), and the taper
// only fires for DOFs that a rule MATCHES.  Both tests below therefore use a
// whole-fault `kind="depth"` rule (default +/-inf bounds match every DOF).
//
// Standalone: no MPI, no mesh — ResolveSlipWeakening takes coords + attrs only.

#include "mfem.hpp"

#include "../../spatial/code/spatial_friction.hpp"

#include "test_macros.hpp"

#include <iostream>
#include <string>

using namespace mfem;
using namespace mfem::seas::spatial;

namespace
{

// Four DOFs at depths 0 / 2500 / 5000 / 10000 m  (code frame: depth = -z).
constexpr int kN = 4;
const real_t kDepths[kN] = { 0.0, 2500.0, 5000.0, 10000.0 };
// Spec expectations (Pa): 4.0, 2.2, 0.40, 0.40 MPa.
const real_t kExpectC0[kN] = { 4.0e6, 2.2e6, 0.40e6, 0.40e6 };

Vector MakeCoords()
{
   Vector coords(3 * kN);
   for (int i = 0; i < kN; ++i)
   {
      coords(3 * i + 0) = 0.0;
      coords(3 * i + 1) = 0.0;
      coords(3 * i + 2) = -kDepths[i];   // z < 0 is down
   }
   return coords;
}

Array<int> MakeAttrs()
{
   Array<int> attr(kN);
   attr = 101;
   return attr;
}

void CheckCohesion(const SlipWeakeningPerDOFParams& p, const char* what)
{
   for (int i = 0; i < kN; ++i)
   {
      TEST_NEAR(p.cohesion(i), kExpectC0[i], 1e3,
                std::string(what) + ": C0(depth="
                + std::to_string(int(kDepths[i])) + " m) = "
                + std::to_string(kExpectC0[i] / 1e6) + " MPa");
   }
}

}  // namespace

// The taper math, driven by a hand-built SpatialRule (no TOML).
void Test_Cohesion_Taper_Direct()
{
   SlipWeakeningBlock blk;
   blk.mu_s_default     = 0.18;
   blk.mu_d_default     = 0.12;
   blk.d_c_default      = 0.30;
   blk.cohesion_default = 0.0;

   SpatialRule r;                       // whole-fault rule (+/-inf bounds)
   r.kind                   = SpatialRule::Kind::Depth;
   r.cohesion_floor_pa      = 0.40e6;
   r.cohesion_grad_pa_per_m = 720.0;    // 0.00072 MPa/m
   r.cohesion_ref_depth_m   = 5000.0;
   r.cohesion_taper_axis    = 'z';
   blk.spatial.push_back(r);

   SpatialFrictionResolver R;
   const Vector coords = MakeCoords();
   const Array<int> attr = MakeAttrs();
   const SlipWeakeningPerDOFParams p = R.ResolveSlipWeakening(blk, coords, attr);

   TEST_ASSERT(p.cohesion.Size() == kN, "per-DOF cohesion sized to fault DOFs");
   CheckCohesion(p, "direct");

   // Below the reference depth the raw taper goes negative; the floor clamps.
   TEST_ASSERT(p.cohesion(3) >= 0.40e6 - 1.0,
               "deep DOF clamped from below by cohesion_floor_pa (not negative)");
   // mu_s / mu_d / d_c pass through untouched by the taper.
   TEST_NEAR(p.mu_s(0), 0.18, 1e-15, "mu_s unchanged by cohesion taper");
   TEST_NEAR(p.mu_d(0), 0.12, 1e-15, "mu_d unchanged by cohesion taper");
}

// The exact TOML surface a TPV26/27 config uses (validates the parser path,
// including that the four taper keys belong to a [[...spatial]] rule).
void Test_Cohesion_Taper_From_Toml()
{
   // parse_root requires the nine top-level blocks (spatial_friction.cpp:839-876);
   // the ones below the friction block are minimal placeholders.  [stress] uses
   // the Phase-1 tpv2627_depth kind, so this also covers that parser branch.
   const std::string toml = R"TOML(
[meta]
schema_version = 1
law            = "slip_weakening"

[material_constant_fallback]
lambda = 32043759360.0
mu     = 32038120320.0
rho    = 2670.0

[pore_pressure]
P_p_pa            = 0.0
P_p_grad_pa_per_m = 9800.0
min_sigma_n_pa    = 0.0

[mesh]
path  = "unused-by-this-test.msh"
order = 1

[velocity]
use_sidecar = false

[stress]
kind = "tpv2627_depth"

[numerics]
ader_order = 2

[time]
tfinal = "1.0s"

[output]
output_dir = "unused-by-this-test"

[friction.slip_weakening]
mu_s_default     = 0.18
mu_d_default     = 0.12
d_c_default      = 0.30
cohesion_default = 0.0

[[friction.slip_weakening.spatial]]
kind                   = "depth"
cohesion_floor_pa      = 0.40e6
cohesion_grad_pa_per_m = 720.0
cohesion_ref_depth_m   = 5000.0
cohesion_taper_axis    = "z"
)TOML";

   const SpatialFrictionConfig cfg = ParseSpatialFrictionConfigString(toml);
   TEST_ASSERT(cfg.slip_weakening.has_value(),
               "[friction.slip_weakening] parsed");

   SpatialFrictionResolver R;
   const Vector coords = MakeCoords();
   const Array<int> attr = MakeAttrs();
   const SlipWeakeningPerDOFParams p =
      R.ResolveSlipWeakening(cfg.slip_weakening.value(), coords, attr);
   CheckCohesion(p, "toml");
}

// Without a matching spatial rule the taper must NOT fire — every DOF keeps
// cohesion_default.  (Guards the [CORR]: putting the keys in the top-level
// block would silently produce a constant cohesion.)
void Test_No_Rule_Means_No_Taper()
{
   SlipWeakeningBlock blk;
   blk.mu_s_default     = 0.18;
   blk.mu_d_default     = 0.12;
   blk.d_c_default      = 0.30;
   blk.cohesion_default = 0.40e6;       // no spatial rule at all

   SpatialFrictionResolver R;
   const Vector coords = MakeCoords();
   const Array<int> attr = MakeAttrs();
   const SlipWeakeningPerDOFParams p = R.ResolveSlipWeakening(blk, coords, attr);

   for (int i = 0; i < kN; ++i)
   {
      TEST_NEAR(p.cohesion(i), 0.40e6, 1e-6,
                "no spatial rule -> constant cohesion_default (no taper)");
   }
}

int main(int /*argc*/, char * /*argv*/[])
{
   std::cout << "=== test_tpv2627_cohesion_taper (Phase 3) ===\n";
   Test_Cohesion_Taper_Direct();
   Test_Cohesion_Taper_From_Toml();
   Test_No_Rule_Means_No_Taper();
   TEST_PRINT_RESULTS();
   return (num_failed == 0) ? 0 : 1;
}
