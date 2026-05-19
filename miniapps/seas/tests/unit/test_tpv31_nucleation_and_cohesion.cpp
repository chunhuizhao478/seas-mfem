// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_tpv31_nucleation_and_cohesion.cpp — runtime-side coverage for
// the two new pieces of infrastructure added for spec-exact TPV31:
//
//   1. `InstantaneousOverstressCircularSpec` + resolver +
//      one-shot applicator (circular cosine-tapered overstress with
//      per-DOF µ-scaling, applied once at init).
//   2. `SpatialRule` depth-linear cohesion taper
//      (`cohesion_grad_pa_per_m` + `cohesion_ref_depth_m` +
//      `cohesion_floor_pa`).
//
// These are the gaps that previously forced tpv31.toml to document
// "OUT OF SCOPE" + "operational stopgap" approximations.

#include "mfem.hpp"

#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/spatial_nucleation.hpp"
#include "../../spatial/code/spatial_friction.hpp"

#include <cmath>
#include <iostream>

using namespace mfem;
using namespace mfem::seas;
using namespace mfem::seas::spatial;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(c, m) do { num_tests++; if (!(c)) { \
   std::cerr << "FAILED: " << m << " line " << __LINE__ << "\n"; \
   num_failed++; } else { std::cout << "  PASSED: " << m << "\n"; \
   num_passed++; } } while (0)

#define TEST_NEAR(v, e, t, m) do { num_tests++; const real_t _v=(v); \
   const real_t _e=(e); const real_t _t=(t); \
   if (std::abs(_v - _e) > _t) { std::cerr << "FAILED: " << m \
   << " (got " << _v << ", expected " << _e \
   << ", tol " << _t << ")\n"; num_failed++; } \
   else { std::cout << "  PASSED: " << m << "\n"; num_passed++; } \
   } while (0)

// =====================================================================
//  Group I — InstantaneousOverstressCircular
// =====================================================================

static InstantaneousOverstressCircularSpec tpv31_nuc_spec()
{
   InstantaneousOverstressCircularSpec s;
   s.center_x_m         = 0.0;
   s.center_y_m         = 7500.0;
   s.center_z_m         = 0.0;
   s.radius_inner_m     = 1400.0;
   s.radius_outer_m     = 2000.0;
   s.delta_tau_peak_pa  = 4.95e6;
   s.mu_ref_pa          = 32.03812032e9;
   s.dip_fraction       = 0.0;
   s.strike_fraction    = 1.0;
   return s;
}

// I-1  At the hypocenter (r = 0): amplitude = peak · (µ_local / µ_ref).
static void I_1_hypocenter_full_amplitude()
{
   std::cout << "\n[I-1] Hypocenter sees full amplitude · µ/µ_ref\n";
   const auto spec = tpv31_nuc_spec();
   Vector coords(3);
   coords(0) = 0.0; coords(1) = 7500.0; coords(2) = 0.0;
   Vector mu(1);
   mu(0) = spec.mu_ref_pa;   // µ = µ_ref → scale = 1
   const auto p = ResolveInstantaneousOverstressCircular(spec,
                                                          /*enabled=*/true,
                                                          coords, mu);
   TEST_NEAR(p.amplitude_strike(0), spec.delta_tau_peak_pa, 0.0,
             "centre amplitude_strike = delta_tau_peak");
   TEST_NEAR(p.amplitude_dip(0),    0.0,                    0.0,
             "centre amplitude_dip = 0 (pure strike)");
   TEST_NEAR(p.radial(0),           1.0,                    0.0,
             "centre radial = 1");
}

// I-2  At r = r_inner: still full amplitude (the formula is f(r_inner) = 1).
static void I_2_inner_edge_full_amplitude()
{
   std::cout << "\n[I-2] Inner edge (r = r_inner) is still full amplitude\n";
   const auto spec = tpv31_nuc_spec();
   Vector coords(3);
   coords(0) = 1400.0; coords(1) = 7500.0; coords(2) = 0.0;   // r = 1400
   Vector mu(1);  mu(0) = spec.mu_ref_pa;
   const auto p = ResolveInstantaneousOverstressCircular(spec, true,
                                                          coords, mu);
   TEST_NEAR(p.amplitude_strike(0), spec.delta_tau_peak_pa, 1e-9,
             "r = r_inner amplitude_strike = delta_tau_peak");
}

// I-3  At r = midpoint of taper: amplitude = 0.5 · peak · (cos(π/2)+1)/2
//      = 0.5 · peak · 0.5 = 0.25 · peak.  Wait, cos(π/2) = 0 → 0.5(1+0) = 0.5.
//      So at r = (r_inner + r_outer)/2 = 1700: spatial factor = 0.5.
static void I_3_taper_midpoint_half_amplitude()
{
   std::cout << "\n[I-3] Midpoint of taper has spatial factor = 0.5\n";
   const auto spec = tpv31_nuc_spec();
   Vector coords(3);
   coords(0) = 1700.0; coords(1) = 7500.0; coords(2) = 0.0;
   Vector mu(1);  mu(0) = spec.mu_ref_pa;
   const auto p = ResolveInstantaneousOverstressCircular(spec, true,
                                                          coords, mu);
   TEST_NEAR(p.amplitude_strike(0), 0.5 * spec.delta_tau_peak_pa, 1e-3,
             "taper-midpoint amplitude_strike = 0.5 · peak");
   TEST_NEAR(p.radial(0), 0.5, 1e-9,
             "taper-midpoint radial = 0.5");
}

// I-4  At r = r_outer: amplitude = 0 (spatial factor = 0.5 · (1 + cos π) = 0).
static void I_4_outer_edge_zero()
{
   std::cout << "\n[I-4] Outer edge (r = r_outer) has amplitude = 0\n";
   const auto spec = tpv31_nuc_spec();
   Vector coords(3);
   coords(0) = 2000.0; coords(1) = 7500.0; coords(2) = 0.0;   // r = 2000
   Vector mu(1);  mu(0) = spec.mu_ref_pa;
   const auto p = ResolveInstantaneousOverstressCircular(spec, true,
                                                          coords, mu);
   TEST_NEAR(p.amplitude_strike(0), 0.0, 1e-9,
             "r = r_outer amplitude_strike = 0");
   TEST_NEAR(p.radial(0), 0.0, 1e-9,
             "r = r_outer radial = 0");
}

// I-5  Beyond r_outer: amplitude = 0.
static void I_5_beyond_outer_zero()
{
   std::cout << "\n[I-5] Beyond r_outer (r = 5000): amplitude = 0\n";
   const auto spec = tpv31_nuc_spec();
   Vector coords(3);
   coords(0) = 5000.0; coords(1) = 7500.0; coords(2) = 0.0;
   Vector mu(1);  mu(0) = spec.mu_ref_pa;
   const auto p = ResolveInstantaneousOverstressCircular(spec, true,
                                                          coords, mu);
   TEST_NEAR(p.amplitude_strike(0), 0.0, 0.0, "far-field amplitude = 0");
}

// I-6  µ-scaling: with µ_local = 2 · µ_ref the amplitude should double.
static void I_6_mu_scaling()
{
   std::cout << "\n[I-6] Per-DOF µ-scaling: amplitude = peak · µ/µ_ref\n";
   const auto spec = tpv31_nuc_spec();
   Vector coords(3);
   coords(0) = 0.0; coords(1) = 7500.0; coords(2) = 0.0;
   Vector mu(1);  mu(0) = 2.0 * spec.mu_ref_pa;
   const auto p = ResolveInstantaneousOverstressCircular(spec, true,
                                                          coords, mu);
   TEST_NEAR(p.amplitude_strike(0), 2.0 * spec.delta_tau_peak_pa, 1e-3,
             "amplitude doubles with µ doubling");
}

// I-7  ApplyInstantaneousOverstressCircular writes the FULL amplitude
//      into DOFData::tau{1,2}_nuc in a single call (no smoothStep ramp).
static void I_7_one_shot_applicator()
{
   std::cout << "\n[I-7] One-shot applicator writes full amplitude immediately\n";
   InstantaneousOverstressCircularPerDOFParams p;
   p.amplitude_dip.SetSize(1);    p.amplitude_dip    = 7.0;
   p.amplitude_strike.SetSize(1); p.amplitude_strike = 11.0;
   p.radial.SetSize(1);           p.radial           = 1.0;
   std::vector<DOFData> dof(1);
   dof[0].tau1_nuc = 0.0;  dof[0].tau2_nuc = 0.0;
   ApplyInstantaneousOverstressCircular(dof, p);
   TEST_NEAR(dof[0].tau1_nuc,  7.0, 0.0, "tau1_nuc = amplitude_dip (one-shot)");
   TEST_NEAR(dof[0].tau2_nuc, 11.0, 0.0, "tau2_nuc = amplitude_strike (one-shot)");
   // A second call would double the amplitudes — the contract is that
   // the driver calls this EXACTLY ONCE during init.  Verify the math:
   ApplyInstantaneousOverstressCircular(dof, p);
   TEST_NEAR(dof[0].tau1_nuc, 14.0, 0.0, "second call doubles (additive — driver must not repeat)");
}

// I-8  Disabled mode returns zero-sized Vectors → applicator is no-op.
static void I_8_disabled_no_op()
{
   std::cout << "\n[I-8] Disabled mode produces no-op applicator\n";
   const auto spec = tpv31_nuc_spec();
   Vector coords(0), mu(0);
   const auto p = ResolveInstantaneousOverstressCircular(spec,
                                                          /*enabled=*/false,
                                                          coords, mu);
   TEST_ASSERT(p.amplitude_dip.Size() == 0, "disabled → empty Vectors");
   std::vector<DOFData> dof(1);
   dof[0].tau1_nuc = 42.0;  dof[0].tau2_nuc = 99.0;
   ApplyInstantaneousOverstressCircular(dof, p);  // must not crash
   TEST_NEAR(dof[0].tau1_nuc, 42.0, 0.0, "tau1_nuc unchanged");
   TEST_NEAR(dof[0].tau2_nuc, 99.0, 0.0, "tau2_nuc unchanged");
}

// =====================================================================
//  Group C — SpatialRule depth-linear cohesion taper
// =====================================================================

// C-1  TPV31 spec ramp: C_0(y) = max(0, 425 · (2400 - y)).  At y = 0:
//      1.02 MPa.  At y = 1200: 0.51 MPa.  At y = 2400: 0.  At y > 2400: 0.
static void C_1_tpv31_cohesion_ramp()
{
   std::cout << "\n[C-1] TPV31 cohesion ramp matches spec at sample depths\n";
   SlipWeakeningBlock cfg;
   cfg.mu_s_default     = 0.580;
   cfg.mu_d_default     = 0.450;
   cfg.d_c_default      = 0.18;
   cfg.cohesion_default = 0.0;

   SpatialRule r;
   r.kind = SpatialRule::Kind::Depth;
   r.y_min_m = 0.0;
   r.y_max_m = 2400.0;
   r.cohesion_taper_axis    = 'y';
   r.cohesion_ref_depth_m   = 2400.0;
   r.cohesion_grad_pa_per_m = 425.0;
   r.cohesion_floor_pa      = 0.0;
   cfg.spatial.push_back(r);

   // 4 DOFs at y = 0, 1200, 2400, 5000 (all on the fault plane z = 0,
   // x = 0); the rule.matches() filter uses only y for Depth-kind, so
   // x/z values don't matter as long as we pass them.
   Vector coords(12);
   for (int i = 0; i < 4; ++i)
   {
      coords(3*i + 0) = 0.0;
      coords(3*i + 2) = 0.0;
   }
   coords(1)  =    0.0;
   coords(4)  = 1200.0;
   coords(7)  = 2400.0;
   coords(10) = 5000.0;
   Array<int> attr(4); attr = 1;

   SpatialFrictionResolver resolver;
   const auto p = resolver.ResolveSlipWeakening(cfg, coords, attr);

   TEST_NEAR(p.cohesion(0), 1.02e6, 1e-3,  "y=0    → 1.02 MPa");
   TEST_NEAR(p.cohesion(1), 0.51e6, 1e-3,  "y=1200 → 0.51 MPa");
   TEST_NEAR(p.cohesion(2), 0.0,    0.0,   "y=2400 → 0");
   // At y = 5000, the rule.matches() Depth-kind y_max_m = 2400 filter
   // rejects this DOF (y > y_max), so cohesion stays at the default 0.
   TEST_NEAR(p.cohesion(3), 0.0,    0.0,   "y=5000 → default 0 (outside rule)");
}

// C-2  Floor clamp: with floor = 200 kPa, y > 2400 should be clamped to
//      200 kPa, not run negative.  We have to enlarge y_max_m so the
//      rule matches y > 2400, otherwise the rule.matches() filter
//      short-circuits.
static void C_2_floor_clamp()
{
   std::cout << "\n[C-2] Floor clamp prevents negative cohesion past ref_depth\n";
   SlipWeakeningBlock cfg;
   cfg.mu_s_default     = 0.580;
   cfg.mu_d_default     = 0.450;
   cfg.d_c_default      = 0.18;
   cfg.cohesion_default = 0.0;

   SpatialRule r;
   r.kind = SpatialRule::Kind::Depth;
   // Loose bounds so the rule matches at every test DOF.
   r.y_min_m = -1.0e9;
   r.y_max_m =  1.0e9;
   r.cohesion_taper_axis    = 'y';
   r.cohesion_ref_depth_m   = 2400.0;
   r.cohesion_grad_pa_per_m = 425.0;
   r.cohesion_floor_pa      = 200.0e3;
   cfg.spatial.push_back(r);

   Vector coords(6);
   coords(0) = 0.0; coords(1) =    0.0; coords(2) = 0.0;
   coords(3) = 0.0; coords(4) = 5000.0; coords(5) = 0.0;
   Array<int> attr(2); attr = 1;

   SpatialFrictionResolver resolver;
   const auto p = resolver.ResolveSlipWeakening(cfg, coords, attr);

   // At y=0: raw = 200e3 + 425·2400 = 200e3 + 1.02e6 = 1.22e6.  Above floor.
   TEST_NEAR(p.cohesion(0), 1.22e6, 1e-3, "y=0    → 1.22 MPa (above floor)");
   // At y=5000: raw = 200e3 + 425·(2400-5000) = 200e3 - 1.105e6 = -905e3.
   // Clamped to floor = 200e3.
   TEST_NEAR(p.cohesion(1), 200.0e3, 1e-3, "y=5000 → clamped to floor 200 kPa");
}

// C-3  Axis 'x' / 'z' selection: the taper axis is independent of the
//      rule's coordinate bounds (which use y by default for Depth kind).
//      Verify the axis-x taper uses x_dof for the ramp.
static void C_3_axis_selection()
{
   std::cout << "\n[C-3] cohesion_taper_axis selects which coord to ramp on\n";
   SlipWeakeningBlock cfg;
   cfg.mu_s_default     = 0.580;
   cfg.mu_d_default     = 0.450;
   cfg.d_c_default      = 0.18;
   cfg.cohesion_default = 0.0;

   SpatialRule r;
   r.kind = SpatialRule::Kind::Box;        // match by all 3 axes
   r.x_min_m = -1e9;  r.x_max_m = 1e9;
   r.y_min_m = -1e9;  r.y_max_m = 1e9;
   r.z_min_m = -1e9;  r.z_max_m = 1e9;
   r.cohesion_taper_axis    = 'x';
   r.cohesion_ref_depth_m   = 2400.0;
   r.cohesion_grad_pa_per_m = 425.0;
   r.cohesion_floor_pa      = 0.0;
   cfg.spatial.push_back(r);

   // Hold y fixed, vary x.
   Vector coords(6);
   coords(0) =    0.0; coords(1) = 10000.0; coords(2) = 0.0;
   coords(3) = 1200.0; coords(4) = 10000.0; coords(5) = 0.0;
   Array<int> attr(2); attr = 1;

   SpatialFrictionResolver resolver;
   const auto p = resolver.ResolveSlipWeakening(cfg, coords, attr);
   TEST_NEAR(p.cohesion(0), 1.02e6, 1e-3, "x=0    → 1.02 MPa (axis-x taper)");
   TEST_NEAR(p.cohesion(1), 0.51e6, 1e-3, "x=1200 → 0.51 MPa (axis-x taper)");
}

int main(int, char**)
{
   std::cout << "Running test_tpv31_nucleation_and_cohesion\n";
   I_1_hypocenter_full_amplitude();
   I_2_inner_edge_full_amplitude();
   I_3_taper_midpoint_half_amplitude();
   I_4_outer_edge_zero();
   I_5_beyond_outer_zero();
   I_6_mu_scaling();
   I_7_one_shot_applicator();
   I_8_disabled_no_op();
   C_1_tpv31_cohesion_ramp();
   C_2_floor_clamp();
   C_3_axis_selection();
   std::cout << "\n========================================\n";
   std::cout << "test_tpv31_nucleation_and_cohesion: "
             << num_passed << " / " << num_tests
             << " passed, " << num_failed << " failed\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
