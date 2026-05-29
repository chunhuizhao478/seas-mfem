// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_tpv31_canonical_rotation.cpp — verify that the TPV31 TOML
// (rotated into the canonical SEAS frame; see CLAUDE.md "Canonical
// Coordinate System") still encodes spec-exact physics.
//
// The spec describes a fault on the z=0 plane with y as the depth axis.
// The canonical SEAS frame puts the fault on y=0 with z as the depth
// axis (z < 0 below surface).  The rotation R: (x_s, y_s, z_s) →
// (x_s, z_s, -y_s) leaves the physics invariant; this test asserts the
// rotated TOML values match the spec under R.

#include "mfem.hpp"

#include "../../spatial/code/spatial_friction.hpp"
#include "../../spatial/code/spatial_stress.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

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

static std::string find_tpv31_toml()
{
   const char* candidates[] = {
      "miniapps/seas/tpv31/configs/tpv31.toml",
      "tpv31/configs/tpv31.toml",
      "../tpv31/configs/tpv31.toml",
      "../../tpv31/configs/tpv31.toml",
   };
   for (const char* p : candidates)
   {
      std::ifstream f(p);
      if (f.good()) { return std::string(p); }
   }
   return std::string();
}

// ---------------------------------------------------------------------
// T-1  Canonical fault geometry (matches TPV205/102/104 default).
// ---------------------------------------------------------------------
static void T1_canonical_fault_geometry()
{
   std::cout << "\n[T-1] Canonical FaultGeometry: ref_normal (0,-1,0), up (0,0,1)\n";
   const std::string path = find_tpv31_toml();
   TEST_ASSERT(!path.empty(), "tpv31.toml found on disk");
   if (path.empty()) { return; }
   SpatialFrictionConfig cfg = LoadSpatialFrictionConfig(path);
   TEST_NEAR(cfg.fault_geometry.ref_normal[0],  0.0, 1e-12, "ref_normal[0]");
   TEST_NEAR(cfg.fault_geometry.ref_normal[1], -1.0, 1e-12, "ref_normal[1]");
   TEST_NEAR(cfg.fault_geometry.ref_normal[2],  0.0, 1e-12, "ref_normal[2]");
   TEST_NEAR(cfg.fault_geometry.up[0], 0.0, 1e-12, "up[0]");
   TEST_NEAR(cfg.fault_geometry.up[1], 0.0, 1e-12, "up[1]");
   TEST_NEAR(cfg.fault_geometry.up[2], 1.0, 1e-12, "up[2]");
}

// ---------------------------------------------------------------------
// T-2  Hypocenter (0, 0, -7500) — spec (0, 7500, 0) rotated through R.
// ---------------------------------------------------------------------
static void T2_hypocenter_rotated()
{
   std::cout << "\n[T-2] Hypocenter rotated: spec (0, 7500, 0) → canonical (0, 0, -7500)\n";
   const std::string path = find_tpv31_toml();
   if (path.empty()) { std::cerr << "  SKIP\n"; return; }
   SpatialFrictionConfig cfg = LoadSpatialFrictionConfig(path);
   TEST_NEAR(cfg.hypocenter.x_m,  0.0,    1e-12, "hypo x");
   TEST_NEAR(cfg.hypocenter.y_m,  0.0,    1e-12, "hypo y");
   TEST_NEAR(cfg.hypocenter.z_m, -7500.0, 1e-9,  "hypo z = -7500 (mesh frame)");
}

// ---------------------------------------------------------------------
// T-3  Stress tensor: spec components rotated to canonical.
//      Spec:  σ_xx=60, σ_yy=0,  σ_zz=60, σ_xz=30  (per-µ MPa scaling)
//      Canon: σ_xx=60, σ_yy=60, σ_zz=0,  σ_xy=30
// ---------------------------------------------------------------------
static void T3_stress_tensor_rotated()
{
   std::cout << "\n[T-3] Cauchy stress tensor components rotated correctly\n";
   const std::string path = find_tpv31_toml();
   if (path.empty()) { std::cerr << "  SKIP\n"; return; }
   SpatialFrictionConfig cfg = LoadSpatialFrictionConfig(path);
   TEST_ASSERT(cfg.stress.kind == StressSourceKind::DepthProportionalToShearModulus,
               "stress kind = depth_proportional");
   const auto& d = cfg.stress.depth_proportional;
   // Note: the parser scales TOML MPa values by 1e6 internally.
   TEST_NEAR(d.sigma_xx_per_mu, 60.0e6, 1.0, "σ_xx_per_µ = 60 MPa (unchanged)");
   TEST_NEAR(d.sigma_yy_per_mu, 60.0e6, 1.0, "σ_yy_per_µ = 60 MPa (was σ_zz_spec = 60)");
   TEST_NEAR(d.sigma_zz_per_mu,  0.0,   1.0, "σ_zz_per_µ = 0 (was σ_yy_spec = 0)");
   TEST_NEAR(d.sigma_xy_per_mu, 30.0e6, 1.0, "σ_xy_per_µ = 30 MPa (was σ_xz_spec = 30, right-lateral)");
   TEST_NEAR(d.sigma_xz_per_mu,  0.0,   1.0, "σ_xz_per_µ = 0");
   TEST_NEAR(d.sigma_yz_per_mu,  0.0,   1.0, "σ_yz_per_µ = 0");
   // Trace invariance: σ_xx + σ_yy + σ_zz = 120 MPa in both frames.
   const real_t trace_canon = d.sigma_xx_per_mu + d.sigma_yy_per_mu
                              + d.sigma_zz_per_mu;
   TEST_NEAR(trace_canon, 120.0e6, 1.0,
             "trace(σ) = 120 MPa (rotation-invariant)");
}

// ---------------------------------------------------------------------
// T-4  Material: depth_axis = "z" + layer values match spec.
//      The actual depth-axis evaluation (depth = max(0, -z) for axis 'z')
//      is exercised indirectly via the spatial-driver smoke run; here we
//      assert the layer table and axis field parsed correctly.
// ---------------------------------------------------------------------
static void T4_material_depth_axis_z()
{
   std::cout << "\n[T-4] Material depth_axis = 'z' + layer values\n";
   const std::string path = find_tpv31_toml();
   if (path.empty()) { std::cerr << "  SKIP\n"; return; }
   SpatialFrictionConfig cfg = LoadSpatialFrictionConfig(path);
   TEST_ASSERT(cfg.material.kind == MaterialKind::DepthProfile1D,
               "material kind = depth_profile_1d");
   TEST_ASSERT(cfg.material.depth_axis == 'z',
               "material depth_axis = 'z' (was 'y' pre-rotation)");

   const auto& L = cfg.material.profile_layers;
   TEST_ASSERT(L.size() == 5, "5 layers per spec table");
   if (L.size() < 5) { return; }
   // Spot-check the layer ranges and surface/deepest values (spec p. 4).
   TEST_NEAR(L[0].depth_top_m, 0.0,    1e-9, "Layer 0 depth_top");
   TEST_NEAR(L[0].depth_bot_m, 2400.0, 1e-9, "Layer 0 depth_bot");
   TEST_NEAR(L[0].vp_ms,       4050.0, 1e-6, "Layer 0 vp (surface)");
   TEST_NEAR(L[0].vs_ms,       2250.0, 1e-6, "Layer 0 vs (surface)");
   TEST_NEAR(L[0].rho_kgm3,    2580.0, 1e-6, "Layer 0 rho (surface)");
   TEST_NEAR(L[3].vp_ms,       5750.0, 1e-6, "Layer 3 vp (depth 5–10 km)");
   TEST_NEAR(L[3].vs_ms,       3450.0, 1e-6, "Layer 3 vs (depth 5–10 km)");
   TEST_NEAR(L[4].depth_top_m, 10000.0, 1e-9, "Layer 4 depth_top (≥10 km)");
   TEST_NEAR(L[4].vp_ms,       6500.0,  1e-6, "Layer 4 vp (deep)");
   TEST_NEAR(L[4].rho_kgm3,    3000.0,  1e-6, "Layer 4 rho (deep)");
}

// ---------------------------------------------------------------------
// T-5  Cohesion taper resolved correctly under canonical axis='z'.
//      Spec ramp: C₀(depth) = max(0, 425 · (2400 - depth)) Pa.
//        depth=0    → 1.02 MPa
//        depth=1200 → 510 kPa
//        depth=2400 → 0
//        depth=3000 → 0 (clamped)
// ---------------------------------------------------------------------
static void T5_cohesion_taper_canonical()
{
   std::cout << "\n[T-5] Cohesion taper ramps correctly with canonical axis='z'\n";
   const std::string path = find_tpv31_toml();
   if (path.empty()) { std::cerr << "  SKIP\n"; return; }
   SpatialFrictionConfig cfg = LoadSpatialFrictionConfig(path);
   TEST_ASSERT(cfg.slip_weakening.has_value(),
               "has [friction.slip_weakening]");
   if (!cfg.slip_weakening.has_value()) { return; }

   // Five probe DOFs on the fault (y=0) at increasing depth:
   //   z =  0      → depth 0    → C₀ = 1.02 MPa
   //   z = -1200   → depth 1200 → C₀ = 510 kPa
   //   z = -2400   → depth 2400 → C₀ = 0
   //   z = -3000   → depth 3000 → C₀ = 0 (outside rule's z bounds; default)
   //   z = -16000  → depth 16k  → outside rule + outside rupture area: default cohesion
   Vector coords(3 * 5);
   const real_t zs[5] = {0.0, -1200.0, -2400.0, -3000.0, -16000.0};
   for (int i = 0; i < 5; ++i)
   {
      coords(3*i + 0) = 0.0;
      coords(3*i + 1) = 0.0;
      coords(3*i + 2) = zs[i];
   }
   Array<int> attr(5);
   for (int i = 0; i < 5; ++i) { attr[i] = cfg.boundary.fault_attr; }

   SpatialFrictionResolver r;
   const SlipWeakeningPerDOFParams p =
      r.ResolveSlipWeakening(*cfg.slip_weakening, coords, attr);

   TEST_NEAR(p.cohesion(0), 1.02e6, 1.0,
             "cohesion at z=0:    1.02 MPa (depth 0)");
   TEST_NEAR(p.cohesion(1),  510e3, 1.0,
             "cohesion at z=-1200: 510 kPa (depth 1200)");
   TEST_NEAR(p.cohesion(2),    0.0, 1.0,
             "cohesion at z=-2400: 0 (depth 2400 — taper end)");
   TEST_NEAR(p.cohesion(3),    0.0, 1.0,
             "cohesion at z=-3000: 0 (outside taper rule bounds)");
   TEST_NEAR(p.cohesion(4),    0.0, 1.0,
             "cohesion at z=-16000: 0 (deep)");
}

// ---------------------------------------------------------------------
// T-6  Nucleation centre rotated.
// ---------------------------------------------------------------------
static void T6_nucleation_center_rotated()
{
   std::cout << "\n[T-6] Nucleation centre rotated: spec (0, 7500, 0) → canonical (0, 0, -7500)\n";
   const std::string path = find_tpv31_toml();
   if (path.empty()) { std::cerr << "  SKIP\n"; return; }
   SpatialFrictionConfig cfg = LoadSpatialFrictionConfig(path);
   TEST_ASSERT(cfg.nucleation.enabled,
               "nucleation enabled");
   TEST_ASSERT(cfg.nucleation.kind ==
               NucleationKind::InstantaneousOverstressCircular,
               "kind = instantaneous_overstress_circular");
   // Spatial-schema InstantaneousOverstressCircularSpec (radius_m + taper_m
   // reproduce the spec's r=1400→2000 cosine taper; delta_tau_pa is the
   // strike-direction overstress — safs's nucleation is strike-only).
   const auto& n = cfg.nucleation.instantaneous_circular;
   TEST_NEAR(n.center_x_m,  0.0,    1e-12, "nuc center_x");
   TEST_NEAR(n.center_y_m,  0.0,    1e-12, "nuc center_y (on fault)");
   TEST_NEAR(n.center_z_m, -7500.0, 1e-9,  "nuc center_z = -7500");
   TEST_NEAR(n.radius_m, 1400.0, 1e-6, "radius_m (full-amplitude, spec inner)");
   TEST_NEAR(n.taper_m,   600.0, 1e-6, "taper_m (cosine width, spec outer-inner)");
   TEST_NEAR(n.delta_tau_pa, 4.95e6, 1.0, "Δτ = 4.95 MPa (strike)");
}

// ---------------------------------------------------------------------
// T-7  Numerics opt-ins inherited from the TPV205 pattern.
// ---------------------------------------------------------------------
static void T7_numerics_opt_ins()
{
   std::cout << "\n[T-7] Numerics: TPV205-style opt-ins\n";
   const std::string path = find_tpv31_toml();
   if (path.empty()) { std::cerr << "  SKIP\n"; return; }
   SpatialFrictionConfig cfg = LoadSpatialFrictionConfig(path);
   // Spatial schema: cfl_safety / fault_iterator / interior_flux are parsed
   // into enums (not strings).
   TEST_ASSERT(cfg.numerics.cfl_safety == CflSafety::Dg,
               "cfl_safety = dg (matches TPV* native)");
   // The spatial driver always sub-steps ("one-shot" is rejected); the TOML
   // was adapted from the spec's one-shot to the supported "substep".
   TEST_ASSERT(cfg.numerics.fault_iterator == FaultIteratorKind::Substep,
               "fault_iterator = substep (spatial driver always sub-steps)");
   // TPV31's depth_profile_1d material is Mode::Coefficient → requires the
   // matrix (bimaterial) interior-flux path.
   TEST_ASSERT(cfg.numerics.interior_flux == InteriorFlux::Matrix,
               "interior_flux = matrix (TPV31 heterogeneous material)");
}

int main(int, char**)
{
   std::cout << "Running test_tpv31_canonical_rotation "
                "(TPV31 spec → canonical SEAS frame)\n";
   T1_canonical_fault_geometry();
   T2_hypocenter_rotated();
   T3_stress_tensor_rotated();
   T4_material_depth_axis_z();
   T5_cohesion_taper_canonical();
   T6_nucleation_center_rotated();
   T7_numerics_opt_ins();
   std::cout << "\n========================================\n";
   std::cout << "test_tpv31_canonical_rotation: "
             << num_passed << " / " << num_tests
             << " passed, " << num_failed << " failed\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
