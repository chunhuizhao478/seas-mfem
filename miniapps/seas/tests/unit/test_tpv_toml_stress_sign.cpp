// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_tpv_toml_stress_sign.cpp — Phase 8 req 2 / D3.2
// (PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24.md): build each
// TPV FaultLocalPrestress source FROM ITS TOML, seed a canonical-frame DOF
// via FaultGeometry::ComputeParamsFaultLocal, and assert the right-lateral /
// compression POSITIVE sign chain:
//   tau2_0 (strike) == +tau_strike_pa   (== +tau_ini_native, no projection)
//   tau1_0 (dip)    == tau_dip_pa        (== 0, pure strike-slip)
//   sigma_n0        == sigma_n_pa - P_p  (effective normal stress)
// Precondition (D1): the canonical fault-local frame is dip=(0,0,-1),
// strike=(+1,0,0) for ref_normal=(0,-1,0), up=(0,0,1) — checked here via
// the pure FaultBasis::ComputeOrientedFrame kernel (the full mesh/per-QP
// version is test_planar_tpv_basis.cpp).
//
// Native references (config/tpv*_params.hpp): TPV205 tau_ini = 70 MPa
// (background), TPV102 75 MPa, TPV104 40 MPa; sigma_n = 120 MPa all three.
//
// The seeding fixture reuses the BP5 mesh (the seed is frame-agnostic and
// uniform, so it validates the TOML->tau_pre_ flow on any fault mesh); if no
// mesh is present the seeding sub-checks are skipped (compile-only), matching
// test_compute_safs_params.cpp.

#include "mfem.hpp"

#include "../../spatial/code/spatial_friction.hpp"
#include "../../fault/fault_geometry.hpp"
#include "../../fault/fault_geometry_safs.inl"
#include "../../fault/fault_basis.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../domain/boundary_config.hpp"
#include "../../config/bp5_params.hpp"

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <sys/stat.h>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(c, m) do { num_tests++; if (!(c)) { \
   std::cerr << "FAILED: " << m << " line " << __LINE__ << "\n"; \
   num_failed++; } else { std::cout << "  PASSED: " << m << "\n"; \
   num_passed++; } } while (0)

#define TEST_NEAR(v, e, t, m) do { num_tests++; const real_t _v=(v); \
   const real_t _e=(e); const real_t _t=(t); \
   if (std::abs(_v - _e) > _t) { std::cerr << "FAILED: " << m \
   << " (got " << _v << ", expected " << _e << ", tol " << _t \
   << ", line " << __LINE__ << ")\n"; num_failed++; } else { \
   std::cout << "  PASSED: " << m << "\n"; num_passed++; } } while (0)

namespace
{
bool FileExists(const std::string& p) { struct stat s; return stat(p.c_str(), &s) == 0; }

std::string FindConfig(const std::string& rel)
{
   for (const auto& pre : {std::string(""), std::string("../"), std::string("../../")})
   { const std::string p = pre + rel; if (FileExists(p)) { return p; } }
   return "";
}
std::string FindBP5Mesh()
{
   const std::vector<std::string> c = {
      "bp5/mesh/bp5_1000m.msh", "../bp5/mesh/bp5_1000m.msh",
      "/Users/chunhuizhao/projects/seas-mfem-safs/miniapps/seas/bp5/mesh/bp5_1000m.msh",
   };
   for (const auto& p : c) { if (FileExists(p)) { return p; } }
   return "";
}

struct Fixture
{
   BP5Params params;
   std::unique_ptr<Mesh> mesh;
   std::unique_ptr<ElasticityDomainOperator<Mesh>> domain_op;
   int nf = 0;
   bool loaded = false;
};

std::unique_ptr<Fixture> BuildFixture()
{
   auto fix = std::make_unique<Fixture>();
   const std::string mp = FindBP5Mesh();
   if (mp.empty()) { return fix; }
   fix->mesh = std::make_unique<Mesh>(mp.c_str(), 1, 1);
   fix->loaded = true;
   BoundaryConfig bc;
   bc.fault_attr = 100;
   bc.dirichlet_attrs = {1, 2, 3, 4};
   bc.natural_attrs   = {5, 6};
   bc.default_dirichlet_func = MakeBP5DirichletFunc(fix->params.Vp);
   LinearElastic le(fix->params.lambda(), fix->params.mu());
   DomainConfig dc;
   fix->domain_op = std::make_unique<ElasticityDomainOperator<Mesh>>(
      *fix->mesh, 1, le, fix->params.Vp, fix->params.Wf, fix->params.lf,
      bc, DGMethod::IP, SolverType::MUMPS_BLR, dc);
   fix->nf = fix->domain_op->GetNumFaultDOFs();
   return fix;
}

// Drive the seed exactly as the driver does for kind="fault_local_prestress",
// from PARSED config values, and assert the sign chain on every fault DOF.
void CheckSeedFromConfig(Fixture& fix, const spatial::SpatialFrictionConfig& cfg,
                         const char* tag)
{
   const real_t tau_strike = cfg.stress.tau_strike_pa;
   const real_t tau_dip    = cfg.stress.tau_dip_pa;
   const real_t sigma_n    = cfg.stress.sigma_n_pa;
   const real_t P_p        = cfg.stress.pore_pressure.P_p_pa;
   const real_t sigma_n_eff = sigma_n - P_p;

   FaultGeometry<Mesh> fg(*fix.domain_op, fix.params);
   fg.ComputeParamsFaultLocal(tau_strike, tau_dip, sigma_n, P_p);

   TEST_ASSERT(fg.HasParams(), std::string(tag) + ": HasParams() after seed");
   TEST_ASSERT(fg.GetTauPre().Size() == 2 * fix.nf,
               std::string(tag) + ": tau_pre size == 2N");
   const Vector& sn  = fg.sigma_n_per_dof();
   const Vector& tau = fg.GetTauPre();
   real_t d_sn = 0, d_dip = 0, d_strike = 0;
   for (int i = 0; i < fix.nf; ++i)
   {
      d_sn     = std::max(d_sn,     std::abs(sn(i)          - sigma_n_eff));
      d_dip    = std::max(d_dip,    std::abs(tau(2 * i)     - tau_dip));
      d_strike = std::max(d_strike, std::abs(tau(2 * i + 1) - tau_strike));
   }
   TEST_NEAR(d_strike, 0.0, 1e-6,
             std::string(tag) + ": tau2_0(strike) == +tau_strike_pa (right-lateral +)");
   TEST_NEAR(d_dip, 0.0, 1e-6,
             std::string(tag) + ": tau1_0(dip) == tau_dip_pa (0)");
   TEST_NEAR(d_sn, 0.0, 1e-6,
             std::string(tag) + ": sigma_n0 == sigma_n_pa - P_p");
   TEST_ASSERT(fix.nf == 0 || tau(1) > 0.0,
               std::string(tag) + ": positive tau_strike_pa -> positive strike slot");
}
} // anon

// Canonical-frame precondition (D1): dip=(0,0,-1), strike=(+1,0,0).
static void T_canonical_frame_precondition()
{
   std::cout << "\n[precondition] canonical fault-local frame (D1)\n";
   Vector ref_normal(3); ref_normal = 0.0; ref_normal(1) = -1.0;
   Vector up(3);         up = 0.0;         up(2) = 1.0;
   Vector n_raw(3);      n_raw = 0.0;      n_raw(1) = -1.0;
   real_t normal[3], dip[3], strike[3], nl; bool flipped;
   FaultBasis::ComputeOrientedFrame(n_raw, 3, ref_normal, up,
                                    normal, dip, strike, flipped, nl);
   TEST_NEAR(dip[2],   -1.0, 1e-12, "precondition dip == (0,0,-1)");
   TEST_NEAR(strike[0], 1.0, 1e-12, "precondition strike == (+1,0,0)");
}

int main(int, char**)
{
   std::cout << "Running Phase 8 req 2 TPV TOML stress-sign tests\n";
   T_canonical_frame_precondition();

   struct Case { const char* tag; std::string rel; real_t tau_ini; };
   std::vector<Case> cases = {
      {"TPV205", "tpv205/configs/tpv205_spatial.toml", 70.0e6},
      {"TPV102", "tpv102/configs/tpv102_spatial.toml", 75.0e6},
      {"TPV104", "tpv104/configs/tpv104_spatial.toml", 40.0e6},
   };

   // Config-level sign chain (always available: parse only).
   std::vector<spatial::SpatialFrictionConfig> cfgs;
   bool all_found = true;
   for (auto& c : cases)
   {
      const std::string path = FindConfig(c.rel);
      if (path.empty()) { all_found = false; break; }
      std::cout << "\n[" << c.tag << "] " << path << "\n";
      spatial::SpatialFrictionConfig cfg = spatial::LoadSpatialFrictionConfig(path);
      TEST_ASSERT(cfg.stress.kind == spatial::StressSourceKind::FaultLocalPrestress,
                  std::string(c.tag) + ": kind == fault_local_prestress");
      TEST_NEAR(cfg.stress.tau_strike_pa, c.tau_ini, 1.0,
                std::string(c.tag) + ": tau_strike_pa == +tau_ini (right-lateral +)");
      TEST_NEAR(cfg.stress.tau_dip_pa, 0.0, 1.0,
                std::string(c.tag) + ": tau_dip_pa == 0");
      TEST_NEAR(cfg.stress.sigma_n_pa, 120.0e6, 1.0,
                std::string(c.tag) + ": sigma_n_pa == 120 MPa (compression +)");
      cfgs.push_back(cfg);
   }

   if (!all_found)
   {
      std::cout << "(TPV config TOMLs not found; config sign-chain skipped)\n";
   }

   // Mesh-backed seed sign chain (skipped if no BP5 mesh).
   auto fix = BuildFixture();
   if (all_found && fix->loaded && fix->nf > 0)
   {
      for (size_t i = 0; i < cfgs.size(); ++i)
      { CheckSeedFromConfig(*fix, cfgs[i], cases[i].tag); }
   }
   else
   {
      std::cout << "(BP5 mesh / fault DOFs unavailable; ComputeParamsFaultLocal "
                   "seeding sub-checks skipped)\n";
   }

   std::cout << "\n========================================\n";
   std::cout << "Phase 8 stress-sign: " << num_passed << " / " << num_tests
             << " passed, " << num_failed << " failed\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
