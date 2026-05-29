// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_resolve_rate_state_guards.cpp — Phase 3 unit tests for
// SpatialFrictionResolver::ResolveRateState guards.
//
// Plan: PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24, §9.3
// "Acceptance Criteria".
//
//   G_R003  no double pore-pressure: a PorePressureSpec{} (zero) leaves the
//           already-effective sigma_n untouched; a non-zero pp double-
//           subtracts P_p (documents the driver-call-site trap).
//   G_R006  per-DOF b accepted (Phase 11a); per-DOF f_0/V_0 still abort.
//   G_R011  a > b (velocity-strengthening) is allowed; a <= 0 still aborts.
//   G_R001  (Phase 11 review, CRITICAL) per-DOF b under SRW is rejected: a
//           depth_profile or a per-DOF b spatial rule aborts when
//           state_evolution=slip_law_strong_rate_weakening (the SRW iterator
//           evolves psi with the scalar b_default while the equilibrium seed
//           would use per-DOF b -> silent t=0 disequilibrium).  Scalar-b SRW
//           still resolves.

#include "mfem.hpp"

#include "../../spatial/code/spatial_friction.hpp"
#include "../../dynamic/heterogeneous_material.hpp"

#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <functional>
#include <iostream>
#include <vector>

#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

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

namespace
{
// 1-element box mesh so MaterialField::EvalAt has a valid transformation.
class TinyMeshHolder
{
public:
   TinyMeshHolder() : mesh_(mfem::Mesh::MakeCartesian3D(
                              1, 1, 1, mfem::Element::HEXAHEDRON,
                              1.0, 1.0, 1.0)) {}
   mfem::Mesh& mesh() { return mesh_; }
private:
   mfem::Mesh mesh_;
};

void make_dofs(int N, Vector& coords, Array<int>& attr, Array<int>& elem)
{
   coords.SetSize(3 * N);
   attr.SetSize(N);
   elem.SetSize(N);
   for (int i = 0; i < N; ++i)
   {
      coords(3 * i + 0) = 1.0;
      coords(3 * i + 1) = 0.0;
      coords(3 * i + 2) = -1000.0;   // 1 km depth
      attr[i] = 101;
      elem[i] = 0;
   }
}

bool RunInChild(const std::function<void()>& body)
{
   ::fflush(stdout);
   ::fflush(stderr);
   const pid_t pid = ::fork();
   if (pid < 0) { return false; }
   if (pid == 0)
   {
      ::freopen("/dev/null", "w", stderr);
      try { body(); }
      catch (...) { ::_exit(1); }
      ::_exit(0);
   }
   int status = 0;
   while (::waitpid(pid, &status, 0) < 0 && errno == EINTR) { /* retry */ }
   if (WIFEXITED(status))   { return WEXITSTATUS(status) != 0; }
   if (WIFSIGNALED(status)) { return true; }
   return true;
}
}  // namespace

// =====================================================================
// R-003: no double pore-pressure subtraction.
// =====================================================================
static void G_R003_no_double_pp()
{
   std::cout << "\n[R-003] ResolveRateState does not double-subtract P_p\n";
   const int N = 1;
   Vector coords; Array<int> attr, elem;
   make_dofs(N, coords, attr, elem);
   TinyMeshHolder mh;
   SpatialFrictionResolver R;
   RateStateBlock cfg;          // defaults: a=0.010, b=0.015 (a<b)
   auto mat = MaterialField::MakeConstant(32.0e9, 32.0e9, 2670.0);

   Vector sn_total(N);
   sn_total(0) = 49.27e6;       // ALREADY effective (as geom.sigma_n_per_dof())

   // Correct driver call: zero pp -> effective stress unchanged.
   RateStatePerDOFParams p =
      R.ResolveRateState(cfg, coords, elem, attr, mat, mh.mesh(),
                         PorePressureSpec{}, sn_total);
   TEST_NEAR(p.sigma_n_eff(0), 49.27e6, 1.0,
             "PorePressureSpec{} leaves the already-effective sigma_n intact");

   // The trap: a non-zero pp double-subtracts P_p (49.27 - 16 = 33.27 MPa).
   PorePressureSpec pp_bad;
   pp_bad.P_p_pa = 16.0e6;
   RateStatePerDOFParams p_bad =
      R.ResolveRateState(cfg, coords, elem, attr, mat, mh.mesh(),
                         pp_bad, sn_total);
   TEST_NEAR(p_bad.sigma_n_eff(0), 33.27e6, 1.0,
             "non-zero pp double-subtracts P_p (documents the trap)");
}

// =====================================================================
// R-006 (Phase 11a relaxed): per-DOF b is now ACCEPTED (DOFData.b routes it
// through the iterator + seed); per-DOF f_0 / V_0 still abort (scalar aging-law
// globals; V_0 pinned by R-009 to FrictionSolver::V0).
// =====================================================================
static void G_R006_per_dof_b_accepted_f0_v0_rejected()
{
   std::cout << "\n[R-006] per-DOF b accepted; per-DOF f_0/V_0 rejected\n";
   const int N = 1;
   Vector coords; Array<int> attr, elem;
   make_dofs(N, coords, attr, elem);

   // per-DOF f_0 override -> abort.
   const bool f0_aborts = RunInChild([&]() {
      TinyMeshHolder mh;
      SpatialFrictionResolver R;
      RateStateBlock cfg;
      SpatialRule rule;
      rule.kind = SpatialRule::Kind::Depth;
      rule.f_0  = 0.5;                          // per-DOF f_0 override -> reject
      cfg.spatial.push_back(rule);
      auto mat = MaterialField::MakeConstant(32.0e9, 32.0e9, 2670.0);
      Vector sn_total(N); sn_total(0) = 49.27e6;
      (void) R.ResolveRateState(cfg, coords, elem, attr, mat, mh.mesh(),
                                PorePressureSpec{}, sn_total);
   });
   TEST_ASSERT(f0_aborts, "per-DOF f_0 override aborts (R-006)");

   // per-DOF V_0 override -> abort.
   const bool v0_aborts = RunInChild([&]() {
      TinyMeshHolder mh;
      SpatialFrictionResolver R;
      RateStateBlock cfg;
      SpatialRule rule;
      rule.kind = SpatialRule::Kind::Depth;
      rule.V_0  = 2.0e-6;                        // per-DOF V_0 override -> reject
      cfg.spatial.push_back(rule);
      auto mat = MaterialField::MakeConstant(32.0e9, 32.0e9, 2670.0);
      Vector sn_total(N); sn_total(0) = 49.27e6;
      (void) R.ResolveRateState(cfg, coords, elem, attr, mat, mh.mesh(),
                                PorePressureSpec{}, sn_total);
   });
   TEST_ASSERT(v0_aborts, "per-DOF V_0 override aborts (R-006)");

   // Phase 11a: a rule overriding a / b / Dc resolves cleanly and b is applied.
   TinyMeshHolder mh;
   SpatialFrictionResolver R;
   RateStateBlock cfg;
   SpatialRule rule;
   rule.kind = SpatialRule::Kind::Depth;
   rule.a    = 0.012;   // allowed per-DOF override
   rule.b    = 0.020;   // Phase 11a: per-DOF b now allowed (was rejected)
   rule.Dc   = 0.6;     // allowed per-DOF override
   cfg.spatial.push_back(rule);
   auto mat = MaterialField::MakeConstant(32.0e9, 32.0e9, 2670.0);
   Vector sn_total(N); sn_total(0) = 49.27e6;
   RateStatePerDOFParams p =
      R.ResolveRateState(cfg, coords, elem, attr, mat, mh.mesh(),
                         PorePressureSpec{}, sn_total);
   TEST_NEAR(p.a(0),  0.012, 1e-15, "per-DOF a override applied (control)");
   TEST_NEAR(p.b(0),  0.020, 1e-15, "per-DOF b override applied (Phase 11a)");
   TEST_NEAR(p.Dc(0), 0.6,   1e-15, "per-DOF Dc override applied (control)");
}

// =====================================================================
// R-011: a > b (velocity-strengthening) allowed; a <= 0 still aborts.
// =====================================================================
static void G_R011_allows_velocity_strengthening()
{
   std::cout << "\n[R-011] ResolveRateState allows a > b (velocity-strengthening)\n";
   const int N = 1;
   Vector coords; Array<int> attr, elem;
   make_dofs(N, coords, attr, elem);

   // a > b must NOT abort.
   TinyMeshHolder mh;
   SpatialFrictionResolver R;
   RateStateBlock cfg;
   cfg.a_default = 0.020;   // a > b
   cfg.b_default = 0.015;
   auto mat = MaterialField::MakeConstant(32.0e9, 32.0e9, 2670.0);
   Vector sn_total(N); sn_total(0) = 49.27e6;
   RateStatePerDOFParams p =
      R.ResolveRateState(cfg, coords, elem, attr, mat, mh.mesh(),
                         PorePressureSpec{}, sn_total);
   TEST_NEAR(p.a(0), 0.020, 1e-15, "a > b resolves without aborting (R-011)");
   TEST_NEAR(p.b(0), 0.015, 1e-15, "b unchanged under a > b");

   // Control: a <= 0 still aborts (positivity check kept).
   const bool a_nonpos_aborts = RunInChild([&]() {
      TinyMeshHolder mh2;
      SpatialFrictionResolver R2;
      RateStateBlock cfg2;
      cfg2.a_default = -0.001;   // a <= 0
      auto mat2 = MaterialField::MakeConstant(32.0e9, 32.0e9, 2670.0);
      Vector sn2(N); sn2(0) = 49.27e6;
      (void) R2.ResolveRateState(cfg2, coords, elem, attr, mat2, mh2.mesh(),
                                 PorePressureSpec{}, sn2);
   });
   TEST_ASSERT(a_nonpos_aborts, "a <= 0 still aborts (positivity kept)");
}

// =====================================================================
// R-001 (Phase 11 review, CRITICAL): per-DOF b is wired through the AGING
// iterator + equilibrium seed only.  The SRW policy evolves psi with the
// scalar b_default, so a non-scalar b (depth_profile or a per-DOF b spatial
// rule) under state_evolution=slip_law_strong_rate_weakening must be rejected
// loudly — otherwise the seed uses per-DOF b while the dynamics use scalar b
// (silent t=0 disequilibrium).  Scalar-b SRW must still resolve.
// =====================================================================
static void G_R001_srw_rejects_per_dof_b()
{
   std::cout << "\n[R-001] SRW + per-DOF b aborts; scalar-b SRW resolves\n";
   const int N = 1;
   Vector coords; Array<int> attr, elem;
   make_dofs(N, coords, attr, elem);

   // SRW + depth_profile -> abort (the guard fires before any profile eval).
   const bool profile_aborts = RunInChild([&]() {
      TinyMeshHolder mh;
      SpatialFrictionResolver R;
      RateStateBlock cfg;
      cfg.state_evolution = StateEvolutionKind::SlipLawStrongRateWeakening;
      cfg.depth_profile.enabled = true;          // per-DOF b source
      auto mat = MaterialField::MakeConstant(32.0e9, 32.0e9, 2670.0);
      Vector sn_total(N); sn_total(0) = 49.27e6;
      (void) R.ResolveRateState(cfg, coords, elem, attr, mat, mh.mesh(),
                                PorePressureSpec{}, sn_total);
   });
   TEST_ASSERT(profile_aborts, "SRW + depth_profile aborts (R-001)");

   // SRW + a per-DOF b spatial rule -> abort.
   const bool rule_aborts = RunInChild([&]() {
      TinyMeshHolder mh;
      SpatialFrictionResolver R;
      RateStateBlock cfg;
      cfg.state_evolution = StateEvolutionKind::SlipLawStrongRateWeakening;
      SpatialRule rule;
      rule.kind = SpatialRule::Kind::Depth;
      rule.b    = 0.020;                          // per-DOF b override
      cfg.spatial.push_back(rule);
      auto mat = MaterialField::MakeConstant(32.0e9, 32.0e9, 2670.0);
      Vector sn_total(N); sn_total(0) = 49.27e6;
      (void) R.ResolveRateState(cfg, coords, elem, attr, mat, mh.mesh(),
                                PorePressureSpec{}, sn_total);
   });
   TEST_ASSERT(rule_aborts, "SRW + per-DOF b spatial rule aborts (R-001)");

   // Control: SRW with SCALAR b (no profile, no b rule) must still resolve, and
   // a per-DOF a / Dc rule (NOT b) is still allowed under SRW.
   TinyMeshHolder mh;
   SpatialFrictionResolver R;
   RateStateBlock cfg;
   cfg.state_evolution = StateEvolutionKind::SlipLawStrongRateWeakening;
   cfg.V_w_default = 0.1;                          // SRW per-DOF loop needs V_w>0
   SpatialRule rule;
   rule.kind = SpatialRule::Kind::Depth;
   rule.a    = 0.012;                              // per-DOF a is fine under SRW
   cfg.spatial.push_back(rule);
   auto mat = MaterialField::MakeConstant(32.0e9, 32.0e9, 2670.0);
   Vector sn_total(N); sn_total(0) = 49.27e6;
   RateStatePerDOFParams p =
      R.ResolveRateState(cfg, coords, elem, attr, mat, mh.mesh(),
                         PorePressureSpec{}, sn_total);
   TEST_NEAR(p.b(0), cfg.b_default, 1e-15,
             "scalar-b SRW resolves; b == b_default (guard does not over-reject)");
   TEST_NEAR(p.a(0), 0.012, 1e-15,
             "per-DOF a override still allowed under SRW (only b is rejected)");
}

int main(int /*argc*/, char** /*argv*/)
{
   std::cout << "Running Phase 3 test_resolve_rate_state_guards\n";
   G_R003_no_double_pp();
   G_R006_per_dof_b_accepted_f0_v0_rejected();
   G_R011_allows_velocity_strengthening();
   G_R001_srw_rejects_per_dof_b();

   std::cout << "\n========================================\n";
   std::cout << "Phase 3 test_resolve_rate_state_guards: "
             << num_passed << " / " << num_tests
             << " passed, " << num_failed << " failed\n";
   std::cout << "========================================\n";
   return num_failed == 0 ? 0 : 1;
}
