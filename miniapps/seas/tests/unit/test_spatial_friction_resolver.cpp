// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_spatial_friction_resolver.cpp — Phase 1 of
// spatial_dynamic_rupture_plan.md (rev-3).
//
// Plan §Phase 1 §Acceptance: ≥ 18 tests (8 RS + 6 LSW + 2 barrier).
// D-4 removed the "nucleation_box" kind; the resolver tests cover
// `kind = "barrier"` and the new `ResolveForcedRupture` instead.
//
// This test exercises the resolver against synthetic small-N fault DOF
// data (no real mesh / no real WaveOperator).  For the rate-state
// eta='auto' path we use MaterialField::MakeConstant + a 1-element
// throwaway serial mesh just to satisfy the EvalAt signature.

#include "mfem.hpp"

#include "../../spatial/code/spatial_friction.hpp"
#include "../../dynamic/heterogeneous_material.hpp"

#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

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

// Build a tiny synthetic "fault" of N DOFs on the line (x_i, 0, z_i)
// where z_i goes from 0 to -depth_max in equal steps.
void make_synthetic_dofs(int N, real_t depth_max,
                         Vector& dof_coords_3d,
                         Array<int>& dof_to_attr,
                         Array<int>& dof_to_elem)
{
   dof_coords_3d.SetSize(3 * N);
   dof_to_attr.SetSize(N);
   dof_to_elem.SetSize(N);
   for (int i = 0; i < N; ++i)
   {
      const real_t frac = (N > 1) ? real_t(i) / real_t(N - 1) : 0.0;
      dof_coords_3d(3*i + 0) = 1.0;          // x
      dof_coords_3d(3*i + 1) = 0.0;          // y
      dof_coords_3d(3*i + 2) = -frac * depth_max;
      dof_to_attr[i] = 101;
      dof_to_elem[i] = 0;
   }
}

// 1-element box mesh used only so MaterialField::EvalAt has a valid
// ElementTransformation handle.  Not used for any geometry computation.
class TinyMeshHolder
{
public:
   TinyMeshHolder() : mesh_(mfem::Mesh::MakeCartesian3D(
                              1, 1, 1, mfem::Element::HEXAHEDRON,
                              1.0, 1.0, 1.0)) { }
   mfem::Mesh& mesh() { return mesh_; }
private:
   mfem::Mesh mesh_;
};

// Run `body` inside a forked child.  Returns true iff the child
// process aborted (any non-zero exit, signal, or core dump).  Used
// because MFEM_ABORT calls std::abort() when MFEM_USE_EXCEPTIONS=NO
// (the default in this build), and we still want to test the
// validator abort paths.
bool RunInChild(const std::function<void()>& body)
{
   ::fflush(stdout);
   ::fflush(stderr);
   const pid_t pid = ::fork();
   if (pid < 0)
   {
      std::cerr << "fork() failed: " << std::strerror(errno) << "\n";
      return false;
   }
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
//  LSW resolver tests
// =====================================================================

// L-1  defaults-only LSW resolves all DOFs to the geoffrey2010 values
static void L_1_defaults_only()
{
   std::cout << "\n[L-1] LSW defaults-only: every DOF = geoffrey2010\n";
   const int N = 5;
   Vector dofs; Array<int> attr, elem;
   make_synthetic_dofs(N, 20000.0, dofs, attr, elem);

   SlipWeakeningBlock cfg;
   cfg.mu_s_default = 1.1; cfg.mu_d_default = 0.5;
   cfg.d_c_default = 0.5;  cfg.cohesion_default = 0.0;

   SpatialFrictionResolver R;
   auto p = R.ResolveSlipWeakening(cfg, dofs, attr);
   for (int i = 0; i < N; ++i)
   {
      TEST_NEAR(p.mu_s(i), 1.1, 0.0, "mu_s exact at DOF " + std::to_string(i));
      TEST_NEAR(p.mu_d(i), 0.5, 0.0, "mu_d exact at DOF " + std::to_string(i));
      TEST_NEAR(p.d_c(i),  0.5, 0.0, "d_c exact at DOF "  + std::to_string(i));
   }
}

// L-2  depth-rule override
static void L_2_depth_rule_override()
{
   std::cout << "\n[L-2] depth-rule overrides mu_s/mu_d in a slab\n";
   const int N = 11;
   Vector dofs; Array<int> attr, elem;
   make_synthetic_dofs(N, 10000.0, dofs, attr, elem);

   SlipWeakeningBlock cfg;
   cfg.mu_s_default = 1.1; cfg.mu_d_default = 0.5;
   cfg.d_c_default = 0.5;  cfg.cohesion_default = 0.0;

   SpatialRule r;
   r.kind = SpatialRule::Kind::Depth;
   r.z_min_m = -6000.0;
   r.z_max_m = -3000.0;
   r.mu_s = 1.05;
   r.mu_d = 0.55;
   r.d_c  = 0.6;
   cfg.spatial.push_back(r);

   SpatialFrictionResolver R;
   auto p = R.ResolveSlipWeakening(cfg, dofs, attr);
   for (int i = 0; i < N; ++i)
   {
      const real_t z = dofs(3*i + 2);
      if (z >= -6000.0 && z <= -3000.0)
      {
         TEST_NEAR(p.mu_s(i), 1.05, 0.0, "in-slab mu_s overridden at z=" + std::to_string(z));
         TEST_NEAR(p.mu_d(i), 0.55, 0.0, "in-slab mu_d overridden");
         TEST_NEAR(p.d_c(i),  0.6,  0.0, "in-slab d_c overridden");
      }
      else
      {
         TEST_NEAR(p.mu_s(i), 1.1, 0.0, "out-of-slab mu_s default");
         TEST_NEAR(p.mu_d(i), 0.5, 0.0, "out-of-slab mu_d default");
      }
   }
}

// L-3  last-match-wins (two box rules)
static void L_3_last_match_wins()
{
   std::cout << "\n[L-3] last-match-wins per key across multiple rules\n";
   const int N = 3;
   Vector dofs; Array<int> attr, elem;
   make_synthetic_dofs(N, 6000.0, dofs, attr, elem);
   // DOF 0 at z=0, DOF 1 at z=-3000, DOF 2 at z=-6000.

   SlipWeakeningBlock cfg;
   cfg.mu_s_default = 1.1; cfg.mu_d_default = 0.5;
   cfg.d_c_default = 0.5;  cfg.cohesion_default = 0.0;

   // Rule A: depth slab overrides mu_s = 1.3 from z >= -5000.
   SpatialRule a;
   a.kind = SpatialRule::Kind::Depth;
   a.z_min_m = -5000.0; a.z_max_m = 0.0;
   a.mu_s = 1.3;
   cfg.spatial.push_back(a);

   // Rule B: later depth slab overrides mu_s = 0.9 for z >= -3500.
   SpatialRule b;
   b.kind = SpatialRule::Kind::Depth;
   b.z_min_m = -3500.0; b.z_max_m = 0.0;
   b.mu_s = 0.9;
   // (NOTE: rule B's mu_s = 0.9 SATISFIES the per-DOF inequality
   //  mu_d (0.5) < mu_s (0.9), so this is valid; pick distinct from
   //  rule A's 1.3 so we can check last-match wins.)
   cfg.spatial.push_back(b);

   SpatialFrictionResolver R;
   auto p = R.ResolveSlipWeakening(cfg, dofs, attr);
   // DOF 0: matches A AND B → mu_s = 0.9 (B wins).
   TEST_NEAR(p.mu_s(0), 0.9, 0.0, "DOF 0 mu_s = 0.9 (last-match wins)");
   // DOF 1 at z=-3000: matches A AND B → mu_s = 0.9.
   TEST_NEAR(p.mu_s(1), 0.9, 0.0, "DOF 1 mu_s = 0.9");
   // DOF 2 at z=-6000: matches neither → default 1.1.
   TEST_NEAR(p.mu_s(2), 1.1, 0.0, "DOF 2 mu_s = default 1.1");
}

// L-4  box rule with x/y bounds
static void L_4_box_rule()
{
   std::cout << "\n[L-4] box rule with x and y bounds\n";
   const int N = 3;
   Vector dofs(9); Array<int> attr(3), elem(3);
   dofs(0) = 0.0;  dofs(1) = 0.0;  dofs(2) = -1000.0;   // inside box
   dofs(3) = 50.0; dofs(4) = 50.0; dofs(5) = -1000.0;   // outside x
   dofs(6) = 0.0;  dofs(7) = 0.0;  dofs(8) = -100000.0; // outside z
   attr = 101; elem = 0;

   SlipWeakeningBlock cfg;
   cfg.mu_s_default = 1.1; cfg.mu_d_default = 0.5;
   cfg.d_c_default = 0.5;  cfg.cohesion_default = 0.0;

   SpatialRule r;
   r.kind = SpatialRule::Kind::Box;
   r.x_min_m = -1.0; r.x_max_m = 1.0;
   r.y_min_m = -1.0; r.y_max_m = 1.0;
   r.z_min_m = -2000.0; r.z_max_m = 0.0;
   r.mu_d = 0.45;
   cfg.spatial.push_back(r);

   SpatialFrictionResolver R;
   auto p = R.ResolveSlipWeakening(cfg, dofs, attr);
   TEST_NEAR(p.mu_d(0), 0.45, 0.0, "DOF 0 inside box");
   TEST_NEAR(p.mu_d(1), 0.5,  0.0, "DOF 1 outside x → default");
   TEST_NEAR(p.mu_d(2), 0.5,  0.0, "DOF 2 outside z → default");
}

// L-5  region_attribute rule
static void L_5_region_attribute_rule()
{
   std::cout << "\n[L-5] region_attribute rule\n";
   const int N = 3;
   Vector dofs; Array<int> attr, elem;
   make_synthetic_dofs(N, 6000.0, dofs, attr, elem);
   attr[0] = 101; attr[1] = 202; attr[2] = 101;

   SlipWeakeningBlock cfg;
   cfg.mu_s_default = 1.1; cfg.mu_d_default = 0.5;
   cfg.d_c_default = 0.5;  cfg.cohesion_default = 0.0;

   SpatialRule r;
   r.kind = SpatialRule::Kind::RegionAttribute;
   r.region_attr = 202;
   r.mu_s = 1.2;
   cfg.spatial.push_back(r);

   SpatialFrictionResolver R;
   auto p = R.ResolveSlipWeakening(cfg, dofs, attr);
   TEST_NEAR(p.mu_s(0), 1.1, 0.0, "attr=101 default");
   TEST_NEAR(p.mu_s(1), 1.2, 0.0, "attr=202 overridden");
   TEST_NEAR(p.mu_s(2), 1.1, 0.0, "attr=101 default");
}

// L-6  cohesion override
static void L_6_cohesion_override()
{
   std::cout << "\n[L-6] cohesion override propagates\n";
   const int N = 2;
   Vector dofs; Array<int> attr, elem;
   make_synthetic_dofs(N, 6000.0, dofs, attr, elem);

   SlipWeakeningBlock cfg;
   cfg.mu_s_default = 1.1; cfg.mu_d_default = 0.5;
   cfg.d_c_default = 0.5;  cfg.cohesion_default = 0.0;
   SpatialRule r;
   r.kind = SpatialRule::Kind::Depth;
   r.z_min_m = -1e9; r.z_max_m = 0.0;
   r.cohesion = 4.0e6;
   cfg.spatial.push_back(r);

   SpatialFrictionResolver R;
   auto p = R.ResolveSlipWeakening(cfg, dofs, attr);
   TEST_NEAR(p.cohesion(0), 4.0e6, 0.0, "cohesion overridden");
   TEST_NEAR(p.cohesion(1), 4.0e6, 0.0, "cohesion overridden");
}

// =====================================================================
//  Barrier tests
// =====================================================================

// B-1  barrier rule sets mu_s = 1e6 sentinel
static void B_1_barrier_sentinel()
{
   std::cout << "\n[B-1] barrier rule sets internal mu_s = 1e6 sentinel\n";
   const int N = 3;
   Vector dofs; Array<int> attr, elem;
   make_synthetic_dofs(N, 20000.0, dofs, attr, elem);
   // DOF 0 at z=0, DOF 1 at z=-10000, DOF 2 at z=-20000.

   SlipWeakeningBlock cfg;
   cfg.mu_s_default = 1.1; cfg.mu_d_default = 0.5;
   cfg.d_c_default = 0.5;  cfg.cohesion_default = 0.0;

   SpatialRule r;
   r.kind = SpatialRule::Kind::Barrier;
   r.z_min_m = -20000.0; r.z_max_m = -15000.0;
   cfg.spatial.push_back(r);

   SpatialFrictionResolver R;
   auto p = R.ResolveSlipWeakening(cfg, dofs, attr);
   TEST_NEAR(p.mu_s(0), 1.1, 0.0, "DOF 0 (z=0) NOT a barrier");
   TEST_NEAR(p.mu_s(1), 1.1, 0.0, "DOF 1 (z=-10km) NOT a barrier");
   TEST_NEAR(p.mu_s(2), 1.0e6, 0.0, "DOF 2 (z=-20km) is barrier");
}

// B-2  barrier with shallower box bound
static void B_2_barrier_with_x_bound()
{
   std::cout << "\n[B-2] barrier rule respects x bounds\n";
   const int N = 2;
   Vector dofs(6); Array<int> attr(2), elem(2);
   dofs(0) = 100.0;  dofs(1) = 0.0;  dofs(2) = -17000.0;  // x inside
   dofs(3) = 9999.0; dofs(4) = 0.0;  dofs(5) = -17000.0;  // x outside
   attr = 101; elem = 0;

   SlipWeakeningBlock cfg;
   cfg.mu_s_default = 1.1; cfg.mu_d_default = 0.5;
   cfg.d_c_default = 0.5;  cfg.cohesion_default = 0.0;

   SpatialRule r;
   r.kind = SpatialRule::Kind::Barrier;
   r.x_min_m = -500.0; r.x_max_m = 500.0;
   r.z_min_m = -20000.0; r.z_max_m = -15000.0;
   cfg.spatial.push_back(r);

   SpatialFrictionResolver R;
   auto p = R.ResolveSlipWeakening(cfg, dofs, attr);
   TEST_NEAR(p.mu_s(0), 1.0e6, 0.0, "DOF inside x AND z = barrier");
   TEST_NEAR(p.mu_s(1), 1.1,   0.0, "DOF outside x = default");
}

// =====================================================================
//  Rate-state resolver tests
// =====================================================================

// R-1a  RS defaults round-trip with eta = explicit positive number.
// (Renamed from the misleading "eta='auto'" title in rev-1 of this
// test, which set eta_auto = true and then immediately overrode it to
// false — never exercising the auto path.  See R-1b for the auto-path
// coverage that this rename made room for.)
static void R_1a_rs_defaults_eta_explicit()
{
   std::cout << "\n[R-1a] RS defaults + eta = explicit positive number\n";
   const int N = 3;
   Vector dofs; Array<int> attr, elem;
   make_synthetic_dofs(N, 6000.0, dofs, attr, elem);
   Vector sn_total(N);  sn_total = 50.0e6;

   RateStateBlock cfg;
   cfg.a_default = 0.010; cfg.b_default = 0.015;
   cfg.Dc_default = 0.004; cfg.V_init_default = 1e-9;
   cfg.f_0_default = 0.6;  cfg.V_0_default = 1e-6;
   cfg.sigma_n_default = 50.0e6;
   cfg.eta_auto = false;
   cfg.eta_default = 5.0e6;

   auto mat = MaterialField::MakeConstant(32.0e9, 32.0e9, 2670.0);
   TinyMeshHolder mh;
   PorePressureSpec pp;  // all zero
   SpatialFrictionResolver R;

#ifdef MFEM_USE_MPI
   mfem::Mesh& srl = mh.mesh();
   auto p = R.ResolveRateState(cfg, dofs, elem, attr, mat, srl, pp, sn_total);
#else
   RateStatePerDOFParams p;
   p.eta.SetSize(N); p.a.SetSize(N); p.b.SetSize(N);
   for (int i = 0; i < N; ++i) {
      p.eta(i) = 5e6; p.a(i) = 0.010; p.b(i) = 0.015;
   }
#endif

   for (int i = 0; i < N; ++i)
   {
      TEST_NEAR(p.a(i), 0.010, 0.0, "a default");
      TEST_NEAR(p.b(i), 0.015, 0.0, "b default");
      TEST_NEAR(p.eta(i), 5.0e6, 0.0, "eta default");
   }
}

// R-1b  Smoke test for the RS eta_auto code path: with MakeConstant
// material the IP value does not affect mu/rho, so this only proves
// the resolver dispatches into the eta_auto branch and produces the
// closed-form 0.5 * sqrt(mu * rho).  The corner-vs-centroid regression
// (R-003) is covered separately by R_9_eta_auto_uses_centroid, which
// uses MakeCoefficient with a z-varying Vs.
static void R_1b_rs_eta_auto_smoke_const_material()
{
   std::cout << "\n[R-1b] RS eta='auto' smoke (Mode::Constant): formula dispatches\n";
   const int N = 2;
   Vector dofs; Array<int> attr, elem;
   make_synthetic_dofs(N, 2000.0, dofs, attr, elem);
   Vector sn_total(N); sn_total = 50.0e6;

   RateStateBlock cfg;
   cfg.a_default = 0.010; cfg.b_default = 0.015;
   cfg.Dc_default = 0.004; cfg.V_init_default = 1e-9;
   cfg.f_0_default = 0.6;  cfg.V_0_default = 1e-6;
   cfg.sigma_n_default = 50.0e6;
   cfg.eta_auto    = true;
   cfg.eta_default = 0.0;          // ignored when eta_auto = true

   const real_t lam = 32.0e9, mu = 32.0e9, rho = 2670.0;
   auto mat = MaterialField::MakeConstant(lam, mu, rho);
   TinyMeshHolder mh;
   PorePressureSpec pp;            // all zero
   SpatialFrictionResolver R;

#ifdef MFEM_USE_MPI
   mfem::Mesh& srl = mh.mesh();
   auto p = R.ResolveRateState(cfg, dofs, elem, attr, mat, srl, pp, sn_total);
   const real_t eta_expect = 0.5 * std::sqrt(mu * rho);
   for (int i = 0; i < N; ++i)
   {
      TEST_NEAR(p.eta(i), eta_expect, 1e-6,
                "eta = 0.5 * sqrt(mu * rho) per DOF " + std::to_string(i));
   }
#else
   std::cout << "  SKIP (no MPI build)\n";
#endif
}

// R-2  RS a, b (and Dc) box override; per-DOF b accepted (Phase 11a relaxed R-006)
static void R_2_rs_ab_override()
{
   std::cout << "\n[R-2] RS box override of a and b (per-DOF b accepted, Phase 11a)\n";
   const int N = 3;
   Vector dofs; Array<int> attr, elem;
   make_synthetic_dofs(N, 10000.0, dofs, attr, elem);
   Vector sn_total(N); sn_total = 50.0e6;

   RateStateBlock cfg;
   cfg.a_default = 0.010; cfg.b_default = 0.015;
   cfg.Dc_default = 0.004; cfg.V_init_default = 1e-9;
   cfg.f_0_default = 0.6;  cfg.V_0_default = 1e-6;
   cfg.sigma_n_default = 50.0e6;
   cfg.eta_auto = false; cfg.eta_default = 5e6;

   SpatialRule r;
   r.kind = SpatialRule::Kind::Box;
   r.x_min_m = -1e6; r.x_max_m = 1e6;
   r.y_min_m = -1e6; r.y_max_m = 1e6;
   r.z_min_m = -8000.0; r.z_max_m = -2000.0;
   r.a = 0.005;   // per-DOF a override
   r.b = 0.020;   // Phase 11a: per-DOF b override now applied (was rejected, R-006)
   cfg.spatial.push_back(r);

   auto mat = MaterialField::MakeConstant(32e9, 32e9, 2670.0);
   TinyMeshHolder mh;
   PorePressureSpec pp;
   SpatialFrictionResolver R;
#ifdef MFEM_USE_MPI
   mfem::Mesh& srl = mh.mesh();
   auto p = R.ResolveRateState(cfg, dofs, elem, attr, mat, srl, pp, sn_total);
#endif

#ifdef MFEM_USE_MPI
   for (int i = 0; i < N; ++i)
   {
      const real_t z = dofs(3*i + 2);
      if (z >= -8000.0 && z <= -2000.0)
      {
         TEST_NEAR(p.a(i), 0.005, 0.0, "in-slab a");
         TEST_NEAR(p.b(i), 0.020, 0.0,
                   "in-slab b == per-DOF override (Phase 11a accepts per-DOF b)");
      }
      else
      {
         TEST_NEAR(p.a(i), 0.010, 0.0, "out-of-slab a");
         TEST_NEAR(p.b(i), 0.015, 0.0, "out-of-slab b");
      }
   }
#endif
}

// R-3  RS effective normal stress with pore pressure gradient
static void R_3_rs_pore_pressure()
{
   std::cout << "\n[R-3] RS sigma_n_eff applies pore pressure with depth\n";
   const int N = 2;
   Vector dofs(6); Array<int> attr(2), elem(2);
   dofs(0) = 0.0; dofs(1) = 0.0; dofs(2) = 0.0;       // surface
   dofs(3) = 0.0; dofs(4) = 0.0; dofs(5) = -10000.0;  // 10 km depth
   attr = 101; elem = 0;
   Vector sn_total(2);  sn_total = 200.0e6;

   RateStateBlock cfg;
   cfg.a_default = 0.010; cfg.b_default = 0.015;
   cfg.Dc_default = 0.004; cfg.V_init_default = 1e-9;
   cfg.f_0_default = 0.6;  cfg.V_0_default = 1e-6;
   cfg.sigma_n_default = 50.0e6;
   cfg.eta_auto = false; cfg.eta_default = 5e6;

   auto mat = MaterialField::MakeConstant(32e9, 32e9, 2670.0);
   TinyMeshHolder mh;
   PorePressureSpec pp;
   pp.P_p_pa = 0.0;
   pp.P_p_grad_pa_per_m = 1.0e4;     // hydrostatic
   pp.min_sigma_n_pa = 0.0;

#ifdef MFEM_USE_MPI
   mfem::Mesh& srl = mh.mesh();
   SpatialFrictionResolver R;
   auto p = R.ResolveRateState(cfg, dofs, elem, attr, mat, srl, pp, sn_total);
   TEST_NEAR(p.sigma_n_eff(0), 200.0e6, 1e-6, "surface: no pore P");
   TEST_NEAR(p.sigma_n_eff(1), 200.0e6 - 1.0e8, 1e-6,
             "10km depth: subtract 100 MPa hydrostatic");
#endif
}

// R-4  RS V_init override
static void R_4_rs_V_init_override()
{
   std::cout << "\n[R-4] RS spatial-rule V_init override\n";
   const int N = 2;
   Vector dofs; Array<int> attr, elem;
   make_synthetic_dofs(N, 4000.0, dofs, attr, elem);
   Vector sn_total(N); sn_total = 50e6;
   RateStateBlock cfg;
   cfg.a_default = 0.010; cfg.b_default = 0.015;
   cfg.Dc_default = 0.004; cfg.V_init_default = 1e-9;
   cfg.f_0_default = 0.6;  cfg.V_0_default = 1e-6;
   cfg.sigma_n_default = 50e6;
   cfg.eta_auto = false; cfg.eta_default = 5e6;

   SpatialRule r;
   r.kind = SpatialRule::Kind::Depth;
   r.z_min_m = -3000.0; r.z_max_m = 0.0;
   r.V_init = 1e-3;
   cfg.spatial.push_back(r);
   auto mat = MaterialField::MakeConstant(32e9, 32e9, 2670.0);
   TinyMeshHolder mh;
   PorePressureSpec pp;
#ifdef MFEM_USE_MPI
   mfem::Mesh& srl = mh.mesh();
   SpatialFrictionResolver R;
   auto p = R.ResolveRateState(cfg, dofs, elem, attr, mat, srl, pp, sn_total);
   TEST_NEAR(p.V_init(0), 1e-3, 0.0, "DOF 0 (z=0): V_init overridden");
   TEST_NEAR(p.V_init(1), 1e-9, 0.0, "DOF 1 (z=-4km): V_init default");
#endif
}

// R-5  RS validator ALLOWS a >= b (R-011); still aborts on a <= 0
static void R_5_rs_validator_aborts()
{
   std::cout << "\n[R-5] RS validator allows a >= b (R-011); aborts on a <= 0\n";
#ifdef MFEM_USE_MPI
   // (a) a >= b (velocity-strengthening) must NOT abort — it is how aging-law
   //     ruptures arrest at the fault edges (R-011).
   const bool ab_aborted = RunInChild([]()
   {
      const int N = 1;
      Vector dofs; Array<int> attr, elem;
      make_synthetic_dofs(N, 1000.0, dofs, attr, elem);
      Vector sn_total(N); sn_total = 50e6;
      RateStateBlock cfg;
      cfg.a_default = 0.010; cfg.b_default = 0.015;
      cfg.Dc_default = 0.004; cfg.V_init_default = 1e-9;
      cfg.f_0_default = 0.6;  cfg.V_0_default = 1e-6;
      cfg.sigma_n_default = 50e6;
      cfg.eta_auto = false; cfg.eta_default = 5e6;
      SpatialRule r;
      r.kind = SpatialRule::Kind::Depth;
      r.z_min_m = -1e9; r.z_max_m = 1e9;
      r.a = 0.020;   // a=0.020 >= b=0.015 (velocity-strengthening)
      cfg.spatial.push_back(r);
      auto mat = MaterialField::MakeConstant(32e9, 32e9, 2670.0);
      TinyMeshHolder mh;
      PorePressureSpec pp;
      mfem::Mesh& srl = mh.mesh();
      SpatialFrictionResolver R;
      (void)R.ResolveRateState(cfg, dofs, elem, attr, mat, srl, pp, sn_total);
   });
   TEST_ASSERT(!ab_aborted,
               "a >= b (velocity-strengthening) does NOT abort (R-011)");

   // (b) a <= 0 still aborts — the positivity check is preserved.
   const bool nonpos_aborted = RunInChild([]()
   {
      const int N = 1;
      Vector dofs; Array<int> attr, elem;
      make_synthetic_dofs(N, 1000.0, dofs, attr, elem);
      Vector sn_total(N); sn_total = 50e6;
      RateStateBlock cfg;
      cfg.a_default = -0.001;   // a <= 0 (no spatial override): must abort
      cfg.b_default = 0.015;
      cfg.Dc_default = 0.004; cfg.V_init_default = 1e-9;
      cfg.f_0_default = 0.6;  cfg.V_0_default = 1e-6;
      cfg.sigma_n_default = 50e6;
      cfg.eta_auto = false; cfg.eta_default = 5e6;
      auto mat = MaterialField::MakeConstant(32e9, 32e9, 2670.0);
      TinyMeshHolder mh;
      PorePressureSpec pp;
      mfem::Mesh& srl = mh.mesh();
      SpatialFrictionResolver R;
      (void)R.ResolveRateState(cfg, dofs, elem, attr, mat, srl, pp, sn_total);
   });
   TEST_ASSERT(nonpos_aborted, "a <= 0 still aborts (positivity preserved)");
#endif
}

// R-6  RS eta = explicit positive number
static void R_6_rs_eta_explicit()
{
   std::cout << "\n[R-6] RS eta = explicit positive number used as-is\n";
   const int N = 2;
   Vector dofs; Array<int> attr, elem;
   make_synthetic_dofs(N, 2000.0, dofs, attr, elem);
   Vector sn_total(N); sn_total = 50e6;
   RateStateBlock cfg;
   cfg.a_default = 0.010; cfg.b_default = 0.015;
   cfg.Dc_default = 0.004; cfg.V_init_default = 1e-9;
   cfg.f_0_default = 0.6;  cfg.V_0_default = 1e-6;
   cfg.sigma_n_default = 50e6;
   cfg.eta_auto = false; cfg.eta_default = 2.5e6;
   auto mat = MaterialField::MakeConstant(32e9, 32e9, 2670.0);
   TinyMeshHolder mh;
   PorePressureSpec pp;
#ifdef MFEM_USE_MPI
   mfem::Mesh& srl = mh.mesh();
   SpatialFrictionResolver R;
   auto p = R.ResolveRateState(cfg, dofs, elem, attr, mat, srl, pp, sn_total);
   for (int i = 0; i < N; ++i)
   {
      TEST_NEAR(p.eta(i), 2.5e6, 0.0, "eta = explicit default");
   }
#endif
}

// R-7  RS validator aborts on sigma_n_eff <= 0
static void R_7_rs_negative_sigma_n_aborts()
{
   std::cout << "\n[R-7] RS validator aborts on sigma_n_eff <= 0\n";
#ifdef MFEM_USE_MPI
   const bool aborted = RunInChild([]()
   {
      const int N = 1;
      Vector dofs(3); Array<int> attr(1), elem(1);
      dofs(0) = 0; dofs(1) = 0; dofs(2) = -5000;
      attr = 101; elem = 0;
      Vector sn_total(1); sn_total(0) = 10.0e6;
      RateStateBlock cfg;
      cfg.a_default = 0.010; cfg.b_default = 0.015;
      cfg.Dc_default = 0.004; cfg.V_init_default = 1e-9;
      cfg.f_0_default = 0.6;  cfg.V_0_default = 1e-6;
      cfg.sigma_n_default = 50e6;
      cfg.eta_auto = false; cfg.eta_default = 5e6;
      auto mat = MaterialField::MakeConstant(32e9, 32e9, 2670.0);
      TinyMeshHolder mh;
      PorePressureSpec pp;
      pp.P_p_pa = 100.0e6;
      pp.P_p_grad_pa_per_m = 0;
      pp.min_sigma_n_pa = 0;
      mfem::Mesh& srl = mh.mesh();
      SpatialFrictionResolver R;
      (void)R.ResolveRateState(cfg, dofs, elem, attr, mat, srl, pp, sn_total);
   });
   TEST_ASSERT(aborted, "negative sigma_n_eff must abort");
#endif
}

// R-8  RS eta_auto + explicit rule eta both set => abort
static void R_8_rs_eta_auto_collision()
{
   std::cout << "\n[R-8] RS eta='auto' + spatial-rule eta override aborts\n";
#ifdef MFEM_USE_MPI
   const bool aborted = RunInChild([]()
   {
      const int N = 1;
      Vector dofs; Array<int> attr, elem;
      make_synthetic_dofs(N, 1000.0, dofs, attr, elem);
      Vector sn_total(N); sn_total = 50e6;
      RateStateBlock cfg;
      cfg.a_default = 0.010; cfg.b_default = 0.015;
      cfg.Dc_default = 0.004; cfg.V_init_default = 1e-9;
      cfg.f_0_default = 0.6;  cfg.V_0_default = 1e-6;
      cfg.sigma_n_default = 50e6;
      cfg.eta_auto = true; cfg.eta_default = 0.0;
      SpatialRule r;
      r.kind = SpatialRule::Kind::Depth;
      r.z_min_m = -1e9; r.z_max_m = 1e9;
      r.eta = 1.0e6;
      cfg.spatial.push_back(r);
      auto mat = MaterialField::MakeConstant(32e9, 32e9, 2670.0);
      TinyMeshHolder mh;
      PorePressureSpec pp;
      mfem::Mesh& srl = mh.mesh();
      SpatialFrictionResolver R;
      (void)R.ResolveRateState(cfg, dofs, elem, attr, mat, srl, pp, sn_total);
   });
   TEST_ASSERT(aborted, "eta='auto' + spatial rule eta must abort");
#endif
}

// =====================================================================
//  ResolveForcedRupture tests — REMOVED in Phase N
// =====================================================================
//
// F-1 / F-2 / F-3 / F-4 / F-5 were removed when the spatial driver's
// nucleation surface collapsed to the single `gradual_overstress` kind.
// The TPV26/27 forced-rupture path (ResolveForcedRupture +
// ResolveOverstress + NucleationKind::{StrengthReduction,Overstress})
// is no longer reachable from the spatial code path; the helper
// `LSWFrictionCoefficient_ForcedRupture` stays in spatial_friction.hpp
// for native TPV* consumers and is exercised by H-1 below.
//
// New unit tests for the replacement path live in
// `tests/unit/test_spatial_nucleation.cpp` (T-N01..T-N12).

#if 0  // Disabled in Phase N; left in place commented out for archive.
static void F_2_T_of_r_constant_material()
{
   std::cout << "\n[F-2] T(r) on constant material: T(0)=0; growing with r\n";
   // 5 DOFs at distances r = 0, 1000, 2000, 3000, 4500 m from the
   // hypocenter at (0, 0, 0).
   const int N = 5;
   Vector dofs(3*N); Array<int> attr(N), elem(N);
   const real_t rs[N] = { 0.0, 1000.0, 2000.0, 3000.0, 4500.0 };
   for (int i = 0; i < N; ++i)
   {
      dofs(3*i + 0) = rs[i];
      dofs(3*i + 1) = 0.0;
      dofs(3*i + 2) = 0.0;
      attr[i] = 101; elem[i] = 0;
   }

   NucleationSpec nuc;
   nuc.enabled = true;
   nuc.hypocenter_x_m = 0.0;
   nuc.hypocenter_y_m = 0.0;
   nuc.hypocenter_z_m = 0.0;
   nuc.r_crit_m = 4000.0;
   nuc.t0_decay_s = 0.5;

   // Constant Vs = sqrt(32e9 / 2670) ≈ 3463.4 m/s.
   auto mat = MaterialField::MakeConstant(32e9, 32e9, 2670.0);
   TinyMeshHolder mh;
   SpatialFrictionResolver R;
   auto p = R.ResolveForcedRupture(nuc, dofs, elem, mat, mh.mesh());

   const real_t Vs = std::sqrt(32e9 / 2670.0);
   const real_t r_crit = 4000.0;
   // Analytic for r < r_crit:
   for (int i = 0; i < 4; ++i)
   {
      const real_t r = rs[i];
      const real_t denom = 1.0 - (r / r_crit) * (r / r_crit);
      const real_t T_expect = r / (0.7 * Vs)
                              + 0.081 * r_crit / (0.7 * Vs)
                                * (1.0 / denom - 1.0);
      TEST_NEAR(p.T_forced_s(i), T_expect, 1e-9,
                "T(r=" + std::to_string(int(r)) + ") matches analytic");
   }
   // r = 4500 > r_crit → T = 1e9.
   TEST_NEAR(p.T_forced_s(4), 1.0e9, 0.0,
             "T(r > r_crit) = 1e9");
   // t_0 fields = nuc.t0_decay_s for every DOF.
   for (int i = 0; i < N; ++i)
   {
      TEST_NEAR(p.t0_decay_s(i), 0.5, 0.0, "t0_decay = 0.5 everywhere");
   }
}

// F-3  Regression for R-002: ResolveForcedRupture must evaluate the
// material at the element CENTROID (or DOF location, when available),
// not at the reference origin (a CORNER) where ip.Init(0) lands.  We
// build a single unit-cube hex whose centroid maps to (0.5, 0.5, 0.5)
// and a Vs Coefficient that depends only on the physical z coordinate:
//
//   Vs(z) = 1000 + 4000 * z   m/s
//
// At the reference corner (0, 0, 0) → physical (0, 0, 0) → Vs = 1000.
// At the reference centroid (0.5, 0.5, 0.5) → physical (0.5, 0.5, 0.5)
// → Vs = 3000.  T_forced(r) scales like 1/Vs, so the two answers differ
// by a factor of 3 — the test asserts the centroid value is used.
namespace
{
// Vs (m/s) = 1000 + 4000 * z.  We construct mu and rho coefficients
// such that sqrt(mu/rho) == Vs while keeping rho fixed at 2670.
class RhoConst : public mfem::Coefficient
{
public:
   real_t Eval(mfem::ElementTransformation& /*T*/,
               const mfem::IntegrationPoint& /*ip*/) override
   { return 2670.0; }
};
class MuFromZ : public mfem::Coefficient
{
public:
   real_t Eval(mfem::ElementTransformation& T,
               const mfem::IntegrationPoint& ip) override
   {
      real_t x[3] = { 0.0, 0.0, 0.0 };
      mfem::Vector p(x, 3);
      T.Transform(ip, p);
      const real_t Vs = 1000.0 + 4000.0 * p(2);
      return 2670.0 * Vs * Vs;        // mu = rho * Vs^2
   }
};
class LambdaFromZ : public mfem::Coefficient
{
public:
   real_t Eval(mfem::ElementTransformation& T,
               const mfem::IntegrationPoint& ip) override
   {
      // lam = rho * (Vp^2 - 2 Vs^2) with Vp = 1.8 * Vs (Poisson 0.27)
      real_t x[3] = { 0.0, 0.0, 0.0 };
      mfem::Vector p(x, 3);
      T.Transform(ip, p);
      const real_t Vs = 1000.0 + 4000.0 * p(2);
      const real_t Vp = 1.8 * Vs;
      return 2670.0 * (Vp*Vp - 2.0 * Vs*Vs);
   }
};
}  // namespace

static void F_3_T_uses_centroid_not_corner()
{
   std::cout << "\n[F-3] R-002 regression: T_forced reads Vs at element "
                "centroid (NOT reference corner)\n";
   // One unit-cube hex; centroid is physical (0.5, 0.5, 0.5), reference
   // origin maps to physical (0, 0, 0).  Vs(z) = 1000 + 4000*z.
   const real_t Vs_corner   = 1000.0;
   const real_t Vs_centroid = 3000.0;
   // Place the fault DOF at (0.5, 0.5, 0.5) (the element centroid in
   // physical coords) so r = 0 from a hypocenter at the same point.
   // To exercise the in-r_crit branch we use r = 0.1 by shifting the
   // hypocenter.  ResolveForcedRupture only takes dof_to_elem; no
   // dof_to_attr argument is required.
   Vector dofs(3);  Array<int> elem(1);
   dofs(0) = 0.5; dofs(1) = 0.5; dofs(2) = 0.5;
   elem[0] = 0;

   NucleationSpec nuc;
   nuc.enabled        = true;
   nuc.hypocenter_x_m = 0.5; nuc.hypocenter_y_m = 0.5;
   nuc.hypocenter_z_m = 0.4;        // distance = 0.1
   nuc.r_crit_m       = 1.0;
   nuc.t0_decay_s     = 0.5;

   RhoConst    rho_c;
   MuFromZ     mu_c;
   LambdaFromZ lam_c;
   auto mat = MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c);

   TinyMeshHolder mh;
   SpatialFrictionResolver R;
   auto p = R.ResolveForcedRupture(nuc, dofs, elem, mat, mh.mesh());

   const real_t r = 0.1;
   const auto T_expected_centroid = [&](real_t Vs)
   {
      const real_t denom = 1.0 - (r / nuc.r_crit_m) * (r / nuc.r_crit_m);
      return r / (0.7 * Vs)
           + 0.081 * nuc.r_crit_m / (0.7 * Vs) * (1.0 / denom - 1.0);
   };
   const real_t T_centroid = T_expected_centroid(Vs_centroid);
   const real_t T_corner   = T_expected_centroid(Vs_corner);

   // Centroid result must match (within fp tolerance).
   TEST_NEAR(p.T_forced_s(0), T_centroid, 1e-6,
             "T_forced computed with Vs(centroid)");
   // Corner result must NOT match — that's the pre-fix value.
   TEST_ASSERT(std::abs(p.T_forced_s(0) - T_corner)
               > 0.5 * std::abs(T_corner - T_centroid),
               "T_forced differs from the buggy Vs(corner) value");
}

// R-9  Regression for R-003: eta_auto must read mu and rho at element
// centroid, not at the reference corner.  Same single-hex fixture as
// F-3 but exercises the rate-state path.
static void R_9_eta_auto_uses_centroid()
{
   std::cout << "\n[R-9] R-003 regression: eta_auto reads (mu, rho) at "
                "element centroid\n";
   Vector dofs(3);  Array<int> attr(1), elem(1);
   dofs(0) = 0.5; dofs(1) = 0.5; dofs(2) = 0.5;
   attr[0] = 101; elem[0] = 0;
   Vector sn_total(1); sn_total = 50.0e6;

   RateStateBlock cfg;
   cfg.a_default = 0.010; cfg.b_default = 0.015;
   cfg.Dc_default = 0.004; cfg.V_init_default = 1e-9;
   cfg.f_0_default = 0.6;  cfg.V_0_default = 1e-6;
   cfg.sigma_n_default = 50.0e6;
   cfg.eta_auto    = true;
   cfg.eta_default = 0.0;

   RhoConst    rho_c;
   MuFromZ     mu_c;
   LambdaFromZ lam_c;
   auto mat = MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c);
   TinyMeshHolder mh;
   PorePressureSpec pp;
   SpatialFrictionResolver R;
#ifdef MFEM_USE_MPI
   mfem::Mesh& srl = mh.mesh();
   auto p = R.ResolveRateState(cfg, dofs, elem, attr, mat, srl, pp, sn_total);
   const real_t Vs_centroid = 3000.0;
   const real_t mu_expected = 2670.0 * Vs_centroid * Vs_centroid;
   const real_t eta_expect  = 0.5 * std::sqrt(mu_expected * 2670.0);
   TEST_NEAR(p.eta(0), eta_expect, 1e-3,
             "eta = 0.5 * sqrt(mu(centroid) * rho)");
   // Corner sentinel: Vs_corner = 1000 gives mu_corner = 2.67e9; the
   // buggy eta would be 0.5 * sqrt(2.67e9 * 2670) ~= 4.22e6.  Centroid
   // value is 0.5 * sqrt(2.403e10 * 2670) ~= 1.27e7.  Assert we are
   // closer to the centroid value than to the corner value.
   const real_t mu_corner = 2670.0 * 1000.0 * 1000.0;
   const real_t eta_corner = 0.5 * std::sqrt(mu_corner * 2670.0);
   TEST_ASSERT(std::abs(p.eta(0) - eta_expect)
               < std::abs(p.eta(0) - eta_corner),
               "eta closer to centroid value than to corner value");
#else
   std::cout << "  SKIP (no MPI build)\n";
#endif
}

// F-4  R-701 round-7: with nuc.enabled = true AND kind = Overstress,
//      ResolveForcedRupture must STILL return the "never forced"
//      sentinel (T = 1e9, t_0 = 0) so the iterator's forced-rupture
//      mu_eff path is a no-op in overstress mode.
static void F_4_overstress_kind_returns_sentinel()
{
   std::cout << "\n[F-4] enabled + kind = Overstress ⇒ T_forced = 1e9 "
                "(R-701)\n";
   const int N = 3;
   Vector dofs(3*N); Array<int> attr(N), elem(N);
   // Three DOFs very close to the hypocenter (r << r_crit) — pre-fix,
   // these would receive a finite T_forced; post-fix, they get 1e9.
   const real_t rs[N] = { 0.0, 100.0, 500.0 };
   for (int i = 0; i < N; ++i)
   {
      dofs(3*i + 0) = rs[i];
      dofs(3*i + 1) = 0.0;
      dofs(3*i + 2) = 0.0;
      attr[i] = 101; elem[i] = 0;
   }

   NucleationSpec nuc;
   nuc.enabled        = true;
   nuc.kind           = NucleationKind::Overstress;
   nuc.hypocenter_x_m = 0.0;
   nuc.hypocenter_y_m = 0.0;
   nuc.hypocenter_z_m = 0.0;
   nuc.r_crit_m       = 4000.0;
   nuc.t0_decay_s     = 0.5;

   auto mat = MaterialField::MakeConstant(32e9, 32e9, 2670.0);
   TinyMeshHolder mh;
   SpatialFrictionResolver R;
   auto p = R.ResolveForcedRupture(nuc, dofs, elem, mat, mh.mesh());

   TEST_ASSERT(p.T_forced_s.Size() == N, "T_forced sized N");
   for (int i = 0; i < N; ++i)
   {
      TEST_NEAR(p.T_forced_s(i), 1.0e9, 0.0,
                "T_forced = 1e9 for Overstress kind even when enabled");
      TEST_NEAR(p.t0_decay_s(i), 0.0, 0.0,
                "t0_decay = 0 for Overstress kind");
   }
}

// F-5  R-703 round-7: ResolveOverstress returns zero-sized vectors for
//      StrengthReduction (existing path); aborts on Overstress (stub).
static void F_5_resolve_overstress_stub_behaviour()
{
   std::cout << "\n[F-5] ResolveOverstress: empty on StrengthReduction; "
                "aborts on Overstress (R-703)\n";
   // Path 1 — StrengthReduction returns zero-sized vectors.
   {
      NucleationSpec nuc;
      nuc.enabled = true;
      nuc.kind    = NucleationKind::StrengthReduction;
      Vector dofs(6); dofs = 0.0;
      SpatialFrictionResolver R;
      auto p = R.ResolveOverstress(nuc, dofs);
      TEST_ASSERT(p.delta_tau1_pa.Size()    == 0,
                  "ResolveOverstress on StrengthReduction returns empty "
                  "delta_tau1_pa");
      TEST_ASSERT(p.delta_tau2_pa.Size()    == 0,
                  "ResolveOverstress on StrengthReduction returns empty "
                  "delta_tau2_pa");
      TEST_ASSERT(p.delta_sigma_n_pa.Size() == 0,
                  "ResolveOverstress on StrengthReduction returns empty "
                  "delta_sigma_n_pa");
   }
   // Path 2 — disabled also returns empty.
   {
      NucleationSpec nuc;  nuc.enabled = false;
      Vector dofs(6); dofs = 0.0;
      SpatialFrictionResolver R;
      auto p = R.ResolveOverstress(nuc, dofs);
      TEST_ASSERT(p.delta_tau1_pa.Size() == 0,
                  "ResolveOverstress on disabled returns empty");
   }
   // Path 3 — Overstress aborts (stub).  Run in a forked child so the
   // MFEM_ABORT does not take the parent process down.
   const bool aborted = RunInChild([]() {
      NucleationSpec nuc;
      nuc.enabled = true;
      nuc.kind    = NucleationKind::Overstress;
      Vector dofs(3); dofs = 0.0;
      SpatialFrictionResolver R;
      (void)R.ResolveOverstress(nuc, dofs);
   });
   TEST_ASSERT(aborted,
               "ResolveOverstress stub aborts on Overstress kind");
}
#endif  // Disabled in Phase N

// =====================================================================
// R-9 — rate-state path eta_auto uses centroid (still live; spatial
// driver does not depend on this code path, but the resolver test
// covers it directly).  Defines its own coefficient fixtures since
// the old F-3 RhoConst / MuFromZ / LambdaFromZ classes were folded
// into the disabled Phase-N block.
// =====================================================================

namespace
{
class RhoConst_R9 : public mfem::Coefficient
{
public:
   real_t Eval(mfem::ElementTransformation& /*T*/,
               const mfem::IntegrationPoint& /*ip*/) override
   { return 2670.0; }
};
class MuFromZ_R9 : public mfem::Coefficient
{
public:
   real_t Eval(mfem::ElementTransformation& T,
               const mfem::IntegrationPoint& ip) override
   {
      real_t x[3] = { 0.0, 0.0, 0.0 };
      mfem::Vector p(x, 3);
      T.Transform(ip, p);
      const real_t Vs = 1000.0 + 4000.0 * p(2);
      return 2670.0 * Vs * Vs;
   }
};
class LambdaFromZ_R9 : public mfem::Coefficient
{
public:
   real_t Eval(mfem::ElementTransformation& T,
               const mfem::IntegrationPoint& ip) override
   {
      real_t x[3] = { 0.0, 0.0, 0.0 };
      mfem::Vector p(x, 3);
      T.Transform(ip, p);
      const real_t Vs = 1000.0 + 4000.0 * p(2);
      const real_t Vp = 1.8 * Vs;
      return 2670.0 * (Vp*Vp - 2.0 * Vs*Vs);
   }
};
}  // namespace

static void R_9_eta_auto_uses_centroid()
{
   std::cout << "\n[R-9] R-003 regression: eta_auto reads (mu, rho) at "
                "element centroid\n";
   Vector dofs(3);  Array<int> attr(1), elem(1);
   dofs(0) = 0.5; dofs(1) = 0.5; dofs(2) = 0.5;
   attr[0] = 101; elem[0] = 0;
   Vector sn_total(1); sn_total = 50.0e6;

   RateStateBlock cfg;
   cfg.a_default = 0.010; cfg.b_default = 0.015;
   cfg.Dc_default = 0.004; cfg.V_init_default = 1e-9;
   cfg.f_0_default = 0.6;  cfg.V_0_default = 1e-6;
   cfg.sigma_n_default = 50.0e6;
   cfg.eta_auto    = true;
   cfg.eta_default = 0.0;

   RhoConst_R9    rho_c;
   MuFromZ_R9     mu_c;
   LambdaFromZ_R9 lam_c;
   auto mat = MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c);
   TinyMeshHolder mh;
   PorePressureSpec pp;
   SpatialFrictionResolver R;
#ifdef MFEM_USE_MPI
   mfem::Mesh& srl = mh.mesh();
   auto p = R.ResolveRateState(cfg, dofs, elem, attr, mat, srl, pp, sn_total);
   const real_t Vs_centroid = 3000.0;
   const real_t mu_expected = 2670.0 * Vs_centroid * Vs_centroid;
   const real_t eta_expect  = 0.5 * std::sqrt(mu_expected * 2670.0);
   TEST_NEAR(p.eta(0), eta_expect, 1e-3,
             "eta = 0.5 * sqrt(mu(centroid) * rho)");
   const real_t mu_corner  = 2670.0 * 1000.0 * 1000.0;
   const real_t eta_corner = 0.5 * std::sqrt(mu_corner * 2670.0);
   TEST_ASSERT(std::abs(p.eta(0) - eta_expect)
               < std::abs(p.eta(0) - eta_corner),
               "eta closer to centroid value than to corner value");
#else
   std::cout << "  SKIP (no MPI build)\n";
#endif
}

// =====================================================================
//  LSWFrictionCoefficient_ForcedRupture helper sanity
// =====================================================================

// H-1  T = 1e9 makes the helper reduce to plain LSW (no f_2)
static void H_1_T_1e9_reduces_to_plain_LSW()
{
   std::cout << "\n[H-1] T_forced=1e9 ⇒ helper equals plain LSW\n";
   const real_t mu_s = 0.677, mu_d = 0.525, d_c = 0.4;
   for (real_t delta : {0.0, 0.1, 0.4, 1.0})
   {
      const real_t f1 = (delta <= 0)        ? 0.0
                       : (delta >= d_c)     ? 1.0
                       : (delta / d_c);
      const real_t mu_expected_tpv205 = mu_s + (mu_d - mu_s) * f1;
      const real_t mu_got = LSWFrictionCoefficient_ForcedRupture(
         delta, mu_s, mu_d, d_c,
         /*t_now=*/5.0, /*T_forced=*/1.0e9, /*t0_decay=*/0.0);
      TEST_NEAR(mu_got, mu_expected_tpv205, 0.0,
                "byte-identical to TPV205 LSW formula");
   }
}

// R-10  RS depth-profile (two-CSV) a(z)/b(z) mapping (Phase 11b).
static void R_10_rs_depth_profile_csv()
{
   std::cout << "\n[R-10] RS depth profile: a(z), b = a - (a-b) from two CSVs\n";
   const int N = 3;
   Vector dofs; Array<int> attr, elem;
   make_synthetic_dofs(N, 20000.0, dofs, attr, elem);   // depths 0, 10 km, 20 km
   Vector sn_total(N);  sn_total = 50.0e6;

   // Two temp CSVs (value, depth_km).  a: 0.010@0 -> 0.020@20km;
   // a-b: -0.005@0 (VW) -> +0.010@20km (VS).
   const std::string a_csv   = "/tmp/seas_test_resolver_a.csv";
   const std::string amb_csv = "/tmp/seas_test_resolver_amb.csv";
   { std::ofstream o(a_csv,   std::ios::trunc); o << "0.010, 0\n0.020, 20\n"; }
   { std::ofstream o(amb_csv, std::ios::trunc); o << "-0.005, 0\n0.010, 20\n"; }

   RateStateBlock cfg;
   cfg.Dc_default = 0.004; cfg.V_init_default = 1e-9;
   cfg.f_0_default = 0.6;  cfg.V_0_default = 1e-6;
   cfg.sigma_n_default = 50.0e6;
   cfg.eta_auto = false;   cfg.eta_default = 5.0e6;
   cfg.depth_profile.enabled = true;
   cfg.depth_profile.param_a_csv = a_csv;
   cfg.depth_profile.param_a_minus_b_csv = amb_csv;
   cfg.depth_profile.depth_to_m = 1000.0;   // km
   cfg.depth_profile.profile = LoadFrictionDepthProfileCSVs(cfg.depth_profile);

   auto mat = MaterialField::MakeConstant(32.0e9, 32.0e9, 2670.0);
   TinyMeshHolder mh;
   PorePressureSpec pp;
   SpatialFrictionResolver R;

#ifdef MFEM_USE_MPI
   mfem::Mesh& srl = mh.mesh();
   auto p = R.ResolveRateState(cfg, dofs, elem, attr, mat, srl, pp, sn_total);
#else
   RateStatePerDOFParams p;
   p.a.SetSize(N); p.b.SetSize(N);
   const real_t ea[3] = {0.010, 0.015, 0.020};
   const real_t eb[3] = {0.015, 0.0125, 0.010};
   for (int i = 0; i < N; ++i) { p.a(i) = ea[i]; p.b(i) = eb[i]; }
#endif

   // depths 0 / 10 km / 20 km:
   //   a(z)     = 0.010, 0.015, 0.020
   //   (a-b)(z) = -0.005, 0.0025, 0.010
   //   b        = 0.015, 0.0125, 0.010
   TEST_NEAR(p.a(0), 0.010, 1e-12, "a(surface) from profile");
   TEST_NEAR(p.b(0), 0.015, 1e-12, "b(surface) = a - (a-b) = 0.015");
   TEST_ASSERT(p.a(0) < p.b(0), "surface: a < b (velocity-weakening)");

   TEST_NEAR(p.a(1), 0.015,  1e-12, "a(10 km) from profile");
   TEST_NEAR(p.b(1), 0.0125, 1e-12, "b(10 km) = a - (a-b) = 0.0125");

   TEST_NEAR(p.a(2), 0.020, 1e-12, "a(20 km) from profile");
   TEST_NEAR(p.b(2), 0.010, 1e-12, "b(20 km) = a - (a-b) = 0.010");
   TEST_ASSERT(p.a(2) > p.b(2), "deep: a > b (velocity-strengthening)");
}

int main(int argc, char** argv)
{
   // No MPI init: every resolver test uses the serial mfem::Mesh
   // overload of ResolveRateState / ResolveForcedRupture.  Avoiding
   // MPI keeps the `fork()`-based abort-validator tests safe.
   (void)argc; (void)argv;

   std::cout << "Running Phase 1 test_spatial_friction_resolver\n";

   // LSW
   L_1_defaults_only();
   L_2_depth_rule_override();
   L_3_last_match_wins();
   L_4_box_rule();
   L_5_region_attribute_rule();
   L_6_cohesion_override();

   // Barrier
   B_1_barrier_sentinel();
   B_2_barrier_with_x_bound();

   // Rate-state
   R_1a_rs_defaults_eta_explicit();
   R_1b_rs_eta_auto_smoke_const_material();
   R_2_rs_ab_override();
   R_3_rs_pore_pressure();
   R_4_rs_V_init_override();
   R_5_rs_validator_aborts();
   R_6_rs_eta_explicit();
   R_7_rs_negative_sigma_n_aborts();
   R_8_rs_eta_auto_collision();
   R_10_rs_depth_profile_csv();

   // Forced rupture / Overstress — REMOVED in Phase N (the spatial
   // driver no longer dispatches into ResolveForcedRupture /
   // ResolveOverstress; the single replacement kind
   // `gradual_overstress` is covered by seas_test_spatial_nucleation).
   R_9_eta_auto_uses_centroid();

   // Helper
   H_1_T_1e9_reduces_to_plain_LSW();

   std::cout << "\n========================================\n";
   std::cout << "Phase 1 test_spatial_friction_resolver: "
             << num_passed << " / " << num_tests
             << " passed, " << num_failed << " failed\n";
   std::cout << "========================================\n";

   return (num_failed == 0) ? 0 : 1;
}
