// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_spatial_constant_stress_source.cpp — Phase 3b of
// spatial_dynamic_rupture_plan.md (rev-3).
//
// Plan §Phase 3b §Acceptance: 5 tests including a planar n·σ·n
// exactness check and a zero-stress sanity check.  The 64-DOF
// curvilinear cross-check vs the Python reference
// `verify_constant_tensor_projection.py` runs as an offline acceptance
// test in Phase 4 (it needs the real SAFS mesh + FaultGeometry which
// the unit-test stage does not have).

#include "mfem.hpp"

#include "../../spatial/code/spatial_stress.hpp"

#include <array>
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

namespace
{

// Resolve sigma_n = n . S . n analytically for a 3x3 symmetric tensor.
real_t analytic_sigma_n(const std::array<real_t, 6>& s,
                        const real_t n[3])
{
   // s = (sxx, syy, szz, sxy, syz, sxz)
   const real_t Sn0 = s[0]*n[0] + s[3]*n[1] + s[5]*n[2];
   const real_t Sn1 = s[3]*n[0] + s[1]*n[1] + s[4]*n[2];
   const real_t Sn2 = s[5]*n[0] + s[4]*n[1] + s[2]*n[2];
   return n[0]*Sn0 + n[1]*Sn1 + n[2]*Sn2;
}

}  // namespace

// C-1  Evaluate(x, y, z) returns the same tensor at every (x, y, z).
static void C_1_evaluate_constant()
{
   std::cout << "\n[C-1] Evaluate(x, y, z) returns the constant tensor everywhere\n";
   ConstantTensorStressSource src(100e6, 60e6, 80e6,
                                  20e6,  0,    5e6);
   const auto S1 = src.Evaluate(0,   0,   0);
   const auto S2 = src.Evaluate(1e5, 1e6, -3000);
   const auto S3 = src.Evaluate(-2e6, 5e5, 1e3);
   TEST_NEAR(S1(0, 0), 100e6, 0.0, "S(0,0,0)[xx]");
   for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
      {
         TEST_NEAR(S2(i, j), S1(i, j), 0.0, "S2 == S1 byte");
         TEST_NEAR(S3(i, j), S1(i, j), 0.0, "S3 == S1 byte");
      }
}

// C-2  Symmetry: ctor takes six components and produces a symmetric
//      DenseMatrix.
static void C_2_symmetric()
{
   std::cout << "\n[C-2] Source is symmetric by construction\n";
   ConstantTensorStressSource src(1e6, 2e6, 3e6, 4e6, 5e6, 6e6);
   const auto& S = src.Tensor();
   TEST_NEAR(S(0, 1), S(1, 0), 0.0, "S_xy == S_yx");
   TEST_NEAR(S(1, 2), S(2, 1), 0.0, "S_yz == S_zy");
   TEST_NEAR(S(0, 2), S(2, 0), 0.0, "S_xz == S_zx");
   TEST_NEAR(S(0, 0), 1e6, 0.0, "S_xx");
   TEST_NEAR(S(1, 1), 2e6, 0.0, "S_yy");
   TEST_NEAR(S(2, 2), 3e6, 0.0, "S_zz");
   TEST_NEAR(S(0, 1), 4e6, 0.0, "S_xy");
   TEST_NEAR(S(1, 2), 5e6, 0.0, "S_yz");
   TEST_NEAR(S(0, 2), 6e6, 0.0, "S_xz");
}

// C-3  BBox is "infinite" and ContainsBBox returns true for any input.
static void C_3_bbox_infinite()
{
   std::cout << "\n[C-3] BBox is infinite; ContainsBBox always true\n";
   ConstantTensorStressSource src(1, 2, 3, 0, 0, 0);
   const auto& bb = src.BBox();
   const real_t inf = std::numeric_limits<real_t>::infinity();
   TEST_ASSERT(bb[0] == -inf && bb[1] == inf, "x bbox infinite");
   TEST_ASSERT(bb[2] == -inf && bb[3] == inf, "y bbox infinite");
   TEST_ASSERT(bb[4] == -inf && bb[5] == inf, "z bbox infinite");
   TEST_ASSERT(src.ContainsBBox(-1e9, 1e9, -1e9, 1e9, -1e9, 1e9),
               "huge bbox contained");
   TEST_ASSERT(src.ContainsBBox(0, 0, 0, 0, 0, 0),
               "trivial bbox contained");
}

// C-4  Analytic n·σ·n match on a small synthetic set of fault normals.
//      This is the L∞ ≤ 1e-10 planar exactness gate; the integration
//      with FaultGeometry happens in Phase 4 acceptance.
static void C_4_analytic_n_sigma_n()
{
   std::cout << "\n[C-4] Analytic n . S . n on three test normals\n";
   const real_t sxx = 100e6, syy = 60e6, szz = 80e6;
   const real_t sxy = 20e6,  syz =  0,   sxz = 5e6;
   ConstantTensorStressSource src(sxx, syy, szz, sxy, syz, sxz);

   const std::array<std::array<real_t, 3>, 3> normals = {{
      {{ 1.0, 0.0, 0.0 }},      // E-normal fault
      {{ 0.0, 1.0, 0.0 }},      // N-normal fault
      {{ 1.0/std::sqrt(2.0), 1.0/std::sqrt(2.0), 0.0 }}, // NE-normal
   }};

   for (const auto& n : normals)
   {
      const auto S = src.Evaluate(0, 0, 0);
      real_t Sn[3];
      for (int r = 0; r < 3; ++r)
      {
         Sn[r] = S(r, 0) * n[0] + S(r, 1) * n[1] + S(r, 2) * n[2];
      }
      const real_t sn_got = n[0]*Sn[0] + n[1]*Sn[1] + n[2]*Sn[2];
      const std::array<real_t, 6> s6 = { sxx, syy, szz, sxy, syz, sxz };
      const real_t sn_expect = analytic_sigma_n(s6, n.data());
      const real_t tol = std::abs(sn_expect) * 1e-10 + 1e-6;
      TEST_NEAR(sn_got, sn_expect, tol,
                std::string("n · S · n analytic match"));
   }
}

// C-5  Zero-stress sanity: all six scalars = 0 produces a zero tensor.
static void C_5_zero_stress()
{
   std::cout << "\n[C-5] Zero scalars -> zero tensor\n";
   ConstantTensorStressSource src(0, 0, 0, 0, 0, 0);
   const auto S = src.Evaluate(123.4, -56.7, -1000.0);
   for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
      {
         TEST_NEAR(S(i, j), 0.0, 0.0, "zero tensor");
      }
}

int main(int, char**)
{
   std::cout << "Running Phase 3b test_spatial_constant_stress_source\n";
   C_1_evaluate_constant();
   C_2_symmetric();
   C_3_bbox_infinite();
   C_4_analytic_n_sigma_n();
   C_5_zero_stress();
   std::cout << "\n========================================\n";
   std::cout << "Phase 3b test_spatial_constant_stress_source: "
             << num_passed << " / " << num_tests
             << " passed, " << num_failed << " failed\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
