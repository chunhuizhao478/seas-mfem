// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_planar_tpv_basis.cpp — Phase 8 / Decision D1
// (PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24.md): for the
// planar y=0 TPV fault with ref_normal=(0,-1,0), up=(0,0,1), the
// driver-local FaultBasis yields the canonical fault-local frame
//   dip    (tangent1) = (0, 0, -1)
//   strike (tangent2) = (+1, 0, 0)
// This is the "test if planar TPV is correct" the user asked for and the
// frame CLAUDE.md documents (can_t1 = down-dip, can_t2 = along-strike).
//
// Two layers:
//   (A) Pure-function test of FaultBasis::ComputeOrientedFrame — pins the
//       EXACT canonical frame for a ref-aligned normal, and documents the
//       intentional v61 sign-flip convention for an anti-aligned normal
//       (the whole frame negates, consistent with the DG face-normal flip
//       between elements/ranks — see fault_basis.hpp).
//   (B) Mesh-backed test: build the basis on a tiny inline mesh with a
//       single y-normal interior face and assert the per-face / per-QP
//       frame is canonical up to that documented global sign.
//
// Framework: the project's custom TEST_ASSERT / TEST_NEAR harness.

#include "mfem.hpp"

#include "../../fault/fault_basis.hpp"

#include <cmath>
#include <iostream>

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
constexpr real_t kTol = 1e-12;

real_t Dot3(const real_t* a, const real_t* b)
{ return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }

// Assert {dip, strike, normal} is an orthonormal frame with normal along
// the y-axis, strike along x, dip along z (the planar-y-fault geometry),
// and equal to the canonical convention UP TO the documented global sign:
//   either (dip=(0,0,-1), strike=(+1,0,0))   [sign_flipped == false]
//   or     (dip=(0,0,+1), strike=(-1,0,0))   [sign_flipped == true].
void AssertCanonicalUpToSign(const real_t* dip, const real_t* strike,
                             const real_t* normal, const char* tag)
{
   // Orthonormal.
   TEST_NEAR(Dot3(dip, dip),       1.0, 1e-10, std::string(tag) + ": |dip| == 1");
   TEST_NEAR(Dot3(strike, strike), 1.0, 1e-10, std::string(tag) + ": |strike| == 1");
   TEST_NEAR(Dot3(normal, normal), 1.0, 1e-10, std::string(tag) + ": |normal| == 1");
   TEST_NEAR(Dot3(dip, strike),    0.0, 1e-10, std::string(tag) + ": dip . strike == 0");
   TEST_NEAR(Dot3(dip, normal),    0.0, 1e-10, std::string(tag) + ": dip . normal == 0");
   TEST_NEAR(Dot3(strike, normal), 0.0, 1e-10, std::string(tag) + ": strike . normal == 0");
   // Axis alignment: strike on x, dip on z, normal on y.
   TEST_NEAR(std::abs(strike[0]), 1.0, 1e-10, std::string(tag) + ": strike along x");
   TEST_NEAR(std::abs(dip[2]),    1.0, 1e-10, std::string(tag) + ": dip along z");
   TEST_NEAR(std::abs(normal[1]), 1.0, 1e-10, std::string(tag) + ": normal along y");
   // Canonical up to the single global sign: strike_x and dip_z share the
   // relation strike=(+1,0,0)<->dip=(0,0,-1) (and the negation of both).
   TEST_NEAR(strike[0], -dip[2], 1e-10,
             std::string(tag) + ": strike_x == -dip_z (canonical handedness)");
}
} // anon

// ---------------------------------------------------------------------
// (A) Pure-function: exact canonical frame for the ref-aligned normal.
// ---------------------------------------------------------------------
static void T_oriented_frame_ref_aligned()
{
   std::cout << "\n[A1] ComputeOrientedFrame, ref-aligned normal (0,-1,0)\n";
   Vector ref_normal(3); ref_normal = 0.0; ref_normal(1) = -1.0;
   Vector up(3);         up = 0.0;         up(2) = 1.0;
   Vector n_raw(3);      n_raw = 0.0;      n_raw(1) = -1.0;  // == ref_normal

   real_t normal[3], dip[3], strike[3], nl; bool flipped;
   FaultBasis::ComputeOrientedFrame(n_raw, 3, ref_normal, up,
                                    normal, dip, strike, flipped, nl);

   TEST_ASSERT(!flipped, "ref-aligned: sign_flipped == false");
   // EXACT canonical frame (the convention CLAUDE.md / D1 document).
   TEST_NEAR(dip[0],    0.0, kTol, "dip == (0,0,-1): x");
   TEST_NEAR(dip[1],    0.0, kTol, "dip == (0,0,-1): y");
   TEST_NEAR(dip[2],   -1.0, kTol, "dip == (0,0,-1): z");
   TEST_NEAR(strike[0], 1.0, kTol, "strike == (+1,0,0): x");
   TEST_NEAR(strike[1], 0.0, kTol, "strike == (+1,0,0): y");
   TEST_NEAR(strike[2], 0.0, kTol, "strike == (+1,0,0): z");
   TEST_NEAR(normal[0], 0.0, kTol, "normal == (0,-1,0): x");
   TEST_NEAR(normal[1],-1.0, kTol, "normal == (0,-1,0): y");
   TEST_NEAR(normal[2], 0.0, kTol, "normal == (0,-1,0): z");
}

// ---------------------------------------------------------------------
// (A) Pure-function: anti-aligned normal negates the whole frame (v61).
// ---------------------------------------------------------------------
static void T_oriented_frame_anti_aligned()
{
   std::cout << "\n[A2] ComputeOrientedFrame, anti-aligned normal (0,+1,0)\n";
   Vector ref_normal(3); ref_normal = 0.0; ref_normal(1) = -1.0;
   Vector up(3);         up = 0.0;         up(2) = 1.0;
   Vector n_raw(3);      n_raw = 0.0;      n_raw(1) = 1.0;   // opposite ref

   real_t normal[3], dip[3], strike[3], nl; bool flipped;
   FaultBasis::ComputeOrientedFrame(n_raw, 3, ref_normal, up,
                                    normal, dip, strike, flipped, nl);

   TEST_ASSERT(flipped, "anti-aligned: sign_flipped == true");
   // The entire frame is the global-negated canonical (documented, correct:
   // the DG face normal also flips, so the embedded slip flips with it).
   TEST_NEAR(dip[2],    1.0, kTol, "anti-aligned dip == (0,0,+1)");
   TEST_NEAR(strike[0],-1.0, kTol, "anti-aligned strike == (-1,0,0)");
   TEST_NEAR(normal[1], 1.0, kTol, "anti-aligned normal == (0,+1,0)");
   AssertCanonicalUpToSign(dip, strike, normal, "A2");
}

// ---------------------------------------------------------------------
// (B) Mesh-backed: tiny inline mesh, single y-normal interior face.
// ---------------------------------------------------------------------
static void T_mesh_backed_basis()
{
   std::cout << "\n[B] FaultBasis on a tiny inline mesh (y-normal interior face)\n";
   // 1 x 2 x 1 hex box -> exactly one interior face, normal along y.
   Mesh mesh = Mesh::MakeCartesian3D(1, 2, 1, Element::HEXAHEDRON,
                                     2.0, 2.0, 2.0);

   Array<int> fault_faces;
   for (int f = 0; f < mesh.GetNumFaces(); ++f)
   {
      if (mesh.GetInteriorFaceTransformations(f) != nullptr)
      { fault_faces.Append(f); }
   }
   TEST_ASSERT(fault_faces.Size() == 1,
               "tiny mesh has exactly one interior (fault) face");
   if (fault_faces.Size() == 0) { return; }

   Vector ref_normal(3); ref_normal = 0.0; ref_normal(1) = -1.0;
   Vector up(3);         up = 0.0;         up(2) = 1.0;

   FaultBasis fb;
   fb.Compute(mesh, fault_faces, ref_normal, up);
   TEST_ASSERT(fb.NumFaces() == 1, "FaultBasis::NumFaces() == 1");

   const FaultBasisData& b = fb.GetBasis(0);
   AssertCanonicalUpToSign(b.tangent1, b.tangent2, b.normal, "B-centroid");

   // Per-QP basis: every QP must carry the same canonical-up-to-sign frame.
   auto* ftr = mesh.GetInteriorFaceTransformations(fault_faces[0]);
   const Geometry::Type face_geom = ftr->GetGeometryType();
   const IntegrationRule& ir = IntRules.Get(face_geom, 2);  // order-2 face rule
   fb.ComputeQPBasis(mesh, fault_faces, ref_normal, up, ir);
   const FaultBasisData& bq = fb.GetBasis(0);
   TEST_ASSERT(!bq.qp_data.empty(), "per-QP basis populated");
   for (size_t q = 0; q < bq.qp_data.size(); ++q)
   {
      AssertCanonicalUpToSign(bq.qp_data[q].tangent1,
                              bq.qp_data[q].tangent2,
                              bq.qp_data[q].normal, "B-qp");
   }
}

int main(int, char**)
{
   std::cout << "Running Phase 8 / D1 planar-TPV basis tests\n";
   T_oriented_frame_ref_aligned();
   T_oriented_frame_anti_aligned();
   T_mesh_backed_basis();

   std::cout << "\n========================================\n";
   std::cout << "Phase 8 planar-basis: " << num_passed << " / " << num_tests
             << " passed, " << num_failed << " failed\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
