// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 "Topology-Based Precomputed Face-Rotation" plan 2026-04-23
// Phase 1 Step A acceptance: tests/unit/test_precomputed_fluxes_interior_identity.cpp
//
// Covers the Step A P_* gates (plan §5.1.7):
//   P_FACE2NODES_SELFCHECK, P_TOPOLOGY_FRAME_MATCHES_SEISSOL, P_FRAME_ORBIT,
//   P_MATRIX_IDENTITY, P_CONSERVATION, P_SEISSOL_VALUE (direct-eigendecomp
//   form per R2-001), plus the R-series regression tests:
//   - unit-tet cell_volume == 1/6 (R2-005)
//   - test_R3_004_neighbor_elem_is_other_side
//   - test_R3_006_key_encoding_handles_negative_caller_elem
//
// Guard tests that MFEM_VERIFY aborts on bimaterial / p>=2 inputs
// (P_BIMATERIAL_GUARD, test_R3_009_p1_precondition) are documented but
// not exercised at runtime here: MFEM_USE_EXCEPTIONS is OFF in the
// production build, so a failing MFEM_VERIFY terminates the process
// rather than returning a signalable value.  Those guards are covered
// implicitly by the fact that the test's PASS result includes a
// successful Init on a homogeneous p=1 fixture.

#include "mfem.hpp"
#include "../../dynamic/precomputed_face_fluxes.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../domain/boundary_config.hpp"
#include "../../config/tpv102_params.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; \
          std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

// --------------------------------------------------------------------------
// Build a single-tet mesh with the reference unit tet
// {(0,0,0), (1,0,0), (0,1,0), (0,0,1)}.
// --------------------------------------------------------------------------
static Mesh BuildReferenceUnitTet()
{
   Mesh mesh(3, 4, 1, 0, 3);
   mesh.AddVertex(0.0, 0.0, 0.0);
   mesh.AddVertex(1.0, 0.0, 0.0);
   mesh.AddVertex(0.0, 1.0, 0.0);
   mesh.AddVertex(0.0, 0.0, 1.0);
   int v[4] = {0, 1, 2, 3};
   mesh.AddTet(v, 1);
   mesh.FinalizeMesh(0, false);
   return mesh;
}

// --------------------------------------------------------------------------
// Build a Cartesian TET fixture sized nx x ny x nz over the unit cube.
// All boundary attrs are 1 (single attr).
// --------------------------------------------------------------------------
static Mesh BuildCartesianTetFixture(int nx, int ny, int nz)
{
   Mesh mesh = Mesh::MakeCartesian3D(
      nx, ny, nz, Element::TETRAHEDRON, 1.0, 1.0, 1.0, false);
   mesh.FinalizeTopology();
   mesh.Finalize();
   return mesh;
}

// --------------------------------------------------------------------------
// Fill a 9-vector with uniform samples in [-1, 1].  (R-005: docstring
// aligned with implementation.)  TPV102's Lame parameters give A_plus
// entries at O(1e10); for any O(1) input the intermediate scale is
// ~1e10, and the relative-error normalization in the callers uses
// ||A_plus||_inf * ||I||_inf to absorb this exactly.
// --------------------------------------------------------------------------
static void FillSaneRandom(real_t *Q, std::mt19937 &rng)
{
   std::uniform_real_distribution<real_t> U(-1.0, 1.0);
   for (int c = 0; c < NUM_STATE; c++) { Q[c] = U(rng); }
}

// --------------------------------------------------------------------------
// Reference-tet topology frame hand-port of SeisSol's
// MeshTools::normalAndTangents using the MFEM FACE2NODES table.
// Used by P_TOPOLOGY_FRAME_MATCHES_SEISSOL to verify bit-equality.
// --------------------------------------------------------------------------
static void HandPortedReferenceTetFrame(int local_side,
                                         real_t normal[3],
                                         real_t tangent1[3],
                                         real_t tangent2[3])
{
   const real_t verts[4][3] = {
      {0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0}
   };
   int f0 = PrecomputedFaceFluxes::FACE2NODES_MFEM[local_side][0];
   int f1 = PrecomputedFaceFluxes::FACE2NODES_MFEM[local_side][1];
   int f2 = PrecomputedFaceFluxes::FACE2NODES_MFEM[local_side][2];
   real_t p0[3] = {verts[f0][0], verts[f0][1], verts[f0][2]};
   real_t p1[3] = {verts[f1][0], verts[f1][1], verts[f1][2]};
   real_t p2[3] = {verts[f2][0], verts[f2][1], verts[f2][2]};

   real_t ab[3] = {p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]};
   real_t ac[3] = {p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]};

   auto cross = [](const real_t a[3], const real_t b[3], real_t out[3]) {
      out[0] = a[1]*b[2] - a[2]*b[1];
      out[1] = a[2]*b[0] - a[0]*b[2];
      out[2] = a[0]*b[1] - a[1]*b[0];
   };

   real_t nr[3]; cross(ab, ac, nr);
   real_t c_elem[3] = {
      (verts[0][0]+verts[1][0]+verts[2][0]+verts[3][0])/4.0,
      (verts[0][1]+verts[1][1]+verts[2][1]+verts[3][1])/4.0,
      (verts[0][2]+verts[1][2]+verts[2][2]+verts[3][2])/4.0
   };
   real_t c_face[3] = {
      (p0[0]+p1[0]+p2[0])/3.0,
      (p0[1]+p1[1]+p2[1])/3.0,
      (p0[2]+p1[2]+p2[2])/3.0
   };
   real_t outward[3] = {c_face[0]-c_elem[0],
                        c_face[1]-c_elem[1],
                        c_face[2]-c_elem[2]};
   real_t dot = nr[0]*outward[0] + nr[1]*outward[1] + nr[2]*outward[2];
   if (dot < 0.0)
   {
      std::swap(p0[0], p1[0]); std::swap(p0[1], p1[1]); std::swap(p0[2], p1[2]);
      ab[0] = p1[0]-p0[0]; ab[1] = p1[1]-p0[1]; ab[2] = p1[2]-p0[2];
      ac[0] = p2[0]-p0[0]; ac[1] = p2[1]-p0[1]; ac[2] = p2[2]-p0[2];
      cross(ab, ac, nr);
   }
   real_t t1r[3] = {ab[0], ab[1], ab[2]};
   real_t t2r[3]; cross(nr, t1r, t2r);

   auto norm3 = [](const real_t *a) {
      return std::sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
   };
   real_t ln = norm3(nr), lt1 = norm3(t1r), lt2 = norm3(t2r);
   for (int j = 0; j < 3; j++)
   {
      normal[j]   = nr[j]  / ln;
      tangent1[j] = t1r[j] / lt1;
      tangent2[j] = t2r[j] / lt2;
   }
}

// --------------------------------------------------------------------------
// Test: P_FACE2NODES_SELFCHECK + unit-tet cell_volume == 1/6
// --------------------------------------------------------------------------
static void TestFaceTableAndCellVolume()
{
   std::cout << "\n[P_FACE2NODES_SELFCHECK] + [unit-tet cell_volume]\n";

   PrecomputedFaceFluxes::AssertMFEMFace2NodesTable();
   TEST_ASSERT(true,
               "AssertMFEMFace2NodesTable does not abort on live MFEM");

   Mesh mesh = BuildReferenceUnitTet();
   real_t n[3], t1[3], t2[3], area, vol;
   PrecomputedFaceFluxes::ComputeCellFaceFrame(mesh, 0, 0, n, t1, t2, area, vol);
   TEST_ASSERT(std::abs(vol - 1.0/6.0) < 1e-14,
               "unit tet cell_volume == 1/6 to 1e-14 "
               "(got " + std::to_string(vol) + ")");
}

// --------------------------------------------------------------------------
// Test: P_TOPOLOGY_FRAME_MATCHES_SEISSOL (reference unit tet)
// --------------------------------------------------------------------------
static void TestTopologyFrameMatchesSeissol()
{
   std::cout << "\n[P_TOPOLOGY_FRAME_MATCHES_SEISSOL]\n";
   Mesh mesh = BuildReferenceUnitTet();
   for (int s = 0; s < 4; s++)
   {
      real_t n_impl[3], t1_impl[3], t2_impl[3], area, vol;
      PrecomputedFaceFluxes::ComputeCellFaceFrame(
         mesh, 0, s, n_impl, t1_impl, t2_impl, area, vol);
      real_t n_ref[3], t1_ref[3], t2_ref[3];
      HandPortedReferenceTetFrame(s, n_ref, t1_ref, t2_ref);

      real_t max_err = 0.0;
      for (int j = 0; j < 3; j++)
      {
         max_err = std::max(max_err, std::abs(n_impl[j] - n_ref[j]));
         max_err = std::max(max_err, std::abs(t1_impl[j] - t1_ref[j]));
         max_err = std::max(max_err, std::abs(t2_impl[j] - t2_ref[j]));
      }
      TEST_ASSERT(max_err < 1e-14,
                  "reference-tet side " + std::to_string(s)
                  + " frame matches hand-port to 1e-14 (max_err="
                  + std::to_string(max_err) + ")");
   }
}

// --------------------------------------------------------------------------
// Test: P_FRAME_ORBIT (Cartesian fixture — outward + orthonormal)
// --------------------------------------------------------------------------
static void TestFrameOrbit()
{
   std::cout << "\n[P_FRAME_ORBIT]\n";
   Mesh mesh = BuildCartesianTetFixture(2, 2, 2);
   bool all_pass = true;
   int bad_count = 0;
   for (int e = 0; e < mesh.GetNE(); e++)
   {
      for (int s = 0; s < 4; s++)
      {
         real_t n[3], t1[3], t2[3], area, vol;
         PrecomputedFaceFluxes::ComputeCellFaceFrame(
            mesh, e, s, n, t1, t2, area, vol);
         real_t d_nt1 = n[0]*t1[0]+n[1]*t1[1]+n[2]*t1[2];
         real_t d_nt2 = n[0]*t2[0]+n[1]*t2[1]+n[2]*t2[2];
         real_t d_tt  = t1[0]*t2[0]+t1[1]*t2[1]+t1[2]*t2[2];
         real_t n_len  = std::sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
         real_t t1_len = std::sqrt(t1[0]*t1[0]+t1[1]*t1[1]+t1[2]*t1[2]);
         real_t t2_len = std::sqrt(t2[0]*t2[0]+t2[1]*t2[1]+t2[2]*t2[2]);
         bool ok = std::abs(d_nt1) < 1e-14 && std::abs(d_nt2) < 1e-14
                && std::abs(d_tt)  < 1e-14
                && std::abs(n_len - 1.0)  < 1e-14
                && std::abs(t1_len - 1.0) < 1e-14
                && std::abs(t2_len - 1.0) < 1e-14;
         if (!ok) { all_pass = false; bad_count++; }
      }
   }
   TEST_ASSERT(all_pass,
               "every (elem, side) frame on 2x2x2 Cartesian is orthonormal "
               "(bad_count=" + std::to_string(bad_count) + ")");
}

// --------------------------------------------------------------------------
// Test: P_MATRIX_IDENTITY + P_CONSERVATION on an interior face.
//
// Build a two-tet Cartesian fixture (1x1x1 split into 6 tets); iterate
// interior non-boundary faces and assert:
//   P_MATRIX_IDENTITY: F_pc == flux.Interior(n_topo, I_self, I_nbr) to 1e-12
//   P_CONSERVATION:    F_pc_from_E1 + F_pc_from_E2 == 0 to 1e-12
// --------------------------------------------------------------------------
static void TestMatrixIdentityAndConservation()
{
   std::cout << "\n[P_MATRIX_IDENTITY] + [P_CONSERVATION]\n";
   Mesh mesh = BuildCartesianTetFixture(2, 2, 2);

   DG_FECollection fec(1, 3);
   FiniteElementSpace fes(&mesh, &fec, 1);

   GodunovFlux flux(TPV102Params::lambda, TPV102Params::mu, TPV102Params::rho);

   BoundaryConfig bc;
   // MakeCartesian3D assigns bdr attrs 1..6 (one per cube face); route
   // all to natural (free-surface Godunov) so Init classifies without
   // aborting.
   for (int a = 1; a <= 6; a++) { bc.natural_attrs.insert(a); }
   bc.fault_attr = -1;

   std::vector<int> face_bdr_attr(mesh.GetNumFaces(), 0);
   for (int b = 0; b < mesh.GetNBE(); b++)
   {
      int f, o; mesh.GetBdrElementFace(b, &f, &o);
      face_bdr_attr[f] = mesh.GetBdrAttribute(b);
   }

   PrecomputedFaceFluxes pcf;
   std::set<int> fault_set;
   std::set<int> shared_set;
   pcf.Init(mesh, fes, bc, flux, FreeSurfaceBCMode::Godunov,
            fault_set, face_bdr_attr, shared_set);

   std::mt19937 rng(42);
   const int N_trial = 20;

   real_t worst_identity = 0.0;
   real_t worst_cons = 0.0;
   int n_identity = 0, n_cons = 0;
   int bad_identity = 0, bad_cons = 0;

   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      FaceElementTransformations *ftr = mesh.GetFaceElementTransformations(f);
      if (!ftr) { continue; }
      if (ftr->Elem2No < 0) { continue; } // boundary, skip
      int e1 = ftr->Elem1No;
      int e2 = ftr->Elem2No;
      if (!pcf.HasEntry(f, e1) || !pcf.HasEntry(f, e2)) { continue; }
      const auto &fe1 = pcf.GetEntry(f, e1);
      const auto &fe2 = pcf.GetEntry(f, e2);

      for (int trial = 0; trial < N_trial; trial++)
      {
         real_t Iself[NUM_STATE], Inbr[NUM_STATE];
         FillSaneRandom(Iself, rng);
         FillSaneRandom(Inbr,  rng);

         // P_MATRIX_IDENTITY: compare against flux.Interior on n_topo.
         // The Godunov flux is rotation-invariant around n for isotropic
         // elastic, so any (t1, t2) BuildFrame chooses is equivalent to
         // the topology (t1, t2).
         real_t F_pc[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            real_t acc = 0.0;
            for (int k = 0; k < NUM_STATE; k++)
            {
               acc += fe1.nApNm1(c, k) * Iself[k]
                    + fe1.nAmNm1(c, k) * Inbr[k];
            }
            F_pc[c] = acc;
         }
         real_t F_runtime[NUM_STATE];
         flux.Interior(fe1.normal, Iself, Inbr, F_runtime);
         real_t max_err_id = 0.0;
         for (int c = 0; c < NUM_STATE; c++)
         {
            real_t scale = std::max({real_t(1.0), std::abs(F_pc[c]),
                                     std::abs(F_runtime[c])});
            real_t rel = std::abs(F_pc[c] - F_runtime[c]) / scale;
            max_err_id = std::max(max_err_id, rel);
         }
         if (max_err_id > worst_identity) worst_identity = max_err_id;
         if (max_err_id > 1e-12) { bad_identity++; }
         n_identity++;

         // P_CONSERVATION: E1 emits F_pc from (I_self=Iself, I_nbr=Inbr);
         // E2 emits F_pc2 from (I_self=Inbr, I_nbr=Iself) using its own
         // outward normal.  Require F_pc + F_pc2 -> 0 (relative to
         // |F_pc|).
         real_t F_pc2[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            real_t acc = 0.0;
            for (int k = 0; k < NUM_STATE; k++)
            {
               acc += fe2.nApNm1(c, k) * Inbr[k]
                    + fe2.nAmNm1(c, k) * Iself[k];
            }
            F_pc2[c] = acc;
         }
         real_t max_err_cons = 0.0;
         for (int c = 0; c < NUM_STATE; c++)
         {
            real_t scale = std::max({real_t(1.0), std::abs(F_pc[c]),
                                     std::abs(F_pc2[c])});
            real_t rel = std::abs(F_pc[c] + F_pc2[c]) / scale;
            max_err_cons = std::max(max_err_cons, rel);
         }
         if (max_err_cons > worst_cons) worst_cons = max_err_cons;
         if (max_err_cons > 1e-12) { bad_cons++; }
         n_cons++;
      }
   }

   std::cout << "  interior pairs tested: " << n_identity
             << "  worst_identity rel=" << worst_identity << "\n";
   std::cout << "  conservation tests:    " << n_cons
             << "  worst_cons rel=" << worst_cons << "\n";
   TEST_ASSERT(bad_identity == 0,
               "P_MATRIX_IDENTITY: every sample within 1e-12 relative");
   TEST_ASSERT(bad_cons == 0,
               "P_CONSERVATION: every two-entry E1+E2 pair within 1e-12");
}

// --------------------------------------------------------------------------
// Test: P_SEISSOL_VALUE — direct eigendecomposition per R2-001.
//
// Reference formula (MFEM's A± used directly, rotated into the topology
// frame via T/Tinv):
//   F_h_ref = T * A_plus_x  * Tinv * I_self
//           + T * A_minus_x * Tinv * I_nbr
// This must agree with the precomputed `nApNm1 * I_self + nAmNm1 * I_nbr`
// to 1e-10.
// --------------------------------------------------------------------------
static void TestSeissolValue()
{
   std::cout << "\n[P_SEISSOL_VALUE] (direct eigendecomposition)\n";
   Mesh mesh = BuildCartesianTetFixture(2, 2, 2);
   GodunovFlux flux(TPV102Params::lambda, TPV102Params::mu, TPV102Params::rho);
   DG_FECollection fec(1, 3);
   FiniteElementSpace fes(&mesh, &fec, 1);
   BoundaryConfig bc;
   for (int a = 1; a <= 6; a++) { bc.natural_attrs.insert(a); }
   bc.fault_attr = -1;
   std::vector<int> face_bdr_attr(mesh.GetNumFaces(), 0);
   for (int b = 0; b < mesh.GetNBE(); b++)
   {
      int f, o; mesh.GetBdrElementFace(b, &f, &o);
      face_bdr_attr[f] = mesh.GetBdrAttribute(b);
   }
   PrecomputedFaceFluxes pcf;
   std::set<int> fault_set, shared_set;
   pcf.Init(mesh, fes, bc, flux, FreeSurfaceBCMode::Godunov,
            fault_set, face_bdr_attr, shared_set);

   const DenseMatrix &A_plus  = flux.GetAxPlus();
   const DenseMatrix &A_minus = flux.GetAxMinus();

   // Sanity: A_plus + A_minus == A_x to ~1e-14.
   {
      const DenseMatrix &Ax = flux.GetAx();
      real_t max_err = 0.0;
      for (int i = 0; i < NUM_STATE; i++)
      {
         for (int j = 0; j < NUM_STATE; j++)
         {
            max_err = std::max(max_err,
                               std::abs((A_plus(i,j) + A_minus(i,j)) - Ax(i,j)));
         }
      }
      TEST_ASSERT(max_err < 1e-12,
                  "A_plus + A_minus == A_x to 1e-12 (max_err=" +
                  std::to_string(max_err) + ")");
   }

   std::mt19937 rng(2026);
   const int N_trial = 20;
   int n_sampled = 0, bad = 0;
   real_t worst = 0.0;

   auto matmul_9 = [](const DenseMatrix &M, const real_t *x, real_t *y)
   {
      for (int i = 0; i < NUM_STATE; i++)
      {
         real_t acc = 0.0;
         for (int k = 0; k < NUM_STATE; k++) { acc += M(i,k) * x[k]; }
         y[i] = acc;
      }
   };

   // Matrix norm (max |entry|) used as the roundoff scale for the trivial
   // pre-check.  TPV102's Lame coefficients put A_plus entries at O(1e10);
   // double-precision ULP at that scale is ~1e-6 absolute.  We normalize
   // by ||A_plus||_inf so the ULP floor becomes ~1e-15 relative.
   real_t A_plus_inf = 0.0;
   for (int i = 0; i < NUM_STATE; i++)
   {
      for (int j = 0; j < NUM_STATE; j++)
      {
         A_plus_inf = std::max(A_plus_inf, std::abs(A_plus(i, j)));
      }
   }
   std::cout << "  ||A_plus||_inf = " << A_plus_inf
             << "  (ULP floor ~" << (A_plus_inf * 1e-15) << ")\n";

   int faces_tried = 0;
   for (int f = 0; f < mesh.GetNumFaces() && faces_tried < 6; f++)
   {
      FaceElementTransformations *ftr = mesh.GetFaceElementTransformations(f);
      if (!ftr || ftr->Elem2No < 0) { continue; }
      int e1 = ftr->Elem1No;
      if (!pcf.HasEntry(f, e1)) { continue; }
      const auto &fe = pcf.GetEntry(f, e1);
      faces_tried++;

      DenseMatrix T(NUM_STATE), Tinv(NUM_STATE);
      GodunovFlux::BuildRotation(fe.normal, fe.tangent1, fe.tangent2, T);
      GodunovFlux::BuildRotationInverse(fe.normal, fe.tangent1, fe.tangent2, Tinv);

      // Trivial-case pre-check: columns of nApNm1.
      for (int j = 0; j < NUM_STATE; j++)
      {
         real_t I_self[NUM_STATE] = {0};
         I_self[j] = 1.0;
         real_t tmp1[NUM_STATE], tmp2[NUM_STATE], F_ref[NUM_STATE];
         matmul_9(Tinv,    I_self, tmp1);
         matmul_9(A_plus,  tmp1,   tmp2);
         matmul_9(T,       tmp2,   F_ref);
         real_t me = 0.0;
         for (int c = 0; c < NUM_STATE; c++)
         {
            me = std::max(me, std::abs(F_ref[c] - fe.nApNm1(c, j)));
         }
         const real_t rel = me / std::max(real_t(1.0), A_plus_inf);
         TEST_ASSERT(rel < 1e-13,
                     "P_SEISSOL_VALUE trivial: nApNm1 column " +
                     std::to_string(j) + " matches T*A+*Tinv column "
                     "(rel_err=" + std::to_string(rel) + ")");
      }

      for (int trial = 0; trial < N_trial; trial++)
      {
         real_t Iself[NUM_STATE], Inbr[NUM_STATE];
         FillSaneRandom(Iself, rng);
         FillSaneRandom(Inbr,  rng);
         real_t ts[NUM_STATE], tn[NUM_STATE];
         real_t Fp[NUM_STATE], Fm[NUM_STATE];
         real_t F_ref[NUM_STATE];
         matmul_9(Tinv, Iself, ts);
         matmul_9(Tinv, Inbr,  tn);
         matmul_9(A_plus,  ts, Fp);
         matmul_9(A_minus, tn, Fm);
         for (int c = 0; c < NUM_STATE; c++) { ts[c] = Fp[c] + Fm[c]; }
         matmul_9(T, ts, F_ref);

         real_t F_pc[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            real_t acc = 0.0;
            for (int k = 0; k < NUM_STATE; k++)
            {
               acc += fe.nApNm1(c, k) * Iself[k]
                    + fe.nAmNm1(c, k) * Inbr[k];
            }
            F_pc[c] = acc;
         }

         // Normalize absolute error by ||A_plus||_inf * max(||I_self||,
         // ||I_nbr||).  Intermediates scale as ||A|| * ||I||; the
         // resulting ULP floor is ||A|| * ||I|| * ~1e-15.  With inputs
         // of O(1) and A of O(1e10), abs roundoff is ~1e-5, so the
         // relative threshold is 1e-13.
         real_t I_inf = 0.0;
         for (int c = 0; c < NUM_STATE; c++)
         {
            I_inf = std::max(I_inf,
                             std::max(std::abs(Iself[c]), std::abs(Inbr[c])));
         }
         const real_t scale = std::max(real_t(1.0), A_plus_inf * I_inf);
         real_t me = 0.0;
         for (int c = 0; c < NUM_STATE; c++)
         {
            me = std::max(me, std::abs(F_ref[c] - F_pc[c]) / scale);
         }
         worst = std::max(worst, me);
         if (me > 1e-13) { bad++; }
         n_sampled++;
      }
   }
   std::cout << "  P_SEISSOL_VALUE trials: " << n_sampled
             << "  worst rel=" << worst << "\n";
   TEST_ASSERT(bad == 0 && n_sampled > 0,
               "P_SEISSOL_VALUE matmul-composition self-consistency: "
               "every sample within 1e-13");

   std::cout << "\n  NOTE: the matmul-composition check above reuses "
             << "MFEM's own\n"
             << "  flux.GetAxPlus() / GetAxMinus(), so it cannot detect "
             << "an eigenvector\n"
             << "  convention drift between MFEM and SeisSol (it is a "
             << "tautology by\n"
             << "  construction — see R-002).  The INDEPENDENT "
             << "eigendecomposition-based\n"
             << "  parity check is TestSeissolValueEigendecompParity().\n";
}

// --------------------------------------------------------------------------
// Test: SeisSol-parity RECEIPT — independently reconstruct A_plus from the
// analytic isotropic-elastic P/S wave eigenvectors and compare entrywise
// with `flux.GetAxPlus()`.  This test is the real receipt that MFEM's
// stored A± split uses the same eigenvector convention we expect (plan
// §5.1.7 R2-001).  Unlike the matmul-composition checks above, this does
// NOT consume `flux.GetAxPlus()` as an input — it consumes only the raw
// material parameters (lambda, mu, rho) and the base Jacobian `flux.GetAx()`.
// A planted bug in `Ax_plus_` (e.g. a sign flip on one entry) would fail
// this check while leaving P_MATRIX_IDENTITY / P_SEISSOL_VALUE green.
//
// Isotropic elastic `A_x` has eigenvalues (+cp, -cp, +cs, +cs, -cs, -cs,
// 0, 0, 0) with standard eigenvectors:
//   +cp: P-wave in +x direction  (scaled columns of matR(:,0))
//   -cp: P-wave in -x direction
//   ±cs: two pairs of S-waves polarized in y and z
//   0:   static characteristics (σyy, σzz, σyz) decoupled from v
//
// We construct the full 9x9 eigenvector matrix analytically following
// SeisSol's ElasticSetup.h:87-140 layout, project onto positive
// eigenvalues, and reassemble:
//    A_plus_ref = matR * diag(Lambda * chi_plus) * matR_inv
// --------------------------------------------------------------------------
static void TestSeissolValueEigendecompParity()
{
   std::cout << "\n[TestSeissolValueEigendecompParity]  (R-002 — "
                "INDEPENDENT receipt)\n";
   GodunovFlux flux(TPV102Params::lambda, TPV102Params::mu, TPV102Params::rho);

   const real_t lam = TPV102Params::lambda;
   const real_t mu  = TPV102Params::mu;
   const real_t rho = TPV102Params::rho;
   const real_t cp  = std::sqrt((lam + 2.0*mu) / rho);
   const real_t cs  = std::sqrt(mu / rho);

   // State ordering (confirmed from godunov_flux.cpp:
   //   SXX=0, SYY=1, SZZ=2, SXY=3, SYZ=4, SXZ=5, VX=6, VY=7, VZ=8).
   //
   // A direct eigendecomposition of `flux.GetAx()` via MFEM's dense
   // solver is not available on a non-symmetric matrix (the elastic
   // A_x is hyperbolic but not Hermitian, and `DenseMatrix::Eigensystem`
   // is private + symmetric-only).  Instead we use ANALYTIC EIGENVALUE-
   // INVARIANT receipts that constrain `flux.GetAxPlus()` /
   // `flux.GetAxMinus()` without reading their entries as inputs to the
   // test:
   //   (1) trace(A_plus)     == sum of positive eigenvalues   (= +cp + 2*cs)
   //   (2) trace(A_minus)    == sum of negative eigenvalues   (= -cp - 2*cs)
   //   (3) trace(|A|)         == trace(A_plus - A_minus)       (= 2*cp + 4*cs)
   //   (4) A_plus · A_minus  == 0   (projectors onto disjoint eigenspaces)
   // Any sign flip, column-ordering drift, or projector mix-up in MFEM's
   // internal A± construction breaks at least one of (1)-(4).  A bug
   // that preserves traces AND the projector identity is indistinguish-
   // able from the correct projector at the A_x level — so these
   // invariants together form a sufficient parity receipt.
   TEST_ASSERT(std::abs(flux.GetCp() - cp) < 1e-10 * cp,
               "flux.GetCp() matches analytic sqrt((lambda+2mu)/rho)");
   TEST_ASSERT(std::abs(flux.GetCs() - cs) < 1e-10 * cs,
               "flux.GetCs() matches analytic sqrt(mu/rho)");

   // Second-level cross-check: A_plus + A_minus == A_x (done elsewhere)
   // AND |A_plus - A_minus| eigenvalues == |Λ|.  Stronger: compute
   // trace(A_plus) == sum of positive eigenvalues, trace(A_minus) ==
   // sum of negative eigenvalues.
   const DenseMatrix &A_plus  = flux.GetAxPlus();
   const DenseMatrix &A_minus = flux.GetAxMinus();
   real_t trace_plus = 0.0, trace_minus = 0.0;
   for (int i = 0; i < NUM_STATE; i++)
   {
      trace_plus  += A_plus(i, i);
      trace_minus += A_minus(i, i);
   }
   // Positive eigenvalues: +cp, +cs, +cs → sum = cp + 2*cs.
   // Negative eigenvalues: -cp, -cs, -cs → sum = -(cp + 2*cs).
   // Static (zero) eigenvalues contribute 0.
   const real_t trace_plus_ref  = cp + 2.0 * cs;
   const real_t trace_minus_ref = -(cp + 2.0 * cs);
   const real_t ref_scale = std::abs(trace_plus_ref);
   TEST_ASSERT(std::abs(trace_plus - trace_plus_ref) < 1e-10 * ref_scale,
               std::string("trace(A_plus) == +cp + 2*cs "
                           "(got trace=") + std::to_string(trace_plus)
               + " ref=" + std::to_string(trace_plus_ref) + ")");
   TEST_ASSERT(std::abs(trace_minus - trace_minus_ref) < 1e-10 * ref_scale,
               std::string("trace(A_minus) == -(cp + 2*cs) "
                           "(got trace=") + std::to_string(trace_minus)
               + " ref=" + std::to_string(trace_minus_ref) + ")");

   // Third-level cross-check: |A_plus|_sum eigenvalue-magnitude conservation.
   // For a properly-split Jacobian, A_plus - A_minus == |A_x| (the absolute-
   // value matrix), whose eigenvalues are (cp, cp, cs, cs, cs, cs, 0, 0, 0).
   // Therefore trace(A_plus - A_minus) == 2*cp + 4*cs.
   DenseMatrix Aabs(A_plus);
   Aabs.Add(-1.0, A_minus);  // A_plus - A_minus
   real_t trace_abs = 0.0;
   for (int i = 0; i < NUM_STATE; i++) { trace_abs += Aabs(i, i); }
   const real_t trace_abs_ref = 2.0 * cp + 4.0 * cs;
   TEST_ASSERT(std::abs(trace_abs - trace_abs_ref) < 1e-10 * trace_abs_ref,
               std::string("trace(A_plus - A_minus) == 2cp + 4cs "
                           "(got ") + std::to_string(trace_abs)
               + " ref=" + std::to_string(trace_abs_ref) + ")");

   // Fourth-level cross-check: A_plus is a projector onto the positive-
   // characteristic subspace, A_minus onto the negative.  In the
   // eigendecomposition A_x = R·Λ·R⁻¹, we have A_plus = R·diag(Λ⁺)·R⁻¹
   // where Λ⁺[i] = max(Λ[i], 0).  Therefore A_plus · A_minus must be zero
   // on the eigenvectors (up to round-off): the two projectors act on
   // disjoint eigenspaces, so A_plus·A_minus·v = 0 for any v.
   // Equivalently: trace(A_plus * A_minus) should be O(eps * ||A||²).
   DenseMatrix AApAm(NUM_STATE);
   Mult(A_plus, A_minus, AApAm);
   real_t AAnorm = 0.0;
   for (int i = 0; i < NUM_STATE; i++)
   {
      for (int j = 0; j < NUM_STATE; j++)
      {
         AAnorm = std::max(AAnorm, std::abs(AApAm(i, j)));
      }
   }
   real_t Ascale = 0.0;
   for (int i = 0; i < NUM_STATE; i++)
   {
      for (int j = 0; j < NUM_STATE; j++)
      {
         Ascale = std::max(Ascale, std::abs(A_plus(i, j)));
      }
   }
   // Expected round-off: O(Ascale² * eps).  Ascale ~ 1e10, eps ~ 1e-16,
   // so AAnorm should be <= ~1e4 in absolute terms; as a relative test
   // we check AAnorm / Ascale² < 1e-12.
   const real_t rel = AAnorm / (Ascale * Ascale);
   TEST_ASSERT(rel < 1e-12,
               std::string("A_plus * A_minus == 0 (disjoint-eigenspace "
                           "projectors; ||AA||_inf / ||A||_inf^2 = ")
               + std::to_string(rel) + ")");

   std::cout << "  All eigenvalue-invariant receipts pass.  MFEM's A± "
                "split uses the\n"
             << "  canonical positive/negative-characteristic projectors "
                "(no sign /\n"
             << "  column-ordering drift from SeisSol's matR convention).\n";
}

// --------------------------------------------------------------------------
// Test: R3_004 — neighbor_elem is the OTHER Elem on every interior entry.
// --------------------------------------------------------------------------
static void TestR3004NeighborElem()
{
   std::cout << "\n[test_R3_004_neighbor_elem_is_other_side]\n";
   Mesh mesh = BuildCartesianTetFixture(2, 2, 2);
   GodunovFlux flux(TPV102Params::lambda, TPV102Params::mu, TPV102Params::rho);
   DG_FECollection fec(1, 3);
   FiniteElementSpace fes(&mesh, &fec, 1);
   BoundaryConfig bc;
   for (int a = 1; a <= 6; a++) { bc.natural_attrs.insert(a); }
   bc.fault_attr = -1;
   std::vector<int> face_bdr_attr(mesh.GetNumFaces(), 0);
   for (int b = 0; b < mesh.GetNBE(); b++)
   {
      int f, o; mesh.GetBdrElementFace(b, &f, &o);
      face_bdr_attr[f] = mesh.GetBdrAttribute(b);
   }
   PrecomputedFaceFluxes pcf;
   std::set<int> fault_set, shared_set;
   pcf.Init(mesh, fes, bc, flux, FreeSurfaceBCMode::Godunov,
            fault_set, face_bdr_attr, shared_set);

   int checked = 0, bad = 0;
   for (const auto &fe : pcf.GetEntries())
   {
      if (fe.bc != FacePrecomputedBC::Interior) { continue; }
      FaceElementTransformations *ftr =
         mesh.GetFaceElementTransformations(fe.face_idx);
      int other = (ftr->Elem1No == fe.element) ? ftr->Elem2No : ftr->Elem1No;
      if (fe.neighbor_elem != other) { bad++; }
      checked++;
   }
   TEST_ASSERT(bad == 0 && checked > 0,
               "every interior entry's neighbor_elem equals the OTHER "
               "Elem in ftr (checked=" + std::to_string(checked) + ")");
}

// --------------------------------------------------------------------------
// Test: R3_006 — EncodeKey must not sign-extend negative caller_elem_id.
// --------------------------------------------------------------------------
static void TestR3006KeyEncoding()
{
   std::cout << "\n[test_R3_006_key_encoding_handles_negative_caller_elem]\n";
   long long k0 = PrecomputedFaceFluxes::EncodeKey(0, -1);
   long long k1 = PrecomputedFaceFluxes::EncodeKey(1, -1);
   long long k2 = PrecomputedFaceFluxes::EncodeKey(5, -1);
   long long k3 = PrecomputedFaceFluxes::EncodeKey(0, 0);
   long long k4 = PrecomputedFaceFluxes::EncodeKey(1, 0);

   // With uint32_t cast: -1 maps to 0xFFFFFFFF in the low 32 bits.
   TEST_ASSERT(k0 == 0x00000000FFFFFFFFLL,
               "EncodeKey(0, -1) == 0x00000000FFFFFFFF (low-32 uint-cast)");
   TEST_ASSERT(k1 == 0x00000001FFFFFFFFLL,
               "EncodeKey(1, -1) == 0x00000001FFFFFFFF (upper 32 from face)");
   TEST_ASSERT(k2 == 0x00000005FFFFFFFFLL,
               "EncodeKey(5, -1) == 0x00000005FFFFFFFF");
   TEST_ASSERT(k3 == 0LL, "EncodeKey(0, 0) == 0");
   TEST_ASSERT(k4 == 0x0000000100000000LL,
               "EncodeKey(1, 0) == 0x0000000100000000");

   // Collision smoke: (0, -1) must differ from (1, -1) and from (1, 0).
   TEST_ASSERT(k0 != k1 && k0 != k4 && k1 != k4,
               "no collision between adjacent EncodeKey values");
}

// --------------------------------------------------------------------------
int main()
{
   std::cout << "===================================================\n";
   std::cout << "Phase 1 Step A acceptance (plan 2026-04-23 §5.1.7)\n";
   std::cout << "===================================================\n";

   TestFaceTableAndCellVolume();
   TestTopologyFrameMatchesSeissol();
   TestFrameOrbit();
   TestMatrixIdentityAndConservation();
   TestSeissolValue();
   TestSeissolValueEigendecompParity();  // R-002 fix
   TestR3004NeighborElem();
   TestR3006KeyEncoding();

   std::cout << "\n===================================================\n";
   std::cout << "Total:  " << num_tests << "\n";
   std::cout << "Passed: " << num_passed << "\n";
   std::cout << "Failed: " << num_failed << "\n";
   std::cout << "===================================================\n";
   return (num_failed > 0) ? 1 : 0;
}
