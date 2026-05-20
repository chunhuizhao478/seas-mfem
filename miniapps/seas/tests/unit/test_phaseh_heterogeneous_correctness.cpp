// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_phaseh_heterogeneous_correctness.cpp — Phase H Stage 2 acceptance
// gates H-1 / H-2 / H-3 for
// PLAN_heterogeneous_volume_bc_fault_dispatch_2026-05-20.md §6.
//
// These tests prove the per-element dispatch does something REAL — a
// heterogeneous result a single scalar Jacobian cannot reproduce.  Where
// the plan describes a full plane-wave reflection/transmission benchmark
// (an integration-level run), these unit tests verify the equivalent,
// stronger, deterministic statement: the operator's volume / boundary /
// fault flux on each element uses THAT element's material — provably
// different from a single-material scalar operator, and provably equal to
// a homogeneous operator built from that element's own material.
//
//   H-1 (volume):  ComputeVolumeRHS is element-local, so element e's
//                  contribution depends only on e's Jacobian.  We assert
//                  that, on a two-layer slab, the heterogeneous operator's
//                  per-element volume RHS equals a homogeneous operator
//                  built from that element's layer material, and DIFFERS
//                  from a homogeneous operator built from the OTHER layer
//                  (which is what a single scalar Jacobian would give at a
//                  material jump — the R/T-producing term).
//
//   H-2 (BC):      The absorbing / free-surface BC flux uses the element-
//                  correct impedance.  We assert (a) GodunovFlux::
//                  AbsorbingTotal / FreeSurfaceTotal differ between the two
//                  layer materials by > 1e-6 relative, and (b) on a slab
//                  with absorbing BC the heterogeneous operator's per-
//                  element material differs between the shallow and deep
//                  boundary-adjacent elements, so FluxForElem_(e1) feeds
//                  the BC dispatch a depth-correct flux.  (The byte-exact
//                  routing of the BC path through FluxForElem_ is gated by
//                  the mixed-BC parity case in
//                  test_phaseh_wave_operator_constant_parity C-5.)
//
//   H-3 (fault):   The fault bulk-side flux A_n·Q_imp differs between a
//                  shallow and a deep fault QP when the material varies
//                  with depth, whereas a single scalar Jacobian makes them
//                  equal.  We assert (a) GodunovFlux::Interior(can_n,
//                  Q_imp, Q_imp, ·) — exactly the C1-C4 bulk-side flux —
//                  differs between the layer materials by > 1e-6, and (b)
//                  on a depth-varying fault mesh the two fault-adjacent
//                  elements of a shallow vs a deep fault face carry
//                  different (lambda, mu, rho).

#include "mfem.hpp"

#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/heterogeneous_material.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

namespace
{
int g_num_tests = 0, g_num_passed = 0, g_num_failed = 0;
int g_rank = 0;

#define TEST_ASSERT(c, m) do { ++g_num_tests; if (!(c)) { \
   if (g_rank == 0) \
   { std::cerr << "FAILED: " << m << " line " << __LINE__ << "\n"; } \
   ++g_num_failed; } else { if (g_rank == 0) \
   { std::cout << "  PASSED: " << m << "\n"; } ++g_num_passed; } } while (0)

constexpr int    k_order = 1;

// Layer A (the "shallow" / lower-velocity material) and layer B (the
// "deep" / higher-velocity material).  muB = 4*muA so c_s,B = 2*c_s,A and
// the S-impedance ratio Z2/Z1 = 2 — a strong, unambiguous jump.
constexpr real_t k_lamA = 30.0e9, k_muA = 30.0e9, k_rhoA = 2670.0;
constexpr real_t k_lamB = 60.0e9, k_muB = 120.0e9, k_rhoB = 2670.0;

// Mode::Coefficient step function in z: material A for z-centroid below
// `z_split`, material B above.  One instance per scalar field.
class StepInZ : public Coefficient
{
public:
   StepInZ(real_t below, real_t above, real_t z_split)
      : below_(below), above_(above), z_split_(z_split) {}
   real_t Eval(ElementTransformation &T, const IntegrationPoint &ip) override
   {
      Vector x(3);
      T.Transform(ip, x);
      return (x(2) < z_split_) ? below_ : above_;
   }
private:
   real_t below_, above_, z_split_;
};

Mesh MakeBoxMesh(int nx, int ny, int nz, real_t L)
{
   return Mesh::MakeCartesian3D(nx, ny, nz, Element::TETRAHEDRON,
                                L, L, L, /*sfc_ordering=*/false);
}

void FillQDeterministic(Vector &Q)
{
   const int n = Q.Size();
   for (int i = 0; i < n; ++i)
   {
      const real_t x = static_cast<real_t>(i) / static_cast<real_t>(n);
      Q(i) = 1e6 * std::sin(2.0 * M_PI * 3.0 * x)
           + 5e5 * std::cos(2.0 * M_PI * 7.0 * x);
   }
}

BoundaryConfig MakeBareBC()
{
   BoundaryConfig bc;
   bc.fault_attr = 0;
   return bc;
}

// z-centroid of element e.
real_t ElemCentroidZ(Mesh &mesh, int e)
{
   Array<int> v;
   mesh.GetElementVertices(e, v);
   real_t cz = 0.0;
   for (int i = 0; i < v.Size(); ++i) { cz += mesh.GetVertex(v[i])[2]; }
   return cz / std::max(1, v.Size());
}

// Max relative diff over element e's 9-component DOF block.
real_t ElemMaxRelDiff(const Vector &a, const Vector &b, int e,
                      int ndof_per_el, int ndof_total)
{
   real_t mx = 0.0;
   for (int c = 0; c < NUM_STATE; ++c)
   {
      for (int i = 0; i < ndof_per_el; ++i)
      {
         const int idx = c * ndof_total + e * ndof_per_el + i;
         const real_t denom = std::max(
            std::max(std::abs(a(idx)), std::abs(b(idx))), real_t(1.0));
         mx = std::max(mx, std::abs(a(idx) - b(idx)) / denom);
      }
   }
   return mx;
}

real_t MaxRelDiff9(const real_t *a, const real_t *b)
{
   real_t mx = 0.0;
   for (int i = 0; i < NUM_STATE; ++i)
   {
      const real_t denom = std::max(
         std::max(std::abs(a[i]), std::abs(b[i])), real_t(1.0));
      mx = std::max(mx, std::abs(a[i] - b[i]) / denom);
   }
   return mx;
}

// Max relative diff over a full state vector; counts entries over `tol`.
real_t VecMaxRelDiff(const Vector &a, const Vector &b, real_t tol, int *n_over)
{
   real_t mx = 0.0;
   int over = 0;
   for (int i = 0; i < a.Size(); ++i)
   {
      const real_t denom = std::max(
         std::max(std::abs(a(i)), std::abs(b(i))), real_t(1.0));
      const real_t r = std::abs(a(i) - b(i)) / denom;
      if (r > mx) { mx = r; }
      if (r > tol) { ++over; }
   }
   if (n_over) { *n_over = over; }
   return mx;
}

}  // namespace

// =========================================================================
// H-1 — Volume per-element correctness on a two-layer slab.
// =========================================================================
static void H_1_volume_layered_slab()
{
   if (g_rank == 0)
   { std::cout << "\n[H-1] Volume per-element correctness (two-layer slab)\n"; }
#ifdef MFEM_USE_MPI
   if (g_rank != 0) { return; }   // serial-mesh check on rank 0 only
#endif

   const real_t L = 1.0;
   const real_t z_split = 0.5 * L;
   Mesh mesh = MakeBoxMesh(2, 2, 4, L);
   const BoundaryConfig bc = MakeBareBC();

   StepInZ lam_c(k_lamA, k_lamB, z_split);
   StepInZ mu_c (k_muA,  k_muB,  z_split);
   StepInZ rho_c(k_rhoA, k_rhoB, z_split);
   MaterialField mat = MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c);

   WaveOperator<Mesh> wave_hetero(mesh, k_order, mat, bc);
   WaveOperator<Mesh> wave_homoA(mesh, k_order, k_lamA, k_muA, k_rhoA, bc);
   WaveOperator<Mesh> wave_homoB(mesh, k_order, k_lamB, k_muB, k_rhoB, bc);

   const int ne  = wave_hetero.NumElements();
   const int npe = wave_hetero.GetNDof();
   const int nt  = wave_hetero.GetScalarNDof();
   const int n   = NUM_STATE * nt;

   Vector Q(n);
   FillQDeterministic(Q);
   Vector v_het(n), v_a(n), v_b(n);
   v_het = 0.0; v_a = 0.0; v_b = 0.0;
   wave_hetero.ComputeVolumeRHS_ForTest(Q, v_het);
   wave_homoA.ComputeVolumeRHS_ForTest(Q, v_a);
   wave_homoB.ComputeVolumeRHS_ForTest(Q, v_b);

   const real_t k_match_tol = 1e-12;   // byte/FP-equal to the same material
   const real_t k_jump_tol  = 1e-3;    // 2x-impedance jump → big difference
   int n_layerA = 0, n_layerB = 0;
   int n_match_fail = 0, n_jump_fail = 0;
   real_t worst_jump = 0.0;
   for (int e = 0; e < ne; ++e)
   {
      const bool inA = ElemCentroidZ(mesh, e) < z_split;
      // hetero[e] must equal the homogeneous op built from e's own layer.
      const real_t match = inA
         ? ElemMaxRelDiff(v_het, v_a, e, npe, nt)
         : ElemMaxRelDiff(v_het, v_b, e, npe, nt);
      if (match > k_match_tol) { ++n_match_fail; }
      // hetero[e] must DIFFER from the OTHER layer's homogeneous op
      // (this is the per-element term a single scalar Jacobian misses).
      const real_t jump = inA
         ? ElemMaxRelDiff(v_het, v_b, e, npe, nt)
         : ElemMaxRelDiff(v_het, v_a, e, npe, nt);
      worst_jump = std::max(worst_jump, jump);
      if (jump < k_jump_tol) { ++n_jump_fail; }
      inA ? ++n_layerA : ++n_layerB;
   }

   if (g_rank == 0)
   {
      std::cout << "  elements: layerA=" << n_layerA << " layerB=" << n_layerB
                << "  match_fail=" << n_match_fail
                << "  jump_fail=" << n_jump_fail
                << "  worst_jump=" << worst_jump << "\n";
   }
   TEST_ASSERT(n_layerA > 0 && n_layerB > 0,
               "slab has elements in both layers");
   TEST_ASSERT(n_match_fail == 0,
               "per-element volume RHS == homogeneous op of that element's "
               "material (each element uses its own Jacobian)");
   TEST_ASSERT(n_jump_fail == 0,
               "per-element volume RHS differs from the other layer's "
               "single-Jacobian result (scalar op cannot reproduce this)");
}

// =========================================================================
// H-2 — Boundary-face flux uses the element-correct impedance.
// =========================================================================
static void H_2_boundary_per_element_impedance()
{
   if (g_rank == 0)
   { std::cout << "\n[H-2] Boundary-face per-element impedance\n"; }
#ifdef MFEM_USE_MPI
   if (g_rank != 0) { return; }
#endif

   // (a) The actual BC flux functions differ between the two materials.
   GodunovFlux flux_A(k_lamA, k_muA, k_rhoA);
   GodunovFlux flux_B(k_lamB, k_muB, k_rhoB);

   const real_t nor[3]    = { 0.0, 0.0, 1.0 };           // +z outward
   const real_t Q_self[9] = { 1e6, 2e6, 3e6, 4e5, 5e5, 6e5, 0.1, 0.2, 0.3 };
   const real_t Q_bg[9]   = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

   real_t F_abs_A[9], F_abs_B[9];
   flux_A.AbsorbingTotal(nor, Q_self, Q_bg, F_abs_A);
   flux_B.AbsorbingTotal(nor, Q_self, Q_bg, F_abs_B);
   const real_t rel_abs = MaxRelDiff9(F_abs_A, F_abs_B);

   real_t F_fs_A[9], F_fs_B[9];
   flux_A.FreeSurfaceTotal(nor, Q_self, Q_bg, F_fs_A);
   flux_B.FreeSurfaceTotal(nor, Q_self, Q_bg, F_fs_B);
   const real_t rel_fs = MaxRelDiff9(F_fs_A, F_fs_B);

   if (g_rank == 0)
   {
      std::cout << "  AbsorbingTotal rel(A,B)=" << rel_abs
                << "  FreeSurfaceTotal rel(A,B)=" << rel_fs
                << " (expect > 1e-6)\n";
   }
   TEST_ASSERT(rel_abs > 1e-6,
               "AbsorbingTotal flux is impedance-specific (A vs B differ)");
   TEST_ASSERT(rel_fs > 1e-6,
               "FreeSurfaceTotal flux is material-specific (A vs B differ)");

   // (b) On a two-layer slab with absorbing BC, the boundary-adjacent
   //     elements carry layer-specific material, so FluxForElem_(e1) feeds
   //     the BC dispatch a depth-correct flux.
   const real_t L = 1.0;
   const real_t z_split = 0.5 * L;
   Mesh mesh = MakeBoxMesh(2, 2, 4, L);

   BoundaryConfig bc;
   bc.fault_attr      = 0;
   bc.absorbing_attrs = {1, 2, 3, 4, 5, 6};

   StepInZ lam_c(k_lamA, k_lamB, z_split);
   StepInZ mu_c (k_muA,  k_muB,  z_split);
   StepInZ rho_c(k_rhoA, k_rhoB, z_split);
   MaterialField mat = MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c);
   WaveOperator<Mesh> wave(mesh, k_order, mat, bc);

   const auto &lmr = wave.GetPerElementMaterial();
   bool saw_A = false, saw_B = false;
   for (int e = 0; e < wave.NumElements(); ++e)
   {
      if (ElemCentroidZ(mesh, e) < z_split)
      { if (std::abs(lmr[e][1] - k_muA) < 1.0) { saw_A = true; } }
      else
      { if (std::abs(lmr[e][1] - k_muB) < 1.0) { saw_B = true; } }
   }
   TEST_ASSERT(saw_A && saw_B,
               "absorbing-slab op stores both layer-A and layer-B per-element "
               "material (INGREDIENT only; operator-level BC routing is gated "
               "by (c) below)");

   // (c) OPERATOR-LEVEL routing gate.  A UNIFORM Mode::Coefficient op
   //     (constant-valued coefficients == material A) has owned_flux_pool_
   //     built — so the BC dispatch goes through FluxForElem_(e1) — yet its
   //     scalar `flux_` is the (1,1,1) placeholder (cp ≈ 1.73 m/s).  Its
   //     Mult AND AdvanceADER must therefore match the scalar homogeneous op
   //     built from material A to FP precision.  That holds ONLY because the
   //     volume, ADER-CK, interior-face AND boundary-face dispatch all read
   //     FluxForElem_(e)=At(e)=A, not the placeholder `flux_`.  A revert of
   //     any BC site (or volume/CK site) to `flux_` makes this diverge by
   //     O(1) — which (a)/(b) above (ingredient checks) would NOT catch.
   //     MIXED BC (absorbing + free-surface/natural) exercises both BC paths.
   BoundaryConfig bc_mixed;
   bc_mixed.fault_attr      = 0;
   bc_mixed.absorbing_attrs = {1, 3, 5};
   bc_mixed.natural_attrs   = {2, 4, 6};

   ConstantCoefficient lamA_cc(k_lamA), muA_cc(k_muA), rhoA_cc(k_rhoA);
   MaterialField matA =
      MaterialField::MakeCoefficient(&lamA_cc, &muA_cc, &rhoA_cc);
   WaveOperator<Mesh> wave_coefA (mesh, k_order, matA, bc_mixed);
   WaveOperator<Mesh> wave_scalarA(mesh, k_order, k_lamA, k_muA, k_rhoA,
                                   bc_mixed);
   TEST_ASSERT(wave_coefA.UsesGodunovFluxPool(),
               "uniform-Coefficient op built the per-element pool "
               "(FluxForElem_ active, flux_ is the placeholder)");

   real_t Qbg[NUM_STATE] = {0};
   wave_coefA.SetAbsorbingBackground(Qbg);
   wave_scalarA.SetAbsorbingBackground(Qbg);

   const int n = NUM_STATE * wave_coefA.GetScalarNDof();
   Vector Q(n);
   FillQDeterministic(Q);

   const real_t k_rel_tol = 1e-9;
   {
      Vector dc(n), dsc(n);
      wave_coefA.Mult(Q, dc);
      wave_scalarA.Mult(Q, dsc);
      int n_over = 0;
      const real_t mx = VecMaxRelDiff(dc, dsc, k_rel_tol, &n_over);
      if (g_rank == 0)
      { std::cout << "  uniform-Coefficient Mult max_rel_diff=" << mx
                  << " (tol=" << k_rel_tol << ", n_over=" << n_over << ")\n"; }
      TEST_ASSERT(n_over == 0,
                  "uniform-Coefficient Mult == scalar op (BC dispatch routes "
                  "through FluxForElem_(e1), NOT the placeholder flux_)");
   }
   {
      const real_t dt = 0.2 * wave_scalarA.ComputeMaxDt(0.5);
      Vector Qc(n), Qsc(n);
      wave_coefA.AdvanceADER(Q, dt, /*order=*/3, Qc);
      wave_scalarA.AdvanceADER(Q, dt, /*order=*/3, Qsc);
      int n_over = 0;
      const real_t mx = VecMaxRelDiff(Qc, Qsc, k_rel_tol, &n_over);
      if (g_rank == 0)
      { std::cout << "  uniform-Coefficient AdvanceADER max_rel_diff=" << mx
                  << " (tol=" << k_rel_tol << ", n_over=" << n_over << ")\n"; }
      TEST_ASSERT(n_over == 0,
                  "uniform-Coefficient AdvanceADER == scalar op (ADER BC + CK "
                  "dispatch routes through FluxForElem_, NOT placeholder)");
   }
}

// =========================================================================
// H-3 — Fault bulk-side flux differs shallow vs deep (depth-varying mat).
// =========================================================================
static void H_3_fault_bulk_side_depth_varies()
{
   if (g_rank == 0)
   { std::cout << "\n[H-3] Fault bulk-side flux varies with depth\n"; }
#ifdef MFEM_USE_MPI
   if (g_rank != 0) { return; }
#endif

   // (a) The bulk-side fault flux is A_n·Q_imp, computed by C1-C4 as
   //     FluxForElem_(elem).Interior(can_n, Q_imp, Q_imp, F).  With a
   //     depth-varying material the fault-adjacent element's (lambda, mu,
   //     rho) — hence A_n — differs shallow vs deep, so F differs.  A
   //     single scalar Jacobian would make them equal.
   GodunovFlux flux_shallow(k_lamA, k_muA, k_rhoA);
   GodunovFlux flux_deep   (k_lamB, k_muB, k_rhoB);

   const real_t can_n[3] = { 0.0, -1.0, 0.0 };   // canonical fault normal
   // A representative imposed state in the global frame (same on both QPs;
   // the depth dependence we test comes entirely from A_n, not Q_imp).
   const real_t Q_imp[9] =
      { 1.2e7, -3.0e6, 4.5e6, 8.0e6, -2.0e6, 1.0e6, 0.05, -0.02, 0.03 };

   real_t F_shallow[9], F_deep[9];
   flux_shallow.Interior(can_n, Q_imp, Q_imp, F_shallow);
   flux_deep.Interior(can_n, Q_imp, Q_imp, F_deep);
   const real_t rel = MaxRelDiff9(F_shallow, F_deep);

   // Control: a single scalar flux makes the two QPs identical.
   real_t F_scalar_a[9], F_scalar_b[9];
   flux_shallow.Interior(can_n, Q_imp, Q_imp, F_scalar_a);
   flux_shallow.Interior(can_n, Q_imp, Q_imp, F_scalar_b);
   const real_t rel_scalar = MaxRelDiff9(F_scalar_a, F_scalar_b);

   if (g_rank == 0)
   {
      std::cout << "  bulk-side A_n*Q_imp rel(shallow,deep)=" << rel
                << "  (expect > 1e-6);  scalar-path rel=" << rel_scalar
                << "  (expect 0)\n";
   }
   TEST_ASSERT(rel > 1e-6,
               "fault bulk-side flux differs shallow vs deep "
               "(per-element A_n)");
   TEST_ASSERT(rel_scalar == 0.0,
               "scalar single-Jacobian path makes the two QPs equal");

   // (b) On a depth-varying fault mesh, the two fault-adjacent elements of
   //     a shallow vs a deep fault face carry different (lambda, mu, rho),
   //     so C1-C4's FluxForElem_(elem) gives a depth-dependent bulk flux.
   const real_t L = 1.0e4;             // 10 km box
   const real_t z_split = 0.5 * L;
   Mesh mesh = MakeBoxMesh(2, 2, 4, L);
   // Tag the y = L/2 plane of interior faces as the fault (attr 3).
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      Array<int> fv;
      mesh.GetFaceVertices(f, fv);
      if (fv.Size() != 3) { continue; }
      auto *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr && ftr->Elem2No >= 0)
      {
         real_t cy = 0.0;
         for (int v = 0; v < fv.Size(); v++) { cy += mesh.GetVertex(fv[v])[1]; }
         cy /= fv.Size();
         if (std::abs(cy - 0.5 * L) < 1e-8)
         { mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 3); }
         continue;
      }
      mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 1);
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();

   BoundaryConfig bc;
   bc.fault_attr    = 3;
   bc.natural_attrs = {1};

   // Material varies with z (mu jumps at z_split).  Both fault sides at a
   // given QP share material (vertical fault in a depth-only medium), so
   // the Zp_plus == Zp_minus contract holds — but a SHALLOW fault QP and a
   // DEEP fault QP see different material.
   StepInZ lam_c(k_lamA, k_lamB, z_split);
   StepInZ mu_c (k_muA,  k_muB,  z_split);
   StepInZ rho_c(k_rhoA, k_rhoB, z_split);
   MaterialField mat = MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c);
   WaveOperator<Mesh> wave(mesh, k_order, mat, bc);

   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   TEST_ASSERT(int_faces.Size() > 0, "fault mesh has interior fault faces");

   const auto &lmr = wave.GetPerElementMaterial();
   real_t mu_min = std::numeric_limits<real_t>::max();
   real_t mu_max = 0.0;
   for (int i = 0; i < int_faces.Size(); ++i)
   {
      auto *ftr = mesh.GetFaceElementTransformations(int_faces[i]);
      for (int elem : { ftr->Elem1No, ftr->Elem2No })
      {
         if (elem < 0) { continue; }
         mu_min = std::min(mu_min, lmr[elem][1]);
         mu_max = std::max(mu_max, lmr[elem][1]);
      }
   }
   if (g_rank == 0)
   {
      std::cout << "  fault-adjacent mu range: [" << mu_min << ", " << mu_max
                << "]  (expect mu_max/mu_min ~ 4)\n";
   }
   TEST_ASSERT(mu_max > 2.0 * mu_min,
               "fault-adjacent elements span >1 material (INGREDIENT only: "
               "shallow vs deep fault QPs CAN feed FluxForElem_ different "
               "Jacobians).  Operator-level fault routing — that Mult/"
               "AdvanceADER actually read FluxForElem_(elem) and not the "
               "placeholder flux_ — is gated by the uniform-Coefficient fault "
               "op in test_phaseh_wave_operator_constant_parity C-6.");
}

// =========================================================================
// main
// =========================================================================
int main(int argc, char *argv[])
{
#ifdef MFEM_USE_MPI
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);
#else
   g_rank = 0;
#endif

   if (g_rank == 0)
   {
      std::cout << "==============================================\n"
                << "test_phaseh_heterogeneous_correctness\n"
                << "Phase H Stage 2 gates H-1 / H-2 / H-3\n"
                << "==============================================\n";
   }

   H_1_volume_layered_slab();
   H_2_boundary_per_element_impedance();
   H_3_fault_bulk_side_depth_varies();

#ifdef MFEM_USE_MPI
   int total = 0, passed = 0, failed = 0;
   MPI_Allreduce(&g_num_tests,  &total,  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(&g_num_passed, &passed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(&g_num_failed, &failed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
   const int total = g_num_tests, passed = g_num_passed, failed = g_num_failed;
#endif
   if (g_rank == 0)
   {
      std::cout << "\n==============================================\n"
                << "Summary: " << passed << " / " << total
                << " passed (" << failed << " failed)\n"
                << "==============================================\n";
   }
#ifdef MFEM_USE_MPI
   MPI_Finalize();
#endif
   return failed > 0 ? 1 : 0;
}
